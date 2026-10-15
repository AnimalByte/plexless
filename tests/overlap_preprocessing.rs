mod common;
use common::{TestDir, read_gzip_text};
use std::process::Command;
const ADAPTER1: &str = "ACGTTGCACTGATCGAGGTACCTATG";
const ADAPTER2: &str = "GCTACAGTTCGATACCGTGACATGAC";
fn dna(n: usize) -> String {
    let mut x = 84729u64;
    (0..n)
        .map(|_| {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1);
            b"ACGT"[(x >> 32) as usize % 4] as char
        })
        .collect()
}
fn rc(s: &str) -> String {
    s.bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            _ => 'N',
        })
        .collect()
}
fn fq(id: &str, s: &str) -> String {
    format!("@{id}\n{s}\n+\n{}\n", "I".repeat(s.len()))
}
fn run(
    test: &TestDir,
    label: &str,
    a: &str,
    b: &str,
    extra: &[&str],
    threads: usize,
) -> std::path::PathBuf {
    let out = test.child(&format!("{label}-{threads}"));
    let result = Command::new(env!("CARGO_BIN_EXE_plexless"))
        .args(["trim", "--r1"])
        .arg(test.write(&format!("{label}-1.fq"), &fq("r/1", a)))
        .arg("--r2")
        .arg(test.write(&format!("{label}-2.fq"), &fq("r/2", b)))
        .arg("--output")
        .arg(&out)
        .arg("--threads")
        .arg(threads.to_string())
        .args([
            "--r1-distal-sequence",
            ADAPTER1,
            "--r2-distal-sequence",
            ADAPTER2,
        ])
        .args(extra)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    out
}
#[test]
fn default_confirms_real_adapter_boundaries_but_preserves_internal_biological_matches() {
    let test = TestDir::new("overlap-real-boundaries");
    for threads in [1, 4, 8] {
        let insert = dna(120);
        let a = format!("{insert}{ADAPTER1}");
        let b = format!("{}{ADAPTER2}", rc(&insert));
        let out = run(&test, "full", &a, &b, &[], threads);
        assert_eq!(
            read_gzip_text(&out.join("trimmed_R1.fastq.gz")),
            fq("r/1", &insert)
        );
        assert_eq!(
            read_gzip_text(&out.join("trimmed_R2.fastq.gz")),
            fq("r/2", &rc(&insert))
        );
        let out = run(
            &test,
            "partial",
            &format!("{insert}{}", &ADAPTER1[..20]),
            &format!("{}{}", rc(&insert), &ADAPTER2[..20]),
            &["--partial-pass-through"],
            threads,
        );
        assert_eq!(
            read_gzip_text(&out.join("trimmed_R1.fastq.gz")),
            fq("r/1", &insert)
        );
        let mut insert = dna(220);
        insert.replace_range(60..60 + ADAPTER1.len(), ADAPTER1);
        let a = &insert[..150];
        let b = rc(&insert[70..]);
        let out = run(&test, "internal", a, &b, &[], threads);
        assert_eq!(
            read_gzip_text(&out.join("trimmed_R1.fastq.gz")),
            fq("r/1", a)
        );
        let out = run(&test, "sequence-only", a, &b, &["--no-overlap"], threads);
        assert_eq!(
            read_gzip_text(&out.join("trimmed_R1.fastq.gz")),
            fq("r/1", &a[..60])
        );
    }
}
#[test]
fn rejected_early_full_candidate_does_not_hide_later_confirmed_full_boundary() {
    let test = TestDir::new("overlap-later-boundary");
    let mut insert = dna(120);
    insert.replace_range(40..40 + ADAPTER1.len(), ADAPTER1);
    let out = run(
        &test,
        "later",
        &format!("{insert}{ADAPTER1}"),
        &format!("{}{ADAPTER2}", rc(&insert)),
        &[],
        4,
    );
    assert_eq!(
        read_gzip_text(&out.join("trimmed_R1.fastq.gz")),
        fq("r/1", &insert)
    );
}
#[test]
fn unequal_mates_and_explicit_front_clipping_use_correct_physical_coordinates() {
    let test = TestDir::new("overlap-asymmetric");
    let insert = dna(120);
    let a = format!("{insert}{ADAPTER1}");
    let b = rc(&insert[20..]);
    let out = run(&test, "unequal", &a, &b, &[], 4);
    assert_eq!(
        read_gzip_text(&out.join("trimmed_R1.fastq.gz")),
        fq("r/1", &insert)
    );
    assert_eq!(
        read_gzip_text(&out.join("trimmed_R2.fastq.gz")),
        fq("r/2", &b)
    );
    let out = run(
        &test,
        "front",
        &a,
        &format!("{}{ADAPTER2}", rc(&insert)),
        &["--trim-front1", "3", "--trim-front2", "7"],
        4,
    );
    assert_eq!(
        read_gzip_text(&out.join("trimmed_R1.fastq.gz")),
        fq("r/1", &insert[3..])
    );
    assert_eq!(
        read_gzip_text(&out.join("trimmed_R2.fastq.gz")),
        fq("r/2", &rc(&insert)[7..])
    );
}
#[test]
fn short_unconfirmed_mates_stay_intact_and_full_known_base_floor_is_configurable() {
    let test = TestDir::new("overlap-short");
    let insert = dna(20);
    let a = format!("{insert}{ADAPTER1}");
    let b = format!("{}{ADAPTER2}", rc(&insert));
    let out = run(&test, "short", &a, &b, &[], 4);
    assert_eq!(
        read_gzip_text(&out.join("trimmed_R1.fastq.gz")),
        fq("r/1", &a)
    );
    let adapter = &ADAPTER1[..16];
    let reads = test.write("weak.fq", &fq("read", &format!("GATTACA{adapter}")));
    for (threshold, expected) in [(20, format!("GATTACA{adapter}")), (16, "GATTACA".into())] {
        let out = test.child(&format!("floor-{threshold}"));
        let result = Command::new(env!("CARGO_BIN_EXE_plexless"))
            .args(["trim", "--reads"])
            .arg(&reads)
            .arg("--output")
            .arg(&out)
            .args([
                "--r1-distal-sequence",
                adapter,
                "--min-full-adapter-known-bases",
                &threshold.to_string(),
            ])
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        assert_eq!(
            read_gzip_text(&out.join("trimmed.fastq.gz")),
            fq("read", &expected)
        );
    }
}
