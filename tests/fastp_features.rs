mod common;
use common::{TestDir, read_gzip_text};
use std::io::Write;
use std::process::{Command, Stdio};
fn fq(name: &str, seq: &str, quality: Option<&str>) -> String {
    format!(
        "@{name}\n{seq}\n+\n{}\n",
        quality
            .map(str::to_owned)
            .unwrap_or_else(|| "I".repeat(seq.len()))
    )
}
fn dna(n: usize) -> String {
    let mut x = 783392u64;
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
fn se(test: &TestDir, label: &str, text: &str, flags: &[&str]) -> std::path::PathBuf {
    let out = test.child(label);
    let result = Command::new(env!("CARGO_BIN_EXE_plexless"))
        .args(["trim", "--reads"])
        .arg(test.write(&format!("{label}.fq"), text))
        .arg("--output")
        .arg(&out)
        .args(flags)
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
fn quality_polyg_polyx_and_plain_output() {
    let test = TestDir::new("fastp-quality");
    let out = se(
        &test,
        "quality",
        &fq("r", "GATTACAACGT", Some("IIIIIII!!!!")),
        &["--cut-tail", "--cut-window-size", "1"],
    );
    assert_eq!(
        read_gzip_text(&out.join("trimmed.fastq.gz")),
        fq("r", "GATTACA", None)
    );
    let out = se(
        &test,
        "poly",
        &fq("r", "ACCTGGGGGGGGGG", None),
        &["--trim-poly-g"],
    );
    assert_eq!(
        read_gzip_text(&out.join("trimmed.fastq.gz")),
        fq("r", "ACCT", None)
    );
    let path = test.child("custom.fq");
    let out = se(
        &test,
        "plain",
        &fq("r", "GATTACA", None),
        &["--out1", path.to_str().unwrap()],
    );
    assert_eq!(
        std::fs::read_to_string(path).unwrap(),
        fq("r", "GATTACA", None)
    );
    assert!(!out.join("sample_metrics.tsv").exists());
}
#[test]
fn umi_headers_and_phred64_are_processed_in_shared_pipeline() {
    let test = TestDir::new("fastp-umi");
    let out = se(
        &test,
        "umi",
        &fq("r/1 extra", "ACGATTACA", None),
        &[
            "--umi",
            "--umi-loc",
            "read1",
            "--umi-len",
            "2",
            "--umi-prefix",
            "UMI",
            "--threads",
            "4",
        ],
    );
    assert_eq!(
        read_gzip_text(&out.join("trimmed.fastq.gz")),
        fq("r:UMI_AC/1 extra", "GATTACA", None)
    );
    let out = se(
        &test,
        "encoding",
        &fq("r", "ACGT", Some("hhhh")),
        &["--phred64", "--fastq-stats"],
    );
    assert_eq!(
        read_gzip_text(&out.join("trimmed.fastq.gz")),
        fq("r", "ACGT", None)
    );
}
#[test]
fn dedup_is_exact_deterministic_and_applied_after_assignment() {
    let test = TestDir::new("fastp-dedup");
    let text = (0..2500)
        .map(|i| {
            fq(
                &format!("r{i}"),
                if i % 2 == 0 { "ACGT" } else { "GATTACA" },
                None,
            )
        })
        .collect::<String>();
    for t in [1, 4, 8] {
        let out = se(
            &test,
            &format!("dedup-{t}"),
            &text,
            &["--dedup", "--threads", &t.to_string(), "--fastq-stats"],
        );
        assert_eq!(
            read_gzip_text(&out.join("trimmed.fastq.gz")),
            fq("r0", "ACGT", None) + &fq("r1", "GATTACA", None)
        );
        let metrics = std::fs::read_to_string(out.join("preprocessing.tsv")).unwrap();
        assert!(metrics.contains("duplicate_fragments\t2498\n"));
        assert!(metrics.contains("fragments_rejected\t2498\n"));
    }
}
#[test]
fn interleaved_stdin_preserves_pair_order_and_processing_limit() {
    let test = TestDir::new("fastp-stdin");
    let out = test.child("out");
    let text = (0..4)
        .map(|i| fq(&format!("r{i}/1"), "GATTACA", None) + &fq(&format!("r{i}/2"), "TGTAATC", None))
        .collect::<String>();
    let mut child = Command::new(env!("CARGO_BIN_EXE_plexless"))
        .args([
            "trim",
            "--stdin",
            "--interleaved-in",
            "--reads-to-process",
            "2",
            "--threads",
            "4",
            "--output",
        ])
        .arg(&out)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(text.as_bytes())
        .unwrap();
    let result = child.wait_with_output().unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_eq!(
        read_gzip_text(&out.join("trimmed_R1.fastq.gz")),
        fq("r0/1", "GATTACA", None) + &fq("r1/1", "GATTACA", None)
    );
}
#[test]
fn merging_uses_unique_biological_overlap_and_parallel_output() {
    let test = TestDir::new("fastp-merge");
    let insert = dna(220);
    let a = test.write("a.fq", &fq("r/1", &insert[..150], None));
    let b = test.write("b.fq", &fq("r/2", &rc(&insert[70..]), None));
    for t in [1, 4, 8] {
        let out = test.child(&format!("merge-{t}"));
        let result = Command::new(env!("CARGO_BIN_EXE_plexless"))
            .arg("trim")
            .arg("--r1")
            .arg(&a)
            .arg("--r2")
            .arg(&b)
            .arg("--output")
            .arg(&out)
            .args(["--merge", "--threads", &t.to_string(), "--fastq-stats"])
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        assert_eq!(
            read_gzip_text(&out.join("trimmed_merged.fastq.gz")),
            fq("r merged", &insert, None)
        );
    }
}

fn pe(
    test: &TestDir,
    label: &str,
    first: &str,
    second: &str,
    flags: &[&str],
) -> std::path::PathBuf {
    let out = test.child(label);
    let result = Command::new(env!("CARGO_BIN_EXE_plexless"))
        .args(["trim", "--r1"])
        .arg(test.write(&format!("{label}.r1.fq"), first))
        .arg("--r2")
        .arg(test.write(&format!("{label}.r2.fq"), second))
        .arg("--output")
        .arg(&out)
        .args(flags)
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
fn default_overlap_trims_without_an_adapter_and_gap_fallback_is_optional() {
    let test = TestDir::new("overlap-only-gap");
    let insert = dna(100);
    let a = format!("{insert}ACCTTGAGTCACGATAGGCT");
    let b = format!("{}GTCAAGCTGTACCTAGGTCA", rc(&insert));
    let out = pe(
        &test,
        "overlap",
        &fq("r/1", &a, None),
        &fq("r/2", &b, None),
        &[],
    );
    assert_eq!(
        read_gzip_text(&out.join("trimmed_R1.fastq.gz")),
        fq("r/1", &insert, None)
    );
    let out = pe(
        &test,
        "strict",
        &fq("r/1", &a, None),
        &fq("r/2", &b, None),
        &["--require-adapter-for-overlap"],
    );
    assert_eq!(
        read_gzip_text(&out.join("trimmed_R1.fastq.gz")),
        fq("r/1", &a, None)
    );
    let mut gap_insert = insert.clone();
    gap_insert.insert(43, 'A');
    let a = format!("{gap_insert}ACCTTGAGTCACGATAGGCT");
    let flags = ["--allow-gap-overlap-trimming", "--overlap-diff-limit", "0"];
    let out = pe(
        &test,
        "gap",
        &fq("r/1", &a, None),
        &fq("r/2", &b, None),
        &flags,
    );
    assert_eq!(
        read_gzip_text(&out.join("trimmed_R1.fastq.gz")),
        fq("r/1", &gap_insert, None)
    );
    assert_eq!(
        read_gzip_text(&out.join("trimmed_R2.fastq.gz")),
        fq("r/2", &rc(&insert), None)
    );
    let out = pe(
        &test,
        "gap-disabled",
        &fq("r/1", &a, None),
        &fq("r/2", &b, None),
        &["--overlap-diff-limit", "0"],
    );
    assert_eq!(
        read_gzip_text(&out.join("trimmed_R1.fastq.gz")),
        fq("r/1", &a, None)
    );
}
#[test]
fn overlap_correction_requires_asymmetric_quality_support() {
    let test = TestDir::new("overlap-correction");
    let insert = dna(100);
    let mut second = rc(&insert).into_bytes();
    second[17] = if second[17] == b'A' { b'C' } else { b'A' };
    let bad = String::from_utf8(second).unwrap();
    let mut quality = vec![b'I'; 100];
    quality[17] = b'!';
    let quality = String::from_utf8(quality).unwrap();
    let out = pe(
        &test,
        "corrected",
        &fq("r/1", &insert, None),
        &fq("r/2", &bad, Some(&quality)),
        &["--correction", "--fastq-stats"],
    );
    assert_eq!(
        read_gzip_text(&out.join("trimmed_R2.fastq.gz")),
        fq("r/2", &rc(&insert), None)
    );
    let metrics = std::fs::read_to_string(out.join("preprocessing.tsv")).unwrap();
    assert!(metrics.contains("corrected_bases\t1\n"));
    let out = pe(
        &test,
        "high-quality-disagreement",
        &fq("r/1", &insert, None),
        &fq("r/2", &bad, None),
        &["--correction"],
    );
    assert_eq!(
        read_gzip_text(&out.join("trimmed_R2.fastq.gz")),
        fq("r/2", &bad, None)
    );
}
#[test]
fn reports_are_standalone_escaped_and_do_not_change_retention() {
    let test = TestDir::new("fastp-reports");
    let json = test.child("report.json");
    let html = test.child("report.html");
    let sequence = dna(70);
    let text = (0..25)
        .map(|i| fq(&format!("r{i}"), &sequence, None))
        .collect::<String>();
    let out = se(
        &test,
        "reports",
        &text,
        &[
            "--json",
            json.to_str().unwrap(),
            "--html",
            html.to_str().unwrap(),
            "--report-title",
            "</script><b>title</b>",
            "--overrepresentation-analysis",
            "--overrepresentation-sampling",
            "1",
            "--evaluate-duplication",
            "--duplication-memory-mib",
            "1",
            "--threads",
            "4",
        ],
    );
    assert_eq!(read_gzip_text(&out.join("trimmed.fastq.gz")), text);
    let report: serde_json::Value = serde_json::from_slice(&std::fs::read(json).unwrap()).unwrap();
    assert_eq!(report["processing"]["fragments_retained"], 25);
    assert_eq!(report["duplication"]["duplicate_fragments"], 24);
    assert_eq!(
        report["read1_before_filtering"]["summary"]["total_bases"],
        1750
    );
    assert_eq!(
        report["read1_after_filtering"]["summary"]["total_bases"],
        1750
    );
    assert!(
        !report["overrepresented_sequences"]["R1_after"]["sequences"]
            .as_array()
            .unwrap()
            .is_empty()
    );
    let html = std::fs::read_to_string(html).unwrap();
    assert!(html.contains("&lt;/script&gt;&lt;b&gt;title&lt;/b&gt;"));
    assert!(!html.contains("</script><b>title"));
    assert!(!html.contains("src=\"http"));
}
#[test]
fn split_output_preserves_records_pairing_and_valid_gzip() {
    let test = TestDir::new("fastp-split");
    let text = (0..601)
        .map(|i| fq(&format!("r{i}"), "GATTACA", None))
        .collect::<String>();
    for flags in [vec!["--split-by-lines", "1000"], vec!["--split", "3"]] {
        let label = if flags[0] == "--split" {
            "count"
        } else {
            "lines"
        };
        let out = se(&test, label, &text, &flags);
        let mut combined = String::new();
        let mut paths: Vec<_> = std::fs::read_dir(&out)
            .unwrap()
            .map(|e| e.unwrap().path())
            .filter(|p| p.extension().is_some_and(|e| e == "gz"))
            .collect();
        paths.sort();
        assert_eq!(paths.len(), 3);
        for path in paths {
            let part = read_gzip_text(&path);
            assert_eq!(part.lines().count() % 4, 0);
            combined.push_str(&part);
        }
        assert_eq!(combined, text);
    }
    let first = (0..520)
        .map(|i| fq(&format!("r{i}/1"), "GATTACA", None))
        .collect::<String>();
    let second = (0..520)
        .map(|i| fq(&format!("r{i}/2"), "TGCATGT", None))
        .collect::<String>();
    let out = pe(
        &test,
        "paired",
        &first,
        &second,
        &["--split-by-lines", "1000", "--threads", "4"],
    );
    for part in 1..=3 {
        let a = read_gzip_text(&out.join(format!("{part:04}.trimmed_R1.fastq.gz")));
        let b = read_gzip_text(&out.join(format!("{part:04}.trimmed_R2.fastq.gz")));
        assert_eq!(a.lines().count(), b.lines().count());
    }
}
#[test]
fn index_filtering_index_umi_and_merged_unmerged_outputs() {
    let test = TestDir::new("fastp-index");
    let blacklist = test.write("blacklist.txt", "ACGT\n");
    let out = se(
        &test,
        "index",
        &(fq("r1 1:N:0:ACGT+TGCA", "GATTACA", None) + &fq("r2 1:N:0:TTTT+AAAA", "GATTACA", None)),
        &["--filter-by-index1", blacklist.to_str().unwrap()],
    );
    assert_eq!(
        read_gzip_text(&out.join("trimmed.fastq.gz")),
        fq("r2 1:N:0:TTTT+AAAA", "GATTACA", None)
    );
    let out = pe(
        &test,
        "index-umi",
        &fq("r/1 1:N:0:ACGT+TGCA", "GATTACA", None),
        &fq("r/2 2:N:0:ACGT+TGCA", "TGCATGT", None),
        &["--umi", "--umi-loc", "per-index"],
    );
    assert!(read_gzip_text(&out.join("trimmed_R1.fastq.gz")).starts_with("@r:ACGT_TGCA/1 "));
    let merged = test.child("merged.fq");
    let out = pe(
        &test,
        "unmerged",
        &fq("r/1", "GATTACA", None),
        &fq("r/2", "TGCATGT", None),
        &[
            "--merge",
            "--include-unmerged",
            "--merged-out",
            merged.to_str().unwrap(),
        ],
    );
    assert_eq!(
        std::fs::read_to_string(merged).unwrap(),
        fq("r/1", "GATTACA", None) + &fq("r/2", "TGCATGT", None)
    );
    assert!(read_gzip_text(&out.join("trimmed_R1.fastq.gz")).is_empty());
}
#[test]
fn failed_and_shared_unpaired_outputs_keep_input_order() {
    let test = TestDir::new("fastp-unpaired");
    let passed = test.child("unpaired.fq");
    let failed = test.child("failed.fq");
    let first = fq("a/1", "GATTACA", Some("!!!!!!!")) + &fq("b/1", "GATTACA", None);
    let second = fq("a/2", "TGCATGT", None) + &fq("b/2", "TGCATGT", Some("!!!!!!!"));
    let out = pe(
        &test,
        "filter",
        &first,
        &second,
        &[
            "--min-average-quality",
            "20",
            "--unpaired1",
            passed.to_str().unwrap(),
            "--unpaired2",
            passed.to_str().unwrap(),
            "--failed-out",
            failed.to_str().unwrap(),
            "--threads",
            "4",
        ],
    );
    assert_eq!(
        std::fs::read_to_string(passed).unwrap(),
        fq("a/2", "TGCATGT", None) + &fq("b/1", "GATTACA", None)
    );
    assert_eq!(
        std::fs::read_to_string(failed).unwrap(),
        fq("a/1", "GATTACA", Some("!!!!!!!")) + &fq("b/2", "TGCATGT", Some("!!!!!!!"))
    );
    assert!(read_gzip_text(&out.join("trimmed_R1.fastq.gz")).is_empty());
    assert!(read_gzip_text(&out.join("trimmed_R2.fastq.gz")).is_empty());
}
#[test]
fn adapter_fasta_profile_and_stdout_are_effective() {
    let test = TestDir::new("fastp-profile");
    let adapter = "ACGTTGCACTGATCGAGGTACCTATGCCAGTA";
    let fasta = test.write("adapter.fa", &format!(">custom\n{adapter}\n"));
    let out = se(
        &test,
        "fasta",
        &fq("r", &format!("GATTACA{adapter}"), None),
        &["--adapter-fasta", fasta.to_str().unwrap()],
    );
    assert_eq!(
        read_gzip_text(&out.join("trimmed.fastq.gz")),
        fq("r", "GATTACA", None)
    );
    let out = se(
        &test,
        "profile",
        &fq("NS500:1:1:1", &format!("{}CCGGGGGGGGGG", dna(48)), None),
        &["--fastp-defaults", "--duplication-memory-mib", "1"],
    );
    assert!(out.join("plexless.json").exists());
    assert!(out.join("plexless.html").exists());
    let metrics = std::fs::read_to_string(out.join("preprocessing.tsv")).unwrap();
    assert!(metrics.contains("R1_poly_g_bases_removed\t10\n"));
    let out = test.child("stream");
    let input = test.write("stream.fq", &fq("r", "GATTACA", None));
    let result = Command::new(env!("CARGO_BIN_EXE_plexless"))
        .args(["trim", "--reads"])
        .arg(input)
        .arg("--output")
        .arg(out)
        .arg("--stdout")
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_eq!(
        String::from_utf8(result.stdout).unwrap(),
        fq("r", "GATTACA", None)
    );
}

#[test]
fn fastp_1_4_0_sequence_and_quality_goldens() {
    let test = TestDir::new("fastp-goldens");
    let golden: serde_json::Value =
        serde_json::from_str(include_str!("fixtures/fastp_v1_4_0.json")).unwrap();
    for case in golden["cases"].as_array().unwrap() {
        let name = case["name"].as_str().unwrap();
        let first = test.write(&format!("{name}.r1.fq"), case["r1"].as_str().unwrap());
        let out = test.child(name);
        let output1 = test.child(&format!("{name}.out1.fq"));
        let output2 = test.child(&format!("{name}.out2.fq"));
        let merged = test.child(&format!("{name}.merged.fq"));
        let mut cmd = Command::new(env!("CARGO_BIN_EXE_plexless"));
        cmd.arg("trim")
            .arg("--output")
            .arg(out)
            .arg("--out1")
            .arg(&output1);
        if let Some(text) = case["r2"].as_str() {
            cmd.arg("--r1")
                .arg(first)
                .arg("--r2")
                .arg(test.write(&format!("{name}.r2.fq"), text))
                .arg("--out2")
                .arg(&output2);
        } else {
            cmd.arg("--reads").arg(first);
        }
        for flag in case["args"].as_array().unwrap() {
            cmd.arg(flag.as_str().unwrap());
        }
        if name == "merge" {
            cmd.arg("--merged-out").arg(&merged);
        }
        let run = cmd.output().unwrap();
        assert!(
            run.status.success(),
            "{name}: {}",
            String::from_utf8_lossy(&run.stderr)
        );
        let paths = if name == "merge" {
            vec![merged]
        } else if case["r2"].is_string() {
            vec![output1, output2]
        } else {
            vec![output1]
        };
        for (path, expected) in paths.iter().zip(case["outputs"].as_array().unwrap()) {
            let text = std::fs::read_to_string(path).unwrap();
            let lines: Vec<_> = text.lines().collect();
            let records: Vec<_> = lines
                .as_chunks::<4>()
                .0
                .iter()
                .map(|r| serde_json::json!([r[1], r[3]]))
                .collect();
            assert_eq!(serde_json::json!(records), *expected, "{name}");
        }
    }
}
#[test]
fn invalid_feature_combinations_fail_before_writing_fastq() {
    let test = TestDir::new("fastp-invalid");
    let input = test.write("in.fq", &fq("read", "GATTACA", None));
    let collision = test.child("same.fq");
    for (name, flags) in [
        (
            "collision",
            vec![
                "--out1",
                collision.to_str().unwrap(),
                "--json",
                collision.to_str().unwrap(),
            ],
        ),
        ("gap", vec!["--allow-gap-overlap-trimming"]),
        ("split", vec!["--split-by-lines", "1001"]),
        (
            "umi",
            vec!["--umi", "--umi-loc", "read1", "--umi-len", "100"],
        ),
        ("index", vec!["--umi", "--umi-loc", "index1"]),
    ] {
        let out = test.child(name);
        let result = Command::new(env!("CARGO_BIN_EXE_plexless"))
            .args(["trim", "--reads"])
            .arg(&input)
            .arg("--output")
            .arg(&out)
            .args(flags)
            .output()
            .unwrap();
        assert!(
            !result.status.success(),
            "{name} should reject invalid config/input"
        );
        assert!(
            !out.join("trimmed.fastq.gz").exists(),
            "{name} emitted invalid output"
        );
    }
}

#[test]
fn output_aliases_default_names_and_reserved_reports_cannot_collide() {
    let test = TestDir::new("fastp-output-collision");
    let input = test.write("in.fq", &fq("r", &dna(70), None));
    for label in ["default-mate", "profile-report", "reserved"] {
        let out = test.child(label);
        let destination = out.join(match label {
            "default-mate" => "trimmed_R1.fastq.gz",
            "profile-report" => "plexless.json",
            _ => "preprocessing.tsv",
        });
        let mut cmd = Command::new(env!("CARGO_BIN_EXE_plexless"));
        cmd.arg("trim").arg("--output").arg(&out);
        if label == "default-mate" {
            cmd.arg("--r1")
                .arg(&input)
                .arg("--r2")
                .arg(&input)
                .arg("--out2")
                .arg(&destination);
        } else {
            cmd.arg("--reads")
                .arg(&input)
                .arg("--out1")
                .arg(&destination);
            if label == "profile-report" {
                cmd.arg("--fastp-defaults");
            }
        }
        let result = cmd.output().unwrap();
        assert!(
            !result.status.success(),
            "{label} should reject a collision"
        );
        assert!(!destination.exists(), "{label} wrote a colliding output");
    }
}
#[test]
fn splitting_never_appends_to_an_existing_later_part() {
    let test = TestDir::new("fastp-split-existing");
    let destination = test.child("custom.fq");
    let old = test.write("0002.custom.fq", "previous content\n");
    let text = (0..300)
        .map(|i| fq(&format!("r{i}"), "GATTACA", None))
        .collect::<String>();
    let input = test.write("in.fq", &text);
    let out = test.child("out");
    let result = Command::new(env!("CARGO_BIN_EXE_plexless"))
        .args(["trim", "--reads"])
        .arg(input)
        .arg("--output")
        .arg(&out)
        .arg("--out1")
        .arg(destination)
        .args(["--split-by-lines", "1000"])
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert_eq!(std::fs::read_to_string(old).unwrap(), "previous content\n");
    assert!(out.join("PLEXLESS_INCOMPLETE").exists());
}

#[test]
fn demux_assignment_precedes_umi_dedup_adapter_and_filtering() {
    let test = TestDir::new("demux-combined-preprocessing");
    let adapter = "ACGTTGCACTGATCGAGGTACCTATGCCAGTATCGGATCAACGTGTCACAGTCTAGGACTTCAGAT";
    let insert = dna(80);
    let sequence = format!("ACGTGG{insert}{adapter}");
    let input = test.write(
        "reads.fq",
        &(fq("first/1", &sequence, None)
            + &fq("duplicate/1", &sequence, None)
            + &fq(
                "short/1",
                &format!("ACGTTT{}{adapter}", &insert[..10]),
                None,
            )
            + &fq(
                "partial/1",
                &format!("ACGTGT{insert}{}", &adapter[..30]),
                None,
            )
            + &fq("unmatched/1", &format!("TTTTGG{insert}{adapter}"), None)),
    );
    let barcodes = test.write("barcodes.tsv", "Set\tID\tSequence\nA\ta\tACGT\n");
    let samples = test.write("samples.tsv", "Sample\tA\nsample\ta\n");
    for threads in [1, 4, 8] {
        let out = test.child(&format!("out-{threads}"));
        let result = Command::new(env!("CARGO_BIN_EXE_plexless"))
            .args(["demux", "--reads"])
            .arg(&input)
            .args(["--structure", "R1_4A", "--barcodes"])
            .arg(&barcodes)
            .arg("--samples")
            .arg(&samples)
            .arg("--output")
            .arg(&out)
            .arg("--threads")
            .arg(threads.to_string())
            .args([
                "--umi",
                "--umi-loc",
                "read1",
                "--umi-len",
                "2",
                "--dedup",
                "--duplication-memory-mib",
                "1",
                "--r1-distal-sequence",
                adapter,
                "--partial-pass-through",
                "--min-length",
                "20",
                "--write-unassigned",
                "--fastq-stats",
            ])
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        assert_eq!(
            read_gzip_text(&out.join("sample.fastq.gz")),
            fq("first:GG/1", &insert, None) + &fq("partial:GT/1", &insert, None)
        );
        let metrics = std::fs::read_to_string(out.join("preprocessing.tsv")).unwrap();
        for (key, count) in [
            ("fragments_processed", 5),
            ("fragments_retained", 3),
            ("fragments_rejected", 2),
            ("duplicate_fragments", 1),
            ("rejected_min_length", 1),
            ("R1_full_adapter_reads", 3),
            ("R1_partial_adapter_reads", 1),
        ] {
            assert!(metrics.contains(&format!("{key}\t{count}\n")), "{metrics}");
        }
    }
}
#[test]
fn demux_merging_keeps_per_sample_fanout_and_reconciles_qc() {
    let test = TestDir::new("demux-merge-fanout");
    let insert = dna(220);
    let r1 = test.write(
        "r1.fq",
        &(fq("a/1", &format!("ACGT{}", &insert[..150]), None)
            + &fq("b/1", &format!("GGGG{}", &insert[..150]), None)),
    );
    let r2 = test.write(
        "r2.fq",
        &(fq("a/2", &rc(&insert[70..]), None) + &fq("b/2", &rc(&insert[70..]), None)),
    );
    let barcodes = test.write(
        "barcodes.tsv",
        "Set\tID\tSequence\nA\ta\tACGT\nA\tb\tGGGG\n",
    );
    let samples = test.write("samples.tsv", "Sample\tA\nfirst\ta\nsecond\tb\n");
    for threads in [1, 4] {
        let out = test.child(&format!("out-{threads}"));
        let result = Command::new(env!("CARGO_BIN_EXE_plexless"))
            .args(["demux", "--r1"])
            .arg(&r1)
            .arg("--r2")
            .arg(&r2)
            .args(["--r1-structure", "R1_4A", "--barcodes"])
            .arg(&barcodes)
            .arg("--samples")
            .arg(&samples)
            .arg("--output")
            .arg(&out)
            .arg("--threads")
            .arg(threads.to_string())
            .args(["--merge", "--fastq-stats"])
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        assert_eq!(
            read_gzip_text(&out.join("first_merged.fastq.gz")),
            fq("a merged", &insert, None)
        );
        assert_eq!(
            read_gzip_text(&out.join("second_merged.fastq.gz")),
            fq("b merged", &insert, None)
        );
        let metrics = std::fs::read_to_string(out.join("preprocessing.tsv")).unwrap();
        assert!(metrics.contains("merged_fragments\t2\n"));
        assert!(metrics.contains("fragments_retained\t2\n"));
    }
}

#[test]
fn merging_after_asymmetric_front_clipping_preserves_mate_supported_biological_overhangs() {
    let test = TestDir::new("merge-front-overhangs");
    let insert = dna(100);
    let a = format!("{insert}ACCTTGAGTCACGATAGGCT");
    let b = format!("{}GTCAAGCTGTACCTAGGTCA", rc(&insert));
    let merged = test.child("merged.fq");
    let out = pe(
        &test,
        "merged",
        &fq("r/1", &a, None),
        &fq("r/2", &b, None),
        &[
            "--merge",
            "--trim-front1",
            "3",
            "--trim-front2",
            "7",
            "--merged-out",
            merged.to_str().unwrap(),
            "--fastq-stats",
        ],
    );
    assert_eq!(
        std::fs::read_to_string(merged).unwrap(),
        fq("r merged", &insert, None)
    );
    let qc = std::fs::read_to_string(out.join("retained_fastq_stats.tsv")).unwrap();
    assert!(qc.contains("Merged\t1\t100\t"), "{qc}");
}

#[test]
fn inferred_distal_signatures_include_opposite_read_umi_cycles_without_counting_them_as_evidence() {
    let test = TestDir::new("demux-umi-inferred-distal");
    let t1 = "ACGTTGCACTGATCGAGGTACCTA";
    let t2 = "GATCCGTAGTCAACGCTTGACGTA";
    let insert = dna(100);
    let first = format!("ACGT{t1}GG{insert}GG{}", rc(t2));
    let second = format!("{t2}CC{}CC{}ACGT", rc(&insert), rc(t1));
    let r1 = test.write("r1.fq", &fq("r/1", &first, None));
    let r2 = test.write("r2.fq", &fq("r/2", &second, None));
    let barcodes = test.write("barcodes.tsv", "Set\tID\tSequence\nA\ta\tACGT\n");
    let samples = test.write("samples.tsv", "Sample\tA\nsample\ta\n");
    let technical = test.write(
        "technical.tsv",
        &format!("Mate\tSequence\nR1\t{t1}\nR2\t{t2}\n"),
    );
    for threads in [1, 4] {
        let out = test.child(&format!("out-{threads}"));
        let result = Command::new(env!("CARGO_BIN_EXE_plexless"))
            .args(["demux", "--r1"])
            .arg(&r1)
            .arg("--r2")
            .arg(&r2)
            .args([
                "--r1-structure",
                "R1_4A24T",
                "--r2-structure",
                "R2_24T",
                "--barcodes",
            ])
            .arg(&barcodes)
            .arg("--samples")
            .arg(&samples)
            .arg("--technical-sequences")
            .arg(&technical)
            .arg("--output")
            .arg(&out)
            .arg("--threads")
            .arg(threads.to_string())
            .args(["--umi", "--umi-loc", "per-read", "--umi-len", "2"])
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        assert_eq!(
            read_gzip_text(&out.join("sample_R1.fastq.gz")),
            fq("r:GG_CC/1", &insert, None)
        );
        assert_eq!(
            read_gzip_text(&out.join("sample_R2.fastq.gz")),
            fq("r:GG_CC/2", &rc(&insert), None)
        );
        let metrics = std::fs::read_to_string(out.join("preprocessing.tsv")).unwrap();
        assert!(metrics.contains("R1_full_adapter_reads\t1\n"));
        assert!(metrics.contains("R2_full_adapter_reads\t1\n"));
    }
}

#[test]
fn paired_phred64_normalizes_both_parser_paths_and_preserves_pair_quality() {
    let test = TestDir::new("paired-phred64");
    let insert = dna(80);
    let q = "h".repeat(80);
    for threads in [1, 4, 8] {
        let out = pe(
            &test,
            &format!("out-{threads}"),
            &fq("r/1", &insert, Some(&q)),
            &fq("r/2", &rc(&insert), Some(&q)),
            &[
                "--phred64",
                "--threads",
                &threads.to_string(),
                "--fastq-stats",
            ],
        );
        assert_eq!(
            read_gzip_text(&out.join("trimmed_R1.fastq.gz")),
            fq("r/1", &insert, None)
        );
        assert_eq!(
            read_gzip_text(&out.join("trimmed_R2.fastq.gz")),
            fq("r/2", &rc(&insert), None)
        );
    }
}
