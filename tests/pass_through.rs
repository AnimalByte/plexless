mod common;

use common::{TestDir, read_gzip_text};
use std::process::Command;

const A: &str = "ACGTTGCACTGATCGAGGTACCTA";
const B: &str = "GCTACAGTTCGATACCGTGA";

fn rc(sequence: &str) -> String {
    sequence
        .bytes()
        .rev()
        .map(|base| match base {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            _ => 'N',
        })
        .collect()
}

fn fastq(name: &str, sequence: &str) -> String {
    let quality: String = (0..sequence.len())
        .map(|i| (33 + (i % 60) as u8) as char)
        .collect();
    format!("@{name}\n{sequence}\n+\n{quality}\n")
}

fn expected_record(name: &str, raw: &str, start: usize, end: usize) -> String {
    let quality: String = (start..end)
        .map(|i| (33 + (i % 60) as u8) as char)
        .collect();
    format!("@{name}\n{}\n+\n{quality}\n", &raw[start..end])
}

#[test]
fn se_default_symmetric_inference_and_disable_across_writers() {
    let test = TestDir::new("pass-through-se");
    let barcodes = test.write("barcodes.tsv", &format!("Set\tID\tSequence\nA\ta\t{A}\n"));
    let samples = test.write("samples.tsv", "Sample\tA\ns\ta\n");
    // T is unknown and reverses to the beginning of the distal structure.
    let proximal = format!("{A}CGN");
    let distal = format!("TNN{}", rc(A));
    let negative = format!("GGGG{}GG", &rc(A)[..12]);
    let inputs = [
        format!("{proximal}GATTACA{distal}AAAA"),
        format!("{proximal}{distal}"),
        format!("{proximal}G{distal}"),
        format!("{proximal}{negative}"),
        format!("{proximal}{}", rc(B)),
    ];
    let ends = [
        proximal.len() + 7,
        proximal.len(),
        proximal.len() + 1,
        inputs[3].len(),
        inputs[4].len(),
    ];
    let reads = test.write(
        "reads.fastq",
        &inputs
            .iter()
            .enumerate()
            .map(|(i, seq)| fastq(&format!("read{i}"), seq))
            .collect::<String>(),
    );
    for threads in [1, 4] {
        for mode in ["direct", "buffered"] {
            for disabled in [false, true] {
                let out = test.child(&format!("out-{threads}-{mode}-{disabled}"));
                let mut command = Command::new(env!("CARGO_BIN_EXE_plexless"));
                command
                    .args([
                        "--threads",
                        &threads.to_string(),
                        "demux",
                        "--structure",
                        "R1_24A3T",
                        "--output-mode",
                        mode,
                    ])
                    .arg("--reads")
                    .arg(&reads)
                    .arg("--barcodes")
                    .arg(&barcodes)
                    .arg("--samples")
                    .arg(&samples)
                    .arg("--output")
                    .arg(&out);
                if disabled {
                    command.arg("--no-pass-through");
                }
                let result = command.output().unwrap();
                assert!(
                    result.status.success(),
                    "{}",
                    String::from_utf8_lossy(&result.stderr)
                );
                let expected: String = inputs
                    .iter()
                    .enumerate()
                    .map(|(i, seq)| {
                        expected_record(
                            &format!("read{i}"),
                            seq,
                            proximal.len(),
                            if disabled { seq.len() } else { ends[i] },
                        )
                    })
                    .collect();
                assert_eq!(read_gzip_text(&out.join("s.fastq.gz")), expected);
                let summary = String::from_utf8_lossy(&result.stderr);
                assert!(summary.contains(if disabled {
                    "pass-through R1: 0"
                } else {
                    "pass-through R1: 3"
                }));
                if disabled {
                    assert!(summary.contains("no pass-through (scanned mates): 0"));
                }
            }
        }
    }
}

#[test]
fn pe_automatic_inference_repeated_pieces_rc_and_independent_mates() {
    let test = TestDir::new("pass-through-pe");
    let b = "TCGAGGTACCTAGCTACAGTTCGA";
    let barcodes = test.write(
        "barcodes.tsv",
        &format!("Set\tID\tSequence\nA\ta\t{A}\nB\tb\t{b}\n"),
    );
    let samples = test.write("samples.tsv", "Sample\tA\tB\ns\ta\tb\n");
    let p1 = format!("{}NN{}", &A[..12], rc(&b[..12]));
    let p2 = format!("{}{}", rc(&A[12..]), &b[12..]);
    let d1 = rc(&p2);
    let d2 = rc(&p1);
    let mut r1 = String::new();
    let mut r2 = String::new();
    let mut expected1 = String::new();
    let mut expected2 = String::new();
    for index in 0..4 {
        let bio1 = "G".repeat(90);
        let bio2 = "C".repeat(173);
        let seq1 = format!("{p1}{bio1}{}", if index & 1 != 0 { &d1 } else { "" });
        let seq2 = format!("{p2}{bio2}{}", if index & 2 != 0 { &d2 } else { "" });
        r1.push_str(&fastq(&format!("pair{index}/1"), &seq1));
        r2.push_str(&fastq(&format!("pair{index}/2"), &seq2));
        expected1.push_str(&expected_record(
            &format!("pair{index}/1"),
            &seq1,
            p1.len(),
            p1.len() + bio1.len(),
        ));
        expected2.push_str(&expected_record(
            &format!("pair{index}/2"),
            &seq2,
            p2.len(),
            p2.len() + bio2.len(),
        ));
    }
    let r1 = test.write("r1.fastq", &r1);
    let r2 = test.write("r2.fastq", &r2);
    for threads in [1, 4] {
        for mode in ["direct", "buffered"] {
            let out = test.child(&format!("out-{threads}-{mode}"));
            let result = Command::new(env!("CARGO_BIN_EXE_plexless"))
                .args([
                    "--threads",
                    &threads.to_string(),
                    "demux",
                    "--r1-structure",
                    "R1_12A2T12B(rc)",
                    "--r2-structure",
                    "R2_12A(rc)12B",
                    "--output-mode",
                    mode,
                ])
                .arg("--r1")
                .arg(&r1)
                .arg("--r2")
                .arg(&r2)
                .arg("--barcodes")
                .arg(&barcodes)
                .arg("--samples")
                .arg(&samples)
                .arg("--output")
                .arg(&out)
                .output()
                .unwrap();
            assert!(
                result.status.success(),
                "{}",
                String::from_utf8_lossy(&result.stderr)
            );
            assert_eq!(read_gzip_text(&out.join("s_R1.fastq.gz")), expected1);
            assert_eq!(read_gzip_text(&out.join("s_R2.fastq.gz")), expected2);
        }
    }
}

#[test]
fn explicit_override_long_structure_orientation_and_sample_specificity() {
    let test = TestDir::new("pass-through-override");
    let other = "TTGACCGATAGCTACGTCATGGAC";
    let barcodes = test.write(
        "barcodes.tsv",
        &format!("Set\tID\tSequence\nA\ta\t{A}\nA\tother\t{other}\nB\tb\t{B}\n"),
    );
    let samples = test.write("samples.tsv", "Sample\tA\tB\ns\ta\tb\nother\tother\tb\n");
    let proximal = format!("{}{}{}", &A[..12], rc(&A[12..]), B);
    let automatic = rc(&proximal);
    // Override is expressed in observed orientation; repeated A consumes canonical pieces.
    let overridden = format!("{B}NN{}{}", rc(&A[..12]), &A[12..]);
    let wrong_sample = format!("{}{}{}", rc(B), &other[12..], rc(&other[..12]));
    let inputs = [
        format!("{proximal}GATTACA{automatic}"),
        format!("{proximal}GG{overridden}"),
        format!("{proximal}CC{wrong_sample}"),
    ];
    let reads = test.write(
        "reads.fastq",
        &inputs
            .iter()
            .enumerate()
            .map(|(i, seq)| fastq(&format!("r{i}"), seq))
            .collect::<String>(),
    );
    for override_enabled in [false, true] {
        let out = test.child(&format!("out-{override_enabled}"));
        let mut command = Command::new(env!("CARGO_BIN_EXE_plexless"));
        command
            .args([
                "demux",
                "--structure",
                "R1_12A12A(rc)20B",
                "--max-mismatches",
                "0",
            ])
            .arg("--reads")
            .arg(&reads)
            .arg("--barcodes")
            .arg(&barcodes)
            .arg("--samples")
            .arg(&samples)
            .arg("--output")
            .arg(&out);
        if override_enabled {
            command.args(["--r1-distal-structure", "R1_20B2T12A(rc)12A"]);
        }
        let result = command.output().unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        let expected: String = inputs
            .iter()
            .enumerate()
            .map(|(i, seq)| {
                let end = match (override_enabled, i) {
                    (false, 0) => proximal.len() + 7,
                    (true, 1) => proximal.len() + 2,
                    _ => seq.len(),
                };
                expected_record(&format!("r{i}"), seq, proximal.len(), end)
            })
            .collect();
        assert_eq!(read_gzip_text(&out.join("s.fastq.gz")), expected);
    }
}

#[test]
fn unavailable_anchor_diagnostic_and_override_validation() {
    let test = TestDir::new("pass-through-config");
    let reads = test.write("reads.fastq", &fastq("r", "ACGTGATTACAACGT"));
    let barcodes = test.write("barcodes.tsv", "Set\tID\tSequence\nA\ta\tACGT\n");
    let samples = test.write("samples.tsv", "Sample\tA\ns\ta\n");
    for (index, extra) in [
        vec![],
        vec!["--r1-distal-structure", "R1_4B"],
        vec!["--r1-distal-structure", "R2_4A"],
        vec!["--r2-distal-structure", "R2_4A"],
        vec!["--r1-distal-structure", "R1_3A"],
        vec!["--r1-distal-structure", "R1_5A"],
    ]
    .iter()
    .enumerate()
    {
        let out = test.child(&format!("out-{index}"));
        let result = Command::new(env!("CARGO_BIN_EXE_plexless"))
            .args(["demux", "--structure", "R1_4A"])
            .args(extra)
            .arg("--reads")
            .arg(&reads)
            .arg("--barcodes")
            .arg(&barcodes)
            .arg("--samples")
            .arg(&samples)
            .arg("--output")
            .arg(&out)
            .output()
            .unwrap();
        assert_eq!(result.status.success(), index == 0);
        if index == 0 {
            assert!(String::from_utf8_lossy(&result.stderr).contains("R1 unavailable"));
            assert!(read_gzip_text(&out.join("s.fastq.gz")).contains("GATTACAACGT"));
        } else {
            assert!(!out.exists());
        }
    }
}

#[test]
fn corrected_assignment_uses_canonical_distal_and_unassigned_is_not_scanned() {
    let test = TestDir::new("pass-through-assignment");
    let barcodes = test.write("barcodes.tsv", &format!("Set\tID\tSequence\nA\ta\t{A}\n"));
    let samples = test.write("samples.tsv", "Sample\tA\ns\ta\n");
    let corrected = format!("C{}GATTACA{}", &A[1..], rc(A));
    let unmatched = format!("{}GATTACA{}", "C".repeat(24), rc(A));
    let reads = test.write(
        "reads.fastq",
        &(fastq("corrected", &corrected) + &fastq("unmatched", &unmatched)),
    );
    for threads in [1, 4] {
        let out = test.child(&format!("out-{threads}"));
        let result = Command::new(env!("CARGO_BIN_EXE_plexless"))
            .args([
                "--threads",
                &threads.to_string(),
                "demux",
                "--structure",
                "R1_24A",
                "--write-unassigned",
            ])
            .arg("--reads")
            .arg(&reads)
            .arg("--barcodes")
            .arg(&barcodes)
            .arg("--samples")
            .arg(&samples)
            .arg("--output")
            .arg(&out)
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        assert_eq!(
            read_gzip_text(&out.join("s.fastq.gz")),
            expected_record("corrected", &corrected, 24, 31)
        );
        assert_eq!(
            read_gzip_text(&out.join("unassigned.fastq.gz")),
            fastq("unmatched", &unmatched)
        );
        assert!(String::from_utf8_lossy(&result.stderr).contains("pass-through R1: 1"));
    }
}

#[test]
fn pe_false_anchor_candidates_do_not_hide_later_distal_or_affect_other_mate() {
    let test = TestDir::new("pass-through-pe-false-anchors");
    let b = "GCTACAGTTCGATACCGTGATCAG";
    let barcodes = test.write(
        "barcodes.tsv",
        &format!("Set\tID\tSequence\nA\ta\t{A}\nB\tb\t{b}\n"),
    );
    let samples = test.write("samples.tsv", "Sample\tA\tB\ns\ta\tb\n");
    let distal1 = format!("NNN{}", rc(b));
    let distal2 = format!("NNN{}", rc(A));
    let mut bio1 = "GATTACA".to_string();
    // Move one observed N through each known position. Some variants break
    // the chosen anchor; the others keep it exact but must fail verification.
    // This exercises the behavior without depending on which anchor is chosen.
    for position in 3..distal1.len() {
        let mut fake = distal1.clone().into_bytes();
        fake[position] = b'N';
        bio1.push_str(&String::from_utf8(fake).unwrap());
        bio1.push_str("GGGGGGG");
    }
    let p1 = format!("{A}NNN");
    let p2 = format!("{b}NNN");
    let inputs1 = [format!("{p1}{bio1}"), format!("{p1}{bio1}{distal1}AAAA")];
    let inputs2 = [format!("{p2}CCTAGG{distal2}AAAA"), format!("{p2}CCTAGG")];
    let r1 = test.write(
        "r1.fastq",
        &(fastq("pair0/1", &inputs1[0]) + &fastq("pair1/1", &inputs1[1])),
    );
    let r2 = test.write(
        "r2.fastq",
        &(fastq("pair0/2", &inputs2[0]) + &fastq("pair1/2", &inputs2[1])),
    );
    for threads in [1, 4] {
        for mode in ["direct", "buffered"] {
            let out = test.child(&format!("out-{threads}-{mode}"));
            let result = Command::new(env!("CARGO_BIN_EXE_plexless"))
                .args([
                    "--threads",
                    &threads.to_string(),
                    "demux",
                    "--r1-structure",
                    "R1_24A3T",
                    "--r2-structure",
                    "R2_24B3T",
                    "--output-mode",
                    mode,
                ])
                .arg("--r1")
                .arg(&r1)
                .arg("--r2")
                .arg(&r2)
                .arg("--barcodes")
                .arg(&barcodes)
                .arg("--samples")
                .arg(&samples)
                .arg("--output")
                .arg(&out)
                .output()
                .unwrap();
            assert!(
                result.status.success(),
                "{}",
                String::from_utf8_lossy(&result.stderr)
            );
            let expected1: String = inputs1
                .iter()
                .enumerate()
                .map(|(i, raw)| {
                    expected_record(&format!("pair{i}/1"), raw, p1.len(), p1.len() + bio1.len())
                })
                .collect();
            let expected2: String = inputs2
                .iter()
                .enumerate()
                .map(|(i, raw)| expected_record(&format!("pair{i}/2"), raw, p2.len(), p2.len() + 6))
                .collect();
            assert_eq!(read_gzip_text(&out.join("s_R1.fastq.gz")), expected1);
            assert_eq!(read_gzip_text(&out.join("s_R2.fastq.gz")), expected2);
        }
    }
}
