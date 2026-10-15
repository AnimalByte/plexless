mod common;
use common::{TestDir, read_gzip_text};
use std::io::Write;
use std::path::Path;
use std::process::{Command, Output};

// Synthetic, non-manufacturer sequences used solely as test fixtures.
const ADAPTER: &str = "ACGTTGCACTGATCGAGGTACCTATGCCAGTATCGGATCAACGTGTCACAGTCTAGGACTTCAGATGCCTAGCATCGTACGATTCGACCTGATGACCTGACGTTGCACTGATCGAGGTACCTATGCCAGTA";
const BARCODE: &str = "GCTACAGTTCGATACCGTGA";
fn fq(id: &str, seq: &str, quality: Option<&str>) -> String {
    format!(
        "@{id}\n{seq}\n+\n{}\n",
        quality
            .map(str::to_owned)
            .unwrap_or_else(|| "I".repeat(seq.len()))
    )
}
fn rc(seq: &str) -> String {
    seq.bytes()
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
fn success(output: Output) {
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}
fn run_trim(reads: &Path, out: &Path, threads: usize, extra: &[&str]) {
    success(
        Command::new(env!("CARGO_BIN_EXE_plexless"))
            .args(["trim", "--reads"])
            .arg(reads)
            .arg("--output")
            .arg(out)
            .arg("--threads")
            .arg(threads.to_string())
            .args(extra)
            .output()
            .unwrap(),
    );
}
fn sequences(path: &Path) -> Vec<String> {
    read_gzip_text(path)
        .lines()
        .skip(1)
        .step_by(4)
        .map(str::to_owned)
        .collect()
}
fn metrics(path: &Path, key: &str) -> u64 {
    std::fs::read_to_string(path.join("preprocessing.tsv"))
        .unwrap()
        .lines()
        .find_map(|l| l.strip_prefix(&format!("{key}\t")))
        .unwrap()
        .parse()
        .unwrap()
}
#[test]
fn trim_se_full_partial_long_adapters_no_proximal_clipping_and_qc() {
    let test = TestDir::new("trim-long");
    for length in [42, 64, 100, 128] {
        let adapter = &ADAPTER[..length];
        let full = format!("GATTACA{adapter}AAAA");
        let partial = format!("GATTACA{}", &adapter[..30]);
        let short = format!("GATTACA{}", &adapter[..8]);
        let text =
            fq("full", &full, None) + &fq("partial", &partial, None) + &fq("short", &short, None);
        let reads = test.write(&format!("{length}.fastq"), &text);
        for threads in [1, 4] {
            let output = test.child(&format!("{length}-{threads}-full"));
            run_trim(
                &reads,
                &output,
                threads,
                &["--r1-distal-sequence", adapter, "--fastq-stats"],
            );
            assert_eq!(
                sequences(&output.join("trimmed.fastq.gz")),
                ["GATTACA", partial.as_str(), short.as_str()]
            );
            assert_eq!(metrics(&output, "R1_full_adapter_reads"), 1);
            assert_eq!(metrics(&output, "R1_partial_adapter_reads"), 0);
            let output = test.child(&format!("{length}-{threads}-partial"));
            run_trim(
                &reads,
                &output,
                threads,
                &[
                    "--r1-distal-sequence",
                    adapter,
                    "--partial-pass-through",
                    "--fastq-stats",
                ],
            );
            assert_eq!(
                sequences(&output.join("trimmed.fastq.gz")),
                ["GATTACA", "GATTACA", short.as_str()]
            );
            assert_eq!(metrics(&output, "R1_full_adapter_reads"), 1);
            assert_eq!(metrics(&output, "R1_partial_adapter_reads"), 1);
            assert_eq!(metrics(&output, "R1_partial_adapter_bases_removed"), 30);
            assert_eq!(metrics(&output, "R1_proximal_bases_removed"), 0);
            assert_eq!(metrics(&output, "fragments_retained"), 3);
            let stats = std::fs::read_to_string(output.join("retained_fastq_stats.tsv")).unwrap();
            assert!(stats.contains("R1\t3\t29\t7\t15\t"), "{stats}");
            let contents = read_gzip_text(&output.join("trimmed.fastq.gz"));
            for record in contents.lines().collect::<Vec<_>>().chunks(4) {
                assert_eq!(record[1].len(), record[3].len());
            }
        }
    }
}
#[test]
fn trim_pair_independent_signatures_filters_and_order_with_parallel_gzip() {
    let test = TestDir::new("trim-pair");
    let adapter2 = rc(ADAPTER);
    let mut r1 = String::new();
    let mut r2 = String::new();
    let mut expected1 = Vec::new();
    let mut expected2 = Vec::new();
    for i in 0..2100 {
        let insert1 = if i % 3 == 0 { "GAT" } else { "GATTACA" };
        let insert2 = if i % 5 == 0 { "CC" } else { "GGTT" };
        r1.push_str(&fq(
            &format!("read{i}/1"),
            &format!("{insert1}{}", &ADAPTER[..30]),
            None,
        ));
        r2.push_str(&fq(
            &format!("read{i}/2"),
            &format!("{insert2}{}", &adapter2[..20]),
            None,
        ));
        if insert1.len() >= 4 && insert2.len() >= 4 {
            expected1.push(fq(&format!("read{i}/1"), insert1, None));
            expected2.push(fq(&format!("read{i}/2"), insert2, None));
        }
    }
    for (name, text) in [("R1", r1), ("R2", r2)] {
        let mut gzip = flate2::write::GzEncoder::new(
            std::fs::File::create(test.child(&format!("{name}.fastq.gz"))).unwrap(),
            flate2::Compression::new(2),
        );
        gzip.write_all(text.as_bytes()).unwrap();
        gzip.finish().unwrap();
    }
    for threads in [1, 4, 8] {
        let out = test.child(&format!("out-{threads}"));
        success(
            Command::new(env!("CARGO_BIN_EXE_plexless"))
                .args(["trim", "--r1"])
                .arg(test.child("R1.fastq.gz"))
                .arg("--r2")
                .arg(test.child("R2.fastq.gz"))
                .arg("--output")
                .arg(&out)
                .arg("--threads")
                .arg(threads.to_string())
                .args([
                    "--r1-distal-sequence",
                    ADAPTER,
                    "--r2-distal-sequence",
                    &adapter2,
                    "--no-overlap",
                    "--partial-pass-through",
                    "--min-length",
                    "4",
                    "--fastq-stats",
                ])
                .output()
                .unwrap(),
        );
        assert_eq!(
            read_gzip_text(&out.join("trimmed_R1.fastq.gz")),
            expected1.concat()
        );
        assert_eq!(
            read_gzip_text(&out.join("trimmed_R2.fastq.gz")),
            expected2.concat()
        );
        assert_eq!(metrics(&out, "fragments_processed"), 2100);
        assert_eq!(metrics(&out, "fragments_retained"), expected1.len() as u64);
        assert_eq!(
            metrics(&out, "fragments_rejected"),
            2100 - expected1.len() as u64
        );
        assert_eq!(metrics(&out, "R1_partial_adapter_reads"), 2100);
        assert_eq!(metrics(&out, "R2_partial_adapter_reads"), 2100);
    }
}
#[test]
fn filters_are_opt_in_and_respect_exact_thresholds() {
    let test = TestDir::new("trim-filters");
    let reads = test.write(
        "reads.fastq",
        &(fq("a", "ACGT", Some("I!II"))
            + &fq("b", "ACGN", None)
            + &fq("c", "AC", Some("!!"))
            + &fq("d", "ACGTN", None)),
    );
    let cases: Vec<(Vec<&str>, Vec<&str>)> = vec![
        (vec![], vec!["ACGT", "ACGN", "AC", "ACGTN"]),
        (vec!["--fastq-stats"], vec!["ACGT", "ACGN", "AC", "ACGTN"]),
        (vec!["--max-n", "0"], vec!["ACGT", "AC"]),
        (vec!["--max-n", "1"], vec!["ACGT", "ACGN", "AC", "ACGTN"]),
        (vec!["--min-length", "4"], vec!["ACGT", "ACGN", "ACGTN"]),
        (
            vec![
                "--qualified-quality",
                "40",
                "--max-low-quality-percent",
                "25",
            ],
            vec!["ACGT", "ACGN", "ACGTN"],
        ),
        (
            vec!["--max-low-quality-percent", "24"],
            vec!["ACGN", "ACGTN"],
        ),
        (
            vec!["--min-average-quality", "30"],
            vec!["ACGT", "ACGN", "ACGTN"],
        ),
        (vec!["--min-average-quality", "30.1"], vec!["ACGN", "ACGTN"]),
        (
            vec![
                "--qualified-quality",
                "41",
                "--max-low-quality-percent",
                "100",
            ],
            vec!["ACGT", "ACGN", "AC", "ACGTN"],
        ),
    ];
    for (i, (options, expected)) in cases.iter().enumerate() {
        for threads in [1, 4] {
            let out = test.child(&format!("out-{i}-{threads}"));
            run_trim(&reads, &out, threads, options);
            assert_eq!(
                sequences(&out.join("trimmed.fastq.gz")),
                *expected,
                "case {i}"
            );
            assert_eq!(metrics(&out, "fragments_retained"), expected.len() as u64);
            assert_eq!(
                metrics(&out, "fragments_rejected"),
                (4 - expected.len()) as u64
            );
        }
    }
}
#[allow(clippy::too_many_arguments)]
fn demux(
    test: &TestDir,
    structure1: &str,
    structure2: Option<&str>,
    r1: &str,
    r2: Option<&str>,
    technical: Option<&str>,
    threads: usize,
    label: &str,
    extra: &[&str],
) -> std::path::PathBuf {
    let barcodes = test.write(
        &format!("barcodes-{label}.tsv"),
        &format!("Set\tID\tSequence\nA\ta\t{BARCODE}\n"),
    );
    let samples = test.write(&format!("samples-{label}.tsv"), "Sample\tA\ns\ta\n");
    let reads = test.write(&format!("r1-{label}.fastq"), r1);
    let out = test.child(&format!("out-{label}-{threads}"));
    let mut command = Command::new(env!("CARGO_BIN_EXE_plexless"));
    command
        .args(["demux", "--barcodes"])
        .arg(barcodes)
        .arg("--samples")
        .arg(samples)
        .arg("--output")
        .arg(&out)
        .arg("--threads")
        .arg(threads.to_string());
    if let Some(r2) = r2 {
        command
            .arg("--r1")
            .arg(reads)
            .arg("--r2")
            .arg(test.write(&format!("r2-{label}.fastq"), r2))
            .args([
                "--r1-structure",
                structure1,
                "--r2-structure",
                structure2.unwrap(),
            ]);
    } else {
        command
            .arg("--reads")
            .arg(reads)
            .args(["--structure", structure1]);
    }
    if let Some(technical) = technical {
        command
            .arg("--technical-sequences")
            .arg(test.write(&format!("tech-{label}.tsv"), technical));
    }
    success(command.args(extra).arg("--no-overlap").output().unwrap());
    out
}
#[test]
fn known_t_multiple_segments_long_signatures_and_physical_orientation() {
    let test = TestDir::new("known-t-orientation");
    for threads in [1, 4] {
        for tech_length in [22, 108] {
            let tech = &ADAPTER[..tech_length];
            // Proximal order: first half T, forward A, remainder T. Distal
            // reverses physical segment order and reverse-complements every T.
            let half = tech_length / 2;
            let proximal = format!("{}{BARCODE}{}", &tech[..half], &tech[half..]);
            let distal = rc(&proximal);
            let raw = format!("{proximal}GATTACA{distal}");
            let structure = format!("R1_{half}T20A{}T", tech_length - half);
            let techfile = format!("Mate\tSequence\nR1\t{}\n", tech.to_ascii_lowercase());
            let label = format!("long-{tech_length}");
            let out = demux(
                &test,
                &structure,
                None,
                &fq("read", &raw, None),
                None,
                Some(&techfile),
                threads,
                &label,
                &["--pass-through-mismatches", "0"],
            );
            assert_eq!(sequences(&out.join("s.fastq.gz")), ["GATTACA"]);
            let mut wrong = distal.clone().into_bytes();
            wrong[0] = if wrong[0] == b'A' { b'C' } else { b'A' };
            let raw = format!("{proximal}GATTACA{}", String::from_utf8(wrong).unwrap());
            let label = format!("wrong-{tech_length}");
            let out = demux(
                &test,
                &structure,
                None,
                &fq("read", &raw, None),
                None,
                Some(&techfile),
                threads,
                &label,
                &["--pass-through-mismatches", "0"],
            );
            assert_eq!(sequences(&out.join("s.fastq.gz")), [&raw[proximal.len()..]]);
            let out = demux(
                &test,
                &structure,
                None,
                &fq("read", &raw, None),
                None,
                None,
                threads,
                &format!("masked-{tech_length}"),
                &["--pass-through-mismatches", "0"],
            );
            assert_eq!(sequences(&out.join("s.fastq.gz")), ["GATTACA"]);
        }
        // A(rc) is already physically complemented; opposite compilation restores A.
        let tech1 = &ADAPTER[..25];
        let tech2 = &ADAPTER[40..82];
        let proximal1 = format!(
            "{}{}{}{}",
            rc(&BARCODE[..10]),
            &tech1[..9],
            rc(&BARCODE[10..]),
            &tech1[9..]
        );
        let proximal2 = tech2.to_owned();
        let raw1 = format!("{proximal1}GATTACA{}", rc(&proximal2));
        let raw2 = format!("{proximal2}CAGT{}", rc(&proximal1));
        let techfile = format!("Mate\tSequence\nR1\t{tech1}\nR2\t{tech2}\n");
        let out = demux(
            &test,
            "R1_10A(rc)9T10A(rc)16T",
            Some("R2_42T"),
            &fq("read/1", &raw1, None),
            Some(&fq("read/2", &raw2, None)),
            Some(&techfile),
            threads,
            "paired",
            &["--pass-through-mismatches", "0"],
        );
        assert_eq!(sequences(&out.join("s_R1.fastq.gz")), ["GATTACA"]);
        assert_eq!(sequences(&out.join("s_R2.fastq.gz")), ["CAGT"]);
    }
}

#[test]
fn demux_filters_after_assignment_and_adapter_detection_serial_and_parallel() {
    let test = TestDir::new("demux-preprocessing");
    let read1 = fq(
        "full",
        &format!("{BARCODE}GATTACA{ADAPTER}"),
        Some(&format!(
            "{}{}{}",
            "!".repeat(20),
            "I".repeat(7),
            "!".repeat(ADAPTER.len())
        )),
    );
    let read2 = fq(
        "partial",
        &format!("{BARCODE}GATTACA{}", &ADAPTER[..30]),
        None,
    );
    let read3 = fq("short", &format!("{BARCODE}GA{ADAPTER}"), None);
    let unassigned = fq("unassigned", "NNNNNNNNNNNNNNNNNNNNACGTN", None);
    for threads in [1, 4] {
        for mode in ["direct", "buffered"] {
            let out = demux(
                &test,
                "R1_20A",
                None,
                &(read1.clone() + &read2 + &read3 + &unassigned),
                None,
                None,
                threads,
                mode,
                &[
                    "--r1-distal-sequence",
                    ADAPTER,
                    "--partial-pass-through",
                    "--min-length",
                    "5",
                    "--min-average-quality",
                    "30",
                    "--max-n",
                    "0",
                    "--fastq-stats",
                    "--write-unassigned",
                    "--output-mode",
                    mode,
                ],
            );
            assert_eq!(sequences(&out.join("s.fastq.gz")), ["GATTACA", "GATTACA"]);
            assert!(!out.join("unassigned.fastq.gz").exists());
            assert_eq!(metrics(&out, "fragments_processed"), 4);
            assert_eq!(metrics(&out, "fragments_rejected"), 2);
            assert_eq!(metrics(&out, "R1_full_adapter_reads"), 2);
            assert_eq!(metrics(&out, "R1_partial_adapter_reads"), 1);
        }
    }
}
#[test]
fn full_detection_wins_over_stronger_partial_at_another_position() {
    let test = TestDir::new("full-priority");
    // Both complete and incomplete copies are present. A later partial match
    // cannot replace the complete match at the first exact adapter boundary.
    let raw = format!("GATTACA{ADAPTER}GGG{}", &ADAPTER[..30]);
    let reads = test.write("reads.fastq", &fq("read", &raw, None));
    let out = test.child("out");
    run_trim(
        &reads,
        &out,
        4,
        &["--r1-distal-sequence", ADAPTER, "--partial-pass-through"],
    );
    assert_eq!(sequences(&out.join("trimmed.fastq.gz")), ["GATTACA"]);
    assert_eq!(metrics(&out, "R1_full_adapter_reads"), 1);
    assert_eq!(metrics(&out, "R1_partial_adapter_reads"), 0);
}
#[test]
fn explicit_proximal_clipping_and_distal_t_structure_in_trim() {
    let test = TestDir::new("trim-proximal");
    let adapter = &ADAPTER[..42];
    let reads = test.write(
        "reads.fastq",
        &fq("read", &format!("{adapter}GATTACA{adapter}"), None),
    );
    let tech = test.write("tech.tsv", &format!("Mate\tSequence\nR1\t{adapter}\n"));
    let out = test.child("out");
    run_trim(
        &reads,
        &out,
        4,
        &[
            "--structure",
            "R1_42T",
            "--technical-sequences",
            tech.to_str().unwrap(),
            "--r1-distal-structure",
            "R1_42T",
        ],
    );
    assert_eq!(sequences(&out.join("trimmed.fastq.gz")), ["GATTACA"]);
    assert_eq!(metrics(&out, "R1_proximal_bases_removed"), 42);
}
#[test]
fn malformed_technical_files_and_conflicting_cli_settings_are_rejected() {
    let test = TestDir::new("preprocessing-validation");
    let reads = test.write("reads.fastq", &fq("read", "GATTACA", None));
    for (i, text) in [
        "",
        "Mate Sequence\nR1 ACGT\n",
        "Mate\tSequence\nR3\tACGT\n",
        "Mate\tSequence\nR1\tACGN\n",
        "Mate\tSequence\nR1\tACGT\nR1\tACGT\n",
        "Mate\tSequence\nR1\tACG\n",
        "Mate\tSequence\nR2\tACGT\n",
        "Mate\tSequence\nR1\tACGT\textra\n",
        "Mate\tSequence\nR1\t\n",
        "Mate\tSequence\n",
    ]
    .iter()
    .enumerate()
    {
        let tech = test.write(&format!("tech{i}.tsv"), text);
        let out = test.child(&format!("out{i}"));
        let output = Command::new(env!("CARGO_BIN_EXE_plexless"))
            .args(["trim", "--reads"])
            .arg(&reads)
            .args(["--structure", "R1_4T", "--technical-sequences"])
            .arg(tech)
            .arg("--output")
            .arg(out)
            .output()
            .unwrap();
        assert!(!output.status.success(), "accepted {text:?}");
    }
    for (i, args) in [
        vec!["--r1-distal-sequence", "ACGN"],
        vec!["--r2-distal-sequence", ADAPTER],
        vec![
            "--r1-distal-sequence",
            ADAPTER,
            "--r1-distal-structure",
            "R1_20T",
        ],
        vec!["--r1-distal-sequence", ADAPTER, "--no-pass-through"],
        vec![
            "--partial-pass-through",
            "--min-partial-adapter-overlap",
            "0",
        ],
        vec!["--max-low-quality-percent", "101"],
        vec!["--min-average-quality", "NaN"],
        vec!["--r1-distal-structure", "R1_20A"],
        vec!["--r1-structure", "R1_4T"],
    ]
    .iter()
    .enumerate()
    {
        let output = Command::new(env!("CARGO_BIN_EXE_plexless"))
            .args(["trim", "--reads"])
            .arg(&reads)
            .arg("--output")
            .arg(test.child(&format!("invalid{i}")))
            .args(args)
            .output()
            .unwrap();
        assert!(!output.status.success(), "accepted {args:?}");
    }
}

#[test]
fn trim_all_rejected_produces_valid_empty_gzip_and_has_no_sample_reports() {
    let test = TestDir::new("trim-all-rejected");
    let reads = test.write("reads.fastq", &fq("read", "GATTACA", None));
    let out = test.child("out");
    run_trim(&reads, &out, 4, &["--min-length", "8"]);
    assert_eq!(read_gzip_text(&out.join("trimmed.fastq.gz")), "");
    assert_eq!(metrics(&out, "fragments_rejected"), 1);
    assert!(!out.join("sample_metrics.tsv").exists());
    assert!(!out.join("PLEXLESS_INCOMPLETE").exists());
}

#[test]
fn demux_paired_filter_rejects_both_mates_after_independent_full_partial_trimming() {
    let test = TestDir::new("demux-paired-filters");
    let adapter2 = rc(ADAPTER);
    let prefix2 = &ADAPTER[40..82];
    let r1 = fq(
        "keep/1",
        &format!("{BARCODE}GATTACA{}", &ADAPTER[..30]),
        None,
    ) + &fq("reject/1", &format!("{BARCODE}GA{ADAPTER}"), None);
    let r2 = fq("keep/2", &format!("{prefix2}CAGT{}", &adapter2[..20]), None)
        + &fq("reject/2", &format!("{prefix2}CAGT{adapter2}"), None);
    for threads in [1, 4] {
        for mode in ["direct", "buffered"] {
            let out = demux(
                &test,
                "R1_20A",
                Some("R2_42T"),
                &r1,
                Some(&r2),
                None,
                threads,
                mode,
                &[
                    "--r1-distal-sequence",
                    ADAPTER,
                    "--r2-distal-sequence",
                    &adapter2,
                    "--partial-pass-through",
                    "--min-length",
                    "4",
                    "--fastq-stats",
                    "--output-mode",
                    mode,
                ],
            );
            assert_eq!(
                read_gzip_text(&out.join("s_R1.fastq.gz")),
                fq("keep/1", "GATTACA", None)
            );
            assert_eq!(
                read_gzip_text(&out.join("s_R2.fastq.gz")),
                fq("keep/2", "CAGT", None)
            );
            assert_eq!(metrics(&out, "fragments_retained"), 1);
            assert_eq!(metrics(&out, "fragments_rejected"), 1);
            for mate in [1, 2] {
                assert_eq!(metrics(&out, &format!("R{mate}_full_adapter_reads")), 1);
                assert_eq!(metrics(&out, &format!("R{mate}_partial_adapter_reads")), 1);
            }
        }
    }
}
