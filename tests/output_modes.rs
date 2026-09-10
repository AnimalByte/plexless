mod common;

use std::fs;
use std::path::{Path, PathBuf};

use plexless::cli::{ByteSizeSetting, DemuxArgs, OutputMode};

use common::{TestDir, read_gzip_text};

fn args(
    reads: Option<PathBuf>,
    r1: Option<PathBuf>,
    r2: Option<PathBuf>,
    barcodes: PathBuf,
    samples: PathBuf,
    output: PathBuf,
    mode: OutputMode,
) -> DemuxArgs {
    DemuxArgs {
        reads,
        r1,
        r2,
        structure: None,
        r1_structure: None,
        r2_structure: None,
        barcodes,
        samples,
        output,
        compression_level: 2,
        output_mode: mode,
        output_chunk_size: ByteSizeSetting::Bytes(32 * 1024),
        output_buffer_memory: ByteSizeSetting::Bytes(4 * 1024 * 1024),
        max_open_files: Some(2),
        max_mismatches: 0,
        fastq_stats: true,
        write_unassigned: true,
        low_sample_fraction: 0.05,
    }
}

fn compare_outputs(direct: &Path, buffered: &Path, fastq_files: &[&str]) {
    for output in [direct, buffered] {
        assert!(
            !output.join("PLEXLESS_INCOMPLETE").exists(),
            "successful output retained its incomplete marker"
        );
    }
    for filename in fastq_files {
        assert_eq!(
            read_gzip_text(&direct.join(filename)),
            read_gzip_text(&buffered.join(filename)),
            "decompressed output differs for {filename}"
        );
    }
    for filename in ["sample_metrics.tsv", "fastq_stats.tsv", "barcode_stats.tsv"] {
        assert_eq!(
            fs::read_to_string(direct.join(filename)).unwrap(),
            fs::read_to_string(buffered.join(filename)).unwrap(),
            "metrics differ for {filename}"
        );
    }
}

#[test]
fn direct_and_buffered_single_end_outputs_and_qc_are_equivalent() {
    let test = TestDir::new("output-modes-se");
    let barcodes = test.write(
        "barcodes.tsv",
        "Set\tID\tSequence\nA\tA0\tACGT\nA\tA1\tTGCA\n",
    );
    let samples = test.write("samples.tsv", "Sample\tA\nsample_0\tA0\nsample_1\tA1\n");
    let mut fastq = String::new();
    for index in 0..2_500usize {
        let barcode = match index % 5 {
            0 => "GGGG",
            value if value % 2 == 0 => "ACGT",
            _ => "TGCA",
        };
        let sequence = format!("{barcode}TTACGTACGT");
        fastq.push_str(&format!(
            "@read{index} instrument metadata\n{sequence}\n+source metadata\n{}\n",
            "I".repeat(sequence.len())
        ));
    }
    let reads = test.write("reads.fastq", &fastq);
    for threads in [1, 4] {
        let direct = test.child(&format!("direct-{threads}"));
        let buffered = test.child(&format!("buffered-{threads}"));

        let mut direct_args = args(
            Some(reads.clone()),
            None,
            None,
            barcodes.clone(),
            samples.clone(),
            direct.clone(),
            OutputMode::Direct,
        );
        direct_args.structure = Some("R1_4A2T".into());
        plexless::demux::run_with_threads(direct_args, threads).unwrap();

        let mut buffered_args = args(
            Some(reads.clone()),
            None,
            None,
            barcodes.clone(),
            samples.clone(),
            buffered.clone(),
            OutputMode::Buffered,
        );
        buffered_args.structure = Some("R1_4A2T".into());
        plexless::demux::run_with_threads(buffered_args, threads).unwrap();

        compare_outputs(
            &direct,
            &buffered,
            &[
                "sample_0.fastq.gz",
                "sample_1.fastq.gz",
                "unassigned.fastq.gz",
            ],
        );
    }
}

#[test]
fn direct_and_buffered_paired_outputs_and_qc_are_equivalent() {
    let test = TestDir::new("output-modes-pe");
    let barcodes = test.write(
        "barcodes.tsv",
        "Set\tID\tSequence\nA\tA0\tACGT\nA\tA1\tTGCA\n",
    );
    let samples = test.write("samples.tsv", "Sample\tA\nsample_0\tA0\nsample_1\tA1\n");
    let mut r1_fastq = String::new();
    let mut r2_fastq = String::new();
    for index in 0..2_200usize {
        if index == 1_100 {
            let sequence = "ACGTTTGGGG";
            r1_fastq.push_str(&format!(
                "@orphan/1 metadata\n{sequence}\n+\n{}\n",
                "I".repeat(sequence.len())
            ));
        }
        let barcode = match index % 5 {
            0 => "GGGG",
            value if value % 2 == 0 => "ACGT",
            _ => "TGCA",
        };
        let r1_sequence = format!("{barcode}TTACGTACGT");
        let r2_sequence = "TTAACCGGTT";
        r1_fastq.push_str(&format!(
            "@read{index}/1 metadata\n{r1_sequence}\n+\n{}\n",
            "I".repeat(r1_sequence.len())
        ));
        r2_fastq.push_str(&format!(
            "@read{index}/2 metadata\n{r2_sequence}\n+\n{}\n",
            "I".repeat(r2_sequence.len())
        ));
    }
    let r1 = test.write("R1.fastq", &r1_fastq);
    let r2 = test.write("R2.fastq", &r2_fastq);
    for threads in [1, 4] {
        let direct = test.child(&format!("direct-{threads}"));
        let buffered = test.child(&format!("buffered-{threads}"));

        let mut direct_args = args(
            None,
            Some(r1.clone()),
            Some(r2.clone()),
            barcodes.clone(),
            samples.clone(),
            direct.clone(),
            OutputMode::Direct,
        );
        direct_args.r1_structure = Some("R1_4A2T".into());
        direct_args.r2_structure = Some("R2_2T".into());
        plexless::demux::run_with_threads(direct_args, threads).unwrap();

        let mut buffered_args = args(
            None,
            Some(r1.clone()),
            Some(r2.clone()),
            barcodes.clone(),
            samples.clone(),
            buffered.clone(),
            OutputMode::Buffered,
        );
        buffered_args.r1_structure = Some("R1_4A2T".into());
        buffered_args.r2_structure = Some("R2_2T".into());
        plexless::demux::run_with_threads(buffered_args, threads).unwrap();

        compare_outputs(
            &direct,
            &buffered,
            &[
                "sample_0_R1.fastq.gz",
                "sample_0_R2.fastq.gz",
                "sample_1_R1.fastq.gz",
                "sample_1_R2.fastq.gz",
                "unassigned_R1.fastq.gz",
                "unassigned_R2.fastq.gz",
            ],
        );
    }
}

#[test]
fn operational_failures_leave_incomplete_markers_in_every_output_path() {
    let test = TestDir::new("output-mode-markers");
    let barcodes = test.write("barcodes.tsv", "Set\tID\tSequence\nA\tA0\tACGT\n");
    let samples = test.write("samples.tsv", "Sample\tA\nsample_0\tA0\n");
    let reads = test.write(
        "malformed.fastq",
        "@valid metadata\nACGTTT\n+\nIIIIII\n@truncated\nACGTTT\n+\nIII\n",
    );

    for mode in [OutputMode::Direct, OutputMode::Buffered] {
        for threads in [1, 4] {
            let output = test.child(&format!("failed-{mode:?}-{threads}"));
            let mut run_args = args(
                Some(reads.clone()),
                None,
                None,
                barcodes.clone(),
                samples.clone(),
                output.clone(),
                mode,
            );
            run_args.structure = Some("R1_4A2T".into());
            let error = plexless::demux::run_with_threads(run_args, threads).unwrap_err();
            assert!(
                error.contains("FASTQ parse error"),
                "unexpected failure for {mode:?}/{threads}: {error}"
            );
            assert!(output.join("PLEXLESS_INCOMPLETE").is_file());
            assert!(!output.join("sample_metrics.tsv").exists());
            assert!(!output.join("fastq_stats.tsv").exists());
            assert!(!output.join("barcode_stats.tsv").exists());

            let mut retry_args = args(
                Some(reads.clone()),
                None,
                None,
                barcodes.clone(),
                samples.clone(),
                output,
                mode,
            );
            retry_args.structure = Some("R1_4A2T".into());
            let retry_error = plexless::demux::run_with_threads(retry_args, threads).unwrap_err();
            assert!(retry_error.contains("PLEXLESS_INCOMPLETE"));
        }
    }
}
