mod common;

use std::fs;

use plexless::cli::DemuxArgs;

use common::{TestDir, read_gzip_text};

#[allow(clippy::too_many_arguments)]
fn run_single(
    test: &TestDir,
    name: &str,
    reads_text: &str,
    barcodes_text: &str,
    samples_text: &str,
    max_mismatches: u8,
    fastq_stats: bool,
    write_unassigned: bool,
) -> std::path::PathBuf {
    let reads = test.write(&format!("{name}.fastq"), reads_text);

    let barcodes = test.write(&format!("{name}_barcodes.tsv"), barcodes_text);

    let samples = test.write(&format!("{name}_samples.tsv"), samples_text);

    let output = test.child(&format!("{name}_output"));

    let args = DemuxArgs {
        reads: Some(reads),
        r1: None,
        r2: None,
        structure: Some("R1_4A4B2T".to_string()),
        r1_structure: None,
        r2_structure: None,
        barcodes,
        samples,
        output: output.clone(),
        compression_level: 2,
        max_mismatches,
        fastq_stats,
        write_unassigned,
    };

    args.validate().expect("CLI arguments should be valid");

    plexless::demux::run(args).expect("Demultiplexing should succeed");

    output
}

#[test]
fn one_mismatch_is_corrected_and_routed() {
    let test = TestDir::new("one-mismatch");

    let output = run_single(
        &test,
        "case",
        "@read1\n\
         ACGATGCAAAGATTACA\n\
         +\n\
         IIIIIIIIIIIIIIIII\n",
        "Set\tID\tSequence\n\
         A\tA01\tACGT\n\
         B\tB01\tTGCA\n",
        "Sample\tA\tB\n\
         sample_1\tA01\tB01\n",
        1,
        false,
        false,
    );

    let observed = read_gzip_text(&output.join("sample_1.fastq.gz"));

    assert_eq!(observed, "@read1\nGATTACA\n+\nIIIIIII\n");
}

#[test]
fn one_n_is_rescued_but_two_ns_are_unmatched() {
    let test = TestDir::new("n-handling");

    let output = run_single(
        &test,
        "case",
        "@one_n\n\
         ACNTTGCAAAGATTACA\n\
         +\n\
         IIIIIIIIIIIIIIIII\n\
         @two_ns\n\
         ANNTTGCAAACCCCCCC\n\
         +\n\
         IIIIIIIIIIIIIIIII\n",
        "Set\tID\tSequence\n\
         A\tA01\tACGT\n\
         B\tB01\tTGCA\n",
        "Sample\tA\tB\n\
         sample_1\tA01\tB01\n",
        1,
        false,
        true,
    );

    let assigned = read_gzip_text(&output.join("sample_1.fastq.gz"));

    assert_eq!(assigned, "@one_n\nGATTACA\n+\nIIIIIII\n");

    let unassigned = read_gzip_text(&output.join("unassigned.fastq.gz"));

    assert_eq!(
        unassigned,
        "@two_ns\nANNTTGCAAACCCCCCC\n+\nIIIIIIIIIIIIIIIII\n"
    );
}

#[test]
fn barcode_beyond_mismatch_limit_is_unmatched_and_written_raw() {
    let test = TestDir::new("unmatched");

    let output = run_single(
        &test,
        "case",
        "@read1\n\
         TTTTTGCAAAGATTACA\n\
         +\n\
         IIIIIIIIIIIIIIIII\n",
        "Set\tID\tSequence\n\
         A\tA01\tACGT\n\
         B\tB01\tTGCA\n",
        "Sample\tA\tB\n\
         sample_1\tA01\tB01\n",
        1,
        false,
        true,
    );

    assert!(!output.join("sample_1.fastq.gz").exists());

    let unassigned = read_gzip_text(&output.join("unassigned.fastq.gz"));

    assert_eq!(
        unassigned,
        "@read1\nTTTTTGCAAAGATTACA\n+\nIIIIIIIIIIIIIIIII\n"
    );
}

#[test]
fn valid_barcodes_with_missing_sample_combination_are_unrouted() {
    let test = TestDir::new("unrouted");

    let output = run_single(
        &test,
        "case",
        "@read1\n\
         TGC ACTAG AAGATTACA\n\
         +\n\
         IIIIIIIIIIIIIIIII\n"
            .replace(' ', "")
            .as_str(),
        "Set\tID\tSequence\n\
         A\tA01\tACGT\n\
         A\tA02\tTGCA\n\
         B\tB01\tGATC\n\
         B\tB02\tCTAG\n",
        "Sample\tA\tB\n\
         sample_1\tA01\tB01\n",
        1,
        false,
        true,
    );

    assert!(!output.join("sample_1.fastq.gz").exists());

    let unassigned = read_gzip_text(&output.join("unassigned.fastq.gz"));

    assert_eq!(
        unassigned,
        "@read1\nTGCACTAGAAGATTACA\n+\nIIIIIIIIIIIIIIIII\n"
    );
}

#[test]
fn read_shorter_than_structure_is_unmatched_without_panicking() {
    let test = TestDir::new("short-read");

    let output = run_single(
        &test,
        "case",
        "@short\n\
         ACGTTGCA\n\
         +\n\
         IIIIIIII\n",
        "Set\tID\tSequence\n\
         A\tA01\tACGT\n\
         B\tB01\tTGCA\n",
        "Sample\tA\tB\n\
         sample_1\tA01\tB01\n",
        1,
        false,
        true,
    );

    let unassigned = read_gzip_text(&output.join("unassigned.fastq.gz"));

    assert_eq!(unassigned, "@short\nACGTTGCA\n+\nIIIIIIII\n");
}

#[test]
fn fastq_stats_are_calculated_from_raw_untrimmed_input() {
    let test = TestDir::new("fastq-stats");

    let output = run_single(
        &test,
        "case",
        "@read1\n\
         ACGTTGCAAAGATTACA\n\
         +\n\
         IIIIIIIIIIIIIIIII\n",
        "Set\tID\tSequence\n\
         A\tA01\tACGT\n\
         B\tB01\tTGCA\n",
        "Sample\tA\tB\n\
         sample_1\tA01\tB01\n",
        1,
        true,
        false,
    );

    let stats =
        fs::read_to_string(output.join("fastq_stats.tsv")).expect("Could not read FASTQ stats");

    let expected = "\
Mate\tReads\tBases\tMinLength\tMaxLength\tMeanLength\tGCPercent\tNPercent\tMeanQuality\tQ20Percent\tQ30Percent
R1\t1\t17\t17\t17\t17.00\t35.2941\t0.0000\t40.00\t100.0000\t100.0000
";

    assert_eq!(stats, expected);
}
