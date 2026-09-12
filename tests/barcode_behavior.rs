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
        cram: None,
        read_mode: None,
        output_format: None,
        structure: Some("R1_4A4B2T".to_string()),
        r1_structure: None,
        r2_structure: None,
        barcodes,
        samples,
        output: output.clone(),
        compression_level: 2,
        output_mode: plexless::cli::OutputMode::Buffered,
        output_chunk_size: plexless::cli::ByteSizeSetting::Auto,
        output_buffer_memory: plexless::cli::ByteSizeSetting::Auto,
        max_open_files: None,
        max_mismatches,
        fastq_stats,
        write_unassigned,
        low_sample_fraction: 0.05,
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
    assert!(!output.join("fastq_stats.tsv").exists());
    assert!(!output.join("barcode_stats.tsv").exists());
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
fn fastq_stats_separate_biological_and_barcode_regions() {
    let test = TestDir::new("fastq-stats");

    let output = run_single(
        &test,
        "case",
        "@read1\n\
         ACNTTGCAAAACTA\n\
         +\n\
         IIII5555!!5?I!\n\
         @read2\n\
         TGCAACGTTTGCNN\n\
         +\n\
         ++++????!!II55\n",
        "Set\tID\tSequence\n\
         A\tA01\tACGT\n\
         A\tA02\tTGCA\n\
         B\tB01\tTGCA\n\
         B\tB02\tACGT\n",
        "Sample\tA\tB\n\
         sample_1\tA01\tB01\n\
         sample_2\tA02\tB02\n",
        1,
        true,
        false,
    );

    let biological_stats =
        fs::read_to_string(output.join("fastq_stats.tsv")).expect("Could not read FASTQ stats");
    let expected_biological = "\
Mate\tReads\tBases\tMinLength\tMaxLength\tMeanLength\tGCPercent\tNPercent\tMeanQuality\tQ20Percent\tQ30Percent
R1\t2\t8\t4\t4\t4.00\t37.5000\t25.0000\t26.25\t87.5000\t50.0000
";
    assert_eq!(biological_stats, expected_biological);

    let barcode_stats =
        fs::read_to_string(output.join("barcode_stats.tsv")).expect("Could not read barcode stats");
    let expected_barcodes = "\
Mate\tBarcodeSymbol\tMatePiece\tStartCycle\tEndCycle\tOrientation\tReads\tBases\tMinLength\tMaxLength\tMeanLength\tGCPercent\tNPercent\tMeanQuality\tQ20Percent\tQ30Percent
R1\tA\t1\t1\t4\tforward\t2\t8\t4\t4\t4.00\t37.5000\t12.5000\t25.00\t50.0000\t50.0000
R1\tB\t1\t5\t8\tforward\t2\t8\t4\t4\t4.00\t50.0000\t0.0000\t25.00\t100.0000\t50.0000
";
    assert_eq!(barcode_stats, expected_barcodes);
}

#[test]
fn paired_fastq_stats_keep_barcode_mate_cycles_and_orientation() {
    let test = TestDir::new("paired-fastq-stats");
    let r1 = test.write("R1.fastq", "@pair/1 metadata\nACNGATT\n+\nII!????\n");
    let r2 = test.write("R2.fastq", "@pair/2 metadata\nACNTTNN\n+\n55+IIII\n");
    let barcodes = test.write("barcodes.tsv", "Set\tID\tSequence\nA\tA01\tACGT\n");
    let samples = test.write("samples.tsv", "Sample\tA\nsample_1\tA01\n");
    let output = test.child("output");

    let args = DemuxArgs {
        reads: None,
        r1: Some(r1),
        r2: Some(r2),
        cram: None,
        read_mode: None,
        output_format: None,
        structure: None,
        r1_structure: Some("R1_2A1T".into()),
        r2_structure: Some("R2_2A(rc)1T".into()),
        barcodes,
        samples,
        output: output.clone(),
        compression_level: 2,
        output_mode: plexless::cli::OutputMode::Buffered,
        output_chunk_size: plexless::cli::ByteSizeSetting::Auto,
        output_buffer_memory: plexless::cli::ByteSizeSetting::Auto,
        max_open_files: None,
        max_mismatches: 0,
        fastq_stats: true,
        write_unassigned: false,
        low_sample_fraction: 0.05,
    };

    plexless::demux::run_with_threads(args, 4).expect("Demultiplexing should succeed");

    let biological_stats = fs::read_to_string(output.join("fastq_stats.tsv")).unwrap();
    assert_eq!(
        biological_stats,
        "Mate\tReads\tBases\tMinLength\tMaxLength\tMeanLength\tGCPercent\tNPercent\tMeanQuality\tQ20Percent\tQ30Percent\n\
R1\t1\t4\t4\t4\t4.00\t25.0000\t0.0000\t30.00\t100.0000\t100.0000\n\
R2\t1\t4\t4\t4\t4.00\t0.0000\t50.0000\t40.00\t100.0000\t100.0000\n"
    );

    let barcode_stats = fs::read_to_string(output.join("barcode_stats.tsv")).unwrap();
    assert_eq!(
        barcode_stats,
        "Mate\tBarcodeSymbol\tMatePiece\tStartCycle\tEndCycle\tOrientation\tReads\tBases\tMinLength\tMaxLength\tMeanLength\tGCPercent\tNPercent\tMeanQuality\tQ20Percent\tQ30Percent\n\
R1\tA\t1\t1\t2\tforward\t1\t2\t2\t2\t2.00\t50.0000\t0.0000\t40.00\t100.0000\t100.0000\n\
R2\tA\t1\t1\t2\treverse-complement\t1\t2\t2\t2\t2.00\t50.0000\t0.0000\t20.00\t100.0000\t0.0000\n"
    );
}
