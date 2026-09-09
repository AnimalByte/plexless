mod common;

use plexless::cli::DemuxArgs;

use common::{TestDir, read_gzip_text};

fn paired_args(
    test: &TestDir,
    name: &str,
    r1_text: &str,
    r2_text: &str,
    write_unassigned: bool,
) -> DemuxArgs {
    let r1 = test.write(&format!("{name}_R1.fastq"), r1_text);
    let r2 = test.write(&format!("{name}_R2.fastq"), r2_text);

    let barcodes = test.write(
        &format!("{name}_barcodes.tsv"),
        "Set\tID\tSequence\n\
         A\tA01\tACGT\n\
         B\tB01\tTGCA\n",
    );

    let samples = test.write(
        &format!("{name}_samples.tsv"),
        "Sample\tA\tB\n\
         sample_1\tA01\tB01\n",
    );

    DemuxArgs {
        reads: None,
        r1: Some(r1),
        r2: Some(r2),
        structure: None,
        r1_structure: Some("R1_2A2B2T".to_string()),
        r2_structure: Some("R2_2A2B2T".to_string()),
        barcodes,
        samples,
        output: test.child(&format!("{name}_output")),
        compression_level: 2,
        output_mode: plexless::cli::OutputMode::Buffered,
        output_chunk_size: plexless::cli::ByteSizeSetting::Auto,
        output_buffer_memory: plexless::cli::ByteSizeSetting::Auto,
        max_open_files: None,
        max_mismatches: 1,
        fastq_stats: false,
        write_unassigned,
        low_sample_fraction: 0.05,
    }
}

#[test]
fn missing_r2_mate_becomes_r1_orphan_and_processing_resynchronizes() {
    let test = TestDir::new("r1-orphan");

    let args = paired_args(
        &test,
        "case",
        "@read1/1\n\
         ACTGAAGATTACA\n\
         +\n\
         IIIIIIIIIIIII\n\
         @orphan/1\n\
         AAAAAACCCCCCC\n\
         +\n\
         IIIIIIIIIIIII\n\
         @read3/1\n\
         ACTGAATTTTTTT\n\
         +\n\
         IIIIIIIIIIIII\n",
        "@read1/2\n\
         GTCACCACAC\n\
         +\n\
         IIIIIIIIII\n\
         @read3/2\n\
         GTCACCGGGG\n\
         +\n\
         IIIIIIIIII\n",
        true,
    );

    let output = args.output.clone();

    plexless::demux::run(args).expect("R1 orphan should not abort paired demultiplexing");

    assert_eq!(
        read_gzip_text(&output.join("sample_1_R1.fastq.gz")),
        "@read1/1\nGATTACA\n+\nIIIIIII\n\
         @read3/1\nTTTTTTT\n+\nIIIIIII\n"
    );

    assert_eq!(
        read_gzip_text(&output.join("sample_1_R2.fastq.gz")),
        "@read1/2\nACAC\n+\nIIII\n\
         @read3/2\nGGGG\n+\nIIII\n"
    );

    assert_eq!(
        read_gzip_text(&output.join("unassigned_R1.fastq.gz")),
        "@orphan/1\nAAAAAACCCCCCC\n+\nIIIIIIIIIIIII\n"
    );

    assert!(!output.join("unassigned_R2.fastq.gz").exists());
}

#[test]
fn missing_r1_mate_becomes_r2_orphan_and_processing_resynchronizes() {
    let test = TestDir::new("r2-orphan");

    let args = paired_args(
        &test,
        "case",
        "@read1/1\n\
         ACTGAAGATTACA\n\
         +\n\
         IIIIIIIIIIIII\n\
         @read3/1\n\
         ACTGAATTTTTTT\n\
         +\n\
         IIIIIIIIIIIII\n",
        "@read1/2\n\
         GTCACCACAC\n\
         +\n\
         IIIIIIIIII\n\
         @orphan/2\n\
         AAAAAACCCC\n\
         +\n\
         IIIIIIIIII\n\
         @read3/2\n\
         GTCACCGGGG\n\
         +\n\
         IIIIIIIIII\n",
        true,
    );

    let output = args.output.clone();

    plexless::demux::run(args).expect("R2 orphan should not abort paired demultiplexing");

    assert_eq!(
        read_gzip_text(&output.join("unassigned_R2.fastq.gz")),
        "@orphan/2\nAAAAAACCCC\n+\nIIIIIIIIII\n"
    );

    assert_eq!(
        read_gzip_text(&output.join("sample_1_R1.fastq.gz")),
        "@read1/1\nGATTACA\n+\nIIIIIII\n\
         @read3/1\nTTTTTTT\n+\nIIIIIII\n"
    );
}

#[test]
fn trailing_r1_record_is_orphan_not_fatal_error() {
    let test = TestDir::new("trailing-r1");

    let args = paired_args(
        &test,
        "case",
        "@read1/1\n\
         ACTGAAGATTACA\n\
         +\n\
         IIIIIIIIIIIII\n\
         @trailing/1\n\
         AAAAAACCCCCCC\n\
         +\n\
         IIIIIIIIIIIII\n",
        "@read1/2\n\
         GTCACCACAC\n\
         +\n\
         IIIIIIIIII\n",
        true,
    );

    let output = args.output.clone();

    plexless::demux::run(args).expect("Trailing R1 should be treated as an orphan");

    assert_eq!(
        read_gzip_text(&output.join("unassigned_R1.fastq.gz")),
        "@trailing/1\nAAAAAACCCCCCC\n+\nIIIIIIIIIIIII\n"
    );
}

#[test]
fn trailing_r2_record_is_orphan_not_fatal_error() {
    let test = TestDir::new("trailing-r2");

    let args = paired_args(
        &test,
        "case",
        "@read1/1\n\
         ACTGAAGATTACA\n\
         +\n\
         IIIIIIIIIIIII\n",
        "@read1/2\n\
         GTCACCACAC\n\
         +\n\
         IIIIIIIIII\n\
         @trailing/2\n\
         AAAAAACCCC\n\
         +\n\
         IIIIIIIIII\n",
        true,
    );

    let output = args.output.clone();

    plexless::demux::run(args).expect("Trailing R2 should be treated as an orphan");

    assert_eq!(
        read_gzip_text(&output.join("unassigned_R2.fastq.gz")),
        "@trailing/2\nAAAAAACCCC\n+\nIIIIIIIIII\n"
    );
}

#[test]
fn orphan_recovery_does_not_require_unassigned_output() {
    let test = TestDir::new("orphan-no-output");

    let args = paired_args(
        &test,
        "case",
        "@read1/1\n\
         ACTGAAGATTACA\n\
         +\n\
         IIIIIIIIIIIII\n\
         @orphan/1\n\
         AAAAAACCCCCCC\n\
         +\n\
         IIIIIIIIIIIII\n\
         @read3/1\n\
         ACTGAATTTTTTT\n\
         +\n\
         IIIIIIIIIIIII\n",
        "@read1/2\n\
         GTCACCACAC\n\
         +\n\
         IIIIIIIIII\n\
         @read3/2\n\
         GTCACCGGGG\n\
         +\n\
         IIIIIIIIII\n",
        false,
    );

    let output = args.output.clone();

    plexless::demux::run(args)
        .expect("Orphan recovery should work even when unassigned output is disabled");

    assert!(!output.join("unassigned_R1.fastq.gz").exists());

    assert_eq!(
        read_gzip_text(&output.join("sample_1_R1.fastq.gz")),
        "@read1/1\nGATTACA\n+\nIIIIIII\n\
         @read3/1\nTTTTTTT\n+\nIIIIIII\n"
    );
}
