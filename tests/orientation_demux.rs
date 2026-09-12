mod common;

use plexless::cli::DemuxArgs;

use common::{TestDir, read_gzip_text};

#[test]
fn paired_reverse_complement_piece_routes_and_trims_each_mate() {
    let test = TestDir::new("paired-rc");
    let barcodes = test.write("barcodes.tsv", "Set\tID\tSequence\nA\tA1\tACGTATGG\n");
    let samples = test.write("samples.tsv", "Sample\tA\nsample1\tA1\n");
    let r1 = test.write("R1.fastq", "@read/1\nACGTAAGATT\n+\nIIIIIIIIII\n");
    let r2 = test.write("R2.fastq", "@read/2\nCCATCCACAC\n+\nIIIIIIIIII\n");
    let output = test.child("output");

    plexless::demux::run(DemuxArgs {
        reads: None,
        r1: Some(r1),
        r2: Some(r2),
        cram: None,
        read_mode: None,
        output_format: None,
        structure: None,
        r1_structure: Some("R1_4A2T".into()),
        r2_structure: Some("R2_4A(rc)2T".into()),
        barcodes,
        samples,
        output: output.clone(),
        compression_level: 2,
        output_mode: plexless::cli::OutputMode::Buffered,
        output_chunk_size: plexless::cli::ByteSizeSetting::Auto,
        output_buffer_memory: plexless::cli::ByteSizeSetting::Auto,
        max_open_files: None,
        max_mismatches: 0,
        fastq_stats: false,
        write_unassigned: false,
        low_sample_fraction: 0.05,
    })
    .expect("Reverse-complement demultiplexing should succeed");

    assert_eq!(
        read_gzip_text(&output.join("sample1_R1.fastq.gz")),
        "@read/1\nGATT\n+\nIIII\n"
    );
    assert_eq!(
        read_gzip_text(&output.join("sample1_R2.fastq.gz")),
        "@read/2\nACAC\n+\nIIII\n"
    );
}

#[test]
fn asymmetric_paired_layout_assembles_repeated_a_and_trims_independently() {
    let test = TestDir::new("asymmetric-paired");
    let barcodes = test.write(
        "barcodes.tsv",
        "Set\tID\tSequence\nA\tA1\tACTG\nB\tB1\tGT\n",
    );
    let samples = test.write("samples.tsv", "Sample\tA\tB\nsample1\tA1\tB1\n");
    let r1 = test.write("R1.fastq", "@read/1\nACGTGATT\n+\nIIIIIIII\n");
    // reverse-complement(CA) = TG, the R2 half of logical A.
    let r2 = test.write("R2.fastq", "@read/2\nCAACAC\n+\nIIIIII\n");
    let output = test.child("output");

    plexless::demux::run(DemuxArgs {
        reads: None,
        r1: Some(r1),
        r2: Some(r2),
        cram: None,
        read_mode: None,
        output_format: None,
        structure: None,
        r1_structure: Some("R1_2A2B".into()),
        r2_structure: Some("R2_2A(rc)".into()),
        barcodes,
        samples,
        output: output.clone(),
        compression_level: 2,
        output_mode: plexless::cli::OutputMode::Buffered,
        output_chunk_size: plexless::cli::ByteSizeSetting::Auto,
        output_buffer_memory: plexless::cli::ByteSizeSetting::Auto,
        max_open_files: None,
        max_mismatches: 0,
        fastq_stats: false,
        write_unassigned: false,
        low_sample_fraction: 0.05,
    })
    .expect("Asymmetric paired demultiplexing should succeed");

    assert_eq!(
        read_gzip_text(&output.join("sample1_R1.fastq.gz")),
        "@read/1\nGATT\n+\nIIII\n"
    );
    assert_eq!(
        read_gzip_text(&output.join("sample1_R2.fastq.gz")),
        "@read/2\nACAC\n+\nIIII\n"
    );
}
