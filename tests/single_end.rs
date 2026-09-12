mod common;

use plexless::cli::DemuxArgs;

use common::{TestDir, read_gzip_text};

#[test]
fn single_end_demux_decodes_a_and_b_and_trims_t() {
    let test = TestDir::new("single-abt");

    let barcodes = test.write(
        "barcodes.tsv",
        "Set\tID\tSequence\n\
         A\tA01\tACGT\n\
         B\tB01\tTGCA\n",
    );

    let samples = test.write(
        "samples.tsv",
        "Sample\tA\tB\n\
         sample_1\tA01\tB01\n",
    );

    let reads = test.write(
        "reads.fastq",
        "@read1\n\
         ACGTTGCAAAGATTACA\n\
         +\n\
         IIIIIIIIIIIIIIIII\n",
    );

    let output = test.child("output");

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
        max_mismatches: 1,
        fastq_stats: false,
        write_unassigned: false,
        low_sample_fraction: 0.05,
    };

    args.validate().expect("CLI arguments should be valid");

    plexless::demux::run(args).expect("Demultiplexing should succeed");

    let observed = read_gzip_text(&output.join("sample_1.fastq.gz"));

    let expected = "\
@read1
GATTACA
+
IIIIIII
";

    assert_eq!(observed, expected);
}
