mod common;

use plexless::cli::DemuxArgs;

use common::{TestDir, read_gzip_text};

#[test]
fn paired_end_demux_combines_a_and_b_across_mates_and_trims_t_from_each_mate() {
    let test = TestDir::new("paired-abt");

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

    let r1 = test.write(
        "R1.fastq",
        "@read1/1\n\
         ACTGAAGATTACA\n\
         +\n\
         IIIIIIIIIIIII\n",
    );

    let r2 = test.write(
        "R2.fastq",
        "@read1/2\n\
         GTCACCACAC\n\
         +\n\
         IIIIIIIIII\n",
    );

    let output = test.child("output");

    let args = DemuxArgs {
        reads: None,
        r1: Some(r1),
        r2: Some(r2),
        structure: None,
        r1_structure: Some("R1_2A2B2T".to_string()),
        r2_structure: Some("R2_2A2B2T".to_string()),
        barcodes,
        samples,
        output: output.clone(),
        compression_level: 2,
        max_mismatches: 1,
        fastq_stats: false,
        write_unassigned: false,
    };

    args.validate().expect("CLI arguments should be valid");

    plexless::demux::run(args).expect("Demultiplexing should succeed");

    let observed_r1 = read_gzip_text(&output.join("sample_1_R1.fastq.gz"));

    let observed_r2 = read_gzip_text(&output.join("sample_1_R2.fastq.gz"));

    let expected_r1 = "\
@read1/1
GATTACA
+
IIIIIII
";

    let expected_r2 = "\
@read1/2
ACAC
+
IIII
";

    assert_eq!(observed_r1, expected_r1);
    assert_eq!(observed_r2, expected_r2);
}
