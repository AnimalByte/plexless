mod common;

use simplex::cli::DemuxArgs;

use common::{TestDir, read_gzip_text};

fn make_args(
    reads: std::path::PathBuf,
    barcodes: std::path::PathBuf,
    samples: std::path::PathBuf,
    output: std::path::PathBuf,
) -> DemuxArgs {
    DemuxArgs {
        reads: Some(reads),
        r1: None,
        r2: None,
        structure: Some("R1_4A2T".to_string()),
        r1_structure: None,
        r2_structure: None,
        barcodes,
        samples,
        output,
        compression_level: 2,
        max_mismatches: 1,
        fastq_stats: false,
        write_unassigned: false,
    }
}

#[test]
fn parallel_single_end_matches_serial_and_preserves_order_across_batches() {
    let test = TestDir::new("parallel-single-equivalence");

    let barcodes = test.write(
        "barcodes.tsv",
        "Set\tID\tSequence\n\
         A\tA01\tACGT\n\
         A\tA02\tTGCA\n",
    );

    let samples = test.write(
        "samples.tsv",
        "Sample\tA\n\
         sample_1\tA01\n\
         sample_2\tA02\n",
    );

    let mut fastq = String::new();

    for index in 0..3_000usize {
        let barcode = if index % 2 == 0 { "ACGT" } else { "TGCA" };
        let insert = format!("{:06}", index);
        let insert = insert
            .bytes()
            .map(|value| match value % 4 {
                0 => 'A',
                1 => 'C',
                2 => 'G',
                _ => 'T',
            })
            .collect::<String>();

        let seq = format!("{barcode}AA{insert}");
        let qual = "I".repeat(seq.len());
        fastq.push_str(&format!("@read{index}\n{seq}\n+\n{qual}\n"));
    }

    let reads = test.write("reads.fastq", &fastq);
    let serial_output = test.child("serial_output");
    let parallel_output = test.child("parallel_output");

    let serial_args = make_args(
        reads.clone(),
        barcodes.clone(),
        samples.clone(),
        serial_output.clone(),
    );
    let parallel_args = make_args(reads, barcodes, samples, parallel_output.clone());

    simplex::demux::run(serial_args).expect("Serial demultiplexing should succeed");
    simplex::demux::run_with_threads(parallel_args, 4)
        .expect("Parallel demultiplexing should succeed");

    for sample in ["sample_1.fastq.gz", "sample_2.fastq.gz"] {
        assert_eq!(
            read_gzip_text(&serial_output.join(sample)),
            read_gzip_text(&parallel_output.join(sample)),
            "parallel output differs for {sample}"
        );
    }
}
