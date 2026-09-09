mod common;

use plexless::cli::DemuxArgs;

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
        output_mode: plexless::cli::OutputMode::Buffered,
        output_chunk_size: plexless::cli::ByteSizeSetting::Auto,
        output_buffer_memory: plexless::cli::ByteSizeSetting::Auto,
        max_open_files: None,
        max_mismatches: 1,
        fastq_stats: false,
        write_unassigned: false,
        low_sample_fraction: 0.05,
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
    let outputs: Vec<_> = [1usize, 2, 8]
        .into_iter()
        .map(|threads| {
            let output = test.child(&format!("threads_{threads}"));
            let args = make_args(
                reads.clone(),
                barcodes.clone(),
                samples.clone(),
                output.clone(),
            );
            plexless::demux::run_with_threads(args, threads)
                .unwrap_or_else(|error| panic!("{threads}-thread demultiplexing failed: {error}"));
            output
        })
        .collect();

    for sample in ["sample_1.fastq.gz", "sample_2.fastq.gz"] {
        let expected = read_gzip_text(&outputs[0].join(sample));
        for (index, output) in outputs.iter().enumerate().skip(1) {
            assert_eq!(
                expected,
                read_gzip_text(&output.join(sample)),
                "{}-thread output differs for {sample}",
                [1usize, 2, 8][index]
            );
        }
    }
}
