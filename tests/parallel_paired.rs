mod common;

use simplex::cli::DemuxArgs;

use common::{TestDir, read_gzip_text};

fn make_args(
    r1: std::path::PathBuf,
    r2: std::path::PathBuf,
    barcodes: std::path::PathBuf,
    samples: std::path::PathBuf,
    output: std::path::PathBuf,
) -> DemuxArgs {
    DemuxArgs {
        reads: None,
        r1: Some(r1),
        r2: Some(r2),
        structure: None,
        r1_structure: Some("R1_2A2B2T".to_string()),
        r2_structure: Some("R2_2A2B2T".to_string()),
        barcodes,
        samples,
        output,
        compression_level: 2,
        max_mismatches: 1,
        fastq_stats: true,
        write_unassigned: true,
    }
}

#[test]
fn parallel_paired_matches_serial_with_orphan_resynchronization() {
    let test = TestDir::new("parallel-paired-equivalence");

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

    let mut r1_text = String::new();
    let mut r2_text = String::new();

    for index in 0..2_200usize {
        if index == 1_100 {
            let orphan_seq = "AAAAAACCCC";
            r1_text.push_str(&format!(
                "@orphan/1\n{orphan_seq}\n+\n{}\n",
                "I".repeat(orphan_seq.len())
            ));
        }

        let insert_r1 = if index % 2 == 0 { "GATT" } else { "TTTT" };
        let insert_r2 = if index % 2 == 0 { "ACAC" } else { "GGGG" };
        let seq_r1 = format!("ACTGAA{insert_r1}");
        let seq_r2 = format!("GTCACC{insert_r2}");

        r1_text.push_str(&format!(
            "@read{index}/1\n{seq_r1}\n+\n{}\n",
            "I".repeat(seq_r1.len())
        ));
        r2_text.push_str(&format!(
            "@read{index}/2\n{seq_r2}\n+\n{}\n",
            "I".repeat(seq_r2.len())
        ));
    }

    let r1 = test.write("R1.fastq", &r1_text);
    let r2 = test.write("R2.fastq", &r2_text);
    let serial_output = test.child("serial_output");
    let parallel_output = test.child("parallel_output");

    let serial_args = make_args(
        r1.clone(),
        r2.clone(),
        barcodes.clone(),
        samples.clone(),
        serial_output.clone(),
    );
    let parallel_args = make_args(r1, r2, barcodes, samples, parallel_output.clone());

    simplex::demux::run(serial_args).expect("Serial paired demultiplexing should succeed");
    simplex::demux::run_with_threads(parallel_args, 4)
        .expect("Parallel paired demultiplexing should succeed");

    for filename in [
        "sample_1_R1.fastq.gz",
        "sample_1_R2.fastq.gz",
        "unassigned_R1.fastq.gz",
        "fastq_stats.tsv",
    ] {
        let serial_path = serial_output.join(filename);
        let parallel_path = parallel_output.join(filename);

        if filename.ends_with(".gz") {
            assert_eq!(
                read_gzip_text(&serial_path),
                read_gzip_text(&parallel_path),
                "parallel output differs for {filename}"
            );
        } else {
            assert_eq!(
                std::fs::read_to_string(&serial_path).expect("Could not read serial stats"),
                std::fs::read_to_string(&parallel_path).expect("Could not read parallel stats"),
                "parallel stats differ"
            );
        }
    }

    assert!(!serial_output.join("unassigned_R2.fastq.gz").exists());
    assert!(!parallel_output.join("unassigned_R2.fastq.gz").exists());
}
