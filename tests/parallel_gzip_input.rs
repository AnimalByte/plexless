mod common;

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

use flate2::Compression;
use flate2::write::GzEncoder;
use plexless::cli::DemuxArgs;

use common::{TestDir, read_gzip_text};

fn write_gzip(test: &TestDir, name: &str, text: &str) -> PathBuf {
    let path = test.child(name);
    let file = fs::File::create(&path).expect("Could not create gzip fixture");
    let mut encoder = GzEncoder::new(file, Compression::new(2));
    encoder
        .write_all(text.as_bytes())
        .expect("Could not write gzip fixture");
    encoder.finish().expect("Could not finish gzip fixture");
    path
}

fn single_args(reads: PathBuf, barcodes: PathBuf, samples: PathBuf, output: PathBuf) -> DemuxArgs {
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
        fastq_stats: true,
        write_unassigned: true,
    }
}

fn paired_args(
    r1: PathBuf,
    r2: PathBuf,
    barcodes: PathBuf,
    samples: PathBuf,
    output: PathBuf,
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

fn assert_same_text_file(left: &Path, right: &Path) {
    assert_eq!(
        fs::read_to_string(left).expect("Could not read left text file"),
        fs::read_to_string(right).expect("Could not read right text file"),
    );
}

#[test]
fn automatic_parallel_gzip_single_end_matches_serial_reference() {
    let test = TestDir::new("automatic-parallel-gzip-se");

    let barcodes = test.write(
        "barcodes.tsv",
        "Set\tID\tSequence\n\
         A\tA01\tACGT\n",
    );
    let samples = test.write(
        "samples.tsv",
        "Sample\tA\n\
         sample_1\tA01\n",
    );

    let mut fastq = String::new();

    for index in 0..3_000usize {
        let insert = if index % 2 == 0 { "GATTACA" } else { "CCCCCCC" };
        let seq = format!("ACGTAA{insert}");
        let qual = "I".repeat(seq.len());
        fastq.push_str(&format!("@read{index}\n{seq}\n+\n{qual}\n"));
    }

    let reads = write_gzip(&test, "reads.fastq.gz", &fastq);
    let serial_output = test.child("serial");
    let parallel_output = test.child("parallel");

    plexless::demux::run(single_args(
        reads.clone(),
        barcodes.clone(),
        samples.clone(),
        serial_output.clone(),
    ))
    .expect("Serial gzip demultiplexing should succeed");

    plexless::demux::run_with_threads(
        single_args(reads, barcodes, samples, parallel_output.clone()),
        8,
    )
    .expect("Automatic parallel gzip demultiplexing should succeed");

    assert_eq!(
        read_gzip_text(&serial_output.join("sample_1.fastq.gz")),
        read_gzip_text(&parallel_output.join("sample_1.fastq.gz")),
    );

    assert_same_text_file(
        &serial_output.join("fastq_stats.tsv"),
        &parallel_output.join("fastq_stats.tsv"),
    );
}

#[test]
fn automatic_parallel_gzip_paired_end_preserves_pair_resync() {
    let test = TestDir::new("automatic-parallel-gzip-pe");

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

        let r1_insert = if index % 2 == 0 { "GATT" } else { "TTTT" };
        let r2_insert = if index % 2 == 0 { "ACAC" } else { "GGGG" };
        let r1_seq = format!("ACTGAA{r1_insert}");
        let r2_seq = format!("GTCACC{r2_insert}");

        r1_text.push_str(&format!(
            "@read{index}/1\n{r1_seq}\n+\n{}\n",
            "I".repeat(r1_seq.len())
        ));
        r2_text.push_str(&format!(
            "@read{index}/2\n{r2_seq}\n+\n{}\n",
            "I".repeat(r2_seq.len())
        ));
    }

    let r1 = write_gzip(&test, "R1.fastq.gz", &r1_text);
    let r2 = write_gzip(&test, "R2.fastq.gz", &r2_text);
    let serial_output = test.child("serial");
    let parallel_output = test.child("parallel");

    plexless::demux::run(paired_args(
        r1.clone(),
        r2.clone(),
        barcodes.clone(),
        samples.clone(),
        serial_output.clone(),
    ))
    .expect("Serial paired gzip demultiplexing should succeed");

    plexless::demux::run_with_threads(
        paired_args(r1, r2, barcodes, samples, parallel_output.clone()),
        8,
    )
    .expect("Automatic parallel paired gzip demultiplexing should succeed");

    for filename in [
        "sample_1_R1.fastq.gz",
        "sample_1_R2.fastq.gz",
        "unassigned_R1.fastq.gz",
    ] {
        assert_eq!(
            read_gzip_text(&serial_output.join(filename)),
            read_gzip_text(&parallel_output.join(filename)),
            "parallel gzip output differs for {filename}"
        );
    }

    assert!(!serial_output.join("unassigned_R2.fastq.gz").exists());
    assert!(!parallel_output.join("unassigned_R2.fastq.gz").exists());

    assert_same_text_file(
        &serial_output.join("fastq_stats.tsv"),
        &parallel_output.join("fastq_stats.tsv"),
    );
}

#[test]
fn automatic_parallel_gzip_rejects_corrupt_footer() {
    let test = TestDir::new("automatic-parallel-gzip-corrupt");

    let barcodes = test.write(
        "barcodes.tsv",
        "Set\tID\tSequence\n\
         A\tA01\tACGT\n",
    );
    let samples = test.write(
        "samples.tsv",
        "Sample\tA\n\
         sample_1\tA01\n",
    );

    let mut fastq = String::new();
    for index in 0..3_000usize {
        let seq = "ACGTAAGATTACA";
        let qual = "I".repeat(seq.len());
        fastq.push_str(&format!("@read{index}\n{seq}\n+\n{qual}\n"));
    }

    let reads = write_gzip(&test, "reads.fastq.gz", &fastq);
    let mut compressed = fs::read(&reads).expect("Could not read gzip fixture");

    assert!(compressed.len() >= 8, "gzip fixture is unexpectedly short");
    let crc_index = compressed.len() - 8;
    compressed[crc_index] ^= 0xff;
    fs::write(&reads, compressed).expect("Could not corrupt gzip fixture");

    let error = plexless::demux::run_with_threads(
        single_args(reads, barcodes, samples, test.child("parallel")),
        8,
    )
    .expect_err("Corrupt gzip input must fail");

    assert!(
        !error.trim().is_empty(),
        "Corrupt gzip failure should include an error message"
    );
}
