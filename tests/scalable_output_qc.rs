mod common;

use std::fs;
use std::path::PathBuf;

use plexless::cli::{ByteSizeSetting, DemuxArgs};

use common::{TestDir, read_gzip_text};

fn barcode(mut value: usize, length: usize) -> String {
    let mut result = vec![b'A'; length];
    for base in result.iter_mut().rev() {
        *base = match value % 4 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        };
        value /= 4;
    }
    String::from_utf8(result).unwrap()
}

fn common_args(
    reads: Option<PathBuf>,
    r1: Option<PathBuf>,
    r2: Option<PathBuf>,
    barcodes: PathBuf,
    samples: PathBuf,
    output: PathBuf,
) -> DemuxArgs {
    DemuxArgs {
        reads,
        r1,
        r2,
        structure: None,
        r1_structure: None,
        r2_structure: None,
        barcodes,
        samples,
        output,
        compression_level: 2,
        output_mode: plexless::cli::OutputMode::Buffered,
        output_chunk_size: ByteSizeSetting::Bytes(32 * 1024),
        output_buffer_memory: ByteSizeSetting::Bytes(4 * 1024 * 1024),
        max_open_files: Some(8),
        max_mismatches: 0,
        fastq_stats: false,
        write_unassigned: true,
        low_sample_fraction: 0.05,
    }
}

#[test]
fn metrics_include_missing_and_low_samples_and_preserve_full_headers() {
    let test = TestDir::new("sample-qc");
    let barcodes = test.write(
        "barcodes.tsv",
        "Set\tID\tSequence\nA\tA0\tAAAA\nA\tA1\tCCCC\nA\tA2\tGGGG\nA\tA3\tTTTT\n",
    );
    let samples = test.write(
        "samples.tsv",
        "Sample\tA\nsample_0\tA0\nsample_1\tA1\nsample_2\tA2\nsample_3\tA3\n",
    );
    let mut fastq = String::new();
    for (sample, count) in [(0, 100), (1, 100), (2, 4)] {
        for index in 0..count {
            let id = if sample == 0 && index == 0 {
                "instrument:1:flowcell:2:3:4 1:N:0:INDEX".to_string()
            } else {
                format!("sample{sample}-read{index} metadata=value")
            };
            let seq = format!("{}TAACCGG", ["AAAA", "CCCC", "GGGG"][sample]);
            fastq.push_str(&format!(
                "@{id}\n{seq}\n+source metadata\n{}\n",
                "I".repeat(seq.len())
            ));
        }
    }
    let reads = test.write("reads.fastq", &fastq);
    let output = test.child("output");
    let mut args = common_args(Some(reads), None, None, barcodes, samples, output.clone());
    args.structure = Some("R1_4A1T".into());
    plexless::demux::run_with_threads(args, 4).unwrap();

    let sample_0 = read_gzip_text(&output.join("sample_0.fastq.gz"));
    assert!(sample_0.starts_with("@instrument:1:flowcell:2:3:4 1:N:0:INDEX\nAACCGG\n+\n"));
    assert!(!sample_0.contains("+source metadata"));

    let metrics = fs::read_to_string(output.join("sample_metrics.tsv")).unwrap();
    assert!(metrics.contains("sample_2\t4\t4\t0\t0\t24\t0\t0\t"));
    assert!(
        metrics
            .lines()
            .any(|line| line.starts_with("sample_2\t") && line.ends_with("LOW_REPRESENTATION"))
    );
    assert!(
        metrics
            .lines()
            .any(|line| line.starts_with("sample_3\t0\t") && line.ends_with("MISSING"))
    );
    let assigned: u64 = metrics
        .lines()
        .skip(1)
        .map(|line| line.split('\t').nth(1).unwrap().parse::<u64>().unwrap())
        .sum();
    assert_eq!(assigned, 204);
}

#[test]
fn sparse_reads_across_384_single_end_outputs_remain_complete() {
    let test = TestDir::new("high-multiplex-se");
    let mut barcode_tsv = String::from("Set\tID\tSequence\n");
    let mut sample_tsv = String::from("Sample\tA\n");
    let mut fastq = String::new();
    for sample in 0..384usize {
        let sequence = barcode(sample, 5);
        barcode_tsv.push_str(&format!("A\tA{sample:03}\t{sequence}\n"));
        sample_tsv.push_str(&format!("sample_{sample:03}\tA{sample:03}\n"));
    }
    for round in 0..8usize {
        for sample in 0..384usize {
            let sequence = barcode(sample, 5);
            let seq = format!("{sequence}TACGTACGT");
            fastq.push_str(&format!(
                "@sample{sample:03}-round{round}\n{seq}\n+\n{}\n",
                "I".repeat(seq.len())
            ));
        }
    }
    let barcodes = test.write("barcodes.tsv", &barcode_tsv);
    let samples = test.write("samples.tsv", &sample_tsv);
    let reads = test.write("reads.fastq", &fastq);
    let output = test.child("output");
    let mut args = common_args(Some(reads), None, None, barcodes, samples, output.clone());
    args.structure = Some("R1_5A1T".into());
    plexless::demux::run_with_threads(args, 8).unwrap();

    for sample in 0..384usize {
        let text = read_gzip_text(&output.join(format!("sample_{sample:03}.fastq.gz")));
        assert_eq!(text.matches("\n+\n").count(), 8);
        for round in 0..8usize {
            assert!(text.contains(&format!("@sample{sample:03}-round{round}\n")));
        }
    }
}

#[test]
fn sparse_reads_across_paired_outputs_preserve_r1_r2_correspondence() {
    let test = TestDir::new("high-multiplex-pe");
    let mut barcode_tsv = String::from("Set\tID\tSequence\n");
    let mut sample_tsv = String::from("Sample\tA\n");
    let mut r1_fastq = String::new();
    let mut r2_fastq = String::new();
    for sample in 0..96usize {
        let sequence = barcode(sample, 4);
        barcode_tsv.push_str(&format!("A\tA{sample:03}\t{sequence}\n"));
        sample_tsv.push_str(&format!("sample_{sample:03}\tA{sample:03}\n"));
    }
    for round in 0..12usize {
        for sample in 0..96usize {
            let sequence = barcode(sample, 4);
            let r1_seq = format!("{sequence}TACGT");
            let r2_seq = "TGCATGCA";
            r1_fastq.push_str(&format!(
                "@sample{sample:03}-round{round}/1 meta\n{r1_seq}\n+\n{}\n",
                "I".repeat(r1_seq.len())
            ));
            r2_fastq.push_str(&format!(
                "@sample{sample:03}-round{round}/2 meta\n{r2_seq}\n+\n{}\n",
                "I".repeat(r2_seq.len())
            ));
        }
    }
    let barcodes = test.write("barcodes.tsv", &barcode_tsv);
    let samples = test.write("samples.tsv", &sample_tsv);
    let r1 = test.write("R1.fastq", &r1_fastq);
    let r2 = test.write("R2.fastq", &r2_fastq);
    let output = test.child("output");
    let mut args = common_args(None, Some(r1), Some(r2), barcodes, samples, output.clone());
    args.r1_structure = Some("R1_4A1T".into());
    plexless::demux::run_with_threads(args, 8).unwrap();

    for sample in 0..96usize {
        let r1 = read_gzip_text(&output.join(format!("sample_{sample:03}_R1.fastq.gz")));
        let r2 = read_gzip_text(&output.join(format!("sample_{sample:03}_R2.fastq.gz")));
        let r1_ids: Vec<_> = r1.lines().step_by(4).collect();
        let r2_ids: Vec<_> = r2.lines().step_by(4).collect();
        assert_eq!(r1_ids.len(), 12);
        assert_eq!(r2_ids.len(), 12);
        for (r1_id, r2_id) in r1_ids.iter().zip(r2_ids) {
            assert_eq!(r1_id.replace("/1", "/X"), r2_id.replace("/2", "/X"));
        }
    }
}
