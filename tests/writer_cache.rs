mod common;

use std::io::Write;

use flate2::Compression;
use flate2::write::GzEncoder;
use plexless::samples::{Sample, SampleSheet};
use plexless::writer::{CompressedWriterManager, OutputMate};

use common::{TestDir, read_gzip_text};

fn gzip_member(fastq: &str) -> Vec<u8> {
    let mut encoder = GzEncoder::new(Vec::new(), Compression::new(2));
    encoder.write_all(fastq.as_bytes()).unwrap();
    encoder.finish().unwrap()
}

#[test]
fn writer_cache_reopen_produces_valid_concatenated_gzip_members() {
    let test = TestDir::new("writer-cache");

    let samples = SampleSheet {
        samples: vec![
            Sample {
                name: "sample_0".to_string(),
                barcode_ids: Vec::new(),
            },
            Sample {
                name: "sample_1".to_string(),
                barcode_ids: Vec::new(),
            },
            Sample {
                name: "sample_2".to_string(),
                barcode_ids: Vec::new(),
            },
        ],
    };

    let output = test.child("output");

    let mut writer = CompressedWriterManager::new(output.clone(), &samples, false, false, 2, 2)
        .expect("Writer manager should initialize");

    for round in 0..10 {
        for sample in 0..3 {
            writer
                .append_sample_member(
                    sample,
                    OutputMate::Single,
                    &gzip_member(&format!("@sample{sample}-round{round}\nACGT\n+\nIIII\n")),
                )
                .expect("Evicted output should reopen for append");
        }
    }

    writer.finish().expect("Writers should finish successfully");

    for sample in 0..3 {
        let observed = read_gzip_text(&output.join(format!("sample_{sample}.fastq.gz")));
        let expected: String = (0..10)
            .map(|round| format!("@sample{sample}-round{round}\nACGT\n+\nIIII\n"))
            .collect();
        assert_eq!(observed, expected);
    }
}

#[test]
fn reserved_unassigned_sample_name_is_rejected_before_output_creation() {
    let test = TestDir::new("writer-reserved-name");
    let samples = SampleSheet {
        samples: vec![Sample {
            name: "unassigned".to_string(),
            barcode_ids: Vec::new(),
        }],
    };
    let output = test.child("output");

    let error = CompressedWriterManager::new(output.clone(), &samples, false, true, 2, 2)
        .err()
        .expect("Reserved output name should be rejected");

    assert!(error.contains("reserved unassigned output filename"));
    assert!(!output.exists());
}
