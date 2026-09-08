mod common;

use simplex::samples::{Sample, SampleSheet};
use simplex::writer::{OutputMate, WriterManager};

use common::{TestDir, read_gzip_text};

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

    let mut writer = WriterManager::new(output.clone(), &samples, false, false, 2, 2)
        .expect("Writer manager should initialize");

    writer
        .write_sample(0, OutputMate::Single, b"first", b"ACGT", b"IIII")
        .expect("First sample_0 write should succeed");

    writer
        .write_sample(1, OutputMate::Single, b"other1", b"AAAA", b"IIII")
        .expect("sample_1 write should succeed");

    writer
        .write_sample(2, OutputMate::Single, b"other2", b"CCCC", b"IIII")
        .expect("sample_2 write should succeed");

    writer
        .write_sample(0, OutputMate::Single, b"second", b"TGCA", b"IIII")
        .expect("Reopened sample_0 write should succeed");

    writer.finish().expect("Writers should finish successfully");

    let observed = read_gzip_text(&output.join("sample_0.fastq.gz"));

    let expected = "\
@first
ACGT
+
IIII
@second
TGCA
+
IIII
";

    assert_eq!(observed, expected);
}
