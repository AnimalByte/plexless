mod common;

use plexless::barcodes::BarcodeCatalog;
use plexless::decoder::Decoder;
use plexless::samples::SampleSheet;
use plexless::structure::ReadLayout;

use common::TestDir;

fn abt_layout() -> ReadLayout {
    ReadLayout::single(Some("R1_4A4B2T")).expect("Test read structure should be valid")
}

fn valid_catalog(test: &TestDir, name: &str) -> BarcodeCatalog {
    let layout = abt_layout();

    let path = test.write(
        name,
        "Set\tID\tSequence\n\
         A\tA01\tACGT\n\
         A\tA02\tTGCA\n\
         B\tB01\tGATC\n\
         B\tB02\tCTAG\n",
    );

    BarcodeCatalog::load(&path, &layout).expect("Test barcode catalog should be valid")
}

fn assert_error<T>(result: Result<T, String>, context: &str) {
    match result {
        Ok(_) => panic!("{context}: expected an error"),
        Err(error) => {
            assert!(
                !error.trim().is_empty(),
                "{context}: error message should not be empty"
            );
        }
    }
}

#[test]
fn duplicate_barcode_id_is_rejected() {
    let test = TestDir::new("duplicate-barcode-id");
    let layout = abt_layout();

    let path = test.write(
        "barcodes.tsv",
        "Set\tID\tSequence\n\
         A\tA01\tACGT\n\
         A\tA01\tTGCA\n\
         B\tB01\tGATC\n",
    );

    assert_error(BarcodeCatalog::load(&path, &layout), "duplicate barcode ID");
}

#[test]
fn duplicate_barcode_sequence_is_rejected() {
    let test = TestDir::new("duplicate-barcode-sequence");
    let layout = abt_layout();

    let path = test.write(
        "barcodes.tsv",
        "Set\tID\tSequence\n\
         A\tA01\tACGT\n\
         A\tA02\tACGT\n\
         B\tB01\tGATC\n",
    );

    assert_error(
        BarcodeCatalog::load(&path, &layout),
        "duplicate barcode sequence",
    );
}

#[test]
fn barcode_length_must_match_read_structure() {
    let test = TestDir::new("barcode-length-mismatch");
    let layout = abt_layout();

    let path = test.write(
        "barcodes.tsv",
        "Set\tID\tSequence\n\
         A\tA01\tACGTA\n\
         B\tB01\tGATC\n",
    );

    assert_error(
        BarcodeCatalog::load(&path, &layout),
        "barcode length mismatch",
    );
}

#[test]
fn duplicate_sample_name_is_rejected() {
    let test = TestDir::new("duplicate-sample-name");
    let layout = abt_layout();
    let catalog = valid_catalog(&test, "barcodes.tsv");

    let samples = test.write(
        "samples.tsv",
        "Sample\tA\tB\n\
         sample_1\tA01\tB01\n\
         sample_1\tA02\tB02\n",
    );

    assert_error(
        SampleSheet::load(&samples, &layout, &catalog),
        "duplicate sample name",
    );
}

#[test]
fn duplicate_sample_barcode_combination_is_rejected() {
    let test = TestDir::new("duplicate-sample-combination");
    let layout = abt_layout();
    let catalog = valid_catalog(&test, "barcodes.tsv");

    let samples = test.write(
        "samples.tsv",
        "Sample\tA\tB\n\
         sample_1\tA01\tB01\n\
         sample_2\tA01\tB01\n",
    );

    assert_error(
        SampleSheet::load(&samples, &layout, &catalog),
        "duplicate sample barcode combination",
    );
}

#[test]
fn unknown_sample_barcode_id_is_rejected() {
    let test = TestDir::new("unknown-sample-barcode");
    let layout = abt_layout();
    let catalog = valid_catalog(&test, "barcodes.tsv");

    let samples = test.write(
        "samples.tsv",
        "Sample\tA\tB\n\
         sample_1\tA99\tB01\n",
    );

    assert_error(
        SampleSheet::load(&samples, &layout, &catalog),
        "unknown sample barcode ID",
    );
}

#[test]
fn invalid_structure_symbol_is_rejected() {
    assert_error(
        ReadLayout::single(Some("R1_4A4x2T")),
        "invalid read-structure symbol",
    );
}

#[test]
fn logical_barcode_longer_than_32_bases_is_rejected() {
    let test = TestDir::new("logical-barcode-too-long");

    let layout =
        ReadLayout::single(Some("R1_20A20A2T")).expect("Structure syntax itself should parse");

    let path = test.write(
        "barcodes.tsv",
        "Set\tID\tSequence\n\
         A\tA01\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n",
    );

    assert_error(
        BarcodeCatalog::load(&path, &layout),
        "logical barcode longer than 32 bases",
    );
}

#[test]
fn unsafe_one_mismatch_whitelist_geometry_is_rejected() {
    let test = TestDir::new("unsafe-whitelist-geometry");

    let layout = ReadLayout::single(Some("R1_4A2T")).expect("Test read structure should be valid");

    let path = test.write(
        "barcodes.tsv",
        "Set\tID\tSequence\n\
         A\tA01\tACGT\n\
         A\tA02\tACGA\n",
    );

    let catalog = BarcodeCatalog::load(&path, &layout)
        .expect("Catalog loading should succeed before distance validation");

    let set = catalog.set(b'A').expect("Barcode set A should exist");

    assert_error(Decoder::new(set, 1), "unsafe one-mismatch barcode geometry");
}

#[test]
fn close_whitelist_is_allowed_when_correction_is_disabled() {
    let test = TestDir::new("exact-only-close-whitelist");

    let layout = ReadLayout::single(Some("R1_4A2T")).expect("Test read structure should be valid");

    let path = test.write(
        "barcodes.tsv",
        "Set\tID\tSequence\n\
         A\tA01\tACGT\n\
         A\tA02\tACGA\n",
    );

    let catalog = BarcodeCatalog::load(&path, &layout).expect("Catalog should load");

    let set = catalog.set(b'A').expect("Barcode set A should exist");

    Decoder::new(set, 0).expect("Close barcodes should be valid with exact-only matching");
}
