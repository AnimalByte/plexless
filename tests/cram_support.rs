mod common;

use std::path::{Path, PathBuf};

use plexless::cli::{ByteSizeSetting, DemuxArgs, OutputFormat, OutputMode, ReadMode};
use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::record::{Aux, AuxArray, Cigar, CigarString};
use rust_htslib::bam::{self, Format, Read};

use common::{TestDir, read_gzip_text};

const SINGLE_UNMAPPED: u16 = 0x4;
const PAIRED_R1_UNMAPPED: u16 = 0x1 | 0x4 | 0x8 | 0x40;
const PAIRED_R2_UNMAPPED: u16 = 0x1 | 0x4 | 0x8 | 0x80;

fn resources(test: &TestDir) -> (PathBuf, PathBuf) {
    let barcodes = test.write(
        "barcodes.tsv",
        "Set\tID\tSequence\nA\tA01\tACGT\nA\tA02\tGGGG\n",
    );
    let samples = test.write("samples.tsv", "Sample\tA\nsample_1\tA01\n");
    (barcodes, samples)
}

fn args(
    cram: PathBuf,
    output: PathBuf,
    barcodes: PathBuf,
    samples: PathBuf,
    mode: ReadMode,
    format: OutputFormat,
) -> DemuxArgs {
    DemuxArgs {
        reads: None,
        r1: None,
        r2: None,
        cram: Some(cram),
        read_mode: Some(mode),
        output_format: Some(format),
        structure: (mode == ReadMode::Single).then(|| "R1_4A2T".into()),
        r1_structure: (mode == ReadMode::Paired).then(|| "R1_4A2T".into()),
        r2_structure: None,
        barcodes,
        samples,
        output,
        compression_level: 2,
        output_mode: OutputMode::Direct,
        output_chunk_size: ByteSizeSetting::Auto,
        output_buffer_memory: ByteSizeSetting::Auto,
        max_open_files: None,
        max_mismatches: 1,
        fastq_stats: true,
        write_unassigned: true,
        low_sample_fraction: 0.05,
    }
}

fn header(paired: bool) -> bam::Header {
    let mut header = bam::Header::new();
    let mut hd = HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", "1.6");
    hd.push_tag(b"SO", if paired { "queryname" } else { "unsorted" });
    header.push_record(&hd);
    header.push_record(
        HeaderRecord::new(b"RG")
            .push_tag(b"ID", "rg1")
            .push_tag(b"SM", "upstream"),
    );
    header.push_record(
        HeaderRecord::new(b"PG")
            .push_tag(b"ID", "upstream")
            .push_tag(b"PN", "fixture"),
    );
    header.push_comment(b"preserved fixture comment");
    header
}

fn record(name: &[u8], flags: u16, seq: &[u8], qual: &[u8]) -> bam::Record {
    assert_eq!(seq.len(), qual.len());
    let mut record = bam::Record::new();
    record.set(name, None, seq, qual);
    record.set_flags(flags);
    record.set_tid(-1);
    record.set_pos(-1);
    record.set_mtid(-1);
    record.set_mpos(-1);
    record.set_mapq(0);
    record.set_insert_size(0);
    record
}

fn write_cram(path: &Path, header: &bam::Header, records: &[bam::Record]) {
    let mut writer = bam::Writer::from_path(path, header, Format::Cram).unwrap();
    for record in records {
        writer.write(record).unwrap();
    }
    drop(writer);
}

fn read_cram(path: &Path) -> (String, Vec<bam::Record>) {
    let mut reader = bam::Reader::from_path(path).unwrap();
    let header = String::from_utf8(reader.header().as_bytes().to_vec()).unwrap();
    let records = reader.records().map(Result::unwrap).collect();
    (header, records)
}

fn write_fastq(path: &Path, records: &[&bam::Record]) {
    let mut text = Vec::new();
    for record in records {
        text.extend_from_slice(b"@");
        text.extend_from_slice(record.qname());
        text.push(b'\n');
        text.extend_from_slice(&record.seq().as_bytes());
        text.extend_from_slice(b"\n+\n");
        text.extend(record.qual().iter().map(|quality| quality + 33));
        text.push(b'\n');
    }
    std::fs::write(path, text).unwrap();
}

fn fastq_args_from_cram(
    mut args: DemuxArgs,
    reads: Option<PathBuf>,
    r1: Option<PathBuf>,
    r2: Option<PathBuf>,
) -> DemuxArgs {
    args.reads = reads;
    args.r1 = r1;
    args.r2 = r2;
    args.cram = None;
    args.read_mode = None;
    args.output_format = None;
    args
}

#[test]
fn single_cram_to_fastq_matches_fastq_routing_and_quality() {
    let test = TestDir::new("cram-se-fastq");
    let (barcodes, samples) = resources(&test);
    let records = vec![
        record(
            b"exact",
            SINGLE_UNMAPPED,
            b"ACGTAAGATTACA",
            &[0, 1, 2, 3, 40, 41, 10, 11, 12, 13, 14, 15, 16],
        ),
        record(b"correct", SINGLE_UNMAPPED, b"ACGAAACCCCC", &[30; 11]),
        record(b"observed_n", SINGLE_UNMAPPED, b"ACGNAATTTTT", &[31; 11]),
        record(b"unmatched", SINGLE_UNMAPPED, b"TTTTAAGGGGG", &[32; 11]),
        record(b"unrouted", SINGLE_UNMAPPED, b"GGGGAACACAC", &[33; 11]),
        record(b"short", SINGLE_UNMAPPED, b"ACGT", &[34; 4]),
        record(
            b"quality_bounds",
            SINGLE_UNMAPPED,
            b"ACGTAAAC",
            &[30, 30, 30, 30, 30, 30, 0, 93],
        ),
    ];
    let cram = test.child("input.cram");
    write_cram(&cram, &header(false), &records);

    let fastq = test.child("input.fastq");
    write_fastq(&fastq, &records.iter().collect::<Vec<_>>());
    let fastq_output = test.child("fastq-output");
    let fastq_args = fastq_args_from_cram(
        args(
            cram.clone(),
            fastq_output.clone(),
            barcodes.clone(),
            samples.clone(),
            ReadMode::Single,
            OutputFormat::Fastq,
        ),
        Some(fastq),
        None,
        None,
    );
    plexless::demux::run_with_threads(fastq_args, 1).unwrap();

    let output = test.child("output");
    plexless::demux::run_with_threads(
        args(
            cram,
            output.clone(),
            barcodes,
            samples,
            ReadMode::Single,
            OutputFormat::Fastq,
        ),
        4,
    )
    .unwrap();

    let assigned = read_gzip_text(&output.join("sample_1.fastq.gz"));
    assert!(assigned.contains("@exact\nGATTACA\n+\n+,-./01\n"));
    assert!(assigned.contains("@correct\nCCCCC\n+\n?????\n"));
    assert!(assigned.contains("@observed_n\nTTTTT\n+\n@@@@@\n"));
    assert!(assigned.contains("@quality_bounds\nAC\n+\n!~\n"));
    let unassigned = read_gzip_text(&output.join("unassigned.fastq.gz"));
    assert!(unassigned.contains("@unmatched\nTTTTAAGGGGG"));
    assert!(unassigned.contains("@unrouted\nGGGGAACACAC"));
    assert!(unassigned.contains("@short\nACGT"));
    for filename in [
        "sample_1.fastq.gz",
        "unassigned.fastq.gz",
        "sample_metrics.tsv",
        "fastq_stats.tsv",
        "barcode_stats.tsv",
    ] {
        let cram_bytes = if filename.ends_with(".gz") {
            read_gzip_text(&output.join(filename)).into_bytes()
        } else {
            std::fs::read(output.join(filename)).unwrap()
        };
        let fastq_bytes = if filename.ends_with(".gz") {
            read_gzip_text(&fastq_output.join(filename)).into_bytes()
        } else {
            std::fs::read(fastq_output.join(filename)).unwrap()
        };
        assert_eq!(
            cram_bytes, fastq_bytes,
            "differential mismatch in {filename}"
        );
    }
    assert!(!output.join("PLEXLESS_INCOMPLETE").exists());
}

#[test]
fn paired_cram_uses_flags_allows_both_orders_and_orphans() {
    let test = TestDir::new("cram-pe-fastq");
    let (barcodes, samples) = resources(&test);
    let records = vec![
        record(b"pair_r2_first", PAIRED_R2_UNMAPPED, b"CCCCCCC", &[20; 7]),
        record(
            b"pair_r2_first",
            PAIRED_R1_UNMAPPED,
            b"ACGTAAGGGGG",
            &[21; 11],
        ),
        record(
            b"pair_r1_first",
            PAIRED_R1_UNMAPPED,
            b"ACGTAATTTTT",
            &[22; 11],
        ),
        record(b"pair_r1_first", PAIRED_R2_UNMAPPED, b"AAAAAAA", &[23; 7]),
        record(b"orphan_r1", PAIRED_R1_UNMAPPED, b"ACGTAACCCCC", &[24; 11]),
        record(b"orphan_r2", PAIRED_R2_UNMAPPED, b"GGGGGGG", &[25; 7]),
    ];
    let cram = test.child("input.cram");
    write_cram(&cram, &header(true), &records);

    let fastq_r1 = test.child("input_R1.fastq");
    let fastq_r2 = test.child("input_R2.fastq");
    write_fastq(&fastq_r1, &[&records[1], &records[2], &records[4]]);
    write_fastq(&fastq_r2, &[&records[0], &records[3], &records[5]]);
    let fastq_output = test.child("fastq-output");
    let fastq_args = fastq_args_from_cram(
        args(
            cram.clone(),
            fastq_output.clone(),
            barcodes.clone(),
            samples.clone(),
            ReadMode::Paired,
            OutputFormat::Fastq,
        ),
        None,
        Some(fastq_r1),
        Some(fastq_r2),
    );
    plexless::demux::run_with_threads(fastq_args, 1).unwrap();
    let output = test.child("output");
    plexless::demux::run_with_threads(
        args(
            cram,
            output.clone(),
            barcodes,
            samples,
            ReadMode::Paired,
            OutputFormat::Fastq,
        ),
        4,
    )
    .unwrap();

    let r1 = read_gzip_text(&output.join("sample_1_R1.fastq.gz"));
    let r2 = read_gzip_text(&output.join("sample_1_R2.fastq.gz"));
    assert!(r1.starts_with("@pair_r2_first\nGGGGG"));
    assert!(r1.contains("@pair_r1_first\nTTTTT"));
    assert!(r2.starts_with("@pair_r2_first\nCCCCCCC"));
    assert!(r2.contains("@pair_r1_first\nAAAAAAA"));
    assert!(read_gzip_text(&output.join("unassigned_R1.fastq.gz")).contains("@orphan_r1"));
    assert!(read_gzip_text(&output.join("unassigned_R2.fastq.gz")).contains("@orphan_r2"));
    for filename in [
        "sample_1_R1.fastq.gz",
        "sample_1_R2.fastq.gz",
        "unassigned_R1.fastq.gz",
        "unassigned_R2.fastq.gz",
        "sample_metrics.tsv",
        "fastq_stats.tsv",
        "barcode_stats.tsv",
    ] {
        let cram_bytes = if filename.ends_with(".gz") {
            read_gzip_text(&output.join(filename)).into_bytes()
        } else {
            std::fs::read(output.join(filename)).unwrap()
        };
        let fastq_bytes = if filename.ends_with(".gz") {
            read_gzip_text(&fastq_output.join(filename)).into_bytes()
        } else {
            std::fs::read(fastq_output.join(filename)).unwrap()
        };
        assert_eq!(
            cram_bytes, fastq_bytes,
            "differential mismatch in {filename}"
        );
    }
}

#[test]
fn cram_to_cram_preserves_and_transforms_metadata_by_policy() {
    let test = TestDir::new("cram-metadata");
    let (barcodes, samples) = resources(&test);
    let seq = b"ACGTAAGATTACA";
    let qual = vec![40u8; seq.len()];
    let mut input = record(b"metadata", SINGLE_UNMAPPED | 0x200 | 0x400, seq, &qual);
    for (tag, value) in [
        (b"RG", "rg1"),
        (b"BC", "ACGT"),
        (b"QT", "IIII"),
        (b"RX", "UMI"),
        (b"QX", "HHH"),
        (b"MI", "molecule-1"),
        (b"ZZ", "custom"),
        (b"OQ", "IIIIIIIIIIIII"),
        (b"BQ", "@@@@@@@@@@@@@"),
        (b"E2", "TGCATTCTAATGT"),
        (b"U2", "HHHHHHHHHHHHH"),
        (b"OA", "unused-reference,1,+,13M,60,0;"),
    ] {
        input.push_aux(tag, Aux::String(value)).unwrap();
    }
    input.push_aux(b"MM", Aux::String("C+m,0;")).unwrap();
    let probabilities = [200u8];
    input
        .push_aux(b"ML", Aux::ArrayU8(AuxArray::from(&probabilities[..])))
        .unwrap();
    input.push_aux(b"MN", Aux::I32(seq.len() as i32)).unwrap();
    input.push_aux(b"MD", Aux::String("13")).unwrap();
    input.push_aux(b"NM", Aux::I32(0)).unwrap();
    input
        .push_aux(b"SA", Aux::String("chr1,1,+,13M,60,0;"))
        .unwrap();
    input.push_aux(b"PT", Aux::String("1;4;+;note")).unwrap();
    input.push_aux(b"CC", Aux::String("chr1")).unwrap();
    input.push_aux(b"CP", Aux::I32(1)).unwrap();

    let cram = test.child("input.cram");
    let mut metadata_header = header(false);
    metadata_header.push_record(
        HeaderRecord::new(b"SQ")
            .push_tag(b"SN", "unused-reference")
            .push_tag(b"LN", 1000),
    );
    write_cram(&cram, &metadata_header, &[input]);
    let output = test.child("output");
    plexless::demux::run_with_threads(
        args(
            cram,
            output.clone(),
            barcodes,
            samples,
            ReadMode::Single,
            OutputFormat::Cram,
        ),
        8,
    )
    .unwrap();

    let (out_header, records) = read_cram(&output.join("sample_1.cram"));
    assert!(out_header.contains("@RG\tID:rg1\tSM:upstream"));
    assert!(out_header.contains("@SQ\tSN:unused-reference\tLN:1000"));
    assert!(out_header.contains("@PG\tID:upstream\tPN:fixture"));
    assert!(out_header.contains("@PG\tID:plexless\tPN:plexless\tVN:"));
    assert!(out_header.contains("@CO\tpreserved fixture comment"));
    assert_eq!(records.len(), 1);
    let record = &records[0];
    assert_eq!(record.qname(), b"metadata");
    assert_eq!(record.flags(), SINGLE_UNMAPPED | 0x200 | 0x400);
    assert_eq!(record.seq().as_bytes(), b"GATTACA");
    assert_eq!(record.qual(), &[40; 7]);
    for (tag, expected) in [
        (b"RG", "rg1"),
        (b"BC", "ACGT"),
        (b"QT", "IIII"),
        (b"RX", "UMI"),
        (b"QX", "HHH"),
        (b"MI", "molecule-1"),
        (b"ZZ", "custom"),
        (b"OQ", "IIIIIII"),
        (b"BQ", "@@@@@@@"),
        (b"E2", "CTAATGT"),
        (b"U2", "HHHHHHH"),
        (b"OA", "unused-reference,1,+,13M,60,0;"),
    ] {
        assert_eq!(record.aux(tag).unwrap(), Aux::String(expected));
    }
    for tag in [
        b"MM", b"ML", b"MN", b"MD", b"NM", b"SA", b"PT", b"CC", b"CP",
    ] {
        assert!(record.aux(tag).is_err());
    }
    let metadata_report = std::fs::read_to_string(output.join("cram_metadata.tsv")).unwrap();
    assert!(metadata_report.contains("rewritten\tOQ\t1"));
    assert!(metadata_report.contains("rewritten\tE2\t1"));
    assert!(metadata_report.contains("removed\tMM\t1"));
    assert!(metadata_report.contains("removed\tMD\t1"));
    assert!(metadata_report.contains("preserved\tZZ\t1"));
    assert!(!output.join("PLEXLESS_INCOMPLETE").exists());
}

#[test]
fn paired_cram_rewrites_mate_sequence_and_quality_tags() {
    let test = TestDir::new("cram-mate-tags");
    let (barcodes, samples) = resources(&test);
    let mut r1 = record(b"pair", PAIRED_R1_UNMAPPED, b"ACGTAAGGGGG", &[30; 11]);
    let mut r2 = record(b"pair", PAIRED_R2_UNMAPPED, b"TTTTTTT", &[31; 7]);
    r1.push_aux(b"R2", Aux::String("TTTTTTT")).unwrap();
    r1.push_aux(b"Q2", Aux::String("@@@@@@@")).unwrap();
    r2.push_aux(b"R2", Aux::String("ACGTAAGGGGG")).unwrap();
    r2.push_aux(b"Q2", Aux::String("???????????")).unwrap();
    let cram = test.child("input.cram");
    write_cram(&cram, &header(true), &[r1, r2]);
    let output = test.child("output");
    plexless::demux::run_with_threads(
        args(
            cram,
            output.clone(),
            barcodes,
            samples,
            ReadMode::Paired,
            OutputFormat::Cram,
        ),
        8,
    )
    .unwrap();
    let (_, records) = read_cram(&output.join("sample_1.cram"));
    assert_eq!(records[0].aux(b"R2").unwrap(), Aux::String("TTTTTTT"));
    assert_eq!(records[0].aux(b"Q2").unwrap(), Aux::String("@@@@@@@"));
    assert_eq!(records[1].aux(b"R2").unwrap(), Aux::String("GGGGG"));
    assert_eq!(records[1].aux(b"Q2").unwrap(), Aux::String("?????"));
}

#[test]
fn invalid_metadata_rewrite_type_fails_and_retains_incomplete_marker() {
    let test = TestDir::new("cram-invalid-rewrite-type");
    let (barcodes, samples) = resources(&test);
    let mut input = record(b"invalid-oq", SINGLE_UNMAPPED, b"ACGTAAGATTACA", &[30; 13]);
    input.push_aux(b"OQ", Aux::I32(30)).unwrap();
    let cram = test.child("input.cram");
    write_cram(&cram, &header(false), &[input]);

    let output = test.child("output");
    let error = plexless::demux::run_with_threads(
        args(
            cram,
            output.clone(),
            barcodes,
            samples,
            ReadMode::Single,
            OutputFormat::Cram,
        ),
        8,
    )
    .unwrap_err();

    assert!(
        error.contains("CRAM tag OQ must have SAM type Z"),
        "{error}"
    );
    assert!(output.join("PLEXLESS_INCOMPLETE").is_file());
    assert!(!output.join("sample_metrics.tsv").exists());
}

#[test]
fn rejects_aligned_and_malformed_paired_cram() {
    let test = TestDir::new("cram-invalid");
    let (barcodes, samples) = resources(&test);
    let aligned_seq = b"ACGTAAGGGGG";
    let mut aligned = record(b"aligned", 0, aligned_seq, &[30; 11]);
    aligned.set_tid(0);
    aligned.set_pos(10);
    aligned.set(
        b"aligned",
        Some(&CigarString(vec![Cigar::Match(aligned_seq.len() as u32)])),
        aligned_seq,
        &[30; 11],
    );
    let aligned_path = test.child("aligned.cram");
    let mut aligned_header = header(false);
    aligned_header.push_record(
        HeaderRecord::new(b"SQ")
            .push_tag(b"SN", "chr1")
            .push_tag(b"LN", 1000),
    );
    write_cram(&aligned_path, &aligned_header, &[aligned]);
    let error = plexless::demux::run_with_threads(
        args(
            aligned_path,
            test.child("aligned-output"),
            barcodes.clone(),
            samples.clone(),
            ReadMode::Single,
            OutputFormat::Fastq,
        ),
        2,
    )
    .unwrap_err();
    assert!(error.contains("aligned CRAM is not supported"));

    let neither = record(b"neither", 0x1 | 0x4 | 0x8, b"ACGTAAGGGGG", &[30; 11]);
    let neither_path = test.child("neither.cram");
    write_cram(&neither_path, &header(true), &[neither]);
    let error = plexless::demux::run_with_threads(
        args(
            neither_path,
            test.child("neither-output"),
            barcodes,
            samples,
            ReadMode::Paired,
            OutputFormat::Fastq,
        ),
        2,
    )
    .unwrap_err();
    assert!(error.contains("exactly one of READ1 or READ2"));
}

#[test]
fn rejects_all_v1_record_contract_violations_and_retains_runtime_marker() {
    let test = TestDir::new("cram-invalid-contract");
    let (barcodes, samples) = resources(&test);

    let cases = vec![
        (
            "secondary",
            ReadMode::Single,
            header(false),
            vec![record(
                b"secondary",
                SINGLE_UNMAPPED | 0x100,
                b"ACGTAAG",
                &[30; 7],
            )],
            "secondary",
        ),
        (
            "supplementary",
            ReadMode::Single,
            header(false),
            vec![record(
                b"supplementary",
                SINGLE_UNMAPPED | 0x800,
                b"ACGTAAG",
                &[30; 7],
            )],
            "supplementary",
        ),
        (
            "both-read-flags",
            ReadMode::Paired,
            header(true),
            vec![record(
                b"both",
                PAIRED_R1_UNMAPPED | 0x80,
                b"ACGTAAG",
                &[30; 7],
            )],
            "exactly one of READ1 or READ2",
        ),
        (
            "duplicate-r1",
            ReadMode::Paired,
            header(true),
            vec![
                record(b"duplicate", PAIRED_R1_UNMAPPED, b"ACGTAAG", &[30; 7]),
                record(b"duplicate", PAIRED_R1_UNMAPPED, b"ACGTAAT", &[30; 7]),
            ],
            "duplicate primary R1",
        ),
        (
            "duplicate-r2",
            ReadMode::Paired,
            header(true),
            vec![
                record(b"duplicate", PAIRED_R2_UNMAPPED, b"ACGTAAG", &[30; 7]),
                record(b"duplicate", PAIRED_R2_UNMAPPED, b"ACGTAAT", &[30; 7]),
            ],
            "duplicate primary R2",
        ),
        (
            "single-declared-paired-record",
            ReadMode::Single,
            header(false),
            vec![record(b"paired", PAIRED_R1_UNMAPPED, b"ACGTAAG", &[30; 7])],
            "--read-mode single",
        ),
        (
            "paired-declared-single-record",
            ReadMode::Paired,
            header(true),
            vec![record(b"single", SINGLE_UNMAPPED, b"ACGTAAG", &[30; 7])],
            "not flagged paired",
        ),
        (
            "missing-sequence",
            ReadMode::Single,
            header(false),
            vec![record(b"missing-seq", SINGLE_UNMAPPED, b"", &[])],
            "has no sequence",
        ),
        (
            "missing-quality",
            ReadMode::Single,
            header(false),
            vec![record(
                b"missing-qual",
                SINGLE_UNMAPPED,
                b"ACGTAAG",
                &[255; 7],
            )],
            "missing quality",
        ),
        (
            "quality-too-high",
            ReadMode::Single,
            header(false),
            vec![record(b"high-qual", SINGLE_UNMAPPED, b"ACGTAAG", &[94; 7])],
            "supported SAM/FASTQ range is 0..=93",
        ),
    ];

    for (name, mode, fixture_header, records, expected) in cases {
        let path = test.child(&format!("{name}.cram"));
        write_cram(&path, &fixture_header, &records);
        let output = test.child(&format!("{name}-output"));
        let error = plexless::demux::run_with_threads(
            args(
                path,
                output.clone(),
                barcodes.clone(),
                samples.clone(),
                mode,
                OutputFormat::Fastq,
            ),
            2,
        )
        .unwrap_err();
        assert!(error.contains(expected), "{name}: {error}");
        assert!(
            output.join("PLEXLESS_INCOMPLETE").is_file(),
            "{name} must retain the incomplete marker"
        );
    }

    let incompatible_path = test.child("coordinate.cram");
    write_cram(
        &incompatible_path,
        &header(false),
        &[record(b"pair", PAIRED_R1_UNMAPPED, b"ACGTAAG", &[30; 7])],
    );
    let output = test.child("coordinate-output");
    let error = plexless::demux::run_with_threads(
        args(
            incompatible_path,
            output.clone(),
            barcodes,
            samples,
            ReadMode::Paired,
            OutputFormat::Fastq,
        ),
        2,
    )
    .unwrap_err();
    assert!(error.contains("queryname-grouped"));
    assert!(
        !output.exists(),
        "startup header validation precedes output creation"
    );
}

#[test]
fn corrupt_cram_decode_failure_retains_incomplete_marker() {
    let test = TestDir::new("cram-corrupt");
    let (barcodes, samples) = resources(&test);
    let valid_path = test.child("valid.cram");
    let records = (0..256)
        .map(|index| {
            record(
                format!("read-{index:04}").as_bytes(),
                SINGLE_UNMAPPED,
                b"ACGTAAGATTACA",
                &[30; 13],
            )
        })
        .collect::<Vec<_>>();
    write_cram(&valid_path, &header(false), &records);
    let mut corrupt_bytes = std::fs::read(&valid_path).unwrap();
    corrupt_bytes.truncate(corrupt_bytes.len() / 2);
    let corrupt_path = test.child("corrupt.cram");
    std::fs::write(&corrupt_path, corrupt_bytes).unwrap();

    let output = test.child("output");
    let error = plexless::demux::run_with_threads(
        args(
            corrupt_path,
            output.clone(),
            barcodes,
            samples,
            ReadMode::Single,
            OutputFormat::Cram,
        ),
        4,
    )
    .unwrap_err();
    assert!(
        error.contains("CRAM decode error") || error.contains("CRAM input producer"),
        "{error}"
    );
    assert!(output.join("PLEXLESS_INCOMPLETE").is_file());
    assert!(!output.join("sample_metrics.tsv").exists());
}
