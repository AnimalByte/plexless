//! Validation-only FASTQ-to-unmapped-CRAM fixture converter.
//!
//! This is deliberately an example target, not Plexless product functionality.
//! It exists so the deterministic validation suite can create an equivalent
//! CRAM representation without a Python HTS dependency or samtools.

use std::cmp::Ordering;
use std::path::Path;

use needletail::{FastxReader, parse_fastx_file};
use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::{self, Format, Writer};

const SINGLE_UNMAPPED: u16 = 0x4;
const PAIRED_R1_UNMAPPED: u16 = 0x1 | 0x4 | 0x8 | 0x40;
const PAIRED_R2_UNMAPPED: u16 = 0x1 | 0x4 | 0x8 | 0x80;

struct FastqRecord {
    id: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

fn main() {
    if let Err(error) = run() {
        eprintln!("Error: {error}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), String> {
    let arguments = std::env::args().skip(1).collect::<Vec<_>>();
    match arguments.as_slice() {
        [mode, input, output] if mode == "single" => convert_single(input, output),
        [mode, r1, r2, output] if mode == "paired" => convert_paired(r1, r2, output),
        _ => Err(
            "usage: make_unmapped_cram_fixture single INPUT.fastq[.gz] OUTPUT.cram\n       make_unmapped_cram_fixture paired R1.fastq[.gz] R2.fastq[.gz] OUTPUT.cram"
                .into(),
        ),
    }
}

fn validation_header(paired: bool) -> bam::Header {
    let mut header = bam::Header::new();
    header.push_record(
        HeaderRecord::new(b"HD")
            .push_tag(b"VN", "1.6")
            .push_tag(b"SO", if paired { "queryname" } else { "unsorted" }),
    );
    header.push_record(
        HeaderRecord::new(b"RG")
            .push_tag(b"ID", "plexless-validation")
            .push_tag(b"SM", "synthetic"),
    );
    header.push_record(
        HeaderRecord::new(b"PG")
            .push_tag(b"ID", "plexless-cram-fixture")
            .push_tag(b"PN", "make_unmapped_cram_fixture")
            .push_tag(b"VN", env!("CARGO_PKG_VERSION")),
    );
    header.push_comment(b"Deterministic Plexless validation fixture");
    header
}

fn convert_single(input: &str, output: &str) -> Result<(), String> {
    let mut reader = open_fastq(input)?;
    let mut writer = Writer::from_path(output, &validation_header(false), Format::Cram)
        .map_err(|error| format!("Could not create CRAM fixture '{output}': {error}"))?;
    while let Some(record) = next_record(&mut *reader, Path::new(input))? {
        writer
            .write(&to_cram_record(record, SINGLE_UNMAPPED)?)
            .map_err(|error| format!("Could not write CRAM fixture '{output}': {error}"))?;
    }
    Ok(())
}

fn convert_paired(r1: &str, r2: &str, output: &str) -> Result<(), String> {
    let mut r1_reader = open_fastq(r1)?;
    let mut r2_reader = open_fastq(r2)?;
    let mut writer = Writer::from_path(output, &validation_header(true), Format::Cram)
        .map_err(|error| format!("Could not create CRAM fixture '{output}': {error}"))?;
    let mut next_r1 = next_record(&mut *r1_reader, Path::new(r1))?;
    let mut next_r2 = next_record(&mut *r2_reader, Path::new(r2))?;

    while next_r1.is_some() || next_r2.is_some() {
        match (&next_r1, &next_r2) {
            (Some(left), Some(right)) => match core_id(&left.id).cmp(core_id(&right.id)) {
                Ordering::Equal => {
                    let left = next_r1.take().ok_or("Missing R1 fixture record")?;
                    let right = next_r2.take().ok_or("Missing R2 fixture record")?;
                    writer
                        .write(&to_cram_record(left, PAIRED_R1_UNMAPPED)?)
                        .map_err(|error| format!("Could not write CRAM R1: {error}"))?;
                    writer
                        .write(&to_cram_record(right, PAIRED_R2_UNMAPPED)?)
                        .map_err(|error| format!("Could not write CRAM R2: {error}"))?;
                    next_r1 = next_record(&mut *r1_reader, Path::new(r1))?;
                    next_r2 = next_record(&mut *r2_reader, Path::new(r2))?;
                }
                Ordering::Less => {
                    let orphan = next_r1.take().ok_or("Missing R1 fixture orphan")?;
                    writer
                        .write(&to_cram_record(orphan, PAIRED_R1_UNMAPPED)?)
                        .map_err(|error| format!("Could not write orphan R1: {error}"))?;
                    next_r1 = next_record(&mut *r1_reader, Path::new(r1))?;
                }
                Ordering::Greater => {
                    let orphan = next_r2.take().ok_or("Missing R2 fixture orphan")?;
                    writer
                        .write(&to_cram_record(orphan, PAIRED_R2_UNMAPPED)?)
                        .map_err(|error| format!("Could not write orphan R2: {error}"))?;
                    next_r2 = next_record(&mut *r2_reader, Path::new(r2))?;
                }
            },
            (Some(_), None) => {
                let orphan = next_r1.take().ok_or("Missing R1 fixture orphan")?;
                writer
                    .write(&to_cram_record(orphan, PAIRED_R1_UNMAPPED)?)
                    .map_err(|error| format!("Could not write orphan R1: {error}"))?;
                next_r1 = next_record(&mut *r1_reader, Path::new(r1))?;
            }
            (None, Some(_)) => {
                let orphan = next_r2.take().ok_or("Missing R2 fixture orphan")?;
                writer
                    .write(&to_cram_record(orphan, PAIRED_R2_UNMAPPED)?)
                    .map_err(|error| format!("Could not write orphan R2: {error}"))?;
                next_r2 = next_record(&mut *r2_reader, Path::new(r2))?;
            }
            (None, None) => break,
        }
    }
    Ok(())
}

fn open_fastq(path: &str) -> Result<Box<dyn FastxReader>, String> {
    parse_fastx_file(path)
        .map_err(|error| format!("Could not open validation FASTQ '{path}': {error}"))
}

fn next_record(reader: &mut dyn FastxReader, path: &Path) -> Result<Option<FastqRecord>, String> {
    let Some(record) = reader.next() else {
        return Ok(None);
    };
    let record = record.map_err(|error| {
        format!(
            "Could not parse validation FASTQ '{}': {error}",
            path.display()
        )
    })?;
    let qual = record
        .qual()
        .ok_or_else(|| format!("Validation input '{}' has no qualities", path.display()))?;
    Ok(Some(FastqRecord {
        id: record.id().to_vec(),
        seq: record.seq().into_owned(),
        qual: qual.to_vec(),
    }))
}

fn core_id(id: &[u8]) -> &[u8] {
    let id = id.strip_prefix(b"@").unwrap_or(id);
    let id = id
        .split(|byte| byte.is_ascii_whitespace())
        .next()
        .unwrap_or(id);
    id.strip_suffix(b"/1")
        .or_else(|| id.strip_suffix(b"/2"))
        .unwrap_or(id)
}

fn to_cram_record(record: FastqRecord, flags: u16) -> Result<bam::Record, String> {
    if record.seq.len() != record.qual.len() {
        return Err("Validation FASTQ sequence/quality lengths differ".into());
    }
    let qname = core_id(&record.id);
    let numeric_qual = record
        .qual
        .iter()
        .map(|quality| {
            quality
                .checked_sub(33)
                .filter(|value| *value <= 93)
                .ok_or_else(|| format!("Invalid FASTQ quality byte {quality}"))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let mut output = bam::Record::new();
    output.set(qname, None, &record.seq, &numeric_qual);
    output.set_flags(flags);
    output.set_tid(-1);
    output.set_pos(-1);
    output.set_mtid(-1);
    output.set_mpos(-1);
    output.set_mapq(0);
    output.set_insert_size(0);
    output
        .push_aux(b"RG", Aux::String("plexless-validation"))
        .map_err(|error| format!("Could not add validation RG tag: {error}"))?;
    Ok(output)
}
