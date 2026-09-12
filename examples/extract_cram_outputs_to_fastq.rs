//! Validation-only decoder for semantic CRAM output comparison.
//!
//! It converts every CRAM in a Plexless output directory into deterministic
//! gzip FASTQ streams.  The deep validator then compares those streams to the
//! same truth tables used for native FASTQ output.

use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use flate2::Compression;
use flate2::write::GzEncoder;
use rust_htslib::bam::{self, Read};

type FastqWriter = BufWriter<GzEncoder<File>>;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Mode {
    Single,
    Paired,
}

fn main() {
    if let Err(error) = run() {
        eprintln!("Error: {error}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), String> {
    let arguments = std::env::args().skip(1).collect::<Vec<_>>();
    let [mode, input, output] = arguments.as_slice() else {
        return Err("usage: extract_cram_outputs_to_fastq single|paired CRAM_DIR FASTQ_DIR".into());
    };
    let mode = match mode.as_str() {
        "single" => Mode::Single,
        "paired" => Mode::Paired,
        _ => return Err("mode must be single or paired".into()),
    };
    let input = Path::new(input);
    let output = Path::new(output);
    if output.exists() {
        return Err(format!(
            "decoded output already exists: '{}'",
            output.display()
        ));
    }
    fs::create_dir_all(output)
        .map_err(|error| format!("Could not create '{}': {error}", output.display()))?;

    let mut paths = fs::read_dir(input)
        .map_err(|error| format!("Could not read '{}': {error}", input.display()))?
        .map(|entry| entry.map(|value| value.path()))
        .collect::<Result<Vec<_>, _>>()
        .map_err(|error| format!("Could not enumerate CRAM outputs: {error}"))?;
    paths.retain(|path| {
        path.extension()
            .is_some_and(|extension| extension == "cram")
    });
    paths.sort();
    if paths.is_empty() {
        return Err(format!("No CRAM outputs found in '{}'", input.display()));
    }

    let mut writers = HashMap::<PathBuf, FastqWriter>::new();
    for path in paths {
        decode_file(&path, output, mode, &mut writers)?;
    }
    for (path, mut writer) in writers {
        writer
            .flush()
            .map_err(|error| format!("Could not flush '{}': {error}", path.display()))?;
        let encoder = writer
            .into_inner()
            .map_err(|error| format!("Could not finish '{}': {error}", path.display()))?;
        encoder
            .finish()
            .map_err(|error| format!("Could not finish '{}': {error}", path.display()))?;
    }
    Ok(())
}

fn decode_file(
    path: &Path,
    output: &Path,
    mode: Mode,
    writers: &mut HashMap<PathBuf, FastqWriter>,
) -> Result<(), String> {
    let stem = path
        .file_stem()
        .and_then(|value| value.to_str())
        .ok_or_else(|| format!("Invalid CRAM output filename: '{}'", path.display()))?;
    let mut reader = bam::Reader::from_path(path)
        .map_err(|error| format!("Could not open '{}': {error}", path.display()))?;
    for result in reader.records() {
        let record =
            result.map_err(|error| format!("Could not decode '{}': {error}", path.display()))?;
        let filename = match mode {
            Mode::Single => {
                if record.is_paired() {
                    return Err(format!(
                        "Paired record found in single output '{}'",
                        path.display()
                    ));
                }
                format!("{stem}.fastq.gz")
            }
            Mode::Paired => {
                let first = record.is_first_in_template();
                let second = record.is_last_in_template();
                if !record.is_paired() || first == second {
                    return Err(format!(
                        "Malformed paired record found in output '{}'",
                        path.display()
                    ));
                }
                format!("{stem}_R{}.fastq.gz", if first { 1 } else { 2 })
            }
        };
        let output_path = output.join(filename);
        if !writers.contains_key(&output_path) {
            let file = File::create(&output_path).map_err(|error| {
                format!("Could not create '{}': {error}", output_path.display())
            })?;
            let encoder = GzEncoder::new(file, Compression::fast());
            writers.insert(output_path.clone(), BufWriter::new(encoder));
        }
        let writer = writers
            .get_mut(&output_path)
            .ok_or("FASTQ validation writer cache failure")?;
        write_record(writer, &record, path)?;
    }
    Ok(())
}

fn write_record(writer: &mut FastqWriter, record: &bam::Record, path: &Path) -> Result<(), String> {
    let qualities = record
        .qual()
        .iter()
        .map(|quality| {
            quality
                .checked_add(33)
                .ok_or_else(|| format!("Quality overflow while decoding '{}'", path.display()))
        })
        .collect::<Result<Vec<_>, _>>()?;
    if qualities.len() != record.seq_len() {
        return Err(format!(
            "Sequence/quality length mismatch while decoding '{}'",
            path.display()
        ));
    }
    writer
        .write_all(b"@")
        .and_then(|()| writer.write_all(record.qname()))
        .and_then(|()| writer.write_all(b"\n"))
        .and_then(|()| writer.write_all(&record.seq().as_bytes()))
        .and_then(|()| writer.write_all(b"\n+\n"))
        .and_then(|()| writer.write_all(&qualities))
        .and_then(|()| writer.write_all(b"\n"))
        .map_err(|error| format!("Could not write decoded FASTQ: {error}"))
}
