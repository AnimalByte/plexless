use std::collections::{HashMap, HashSet};
use std::fs::{self, File, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use flate2::Compression;
use flate2::write::GzEncoder;

use crate::samples::SampleSheet;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum OutputMate {
    Single,
    R1,
    R2,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum WriterKey {
    Sample { sample_id: u32, mate: OutputMate },

    Unassigned { mate: OutputMate },
}

struct CachedWriter {
    writer: GzEncoder<BufWriter<File>>,
    last_used: u64,
}

struct CachedRawWriter {
    writer: BufWriter<File>,
    last_used: u64,
}

pub struct WriterManager {
    output_dir: PathBuf,
    sample_names: Vec<String>,
    cache: HashMap<WriterKey, CachedWriter>,
    max_open_writers: usize,
    write_unassigned: bool,
    paired: bool,
    compression_level: u32,
    tick: u64,
}

pub struct CompressedWriterManager {
    output_dir: PathBuf,
    sample_names: Vec<String>,
    cache: HashMap<WriterKey, CachedRawWriter>,
    max_open_writers: usize,
    write_unassigned: bool,
    paired: bool,
    tick: u64,
}

impl WriterManager {
    pub fn new(
        output_dir: PathBuf,
        samples: &SampleSheet,
        paired: bool,
        write_unassigned: bool,
        max_open_writers: usize,
        compression_level: u32,
    ) -> Result<Self, String> {
        if max_open_writers == 0 {
            return Err("Maximum open writers must be greater than 0".into());
        }

        if compression_level > 9 {
            return Err("Compression level must be between 0 and 9".into());
        }

        prepare_output_dir(&output_dir)?;
        let sample_names = build_sample_names(samples)?;

        Ok(Self {
            output_dir,
            sample_names,
            cache: HashMap::new(),
            max_open_writers,
            write_unassigned,
            paired,
            compression_level,
            tick: 0,
        })
    }

    pub fn write_sample(
        &mut self,
        sample_id: u32,
        mate: OutputMate,
        id: &[u8],
        seq: &[u8],
        qual: &[u8],
    ) -> Result<(), String> {
        validate_mate(self.paired, mate)?;
        validate_sample_id(&self.sample_names, sample_id)?;

        let key = WriterKey::Sample { sample_id, mate };

        self.write_record(key, id, seq, qual)
    }

    pub fn write_unassigned(
        &mut self,
        mate: OutputMate,
        id: &[u8],
        seq: &[u8],
        qual: &[u8],
    ) -> Result<(), String> {
        if !self.write_unassigned {
            return Ok(());
        }

        validate_mate(self.paired, mate)?;

        let key = WriterKey::Unassigned { mate };

        self.write_record(key, id, seq, qual)
    }

    pub fn finish(mut self) -> Result<(), String> {
        for (_, cached) in self.cache.drain() {
            finish_writer(cached.writer)?;
        }

        Ok(())
    }

    fn write_record(
        &mut self,
        key: WriterKey,
        id: &[u8],
        seq: &[u8],
        qual: &[u8],
    ) -> Result<(), String> {
        if seq.len() != qual.len() {
            return Err("Sequence and quality lengths do not match".into());
        }

        self.tick = self.tick.wrapping_add(1);
        let tick = self.tick;

        if !self.cache.contains_key(&key) {
            self.open_writer(key)?;
        }

        let cached = self.cache.get_mut(&key).ok_or("Writer cache failure")?;
        cached.last_used = tick;

        let id = id.strip_prefix(b"@").unwrap_or(id);

        cached
            .writer
            .write_all(b"@")
            .map_err(|e| format!("Could not write FASTQ: {e}"))?;
        cached
            .writer
            .write_all(id)
            .map_err(|e| format!("Could not write FASTQ: {e}"))?;
        cached
            .writer
            .write_all(b"\n")
            .map_err(|e| format!("Could not write FASTQ: {e}"))?;
        cached
            .writer
            .write_all(seq)
            .map_err(|e| format!("Could not write FASTQ: {e}"))?;
        cached
            .writer
            .write_all(b"\n+\n")
            .map_err(|e| format!("Could not write FASTQ: {e}"))?;
        cached
            .writer
            .write_all(qual)
            .map_err(|e| format!("Could not write FASTQ: {e}"))?;
        cached
            .writer
            .write_all(b"\n")
            .map_err(|e| format!("Could not write FASTQ: {e}"))?;

        Ok(())
    }

    fn open_writer(&mut self, key: WriterKey) -> Result<(), String> {
        if self.cache.len() >= self.max_open_writers {
            self.evict_oldest()?;
        }

        let path = output_path(&self.output_dir, &self.sample_names, key)?;

        let file = OpenOptions::new()
            .create(true)
            .append(true)
            .open(&path)
            .map_err(|e| format!("Could not open output '{}': {e}", path.display()))?;

        let buffered = BufWriter::new(file);
        let writer = GzEncoder::new(buffered, Compression::new(self.compression_level));

        self.cache.insert(
            key,
            CachedWriter {
                writer,
                last_used: self.tick,
            },
        );

        Ok(())
    }

    fn evict_oldest(&mut self) -> Result<(), String> {
        let oldest_key = self
            .cache
            .iter()
            .min_by_key(|(_, cached)| cached.last_used)
            .map(|(&key, _)| key)
            .ok_or("Writer cache is empty")?;

        let cached = self
            .cache
            .remove(&oldest_key)
            .ok_or("Writer cache failure")?;

        finish_writer(cached.writer)
    }
}

impl CompressedWriterManager {
    pub fn new(
        output_dir: PathBuf,
        samples: &SampleSheet,
        paired: bool,
        write_unassigned: bool,
        max_open_writers: usize,
        compression_level: u32,
    ) -> Result<Self, String> {
        if max_open_writers == 0 {
            return Err("Maximum open writers must be greater than 0".into());
        }

        if compression_level > 9 {
            return Err("Compression level must be between 0 and 9".into());
        }

        prepare_output_dir(&output_dir)?;
        let sample_names = build_sample_names(samples)?;

        Ok(Self {
            output_dir,
            sample_names,
            cache: HashMap::new(),
            max_open_writers,
            write_unassigned,
            paired,
            tick: 0,
        })
    }

    pub fn append_sample_member(
        &mut self,
        sample_id: u32,
        mate: OutputMate,
        member: &[u8],
    ) -> Result<(), String> {
        validate_mate(self.paired, mate)?;
        validate_sample_id(&self.sample_names, sample_id)?;

        let key = WriterKey::Sample { sample_id, mate };
        self.append_member(key, member)
    }

    pub fn append_unassigned_member(
        &mut self,
        mate: OutputMate,
        member: &[u8],
    ) -> Result<(), String> {
        if !self.write_unassigned {
            return Ok(());
        }

        validate_mate(self.paired, mate)?;

        let key = WriterKey::Unassigned { mate };
        self.append_member(key, member)
    }

    pub fn finish(mut self) -> Result<(), String> {
        for (_, cached) in self.cache.drain() {
            finish_raw_writer(cached.writer)?;
        }

        Ok(())
    }

    fn append_member(&mut self, key: WriterKey, member: &[u8]) -> Result<(), String> {
        if member.is_empty() {
            return Ok(());
        }

        self.tick = self.tick.wrapping_add(1);
        let tick = self.tick;

        if !self.cache.contains_key(&key) {
            self.open_writer(key)?;
        }

        let cached = self.cache.get_mut(&key).ok_or("Writer cache failure")?;
        cached.last_used = tick;
        cached
            .writer
            .write_all(member)
            .map_err(|e| format!("Could not append compressed FASTQ member: {e}"))
    }

    fn open_writer(&mut self, key: WriterKey) -> Result<(), String> {
        if self.cache.len() >= self.max_open_writers {
            self.evict_oldest()?;
        }

        let path = output_path(&self.output_dir, &self.sample_names, key)?;

        let file = OpenOptions::new()
            .create(true)
            .append(true)
            .open(&path)
            .map_err(|e| format!("Could not open output '{}': {e}", path.display()))?;

        self.cache.insert(
            key,
            CachedRawWriter {
                writer: BufWriter::new(file),
                last_used: self.tick,
            },
        );

        Ok(())
    }

    fn evict_oldest(&mut self) -> Result<(), String> {
        let oldest_key = self
            .cache
            .iter()
            .min_by_key(|(_, cached)| cached.last_used)
            .map(|(&key, _)| key)
            .ok_or("Writer cache is empty")?;

        let cached = self
            .cache
            .remove(&oldest_key)
            .ok_or("Writer cache failure")?;

        finish_raw_writer(cached.writer)
    }
}

fn build_sample_names(samples: &SampleSheet) -> Result<Vec<String>, String> {
    let mut sample_names = Vec::with_capacity(samples.samples.len());
    let mut seen_names = HashSet::new();

    for sample in &samples.samples {
        let name = sanitize_name(&sample.name);

        if name.is_empty() {
            return Err(format!(
                "Sample '{}' cannot be converted to a valid filename",
                sample.name
            ));
        }

        if !seen_names.insert(name.clone()) {
            return Err(format!(
                "Sample filename collision after sanitizing '{}'",
                sample.name
            ));
        }

        sample_names.push(name);
    }

    Ok(sample_names)
}

fn validate_sample_id(sample_names: &[String], sample_id: u32) -> Result<(), String> {
    let sample_index = usize::try_from(sample_id).map_err(|_| "Invalid sample ID")?;

    if sample_index >= sample_names.len() {
        return Err(format!("Unknown sample ID {sample_id}"));
    }

    Ok(())
}

fn validate_mate(paired: bool, mate: OutputMate) -> Result<(), String> {
    match (paired, mate) {
        (false, OutputMate::Single) | (true, OutputMate::R1 | OutputMate::R2) => Ok(()),
        (false, _) => Err("Single-end output requires OutputMate::Single".into()),
        (true, OutputMate::Single) => Err("Paired-end output requires R1 or R2".into()),
    }
}

fn output_path(
    output_dir: &Path,
    sample_names: &[String],
    key: WriterKey,
) -> Result<PathBuf, String> {
    let (name, mate) = match key {
        WriterKey::Sample { sample_id, mate } => {
            let index = usize::try_from(sample_id).map_err(|_| "Invalid sample ID")?;
            let name = sample_names.get(index).ok_or("Unknown sample ID")?;
            (name.as_str(), mate)
        }
        WriterKey::Unassigned { mate } => ("unassigned", mate),
    };

    let filename = match mate {
        OutputMate::Single => format!("{name}.fastq.gz"),
        OutputMate::R1 => format!("{name}_R1.fastq.gz"),
        OutputMate::R2 => format!("{name}_R2.fastq.gz"),
    };

    Ok(output_dir.join(filename))
}

fn finish_writer(writer: GzEncoder<BufWriter<File>>) -> Result<(), String> {
    let mut buffered = writer
        .finish()
        .map_err(|e| format!("Could not finish gzip output: {e}"))?;

    buffered
        .flush()
        .map_err(|e| format!("Could not flush output: {e}"))
}

fn finish_raw_writer(mut writer: BufWriter<File>) -> Result<(), String> {
    writer
        .flush()
        .map_err(|e| format!("Could not flush output: {e}"))
}

fn prepare_output_dir(path: &Path) -> Result<(), String> {
    if path.exists() {
        if !path.is_dir() {
            return Err(format!(
                "Output path '{}' exists but is not a directory",
                path.display()
            ));
        }

        if fs::read_dir(path)
            .map_err(|e| format!("Could not read output directory '{}': {e}", path.display()))?
            .next()
            .is_some()
        {
            return Err(format!(
                "Output directory '{}' is not empty",
                path.display()
            ));
        }
    } else {
        fs::create_dir_all(path).map_err(|e| {
            format!(
                "Could not create output directory '{}': {e}",
                path.display()
            )
        })?;
    }

    Ok(())
}

fn sanitize_name(name: &str) -> String {
    name.chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || matches!(c, '-' | '_' | '.') {
                c
            } else {
                '_'
            }
        })
        .collect()
}
