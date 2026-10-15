//! Trim-only entry point. Reuses the bounded, ordered FASTQ processing pipeline.
use crate::cli::TrimArgs;

pub fn run_with_threads(args: TrimArgs, threads: usize) -> Result<(), String> {
    let mut args = args.into_pipeline_args();
    crate::fastp_options::prepare(&mut args)?;
    let input = crate::input::resolve_input_source(&args)?;
    crate::parallel::run_trim(args, threads, input)
}

pub(crate) fn ensure_outputs(
    dir: &std::path::Path,
    paired: bool,
    level: u32,
    options: &crate::fastp_options::FastpOptions,
) -> Result<(), String> {
    let layout = crate::output::OutputLayout::new(1, paired, false).with_fastp(options);
    for index in 0..layout.stream_count()? {
        let key = layout.key_at(index)?;
        if matches!(
            key.target,
            crate::output::OutputTarget::Auxiliary {
                stream: crate::output::AuxiliaryStream::Stdout
            }
        ) {
            continue;
        }
        let path = crate::writer::configured_path(dir, &["trimmed".into()], key, options, 0)?;
        if !path.exists() {
            let file = std::fs::OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&path)
                .map_err(|e| format!("Could not create '{}': {e}", path.display()))?;
            if layout.gzip(key) {
                flate2::write::GzEncoder::new(file, flate2::Compression::new(level))
                    .finish()
                    .map_err(|e| e.to_string())?;
            }
        }
    }
    Ok(())
}

pub(crate) fn print_summary(counts: &crate::demux::DemuxCounts) {
    eprintln!("Trimming complete");
    eprintln!("  fragments processed: {}", counts.preprocessing.processed);
    eprintln!("  fragments retained: {}", counts.preprocessing.retained);
    eprintln!("  fragments rejected: {}", counts.preprocessing.rejected);
    for mate in 0..2 {
        eprintln!(
            "  R{} full adapter reads: {}",
            mate + 1,
            counts.preprocessing.full[mate]
        );
        eprintln!(
            "  R{} partial adapter reads: {}",
            mate + 1,
            counts.preprocessing.partial[mate]
        );
    }
    if counts.orphan_r1 + counts.orphan_r2 > 0 {
        eprintln!("  orphan R1: {}", counts.orphan_r1);
        eprintln!("  orphan R2: {}", counts.orphan_r2);
    }
}
