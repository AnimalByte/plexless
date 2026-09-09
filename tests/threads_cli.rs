use clap::Parser;
use plexless::cli::{ByteSizeSetting, Cli, Command, OutputMode};

fn parse_demux(extra: &[&str]) -> Cli {
    let mut arguments = vec![
        "plexless",
        "demux",
        "--reads",
        "reads.fastq",
        "--structure",
        "R1_4A",
        "--barcodes",
        "barcodes.tsv",
        "--samples",
        "samples.tsv",
        "--output",
        "out",
    ];
    arguments.extend_from_slice(extra);
    Cli::try_parse_from(arguments).unwrap()
}

#[test]
fn threads_is_a_global_demux_option() {
    let cli = Cli::try_parse_from([
        "simplex",
        "demux",
        "--threads",
        "4",
        "--reads",
        "reads.fastq",
        "--structure",
        "R1_4A",
        "--barcodes",
        "barcodes.tsv",
        "--samples",
        "samples.tsv",
        "--output",
        "out",
    ])
    .expect("--threads should parse after the demux subcommand");

    assert_eq!(cli.threads, 4);
}

#[test]
fn advanced_output_and_qc_options_parse_human_sizes() {
    let cli = Cli::try_parse_from([
        "plexless",
        "demux",
        "--reads",
        "reads.fastq",
        "--structure",
        "R1_4A",
        "--barcodes",
        "barcodes.tsv",
        "--samples",
        "samples.tsv",
        "--output",
        "out",
        "--output-chunk-size",
        "256K",
        "--output-buffer-memory",
        "2G",
        "--max-open-files",
        "128",
        "--low-sample-fraction",
        "0.02",
    ])
    .unwrap();
    let Command::Demux(args) = cli.command;
    assert_eq!(args.output_chunk_size, ByteSizeSetting::Bytes(256 * 1024));
    assert_eq!(
        args.output_buffer_memory,
        ByteSizeSetting::Bytes(2 * 1024 * 1024 * 1024)
    );
    assert_eq!(args.max_open_files, Some(128));
    assert_eq!(args.low_sample_fraction, 0.02);
}

#[test]
fn zero_threads_is_rejected() {
    let result = Cli::try_parse_from([
        "simplex",
        "demux",
        "--threads",
        "0",
        "--reads",
        "reads.fastq",
        "--structure",
        "R1_4A",
        "--barcodes",
        "barcodes.tsv",
        "--samples",
        "samples.tsv",
        "--output",
        "out",
    ]);

    assert!(result.is_err());
}

#[test]
fn output_mode_defaults_to_auto() {
    let Command::Demux(args) = parse_demux(&[]).command;
    assert_eq!(args.output_mode, OutputMode::Auto);
}

#[test]
fn all_output_modes_parse_and_invalid_values_fail() {
    for (value, expected) in [
        ("auto", OutputMode::Auto),
        ("direct", OutputMode::Direct),
        ("buffered", OutputMode::Buffered),
    ] {
        let Command::Demux(args) = parse_demux(&["--output-mode", value]).command;
        assert_eq!(args.output_mode, expected);
    }

    let result = Cli::try_parse_from([
        "plexless",
        "demux",
        "--reads",
        "reads.fastq",
        "--structure",
        "R1_4A",
        "--barcodes",
        "barcodes.tsv",
        "--samples",
        "samples.tsv",
        "--output",
        "out",
        "--output-mode",
        "unknown",
    ]);
    assert!(result.is_err());
}
