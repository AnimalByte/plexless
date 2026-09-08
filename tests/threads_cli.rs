use clap::Parser;
use simplex::cli::Cli;

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
