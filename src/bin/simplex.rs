//! Backward-compatible executable name for workflows migrating from simplex.

use clap::Parser;
use plexless::cli::{Cli, Command};

fn main() {
    if let Err(error) = run() {
        eprintln!("Error: {error}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), String> {
    let cli = Cli::parse();
    let threads = cli.threads;

    match cli.command {
        Command::Demux(args) => {
            args.validate()?;
            plexless::demux::run_with_threads(args, threads)
        }
    }
}
