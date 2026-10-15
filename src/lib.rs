pub mod barcodes;
pub mod cli;
pub(crate) mod cram;
pub mod decoder;
pub mod demux;
pub mod encoder;
pub mod fastq;
pub(crate) mod input;
pub(crate) mod output;
pub(crate) mod parallel;
pub(crate) mod parallel_input;
pub mod pass_through;
pub(crate) mod qc;
pub mod routing;
pub mod samples;
pub mod stats;
pub mod structure;
pub(crate) mod thread_plan;
pub mod writer;

pub(crate) mod preprocessing;
pub(crate) mod technical;
pub mod trim;

pub(crate) mod overlap;

pub mod fastp_options;
pub(crate) mod quality_trim;

pub(crate) mod fastp_runtime;

pub(crate) mod read_ops;

pub(crate) mod adapter_detection;

pub(crate) mod distribution;
pub(crate) mod report;

pub(crate) mod overrepresentation;
