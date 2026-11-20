#![cfg_attr(feature = "no-kmer", allow(dead_code))]

pub mod cli;
pub mod core;
pub mod report;
pub mod simd;

/// Version string embedded in reports. The `##FastQC` header line carries the
/// FastQC format version so downstream parsers (MultiQC) recognise the file.
pub const KIRA_QC_VERSION: &str = env!("CARGO_PKG_VERSION");

/// FastQC data-format version implemented by the text report.
pub const FASTQC_FORMAT_VERSION: &str = "0.12.1";
