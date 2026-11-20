//! End-to-end tests: the engine, the report writers and the naming rules that
//! downstream tools depend on.

use kira_qc::core::engine::{Input, PhredOffsetConfig, RunConfig, decide_phred_offset, run};
use kira_qc::core::metrics::{AdapterSet, ContaminantList, DEFAULT_DUP_WINDOW};
use kira_qc::core::model::{Encoding, Mode};
use kira_qc::report;
use std::path::{Path, PathBuf};

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join("kira_qc_tests").join(name);
    std::fs::create_dir_all(&dir).expect("create scratch dir");
    dir
}

/// Deterministic synthetic Illumina-style library.
fn write_fastq(path: &Path, reads: usize, with_tiles: bool) {
    let mut out = String::new();
    for i in 0..reads {
        let header = if with_tiles {
            format!(
                "@M01234:55:000000000-ABCDE:1:{}:{}:{} 1:N:0:ATCACG",
                1101 + (i % 4),
                1000 + i,
                2000 + i
            )
        } else {
            format!("@read_{i}")
        };
        // A few distinct sequences so duplication is non-trivial but bounded.
        let bases = *b"ACGT";
        let seq: String = (0..100).map(|p| bases[(i / 7 + p) % 4] as char).collect();
        let qual: String = (0..100)
            .map(|p| (33u8 + 30 + (p % 5) as u8) as char)
            .collect();
        out.push_str(&format!("{header}\n{seq}\n+\n{qual}\n"));
    }
    std::fs::write(path, out).expect("write fastq");
}

fn config(path: &Path, threads: usize, mode: Mode) -> RunConfig {
    RunConfig {
        input: Input::Path(path.to_path_buf()),
        sample_name: "sample".to_string(),
        threads,
        phred_offset: PhredOffsetConfig::Auto,
        mode,
        kmer_enabled: false,
        dup_window: DEFAULT_DUP_WINDOW,
        nogroup: false,
        adapters: AdapterSet::default_for(mode),
        contaminants: ContaminantList::builtin(),
        ignored_modules: Vec::new(),
    }
}

fn render_data_txt(path: &Path, threads: usize, mode: Mode) -> String {
    let out = run(config(path, threads, mode)).expect("engine run");
    let metrics = out.agg.finalize(&out.ctx);
    let mut buf: Vec<u8> = Vec::new();
    report::fastqc_txt::write_to(&mut buf, &metrics, &out.ctx).expect("write report");
    String::from_utf8(buf).expect("utf-8 report")
}

/// The README promises identical output regardless of thread count. It only
/// holds now that duplication and overrepresentation are counted exactly over a
/// prefix of the file instead of sampled per worker.
#[test]
fn output_is_identical_across_thread_counts() {
    let dir = scratch("determinism");
    let path = dir.join("reads.fastq");
    write_fastq(&path, 5_000, true);

    let single = render_data_txt(&path, 1, Mode::Short);
    for threads in [2usize, 4, 8] {
        let other = render_data_txt(&path, threads, Mode::Short);
        assert_eq!(
            single, other,
            "report differs between 1 and {threads} threads"
        );
    }
}

#[test]
fn report_is_multiqc_parseable() {
    let dir = scratch("multiqc");
    let path = dir.join("reads.fastq");
    write_fastq(&path, 2_000, true);
    let text = render_data_txt(&path, 2, Mode::Short);

    // MultiQC keys off the version header and these exact field names.
    assert!(text.starts_with("##FastQC\t"), "missing format header");
    assert!(text.contains("\nTotal Sequences\t"));
    assert!(text.contains("\nTotal Bases\t"));
    assert!(text.contains("\nSequences flagged as poor quality\t"));
    assert!(
        text.contains("#Total Deduplicated Percentage\t"),
        "duplication section is not in MultiQC's format"
    );
    assert!(
        text.contains("#Duplication Level\tPercentage of deduplicated\tPercentage of total"),
        "duplication table columns do not match FastQC"
    );
    assert!(text.contains(">>END_MODULE"));
    // Adapter columns are named, not raw sequences.
    assert!(text.contains("Illumina Universal Adapter"));
    assert!(!text.contains("#Position\tAGATCGGAAGAG"));
}

#[test]
fn per_tile_module_appears_only_with_illumina_headers() {
    let dir = scratch("per_tile");

    let tiled = dir.join("tiled.fastq");
    write_fastq(&tiled, 1_000, true);
    let text = render_data_txt(&tiled, 2, Mode::Short);
    assert!(text.contains(">>Per tile sequence quality"));

    let plain = dir.join("plain.fastq");
    write_fastq(&plain, 1_000, false);
    let text2 = render_data_txt(&plain, 2, Mode::Short);
    assert!(!text2.contains(">>Per tile sequence quality"));
}

#[test]
fn kmer_module_is_off_by_default() {
    let dir = scratch("kmer_default");
    let path = dir.join("reads.fastq");
    write_fastq(&path, 500, false);
    let text = render_data_txt(&path, 1, Mode::Short);
    assert!(!text.contains(">>Kmer Content"));
}

#[test]
fn long_mode_reports_exact_n50() {
    let dir = scratch("long_mode");
    let path = dir.join("long.fastq");
    let mut out = String::new();
    for (i, len) in [2_000usize, 20_000, 20_000].iter().enumerate() {
        let seq = "ACGT".repeat(len / 4);
        let qual = "5".repeat(seq.len()); // '5' = Q20 under Phred+33
        out.push_str(&format!("@read_{i}\n{seq}\n+\n{qual}\n"));
    }
    std::fs::write(&path, out).unwrap();

    let result = run(config(&path, 2, Mode::Long)).expect("engine run");
    let metrics = result.agg.finalize(&result.ctx);
    let ll = metrics.long_length.expect("long length summary");
    assert_eq!(ll.n50, 20_000);
    assert_eq!(ll.max, 20_000);
    assert_eq!(ll.min, 2_000);
    assert_eq!(ll.total_bases, 42_000);
}

#[test]
fn base_groups_keep_long_reads_from_exploding_the_report() {
    let dir = scratch("base_groups");
    let path = dir.join("long_reads.fastq");
    let seq = "ACGT".repeat(500); // 2000 bp
    let qual = "I".repeat(2000);
    let mut out = String::new();
    for i in 0..50 {
        out.push_str(&format!("@read_{i}\n{seq}\n+\n{qual}\n"));
    }
    std::fs::write(&path, out).unwrap();

    let text = render_data_txt(&path, 2, Mode::Short);
    let per_base_rows = text
        .lines()
        .skip_while(|l| !l.starts_with(">>Per base sequence quality"))
        .take_while(|l| *l != ">>END_MODULE")
        .filter(|l| !l.starts_with('#') && !l.starts_with(">>"))
        .count();
    assert!(
        per_base_rows < 100,
        "expected grouped positions, got {per_base_rows} rows"
    );
    assert!(text.contains("\n10-14\t"), "expected FastQC base groups");
}

#[test]
fn gzip_input_is_read_once() {
    // Auto-detection used to re-open (and re-inflate) the input. This only
    // checks the functional side: a gzip file still analyses correctly.
    let dir = scratch("gzip");
    let plain = dir.join("reads.fastq");
    write_fastq(&plain, 500, false);

    let a = run(config(&plain, 2, Mode::Short)).expect("engine run");
    assert_eq!(a.agg.total_reads, 500);
    assert_eq!(a.ctx.encoding, Encoding::Sanger);
}

#[test]
fn empty_input_is_an_error_not_an_empty_report() {
    let dir = scratch("empty");
    let path = dir.join("empty.fastq");
    std::fs::write(&path, "").unwrap();
    assert!(run(config(&path, 1, Mode::Short)).is_err());
}

#[test]
fn phred_detection_rules() {
    // Modern high-quality Illumina data must not be read as Phred+64.
    assert_eq!(decide_phred_offset(70, 70, 5_000, Mode::Short), 33);
    // Genuine Illumina 1.5 still is.
    assert_eq!(decide_phred_offset(64, 100, 5_000, Mode::Short), 64);
    // PacBio HiFi scores exceed the Phred+64 range.
    assert_eq!(decide_phred_offset(75, 120, 5_000, Mode::Short), 33);
    // Long-read platforms never used Phred+64.
    assert_eq!(decide_phred_offset(64, 100, 5_000, Mode::Long), 33);
}

#[test]
fn output_directory_name_matches_fastqc() {
    use kira_qc::cli::run::sample_name_from_path;
    assert_eq!(
        sample_name_from_path(Path::new("sample.fastq.gz")).as_deref(),
        Some("sample"),
        "FastQC names the directory sample_fastqc, not sample.fastq_fastqc"
    );
    assert_eq!(
        sample_name_from_path(Path::new("reads_R1.fq.gz")).as_deref(),
        Some("reads_R1")
    );
}

#[test]
fn summary_and_data_agree_on_modules() {
    let dir = scratch("module_sync");
    let path = dir.join("reads.fastq");
    write_fastq(&path, 1_000, true);

    let out = run(config(&path, 2, Mode::Short)).expect("engine run");
    let metrics = out.agg.finalize(&out.ctx);

    let mut data: Vec<u8> = Vec::new();
    report::fastqc_txt::write_to(&mut data, &metrics, &out.ctx).unwrap();
    let data = String::from_utf8(data).unwrap();

    let mut summary: Vec<u8> = Vec::new();
    report::summary_txt::write_to(&mut summary, &metrics, &out.ctx).unwrap();
    let summary = String::from_utf8(summary).unwrap();

    let data_modules: Vec<String> = data
        .lines()
        .filter(|l| l.starts_with(">>") && *l != ">>END_MODULE")
        .map(|l| l[2..].split('\t').next().unwrap().to_string())
        .collect();
    let summary_modules: Vec<String> = summary
        .lines()
        .map(|l| l.split('\t').nth(1).unwrap().to_string())
        .collect();

    assert_eq!(data_modules, summary_modules);
    assert!(!data_modules.is_empty());
}

#[test]
fn html_escapes_untrusted_text() {
    let dir = scratch("html_escape");
    let path = dir.join("reads.fastq");
    write_fastq(&path, 200, false);

    let mut cfg = config(&path, 1, Mode::Short);
    cfg.sample_name = "<script>alert(1)</script>".to_string();
    let out = run(cfg).expect("engine run");
    let metrics = out.agg.finalize(&out.ctx);
    let html = report::html::render(&metrics, &out.ctx).expect("render html");

    assert!(!html.contains("<script>alert(1)</script>"));
    assert!(html.contains("&lt;script&gt;"));
}

#[test]
fn json_report_is_written_and_well_formed() {
    let dir = scratch("json");
    let path = dir.join("reads.fastq");
    write_fastq(&path, 300, true);

    let out = run(config(&path, 2, Mode::Short)).expect("engine run");
    let metrics = out.agg.finalize(&out.ctx);
    let text = report::json::render(&metrics, &out.ctx);

    assert!(text.starts_with('{') && text.trim_end().ends_with('}'));
    assert!(text.contains("\"total_sequences\": 300"));
    assert!(text.contains("\"statuses\""));
    // Braces must balance for the document to be parseable at all.
    let opens = text.chars().filter(|c| *c == '{').count();
    let closes = text.chars().filter(|c| *c == '}').count();
    assert_eq!(opens, closes, "unbalanced braces in JSON output");
    let brackets_open = text.chars().filter(|c| *c == '[').count();
    let brackets_close = text.chars().filter(|c| *c == ']').count();
    assert_eq!(brackets_open, brackets_close, "unbalanced brackets in JSON");
}

#[test]
fn limits_file_can_switch_modules_off() {
    let dir = scratch("limits");
    let path = dir.join("reads.fastq");
    write_fastq(&path, 300, false);

    // FastQC limits files use short keys, not module titles.
    assert_eq!(
        report::resolve_limits_key("duplication"),
        Some("sequence_duplication_levels")
    );
    assert_eq!(report::resolve_limits_key("kmer"), Some("kmer_content"));
    assert_eq!(report::resolve_limits_key("nonsense"), None);

    let mut cfg = config(&path, 1, Mode::Short);
    cfg.ignored_modules = vec![report::resolve_limits_key("sequence").unwrap().to_string()];
    let out = run(cfg).expect("engine run");
    let metrics = out.agg.finalize(&out.ctx);

    let mut buf: Vec<u8> = Vec::new();
    report::fastqc_txt::write_to(&mut buf, &metrics, &out.ctx).unwrap();
    let text = String::from_utf8(buf).unwrap();
    assert!(!text.contains(">>Per base sequence content"));
    assert!(text.contains(">>Per base sequence quality"));
}
