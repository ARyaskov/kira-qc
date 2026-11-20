//! End-to-end tests: the engine, the report writers and the naming rules that
//! downstream tools depend on.

use kira_qc::core::engine::{Input, PhredOffsetConfig, RunConfig, decide_phred_offset, run};
use kira_qc::core::metrics::{AdapterSet, ContaminantList, DEFAULT_DUP_WINDOW};
use kira_qc::core::model::{Encoding, Mode, Status};
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

#[test]
fn output_is_identical_across_thread_counts() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("reads.fastq");
    // More than four 2 MiB chunks of packed sequence and quality.
    write_fastq(&path, 45_000, true);
    let artifacts = |threads: &str| {
        let output = dir.path().join(threads);
        let result = std::process::Command::new(env!("CARGO_BIN_EXE_kira-qc"))
            .args([
                "--quiet",
                "--json",
                "--dup-window",
                "20000",
                "-t",
                threads,
                "-o",
            ])
            .env("SOURCE_DATE_EPOCH", "0")
            .arg(&output)
            .arg(&path)
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        [
            "reads_fastqc/fastqc_data.txt",
            "reads_fastqc/summary.txt",
            "reads_fastqc/fastqc_report.html",
            "reads_fastqc/kira_qc.json",
            "reads_fastqc.html",
            "reads_fastqc.zip",
        ]
        .map(|name| std::fs::read(output.join(name)).unwrap())
    };
    let single = artifacts("1");
    for threads in ["2", "4", "16"] {
        let other = artifacts(threads);
        assert!(
            single == other,
            "reports differ between 1 and {threads} threads"
        );
    }
}

#[test]
fn fastqc_text_uses_the_expected_headers() {
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
fn partial_duplication_window_keeps_fastqc_headers_tabular() {
    let path = scratch("partial_dup_window").join("reads.fastq");
    write_fastq(&path, 11, false);
    let mut cfg = config(&path, 1, Mode::Short);
    cfg.dup_window = 10;
    let output = run(cfg).unwrap();
    assert!(!output.ctx.dup_window_complete);
    let metrics = output.agg.finalize(&output.ctx);
    let mut txt = Vec::new();
    report::fastqc_txt::write_to(&mut txt, &metrics, &output.ctx).unwrap();
    let txt = String::from_utf8(txt).unwrap();
    let start = txt.find(">>Sequence Duplication Levels\t").unwrap();
    let section = txt[start..].split(">>END_MODULE").next().unwrap();
    assert!(section.contains("#Total Deduplicated Percentage\t"));
    for header in section.lines().filter(|line| line.starts_with('#')) {
        assert!(
            header.split('\t').count() >= 2,
            "invalid column header: {header}"
        );
    }
    let json = report::json::render(&metrics, &output.ctx);
    assert!(json.contains("\"reads_analysed\": 10"));
    assert!(json.contains("\"whole_file\": false"));
}

#[test]
fn basic_gc_uses_called_bases_and_truncates_like_fastqc() {
    let path = scratch("basic_gc").join("reads.fastq");
    std::fs::write(&path, "@one\nGCAANN\n+\nIIIIII\n@two\nGNNN\n+\nIIII\n").unwrap();
    let text = render_data_txt(&path, 1, Mode::Short);
    assert!(text.contains("\n%GC\t60\n"), "{text}");
    std::fs::write(&path, "@one\nGCAAAAN\n+\nIIIIIII\n").unwrap();
    let text = render_data_txt(&path, 1, Mode::Short);
    assert!(text.contains("\n%GC\t33\n"), "{text}");
}

#[test]
fn sequence_quality_uses_the_lower_modal_score_and_inclusive_thresholds() {
    let path = scratch("sequence_quality_thresholds").join("reads.fastq");
    for (quality, status) in [("56", "fail"), ("<=", "warn")] {
        let input = format!("@first\nAC\n+\n{quality}\n@second\nAC\n+\nII\n");
        std::fs::write(&path, input).unwrap();
        let text = render_data_txt(&path, 1, Mode::Short);
        assert!(text.contains(&format!(">>Per sequence quality scores\t{status}\n")));
    }
}

#[test]
fn adapter_content_includes_ambiguous_bases_in_position_coverage() {
    let path = scratch("adapter_ambiguous").join("reads.fastq");
    std::fs::write(
        &path,
        "@one\nAGATCGGAAGAGR\n+\nIIIIIIIIIIIII\n@two\nAGATCGGAAGAGA\n+\nIIIIIIIIIIIII\n",
    )
    .unwrap();
    let mut cfg = config(&path, 1, Mode::Short);
    cfg.nogroup = true;
    let output = run(cfg).unwrap();
    let metrics = output.agg.finalize(&output.ctx);
    assert_eq!(metrics.adapter_content.last().unwrap().values[0], 100.0);
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

#[cfg(feature = "no-kmer")]
#[test]
fn requesting_kmers_without_the_module_fails() {
    let path = scratch("kmer_disabled").join("reads.fastq");
    write_fastq(&path, 10, false);
    let mut cfg = config(&path, 1, Mode::Short);
    cfg.kmer_enabled = true;
    let error = run(cfg).err().expect("unavailable analysis must fail");
    assert!(error.to_string().contains("no-kmer"));
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
fn long_mode_n_warning_includes_reads_above_the_higher_threshold() {
    let path = scratch("long_n_thresholds").join("reads.fastq");
    let mut input = String::new();
    for i in 0..100 {
        let n_count = match i {
            0..=3 => 26,
            4..=7 => 11,
            _ => 0,
        };
        input.push_str(&format!(
            "@read{i}\n{}{}\n+\n{}\n",
            "N".repeat(n_count),
            "A".repeat(105 - n_count),
            "I".repeat(105)
        ));
    }
    std::fs::write(&path, input).unwrap();
    let output = run(config(&path, 1, Mode::Long)).unwrap();
    assert_eq!(output.agg.reads_n_gt10, 8);
    assert_eq!(output.agg.reads_n_gt20, 4);
    assert_eq!(
        output.agg.finalize(&output.ctx).statuses.per_seq_n,
        Status::Warn
    );
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
    assert!(text.contains("\n10-49\t"), "expected FastQC base groups");
}

#[test]
fn grouped_quality_uses_per_cycle_statistics_and_coverage() {
    let path = scratch("grouped_quality").join("reads.fastq");
    let mut quality = vec![b'I'; 100];
    quality[9] = b'5';
    quality[10] = b'H';
    let quality = String::from_utf8(quality).unwrap();
    let seq = "ACGT".repeat(25);
    let mut input = format!("@read\n{seq}\n+\n{quality}\n").repeat(101);
    input.push_str(&format!("@tail\n{seq}A\n+\n{quality}!\n"));
    std::fs::write(&path, input).unwrap();

    let output = run(config(&path, 1, Mode::Short)).unwrap();
    let metrics = output.agg.finalize(&output.ctx);
    let grouped = &metrics.per_base_qual[9];
    assert_eq!(grouped.group.label(), "10-11");
    assert_eq!(grouped.median, 29.5);
    let tail = metrics.per_base_qual.last().unwrap();
    assert_eq!(tail.group.label(), "100-101");
    assert_eq!(tail.mean, 20.0);
    assert_eq!(tail.median, 40.0);
    assert_eq!(metrics.statuses.per_base_qual, Status::Pass);
}

#[test]
fn gzip_and_multi_block_bgzf_match_plain_fastq() {
    use flate2::{Compression, write::GzEncoder};
    use std::io::Write;

    fn gzip(bytes: &[u8]) -> Vec<u8> {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::fast());
        encoder.write_all(bytes).unwrap();
        encoder.finish().unwrap()
    }

    fn bgzf_block(bytes: &[u8]) -> Vec<u8> {
        let gzip = gzip(bytes);
        let mut block = gzip[..10].to_vec();
        block[3] |= 4; // FEXTRA
        block.extend_from_slice(&[6, 0, b'B', b'C', 2, 0]);
        let size_minus_one = u16::try_from(gzip.len() + 8 - 1).unwrap();
        block.extend_from_slice(&size_minus_one.to_le_bytes());
        block.extend_from_slice(&gzip[10..]);
        block
    }

    let dir = tempfile::tempdir().unwrap();
    let plain = dir.path().join("reads.fastq");
    let gz = dir.path().join("reads.fastq.gz");
    let concatenated = dir.path().join("concatenated.fastq.gz");
    let bgzf = dir.path().join("reads.fastq.bgzf");
    write_fastq(&plain, 1_000, true);
    let bytes = std::fs::read(&plain).unwrap();
    std::fs::write(&gz, gzip(&bytes)).unwrap();
    let members: Vec<u8> = bytes.chunks(32_768).flat_map(gzip).collect();
    std::fs::write(&concatenated, members).unwrap();
    let mut blocks: Vec<u8> = bytes.chunks(32_768).flat_map(bgzf_block).collect();
    let real_blocks_len = blocks.len();
    // The BGZF EOF marker has a canonical byte representation.
    blocks.extend([
        31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, b'B', b'C', 2, 0, 27, 0, 3, 0, 0, 0, 0, 0, 0, 0,
        0, 0,
    ]);
    std::fs::write(&bgzf, &blocks).unwrap();
    let normalized_report = |path: &Path| {
        render_data_txt(path, 2, Mode::Short)
            .lines()
            .filter(|line| !line.starts_with("Filename\t"))
            .collect::<Vec<_>>()
            .join("\n")
    };
    let expected = normalized_report(&plain);
    assert_eq!(normalized_report(&gz), expected);
    assert_eq!(normalized_report(&concatenated), expected);
    assert_eq!(normalized_report(&bgzf), expected);
    std::fs::write(&bgzf, &blocks[..real_blocks_len]).unwrap();
    let error = run(config(&bgzf, 2, Mode::Short))
        .err()
        .expect("missing BGZF EOF marker must fail");
    assert!(
        error.to_string().contains("BgzfMissingEofMarker"),
        "{error}"
    );
}

#[test]
fn empty_input_is_an_error_not_an_empty_report() {
    let dir = scratch("empty");
    let path = dir.join("empty.fastq");
    std::fs::write(&path, "").unwrap();
    assert!(run(config(&path, 1, Mode::Short)).is_err());
}

#[test]
fn automatic_phred64_applies_to_the_chunk_that_resolves_encoding() {
    let path = scratch("phred64_chunk").join("reads.fastq");
    let record = format!("@read\n{}\n+\n{}\n", "ACGT".repeat(25), "h".repeat(100));
    std::fs::write(&path, record.repeat(110_000)).unwrap();

    let automatic = run(config(&path, 2, Mode::Short)).unwrap();
    let mut fixed_config = config(&path, 2, Mode::Short);
    fixed_config.phred_offset = PhredOffsetConfig::Fixed(64);
    let fixed = run(fixed_config).unwrap();
    assert_eq!(automatic.ctx.encoding, Encoding::Illumina15);
    assert_eq!(automatic.agg.total_reads, 110_000);
    assert_eq!(automatic.agg.per_pos_qual, fixed.agg.per_pos_qual);
    assert_eq!(automatic.agg.per_seq_mean_q_hist[40], 110_000);
}

#[test]
fn malformed_input_after_a_full_chunk_fails_the_run() {
    use std::io::Write;

    let path = scratch("truncated_after_chunk").join("reads.fastq");
    write_fastq(&path, 11_000, false);
    let mut file = std::fs::OpenOptions::new()
        .append(true)
        .open(&path)
        .unwrap();
    file.write_all(b"@broken\nACGT\n+\n!\n").unwrap();
    drop(file);

    let mut cfg = config(&path, 4, Mode::Short);
    cfg.phred_offset = PhredOffsetConfig::Fixed(33);
    let error = run(cfg).err().expect("partial input must fail").to_string();
    assert!(error.contains("read 11001"), "{error}");
}

#[test]
fn noextract_preserves_unrelated_files_and_rejects_path_sample_names() {
    use std::process::Command;

    let root = tempfile::tempdir().unwrap();
    let input = root.path().join("reads.fastq");
    write_fastq(&input, 10, false);
    let existing = root.path().join("reads_fastqc");
    std::fs::create_dir(&existing).unwrap();
    std::fs::write(existing.join("summary.txt"), "previous report").unwrap();

    let run_cli = |extra: &[&str]| {
        Command::new(env!("CARGO_BIN_EXE_kira-qc"))
            .arg(&input)
            .args(["--noextract", "--quiet", "-t", "1", "-o"])
            .arg(root.path())
            .args(extra)
            .output()
            .unwrap()
    };
    let output = run_cli(&[]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        std::fs::read_to_string(existing.join("summary.txt")).unwrap(),
        "previous report"
    );
    let archive = std::fs::File::open(root.path().join("reads_fastqc.zip")).unwrap();
    let mut archive = zip::ZipArchive::new(archive).unwrap();
    assert!(archive.by_name("reads_fastqc/fastqc_data.txt").is_ok());
    let standalone = std::fs::read(root.path().join("reads_fastqc.html")).unwrap();
    let mut bundled = Vec::new();
    std::io::Read::read_to_end(
        &mut archive.by_name("reads_fastqc/fastqc_report.html").unwrap(),
        &mut bundled,
    )
    .unwrap();
    assert_eq!(standalone, bundled);
    assert_eq!(std::fs::read_dir(root.path()).unwrap().count(), 5);

    let output = run_cli(&["--sample-name", "../other"]);
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("filename component"));
}

#[test]
fn colliding_sample_outputs_fail_before_any_reports_are_written() {
    use std::process::Command;

    let root = tempfile::tempdir().unwrap();
    let first = root.path().join("lane1/reads.fastq");
    let second = root.path().join("lane2/reads.fq");
    for path in [&first, &second] {
        std::fs::create_dir_all(path.parent().unwrap()).unwrap();
        write_fastq(path, 10, false);
    }
    let outdir = root.path().join("new/qc");
    let output = Command::new(env!("CARGO_BIN_EXE_kira-qc"))
        .args([&first, &second])
        .args(["--quiet", "-t", "1", "-o"])
        .arg(&outdir)
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("same report"));
    assert!(!outdir.exists());

    let output = Command::new(env!("CARGO_BIN_EXE_kira-qc"))
        .args([&first, &second])
        .args(["--quiet", "--no-zip", "-t", "1"])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    for path in [&first, &second] {
        assert!(
            path.parent()
                .unwrap()
                .join("reads_fastqc/summary.txt")
                .is_file()
        );
    }
}

#[test]
fn named_stdin_reports_expose_distinct_multiqc_filenames() {
    use std::process::{Command, Stdio};

    let root = tempfile::tempdir().unwrap();
    let input = root.path().join("reads.fastq");
    write_fastq(&input, 10, false);
    for sample in ["first", "second"] {
        let output = Command::new(env!("CARGO_BIN_EXE_kira-qc"))
            .args([
                "-",
                "--sample-name",
                sample,
                "--quiet",
                "--no-zip",
                "-t",
                "1",
                "-o",
            ])
            .arg(root.path())
            .stdin(Stdio::from(std::fs::File::open(&input).unwrap()))
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let report = root.path().join(format!("{sample}_fastqc/fastqc_data.txt"));
        assert!(
            std::fs::read_to_string(report)
                .unwrap()
                .contains(&format!("\nFilename\t{sample}\n"))
        );
    }
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

#[test]
fn json_module_data_follows_each_ignore_directive_independently() {
    let path = scratch("json_ignore_directives").join("reads.fastq");
    write_fastq(&path, 10, false);
    let output = run(config(&path, 1, Mode::Short)).unwrap();
    let metrics = output.agg.finalize(&output.ctx);
    let mut ctx = output.ctx;
    for ignored in ["duplication", "overrepresented"] {
        ctx.ignored_modules = vec![report::resolve_limits_key(ignored).unwrap().to_string()];
        let json = report::json::render(&metrics, &ctx);
        assert_eq!(
            json.contains("\"overrepresented_sequences\": ["),
            ignored != "overrepresented"
        );
        assert_eq!(
            json.contains("\"duplication\": {"),
            ignored != "duplication"
        );
    }
    ctx.ignored_modules = vec![
        "adapter_content".to_string(),
        "per_sequence_gc_content".to_string(),
    ];
    let json = report::json::render(&metrics, &ctx);
    assert!(!json.contains("\"adapter_content\": ["));
    assert!(!json.contains("\"per_sequence_gc\": ["));
}
