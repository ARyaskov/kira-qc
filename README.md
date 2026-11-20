# kira-qc

## Project overview

kira-qc is a FastQC-compatible QC tool written in Rust. It reads FASTQ,
FASTQ.GZ and BGZF input (or a stream on stdin) and writes FastQC-style text,
HTML and LaTeX reports with deterministic, bounded-memory processing and
SIMD-accelerated hot paths. It targets short-read Illumina QC and has a
dedicated long-read mode for ONT and PacBio data.

## Key features

- FastQC-compatible outputs: `fastqc_data.txt`, `summary.txt`, HTML and a ZIP bundle
- MultiQC-parseable duplication and basic-statistics sections
- FASTQ, FASTQ.GZ, BGZF and stdin input
- Streaming, bounded memory processing
- SIMD acceleration (AVX2 on x86_64, NEON on aarch64, scalar fallback elsewhere)
- Byte-identical output regardless of thread count
- Short-read (Illumina) and long-read (ONT/PacBio) modes
- Optional machine-readable JSON report

## Supported QC modules

| Module | Short-read mode | Long-read mode |
|---|---|---|
| Basic Statistics | Yes | Yes |
| Per base sequence quality | Yes | No |
| Per tile sequence quality | Yes, when headers carry tile numbers | No |
| Per sequence quality scores | Yes | Yes |
| Per base sequence content | Yes | No |
| Per sequence GC content | Yes | Yes |
| Per base N content | Yes | No |
| Sequence Length Distribution | Yes | Yes (log-binned, exact N50/N90) |
| Sequence Duplication Levels | Yes | No |
| Overrepresented Sequences | Yes | No |
| Adapter Content | Yes | Yes (presence, both strands) |
| K-mer Content (k=7) | Opt-in via `--kmers` | No |
| Per sequence N content | No | Yes |

## Installation

Install from crates.io (Windows / Linux / macOS):

```bash
cargo install kira-qc
```

From source:

```bash
cargo build --release
```

Build options:

- `--no-default-features` drops the LaTeX/PDF export and cuts the dependency
  tree from 107 crates to 46.
- `--features no-kmer` strips the K-mer module at compile time. It is already
  off at runtime unless `--kmers` is given.

The default build needs no C toolchain: the ZIP writer uses only the pure-Rust
deflate backend.

The minimum supported Rust version for kira-qc is 1.95.

## Usage

kira-qc accepts the same invocation shape as `fastqc`:

```bash
kira-qc reads.fastq.gz -o qc/
```

The `run` subcommand from earlier versions still works:

```bash
kira-qc run reads.fastq --out qc/
```

Several files at once:

```bash
kira-qc sample1.fastq.gz sample2.fastq.gz -o qc/ -t 8
```

From a pipe, which is how to handle formats kira-qc does not read natively:

```bash
samtools fastq input.bam | kira-qc - -o qc/ --sample-name input
```

Long reads:

```bash
kira-qc reads.fastq.gz -o qc/ --mode long
```

## CLI options

| Option | Description | Default |
|---|---|---|
| `<FASTQ>...` | One or more input files, or `-` for stdin | Required |
| `-o`, `--outdir <DIR>` | Output directory (alias `--out`) | Directory of the input file |
| `-t`, `--threads <N>` | Worker threads | Logical CPU count |
| `--sample-name <NAME>` | Sample name; single input only | Input file stem |
| `--phred-offset auto\|33\|64` | Quality encoding | `auto` |
| `--mode short\|long` | Sequencing technology profile | `short` |
| `--nogroup` | One row per base instead of base groups | Off |
| `-k`, `--kmers` | Enable the K-mer Content module | Off |
| `--dup-window <N>` | Reads used for duplication estimates, 0 = all | 1000000 |
| `--no-zip` | Do not create the ZIP bundle | Off |
| `--noextract` | Keep ZIP and standalone HTML without extracting the directory | Off |
| `--extract` | Keep the directory (the default) | On |
| `-a`, `--adapters <FILE>` | Custom adapter list | Built-in list |
| `-c`, `--contaminants <FILE>` | Custom contaminant list | Built-in list |
| `-l`, `--limits <FILE>` | FastQC limits file; `ignore` directives only | None |
| `-f`, `--format fastq\|bam\|sam` | Input format | `fastq` |
| `--json` | Also write `kira_qc.json` | Off |
| `--export-latex summary\|supplement` | LaTeX export | Off |
| `-q`, `--quiet` | Suppress warnings | Off |

Adapter and contaminant files use FastQC's format: `name<TAB>sequence`, with
`#` starting a comment.
Custom adapter probes accept A, C, G, T and literal N, case-insensitively.

## Output description

```
<out>/<sample>_fastqc/
  fastqc_data.txt
  summary.txt
  fastqc_report.html
  kira_qc.json      (with --json)
  latex/            (with --export-latex)
<out>/<sample>_fastqc.zip
<out>/<sample>_fastqc.html
```

The sample name strips both the compression suffix and the FASTQ suffix, so
`sample.fastq.gz` produces `sample_fastqc`, matching FastQC and the wrappers
that expect that layout.

## FastQC compatibility

- `fastqc_data.txt` starts with a `##FastQC 0.12.1` header line and uses
  FastQC's module titles, base groups and status vocabulary.
- Basic Statistics includes `Total Bases` and `Sequences flagged as poor
  quality`, both of which MultiQC reads.
- Sequence Duplication Levels emits `#Total Deduplicated Percentage` and the
  `Percentage of deduplicated` / `Percentage of total` columns over FastQC's 16
  duplication levels, which is what MultiQC's duplication plot and the "% Dups"
  General Statistics column need.
- Adapter Content columns are adapter names, matching FastQC's 12 bp probe list
  (Illumina Universal, Small RNA 3' and 5', Nextera, PolyA, PolyG).
- Module statuses use FastQC's thresholds: lower quartile and median for per
  base quality, the distribution mode for per sequence quality, |A-T| and |G-C|
  imbalance for per base content, a fitted normal for per sequence GC, and 0.1 %
  / 1 % for overrepresented sequences.

Known differences:

- Duplication and overrepresentation are estimated from the first
  `--dup-window` reads (one million by default) and counted exactly within that
  window, rather than from the first 100 000 distinct sequences of the whole
  file. The HTML and JSON reports state the window when it covers only part
  of the file. A run fails if the window contains more than two million
  distinct tracked sequences. Use a smaller `--dup-window` in that case.
- `--limits` honours `ignore` directives; threshold overrides are reported as
  unsupported rather than silently dropped.
- BAM, CRAM and SAM are not read natively. Pipe them through `samtools fastq`.
- Long-read mode is deliberately not FastQC-compatible.
- Optional K-mer Content flags positional enrichment when the observed-to-expected
  ratio reaches 3 for a warning or 5 for a failure and the adjusted probability
  bound is at most 0.05 under the random-placement model.

## Determinism and reproducibility

kira-qc produces byte-identical output regardless of thread count. Every
aggregate is order-independent, and the duplication window is a prefix of the
file defined by global read index rather than by whichever worker saw a read.
Set `SOURCE_DATE_EPOCH` to pin the report timestamp; the ZIP bundle already uses
a fixed archive timestamp. A multi-chunk CLI test compares the complete reports
and ZIP bundles across 1, 2, 4 and 16 threads.

## Performance notes

- The input is read once. Short-read auto detection buffers an initial sample
  of up to about 16 MiB of packed reads. Long-read auto uses Phred+33 immediately.
- Short-read processing makes a single fused pass over sequence and quality that
  fills the per-position accumulators, base composition and mean quality.
- Chunks are 2 MiB and the producer-worker channel is bounded, so peak memory
  stays proportional to the thread count rather than growing with it.
- Adapter matching uses one Aho-Corasick automaton with overlapping search.

## Long-read mode notes

FastQC's per-base modules assume uniform read length and are not appropriate for
ONT or PacBio data. Long-read mode disables them and instead reports:

- Log-binned length distribution with exact N50, N90, median and mean
- Per-read quality, averaged over error probabilities rather than Phred scores,
  which is the convention used by NanoPlot and related tools
- Per-read GC content
- Per-read N content
- Presence of ONT ligation, ONT rapid and PacBio SMRTbell adapters, searched on
  both strands

## Limitations

- FASTQ, FASTQ.GZ, BGZF and stdin only; no native BAM/CRAM
- Compression other than gzip and BGZF must be piped in
- No GPU acceleration
- Paired-end files are analysed independently, not as pairs

## License

MIT. See `LICENSE`.
