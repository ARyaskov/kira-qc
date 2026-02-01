# kira-qc

## Project overview

kira-qc is a FastQC-compatible QC tool written in Rust. It generates FastQC-style
text and HTML/LaTeX reports from FASTQ/FASTQ.GZ input with deterministic, bounded-memory
processing and SIMD-accelerated hot paths. It is designed for short-read Illumina
QC with a dedicated long-read mode for ONT/PacBio data.

## Key features

- FastQC-compatible outputs: `fastqc_data.txt`, `summary.txt`, HTML, and ZIP bundle
- FASTQ and FASTQ.GZ support
- Streaming, bounded memory processing
- SIMD acceleration (AVX2 on x86_64, NEON on aarch64)
- Deterministic multithreaded execution
- Short-read (Illumina) and long-read (ONT/PacBio) modes

## Supported QC modules

| Module | Short-read mode | Long-read mode |
|---|---|---|
| Basic Statistics | Yes | Yes |
| Per base sequence quality | Yes | No |
| Per sequence quality scores | Yes | Yes |
| Per base sequence content | Yes | No |
| Per sequence GC content | Yes | Yes |
| Per base N content | Yes | No |
| Sequence Length Distribution | Yes | Yes (log-binned + N50/N90) |
| Sequence Duplication Levels | Yes | No |
| Overrepresented Sequences | Yes | No |
| Adapter Content | Yes | Yes (summary only) |
| K-mer Content (k=7) | Yes | No |
| Per sequence N content | No | Yes |

## Installation

Install from crates.io (Rust 1.91+ / Windows / Linux / macOS):

```bash
cargo install kira-qc
```

Or

From source (Rust 1.91+):

```
# build
cargo build --release

# run
./target/release/kira-qc --help
```

Notes:
- On x86_64, AVX2 is required at runtime.
- On aarch64 (e.g., Apple Silicon), NEON kernels are used.
- Optional build feature to disable K-mer Content (compile-time):
  - `cargo build --release --features no-kmer`

## Usage

Basic FASTQ:

```
kira-qc run reads.fastq --out qc/
```

FASTQ.GZ:

```
kira-qc run reads.fastq.gz --out qc/
```

Short-read (default):

```
kira-qc run reads.fastq.gz --out qc/ --mode short
```

Long-read:

```
kira-qc run reads.fastq.gz --out qc/ --mode long
```

LaTeX export (supplement):

```
kira-qc run reads.fastq.gz --out qc/ --export-latex supplement
```

## CLI options

| Option | Description | Default |
|---|---|---|
| `run` | Run QC on a single FASTQ/FASTQ.GZ file | Required |
| `--out <DIR>` | Output directory | Required |
| `--threads <N>` | Number of worker threads | Logical CPU count |
| `--sample-name <NAME>` | Sample name (used in output folder/ZIP) | Input file stem |
| `--phred-offset auto\|33\|64` | Quality encoding detection or fixed offset | `auto` |
| `--mode short\|long` | QC mode: short-read or long-read | `short` |
| `--no-zip` | Disable ZIP bundle creation | Off (ZIP enabled) |
| `--export-latex summary\|supplement` | Generate LaTeX export | Disabled |

## Output description

Output directory structure:

```
<out>/<sample_name>_fastqc/
  fastqc_data.txt
  summary.txt
  fastqc_report.html
  latex/ (if --export-latex)
<out>/<sample_name>_fastqc.zip
```

- `fastqc_data.txt`: FastQC-style module sections and tabular data.
- `summary.txt`: One-line PASS/WARN/FAIL status per module.
- `fastqc_report.html`: Self-contained HTML report (no external assets).
- `{sample_name}_fastqc.zip`: ZIP bundle of the output directory (unless `--no-zip`).
- `latex/`: LaTeX export (optional), including SVG figures and `kira_qc.tex`.

## LaTeX export for publications

kira-qc can generate an article-ready LaTeX report alongside the HTML output.
Use `--export-latex summary` for a concise figure set or `--export-latex supplement`
for a full module-by-module report. PDF figures are generated directly.

LaTeX artifacts are written under:

```
<out>/<sample_name>_fastqc/latex/
```

For Overleaf upload (Upload Project Overleaf feature), a ready-to-use ZIP is created:

```
<out>/<sample_name>_fastqc/latex/kira_qc_latex.zip
```

## Determinism and reproducibility

kira-qc produces identical results regardless of thread count or scheduling.
All aggregation is deterministic and merged by chunk index. This ensures that
pipelines can rely on stable outputs for regression tests and reproducible QC.

## Performance notes

- Streaming design with bounded memory for large inputs
- SIMD acceleration for base counting and quality processing
- Adapter matching includes a SIMD prefix prefilter before full pattern matching

## Long-read mode notes

FastQC's per-base modules assume uniform read length and are not appropriate for
ONT/PacBio data. Long-read mode disables those plots and instead focuses on:
- Log-binned length distribution with N50/N90
- Per-read quality scores
- Per-read GC content
- Per-read N content
- Adapter presence summary

## FastQC compatibility notes

- Drop-in replacement for FastQC in short-read pipelines
- Output formats and module names follow FastQC conventions
- Some heuristics (e.g., duplication, k-mer content) are approximate by design

## Limitations

- FASTQ/FASTQ.GZ input only (no BAM/CRAM)
- No GPU acceleration
- Long-read mode is not FastQC-compatible by design

## License

MIT,
See `LICENSE`.
