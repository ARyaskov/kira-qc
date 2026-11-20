# kira-qc processing pipeline

## Overview

kira-qc performs FastQC-style quality control for high-throughput sequencing
data. It targets short-read Illumina libraries by default and has a long-read
mode for ONT and PacBio data that avoids the misleading per-base summaries
FastQC would produce there.

## Input formats

- FASTQ (plain text)
- FASTQ.GZ (gzip)
- BGZF
- stdin, selected by passing `-` as the file name

Quality scores are interpreted as Phred values. The offset can be fixed at 33 or
64 or detected automatically.

## Processing pipeline

### 1) Input detection and opening
**Input:** file path or stdin
**Output:** record stream
**Purpose:** pick a backend from the file's magic bytes.

Plain FASTQ is memory-mapped. Gzip and BGZF use streaming decompression with
bounded buffering. A stdin stream is read through a buffered reader.

### 2) Chunking and encoding detection
**Input:** record stream
**Output:** packed record chunks carrying a global read index and the Phred offset
**Purpose:** split the input into bounded chunks and settle the quality encoding.

The producer thread packs `seq||qual` pairs into 2 MiB chunks and parses the
flowcell tile out of the read header when the header format carries one. Tile
parsing stops after the first chunk if no read had one.

Short-read auto detection samples the input up to a chunk boundary, with limits
of 100,000 reads or about 16 MiB of packed reads. Long-read auto mode uses
Phred+33 immediately. The selected offset travels with every chunk, so the
input is read once and workers need no synchronisation. Short-read detection
defaults to Phred+33 and requires positive evidence
before selecting Phred+64: any byte below `;` rules out Phred+64, and any byte
above the Phred+64 ceiling rules it out as well, which keeps high-quality
NovaSeq and PacBio HiFi data from being shifted by 31.

### 3) Parallel aggregation
**Input:** record chunks
**Output:** per-worker metric aggregates
**Purpose:** compute streaming, mergeable summaries without global locks.

Short-read processing makes one fused pass over each read that fills the
per-position quality histograms and base counts, the GC and N totals and the
read's mean quality. Long-read processing uses SIMD base counting and a
per-read quality histogram, from which the mean error probability is derived.

### 4) Specialised modules
**Input:** parsed reads
**Output:** module-specific aggregates

- **Duplication and overrepresentation** share one exact sequence table, as they
  do in FastQC. Sequences longer than 75 bp are truncated to 50 bp before
  hashing. Counting covers the first `--dup-window` reads; because the window is
  a prefix defined by global read index, the result does not depend on how
  chunks were distributed.
- **Adapter Content** matches 12 bp adapter prefixes with a single Aho-Corasick
  automaton in overlapping mode, so an adapter starting inside another adapter's
  match is still reported. Long-read mode uses 20-28 bp probes and searches both
  strands.
- **Per tile sequence quality** accumulates quality sums per tile and position,
  capped at 4096 tiles.
- **K-mer Content** is opt-in and uses a rolling 2-bit encoding over ten
  positional bins.

### 5) Deterministic reduction
**Input:** per-worker aggregates
**Output:** merged aggregate
**Purpose:** combine partial results.

Every accumulator is a commutative, associative sum, so the merge order is
irrelevant and the result is identical for any thread count.

### 6) Finalisation
**Input:** merged aggregate
**Output:** rows and statuses per module

Per-position data is collapsed into FastQC's widening base groups: positions 1
to 9 individually, then intervals that grow with read length. This keeps the
report bounded for long reads, and it is the row format MultiQC expects.

Statuses use FastQC's thresholds. Adapter percentages are divided by the number
of reads that actually reach each position, so variable-length libraries are not
diluted.

### 7) Report writing
**Input:** finalised metrics
**Output:** report files

All writers walk one shared module registry, so `fastqc_data.txt`, `summary.txt`,
the HTML report and the JSON report can never disagree about which modules the
run produced. Every interpolated string is escaped for its target format.

## Output artifacts

- `fastqc_data.txt`: FastQC 0.12-format module sections
- `summary.txt`: one PASS/WARN/FAIL line per module
- `fastqc_report.html`: self-contained HTML report
- `kira_qc.json`: machine-readable report (`--json`)
- `<sample>_fastqc.zip`: bundle of the above (unless `--no-zip`)
- `latex/`: LaTeX report with SVG and PDF figures (`--export-latex`)

## Determinism and reproducibility

Output is byte-identical across thread counts, including the ZIP bundle, which
uses a fixed archive timestamp. `SOURCE_DATE_EPOCH` pins the report timestamp.
A multi-chunk test compares the complete reports and ZIP bundles across 1, 2, 4
and 16 threads. The compatibility job also compares selected metrics against
FastQC 0.12.1 and imports its ZIP and the generated ZIPs with MultiQC 1.33.

## Performance considerations

- The input is read once; encoding detection does not re-open or re-inflate it
- Memory-mapped I/O for plain FASTQ, streaming decompression otherwise
- 2 MiB chunks and a bounded channel keep peak memory proportional to threads
- One fused pass per short read instead of separate counting passes
- Worker-local aggregation with no shared mutable state
