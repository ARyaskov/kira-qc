# kira-qc processing pipeline

## Overview

kira-qc performs FastQC-style quality control for high-throughput sequencing data. It targets short-read Illumina libraries by default and provides a dedicated long-read mode for ONT/PacBio data to avoid misleading per-base summaries.

## Input formats

- FASTQ (plain text)
- FASTQ.GZ (gzip-compressed FASTQ)

Quality scores are interpreted as Phred values. The Phred offset can be fixed (33 or 64) or auto-detected from an initial subset of reads.

## Processing pipeline

### 1) Input detection and opening
**Input:** file path  
**Output:** byte stream source  
**Purpose:** select the appropriate reader based on file extension or magic bytes.  

Plain FASTQ is memory-mapped for zero-copy access. FASTQ.GZ uses streaming decompression with bounded buffering.

### 2) Phred encoding detection
**Input:** initial read subset  
**Output:** selected Phred offset  
**Purpose:** determine whether quality scores are Phred+33 or Phred+64 (unless fixed by CLI).  

### 3) Chunking and parallel parsing
**Input:** byte stream  
**Output:** FASTQ record chunks  
**Purpose:** split input into bounded chunks that align with FASTQ 4-line records and distribute work across threads.  

The producer thread emits ordered chunks to a bounded channel. Worker threads parse reads and accumulate per-chunk metrics without global locks.

### 4) Core QC aggregation
**Input:** parsed reads  
**Output:** per-chunk metric aggregates  
**Purpose:** compute streaming, mergeable summaries of base composition, quality distributions, GC, N content, and length statistics.  

SIMD-accelerated kernels are used on x86_64 (AVX2) and aarch64 (NEON) where applicable.

### 5) Specialized modules
**Input:** parsed reads  
**Output:** module-specific aggregates  
**Purpose:** compute heavier QC modules with bounded memory:

- Sequence Duplication Levels (streaming heavy-hitter estimates)
- Overrepresented Sequences (top-K with source heuristics)
- Adapter Content (multi-pattern matching with SIMD prefilter)
- K-mer Content (short-read only; CMS + heavy-hitter tracking)

### 6) Mode-specific behavior
**Input:** aggregated metrics + mode  
**Output:** final module set  
**Purpose:** enable or disable modules based on sequencing technology:

- Short-read mode: full FastQC-style module set
- Long-read mode: length-focused and per-read summaries; per-base modules disabled

### 7) Deterministic reduction and finalization
**Input:** per-chunk aggregates  
**Output:** final metrics and module statuses  
**Purpose:** merge partial results in chunk index order to ensure deterministic outputs across thread counts.

## Output artifacts

- `fastqc_data.txt`: FastQC-style module sections and tabular data
- `summary.txt`: one-line PASS/WARN/FAIL per module
- `fastqc_report.html`: self-contained HTML report
- `{sample_name}_fastqc.zip`: ZIP bundle (unless `--no-zip`)
- `latex/` (optional): LaTeX report with SVG/PDF figures and tables (`--export-latex`)

## Determinism and reproducibility

All aggregation is mergeable and reduced in a fixed order by chunk index. 
This ensures byte-identical outputs regardless of thread scheduling, which is critical for pipeline reproducibility and regression testing.

## Performance considerations

- Memory-mapped I/O for plain FASTQ
- Streaming decompression for FASTQ.GZ with bounded buffers
- Bounded channels and worker-local aggregation to avoid global contention
- SIMD acceleration for base counting and quality processing where available
