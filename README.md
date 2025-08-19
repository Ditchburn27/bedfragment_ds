# bedfragment_ds

Tool for QC-filtered, randomized downsampling and BigWig tracks from BED or BAM files

---

## Overview

`bedfragment_ds` is a Rust command-line tool for high-throughput downsampling and track generation from chromatin fragment files. It supports both:

- **BED files** (e.g., deduplicated, paired-end cut/fragment files)
- **BAM files** (paired-end, deduplicated, mapped, coordinate-sorted)

The tool excludes low-yield libraries based on Z-score, randomly downsamples all samples to the lowest retained fragment count, and generates signal tracks (BigWig) in 50bp bins for direct visualization and comparison.

---

## Features

- Supports **parallel processing** of multiple files for high throughput
- Works with both BED and BAM (paired-end) inputs
- Automatic or user-defined thread count
- QC filtering to exclude samples with anomalously low fragment counts
- Downsamples to equal depth for all retained samples
- **Optionally** supports blacklist input for BAM mode
- Progress bars for real-time workflow monitoring

---

## Installation

### Prerequisites

- Rust toolchain (install via [rustup.rs](https://rustup.rs/))
- External tools **must be in your $PATH** for relevant steps:
  - [`bedtools`](https://bedtools.readthedocs.io/)
  - [`samtools`](http://www.htslib.org/doc/samtools.html)
  - [`bamCoverage`](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html)
  - [`bedGraphToBigWig`](https://genome.ucsc.edu/goldenPath/help/bigWig.html) (if using BED mode)
  - `awk`, `sort`

### Build
```bash
git clone https://github.com/YOURUSERNAME/fragment_downsample.git
cd fragment_downsample
cargo build --release
```

This will produce `./target/release/fragment_downsample`

---

## Usage

### 1. BED Input Mode

Requires a chromosome sizes file and one or more BED files:


```bash
./fragment_downsample
--input-type bed
--chrom-sizes mm10.chrom.sizes
sample1.bed sample2.bed sample3.bed
--threads 8
```


- **--chrom-sizes**: tab-separated file of `chrom\tlength` per line (UCSC chrom.sizes format)
- Output: One BigWig per sample, downsampled and binned to 50bp

---

### 2. BAM Input Mode (Paired-end)

Works directly from deduped, sorted, paired-end BAM files:

```bash
./fragment_downsample
--input-type bam
sample1.bam sample2.bam
--blacklist mm10-blacklist.bed
--threads 8
```

- **--blacklist** (optional): BED file of regions to exclude in bamCoverage
- Output: One BigWig per sample, from downsampled properly paired fragments

---

### Options (common to both modes)

- `--exclude-sd <float>`: Z-score threshold to exclude low-yield samples (default 1.5)
- `--threads <int>`: Number of parallel threads (default: all CPU cores)
- `--keep-bedgraph`: Keep intermediate .bedGraph files (BED mode only)
- `--keep-tmp-bam`: Keep downsampled intermediate BAMs (BAM mode only)

---

## Example outputs

- `sample1_50bp.bw`, `sample2_50bp.bw`, ... (per-sample BigWig tracks)
- [Optionally] Downsampled intermediates (`*_downsampled.bam` or `*_downsampled.bed`)

---

## Troubleshooting

- Ensure all dependencies (bedtools, samtools, bamCoverage, bedGraphToBigWig) are in your `$PATH`.
- Your BAM files **must be paired-end, indexed, sorted, and deduplicated** for best results.
- For any problems, run with more threads disabled (`--threads 1`) to check serial behavior.
- Check intermediate files and logs for filtering, downsampling, and track generation steps.