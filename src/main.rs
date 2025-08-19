use clap::{Parser, ValueEnum};
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use rayon::prelude::*;
use rand::random;
use regex::Regex;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::process::{Command, Stdio};
use std::sync::Arc;

#[derive(ValueEnum, Clone)]
enum InputType {
    Bed,
    Bam,
}

#[derive(Parser)]
#[clap(name = "bedfragment_ds", version = "6.3")]
struct Args {
    /// Input mode: 'bed' (default) or 'bam'
    #[clap(long, value_enum, default_value_t = InputType::Bed)]
    input_type: InputType,

    /// Chromosome sizes file (required if input_type=bed)
    #[clap(long, required_if_eq("input_type", "bed"))]
    chrom_sizes: Option<PathBuf>,

    /// Optional blacklist BED file (only used if input_type=bam)
    #[clap(long)]
    blacklist: Option<PathBuf>,

    /// Fragment BED or BAM files to process
    files: Vec<PathBuf>,

    /// Z-score threshold for excluding low-yield libraries (default 1.5)
    #[clap(short, long, default_value = "1.5")]
    exclude_sd: f64,

    /// Keep intermediate bedGraph files (only in bed mode)
    #[clap(long)]
    keep_bedgraph: bool,

    /// Whether to keep temporary downsampled BAM files (only for BAM input)
    #[clap(long)]
    keep_tmp_bam: bool,

    /// Number of threads (0 = use all available cores)
    #[clap(short = 't', long, default_value = "0")]
    threads: usize,
}

fn mean(values: &[usize]) -> f64 {
    values.iter().map(|&v| v as f64).sum::<f64>() / values.len() as f64
}

fn std_dev(values: &[usize], mean: f64) -> f64 {
    let var = values.iter().map(|&v| {
        let diff = v as f64 - mean;
        diff * diff
    }).sum::<f64>() / values.len() as f64;
    var.sqrt()
}

fn parse_chrom_order(chrom_sizes: &PathBuf) -> Result<HashMap<String, usize>, Box<dyn Error>> {
    let file = File::open(chrom_sizes)?;
    let reader = BufReader::new(file);
    let mut map = HashMap::new();
    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let chrom = line.split_whitespace().next().unwrap().to_string();
        map.insert(chrom, i);
    }
    Ok(map)
}

fn count_fragments(path: &PathBuf) -> Result<usize, Box<dyn Error>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut count = 0usize;
    for (i, line) in reader.lines().enumerate() {
        if i == 0 {
            continue; // skip header
        }
        if !line?.trim().is_empty() {
            count += 1;
        }
    }
    Ok(count)
}

fn reservoir_sample(
    path: &PathBuf,
    min_count: usize,
) -> Result<(String, Vec<String>), Box<dyn Error>> {
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);
    let mut header = String::new();
    reader.read_line(&mut header)?;
    header = header.trim_end().to_string();

    let mut sample: Vec<String> = Vec::with_capacity(min_count);
    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        if i < min_count {
            sample.push(line);
        } else {
            let j = random::<usize>() % (i + 1);
            if j < min_count {
                sample[j] = line;
            }
        }
    }
    Ok((header, sample))
}

fn create_50bp_bins(chrom_sizes: &PathBuf) -> Result<PathBuf, Box<dyn Error>> {
    let bins_path = PathBuf::from("genome_50bp_bins.bed");
    if bins_path.exists() && bins_path.metadata()?.len() > 0 {
        return Ok(bins_path);
    }
    let status = Command::new("bedtools")
        .args(["makewindows", "-g"])
        .arg(chrom_sizes)
        .args(["-w", "50"])
        .stdout(File::create(&bins_path)?)
        .status()?;
    if !status.success() {
        return Err("bedtools makewindows failed".into());
    }
    if bins_path.metadata()?.len() == 0 {
        return Err("bedtools makewindows produced empty bins file".into());
    }
    Ok(bins_path)
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();

    let nthreads = if args.threads > 0 {
        args.threads
    } else {
        num_cpus::get()
    };

    rayon::ThreadPoolBuilder::new()
        .num_threads(nthreads)
        .build_global()
        .unwrap();

    if args.files.is_empty() {
        eprintln!("No fragment files provided.");
        std::process::exit(1);
    }

    match args.input_type {
        InputType::Bed => {
            // BED pipeline (unchanged)
            let chrom_sizes = args.chrom_sizes.as_ref().unwrap();
            let chrom_order = Arc::new(parse_chrom_order(&chrom_sizes)?);
            let bins_bed = Arc::new(create_50bp_bins(&chrom_sizes)?);

            let mut frag_counts = Vec::new();
            for f in &args.files {
                let c = count_fragments(f)?;
                frag_counts.push((f.clone(), c));
            }
            let counts_only: Vec<_> = frag_counts.iter().map(|(_, c)| *c).collect();
            let mean_val = mean(&counts_only);
            let sd_val = std_dev(&counts_only, mean_val);
            let cutoff = (mean_val - args.exclude_sd * sd_val).max(0.0);
            eprintln!("QC: Mean={}, SD={}, cutoff={}", mean_val, sd_val, cutoff);
            let filtered = frag_counts
                .iter()
                .filter(|(_, c)| (*c as f64) >= cutoff)
                .cloned()
                .collect::<Vec<_>>();
            if filtered.is_empty() {
                eprintln!("No samples pass the QC cutoff");
                std::process::exit(1);
            }
            let excluded = frag_counts
                .iter()
                .filter(|(_, c)| (*c as f64) < cutoff)
                .cloned()
                .collect::<Vec<_>>();
            if !excluded.is_empty() {
                eprintln!("Excluded samples with low fragment counts:");
                for (f, c) in &excluded {
                    eprintln!("  {} => {}", f.display(), c);
                }
            }
            let min_frag_count = filtered.iter().map(|(_, c)| *c).min().unwrap();

            let m = Arc::new(MultiProgress::new());

            filtered.par_iter().for_each(|(file_path, _)| {
                let pb = m.add(ProgressBar::new(6));
                pb.set_style(
                    ProgressStyle::default_bar()
                        .template("{msg} {bar:40.cyan/blue} {pos}/{len} ({eta})")
                        .expect("Progress bar template error")
                        .progress_chars("#>-"),
                );
                let filename = file_path.file_name().unwrap().to_string_lossy().to_string();
                let msg = format!("Processing {}", filename);
                pb.set_message(msg.clone());

                if let Ok((header, mut sample)) = reservoir_sample(file_path, min_frag_count)
                {
                    pb.inc(1);

                    let order_map = chrom_order.clone();
                    sample.retain(|line| {
                        let chrom = line.split('\t').next().unwrap();
                        order_map.contains_key(chrom)
                    });

                    sample.sort_by(|a, b| {
                        let a_parts: Vec<&str> = a.split('\t').collect();
                        let b_parts: Vec<&str> = b.split('\t').collect();
                        let a_rank = *order_map.get(a_parts[0]).unwrap_or(&usize::MAX);
                        let b_rank = *order_map.get(b_parts[0]).unwrap_or(&usize::MAX);
                        if a_rank == b_rank {
                            let a_start = a_parts[1].parse::<u32>().unwrap_or(0);
                            let b_start = b_parts[1].parse::<u32>().unwrap_or(0);
                            return a_start.cmp(&b_start);
                        }
                        a_rank.cmp(&b_rank)
                    });
                    pb.inc(1);

                    let out_bed = file_path.with_file_name(format!(
                        "{}_downsampled.bed",
                        file_path.file_stem().unwrap().to_string_lossy()
                    ));
                    {
                        let out_file = File::create(&out_bed).unwrap();
                        let mut writer = BufWriter::new(out_file);
                        writeln!(writer, "{}", header).unwrap();
                        for line in &sample {
                            writeln!(writer, "{}", line).unwrap();
                        }
                    }
                    pb.inc(1);

                    let sorted_bed = out_bed.with_file_name(format!(
                        "{}_sorted.bed",
                        out_bed.file_stem().unwrap().to_string_lossy()
                    ));
                    let bedtools_sort_status = Command::new("bedtools")
                        .args(["sort", "-faidx"])
                        .arg(&*chrom_sizes)
                        .args(["-i"])
                        .arg(&out_bed)
                        .stdout(File::create(&sorted_bed).unwrap())
                        .status()
                        .expect("bedtools sort failed");
                    if !bedtools_sort_status.success() {
                        eprintln!("bedtools sort failed for {}", out_bed.display());
                        let msg = format!("Sort failed for {}", filename);
                        pb.finish_with_message(msg);
                        return;
                    }
                    pb.inc(1);

                    let coverage_bed =
                        file_path.with_file_name(format!("{}_50bp_counts.bed", filename));
                    let coverage_status = Command::new("bedtools")
                        .args(["coverage", "-a"])
                        .arg(&*bins_bed)
                        .args(["-b"])
                        .arg(&sorted_bed)
                        .args(["-counts"])
                        .stdout(File::create(&coverage_bed).unwrap())
                        .status()
                        .expect("bedtools coverage failed");
                    if !coverage_status.success() {
                        eprintln!("bedtools coverage failed for {}", sorted_bed.display());
                        let msg = format!("Coverage failed for {}", filename);
                        pb.finish_with_message(msg);
                        return;
                    }
                    pb.inc(1);

                    let bedgraph =
                        file_path.with_file_name(format!("{}_50bp.bedGraph", filename));
                    let awk_status = Command::new("awk")
                        .arg(r#"OFS="\t" {print $1, $2, $3, $4}"#)
                        .stdin(File::open(&coverage_bed).unwrap())
                        .stdout(File::create(&bedgraph).unwrap())
                        .status()
                        .expect("awk command failed");
                    if !awk_status.success() {
                        eprintln!("awk conversion failed for {}", coverage_bed.display());
                        let msg = format!("awk failed for {}", filename);
                        pb.finish_with_message(msg);
                        return;
                    }

                    let sorted_bedgraph =
                        file_path.with_file_name(format!("{}_50bp_sorted.bedGraph", filename));
                    let sort_status = Command::new("sort")
                        .args(["--parallel=1", "-k1,1", "-k2,2n"])
                        .arg(&bedgraph)
                        .stdout(File::create(&sorted_bedgraph).unwrap())
                        .status()
                        .expect("sort failed");
                    if !sort_status.success() {
                        eprintln!("Sorting bedGraph failed for {}", bedgraph.display());
                        let msg = format!("bedGraph sort failed {}", filename);
                        pb.finish_with_message(msg);
                        return;
                    }

                    let bigwig = file_path.with_file_name(format!("{}_50bp.bw", filename));
                    let bw_status = Command::new("bedGraphToBigWig")
                        .arg(&sorted_bedgraph)
                        .arg(&*chrom_sizes)
                        .arg(&bigwig)
                        .status()
                        .expect("bedGraphToBigWig failed");
                    if bw_status.success() {
                        eprintln!("Wrote {}", bigwig.display());
                        let msg = format!("Completed {}", filename);
                        pb.finish_with_message(msg);
                    } else {
                        eprintln!("bedGraphToBigWig failed for {}", sorted_bedgraph.display());
                        let msg = format!("BigWig failed {}", filename);
                        pb.finish_with_message(msg);
                    }

                    if !args.keep_bedgraph {
                        let _ = std::fs::remove_file(&coverage_bed);
                        let _ = std::fs::remove_file(&bedgraph);
                        let _ = std::fs::remove_file(&sorted_bedgraph);
                        let _ = std::fs::remove_file(&sorted_bed);
                        let _ = std::fs::remove_file(&out_bed);
                    }
                } else {
                    let msg = format!("Sampling failed {}", file_path.display().to_string());
                    pb.finish_with_message(msg);
                }
            });
        }
        InputType::Bam => {
            if args.files.is_empty() {
                eprintln!("No BAM files provided");
                std::process::exit(1);
            }
            let min_count = {
                let mut counts = Vec::new();
                for f in &args.files {
                    let count_output = Command::new("samtools")
                        .args(&["view", "-c", "-f", "2", "-F", "260", f.to_str().unwrap()])
                        .output()
                        .expect("failed to run samtools count");
                    if !count_output.status.success() {
                        eprintln!("samtools count failed for {}", f.display());
                        std::process::exit(1);
                    }
                    let count_str = String::from_utf8_lossy(&count_output.stdout);
                    let sample_count: usize = count_str.trim().parse().unwrap_or(0);
                    counts.push((f.clone(), sample_count));
                }
                let counts_only: Vec<_> = counts.iter().map(|(_, c)| *c).collect();
                let mean_val = mean(&counts_only);
                let sd_val = std_dev(&counts_only, mean_val);
                let cutoff = (mean_val - args.exclude_sd * sd_val).max(0.0);
                eprintln!("QC: Mean={}, SD={}, cutoff={}", mean_val, sd_val, cutoff);
                let filtered = counts
                    .iter()
                    .filter(|(_, c)| (*c as f64) >= cutoff)
                    .cloned()
                    .collect::<Vec<_>>();
                if filtered.is_empty() {
                    eprintln!("No BAM samples pass the QC cutoff");
                    std::process::exit(1);
                }
                let excluded = counts
                    .iter()
                    .filter(|(_, c)| (*c as f64) < cutoff)
                    .cloned()
                    .collect::<Vec<_>>();
                if !excluded.is_empty() {
                    eprintln!("Excluded BAM samples with low fragment counts:");
                    for (f, c) in &excluded {
                        eprintln!("  {} => {}", f.display(), c);
                    }
                }
                filtered.iter().map(|(_, c)| *c).min().unwrap()
            };

            let m = Arc::new(MultiProgress::new());

            args.files.par_iter().for_each(|file_path| {
                let file_str = file_path.to_str().unwrap();
                let count_output = Command::new("samtools")
                    .args(&["view", "-c", "-f", "2", "-F", "260", file_str])
                    .output()
                    .expect("failed to run samtools count");
                let count_str = String::from_utf8_lossy(&count_output.stdout);
                let sample_count = count_str.trim().parse::<f64>().unwrap_or(0.0);
                if sample_count < min_count as f64 {
                    return;
                }

                let pb = m.add(ProgressBar::new(1));
                pb.set_style(
                    ProgressStyle::default_spinner()
                        .template("{msg} {spinner} {elapsed_precise}")
                        .expect("Progress bar template error"),
                );

                let filename = file_path.file_name().unwrap().to_string_lossy().to_string();
                let msg = format!("Processing BAM {}", filename);
                pb.set_message(msg.clone());

                let fraction = (min_count as f64 / sample_count).min(1.0);
                let seed_fraction = format!("42.{:03}", (fraction * 1000.0) as u32);

                let tmp_bam = file_path.with_file_name(format!("{}_downsampled.bam", filename));
                // Write downsampled BAM to disk
                let samtools_status = Command::new("samtools")
                    .args(&[
                        "view",
                        "-b",
                        "-s",
                        &seed_fraction,
                        "-f",
                        "2",
                        "-F",
                        "260",
                        file_str,
                    ])
                    .stdout(File::create(&tmp_bam).unwrap())
                    .status()
                    .expect("samtools downsampling failed");
                if !samtools_status.success() {
                    eprintln!("samtools downsampling failed for {}", filename);
                    pb.finish_with_message(format!("Failed {}", filename));
                    return;
                }

                // Index the downsampled BAM file
                let samtools_index_status = Command::new("samtools")
                    .args(&["index", tmp_bam.to_str().unwrap()])
                    .status()
                    .expect("samtools index failed for downsampled BAM");
                if !samtools_index_status.success() {
                    eprintln!("samtools index failed for {}", filename);
                    pb.finish_with_message(format!("Failed {}", filename));
                    return;
                }

                let bamcov_out = file_path.with_file_name(format!("{}_50bp.bw", filename));

                let mut bamcov_cmd = Command::new("bamCoverage");
                bamcov_cmd.args(&[
                    "-p", "1",
                    "-b", tmp_bam.to_str().unwrap(),
                    "--binSize", "50",
                    "--normalizeUsing", "None",
                    "-o", bamcov_out.to_str().unwrap(),
                ]);
                if let Some(blacklist_path) = &args.blacklist {
                    bamcov_cmd.args(&["--blackListFileName", blacklist_path.to_str().unwrap()]);
                }

                let bamcov_status = bamcov_cmd.status().unwrap_or_else(|e| {
                    eprintln!("Failed bamCoverage for {}: {}", filename, e);
                    std::process::exit(1);
                });
                if bamcov_status.success() {
                    eprintln!("Wrote {}", bamcov_out.display());
                    pb.finish_with_message(format!("Completed {}", filename));
                } else {
                    eprintln!("bamCoverage failed for {}", filename);
                    pb.finish_with_message(format!("Failed {}", filename));
                }

                if !args.keep_tmp_bam {
                    let _ = std::fs::remove_file(&tmp_bam);
                    let bai_path = tmp_bam.with_extension("bam.bai");
                    let _ = std::fs::remove_file(&bai_path);
                    let bai_path2 = tmp_bam.with_extension("bai");
                    let _ = std::fs::remove_file(&bai_path2);
                }
            });
        }
    }

    Ok(())
}
