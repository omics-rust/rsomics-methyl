use std::fs::File;
use std::io::{self, BufReader, Write};
use std::path::{Path, PathBuf};

use clap::{Args, Parser, Subcommand, ValueEnum};
use noodles::core::Region;
use rsomics_common::{
    OutputArgs, Result, RsomicsError, ToolMeta, reject_output_alias, write_output,
};
use rsomics_methyl::extract::ExtractOptions;
use rsomics_methyl::mbias::MbiasOptions;
use rsomics_methyl::merge_context::{MergeContextStats, merge_context};
use rsomics_methyl::per_read::{PerReadMetric, PerReadOptions, per_read};
use rsomics_methyl::{BisulfiteStrand, ReadBounds, TrimmingOptions};
use serde::Serialize;

use crate::extract_output::extract_to_outputs;
use crate::mbias_output::mbias_to_outputs;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Debug, Parser)]
#[command(
    name = "rsomics-methyl",
    version,
    about = "Extract and inspect bisulfite-sequencing methylation evidence",
    arg_required_else_help = true
)]
pub struct Cli {
    #[command(subcommand)]
    command: Command,

    #[command(flatten)]
    pub output: OutputArgs,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Extract per-cytosine methylation counts from a sorted alignment.
    Extract(ExtractArgs),

    /// Merge strand-specific CpG and CHG cytosine metrics.
    MergeContext(MergeContextArgs),

    /// Report methylation bias by read position and bisulfite strand.
    Mbias(MbiasArgs),

    /// Report CpG methylation evidence for each alignment record.
    PerRead(PerReadArgs),
}

#[derive(Debug, Args)]
struct ExtractArgs {
    /// Indexed reference FASTA.
    reference: PathBuf,

    /// Coordinate-sorted indexed BAM or CRAM.
    input: PathBuf,

    /// Reference region using 1-based inclusive coordinates.
    #[arg(short = 'r', long)]
    region: Option<Region>,

    /// Prefix for transactional extraction outputs.
    #[arg(short, long)]
    output_prefix: PathBuf,

    #[command(flatten)]
    filters: PileupFilterArgs,

    /// Output representation.
    #[arg(long, value_enum, default_value = "standard")]
    format: ExtractFormat,

    /// Write one exhaustive Bismark-style cytosine report.
    #[arg(long, visible_alias = "cytosine_report", conflicts_with = "format")]
    cytosine_report: bool,

    /// Combine complementary CpG and CHG strand calls.
    #[arg(long)]
    merge_context: bool,

    /// Minimum methylated plus unmethylated depth; exhaustive reports include zero depth.
    #[arg(short = 'd', long, default_value_t = 1)]
    minimum_depth: u64,

    #[command(flatten)]
    variants: VariantFilterArgs,
}

#[derive(Debug, Args)]
struct VariantFilterArgs {
    /// Minimum usable opposite-strand depth before testing a site.
    #[arg(long, visible_alias = "minOppositeDepth", default_value_t = 0)]
    minimum_opposite_depth: u64,

    /// Largest allowed non-reference fraction on the opposite strand.
    #[arg(long, visible_alias = "maxVariantFrac", default_value_t = 0.0)]
    maximum_variant_fraction: f64,
}

#[derive(Debug, Args)]
struct PileupFilterArgs {
    #[command(flatten)]
    bed: BedArgs,

    /// Minimum alignment mapping quality.
    #[arg(short = 'q', long, default_value_t = 10)]
    minimum_mapping_quality: u8,

    /// Minimum base quality; must be positive.
    #[arg(short = 'p', long, default_value_t = 5)]
    minimum_base_quality: u8,

    /// Minimum converted fraction among informative non-CpG cytosines.
    #[arg(long, visible_alias = "minConversionEfficiency", default_value_t = 0.0)]
    minimum_conversion_efficiency: f64,

    /// SAM flag bits that exclude a record when any are set.
    #[arg(short = 'F', long, default_value_t = 0x0f00)]
    ignore_flags: u16,

    /// SAM flag bits that must all be set.
    #[arg(short = 'R', long, default_value_t = 0)]
    require_flags: u16,

    /// Include duplicate-marked records.
    #[arg(long)]
    keep_duplicates: bool,

    /// Include paired records whose mate is unmapped.
    #[arg(long)]
    keep_singletons: bool,

    /// Include paired records without the proper-pair flag.
    #[arg(long)]
    keep_discordant: bool,

    /// Do not exclude records with NH greater than one.
    #[arg(long)]
    ignore_nh: bool,

    /// Disable CpG metrics.
    #[arg(long)]
    no_cpg: bool,

    /// Include CHG metrics.
    #[arg(long)]
    chg: bool,

    /// Include CHH metrics.
    #[arg(long)]
    chh: bool,

    #[command(flatten)]
    trimming: TrimmingArgs,
}

#[derive(Debug, Args)]
struct BedArgs {
    /// BED intervals to include; plain text and gzip are accepted.
    #[arg(short = 'l', long, value_name = "FILE")]
    bed: Option<PathBuf>,

    /// Use BED column 6 to select top (+) or bottom (-) evidence.
    #[arg(long, visible_alias = "keepStrand", requires = "bed")]
    keep_bed_strand: bool,
}

#[derive(Debug, Args)]
struct TrimmingArgs {
    /// Inclusive 1-based read 1 start,end and read 2 start,end for OT calls.
    #[arg(long, visible_alias = "OT", value_name = "A,B,C,D")]
    ot: Option<ReadBounds>,

    /// Inclusive 1-based read bounds for OB calls.
    #[arg(long, visible_alias = "OB", value_name = "A,B,C,D")]
    ob: Option<ReadBounds>,

    /// Inclusive 1-based read bounds for CTOT calls.
    #[arg(long, visible_alias = "CTOT", value_name = "A,B,C,D")]
    ctot: Option<ReadBounds>,

    /// Inclusive 1-based read bounds for CTOB calls.
    #[arg(long, visible_alias = "CTOB", value_name = "A,B,C,D")]
    ctob: Option<ReadBounds>,

    /// Fixed counts removed from the left,right ends of OT read 1 and read 2.
    #[arg(long = "trim-ot", visible_alias = "nOT", value_name = "A,B,C,D")]
    trim_ot: Option<ReadBounds>,

    /// Fixed end-removal counts for OB calls.
    #[arg(long = "trim-ob", visible_alias = "nOB", value_name = "A,B,C,D")]
    trim_ob: Option<ReadBounds>,

    /// Fixed end-removal counts for CTOT calls.
    #[arg(long = "trim-ctot", visible_alias = "nCTOT", value_name = "A,B,C,D")]
    trim_ctot: Option<ReadBounds>,

    /// Fixed end-removal counts for CTOB calls.
    #[arg(long = "trim-ctob", visible_alias = "nCTOB", value_name = "A,B,C,D")]
    trim_ctob: Option<ReadBounds>,
}

impl TrimmingArgs {
    fn options(&self) -> Result<TrimmingOptions> {
        let mut options = TrimmingOptions::default();
        for (strand, bounds) in [
            (BisulfiteStrand::Ot, self.ot),
            (BisulfiteStrand::Ob, self.ob),
            (BisulfiteStrand::Ctot, self.ctot),
            (BisulfiteStrand::Ctob, self.ctob),
        ] {
            if let Some(bounds) = bounds {
                options.set_inclusion(strand, bounds)?;
            }
        }
        for (strand, bounds) in [
            (BisulfiteStrand::Ot, self.trim_ot),
            (BisulfiteStrand::Ob, self.trim_ob),
            (BisulfiteStrand::Ctot, self.trim_ctot),
            (BisulfiteStrand::Ctob, self.trim_ctob),
        ] {
            if let Some(bounds) = bounds {
                options.set_fixed_ends(strand, bounds);
            }
        }
        Ok(options)
    }
}

#[derive(Debug, Args)]
struct MbiasArgs {
    /// Indexed reference FASTA.
    reference: PathBuf,

    /// Coordinate-sorted indexed BAM or CRAM.
    input: PathBuf,

    /// Prefix for the transactional TSV and strand-specific SVG outputs.
    #[arg(short, long)]
    output_prefix: PathBuf,

    /// Reference region using 1-based inclusive coordinates.
    #[arg(short = 'r', long)]
    region: Option<Region>,

    #[command(flatten)]
    filters: PileupFilterArgs,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
pub(crate) enum ExtractFormat {
    Standard,
    Fraction,
    Counts,
    Logit,
    MethylKit,
    CytosineReport,
}

#[derive(Debug, Args)]
struct MergeContextArgs {
    /// Indexed reference FASTA.
    reference: PathBuf,

    /// Coordinate-sorted six-column methylation input; use - for standard input.
    input: PathBuf,

    /// Transactional output; omit or use - for standard output.
    #[arg(short, long)]
    output: Option<PathBuf>,
}

#[derive(Debug, Args)]
struct PerReadArgs {
    /// Indexed reference FASTA.
    reference: PathBuf,

    /// Coordinate-sorted indexed BAM or CRAM.
    input: PathBuf,

    /// Reference region using 1-based inclusive coordinates.
    #[arg(short = 'r', long)]
    region: Option<Region>,

    /// Transactional output; omit or use - for standard output.
    #[arg(short, long)]
    output: Option<PathBuf>,

    #[command(flatten)]
    bed: BedArgs,

    /// Minimum alignment mapping quality.
    #[arg(short = 'q', long, default_value_t = 10)]
    minimum_mapping_quality: u8,

    /// Minimum base quality; must be positive.
    #[arg(short = 'p', long, default_value_t = 5)]
    minimum_base_quality: u8,

    /// SAM flag bits that exclude a record when any are set.
    #[arg(short = 'F', long, default_value_t = 0)]
    ignore_flags: u16,

    /// SAM flag bits that must all be set.
    #[arg(short = 'R', long, default_value_t = 0)]
    require_flags: u16,

    /// Include records with NH greater than one.
    #[arg(long)]
    ignore_nh: bool,
}

#[derive(Debug, Serialize)]
pub struct ExecutionReport {
    operation: &'static str,
    input_records: u64,
    output_records: u64,
    filtered_records: u64,
    excluded_variant_sites: u64,
    merged_records: u64,
    outputs: Vec<String>,
}

impl Cli {
    pub fn execute(self) -> Result<ExecutionReport> {
        match self.command {
            Command::Extract(args) => execute_extract(args),
            Command::Mbias(args) => execute_mbias(args),
            Command::MergeContext(args) => execute_merge_context(args, self.output.json),
            Command::PerRead(args) => execute_per_read(args, self.output.json),
        }
    }
}

fn execute_per_read(args: PerReadArgs, json: bool) -> Result<ExecutionReport> {
    if json && is_stdout(args.output.as_deref()) {
        return Err(RsomicsError::ConfigError(
            "--json requires a named per-read output".into(),
        ));
    }
    if let Some(output) = args.output.as_deref() {
        reject_output_alias(output, [args.reference.as_path(), args.input.as_path()])?;
        if let Some(bed) = args.bed.bed.as_deref() {
            reject_output_alias(output, [bed])?;
        }
    }
    let options = PerReadOptions {
        region: args.region,
        bed: args.bed.bed,
        keep_bed_strand: args.bed.keep_bed_strand,
        minimum_mapping_quality: args.minimum_mapping_quality,
        minimum_base_quality: args.minimum_base_quality,
        ignore_flags: args.ignore_flags,
        require_flags: args.require_flags,
        ignore_nh: args.ignore_nh,
    };
    let stats = write_output(args.output.as_deref(), |writer| {
        per_read(&args.input, &args.reference, options, |metric| {
            write_per_read(writer, metric)
        })
    })?;
    Ok(ExecutionReport {
        operation: "per-read",
        input_records: stats.input_records,
        output_records: stats.output_records,
        filtered_records: stats.filtered_records,
        excluded_variant_sites: 0,
        merged_records: 0,
        outputs: vec![
            args.output
                .as_deref()
                .map_or_else(|| "stdout".into(), |path| path.display().to_string()),
        ],
    })
}

fn write_per_read(writer: &mut dyn Write, metric: PerReadMetric) -> Result<()> {
    write_per_read_fields(writer, metric).map_err(RsomicsError::Io)
}

fn write_per_read_fields(writer: &mut dyn Write, metric: PerReadMetric) -> io::Result<()> {
    writer.write_all(metric.name().as_bytes())?;
    writer.write_all(b"\t")?;
    writer.write_all(metric.chromosome().as_bytes())?;
    writer.write_all(b"\t")?;
    write_integer(writer, metric.start())?;
    writer.write_all(b"\t")?;
    let informative = metric.informative_bases();
    if informative == 0 {
        writer.write_all(b"0.0\t0\n")?;
        return Ok(());
    }
    write_percentage(writer, metric.methylated(), informative)?;
    writer.write_all(b"\t")?;
    write_integer(writer, informative)?;
    writer.write_all(b"\n")
}

fn write_integer(writer: &mut dyn Write, value: u64) -> io::Result<()> {
    let mut buffer = itoa::Buffer::new();
    writer.write_all(buffer.format(value).as_bytes())
}

fn write_percentage(writer: &mut dyn Write, methylated: u64, informative: u64) -> io::Result<()> {
    let denominator = u128::from(informative);
    let numerator = u128::from(methylated) * 100_000_000;
    let mut scaled = numerator / denominator;
    let remainder = numerator % denominator;
    if remainder * 2 > denominator || remainder * 2 == denominator && scaled % 2 == 1 {
        scaled += 1;
    }
    let whole = u64::try_from(scaled / 1_000_000).expect("percentage is at most 100");
    write_integer(writer, whole)?;
    writer.write_all(b".")?;
    let mut fraction = u64::try_from(scaled % 1_000_000).expect("fraction has six digits");
    let mut digits = [b'0'; 6];
    for digit in digits.iter_mut().rev() {
        *digit += u8::try_from(fraction % 10).expect("decimal digit fits in u8");
        fraction /= 10;
    }
    writer.write_all(&digits)
}

fn execute_extract(args: ExtractArgs) -> Result<ExecutionReport> {
    let filters = args.filters;
    let trimming = filters.trimming.options()?;
    let format = if args.cytosine_report {
        ExtractFormat::CytosineReport
    } else {
        args.format
    };
    let options = ExtractOptions {
        region: args.region,
        bed: filters.bed.bed,
        keep_bed_strand: filters.bed.keep_bed_strand,
        trimming,
        minimum_mapping_quality: filters.minimum_mapping_quality,
        minimum_base_quality: filters.minimum_base_quality,
        minimum_conversion_efficiency: filters.minimum_conversion_efficiency,
        minimum_opposite_depth: args.variants.minimum_opposite_depth,
        maximum_variant_fraction: args.variants.maximum_variant_fraction,
        ignore_flags: filters.ignore_flags,
        require_flags: filters.require_flags,
        keep_duplicates: filters.keep_duplicates,
        keep_singletons: filters.keep_singletons,
        keep_discordant: filters.keep_discordant,
        ignore_nh: filters.ignore_nh,
        minimum_depth: args.minimum_depth,
        cpg: !filters.no_cpg,
        chg: filters.chg,
        chh: filters.chh,
    };
    let result = extract_to_outputs(
        &args.input,
        &args.reference,
        &args.output_prefix,
        format,
        args.merge_context,
        options,
    )?;
    Ok(ExecutionReport {
        operation: "extract",
        input_records: result.stats.input_records,
        output_records: result.output_records,
        filtered_records: result.stats.filtered_records,
        excluded_variant_sites: result.stats.excluded_variant_sites,
        merged_records: result.merged_records,
        outputs: result
            .outputs
            .iter()
            .map(|path| path.display().to_string())
            .collect(),
    })
}

fn execute_mbias(args: MbiasArgs) -> Result<ExecutionReport> {
    let filters = args.filters;
    let trimming = filters.trimming.options()?;
    let result = mbias_to_outputs(
        &args.input,
        &args.reference,
        &args.output_prefix,
        MbiasOptions {
            region: args.region,
            bed: filters.bed.bed,
            keep_bed_strand: filters.bed.keep_bed_strand,
            trimming,
            minimum_mapping_quality: filters.minimum_mapping_quality,
            minimum_base_quality: filters.minimum_base_quality,
            minimum_conversion_efficiency: filters.minimum_conversion_efficiency,
            ignore_flags: filters.ignore_flags,
            require_flags: filters.require_flags,
            keep_duplicates: filters.keep_duplicates,
            keep_singletons: filters.keep_singletons,
            keep_discordant: filters.keep_discordant,
            ignore_nh: filters.ignore_nh,
            cpg: !filters.no_cpg,
            chg: filters.chg,
            chh: filters.chh,
        },
    )?;
    Ok(ExecutionReport {
        operation: "mbias",
        input_records: result.stats.input_records,
        output_records: result.metrics,
        filtered_records: result.stats.filtered_records,
        excluded_variant_sites: 0,
        merged_records: 0,
        outputs: result
            .outputs
            .iter()
            .map(|path| path.display().to_string())
            .collect(),
    })
}

fn execute_merge_context(args: MergeContextArgs, json: bool) -> Result<ExecutionReport> {
    if json && is_stdout(args.output.as_deref()) {
        return Err(RsomicsError::ConfigError(
            "--json requires a named methylation output".into(),
        ));
    }
    if let Some(output) = args.output.as_deref() {
        reject_output_alias(output, [args.reference.as_path(), args.input.as_path()])?;
    }
    let stats = if args.input == Path::new("-") {
        let stdin = io::stdin();
        write_output(args.output.as_deref(), |writer| {
            merge_context(stdin.lock(), writer, &args.reference)
        })?
    } else {
        let input = File::open(&args.input).map_err(|error| {
            RsomicsError::InvalidInput(format!("opening {}: {error}", args.input.display()))
        })?;
        write_output(args.output.as_deref(), |writer| {
            merge_context(BufReader::new(input), writer, &args.reference)
        })?
    };
    Ok(report(stats, args.output.as_deref()))
}

fn report(stats: MergeContextStats, output: Option<&Path>) -> ExecutionReport {
    ExecutionReport {
        operation: "merge-context",
        input_records: stats.input_records,
        output_records: stats.output_records,
        filtered_records: 0,
        excluded_variant_sites: 0,
        merged_records: stats.merged_records,
        outputs: vec![output.map_or_else(|| "stdout".into(), |path| path.display().to_string())],
    }
}

fn is_stdout(output: Option<&Path>) -> bool {
    output.is_none_or(|path| path == Path::new("-"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn integer_percentage_matches_six_decimal_float_output() {
        for informative in 1..=1024 {
            for methylated in 0..=informative {
                let mut actual = Vec::new();
                write_percentage(&mut actual, methylated, informative).unwrap();
                let expected = format!("{:.6}", methylated as f64 * 100.0 / informative as f64);
                assert_eq!(actual, expected.as_bytes());
            }
        }
    }
}
