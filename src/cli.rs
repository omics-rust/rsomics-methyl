use std::fs::File;
use std::io::{self, BufReader};
use std::path::{Path, PathBuf};

use clap::{Args, Parser, Subcommand, ValueEnum};
use rsomics_common::{
    OutputArgs, Result, RsomicsError, ToolMeta, reject_output_alias, write_output,
};
use rsomics_methyl::extract::ExtractOptions;
use rsomics_methyl::merge_context::{MergeContextStats, merge_context};
use serde::Serialize;

use crate::extract_output::extract_to_standard_outputs;

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
}

#[derive(Debug, Args)]
struct ExtractArgs {
    /// Indexed reference FASTA.
    reference: PathBuf,

    /// Coordinate-sorted indexed BAM or CRAM.
    input: PathBuf,

    /// Prefix for transactional context-specific bedGraph outputs.
    #[arg(short, long)]
    output_prefix: PathBuf,

    /// Output representation.
    #[arg(long, value_enum, default_value = "standard")]
    format: ExtractFormat,

    /// Minimum alignment mapping quality.
    #[arg(short = 'q', long, default_value_t = 10)]
    minimum_mapping_quality: u8,

    /// Minimum base quality; must be positive.
    #[arg(short = 'p', long, default_value_t = 5)]
    minimum_base_quality: u8,

    /// Minimum methylated plus unmethylated depth.
    #[arg(short = 'd', long, default_value_t = 1)]
    minimum_depth: u64,

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

    /// Disable CpG output.
    #[arg(long)]
    no_cpg: bool,

    /// Emit CHG output.
    #[arg(long)]
    chg: bool,

    /// Emit CHH output.
    #[arg(long)]
    chh: bool,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
pub(crate) enum ExtractFormat {
    Standard,
    Fraction,
    Counts,
    Logit,
    MethylKit,
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

#[derive(Debug, Serialize)]
pub struct ExecutionReport {
    operation: &'static str,
    input_records: u64,
    output_records: u64,
    filtered_records: u64,
    merged_records: u64,
    outputs: Vec<String>,
}

impl Cli {
    pub fn execute(self) -> Result<ExecutionReport> {
        match self.command {
            Command::Extract(args) => execute_extract(args),
            Command::MergeContext(args) => execute_merge_context(args, self.output.json),
        }
    }
}

fn execute_extract(args: ExtractArgs) -> Result<ExecutionReport> {
    let options = ExtractOptions {
        minimum_mapping_quality: args.minimum_mapping_quality,
        minimum_base_quality: args.minimum_base_quality,
        ignore_flags: args.ignore_flags,
        require_flags: args.require_flags,
        keep_duplicates: args.keep_duplicates,
        keep_singletons: args.keep_singletons,
        keep_discordant: args.keep_discordant,
        ignore_nh: args.ignore_nh,
        minimum_depth: args.minimum_depth,
        cpg: !args.no_cpg,
        chg: args.chg,
        chh: args.chh,
    };
    let result = extract_to_standard_outputs(
        &args.input,
        &args.reference,
        &args.output_prefix,
        args.format,
        options,
    )?;
    Ok(ExecutionReport {
        operation: "extract",
        input_records: result.stats.input_records,
        output_records: result.stats.emitted_sites,
        filtered_records: result.stats.filtered_records,
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
        merged_records: stats.merged_records,
        outputs: vec![output.map_or_else(|| "stdout".into(), |path| path.display().to_string())],
    }
}

fn is_stdout(output: Option<&Path>) -> bool {
    output.is_none_or(|path| path == Path::new("-"))
}
