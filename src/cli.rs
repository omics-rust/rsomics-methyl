use std::fs::File;
use std::io::{self, BufReader};
use std::path::{Path, PathBuf};

use clap::{Args, Parser, Subcommand};
use rsomics_common::{
    OutputArgs, Result, RsomicsError, ToolMeta, reject_output_alias, write_output,
};
use rsomics_methyl::merge_context::{MergeContextStats, merge_context};
use serde::Serialize;

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
    /// Merge strand-specific CpG and CHG cytosine metrics.
    MergeContext(MergeContextArgs),
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
    merged_records: u64,
    output: String,
}

impl Cli {
    pub fn execute(self) -> Result<ExecutionReport> {
        match self.command {
            Command::MergeContext(args) => execute_merge_context(args, self.output.json),
        }
    }
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
        merged_records: stats.merged_records,
        output: output.map_or_else(|| "stdout".into(), |path| path.display().to_string()),
    }
}

fn is_stdout(output: Option<&Path>) -> bool {
    output.is_none_or(|path| path == Path::new("-"))
}
