mod cli;
mod extract_output;
mod mbias_output;
mod output;

use cli::{Cli, META};

fn main() -> std::process::ExitCode {
    let cli = rsomics_help::parse::<Cli>();
    let output = cli.output.clone();
    rsomics_common::run(&output, META, || cli.execute())
}
