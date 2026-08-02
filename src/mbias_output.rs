use std::ffi::OsString;
use std::io::Write;
use std::path::{Path, PathBuf};

use rsomics_common::{Result, RsomicsError, reject_output_alias};
use rsomics_methyl::BisulfiteStrand;
use rsomics_methyl::mbias::{MbiasMetric, MbiasOptions, MbiasStats, MbiasSuggestion, mbias};

use crate::output::{TransactionalOutput, commit_all};

pub(crate) struct MbiasOutputResult {
    pub(crate) stats: MbiasStats,
    pub(crate) metrics: u64,
    pub(crate) outputs: Vec<PathBuf>,
}

pub(crate) fn mbias_to_outputs(
    input: &Path,
    reference: &Path,
    prefix: &Path,
    options: MbiasOptions,
) -> Result<MbiasOutputResult> {
    let bed = options.bed.clone();
    let result = mbias(input, reference, options)?;
    let mut files = Vec::new();
    let table_path = suffix_path(prefix, "_mbias.tsv");
    reject_output_alias(&table_path, [input, reference])?;
    if let Some(bed) = bed.as_deref() {
        reject_output_alias(&table_path, [bed])?;
    }
    files.push(OutputFile {
        strand: None,
        output: TransactionalOutput::new(&table_path)?,
    });
    for suggestion in result.suggestions() {
        let path = suffix_path(prefix, &format!("_{}.svg", suggestion.strand().label()));
        reject_output_alias(&path, [input, reference])?;
        if let Some(bed) = bed.as_deref() {
            reject_output_alias(&path, [bed])?;
        }
        reject_output_alias(
            &path,
            files.iter().map(|file: &OutputFile| file.output.path()),
        )?;
        files.push(OutputFile {
            strand: Some(suggestion.strand()),
            output: TransactionalOutput::new(&path)?,
        });
    }
    write_table(files[0].output.writer(), result.metrics())?;
    for file in files.iter_mut().skip(1) {
        let strand = file
            .strand
            .expect("SVG output has an associated bisulfite strand");
        let suggestion = result
            .suggestions()
            .iter()
            .find(|suggestion| suggestion.strand() == strand)
            .expect("each SVG strand has an inclusion suggestion");
        write_svg(file.output.writer(), strand, result.metrics(), suggestion)?;
    }
    commit_all(&mut files, |file| &mut file.output)?;
    let metrics = u64::try_from(result.metrics().len())
        .map_err(|error| RsomicsError::InvalidInput(error.to_string()))?;
    Ok(MbiasOutputResult {
        stats: result.stats().clone(),
        metrics,
        outputs: files
            .iter()
            .map(|file| file.output.path().to_owned())
            .collect(),
    })
}

struct OutputFile {
    strand: Option<BisulfiteStrand>,
    output: TransactionalOutput,
}

fn write_table(writer: &mut dyn Write, metrics: &[MbiasMetric]) -> Result<()> {
    writeln!(writer, "Strand\tRead\tPosition\tnMethylated\tnUnmethylated")
        .map_err(RsomicsError::Io)?;
    for metric in metrics {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}",
            metric.strand().label(),
            metric.read(),
            metric.position(),
            metric.methylated(),
            metric.unmethylated()
        )
        .map_err(RsomicsError::Io)?;
    }
    Ok(())
}

fn write_svg(
    writer: &mut dyn Write,
    strand: BisulfiteStrand,
    metrics: &[MbiasMetric],
    suggestion: &MbiasSuggestion,
) -> Result<()> {
    let metrics = metrics
        .iter()
        .filter(|metric| metric.strand() == strand)
        .collect::<Vec<_>>();
    let maximum = metrics
        .iter()
        .map(|metric| metric.position())
        .max()
        .ok_or_else(|| RsomicsError::InvalidInput("M-bias SVG has no metrics".into()))?;
    writeln!(
        writer,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 720 640\">"
    )
    .map_err(RsomicsError::Io)?;
    writeln!(writer, "<title>{} M-bias</title>", strand.label()).map_err(RsomicsError::Io)?;
    writeln!(
        writer,
        "<rect width=\"720\" height=\"640\" fill=\"white\"/>"
    )
    .map_err(RsomicsError::Io)?;
    writeln!(
        writer,
        "<text x=\"360\" y=\"28\" text-anchor=\"middle\">{} strand</text>",
        strand.label()
    )
    .map_err(RsomicsError::Io)?;
    writeln!(
        writer,
        "<path d=\"M 70 50 V 570 H 680\" fill=\"none\" stroke=\"black\"/>"
    )
    .map_err(RsomicsError::Io)?;
    for percentage in [0u32, 25, 50, 75, 100] {
        let y = y(f64::from(percentage) / 100.0);
        writeln!(
            writer,
            "<line x1=\"65\" y1=\"{y:.3}\" x2=\"680\" y2=\"{y:.3}\" stroke=\"#dddddd\"/>"
        )
        .map_err(RsomicsError::Io)?;
        writeln!(
            writer,
            "<text x=\"58\" y=\"{:.3}\" text-anchor=\"end\">{percentage}%</text>",
            y + 4.0
        )
        .map_err(RsomicsError::Io)?;
    }
    for (read, color) in [(1, "#f8766d"), (2, "#00bfc4")] {
        let values = metrics
            .iter()
            .copied()
            .filter(|metric| metric.read() == read)
            .collect::<Vec<_>>();
        if values.is_empty() {
            continue;
        }
        let mut band = values
            .iter()
            .map(|metric| {
                let (lower, _) = metric.confidence_interval();
                format!("{:.3},{:.3}", x(metric.position(), maximum), y(lower))
            })
            .collect::<Vec<_>>();
        band.extend(values.iter().rev().map(|metric| {
            let (_, upper) = metric.confidence_interval();
            format!("{:.3},{:.3}", x(metric.position(), maximum), y(upper))
        }));
        writeln!(
            writer,
            "<polygon points=\"{}\" fill=\"{color}\" fill-opacity=\"0.18\"/>",
            band.join(" ")
        )
        .map_err(RsomicsError::Io)?;
        let points = values
            .iter()
            .map(|metric| {
                format!(
                    "{:.3},{:.3}",
                    x(metric.position(), maximum),
                    y(metric.fraction())
                )
            })
            .collect::<Vec<_>>()
            .join(" ");
        writeln!(writer, "<polyline id=\"read-{read}\" points=\"{points}\" fill=\"none\" stroke=\"{color}\" stroke-width=\"2\"/>")
            .map_err(RsomicsError::Io)?;
    }
    for bound in suggestion.bounds().into_iter().filter(|bound| *bound > 0) {
        let x = x(bound, maximum);
        writeln!(writer, "<line x1=\"{x:.3}\" y1=\"50\" x2=\"{x:.3}\" y2=\"570\" stroke=\"#555555\" stroke-dasharray=\"5 3\"/>")
            .map_err(RsomicsError::Io)?;
    }
    let bounds = suggestion.bounds();
    writeln!(
        writer,
        "<text x=\"680\" y=\"610\" text-anchor=\"end\">--{} {},{},{},{}</text>",
        strand.label(),
        bounds[0],
        bounds[1],
        bounds[2],
        bounds[3]
    )
    .map_err(RsomicsError::Io)?;
    writeln!(
        writer,
        "<text x=\"375\" y=\"632\" text-anchor=\"middle\">Position along mapped read</text>"
    )
    .map_err(RsomicsError::Io)?;
    writeln!(writer, "</svg>").map_err(RsomicsError::Io)
}

fn x(position: u64, maximum: u64) -> f64 {
    70.0 + 610.0 * position as f64 / maximum as f64
}

fn y(fraction: f64) -> f64 {
    570.0 - 520.0 * fraction
}

fn suffix_path(prefix: &Path, suffix: &str) -> PathBuf {
    let mut path = OsString::from(prefix.as_os_str());
    path.push(suffix);
    PathBuf::from(path)
}
