use std::collections::{BTreeMap, BTreeSet};
use std::ffi::OsString;
use std::io::Write;
use std::path::{Path, PathBuf};

use rsomics_common::{Result, RsomicsError, reject_output_alias};
use rsomics_methyl::extract::{
    ExcludedVariantSite, ExtractEvent, ExtractOptions, ExtractStats, SiteMetric, extract,
    extract_all_cytosines, extract_events,
};
use rsomics_methyl::{ReferenceStrand, SequenceContext};

use crate::cli::ExtractFormat;
use crate::output::{TransactionalOutput, commit_all};

pub struct ExtractOutputResult {
    pub stats: ExtractStats,
    pub output_records: u64,
    pub merged_records: u64,
    pub outputs: Vec<PathBuf>,
}

struct CytosineReportOutput {
    output: TransactionalOutput,
    records: u64,
}

impl CytosineReportOutput {
    fn new(prefix: &Path, input: &Path, reference: &Path, bed: Option<&Path>) -> Result<Self> {
        let path = cytosine_report_path(prefix);
        reject_output_alias(&path, [input, reference])?;
        if let Some(bed) = bed {
            reject_output_alias(&path, [bed])?;
        }
        Ok(Self {
            output: TransactionalOutput::new(&path)?,
            records: 0,
        })
    }

    fn write(&mut self, metric: &SiteMetric) -> Result<()> {
        let strand = match metric.strand() {
            ReferenceStrand::Forward => '+',
            ReferenceStrand::Reverse => '-',
        };
        let context = match metric.context() {
            SequenceContext::Cpg => "CG",
            SequenceContext::Chg => "CHG",
            SequenceContext::Chh => "CHH",
        };
        let trinucleotide = metric.trinucleotide();
        let trinucleotide = std::str::from_utf8(&trinucleotide)
            .expect("trinucleotide contexts contain only ASCII bases");
        writeln!(
            self.output.writer(),
            "{}\t{}\t{strand}\t{}\t{}\t{context}\t{trinucleotide}",
            metric.chromosome(),
            metric.end(),
            metric.methylated(),
            metric.unmethylated()
        )
        .map_err(RsomicsError::Io)?;
        self.records = self
            .records
            .checked_add(1)
            .ok_or_else(|| RsomicsError::InvalidInput("output record count overflows".into()))?;
        Ok(())
    }

    fn commit(mut self) -> Result<u64> {
        commit_all(std::slice::from_mut(&mut self.output), |output| output)?;
        Ok(self.records)
    }
}

pub fn extract_to_outputs(
    input: &Path,
    reference: &Path,
    prefix: &Path,
    format: ExtractFormat,
    merge_context: bool,
    mut options: ExtractOptions,
) -> Result<ExtractOutputResult> {
    if merge_context
        && matches!(
            format,
            ExtractFormat::MethylKit | ExtractFormat::CytosineReport
        )
    {
        return Err(RsomicsError::ConfigError(format!(
            "{} output cannot merge complementary contexts",
            format.label()
        )));
    }
    if matches!(format, ExtractFormat::CytosineReport) {
        let mut output =
            CytosineReportOutput::new(prefix, input, reference, options.bed.as_deref())?;
        let path = output.output.path().to_owned();
        let stats =
            extract_all_cytosines(input, reference, options, |metric| output.write(&metric))?;
        let output_records = output.commit()?;
        return Ok(ExtractOutputResult {
            stats,
            output_records,
            merged_records: 0,
            outputs: vec![path],
        });
    }
    let minimum_depth = options.minimum_depth;
    if merge_context {
        options.minimum_depth = 1;
    }
    let mut outputs = ContextOutputs::new(
        prefix,
        input,
        reference,
        format,
        merge_context,
        minimum_depth,
        &options,
    )?;
    let paths = outputs.paths();
    let stats = if merge_context {
        extract_events(input, reference, options, |event| {
            outputs.write_event(&event)
        })?
    } else {
        extract(input, reference, options, |metric| outputs.write(&metric))?
    };
    let output_stats = outputs.commit()?;
    Ok(ExtractOutputResult {
        stats,
        output_records: output_stats.output_records,
        merged_records: output_stats.merged_records,
        outputs: paths,
    })
}

struct ContextOutputs {
    format: ExtractFormat,
    merge_context: bool,
    minimum_depth: u64,
    stats: OutputStats,
    entries: Vec<OutputEntry>,
}

#[derive(Clone, Copy, Debug, Default)]
struct OutputStats {
    output_records: u64,
    merged_records: u64,
}

impl OutputStats {
    fn add(&mut self, other: Self) -> Result<()> {
        self.output_records = self
            .output_records
            .checked_add(other.output_records)
            .ok_or_else(|| RsomicsError::InvalidInput("output record count overflows".into()))?;
        self.merged_records = self
            .merged_records
            .checked_add(other.merged_records)
            .ok_or_else(|| RsomicsError::InvalidInput("merged record count overflows".into()))?;
        Ok(())
    }
}

impl ContextOutputs {
    fn new(
        prefix: &Path,
        input: &Path,
        reference: &Path,
        format: ExtractFormat,
        merge_context: bool,
        minimum_depth: u64,
        options: &ExtractOptions,
    ) -> Result<Self> {
        let mut entries = Vec::new();
        for (enabled, context, label) in [
            (options.cpg, SequenceContext::Cpg, "CpG"),
            (options.chg, SequenceContext::Chg, "CHG"),
            (options.chh, SequenceContext::Chh, "CHH"),
        ] {
            if !enabled {
                continue;
            }
            let path = context_path(prefix, label, format);
            reject_output_alias(&path, [input, reference])?;
            if let Some(bed) = options.bed.as_deref() {
                reject_output_alias(&path, [bed])?;
            }
            reject_output_alias(
                &path,
                entries
                    .iter()
                    .map(|entry: &OutputEntry| entry.output.path()),
            )?;
            entries.push(OutputEntry::new(path, context, label)?);
        }
        if entries.is_empty() {
            return Err(RsomicsError::ConfigError(
                "at least one methylation context must be enabled".into(),
            ));
        }
        let mut outputs = Self {
            format,
            merge_context,
            minimum_depth,
            stats: OutputStats::default(),
            entries,
        };
        for entry in &mut outputs.entries {
            let label = entry.label;
            format.write_header(entry.file(), prefix, label, merge_context)?;
        }
        Ok(outputs)
    }

    fn paths(&self) -> Vec<PathBuf> {
        self.entries
            .iter()
            .map(|entry| entry.output.path().to_owned())
            .collect()
    }

    fn write(&mut self, metric: &SiteMetric) -> Result<()> {
        let format = self.format;
        let merge_context = self.merge_context;
        let minimum_depth = self.minimum_depth;
        let entry = self
            .entries
            .iter_mut()
            .find(|entry| entry.context == metric.context())
            .ok_or_else(|| RsomicsError::ConfigError("methylation context has no output".into()))?;
        self.stats
            .add(entry.push(format, merge_context, minimum_depth, metric)?)
    }

    fn write_event(&mut self, event: &ExtractEvent) -> Result<()> {
        match event {
            ExtractEvent::Site(metric) => self.write(metric),
            ExtractEvent::ExcludedVariant(site) => {
                let format = self.format;
                let minimum_depth = self.minimum_depth;
                let entry = self
                    .entries
                    .iter_mut()
                    .find(|entry| entry.context == site.context())
                    .ok_or_else(|| {
                        RsomicsError::ConfigError("methylation context has no output".into())
                    })?;
                self.stats.add(entry.exclude(format, minimum_depth, site)?)
            }
        }
    }

    fn commit(mut self) -> Result<OutputStats> {
        for entry in &mut self.entries {
            self.stats
                .add(entry.finish(self.format, self.minimum_depth)?)?;
        }
        commit_all(&mut self.entries, |entry| &mut entry.output)?;
        Ok(self.stats)
    }
}

struct OutputEntry {
    context: SequenceContext,
    label: &'static str,
    chromosome: Option<String>,
    pending: BTreeMap<(u64, u64), OutputMetric>,
    excluded: BTreeSet<(u64, u64)>,
    line: Vec<u8>,
    output: TransactionalOutput,
}

impl OutputEntry {
    fn new(path: PathBuf, context: SequenceContext, label: &'static str) -> Result<Self> {
        let output = TransactionalOutput::new(&path)?;
        Ok(Self {
            context,
            label,
            chromosome: None,
            pending: BTreeMap::new(),
            excluded: BTreeSet::new(),
            line: Vec::with_capacity(128),
            output,
        })
    }

    fn file(&mut self) -> &mut dyn Write {
        self.output.writer()
    }

    fn push(
        &mut self,
        format: ExtractFormat,
        merge_context: bool,
        minimum_depth: u64,
        metric: &SiteMetric,
    ) -> Result<OutputStats> {
        let mut stats = OutputStats::default();
        if !merge_context || metric.context() == SequenceContext::Chh {
            stats.add(self.write_site_metric(format, minimum_depth, metric)?)?;
            return Ok(stats);
        }
        let mut metric = OutputMetric::from(metric);
        if self.chromosome.as_deref() != Some(metric.chromosome.as_str()) {
            stats.add(self.flush_all(format, minimum_depth)?)?;
            self.chromosome = Some(metric.chromosome.clone());
        }
        let position = metric.start;
        metric.apply_merged_span()?;
        let key = (metric.start, metric.end);
        if !self.excluded.contains(&key) {
            if let Some(existing) = self.pending.get_mut(&key) {
                existing.merge(&metric)?;
                stats.merged_records = 1;
            } else {
                self.pending.insert(key, metric);
            }
        }
        let settled = position
            .checked_add(1)
            .ok_or_else(|| RsomicsError::InvalidInput("metric coordinate overflows".into()))?;
        stats.add(self.flush_settled(format, minimum_depth, settled)?)?;
        Ok(stats)
    }

    fn exclude(
        &mut self,
        format: ExtractFormat,
        minimum_depth: u64,
        site: &ExcludedVariantSite,
    ) -> Result<OutputStats> {
        if site.context() == SequenceContext::Chh {
            return Ok(OutputStats::default());
        }
        let mut stats = OutputStats::default();
        if self.chromosome.as_deref() != Some(site.chromosome()) {
            stats.add(self.flush_all(format, minimum_depth)?)?;
            self.chromosome = Some(site.chromosome().to_owned());
        }
        let mut metric = OutputMetric::from_variant(site)?;
        let position = metric.start;
        metric.apply_merged_span()?;
        let key = (metric.start, metric.end);
        self.pending.remove(&key);
        self.excluded.insert(key);
        let settled = position
            .checked_add(1)
            .ok_or_else(|| RsomicsError::InvalidInput("metric coordinate overflows".into()))?;
        stats.add(self.flush_settled(format, minimum_depth, settled)?)?;
        Ok(stats)
    }

    fn finish(&mut self, format: ExtractFormat, minimum_depth: u64) -> Result<OutputStats> {
        self.flush_all(format, minimum_depth)
    }

    fn flush_settled(
        &mut self,
        format: ExtractFormat,
        minimum_depth: u64,
        settled: u64,
    ) -> Result<OutputStats> {
        let mut stats = OutputStats::default();
        while self
            .pending
            .first_key_value()
            .is_some_and(|(_, metric)| metric.end <= settled)
        {
            let (_, metric) = self
                .pending
                .pop_first()
                .expect("pending output is present after inspection");
            stats.add(self.write_metric(format, minimum_depth, metric)?)?;
        }
        self.excluded.retain(|(_, end)| *end > settled);
        Ok(stats)
    }

    fn flush_all(&mut self, format: ExtractFormat, minimum_depth: u64) -> Result<OutputStats> {
        let mut stats = OutputStats::default();
        while let Some((_, metric)) = self.pending.pop_first() {
            stats.add(self.write_metric(format, minimum_depth, metric)?)?;
        }
        self.excluded.clear();
        Ok(stats)
    }

    fn write_metric(
        &mut self,
        format: ExtractFormat,
        minimum_depth: u64,
        metric: OutputMetric,
    ) -> Result<OutputStats> {
        let depth = metric.depth()?;
        if depth < minimum_depth {
            return Ok(OutputStats::default());
        }
        let values = MetricValues {
            chromosome: &metric.chromosome,
            start: metric.start,
            end: metric.end,
            strand: metric.strand,
            methylated: metric.methylated,
            unmethylated: metric.unmethylated,
            depth,
        };
        if matches!(format, ExtractFormat::Standard) {
            write_standard(&mut self.output, &mut self.line, values)?;
        } else {
            format.write_values(self.file(), values)?;
        }
        Ok(OutputStats {
            output_records: 1,
            merged_records: 0,
        })
    }

    fn write_site_metric(
        &mut self,
        format: ExtractFormat,
        minimum_depth: u64,
        metric: &SiteMetric,
    ) -> Result<OutputStats> {
        let depth = metric.depth();
        if depth < minimum_depth {
            return Ok(OutputStats::default());
        }
        let values = MetricValues {
            chromosome: metric.chromosome(),
            start: metric.start(),
            end: metric.end(),
            strand: metric.strand(),
            methylated: metric.methylated(),
            unmethylated: metric.unmethylated(),
            depth,
        };
        if matches!(format, ExtractFormat::Standard) {
            write_standard(&mut self.output, &mut self.line, values)?;
        } else {
            format.write_values(self.file(), values)?;
        }
        Ok(OutputStats {
            output_records: 1,
            merged_records: 0,
        })
    }
}

struct OutputMetric {
    chromosome: String,
    start: u64,
    end: u64,
    context: SequenceContext,
    strand: ReferenceStrand,
    methylated: u64,
    unmethylated: u64,
}

#[derive(Clone, Copy)]
struct MetricValues<'a> {
    chromosome: &'a str,
    start: u64,
    end: u64,
    strand: ReferenceStrand,
    methylated: u64,
    unmethylated: u64,
    depth: u64,
}

fn write_standard(
    output: &mut TransactionalOutput,
    line: &mut Vec<u8>,
    metric: MetricValues<'_>,
) -> Result<()> {
    line.clear();
    line.extend_from_slice(metric.chromosome.as_bytes());
    let mut buffer = itoa::Buffer::new();
    for value in [
        metric.start,
        metric.end,
        percentage(metric.methylated, metric.depth),
        metric.methylated,
        metric.unmethylated,
    ] {
        line.push(b'\t');
        line.extend_from_slice(buffer.format(value).as_bytes());
    }
    line.push(b'\n');
    output.writer().write_all(line).map_err(RsomicsError::Io)
}

impl From<&SiteMetric> for OutputMetric {
    fn from(metric: &SiteMetric) -> Self {
        Self {
            chromosome: metric.chromosome().to_owned(),
            start: metric.start(),
            end: metric.end(),
            context: metric.context(),
            strand: metric.strand(),
            methylated: metric.methylated(),
            unmethylated: metric.unmethylated(),
        }
    }
}

impl OutputMetric {
    fn from_variant(site: &ExcludedVariantSite) -> Result<Self> {
        Ok(Self {
            chromosome: site.chromosome().to_owned(),
            start: site.start(),
            end: site
                .start()
                .checked_add(1)
                .ok_or_else(|| RsomicsError::InvalidInput("site end overflows".into()))?,
            context: site.context(),
            strand: site.strand(),
            methylated: 0,
            unmethylated: 0,
        })
    }

    fn depth(&self) -> Result<u64> {
        self.methylated
            .checked_add(self.unmethylated)
            .ok_or_else(|| RsomicsError::InvalidInput("methylation depth overflows".into()))
    }

    fn apply_merged_span(&mut self) -> Result<()> {
        let offset = match self.context {
            SequenceContext::Cpg => 1,
            SequenceContext::Chg => 2,
            SequenceContext::Chh => 0,
        };
        match self.strand {
            ReferenceStrand::Forward => {
                self.end = self
                    .start
                    .checked_add(offset + 1)
                    .ok_or_else(|| RsomicsError::InvalidInput("context end overflows".into()))?;
            }
            ReferenceStrand::Reverse => {
                self.start = self.start.checked_sub(offset).ok_or_else(|| {
                    RsomicsError::InvalidInput("reverse context start underflows".into())
                })?;
            }
        }
        Ok(())
    }

    fn merge(&mut self, other: &Self) -> Result<()> {
        self.methylated = self
            .methylated
            .checked_add(other.methylated)
            .ok_or_else(|| RsomicsError::InvalidInput("methylated count overflows".into()))?;
        self.unmethylated = self
            .unmethylated
            .checked_add(other.unmethylated)
            .ok_or_else(|| RsomicsError::InvalidInput("unmethylated count overflows".into()))?;
        Ok(())
    }
}

impl ExtractFormat {
    fn label(self) -> &'static str {
        match self {
            Self::Standard => "standard",
            Self::Fraction => "fraction",
            Self::Counts => "counts",
            Self::Logit => "logit",
            Self::MethylKit => "methylKit",
            Self::CytosineReport => "cytosine-report",
        }
    }

    fn suffix(self, context: &str) -> String {
        match self {
            Self::Standard => format!("_{context}.bedGraph"),
            Self::Fraction => format!("_{context}.meth.bedGraph"),
            Self::Counts => format!("_{context}.counts.bedGraph"),
            Self::Logit => format!("_{context}.logit.bedGraph"),
            Self::MethylKit => format!("_{context}.methylKit"),
            Self::CytosineReport => {
                unreachable!("cytosine report has one context-independent path")
            }
        }
    }

    fn write_header(
        self,
        writer: &mut dyn Write,
        prefix: &Path,
        context: &str,
        merged: bool,
    ) -> Result<()> {
        if matches!(self, Self::MethylKit) {
            writeln!(writer, "chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT")
                .map_err(RsomicsError::Io)
        } else {
            let description = match self {
                Self::Standard => "methylation levels",
                Self::Fraction => "methylation fractions",
                Self::Counts => "methylation counts",
                Self::Logit => "logit transformed methylation fractions",
                Self::MethylKit => unreachable!(),
                Self::CytosineReport => unreachable!(),
            };
            let merged = if merged { " merged" } else { "" };
            writeln!(
                writer,
                "track type=\"bedGraph\" description=\"{} {context}{merged} {description}\"",
                prefix.display()
            )
            .map_err(RsomicsError::Io)
        }
    }

    fn write_values(self, writer: &mut dyn Write, metric: MetricValues<'_>) -> Result<()> {
        let MetricValues {
            chromosome,
            start,
            end,
            strand,
            methylated,
            unmethylated,
            depth,
        } = metric;
        match self {
            Self::Standard => writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}",
                chromosome,
                start,
                end,
                percentage(methylated, depth),
                methylated,
                unmethylated
            ),
            Self::Fraction => {
                let fraction = methylated as f64 / depth as f64;
                writeln!(writer, "{chromosome}\t{start}\t{end}\t{fraction:.6}")
            }
            Self::Counts => writeln!(writer, "{chromosome}\t{start}\t{end}\t{depth}"),
            Self::Logit => {
                let fraction = methylated as f64 / depth as f64;
                writeln!(
                    writer,
                    "{chromosome}\t{start}\t{end}\t{:.6}",
                    fraction.ln() - (-fraction).ln_1p()
                )
            }
            Self::MethylKit => {
                let fraction = methylated as f64 / depth as f64;
                let position = end;
                let strand = match strand {
                    ReferenceStrand::Forward => 'F',
                    ReferenceStrand::Reverse => 'R',
                };
                writeln!(
                    writer,
                    "{0}.{1}\t{0}\t{1}\t{strand}\t{depth}\t{2:6.2}\t{3:6.2}",
                    chromosome,
                    position,
                    fraction * 100.0,
                    (1.0 - fraction) * 100.0
                )
            }
            Self::CytosineReport => unreachable!("cytosine reports use their dedicated writer"),
        }
        .map_err(RsomicsError::Io)
    }
}

fn percentage(methylated: u64, depth: u64) -> u64 {
    if depth == 0 {
        0
    } else {
        (u128::from(methylated) * 100 / u128::from(depth)) as u64
    }
}

fn context_path(prefix: &Path, context: &str, format: ExtractFormat) -> PathBuf {
    let mut path = OsString::from(prefix.as_os_str());
    path.push(format.suffix(context));
    PathBuf::from(path)
}

fn cytosine_report_path(prefix: &Path) -> PathBuf {
    let mut path = OsString::from(prefix.as_os_str());
    path.push(".cytosine_report.txt");
    PathBuf::from(path)
}

#[cfg(test)]
mod tests {
    use std::fs;

    use super::*;

    #[test]
    fn later_commit_failure_restores_an_earlier_output() {
        let directory = tempfile::tempdir().unwrap();
        let prefix = directory.path().join("result");
        let first = context_path(&prefix, "CpG", ExtractFormat::Standard);
        fs::write(&first, b"old\n").unwrap();
        let mut outputs = ContextOutputs::new(
            &prefix,
            Path::new("input.bam"),
            Path::new("reference.fa"),
            ExtractFormat::Standard,
            false,
            1,
            &ExtractOptions {
                chg: true,
                ..ExtractOptions::default()
            },
        )
        .unwrap();
        for entry in &mut outputs.entries {
            entry.file().write_all(b"new\n").unwrap();
        }
        let second = context_path(&prefix, "CHG", ExtractFormat::Standard);
        fs::create_dir(&second).unwrap();

        assert!(outputs.commit().is_err());
        assert_eq!(fs::read(first).unwrap(), b"old\n");
        assert!(second.is_dir());
    }
}
