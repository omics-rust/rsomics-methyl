use std::collections::BTreeMap;
use std::ffi::OsString;
use std::fs;
use std::io::{self, Write};
#[cfg(unix)]
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};

use rsomics_common::{Context, Result, RsomicsError, reject_output_alias};
use rsomics_methyl::extract::{ExtractOptions, ExtractStats, SiteMetric, extract};
use rsomics_methyl::{ReferenceStrand, SequenceContext};
use tempfile::{Builder, NamedTempFile};

use crate::cli::ExtractFormat;

pub struct ExtractOutputResult {
    pub stats: ExtractStats,
    pub output_records: u64,
    pub merged_records: u64,
    pub outputs: Vec<PathBuf>,
}

pub fn extract_to_standard_outputs(
    input: &Path,
    reference: &Path,
    prefix: &Path,
    format: ExtractFormat,
    merge_context: bool,
    mut options: ExtractOptions,
) -> Result<ExtractOutputResult> {
    if merge_context && matches!(format, ExtractFormat::MethylKit) {
        return Err(RsomicsError::ConfigError(
            "methylKit output cannot merge complementary contexts".into(),
        ));
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
    let stats = extract(input, reference, options, |metric| outputs.write(&metric))?;
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
            reject_output_alias(
                &path,
                entries
                    .iter()
                    .map(|entry: &OutputEntry| entry.path.as_path()),
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
            .map(|entry| entry.path.clone())
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

    fn commit(mut self) -> Result<OutputStats> {
        for entry in &mut self.entries {
            self.stats
                .add(entry.finish(self.format, self.minimum_depth)?)?;
            entry.staged_mut().prepare()?;
        }
        for index in 0..self.entries.len() {
            let staged = self.entries[index]
                .staged
                .take()
                .expect("staged output is present before commit");
            if let Err(error) = staged.commit() {
                return Err(self.restore(index, error));
            }
        }
        Ok(self.stats)
    }

    fn restore(&mut self, through: usize, mut cause: RsomicsError) -> RsomicsError {
        for entry in self.entries[..=through].iter_mut().rev() {
            if let Some(backup) = entry.backup.take() {
                cause = backup.restore(&entry.path, cause);
            }
        }
        cause
    }
}

struct OutputEntry {
    path: PathBuf,
    context: SequenceContext,
    label: &'static str,
    chromosome: Option<String>,
    pending: BTreeMap<(u64, u64), OutputMetric>,
    staged: Option<Staged>,
    backup: Option<Backup>,
}

impl OutputEntry {
    fn new(path: PathBuf, context: SequenceContext, label: &'static str) -> Result<Self> {
        let staged = Staged::new(&path)?;
        let backup = Backup::new(&path)?;
        Ok(Self {
            path,
            context,
            label,
            chromosome: None,
            pending: BTreeMap::new(),
            staged: Some(staged),
            backup: Some(backup),
        })
    }

    fn staged_mut(&mut self) -> &mut Staged {
        self.staged
            .as_mut()
            .expect("staged output is present before commit")
    }

    fn file(&mut self) -> &mut fs::File {
        self.staged_mut().file()
    }

    fn push(
        &mut self,
        format: ExtractFormat,
        merge_context: bool,
        minimum_depth: u64,
        metric: &SiteMetric,
    ) -> Result<OutputStats> {
        let mut stats = OutputStats::default();
        let mut metric = OutputMetric::from(metric);
        if !merge_context || metric.context == SequenceContext::Chh {
            stats.add(self.write_metric(format, minimum_depth, metric)?)?;
            return Ok(stats);
        }
        if self.chromosome.as_deref() != Some(metric.chromosome.as_str()) {
            stats.add(self.flush_all(format, minimum_depth)?)?;
            self.chromosome = Some(metric.chromosome.clone());
        }
        let position = metric.start;
        metric.apply_merged_span()?;
        let key = (metric.start, metric.end);
        if let Some(existing) = self.pending.get_mut(&key) {
            existing.merge(&metric)?;
            stats.merged_records = 1;
        } else {
            self.pending.insert(key, metric);
        }
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
        Ok(stats)
    }

    fn flush_all(&mut self, format: ExtractFormat, minimum_depth: u64) -> Result<OutputStats> {
        let mut stats = OutputStats::default();
        while let Some((_, metric)) = self.pending.pop_first() {
            stats.add(self.write_metric(format, minimum_depth, metric)?)?;
        }
        Ok(stats)
    }

    fn write_metric(
        &mut self,
        format: ExtractFormat,
        minimum_depth: u64,
        metric: OutputMetric,
    ) -> Result<OutputStats> {
        if metric.depth()? < minimum_depth {
            return Ok(OutputStats::default());
        }
        format.write_metric(self.file(), &metric)?;
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
    fn depth(&self) -> Result<u64> {
        self.methylated
            .checked_add(self.unmethylated)
            .ok_or_else(|| RsomicsError::InvalidInput("methylation depth overflows".into()))
    }

    fn percentage(&self) -> Result<u64> {
        let depth = self.depth()?;
        if depth == 0 {
            Ok(0)
        } else {
            Ok((u128::from(self.methylated) * 100 / u128::from(depth)) as u64)
        }
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

struct Staged {
    path: PathBuf,
    parent: PathBuf,
    temporary: NamedTempFile,
}

impl Staged {
    fn new(path: &Path) -> Result<Self> {
        let parent = parent(path);
        let permissions = match fs::metadata(path) {
            Ok(metadata) if metadata.is_file() => Some(metadata.permissions()),
            Ok(_) => {
                return Err(RsomicsError::InvalidInput(format!(
                    "output {} is not a regular file",
                    path.display()
                )));
            }
            Err(error) if error.kind() == io::ErrorKind::NotFound => None,
            Err(error) => return Err(RsomicsError::Io(error)),
        };
        let mut builder = Builder::new();
        builder.prefix(".rsomics-methyl-");
        #[cfg(unix)]
        if permissions.is_none() {
            builder.permissions(fs::Permissions::from_mode(0o666));
        }
        if let Some(existing) = permissions {
            builder.permissions(existing);
        }
        let temporary = builder
            .tempfile_in(parent)
            .rs_with_context(|| format!("creating temporary output beside {}", path.display()))?;
        Ok(Self {
            path: path.to_owned(),
            parent: parent.to_owned(),
            temporary,
        })
    }

    fn file(&mut self) -> &mut fs::File {
        self.temporary.as_file_mut()
    }

    fn prepare(&mut self) -> Result<()> {
        self.file()
            .flush()
            .rs_with_context(|| format!("flushing output {}", self.path.display()))?;
        self.file()
            .sync_all()
            .rs_with_context(|| format!("syncing output {}", self.path.display()))
    }

    fn commit(self) -> Result<()> {
        self.temporary.persist(&self.path).map_err(|error| {
            RsomicsError::Io(io::Error::new(
                error.error.kind(),
                format!("committing output {}: {}", self.path.display(), error.error),
            ))
        })?;
        #[cfg(unix)]
        fs::File::open(&self.parent)
            .and_then(|directory| directory.sync_all())
            .rs_with_context(|| format!("syncing output directory {}", self.parent.display()))?;
        Ok(())
    }
}

enum Backup {
    Absent,
    Existing(NamedTempFile),
}

impl Backup {
    fn new(path: &Path) -> Result<Self> {
        match fs::metadata(path) {
            Ok(metadata) if !metadata.is_file() => Err(RsomicsError::InvalidInput(format!(
                "output {} is not a regular file",
                path.display()
            ))),
            Ok(_) => Builder::new()
                .prefix(".rsomics-methyl-backup-")
                .make_in(parent(path), |backup| {
                    fs::hard_link(path, backup)?;
                    fs::File::open(backup)
                })
                .map(Self::Existing)
                .rs_with_context(|| format!("backing up output {}", path.display())),
            Err(error) if error.kind() == io::ErrorKind::NotFound => Ok(Self::Absent),
            Err(error) => Err(RsomicsError::Io(error)),
        }
    }

    fn restore(self, path: &Path, cause: RsomicsError) -> RsomicsError {
        let restored = match self {
            Self::Absent => match fs::remove_file(path) {
                Ok(()) => Ok(()),
                Err(error) if error.kind() == io::ErrorKind::NotFound => Ok(()),
                Err(error) => Err(error),
            },
            Self::Existing(backup) => backup
                .persist(path)
                .map(|_| ())
                .map_err(|error| error.error),
        };
        match restored {
            Ok(()) => cause,
            Err(error) => RsomicsError::Io(io::Error::new(
                error.kind(),
                format!(
                    "{cause}; also failed to restore output {}: {error}",
                    path.display()
                ),
            )),
        }
    }
}

impl ExtractFormat {
    fn suffix(self, context: &str) -> String {
        match self {
            Self::Standard => format!("_{context}.bedGraph"),
            Self::Fraction => format!("_{context}.meth.bedGraph"),
            Self::Counts => format!("_{context}.counts.bedGraph"),
            Self::Logit => format!("_{context}.logit.bedGraph"),
            Self::MethylKit => format!("_{context}.methylKit"),
        }
    }

    fn write_header(
        self,
        writer: &mut fs::File,
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

    fn write_metric(self, writer: &mut fs::File, metric: &OutputMetric) -> Result<()> {
        let depth = metric.depth()?;
        let fraction = metric.methylated as f64 / depth as f64;
        match self {
            Self::Standard => writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}",
                metric.chromosome,
                metric.start,
                metric.end,
                metric.percentage()?,
                metric.methylated,
                metric.unmethylated
            ),
            Self::Fraction => writeln!(
                writer,
                "{}\t{}\t{}\t{fraction:.6}",
                metric.chromosome, metric.start, metric.end
            ),
            Self::Counts => writeln!(
                writer,
                "{}\t{}\t{}\t{depth}",
                metric.chromosome, metric.start, metric.end
            ),
            Self::Logit => writeln!(
                writer,
                "{}\t{}\t{}\t{:.6}",
                metric.chromosome,
                metric.start,
                metric.end,
                fraction.ln() - (-fraction).ln_1p()
            ),
            Self::MethylKit => {
                let position = metric.end;
                let strand = match metric.strand {
                    ReferenceStrand::Forward => 'F',
                    ReferenceStrand::Reverse => 'R',
                };
                writeln!(
                    writer,
                    "{0}.{1}\t{0}\t{1}\t{strand}\t{depth}\t{2:6.2}\t{3:6.2}",
                    metric.chromosome,
                    position,
                    fraction * 100.0,
                    (1.0 - fraction) * 100.0
                )
            }
        }
        .map_err(RsomicsError::Io)
    }
}

fn context_path(prefix: &Path, context: &str, format: ExtractFormat) -> PathBuf {
    let mut path = OsString::from(prefix.as_os_str());
    path.push(format.suffix(context));
    PathBuf::from(path)
}

fn parent(path: &Path) -> &Path {
    path.parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."))
}

#[cfg(test)]
mod tests {
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
