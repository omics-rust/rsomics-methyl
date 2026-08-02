use std::ffi::OsString;
use std::fs;
use std::io::{self, Write};
#[cfg(unix)]
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};

use rsomics_common::{Context, Result, RsomicsError, reject_output_alias};
use rsomics_methyl::SequenceContext;
use rsomics_methyl::extract::{ExtractOptions, ExtractStats, SiteMetric, extract};
use tempfile::{Builder, NamedTempFile};

pub struct ExtractOutputResult {
    pub stats: ExtractStats,
    pub outputs: Vec<PathBuf>,
}

pub fn extract_to_standard_outputs(
    input: &Path,
    reference: &Path,
    prefix: &Path,
    options: ExtractOptions,
) -> Result<ExtractOutputResult> {
    let mut outputs = ContextOutputs::new(prefix, input, reference, &options)?;
    let paths = outputs.paths();
    let stats = extract(input, reference, options, |metric| outputs.write(&metric))?;
    outputs.commit()?;
    Ok(ExtractOutputResult {
        stats,
        outputs: paths,
    })
}

struct ContextOutputs {
    entries: Vec<OutputEntry>,
}

impl ContextOutputs {
    fn new(
        prefix: &Path,
        input: &Path,
        reference: &Path,
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
            let path = context_path(prefix, label);
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
        let mut outputs = Self { entries };
        for entry in &mut outputs.entries {
            let label = entry.label;
            writeln!(
                entry.file(),
                "track type=\"bedGraph\" description=\"{} {} methylation levels\"",
                prefix.display(),
                label
            )
            .map_err(RsomicsError::Io)?;
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
        let entry = self
            .entries
            .iter_mut()
            .find(|entry| entry.context == metric.context())
            .ok_or_else(|| RsomicsError::ConfigError("methylation context has no output".into()))?;
        writeln!(
            entry.file(),
            "{}\t{}\t{}\t{}\t{}\t{}",
            metric.chromosome(),
            metric.start(),
            metric.end(),
            metric.percentage(),
            metric.methylated(),
            metric.unmethylated()
        )
        .map_err(RsomicsError::Io)
    }

    fn commit(mut self) -> Result<()> {
        for entry in &mut self.entries {
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
        Ok(())
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

fn context_path(prefix: &Path, context: &str) -> PathBuf {
    let mut path = OsString::from(prefix.as_os_str());
    path.push(format!("_{context}.bedGraph"));
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
        let first = context_path(&prefix, "CpG");
        fs::write(&first, b"old\n").unwrap();
        let mut outputs = ContextOutputs::new(
            &prefix,
            Path::new("input.bam"),
            Path::new("reference.fa"),
            &ExtractOptions {
                chg: true,
                ..ExtractOptions::default()
            },
        )
        .unwrap();
        for entry in &mut outputs.entries {
            entry.file().write_all(b"new\n").unwrap();
        }
        let second = context_path(&prefix, "CHG");
        fs::create_dir(&second).unwrap();

        assert!(outputs.commit().is_err());
        assert_eq!(fs::read(first).unwrap(), b"old\n");
        assert!(second.is_dir());
    }
}
