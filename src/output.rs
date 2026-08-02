use std::fs;
use std::io::{self, BufWriter, Write};
#[cfg(unix)]
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};

use rsomics_common::{Context, Result, RsomicsError};
use tempfile::{Builder, NamedTempFile};

pub(crate) struct TransactionalOutput {
    path: PathBuf,
    parent: PathBuf,
    temporary: Option<BufWriter<NamedTempFile>>,
    backup: Option<Backup>,
}

const OUTPUT_BUFFER_SIZE: usize = 1024 * 1024;

pub(crate) fn commit_all<T>(
    items: &mut [T],
    mut output: impl FnMut(&mut T) -> &mut TransactionalOutput,
) -> Result<()> {
    for item in &mut *items {
        output(item).prepare()?;
    }
    for index in 0..items.len() {
        if let Err(error) = output(&mut items[index]).commit() {
            let mut cause = error;
            for item in items[..=index].iter_mut().rev() {
                cause = output(item).restore(cause);
            }
            return Err(cause);
        }
    }
    Ok(())
}

impl TransactionalOutput {
    pub(crate) fn new(path: &Path) -> Result<Self> {
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
            temporary: Some(BufWriter::with_capacity(OUTPUT_BUFFER_SIZE, temporary)),
            backup: Some(Backup::new(path)?),
        })
    }

    pub(crate) fn path(&self) -> &Path {
        &self.path
    }

    pub(crate) fn writer(&mut self) -> &mut dyn Write {
        self.temporary
            .as_mut()
            .expect("temporary output is present before commit")
    }

    pub(crate) fn prepare(&mut self) -> Result<()> {
        self.writer()
            .flush()
            .rs_with_context(|| format!("flushing output {}", self.path.display()))?;
        self.temporary
            .as_ref()
            .expect("temporary output is present before commit")
            .get_ref()
            .as_file()
            .sync_data()
            .rs_with_context(|| format!("syncing output {}", self.path.display()))
    }

    pub(crate) fn commit(&mut self) -> Result<()> {
        let temporary = self
            .temporary
            .take()
            .expect("temporary output is present before commit")
            .into_inner()
            .map_err(|error| RsomicsError::Io(error.into_error()))?;
        temporary.persist(&self.path).map_err(|error| {
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

    pub(crate) fn restore(&mut self, cause: RsomicsError) -> RsomicsError {
        let Some(backup) = self.backup.take() else {
            return cause;
        };
        backup.restore(&self.path, cause)
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

fn parent(path: &Path) -> &Path {
    path.parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."))
}
