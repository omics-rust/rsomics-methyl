use std::collections::HashMap;
use std::fs::File;
use std::ops::Range;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use noodles::core::{Position, Region};
use noodles::fasta;
use noodles::sam;
use rsomics_common::{Result, RsomicsError};

const CHUNK_SIZE: usize = 1024 * 1024;

pub(crate) struct IndexedReference {
    reader: fasta::io::IndexedReader<fasta::io::BufReader<File>>,
    path: PathBuf,
    lengths: HashMap<Vec<u8>, u64>,
    chromosome: Vec<u8>,
    sequence_start: usize,
    sequence: Vec<u8>,
    #[cfg(test)]
    fetch_count: usize,
}

impl IndexedReference {
    pub(crate) fn open(path: &Path) -> Result<Self> {
        let reader = fasta::io::indexed_reader::Builder::default()
            .build_from_path(path)
            .map_err(|error| Self::path_error(path, error))?;
        let mut lengths = HashMap::new();
        for record in reader.index().as_ref() {
            if lengths
                .insert(record.name().to_vec(), record.length())
                .is_some()
            {
                return Err(Self::path_error(
                    path,
                    format!(
                        "duplicate reference {} in FASTA index",
                        String::from_utf8_lossy(record.name())
                    ),
                ));
            }
        }
        Ok(Self {
            reader,
            path: path.to_path_buf(),
            lengths,
            chromosome: Vec::new(),
            sequence_start: 0,
            sequence: Vec::new(),
            #[cfg(test)]
            fetch_count: 0,
        })
    }

    pub(crate) fn length(&self, chromosome: &str) -> Result<u64> {
        self.lengths
            .get(chromosome.as_bytes())
            .copied()
            .ok_or_else(|| self.error(format!("unknown chromosome {chromosome}")))
    }

    pub(crate) fn validate_header(&self, header: &sam::Header) -> Result<Vec<ReferenceSequence>> {
        header
            .reference_sequences()
            .iter()
            .map(|(name, sequence)| {
                let name = String::from_utf8(name.to_vec()).map_err(|_| {
                    self.error("alignment header contains a non-UTF-8 reference name")
                })?;
                let alignment_length = u64::try_from(usize::from(sequence.length()))
                    .map_err(|error| self.error(error))?;
                let reference_length = self.length(&name)?;
                if alignment_length != reference_length {
                    return Err(self.error(format!(
                        "header length {alignment_length} for {name} differs from indexed reference length {reference_length}"
                    )));
                }
                Ok(ReferenceSequence {
                    name: Arc::from(name),
                    length: reference_length,
                })
            })
            .collect()
    }

    pub(crate) fn sequence(&mut self, chromosome: &str, range: Range<usize>) -> Result<&[u8]> {
        let needs_fetch = self.chromosome != chromosome.as_bytes()
            || range.start < self.sequence_start
            || range.end > self.sequence_start + self.sequence.len();
        if needs_fetch {
            let length =
                usize::try_from(self.length(chromosome)?).map_err(|error| self.error(error))?;
            if range.start >= range.end || range.end > length {
                return Err(self.error("requested reference range is invalid"));
            }
            let start = range.start;
            let end = start
                .checked_add(CHUNK_SIZE)
                .map_or(length, |value| value.min(length))
                .max(range.end);
            let interval_start =
                Position::try_from(start + 1).map_err(|error| self.error(error))?;
            let interval_end = Position::try_from(end).map_err(|error| self.error(error))?;
            let record = self
                .reader
                .query(&Region::new(
                    chromosome.as_bytes().to_vec(),
                    interval_start..=interval_end,
                ))
                .map_err(|error| self.error(error))?;
            self.chromosome.clear();
            self.chromosome.extend_from_slice(chromosome.as_bytes());
            self.sequence_start = start;
            self.sequence.clear();
            self.sequence.extend_from_slice(record.sequence().as_ref());
            #[cfg(test)]
            {
                self.fetch_count += 1;
            }
        }
        let start = range
            .start
            .checked_sub(self.sequence_start)
            .ok_or_else(|| self.error("reference cache start is invalid"))?;
        let end = range
            .end
            .checked_sub(self.sequence_start)
            .ok_or_else(|| self.error("reference cache end is invalid"))?;
        self.sequence
            .get(start..end)
            .ok_or_else(|| self.error("reference cache range is invalid"))
    }

    pub(crate) fn error(&self, error: impl std::fmt::Display) -> RsomicsError {
        Self::path_error(&self.path, error)
    }

    fn path_error(path: &Path, error: impl std::fmt::Display) -> RsomicsError {
        RsomicsError::InvalidInput(format!(
            "reading indexed reference {}: {error}",
            path.display()
        ))
    }
}

#[derive(Clone, Debug)]
pub(crate) struct ReferenceSequence {
    pub(crate) name: Arc<str>,
    pub(crate) length: u64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_duplicate_fasta_index_names() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("reference.fa");
        std::fs::write(&path, b">chr1\nA\n").unwrap();
        std::fs::write(
            directory.path().join("reference.fa.fai"),
            b"chr1\t1\t6\t1\t2\nchr1\t1\t6\t1\t2\n",
        )
        .unwrap();

        assert!(IndexedReference::open(&path).is_err());
    }

    #[test]
    fn advances_the_reference_cache_across_chunk_boundaries() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("reference.fa");
        let length = CHUNK_SIZE + 8;
        let sequence = vec![b'A'; length];
        let mut fasta = b">chr1\n".to_vec();
        for line in sequence.chunks(64) {
            fasta.extend_from_slice(line);
            fasta.push(b'\n');
        }
        std::fs::write(&path, fasta).unwrap();
        std::fs::write(
            directory.path().join("reference.fa.fai"),
            format!("chr1\t{length}\t6\t64\t65\n"),
        )
        .unwrap();

        let mut reference = IndexedReference::open(&path).unwrap();
        assert_eq!(reference.sequence("chr1", 0..5).unwrap(), b"AAAAA");
        assert_eq!(
            reference
                .sequence("chr1", CHUNK_SIZE - 2..CHUNK_SIZE + 3)
                .unwrap(),
            b"AAAAA"
        );
        assert_eq!(
            reference
                .sequence("chr1", CHUNK_SIZE - 1..CHUNK_SIZE + 4)
                .unwrap(),
            b"AAAAA"
        );
        assert_eq!(reference.fetch_count, 2);
    }
}
