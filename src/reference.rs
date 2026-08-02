use std::fs::File;
use std::ops::Range;
use std::path::{Path, PathBuf};

use noodles::core::{Position, Region};
use noodles::fasta;
use noodles::sam;
use rsomics_common::{Result, RsomicsError};

const CHUNK_SIZE: usize = 1024 * 1024;

pub(crate) struct IndexedReference {
    reader: fasta::io::IndexedReader<fasta::io::BufReader<File>>,
    path: PathBuf,
    chromosome: Vec<u8>,
    sequence_start: usize,
    sequence: Vec<u8>,
}

impl IndexedReference {
    pub(crate) fn open(path: &Path) -> Result<Self> {
        let reader = fasta::io::indexed_reader::Builder::default()
            .build_from_path(path)
            .map_err(|error| Self::path_error(path, error))?;
        Ok(Self {
            reader,
            path: path.to_path_buf(),
            chromosome: Vec::new(),
            sequence_start: 0,
            sequence: Vec::new(),
        })
    }

    pub(crate) fn length(&self, chromosome: &str) -> Result<u64> {
        self.reader
            .index()
            .as_ref()
            .iter()
            .find(|record| record.name() == chromosome.as_bytes())
            .map(fasta::fai::Record::length)
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
                    name,
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
            let start = range.start / CHUNK_SIZE * CHUNK_SIZE;
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
    pub(crate) name: String,
    pub(crate) length: u64,
}
