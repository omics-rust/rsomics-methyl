use rsomics_bamio::raw::RawRecord;
use rsomics_common::Result;

use crate::alignment::invalid_record;
use crate::context::{ReferenceStrand, SequenceContext, classify};
use crate::reference::{IndexedReference, ReferenceSequence};
use crate::strand::{BisulfiteStrand, bisulfite_strand};

const READ_2: u16 = 0x80;

pub(crate) struct AlignmentCaller {
    reference: IndexedReference,
    references: Vec<ReferenceSequence>,
    minimum_base_quality: u8,
}

pub(crate) struct AlignmentLocation {
    pub(crate) chromosome: String,
    pub(crate) reference_id: usize,
    pub(crate) start: u64,
    pub(crate) end: u64,
    pub(crate) strand: BisulfiteStrand,
}

#[derive(Clone, Copy)]
pub(crate) struct MethylationCall {
    pub(crate) context: SequenceContext,
    pub(crate) reference_id: usize,
    pub(crate) reference_position: u64,
    pub(crate) query_position: u64,
    pub(crate) strand: BisulfiteStrand,
    pub(crate) read: u8,
    pub(crate) methylated: bool,
}

impl AlignmentCaller {
    pub(crate) fn new(
        reference: IndexedReference,
        references: Vec<ReferenceSequence>,
        minimum_base_quality: u8,
    ) -> Self {
        Self {
            reference,
            references,
            minimum_base_quality,
        }
    }

    pub(crate) fn visit(
        &mut self,
        record: &RawRecord,
        mut emit: impl FnMut(MethylationCall) -> Result<()>,
    ) -> Result<AlignmentLocation> {
        let reference_id = usize::try_from(record.reference_sequence_id())
            .map_err(|error| invalid_record(record, error))?;
        let reference = self.references.get(reference_id).ok_or_else(|| {
            invalid_record(record, format!("reference ID {reference_id} is absent"))
        })?;
        let start = usize::try_from(record.alignment_start())
            .map_err(|error| invalid_record(record, error))?;
        let strand = bisulfite_strand(record)?;
        let read = read_number(record);
        let mut query_position = 0usize;
        let mut reference_position = start;
        for (kind, raw_length) in record.decoded_cigar()? {
            let length =
                usize::try_from(raw_length).map_err(|error| invalid_record(record, error))?;
            match kind {
                0 | 7 | 8 => {
                    for _ in 0..length {
                        if query_position >= record.sequence_len() {
                            return Err(invalid_record(
                                record,
                                "CIGAR consumes beyond the sequence",
                            ));
                        }
                        let quality = record
                            .quality_scores()
                            .get(query_position)
                            .copied()
                            .unwrap_or(u8::MAX);
                        if quality >= self.minimum_base_quality
                            && let Some(context) =
                                classify(&mut self.reference, &reference.name, reference_position)?
                            && strand.is_top() == (context.strand == ReferenceStrand::Forward)
                            && let Some(methylated) =
                                methylation_state(strand, record.seq_nibble(query_position))
                        {
                            emit(MethylationCall {
                                context: context.kind,
                                reference_id,
                                reference_position: u64::try_from(reference_position)
                                    .map_err(|error| invalid_record(record, error))?,
                                query_position: u64::try_from(query_position)
                                    .map_err(|error| invalid_record(record, error))?,
                                strand,
                                read,
                                methylated,
                            })?;
                        }
                        query_position = checked_advance(query_position, 1, record)?;
                        reference_position = checked_advance(reference_position, 1, record)?;
                    }
                }
                1 | 4 => {
                    query_position = checked_advance(query_position, length, record)?;
                }
                2 | 3 => {
                    reference_position = checked_advance(reference_position, length, record)?;
                }
                5 | 6 => {}
                _ => {
                    return Err(invalid_record(
                        record,
                        format!("unsupported CIGAR operation {kind}"),
                    ));
                }
            }
        }
        if query_position != record.sequence_len() {
            return Err(invalid_record(
                record,
                format!(
                    "CIGAR consumes {query_position} query bases instead of {}",
                    record.sequence_len()
                ),
            ));
        }
        let reference_length =
            usize::try_from(reference.length).map_err(|error| invalid_record(record, error))?;
        if reference_position > reference_length {
            return Err(invalid_record(record, "CIGAR extends beyond the reference"));
        }
        Ok(AlignmentLocation {
            chromosome: reference.name.clone(),
            reference_id,
            start: u64::try_from(start).map_err(|error| invalid_record(record, error))?,
            end: u64::try_from(reference_position)
                .map_err(|error| invalid_record(record, error))?,
            strand,
        })
    }
}

pub(crate) fn read_number(record: &RawRecord) -> u8 {
    if record.flags() & READ_2 == 0 { 1 } else { 2 }
}

fn methylation_state(strand: BisulfiteStrand, base: u8) -> Option<bool> {
    match (strand.is_top(), base) {
        (true, 2) | (false, 4) => Some(true),
        (true, 8) | (false, 1) => Some(false),
        _ => None,
    }
}

fn checked_advance(position: usize, length: usize, record: &RawRecord) -> Result<usize> {
    position
        .checked_add(length)
        .ok_or_else(|| invalid_record(record, "CIGAR coordinate overflows"))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn raw(flags: u16, converted: &[u8], sequence: &[u8]) -> RawRecord {
        let mut payload = Vec::new();
        payload.extend_from_slice(&0i32.to_le_bytes());
        payload.extend_from_slice(&0i32.to_le_bytes());
        payload.push(5);
        payload.push(60);
        payload.extend_from_slice(&0u16.to_le_bytes());
        payload.extend_from_slice(&1u16.to_le_bytes());
        payload.extend_from_slice(&flags.to_le_bytes());
        payload.extend_from_slice(&(sequence.len() as u32).to_le_bytes());
        payload.extend_from_slice(&(-1i32).to_le_bytes());
        payload.extend_from_slice(&(-1i32).to_le_bytes());
        payload.extend_from_slice(&0i32.to_le_bytes());
        payload.extend_from_slice(b"read\0");
        payload.extend_from_slice(&((sequence.len() as u32) << 4).to_le_bytes());
        for pair in sequence.chunks(2) {
            let high = base_code(pair[0]);
            let low = pair.get(1).copied().map_or(0, base_code);
            payload.push(high << 4 | low);
        }
        payload.extend(std::iter::repeat_n(40, sequence.len()));
        let mut record = RawRecord::try_from(payload).unwrap();
        let mut value = converted.to_vec();
        value.push(0);
        record.append_aux(*b"XG", b'Z', &value).unwrap();
        record
    }

    fn base_code(base: u8) -> u8 {
        match base {
            b'A' => 1,
            b'C' => 2,
            b'G' => 4,
            b'T' => 8,
            _ => 15,
        }
    }

    fn caller() -> (tempfile::TempDir, AlignmentCaller) {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("reference.fa");
        std::fs::write(&path, b">chr1\nCG\n").unwrap();
        std::fs::write(
            directory.path().join("reference.fa.fai"),
            b"chr1\t2\t6\t2\t3\n",
        )
        .unwrap();
        let caller = AlignmentCaller::new(
            IndexedReference::open(&path).unwrap(),
            vec![ReferenceSequence {
                name: "chr1".into(),
                length: 2,
            }],
            5,
        );
        (directory, caller)
    }

    #[test]
    fn reports_nondirectional_read_two_strands() {
        let (_directory, mut caller) = caller();
        let mut calls = Vec::new();
        caller
            .visit(&raw(0x81, b"CT", b"CA"), |call| {
                calls.push(call);
                Ok(())
            })
            .unwrap();
        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].strand, BisulfiteStrand::Ctot);
        assert_eq!(calls[0].read, 2);
        assert_eq!(calls[0].query_position, 0);
        assert!(calls[0].methylated);

        calls.clear();
        caller
            .visit(&raw(0x91, b"GA", b"AG"), |call| {
                calls.push(call);
                Ok(())
            })
            .unwrap();
        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].strand, BisulfiteStrand::Ctob);
        assert_eq!(calls[0].read, 2);
        assert_eq!(calls[0].query_position, 1);
        assert!(calls[0].methylated);
    }
}
