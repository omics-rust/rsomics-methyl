use rsomics_bamio::raw::RawRecord;
use rsomics_common::Result;

use crate::alignment::invalid_record;
use crate::context::{ReferenceStrand, SequenceContext, classify_call, is_cpg_call};
use crate::reference::{IndexedReference, ReferenceSequence};
use crate::strand::{BisulfiteStrand, bisulfite_strand};

const READ_2: u16 = 0x80;

pub(crate) struct AlignmentCaller {
    reference: IndexedReference,
    references: Vec<ReferenceSequence>,
    minimum_base_quality: u8,
}

pub(crate) struct AlignmentLocation {
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
        emit: impl FnMut(MethylationCall) -> Result<()>,
    ) -> Result<AlignmentLocation> {
        self.visit_contexts::<false>(record, emit)
    }

    pub(crate) fn visit_cpg(
        &mut self,
        record: &RawRecord,
        emit: impl FnMut(MethylationCall) -> Result<()>,
    ) -> Result<AlignmentLocation> {
        self.visit_contexts::<true>(record, emit)
    }

    fn visit_contexts<const CPG_ONLY: bool>(
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
        let reference_length =
            usize::try_from(reference.length).map_err(|error| invalid_record(record, error))?;
        let strand = bisulfite_strand(record)?;
        let read = read_number(record);
        let sequence_length = record.sequence_len();
        let sequence = record.seq_bytes_packed();
        let quality_scores = record.quality_scores();
        let mut query_position = 0usize;
        let mut reference_position = start;
        for (kind, raw_length) in cigar_operations(record)? {
            if raw_length == 0 {
                return Err(invalid_record(
                    record,
                    "CIGAR contains a zero-length operation",
                ));
            }
            let length =
                usize::try_from(raw_length).map_err(|error| invalid_record(record, error))?;
            match kind {
                0 | 7 | 8 => {
                    for _ in 0..length {
                        if query_position >= sequence_length {
                            return Err(invalid_record(
                                record,
                                "CIGAR consumes beyond the sequence",
                            ));
                        }
                        let quality = quality_scores
                            .get(query_position)
                            .copied()
                            .unwrap_or(u8::MAX);
                        if quality >= self.minimum_base_quality
                            && let Some(methylated) =
                                methylation_state(strand, packed_base(sequence, query_position))
                            && let Some(context) = Self::call_context::<CPG_ONLY>(
                                &mut self.reference,
                                reference_id,
                                reference,
                                reference_length,
                                reference_position,
                                strand,
                            )?
                        {
                            emit(MethylationCall {
                                context,
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
        if query_position != sequence_length {
            return Err(invalid_record(
                record,
                format!(
                    "CIGAR consumes {query_position} query bases instead of {}",
                    sequence_length
                ),
            ));
        }
        if reference_position > reference_length {
            return Err(invalid_record(record, "CIGAR extends beyond the reference"));
        }
        Ok(AlignmentLocation {
            reference_id,
            start: u64::try_from(start).map_err(|error| invalid_record(record, error))?,
            end: u64::try_from(reference_position)
                .map_err(|error| invalid_record(record, error))?,
            strand,
        })
    }

    #[inline]
    fn call_context<const CPG_ONLY: bool>(
        reference_cache: &mut IndexedReference,
        reference_id: usize,
        reference: &ReferenceSequence,
        reference_length: usize,
        reference_position: usize,
        strand: BisulfiteStrand,
    ) -> Result<Option<SequenceContext>> {
        if CPG_ONLY {
            return is_cpg_call(
                reference_cache,
                reference_id,
                &reference.name,
                reference_length,
                reference_position,
                strand.is_top(),
            )
            .map(|is_cpg| is_cpg.then_some(SequenceContext::Cpg));
        }
        Ok(classify_call(
            reference_cache,
            reference_id,
            &reference.name,
            reference_length,
            reference_position,
        )?
        .and_then(|(context, reference_strand)| {
            (strand.is_top() == (reference_strand == ReferenceStrand::Forward)).then_some(context)
        }))
    }

    pub(crate) fn reference_name(&self, reference_id: usize) -> &str {
        self.references
            .get(reference_id)
            .map(|reference| reference.name.as_ref())
            .expect("alignment location has a validated reference ID")
    }
}

enum CigarOperations {
    One(Option<(u8, u32)>),
    Many(std::vec::IntoIter<(u8, u32)>),
}

impl Iterator for CigarOperations {
    type Item = (u8, u32);

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            Self::One(operation) => operation.take(),
            Self::Many(operations) => operations.next(),
        }
    }
}

fn cigar_operations(record: &RawRecord) -> Result<CigarOperations> {
    let mut operations = record.cigar_ops();
    let first = operations.next();
    let second = operations.next();
    if let (Some(operation), None) = (first, second) {
        return Ok(CigarOperations::One(Some(operation)));
    }
    Ok(CigarOperations::Many(record.decoded_cigar()?.into_iter()))
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

#[inline]
fn packed_base(sequence: &[u8], position: usize) -> u8 {
    let byte = sequence[position / 2];
    if position.is_multiple_of(2) {
        byte >> 4
    } else {
        byte & 0x0f
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
