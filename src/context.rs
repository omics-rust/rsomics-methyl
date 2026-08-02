use rsomics_common::Result;

use crate::reference::IndexedReference;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SequenceContext {
    Cpg,
    Chg,
    Chh,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ReferenceStrand {
    Forward,
    Reverse,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct CytosineContext {
    pub(crate) kind: SequenceContext,
    pub(crate) strand: ReferenceStrand,
    pub(crate) trinucleotide: [u8; 3],
}

pub(crate) fn classify(
    reference: &mut IndexedReference,
    reference_id: usize,
    chromosome: &str,
    length: usize,
    position: usize,
) -> Result<Option<CytosineContext>> {
    let (sequence, offset) = context_window(reference, reference_id, chromosome, length, position)?;
    let Some((kind, strand)) = classify_sequence(sequence, offset) else {
        return Ok(None);
    };
    let trinucleotide = match strand {
        ReferenceStrand::Forward => [
            b'C',
            normalized(sequence.get(offset + 1)),
            normalized(sequence.get(offset + 2)),
        ],
        ReferenceStrand::Reverse => [
            b'C',
            complement(offset.checked_sub(1).and_then(|index| sequence.get(index))),
            complement(offset.checked_sub(2).and_then(|index| sequence.get(index))),
        ],
    };
    Ok(Some(CytosineContext {
        kind,
        strand,
        trinucleotide,
    }))
}

#[inline]
pub(crate) fn classify_call(
    reference: &mut IndexedReference,
    reference_id: usize,
    chromosome: &str,
    length: usize,
    position: usize,
) -> Result<Option<(SequenceContext, ReferenceStrand)>> {
    let (sequence, offset) = context_window(reference, reference_id, chromosome, length, position)?;
    Ok(classify_sequence(sequence, offset))
}

#[inline]
fn context_window<'a>(
    reference: &'a mut IndexedReference,
    reference_id: usize,
    chromosome: &str,
    length: usize,
    position: usize,
) -> Result<(&'a [u8], usize)> {
    if position >= length {
        return Err(reference.error(format!(
            "{chromosome}:{position} is outside reference length {length}"
        )));
    }
    let start = position.saturating_sub(2);
    let end = position.saturating_add(3).min(length);
    let sequence = reference.sequence_by_id(reference_id, chromosome, start..end)?;
    Ok((sequence, position - start))
}

fn classify_sequence(sequence: &[u8], offset: usize) -> Option<(SequenceContext, ReferenceStrand)> {
    let (kind, strand) = match sequence[offset].to_ascii_uppercase() {
        b'C' => (
            if sequence
                .get(offset + 1)
                .is_some_and(|base| base.eq_ignore_ascii_case(&b'G'))
            {
                SequenceContext::Cpg
            } else if sequence
                .get(offset + 2)
                .is_some_and(|base| base.eq_ignore_ascii_case(&b'G'))
            {
                SequenceContext::Chg
            } else {
                SequenceContext::Chh
            },
            ReferenceStrand::Forward,
        ),
        b'G' => (
            if offset >= 1 && sequence[offset - 1].eq_ignore_ascii_case(&b'C') {
                SequenceContext::Cpg
            } else if offset >= 2 && sequence[offset - 2].eq_ignore_ascii_case(&b'C') {
                SequenceContext::Chg
            } else {
                SequenceContext::Chh
            },
            ReferenceStrand::Reverse,
        ),
        _ => return None,
    };
    Some((kind, strand))
}

fn normalized(base: Option<&u8>) -> u8 {
    base.map(u8::to_ascii_uppercase)
        .filter(|base| matches!(base, b'A' | b'C' | b'G' | b'T'))
        .unwrap_or(b'N')
}

fn complement(base: Option<&u8>) -> u8 {
    match normalized(base) {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn classifies_both_strands_and_boundaries() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("reference.fa");
        std::fs::write(&path, b">chr1\nCGCAGCT\n").unwrap();
        std::fs::write(
            directory.path().join("reference.fa.fai"),
            b"chr1\t7\t6\t7\t8\n",
        )
        .unwrap();
        let mut reference = IndexedReference::open(&path).unwrap();
        assert_eq!(
            classify(&mut reference, 0, "chr1", 7, 0)
                .unwrap()
                .unwrap()
                .kind,
            SequenceContext::Cpg
        );
        assert_eq!(
            classify(&mut reference, 0, "chr1", 7, 0)
                .unwrap()
                .unwrap()
                .trinucleotide,
            *b"CGC"
        );
        assert_eq!(
            classify(&mut reference, 0, "chr1", 7, 1)
                .unwrap()
                .unwrap()
                .strand,
            ReferenceStrand::Reverse
        );
        assert_eq!(
            classify(&mut reference, 0, "chr1", 7, 1)
                .unwrap()
                .unwrap()
                .trinucleotide,
            *b"CGN"
        );
        assert_eq!(
            classify(&mut reference, 0, "chr1", 7, 2)
                .unwrap()
                .unwrap()
                .kind,
            SequenceContext::Chg
        );
        assert_eq!(
            classify(&mut reference, 0, "chr1", 7, 5)
                .unwrap()
                .unwrap()
                .kind,
            SequenceContext::Chh
        );
        assert_eq!(classify(&mut reference, 0, "chr1", 7, 6).unwrap(), None);
    }
}
