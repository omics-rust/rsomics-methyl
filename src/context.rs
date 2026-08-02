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
}

pub(crate) fn classify(
    reference: &mut IndexedReference,
    chromosome: &str,
    position: usize,
) -> Result<Option<CytosineContext>> {
    let length =
        usize::try_from(reference.length(chromosome)?).map_err(|error| reference.error(error))?;
    if position >= length {
        return Err(reference.error(format!(
            "{chromosome}:{position} is outside reference length {length}"
        )));
    }
    let start = position.saturating_sub(2);
    let end = position.saturating_add(3).min(length);
    let offset = position - start;
    let sequence = reference.sequence(chromosome, start..end)?;
    let context = match sequence[offset].to_ascii_uppercase() {
        b'C' => CytosineContext {
            kind: if sequence
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
            strand: ReferenceStrand::Forward,
        },
        b'G' => CytosineContext {
            kind: if offset >= 1 && sequence[offset - 1].eq_ignore_ascii_case(&b'C') {
                SequenceContext::Cpg
            } else if offset >= 2 && sequence[offset - 2].eq_ignore_ascii_case(&b'C') {
                SequenceContext::Chg
            } else {
                SequenceContext::Chh
            },
            strand: ReferenceStrand::Reverse,
        },
        _ => return Ok(None),
    };
    Ok(Some(context))
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
            classify(&mut reference, "chr1", 0).unwrap().unwrap().kind,
            SequenceContext::Cpg
        );
        assert_eq!(
            classify(&mut reference, "chr1", 1).unwrap().unwrap().strand,
            ReferenceStrand::Reverse
        );
        assert_eq!(
            classify(&mut reference, "chr1", 2).unwrap().unwrap().kind,
            SequenceContext::Chg
        );
        assert_eq!(
            classify(&mut reference, "chr1", 5).unwrap().unwrap().kind,
            SequenceContext::Chh
        );
        assert_eq!(classify(&mut reference, "chr1", 6).unwrap(), None);
    }
}
