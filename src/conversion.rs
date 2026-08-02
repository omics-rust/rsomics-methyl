use std::path::Path;

use rsomics_bamio::raw::RawRecord;
use rsomics_common::{Result, RsomicsError};

use crate::calling::AlignmentCaller;
use crate::context::SequenceContext;
use crate::reference::{IndexedReference, ReferenceSequence};

pub(crate) struct ConversionFilter {
    caller: AlignmentCaller,
    minimum: f64,
}

impl ConversionFilter {
    pub(crate) fn new(
        reference: &Path,
        references: Vec<ReferenceSequence>,
        minimum_base_quality: u8,
        minimum: f64,
    ) -> Result<Option<Self>> {
        validate_conversion_efficiency(minimum)?;
        if minimum == 0.0 {
            return Ok(None);
        }
        Ok(Some(Self {
            caller: AlignmentCaller::new(
                IndexedReference::open(reference)?,
                references,
                minimum_base_quality,
            ),
            minimum,
        }))
    }

    pub(crate) fn passes(&mut self, record: &RawRecord) -> Result<bool> {
        let mut methylated = 0u64;
        let mut unmethylated = 0u64;
        self.caller.visit(record, |call| {
            if call.context != SequenceContext::Cpg {
                if call.methylated {
                    methylated = increment(methylated)?;
                } else {
                    unmethylated = increment(unmethylated)?;
                }
            }
            Ok(())
        })?;
        let informative = methylated
            .checked_add(unmethylated)
            .ok_or_else(|| RsomicsError::InvalidInput("conversion count overflows".into()))?;
        Ok(informative == 0 || unmethylated as f64 / informative as f64 >= self.minimum)
    }
}

pub(crate) fn validate_conversion_efficiency(value: f64) -> Result<()> {
    if !value.is_finite() || !(0.0..=1.0).contains(&value) {
        return Err(RsomicsError::ConfigError(
            "minimum conversion efficiency must be between 0 and 1".into(),
        ));
    }
    Ok(())
}

fn increment(value: u64) -> Result<u64> {
    value
        .checked_add(1)
        .ok_or_else(|| RsomicsError::InvalidInput("conversion count overflows".into()))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn raw(qualities: &[u8]) -> RawRecord {
        let sequence = b"CAAATTAA";
        let cigar = [(0u8, 3u32), (1, 1), (0, 4)];
        let mut payload = Vec::new();
        payload.extend_from_slice(&0i32.to_le_bytes());
        payload.extend_from_slice(&0i32.to_le_bytes());
        payload.push(5);
        payload.push(60);
        payload.extend_from_slice(&0u16.to_le_bytes());
        payload.extend_from_slice(&(cigar.len() as u16).to_le_bytes());
        payload.extend_from_slice(&0u16.to_le_bytes());
        payload.extend_from_slice(&(sequence.len() as u32).to_le_bytes());
        payload.extend_from_slice(&(-1i32).to_le_bytes());
        payload.extend_from_slice(&(-1i32).to_le_bytes());
        payload.extend_from_slice(&0i32.to_le_bytes());
        payload.extend_from_slice(b"read\0");
        for (kind, length) in cigar {
            payload.extend_from_slice(&((length << 4) | u32::from(kind)).to_le_bytes());
        }
        for pair in sequence.chunks(2) {
            payload.push(base(pair[0]) << 4 | pair.get(1).copied().map_or(0, base));
        }
        payload.extend_from_slice(qualities);
        let mut record = RawRecord::try_from(payload).unwrap();
        record.append_aux(*b"XG", b'Z', b"CT\0").unwrap();
        record
    }

    fn base(value: u8) -> u8 {
        match value {
            b'A' => 1,
            b'C' => 2,
            b'G' => 4,
            b'T' => 8,
            _ => 15,
        }
    }

    fn conversion_filter(minimum: f64) -> (tempfile::TempDir, ConversionFilter) {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("reference.fa");
        std::fs::write(&path, b">chr1\nCAATCAA\n").unwrap();
        std::fs::write(
            directory.path().join("reference.fa.fai"),
            b"chr1\t7\t6\t7\t8\n",
        )
        .unwrap();
        let filter = ConversionFilter::new(
            &path,
            vec![ReferenceSequence {
                name: "chr1".into(),
                length: 7,
            }],
            5,
            minimum,
        )
        .unwrap()
        .unwrap();
        (directory, filter)
    }

    #[test]
    fn follows_the_full_cigar_and_ignores_low_quality_bases() {
        let (_directory, mut filter) = conversion_filter(0.75);
        assert!(!filter.passes(&raw(&[40; 8])).unwrap());

        let (_directory, mut filter) = conversion_filter(1.0);
        let mut qualities = [40; 8];
        qualities[0] = 0;
        assert!(filter.passes(&raw(&qualities)).unwrap());
    }

    #[test]
    fn rejects_invalid_thresholds() {
        for value in [f64::NAN, f64::INFINITY, -0.1, 1.1] {
            assert!(validate_conversion_efficiency(value).is_err());
        }
    }
}
