use std::io;

use noodles::core::{Position, Region};
use noodles::sam;
use rsomics_bamio::IndexedAlignmentReader;
use rsomics_common::{Result, RsomicsError};

use crate::reference::ReferenceSequence;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct ReferenceRange {
    pub(crate) reference_id: usize,
    pub(crate) start: u64,
    pub(crate) end: u64,
}

impl ReferenceRange {
    pub(crate) fn contains(&self, reference_id: usize, position: u64) -> bool {
        reference_id == self.reference_id && position >= self.start && position < self.end
    }
}

pub(crate) struct RegionSelection {
    pub(crate) query: Region,
    pub(crate) range: ReferenceRange,
}

pub(crate) type AlignmentRecordResult = io::Result<Box<dyn sam::alignment::Record>>;

pub(crate) fn alignment_records<'r, 'h: 'r>(
    reader: &'r mut IndexedAlignmentReader,
    header: &'h sam::Header,
    selection: Option<&RegionSelection>,
) -> io::Result<Box<dyn Iterator<Item = AlignmentRecordResult> + 'r>> {
    match selection {
        Some(selection) => reader
            .query(header, &selection.query)
            .map(|records| Box::new(records) as Box<dyn Iterator<Item = AlignmentRecordResult>>),
        None => Ok(Box::new(reader.records(header))),
    }
}

pub(crate) fn resolve_region(
    references: &[ReferenceSequence],
    region: &Region,
) -> Result<RegionSelection> {
    let reference_id = references
        .iter()
        .position(|reference| reference.name.as_bytes() == region.name())
        .ok_or_else(|| {
            RsomicsError::InvalidInput(format!(
                "region reference {} is absent from the alignment header",
                region.name()
            ))
        })?;
    let reference = &references[reference_id];
    let interval = region.interval();
    let start = interval
        .start()
        .map(|position| u64::try_from(position.get()))
        .transpose()
        .map_err(invalid_coordinate)?
        .map_or(0, |position| position - 1);
    let end = interval
        .end()
        .map(|position| u64::try_from(position.get()))
        .transpose()
        .map_err(invalid_coordinate)?
        .unwrap_or(reference.length)
        .min(reference.length);
    if start >= end {
        return Err(RsomicsError::InvalidInput(format!(
            "region {region} is outside reference {} of length {}",
            reference.name, reference.length
        )));
    }
    let query_start = position(start + 1)?;
    let query_end = position(end)?;
    Ok(RegionSelection {
        query: Region::new(region.name().to_vec(), query_start..=query_end),
        range: ReferenceRange {
            reference_id,
            start,
            end,
        },
    })
}

fn position(value: u64) -> Result<Position> {
    let value = usize::try_from(value).map_err(invalid_coordinate)?;
    Position::try_from(value).map_err(invalid_coordinate)
}

fn invalid_coordinate(error: impl std::fmt::Display) -> RsomicsError {
    RsomicsError::InvalidInput(format!("region coordinate is invalid: {error}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn references() -> Vec<ReferenceSequence> {
        vec![ReferenceSequence {
            name: "chr1".into(),
            length: 10,
        }]
    }

    #[test]
    fn clips_region_end_to_the_reference() {
        let selection = resolve_region(&references(), &"chr1:8-20".parse().unwrap()).unwrap();
        assert_eq!(selection.query.to_string(), "chr1:8-10");
        assert_eq!(selection.range.start, 7);
        assert_eq!(selection.range.end, 10);
    }

    #[test]
    fn rejects_unknown_and_outside_regions() {
        assert!(resolve_region(&references(), &"chr2".parse().unwrap()).is_err());
        assert!(resolve_region(&references(), &"chr1:11".parse().unwrap()).is_err());
    }
}
