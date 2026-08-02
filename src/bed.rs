use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::Path;

use flate2::read::MultiGzDecoder;
use rsomics_common::{Result, RsomicsError};
use rsomics_intervals::Interval;

use crate::reference::ReferenceSequence;

pub(crate) struct BedSelection {
    references: Vec<ReferenceIntervals>,
}

#[derive(Default)]
struct ReferenceIntervals {
    top: Vec<Interval<usize>>,
    bottom: Vec<Interval<usize>>,
}

enum StrandSelection {
    Both,
    Top,
    Bottom,
}

impl BedSelection {
    pub(crate) fn load(
        path: &Path,
        references: &[ReferenceSequence],
        keep_strand: bool,
    ) -> Result<Self> {
        let reference_ids = references
            .iter()
            .enumerate()
            .map(|(id, reference)| (reference.name.as_ref(), id))
            .collect::<HashMap<_, _>>();
        let mut selection = Self {
            references: (0..references.len())
                .map(|_| ReferenceIntervals::default())
                .collect(),
        };
        for (index, line) in reader(path)?.lines().enumerate() {
            let line_number = index + 1;
            let line = line.map_err(|error| bed_error(path, line_number, error))?;
            let trimmed = line.trim();
            if trimmed.is_empty()
                || trimmed.starts_with('#')
                || trimmed
                    .split_whitespace()
                    .next()
                    .is_some_and(|field| matches!(field, "track" | "browser"))
            {
                continue;
            }
            let fields = trimmed.split_whitespace().collect::<Vec<_>>();
            if fields.len() < 3 {
                return Err(bed_error(
                    path,
                    line_number,
                    "expected at least three fields",
                ));
            }
            let reference_id = *reference_ids.get(fields[0]).ok_or_else(|| {
                bed_error(
                    path,
                    line_number,
                    format!("unknown reference {}", fields[0]),
                )
            })?;
            let reference_length = references[reference_id].length;
            let start = coordinate(path, line_number, "start", fields[1])?;
            let end = coordinate(path, line_number, "end", fields[2])?;
            if start >= end {
                return Err(bed_error(
                    path,
                    line_number,
                    format!("start {start} must be less than end {end}"),
                ));
            }
            if start >= reference_length {
                return Err(bed_error(
                    path,
                    line_number,
                    format!(
                        "interval {start}-{end} lies outside {} length {reference_length}",
                        fields[0]
                    ),
                ));
            }
            let interval = Interval::new(reference_id, start, end.min(reference_length))
                .map_err(|error| bed_error(path, line_number, error))?;
            let strand = if keep_strand {
                match fields.get(5).copied() {
                    None | Some(".") => StrandSelection::Both,
                    Some("+") => StrandSelection::Top,
                    Some("-") => StrandSelection::Bottom,
                    Some(value) => {
                        return Err(bed_error(
                            path,
                            line_number,
                            format!("invalid BED strand {value}"),
                        ));
                    }
                }
            } else {
                StrandSelection::Both
            };
            let intervals = &mut selection.references[reference_id];
            match strand {
                StrandSelection::Both => {
                    intervals.top.push(interval);
                    intervals.bottom.push(interval);
                }
                StrandSelection::Top => intervals.top.push(interval),
                StrandSelection::Bottom => intervals.bottom.push(interval),
            }
        }
        for (reference_id, intervals) in selection.references.iter_mut().enumerate() {
            merge(reference_id, &mut intervals.top)?;
            merge(reference_id, &mut intervals.bottom)?;
        }
        Ok(selection)
    }

    pub(crate) fn contains(&self, reference_id: usize, position: u64, top: bool) -> bool {
        let Some(intervals) = self.intervals(reference_id, top) else {
            return false;
        };
        let index = intervals.partition_point(|interval| interval.end() <= position);
        intervals
            .get(index)
            .is_some_and(|interval| interval.start() <= position)
    }

    pub(crate) fn overlaps(&self, reference_id: usize, start: u64, end: u64, top: bool) -> bool {
        if start >= end {
            return false;
        }
        let Some(intervals) = self.intervals(reference_id, top) else {
            return false;
        };
        let index = intervals.partition_point(|interval| interval.end() <= start);
        intervals
            .get(index)
            .is_some_and(|interval| interval.start() < end)
    }

    fn intervals(&self, reference_id: usize, top: bool) -> Option<&[Interval<usize>]> {
        let intervals = self.references.get(reference_id)?;
        Some(if top {
            &intervals.top
        } else {
            &intervals.bottom
        })
    }
}

fn reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let mut file = File::open(path).map_err(|error| {
        RsomicsError::InvalidInput(format!("opening BED {}: {error}", path.display()))
    })?;
    let mut magic = [0; 2];
    let count = file.read(&mut magic).map_err(|error| {
        RsomicsError::InvalidInput(format!("reading BED {}: {error}", path.display()))
    })?;
    file.seek(SeekFrom::Start(0)).map_err(|error| {
        RsomicsError::InvalidInput(format!("reading BED {}: {error}", path.display()))
    })?;
    if count == magic.len() && magic == [0x1f, 0x8b] {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

fn coordinate(path: &Path, line: usize, field: &str, value: &str) -> Result<u64> {
    value.parse().map_err(|error| {
        bed_error(
            path,
            line,
            format!("invalid {field} coordinate {value}: {error}"),
        )
    })
}

fn merge(reference_id: usize, intervals: &mut Vec<Interval<usize>>) -> Result<()> {
    intervals.sort_unstable_by_key(|interval| (interval.start(), interval.end()));
    let mut merged: Vec<Interval<usize>> = Vec::with_capacity(intervals.len());
    for interval in intervals.drain(..) {
        if let Some(previous) = merged.last_mut()
            && interval.start() <= previous.end()
        {
            if interval.end() > previous.end() {
                *previous = Interval::new(reference_id, previous.start(), interval.end())
                    .map_err(|error| RsomicsError::InvalidInput(error.to_string()))?;
            }
        } else {
            merged.push(interval);
        }
    }
    *intervals = merged;
    Ok(())
}

fn bed_error(path: &Path, line: usize, error: impl std::fmt::Display) -> RsomicsError {
    RsomicsError::InvalidInput(format!(
        "reading BED {} line {line}: {error}",
        path.display()
    ))
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use flate2::Compression;
    use flate2::write::GzEncoder;

    use super::*;

    fn references() -> Vec<ReferenceSequence> {
        vec![ReferenceSequence {
            name: "chr1".into(),
            length: 20,
        }]
    }

    #[test]
    fn merges_intervals_and_applies_bed_strands() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("selection.bed");
        std::fs::write(
            &path,
            b"track name=x\nchr1 2 5 a 0 +\nchr1\t5\t8\tb\t0\t+\nchr1 6 12 c 0 -\nchr1 10 14 d 0 -\nchr1 15 16\n",
        )
        .unwrap();
        let selection = BedSelection::load(&path, &references(), true).unwrap();

        assert!(selection.contains(0, 2, true));
        assert!(selection.contains(0, 7, true));
        assert!(!selection.contains(0, 8, true));
        assert!(!selection.contains(0, 3, false));
        assert!(selection.contains(0, 7, false));
        assert!(selection.contains(0, 10, false));
        assert!(selection.contains(0, 15, true));
        assert!(selection.contains(0, 15, false));
        assert!(selection.overlaps(0, 7, 11, true));
        assert!(selection.overlaps(0, 7, 11, false));
        assert!(!selection.overlaps(0, 8, 10, true));
    }

    #[test]
    fn reads_gzip_by_magic_and_clips_the_reference_end() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("selection.data");
        let file = File::create(&path).unwrap();
        let mut encoder = GzEncoder::new(file, Compression::default());
        encoder.write_all(b"chr1 18 25\n").unwrap();
        encoder.finish().unwrap();

        let selection = BedSelection::load(&path, &references(), false).unwrap();
        assert!(selection.contains(0, 19, true));
        assert!(!selection.contains(0, 20, false));
    }

    #[test]
    fn rejects_malformed_and_outside_intervals() {
        let directory = tempfile::tempdir().unwrap();
        for (name, value) in [
            ("negative", "chr1 -1 3\n"),
            ("empty", "chr1 3 3\n"),
            ("outside", "chr1 20 21\n"),
            ("unknown", "chr2 1 2\n"),
            ("strand", "chr1 1 2 x 0 ?\n"),
        ] {
            let path = directory.path().join(name);
            std::fs::write(&path, value).unwrap();
            assert!(BedSelection::load(&path, &references(), true).is_err());
        }
    }

    #[test]
    fn empty_bed_matches_nothing() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("empty.bed");
        std::fs::write(&path, b"# no intervals\n").unwrap();
        let selection = BedSelection::load(&path, &references(), false).unwrap();
        assert!(!selection.contains(0, 0, true));
        assert!(!selection.overlaps(0, 0, 20, false));
    }
}
