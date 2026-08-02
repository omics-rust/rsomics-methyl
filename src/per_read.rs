use std::path::Path;

use noodles::core::Region;
use rsomics_bamio::raw::{RawRecord, RawRecordEncoder};
use rsomics_common::{Result, RsomicsError};

use crate::alignment::{AlignmentFilter, invalid_record};
use crate::calling::AlignmentCaller;
use crate::context::SequenceContext;
use crate::reference::IndexedReference;
use crate::selection::{AlignmentRecordResult, alignment_records, resolve_region};

#[derive(Clone, Debug)]
pub struct PerReadOptions {
    pub region: Option<Region>,
    pub minimum_mapping_quality: u8,
    pub minimum_base_quality: u8,
    pub ignore_flags: u16,
    pub require_flags: u16,
    pub ignore_nh: bool,
}

impl Default for PerReadOptions {
    fn default() -> Self {
        Self {
            region: None,
            minimum_mapping_quality: 10,
            minimum_base_quality: 5,
            ignore_flags: 0,
            require_flags: 0,
            ignore_nh: false,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct PerReadMetric {
    name: String,
    chromosome: String,
    start: u64,
    methylated: u64,
    unmethylated: u64,
}

impl PerReadMetric {
    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn chromosome(&self) -> &str {
        &self.chromosome
    }

    pub fn start(&self) -> u64 {
        self.start
    }

    pub fn methylated(&self) -> u64 {
        self.methylated
    }

    pub fn unmethylated(&self) -> u64 {
        self.unmethylated
    }

    pub fn informative_bases(&self) -> u64 {
        self.methylated + self.unmethylated
    }

    pub fn percentage(&self) -> f64 {
        let informative = self.informative_bases();
        if informative == 0 {
            0.0
        } else {
            self.methylated as f64 * 100.0 / informative as f64
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PerReadStats {
    pub input_records: u64,
    pub filtered_records: u64,
    pub output_records: u64,
}

pub fn per_read(
    input: &Path,
    reference: &Path,
    options: PerReadOptions,
    mut emit: impl FnMut(PerReadMetric) -> Result<()>,
) -> Result<PerReadStats> {
    if options.minimum_base_quality == 0 {
        return Err(RsomicsError::ConfigError(
            "minimum base quality must be positive".into(),
        ));
    }
    let mut reader = rsomics_bamio::open_indexed_alignment(input, Some(reference))?;
    let header = reader
        .read_header()
        .map_err(|error| alignment_error(input, error))?;
    let indexed_reference = IndexedReference::open(reference)?;
    let references = indexed_reference.validate_header(&header)?;
    let selection = options
        .region
        .as_ref()
        .map(|region| resolve_region(&references, region))
        .transpose()?;
    let filter = AlignmentFilter {
        minimum_mapping_quality: options.minimum_mapping_quality,
        ignore_flags: options.ignore_flags,
        require_flags: options.require_flags,
        reject_duplicates: false,
        reject_singletons: false,
        reject_discordant: false,
        reject_multimappers: !options.ignore_nh,
    };
    let mut caller = PerReadCaller {
        calls: AlignmentCaller::new(indexed_reference, references, options.minimum_base_quality),
    };
    let mut encoder = RawRecordEncoder::new();
    let mut stats = PerReadStats::default();
    let mut process = |result: AlignmentRecordResult| -> Result<()> {
        let record = result.map_err(|error| alignment_error(input, error))?;
        let record = encoder.encode(&header, record.as_ref())?;
        if let Some(selection) = &selection {
            let reference_id = usize::try_from(record.reference_sequence_id())
                .map_err(|error| invalid_record(&record, error))?;
            let start = u64::try_from(record.alignment_start())
                .map_err(|error| invalid_record(&record, error))?;
            if !selection.range.contains(reference_id, start) {
                return Ok(());
            }
        }
        stats.input_records = checked_increment(stats.input_records, "input record")?;
        if !filter.passes(&record)? {
            stats.filtered_records = checked_increment(stats.filtered_records, "filtered record")?;
            return Ok(());
        }
        emit(caller.metric(&record)?)?;
        stats.output_records = checked_increment(stats.output_records, "output record")?;
        Ok(())
    };
    let records = alignment_records(&mut reader, &header, selection.as_ref())
        .map_err(|error| alignment_error(input, error))?;
    for result in records {
        process(result)?;
    }
    Ok(stats)
}

struct PerReadCaller {
    calls: AlignmentCaller,
}

impl PerReadCaller {
    fn metric(&mut self, record: &RawRecord) -> Result<PerReadMetric> {
        let mut methylated = 0u64;
        let mut unmethylated = 0u64;
        let location = self.calls.visit(record, |call| {
            if call.context == SequenceContext::Cpg {
                if call.methylated {
                    methylated = checked_increment(methylated, "methylated count")?;
                } else {
                    unmethylated = checked_increment(unmethylated, "unmethylated count")?;
                }
            }
            Ok(())
        })?;
        let name = String::from_utf8(record.name().to_vec())
            .map_err(|_| invalid_record(record, "read name is not UTF-8"))?;
        Ok(PerReadMetric {
            name,
            chromosome: location.chromosome,
            start: location.start,
            methylated,
            unmethylated,
        })
    }
}

fn checked_increment(value: u64, field: &str) -> Result<u64> {
    value
        .checked_add(1)
        .ok_or_else(|| RsomicsError::InvalidInput(format!("{field} overflows")))
}

fn alignment_error(path: &Path, error: impl std::fmt::Display) -> RsomicsError {
    RsomicsError::InvalidInput(format!("reading alignment {}: {error}", path.display()))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn raw(cigar: &[(u8, u32)], sequence: &[u8]) -> RawRecord {
        raw_with_qualities(cigar, sequence, &vec![40; sequence.len()])
    }

    fn raw_with_qualities(cigar: &[(u8, u32)], sequence: &[u8], qualities: &[u8]) -> RawRecord {
        assert_eq!(sequence.len(), qualities.len());
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
        for &(kind, length) in cigar {
            payload.extend_from_slice(&((length << 4) | u32::from(kind)).to_le_bytes());
        }
        for pair in sequence.chunks(2) {
            let high = base_code(pair[0]);
            let low = pair.get(1).copied().map_or(0, base_code);
            payload.push(high << 4 | low);
        }
        payload.extend_from_slice(qualities);
        RawRecord::try_from(payload).unwrap()
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

    fn caller(sequence: &[u8]) -> (tempfile::TempDir, PerReadCaller) {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("reference.fa");
        let mut fasta = b">chr1\n".to_vec();
        fasta.extend_from_slice(sequence);
        fasta.push(b'\n');
        std::fs::write(&path, fasta).unwrap();
        std::fs::write(
            directory.path().join("reference.fa.fai"),
            format!(
                "chr1\t{}\t6\t{}\t{}\n",
                sequence.len(),
                sequence.len(),
                sequence.len() + 1
            ),
        )
        .unwrap();
        let caller = PerReadCaller {
            calls: AlignmentCaller::new(
                IndexedReference::open(&path).unwrap(),
                vec![crate::reference::ReferenceSequence {
                    name: "chr1".into(),
                    length: sequence.len() as u64,
                }],
                PerReadOptions::default().minimum_base_quality,
            ),
        };
        (directory, caller)
    }

    #[test]
    fn consumes_the_shared_long_cigar_contract() {
        let (_directory, mut caller) = caller(b"CG");
        let mut record = raw(&[(4, 1), (3, 1)], b"C");
        let count = usize::from(u16::MAX) + 1;
        let mut value = Vec::with_capacity(5 + count * 4);
        value.push(b'I');
        value.extend_from_slice(&(count as u32).to_le_bytes());
        for _ in 0..count - 1 {
            value.extend_from_slice(&((1u32 << 4) | 5).to_le_bytes());
        }
        value.extend_from_slice(&(1u32 << 4).to_le_bytes());
        record.append_aux(*b"CG", b'B', &value).unwrap();

        let metric = caller.metric(&record).unwrap();
        assert_eq!(metric.methylated(), 1);
        assert_eq!(metric.informative_bases(), 1);
    }

    #[test]
    fn reads_cpgs_beyond_the_upstream_ten_kilobase_window() {
        let mut reference = vec![b'A'; 20_003];
        reference[0..2].copy_from_slice(b"CG");
        reference[20_001..20_003].copy_from_slice(b"CG");
        let (_directory, mut caller) = caller(&reference);
        let record = raw(&[(0, 1), (3, 20_000), (0, 1)], b"CC");

        let metric = caller.metric(&record).unwrap();
        assert_eq!(metric.methylated(), 2);
        assert_eq!(metric.informative_bases(), 2);
    }

    #[test]
    fn low_quality_bases_do_not_advance_the_alignment_twice() {
        let (_directory, mut caller) = caller(b"CGCG");
        let record = raw_with_qualities(&[(0, 1), (2, 1), (0, 1)], b"CC", &[0, 40]);

        let metric = caller.metric(&record).unwrap();
        assert_eq!(metric.methylated(), 1);
        assert_eq!(metric.informative_bases(), 1);
    }
}
