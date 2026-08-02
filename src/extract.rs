use std::collections::HashMap;
use std::path::Path;

use rsomics_bamio::raw::{RawRecord, RawRecordEncoder};
use rsomics_common::{Result, RsomicsError};
use rsomics_pileup::{Column, PileupEngine, PileupError, PileupOptions};

use crate::context::{ReferenceStrand, SequenceContext, classify};
use crate::reference::{IndexedReference, ReferenceSequence};
use crate::strand::{BisulfiteStrand, aux_integer, bisulfite_strand};

const PAIRED: u16 = 0x1;
const PROPER_PAIR: u16 = 0x2;
const UNMAPPED: u16 = 0x4;
const MATE_UNMAPPED: u16 = 0x8;
const DUPLICATE: u16 = 0x400;

#[derive(Clone, Debug)]
pub struct ExtractOptions {
    pub minimum_mapping_quality: u8,
    pub minimum_base_quality: u8,
    pub ignore_flags: u16,
    pub require_flags: u16,
    pub keep_duplicates: bool,
    pub keep_singletons: bool,
    pub keep_discordant: bool,
    pub ignore_nh: bool,
    pub minimum_depth: u64,
    pub cpg: bool,
    pub chg: bool,
    pub chh: bool,
}

impl Default for ExtractOptions {
    fn default() -> Self {
        Self {
            minimum_mapping_quality: 10,
            minimum_base_quality: 5,
            ignore_flags: 0x0f00,
            require_flags: 0,
            keep_duplicates: false,
            keep_singletons: false,
            keep_discordant: false,
            ignore_nh: false,
            minimum_depth: 1,
            cpg: true,
            chg: false,
            chh: false,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SiteMetric {
    chromosome: String,
    start: u64,
    end: u64,
    context: SequenceContext,
    strand: ReferenceStrand,
    methylated: u64,
    unmethylated: u64,
}

impl SiteMetric {
    pub fn chromosome(&self) -> &str {
        &self.chromosome
    }

    pub fn start(&self) -> u64 {
        self.start
    }

    pub fn end(&self) -> u64 {
        self.end
    }

    pub fn context(&self) -> SequenceContext {
        self.context
    }

    pub fn strand(&self) -> ReferenceStrand {
        self.strand
    }

    pub fn methylated(&self) -> u64 {
        self.methylated
    }

    pub fn unmethylated(&self) -> u64 {
        self.unmethylated
    }

    pub fn depth(&self) -> u64 {
        self.methylated + self.unmethylated
    }

    pub fn percentage(&self) -> u64 {
        let depth = self.depth();
        if depth == 0 {
            0
        } else {
            (u128::from(self.methylated) * 100 / u128::from(depth)) as u64
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ExtractStats {
    pub input_records: u64,
    pub filtered_records: u64,
    pub examined_columns: u64,
    pub emitted_sites: u64,
}

struct Extractor {
    reference: IndexedReference,
    references: Vec<ReferenceSequence>,
    options: ExtractOptions,
    stats: ExtractStats,
}

impl Extractor {
    fn column(&mut self, column: &Column<'_>) -> Result<Option<SiteMetric>> {
        self.stats.examined_columns = checked_increment(self.stats.examined_columns, "column")?;
        let reference_id = usize::try_from(column.reference_id()).map_err(invalid_coordinate)?;
        let reference = self.references.get(reference_id).ok_or_else(|| {
            RsomicsError::InvalidInput(format!("pileup reference ID {reference_id} is absent"))
        })?;
        let position = usize::try_from(column.position()).map_err(invalid_coordinate)?;
        let Some(context) = classify(&mut self.reference, &reference.name, position)? else {
            return Ok(None);
        };
        if !self.includes(context.kind) {
            return Ok(None);
        }
        let mut evidence = Vec::with_capacity(column.len());
        for entry in column.entries() {
            let projection = entry.projection();
            if projection.is_deletion || projection.is_reference_skip {
                continue;
            }
            let record = entry.record();
            let quality = record
                .quality_scores()
                .get(projection.qpos)
                .copied()
                .unwrap_or(u8::MAX);
            evidence.push(Evidence {
                record,
                base: record.seq_nibble(projection.qpos),
                quality,
                strand: bisulfite_strand(record)?,
            });
        }
        adjust_overlaps(&mut evidence);
        let mut methylated = 0u64;
        let mut unmethylated = 0u64;
        for value in evidence {
            if value.quality < self.options.minimum_base_quality
                || value.strand.is_top() != (context.strand == ReferenceStrand::Forward)
            {
                continue;
            }
            match (value.strand.is_top(), value.base) {
                (true, 2) | (false, 4) => {
                    methylated = checked_increment(methylated, "methylated count")?;
                }
                (true, 8) | (false, 1) => {
                    unmethylated = checked_increment(unmethylated, "unmethylated count")?;
                }
                _ => {}
            }
        }
        let depth = methylated
            .checked_add(unmethylated)
            .ok_or_else(|| RsomicsError::InvalidInput("methylation depth overflows".into()))?;
        if depth < self.options.minimum_depth {
            return Ok(None);
        }
        self.stats.emitted_sites = checked_increment(self.stats.emitted_sites, "site")?;
        let start = u64::try_from(position).map_err(invalid_coordinate)?;
        Ok(Some(SiteMetric {
            chromosome: reference.name.clone(),
            start,
            end: start
                .checked_add(1)
                .ok_or_else(|| RsomicsError::InvalidInput("site end overflows".into()))?,
            context: context.kind,
            strand: context.strand,
            methylated,
            unmethylated,
        }))
    }

    fn includes(&self, context: SequenceContext) -> bool {
        match context {
            SequenceContext::Cpg => self.options.cpg,
            SequenceContext::Chg => self.options.chg,
            SequenceContext::Chh => self.options.chh,
        }
    }
}

struct Evidence<'a> {
    record: &'a RawRecord,
    base: u8,
    quality: u8,
    strand: BisulfiteStrand,
}

pub fn extract(
    input: &Path,
    reference: &Path,
    options: ExtractOptions,
    mut emit: impl FnMut(SiteMetric) -> Result<()>,
) -> Result<ExtractStats> {
    if options.minimum_base_quality == 0 {
        return Err(RsomicsError::ConfigError(
            "minimum base quality must be positive".into(),
        ));
    }
    if !options.cpg && !options.chg && !options.chh {
        return Err(RsomicsError::ConfigError(
            "at least one methylation context must be enabled".into(),
        ));
    }
    let mut reader = rsomics_bamio::open_indexed_alignment(input, Some(reference))?;
    let header = reader
        .read_header()
        .map_err(|error| alignment_error(input, error))?;
    let indexed_reference = IndexedReference::open(reference)?;
    let references = indexed_reference.validate_header(&header)?;
    let lengths = references.iter().map(|reference| reference.length);
    let mut pileup = PileupEngine::new(lengths, PileupOptions::default());
    let mut extractor = Extractor {
        reference: indexed_reference,
        references,
        options,
        stats: ExtractStats::default(),
    };
    let mut encoder = RawRecordEncoder::new();
    for result in reader.records(&header) {
        let record = result.map_err(|error| alignment_error(input, error))?;
        let raw = encoder.encode(&header, record.as_ref())?;
        extractor.stats.input_records =
            checked_increment(extractor.stats.input_records, "input record")?;
        if !passes_filters(&raw, &extractor.options)? {
            extractor.stats.filtered_records =
                checked_increment(extractor.stats.filtered_records, "filtered record")?;
            continue;
        }
        pileup
            .push(raw)
            .map_err(|error| pileup_error(input, error))?;
        pileup.drain(|column| {
            if let Some(metric) = extractor.column(column)? {
                emit(metric)?;
            }
            Ok::<(), RsomicsError>(())
        })?;
    }
    pileup
        .finish()
        .map_err(|error| pileup_error(input, error))?;
    pileup.drain(|column| {
        if let Some(metric) = extractor.column(column)? {
            emit(metric)?;
        }
        Ok::<(), RsomicsError>(())
    })?;
    Ok(extractor.stats)
}

fn passes_filters(record: &RawRecord, options: &ExtractOptions) -> Result<bool> {
    let flags = record.flags();
    let ignored = if options.keep_duplicates {
        options.ignore_flags & !DUPLICATE
    } else {
        options.ignore_flags | DUPLICATE
    };
    if flags & UNMAPPED != 0
        || record.mapping_quality() < options.minimum_mapping_quality
        || flags & ignored != 0
        || (options.require_flags != 0 && flags & options.require_flags != options.require_flags)
        || (!options.keep_singletons
            && flags & (PAIRED | MATE_UNMAPPED) == (PAIRED | MATE_UNMAPPED))
        || (!options.keep_discordant && flags & (PAIRED | PROPER_PAIR) == PAIRED)
    {
        return Ok(false);
    }
    if !options.ignore_nh
        && let Some(value) = aux_integer(record, *b"NH")?
        && value > 1
    {
        return Ok(false);
    }
    Ok(true)
}

fn adjust_overlaps(evidence: &mut [Evidence<'_>]) {
    let mut pending = HashMap::<&[u8], usize>::new();
    for index in 0..evidence.len() {
        let flags = evidence[index].record.flags();
        if flags & PAIRED == 0 || flags & MATE_UNMAPPED != 0 {
            continue;
        }
        let name = evidence[index].record.name();
        if let Some(first) = pending.remove(name) {
            let (left, right) = evidence.split_at_mut(index);
            adjust_pair(&mut left[first], &mut right[0]);
        } else {
            pending.insert(name, index);
        }
    }
}

fn adjust_pair(first: &mut Evidence<'_>, second: &mut Evidence<'_>) {
    if first.strand.is_top() != second.strand.is_top() {
        return;
    }
    if first.base == second.base {
        if first.quality > second.quality {
            first.quality = boosted(first.quality);
            second.quality = 0;
        } else {
            second.quality = boosted(second.quality);
            first.quality = 0;
        }
    } else if first.quality > second.quality && first.base != 15 {
        first.quality -= second.quality;
        second.quality = 0;
    } else if second.quality > first.quality && second.base != 15 {
        second.quality -= first.quality;
        first.quality = 0;
    } else {
        first.quality = 0;
        second.quality = 0;
    }
}

fn boosted(quality: u8) -> u8 {
    quality.saturating_add(quality / 5)
}

fn checked_increment(value: u64, field: &str) -> Result<u64> {
    value
        .checked_add(1)
        .ok_or_else(|| RsomicsError::InvalidInput(format!("{field} count overflows")))
}

fn invalid_coordinate(error: impl std::fmt::Display) -> RsomicsError {
    RsomicsError::InvalidInput(format!("invalid pileup coordinate: {error}"))
}

fn alignment_error(path: &Path, error: impl std::fmt::Display) -> RsomicsError {
    RsomicsError::InvalidInput(format!("reading alignment {}: {error}", path.display()))
}

fn pileup_error(path: &Path, error: PileupError) -> RsomicsError {
    RsomicsError::InvalidInput(format!("building pileup from {}: {error}", path.display()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quality_boost_is_checked() {
        assert_eq!(boosted(40), 48);
        assert_eq!(boosted(250), 255);
    }
}
