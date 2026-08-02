use std::collections::HashMap;
use std::path::Path;

use noodles::core::Region;
use rsomics_bamio::raw::{RawRecord, RawRecordEncoder};
use rsomics_common::{Result, RsomicsError};
use rsomics_pileup::{Column, PileupEngine, PileupError, PileupOptions};

use crate::alignment::{AlignmentFilter, DUPLICATE};
use crate::calling::read_number;
use crate::context::{ReferenceStrand, SequenceContext, classify};
use crate::reference::{IndexedReference, ReferenceSequence};
use crate::selection::{AlignmentRecordResult, ReferenceRange, alignment_records, resolve_region};
use crate::strand::{BisulfiteStrand, bisulfite_strand};
use crate::trimming::TrimmingOptions;

const PAIRED: u16 = 0x1;
const MATE_UNMAPPED: u16 = 0x8;

#[derive(Clone, Debug)]
pub struct ExtractOptions {
    pub region: Option<Region>,
    pub trimming: TrimmingOptions,
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
            region: None,
            trimming: TrimmingOptions::default(),
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
    trinucleotide: [u8; 3],
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

    pub fn trinucleotide(&self) -> [u8; 3] {
        self.trinucleotide
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
    selection: Option<ReferenceRange>,
    options: ExtractOptions,
    exhaustive: bool,
    report_reference_id: usize,
    report_position: u64,
    stats: ExtractStats,
}

impl Extractor {
    fn visit_column(
        &mut self,
        column: &Column<'_>,
        emit: &mut impl FnMut(SiteMetric) -> Result<()>,
    ) -> Result<()> {
        let reference_id = usize::try_from(column.reference_id()).map_err(invalid_coordinate)?;
        let position = u64::try_from(column.position()).map_err(invalid_coordinate)?;
        if self
            .selection
            .is_some_and(|selection| !selection.contains(reference_id, position))
        {
            return Ok(());
        }
        if reference_id >= self.references.len() {
            return Err(RsomicsError::InvalidInput(format!(
                "pileup reference ID {reference_id} is absent"
            )));
        }
        if self.exhaustive {
            self.emit_until(reference_id, position, emit)?;
        }
        if let Some(metric) = self.column(column, reference_id, position)? {
            emit(metric)?;
        }
        if self.exhaustive {
            self.report_reference_id = reference_id;
            self.report_position = position
                .checked_add(1)
                .ok_or_else(|| RsomicsError::InvalidInput("site end overflows".into()))?;
        }
        Ok(())
    }

    fn emit_remaining(&mut self, emit: &mut impl FnMut(SiteMetric) -> Result<()>) -> Result<()> {
        if !self.exhaustive {
            return Ok(());
        }
        if let Some(selection) = self.selection {
            return self.emit_until(selection.reference_id, selection.end, emit);
        }
        while self.report_reference_id < self.references.len() {
            let end = self.references[self.report_reference_id].length;
            self.emit_until(self.report_reference_id, end, emit)?;
            self.report_reference_id += 1;
            self.report_position = 0;
        }
        Ok(())
    }

    fn emit_until(
        &mut self,
        reference_id: usize,
        position: u64,
        emit: &mut impl FnMut(SiteMetric) -> Result<()>,
    ) -> Result<()> {
        while self.report_reference_id < reference_id {
            let end = self
                .references
                .get(self.report_reference_id)
                .ok_or_else(|| RsomicsError::InvalidInput("report reference is absent".into()))?
                .length;
            self.emit_current_range(end, emit)?;
            self.report_reference_id += 1;
            self.report_position = 0;
        }
        if self.report_reference_id != reference_id || self.report_position > position {
            return Err(RsomicsError::InvalidInput(
                "methylation report positions are out of order".into(),
            ));
        }
        self.emit_current_range(position, emit)
    }

    fn emit_current_range(
        &mut self,
        end: u64,
        emit: &mut impl FnMut(SiteMetric) -> Result<()>,
    ) -> Result<()> {
        while self.report_position < end {
            if let Some(metric) = self.zero_metric()? {
                emit(metric)?;
            }
            self.report_position = self
                .report_position
                .checked_add(1)
                .ok_or_else(|| RsomicsError::InvalidInput("report position overflows".into()))?;
        }
        Ok(())
    }

    fn zero_metric(&mut self) -> Result<Option<SiteMetric>> {
        let reference = self
            .references
            .get(self.report_reference_id)
            .ok_or_else(|| RsomicsError::InvalidInput("report reference is absent".into()))?;
        let position = usize::try_from(self.report_position).map_err(invalid_coordinate)?;
        let Some(context) = classify(&mut self.reference, &reference.name, position)? else {
            return Ok(None);
        };
        if !self.includes(context.kind) {
            return Ok(None);
        }
        self.stats.emitted_sites = checked_increment(self.stats.emitted_sites, "site")?;
        Ok(Some(SiteMetric {
            chromosome: reference.name.clone(),
            start: self.report_position,
            end: self
                .report_position
                .checked_add(1)
                .ok_or_else(|| RsomicsError::InvalidInput("site end overflows".into()))?,
            context: context.kind,
            strand: context.strand,
            trinucleotide: context.trinucleotide,
            methylated: 0,
            unmethylated: 0,
        }))
    }

    fn column(
        &mut self,
        column: &Column<'_>,
        reference_id: usize,
        raw_position: u64,
    ) -> Result<Option<SiteMetric>> {
        self.stats.examined_columns = checked_increment(self.stats.examined_columns, "column")?;
        let reference = self.references.get(reference_id).ok_or_else(|| {
            RsomicsError::InvalidInput(format!("pileup reference ID {reference_id} is absent"))
        })?;
        let position = usize::try_from(raw_position).map_err(invalid_coordinate)?;
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
            let strand = bisulfite_strand(record)?;
            if !self.options.trimming.includes(
                strand,
                read_number(record),
                u64::try_from(record.sequence_len()).map_err(invalid_coordinate)?,
                u64::try_from(projection.qpos).map_err(invalid_coordinate)?,
            )? {
                continue;
            }
            let quality = record
                .quality_scores()
                .get(projection.qpos)
                .copied()
                .unwrap_or(u8::MAX);
            evidence.push(Evidence {
                record,
                base: record.seq_nibble(projection.qpos),
                quality,
                strand,
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
        if !self.exhaustive && depth < self.options.minimum_depth {
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
            trinucleotide: context.trinucleotide,
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
    emit: impl FnMut(SiteMetric) -> Result<()>,
) -> Result<ExtractStats> {
    extract_with_mode(input, reference, options, false, emit)
}

pub fn extract_all_cytosines(
    input: &Path,
    reference: &Path,
    options: ExtractOptions,
    emit: impl FnMut(SiteMetric) -> Result<()>,
) -> Result<ExtractStats> {
    extract_with_mode(input, reference, options, true, emit)
}

fn extract_with_mode(
    input: &Path,
    reference: &Path,
    options: ExtractOptions,
    exhaustive: bool,
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
    let selection = options
        .region
        .as_ref()
        .map(|region| resolve_region(&references, region))
        .transpose()?;
    let lengths = references.iter().map(|reference| reference.length);
    let mut pileup = PileupEngine::new(lengths, PileupOptions::default());
    let filter = AlignmentFilter {
        minimum_mapping_quality: options.minimum_mapping_quality,
        ignore_flags: if options.keep_duplicates {
            options.ignore_flags & !DUPLICATE
        } else {
            options.ignore_flags
        },
        require_flags: options.require_flags,
        reject_duplicates: !options.keep_duplicates,
        reject_singletons: !options.keep_singletons,
        reject_discordant: !options.keep_discordant,
        reject_multimappers: !options.ignore_nh,
    };
    let (report_reference_id, report_position) = selection
        .as_ref()
        .map(|selection| (selection.range.reference_id, selection.range.start))
        .unwrap_or((0, 0));
    let mut extractor = Extractor {
        reference: indexed_reference,
        references,
        selection: selection.as_ref().map(|selection| selection.range),
        options,
        exhaustive,
        report_reference_id,
        report_position,
        stats: ExtractStats::default(),
    };
    let mut encoder = RawRecordEncoder::new();
    let mut process = |result: AlignmentRecordResult| -> Result<()> {
        let record = result.map_err(|error| alignment_error(input, error))?;
        let raw = encoder.encode(&header, record.as_ref())?;
        extractor.stats.input_records =
            checked_increment(extractor.stats.input_records, "input record")?;
        if !filter.passes(&raw)? {
            extractor.stats.filtered_records =
                checked_increment(extractor.stats.filtered_records, "filtered record")?;
            return Ok(());
        }
        pileup
            .push(raw)
            .map_err(|error| pileup_error(input, error))?;
        pileup.drain(|column| {
            extractor.visit_column(column, &mut emit)?;
            Ok::<(), RsomicsError>(())
        })?;
        Ok(())
    };
    let records = alignment_records(&mut reader, &header, selection.as_ref())
        .map_err(|error| alignment_error(input, error))?;
    for result in records {
        process(result)?;
    }
    pileup
        .finish()
        .map_err(|error| pileup_error(input, error))?;
    pileup.drain(|column| {
        extractor.visit_column(column, &mut emit)?;
        Ok::<(), RsomicsError>(())
    })?;
    extractor.emit_remaining(&mut emit)?;
    Ok(extractor.stats)
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

    #[test]
    fn exhaustive_scan_crosses_reference_boundaries() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("reference.fa");
        std::fs::write(&path, b">chr1\nCAG\n>chr2\nGTT\n").unwrap();
        std::fs::write(
            directory.path().join("reference.fa.fai"),
            b"chr1\t3\t6\t3\t4\nchr2\t3\t16\t3\t4\n",
        )
        .unwrap();
        let mut extractor = Extractor {
            reference: IndexedReference::open(&path).unwrap(),
            references: vec![
                ReferenceSequence {
                    name: "chr1".into(),
                    length: 3,
                },
                ReferenceSequence {
                    name: "chr2".into(),
                    length: 3,
                },
            ],
            selection: None,
            options: ExtractOptions {
                cpg: false,
                chg: true,
                chh: true,
                ..ExtractOptions::default()
            },
            exhaustive: true,
            report_reference_id: 0,
            report_position: 0,
            stats: ExtractStats::default(),
        };
        let mut observed = Vec::new();
        extractor
            .emit_remaining(&mut |metric| {
                observed.push((
                    metric.chromosome().to_owned(),
                    metric.start(),
                    metric.context(),
                    metric.trinucleotide(),
                ));
                Ok(())
            })
            .unwrap();
        assert_eq!(
            observed,
            [
                ("chr1".into(), 0, SequenceContext::Chg, *b"CAG"),
                ("chr1".into(), 2, SequenceContext::Chg, *b"CTG"),
                ("chr2".into(), 0, SequenceContext::Chh, *b"CNN"),
            ]
        );
    }
}
