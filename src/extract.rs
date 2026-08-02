use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use noodles::core::Region;
use rsomics_bamio::raw::RawRecord;
use rsomics_common::{Result, RsomicsError};
use rsomics_pileup::{Column, ColumnEntry, PileupEngine, PileupError, PileupOptions};
use smallvec::SmallVec;

use crate::alignment::{AlignmentFilter, DUPLICATE};
use crate::bed::BedSelection;
use crate::calling::read_number;
use crate::context::{CytosineContext, ReferenceStrand, SequenceContext, classify};
use crate::conversion::{ConversionFilter, validate_conversion_efficiency};
use crate::reference::{IndexedReference, ReferenceSequence};
use crate::selection::{ReferenceRange, alignment_error, resolve_region, visit_alignment_records};
use crate::strand::{BisulfiteStrand, bisulfite_strand};
use crate::trimming::TrimmingOptions;

const PAIRED: u16 = 0x1;
const MATE_UNMAPPED: u16 = 0x8;

#[derive(Clone, Debug)]
pub struct ExtractOptions {
    pub region: Option<Region>,
    pub bed: Option<PathBuf>,
    pub keep_bed_strand: bool,
    pub trimming: TrimmingOptions,
    pub minimum_mapping_quality: u8,
    pub minimum_base_quality: u8,
    pub minimum_conversion_efficiency: f64,
    pub minimum_opposite_depth: u64,
    pub maximum_variant_fraction: f64,
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
            bed: None,
            keep_bed_strand: false,
            trimming: TrimmingOptions::default(),
            minimum_mapping_quality: 10,
            minimum_base_quality: 5,
            minimum_conversion_efficiency: 0.0,
            minimum_opposite_depth: 0,
            maximum_variant_fraction: 0.0,
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
    chromosome: Arc<str>,
    start: u64,
    end: u64,
    context: SequenceContext,
    strand: ReferenceStrand,
    trinucleotide: [u8; 3],
    methylated: u64,
    unmethylated: u64,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ExcludedVariantSite {
    chromosome: Arc<str>,
    start: u64,
    context: SequenceContext,
    strand: ReferenceStrand,
    opposite_depth: u64,
    variant_bases: u64,
}

impl ExcludedVariantSite {
    pub fn chromosome(&self) -> &str {
        &self.chromosome
    }

    pub fn start(&self) -> u64 {
        self.start
    }

    pub fn context(&self) -> SequenceContext {
        self.context
    }

    pub fn strand(&self) -> ReferenceStrand {
        self.strand
    }

    pub fn opposite_depth(&self) -> u64 {
        self.opposite_depth
    }

    pub fn variant_bases(&self) -> u64 {
        self.variant_bases
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ExtractEvent {
    Site(SiteMetric),
    ExcludedVariant(ExcludedVariantSite),
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
    pub excluded_variant_sites: u64,
}

struct Extractor {
    reference: IndexedReference,
    references: Vec<ReferenceSequence>,
    selection: Option<ReferenceRange>,
    bed: Option<BedSelection>,
    variant_filter: Option<VariantFilter>,
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
        emit: &mut impl FnMut(ExtractEvent) -> Result<()>,
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
        if let Some(event) = self.column(column, reference_id, position)? {
            emit(event)?;
        }
        if self.exhaustive {
            self.report_reference_id = reference_id;
            self.report_position = position
                .checked_add(1)
                .ok_or_else(|| RsomicsError::InvalidInput("site end overflows".into()))?;
        }
        Ok(())
    }

    fn emit_remaining(&mut self, emit: &mut impl FnMut(ExtractEvent) -> Result<()>) -> Result<()> {
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
        emit: &mut impl FnMut(ExtractEvent) -> Result<()>,
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
        emit: &mut impl FnMut(ExtractEvent) -> Result<()>,
    ) -> Result<()> {
        while self.report_position < end {
            if let Some(metric) = self.zero_metric()? {
                emit(ExtractEvent::Site(metric))?;
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
        let length = usize::try_from(reference.length).map_err(invalid_coordinate)?;
        let Some(context) = classify(&mut self.reference, &reference.name, length, position)?
        else {
            return Ok(None);
        };
        if !self.includes(context.kind)
            || self.bed.as_ref().is_some_and(|selection| {
                !selection.contains(
                    self.report_reference_id,
                    self.report_position,
                    context.strand == ReferenceStrand::Forward,
                )
            })
        {
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
    ) -> Result<Option<ExtractEvent>> {
        self.stats.examined_columns = checked_increment(self.stats.examined_columns, "column")?;
        let reference = self.references.get(reference_id).ok_or_else(|| {
            RsomicsError::InvalidInput(format!("pileup reference ID {reference_id} is absent"))
        })?;
        let position = usize::try_from(raw_position).map_err(invalid_coordinate)?;
        let length = usize::try_from(reference.length).map_err(invalid_coordinate)?;
        let Some(context) = classify(&mut self.reference, &reference.name, length, position)?
        else {
            return Ok(None);
        };
        if !self.includes(context.kind)
            || self.bed.as_ref().is_some_and(|selection| {
                !selection.contains(
                    reference_id,
                    raw_position,
                    context.strand == ReferenceStrand::Forward,
                )
            })
        {
            return Ok(None);
        }
        let mut counts = EvidenceCounts::default();
        if column
            .entries()
            .any(|entry| entry.record().flags() & PAIRED != 0)
        {
            let mut evidence = SmallVec::<[Evidence<'_>; 16]>::with_capacity(column.len());
            for entry in column.entries() {
                if let Some(value) = self.evidence(entry)? {
                    evidence.push(value);
                }
            }
            adjust_overlaps(&mut evidence);
            for value in evidence {
                self.observe(value, reference_id, raw_position, context, &mut counts)?;
            }
        } else {
            for entry in column.entries() {
                if let Some(value) = self.evidence(entry)? {
                    self.observe(value, reference_id, raw_position, context, &mut counts)?;
                }
            }
        }
        let depth = counts
            .methylated
            .checked_add(counts.unmethylated)
            .ok_or_else(|| RsomicsError::InvalidInput("methylation depth overflows".into()))?;
        if self
            .variant_filter
            .is_some_and(|filter| filter.excludes(counts.opposite))
        {
            self.stats.excluded_variant_sites =
                checked_increment(self.stats.excluded_variant_sites, "excluded variant site")?;
            return Ok(Some(ExtractEvent::ExcludedVariant(ExcludedVariantSite {
                chromosome: reference.name.clone(),
                start: raw_position,
                context: context.kind,
                strand: context.strand,
                opposite_depth: counts.opposite.depth,
                variant_bases: counts.opposite.variant_bases,
            })));
        }
        if !self.exhaustive && depth < self.options.minimum_depth {
            return Ok(None);
        }
        self.stats.emitted_sites = checked_increment(self.stats.emitted_sites, "site")?;
        let start = u64::try_from(position).map_err(invalid_coordinate)?;
        Ok(Some(ExtractEvent::Site(SiteMetric {
            chromosome: reference.name.clone(),
            start,
            end: start
                .checked_add(1)
                .ok_or_else(|| RsomicsError::InvalidInput("site end overflows".into()))?,
            context: context.kind,
            strand: context.strand,
            trinucleotide: context.trinucleotide,
            methylated: counts.methylated,
            unmethylated: counts.unmethylated,
        })))
    }

    fn evidence<'a>(&self, entry: ColumnEntry<'a>) -> Result<Option<Evidence<'a>>> {
        let projection = entry.projection();
        if projection.is_deletion || projection.is_reference_skip {
            return Ok(None);
        }
        let record = entry.record();
        let strand = bisulfite_strand(record)?;
        if !self.options.trimming.includes(
            strand,
            read_number(record),
            u64::try_from(record.sequence_len()).map_err(invalid_coordinate)?,
            u64::try_from(projection.qpos).map_err(invalid_coordinate)?,
        )? {
            return Ok(None);
        }
        let quality = record
            .quality_scores()
            .get(projection.qpos)
            .copied()
            .unwrap_or(u8::MAX);
        Ok(Some(Evidence {
            record,
            base: record.seq_nibble(projection.qpos),
            quality,
            strand,
            paired_with_mate: false,
        }))
    }

    fn observe(
        &self,
        value: Evidence<'_>,
        reference_id: usize,
        position: u64,
        context: CytosineContext,
        counts: &mut EvidenceCounts,
    ) -> Result<()> {
        if value.quality < self.options.minimum_base_quality
            || self.bed.as_ref().is_some_and(|selection| {
                !selection.contains(reference_id, position, value.strand.is_top())
            })
        {
            return Ok(());
        }
        if value.strand.is_top() != (context.strand == ReferenceStrand::Forward) {
            if self.variant_filter.is_some() {
                counts.opposite.observe(value.strand, value.base)?;
            }
            return Ok(());
        }
        match (value.strand.is_top(), value.base) {
            (true, 2) | (false, 4) => {
                counts.methylated = checked_increment(counts.methylated, "methylated count")?;
            }
            (true, 8) | (false, 1) => {
                counts.unmethylated = checked_increment(counts.unmethylated, "unmethylated count")?;
            }
            _ => {}
        }
        Ok(())
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
    paired_with_mate: bool,
}

#[derive(Default)]
struct EvidenceCounts {
    methylated: u64,
    unmethylated: u64,
    opposite: OppositeEvidence,
}

#[derive(Clone, Copy)]
struct VariantFilter {
    minimum_depth: u64,
    maximum_fraction: f64,
}

impl VariantFilter {
    fn new(minimum_depth: u64, maximum_fraction: f64) -> Result<Option<Self>> {
        if !maximum_fraction.is_finite() || !(0.0..=1.0).contains(&maximum_fraction) {
            return Err(RsomicsError::ConfigError(
                "maximum variant fraction must be between 0 and 1".into(),
            ));
        }
        Ok((minimum_depth > 0).then_some(Self {
            minimum_depth,
            maximum_fraction,
        }))
    }

    fn excludes(self, evidence: OppositeEvidence) -> bool {
        evidence.depth >= self.minimum_depth
            && evidence.variant_bases as f64 / evidence.depth as f64 > self.maximum_fraction
    }
}

#[derive(Clone, Copy, Default)]
struct OppositeEvidence {
    depth: u64,
    variant_bases: u64,
}

impl OppositeEvidence {
    fn observe(&mut self, strand: BisulfiteStrand, base: u8) -> Result<()> {
        if base == 15 {
            return Ok(());
        }
        self.depth = checked_increment(self.depth, "opposite-strand depth")?;
        let expected = if strand.is_top() { 4 } else { 2 };
        if base != expected {
            self.variant_bases = checked_increment(self.variant_bases, "opposite-strand variant")?;
        }
        Ok(())
    }
}

pub fn extract(
    input: &Path,
    reference: &Path,
    options: ExtractOptions,
    emit: impl FnMut(SiteMetric) -> Result<()>,
) -> Result<ExtractStats> {
    let mut emit = emit;
    extract_with_mode(input, reference, options, false, |event| match event {
        ExtractEvent::Site(metric) => emit(metric),
        ExtractEvent::ExcludedVariant(_) => Ok(()),
    })
}

pub fn extract_events(
    input: &Path,
    reference: &Path,
    options: ExtractOptions,
    emit: impl FnMut(ExtractEvent) -> Result<()>,
) -> Result<ExtractStats> {
    extract_with_mode(input, reference, options, false, emit)
}

pub fn extract_all_cytosines(
    input: &Path,
    reference: &Path,
    options: ExtractOptions,
    emit: impl FnMut(SiteMetric) -> Result<()>,
) -> Result<ExtractStats> {
    let mut emit = emit;
    extract_with_mode(input, reference, options, true, |event| match event {
        ExtractEvent::Site(metric) => emit(metric),
        ExtractEvent::ExcludedVariant(_) => Ok(()),
    })
}

fn extract_with_mode(
    input: &Path,
    reference: &Path,
    options: ExtractOptions,
    exhaustive: bool,
    mut emit: impl FnMut(ExtractEvent) -> Result<()>,
) -> Result<ExtractStats> {
    if options.minimum_base_quality == 0 {
        return Err(RsomicsError::ConfigError(
            "minimum base quality must be positive".into(),
        ));
    }
    validate_conversion_efficiency(options.minimum_conversion_efficiency)?;
    let variant_filter = VariantFilter::new(
        options.minimum_opposite_depth,
        options.maximum_variant_fraction,
    )?;
    if !options.cpg && !options.chg && !options.chh {
        return Err(RsomicsError::ConfigError(
            "at least one methylation context must be enabled".into(),
        ));
    }
    if options.keep_bed_strand && options.bed.is_none() {
        return Err(RsomicsError::ConfigError(
            "BED strand filtering requires --bed".into(),
        ));
    }
    let mut reader = rsomics_bamio::open_indexed_alignment(input, Some(reference))?;
    let header = reader
        .read_header()
        .map_err(|error| alignment_error(input, error))?;
    let indexed_reference = IndexedReference::open(reference)?;
    let references = indexed_reference.validate_header(&header)?;
    let mut conversion = ConversionFilter::new(
        reference,
        references.clone(),
        options.minimum_base_quality,
        options.minimum_conversion_efficiency,
    )?;
    let bed = options
        .bed
        .as_deref()
        .map(|path| BedSelection::load(path, &references, options.keep_bed_strand))
        .transpose()?;
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
        bed,
        variant_filter,
        options,
        exhaustive,
        report_reference_id,
        report_position,
        stats: ExtractStats::default(),
    };
    visit_alignment_records(input, &mut reader, &header, selection.as_ref(), |raw| {
        extractor.stats.input_records =
            checked_increment(extractor.stats.input_records, "input record")?;
        let mut passes = filter.passes(&raw)?;
        if passes && let Some(conversion) = conversion.as_mut() {
            passes = conversion.passes(&raw)?;
        }
        if !passes {
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
    })?;
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
    if evidence.len() <= 32 {
        for index in 0..evidence.len() {
            if !eligible_for_overlap(&evidence[index]) {
                continue;
            }
            let name = evidence[index].record.name();
            let Some(first) = (0..index).rev().find(|&candidate| {
                !evidence[candidate].paired_with_mate
                    && eligible_for_overlap(&evidence[candidate])
                    && evidence[candidate].record.name() == name
            }) else {
                continue;
            };
            let (left, right) = evidence.split_at_mut(index);
            left[first].paired_with_mate = true;
            right[0].paired_with_mate = true;
            adjust_pair(&mut left[first], &mut right[0]);
        }
        return;
    }

    let mut pending = HashMap::<&[u8], usize>::new();
    for index in 0..evidence.len() {
        if !eligible_for_overlap(&evidence[index]) {
            continue;
        }
        let name = evidence[index].record.name();
        if let Some(first) = pending.remove(name) {
            let (left, right) = evidence.split_at_mut(index);
            left[first].paired_with_mate = true;
            right[0].paired_with_mate = true;
            adjust_pair(&mut left[first], &mut right[0]);
        } else {
            pending.insert(name, index);
        }
    }
}

fn eligible_for_overlap(evidence: &Evidence<'_>) -> bool {
    let flags = evidence.record.flags();
    flags & PAIRED != 0 && flags & MATE_UNMAPPED == 0
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
    fn rejects_invalid_variant_fractions() {
        for value in [f64::NAN, f64::INFINITY, -0.1, 1.1] {
            assert!(VariantFilter::new(1, value).is_err());
        }
    }

    #[test]
    fn overlap_adjustment_matches_across_depth_paths() {
        for pairs in [1, 17] {
            let records = (0..pairs)
                .flat_map(|pair| {
                    let name = format!("p{pair}");
                    [paired_record(&name), paired_record(&name)]
                })
                .collect::<Vec<_>>();
            let mut evidence = records
                .iter()
                .enumerate()
                .map(|(index, record)| Evidence {
                    record,
                    base: 2,
                    quality: if index % 2 == 0 { 40 } else { 30 },
                    strand: BisulfiteStrand::Ot,
                    paired_with_mate: false,
                })
                .collect::<Vec<_>>();

            adjust_overlaps(&mut evidence);

            for pair in evidence.chunks_exact(2) {
                assert_eq!([pair[0].quality, pair[1].quality], [48, 0]);
                assert!(pair[0].paired_with_mate);
                assert!(pair[1].paired_with_mate);
            }
        }
    }

    fn paired_record(name: &str) -> RawRecord {
        let mut payload = Vec::new();
        payload.extend_from_slice(&0i32.to_le_bytes());
        payload.extend_from_slice(&0i32.to_le_bytes());
        payload.push(u8::try_from(name.len() + 1).unwrap());
        payload.push(60);
        payload.extend_from_slice(&0u16.to_le_bytes());
        payload.extend_from_slice(&1u16.to_le_bytes());
        payload.extend_from_slice(&PAIRED.to_le_bytes());
        payload.extend_from_slice(&1u32.to_le_bytes());
        payload.extend_from_slice(&0i32.to_le_bytes());
        payload.extend_from_slice(&0i32.to_le_bytes());
        payload.extend_from_slice(&0i32.to_le_bytes());
        payload.extend_from_slice(name.as_bytes());
        payload.push(0);
        payload.extend_from_slice(&(1u32 << 4).to_le_bytes());
        payload.push(0x20);
        payload.push(40);
        RawRecord::try_from(payload).unwrap()
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
            bed: None,
            variant_filter: None,
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
            .emit_remaining(&mut |event| {
                let ExtractEvent::Site(metric) = event else {
                    panic!("zero-coverage scan emitted a variant event");
                };
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
