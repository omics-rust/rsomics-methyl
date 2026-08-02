use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use noodles::core::Region;
use rsomics_common::{Result, RsomicsError};

use crate::alignment::{AlignmentFilter, DUPLICATE};
use crate::bed::BedSelection;
use crate::calling::AlignmentCaller;
use crate::context::SequenceContext;
use crate::conversion::{ConversionFilter, validate_conversion_efficiency};
use crate::reference::IndexedReference;
use crate::selection::{alignment_error, resolve_region, visit_alignment_records};
use crate::strand::BisulfiteStrand;
use crate::trimming::TrimmingOptions;

#[derive(Clone, Debug)]
pub struct MbiasOptions {
    pub region: Option<Region>,
    pub bed: Option<PathBuf>,
    pub keep_bed_strand: bool,
    pub trimming: TrimmingOptions,
    pub minimum_mapping_quality: u8,
    pub minimum_base_quality: u8,
    pub minimum_conversion_efficiency: f64,
    pub ignore_flags: u16,
    pub require_flags: u16,
    pub keep_duplicates: bool,
    pub keep_singletons: bool,
    pub keep_discordant: bool,
    pub ignore_nh: bool,
    pub cpg: bool,
    pub chg: bool,
    pub chh: bool,
}

impl Default for MbiasOptions {
    fn default() -> Self {
        Self {
            region: None,
            bed: None,
            keep_bed_strand: false,
            trimming: TrimmingOptions::default(),
            minimum_mapping_quality: 10,
            minimum_base_quality: 5,
            minimum_conversion_efficiency: 0.0,
            ignore_flags: 0x0f00,
            require_flags: 0,
            keep_duplicates: false,
            keep_singletons: false,
            keep_discordant: false,
            ignore_nh: false,
            cpg: true,
            chg: false,
            chh: false,
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MbiasMetric {
    strand: BisulfiteStrand,
    read: u8,
    position: u64,
    methylated: u64,
    unmethylated: u64,
}

impl MbiasMetric {
    pub fn strand(&self) -> BisulfiteStrand {
        self.strand
    }

    pub fn read(&self) -> u8 {
        self.read
    }

    pub fn position(&self) -> u64 {
        self.position
    }

    pub fn methylated(&self) -> u64 {
        self.methylated
    }

    pub fn unmethylated(&self) -> u64 {
        self.unmethylated
    }

    pub fn fraction(&self) -> f64 {
        self.methylated as f64 / (self.methylated + self.unmethylated) as f64
    }

    pub fn confidence_interval(&self) -> (f64, f64) {
        let counts = Counts {
            methylated: self.methylated,
            unmethylated: self.unmethylated,
        };
        (
            confidence_interval(counts, false),
            confidence_interval(counts, true),
        )
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MbiasSuggestion {
    strand: BisulfiteStrand,
    bounds: [u64; 4],
}

impl MbiasSuggestion {
    pub fn strand(&self) -> BisulfiteStrand {
        self.strand
    }

    pub fn bounds(&self) -> [u64; 4] {
        self.bounds
    }
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct MbiasStats {
    pub input_records: u64,
    pub filtered_records: u64,
    pub calls: u64,
}

pub struct MbiasResult {
    metrics: Vec<MbiasMetric>,
    suggestions: Vec<MbiasSuggestion>,
    stats: MbiasStats,
}

impl MbiasResult {
    pub fn metrics(&self) -> &[MbiasMetric] {
        &self.metrics
    }

    pub fn suggestions(&self) -> &[MbiasSuggestion] {
        &self.suggestions
    }

    pub fn stats(&self) -> &MbiasStats {
        &self.stats
    }
}

#[derive(Clone, Copy, Default)]
struct Counts {
    methylated: u64,
    unmethylated: u64,
}

pub fn mbias(input: &Path, reference: &Path, options: MbiasOptions) -> Result<MbiasResult> {
    if options.minimum_base_quality == 0 {
        return Err(RsomicsError::ConfigError(
            "minimum base quality must be positive".into(),
        ));
    }
    validate_conversion_efficiency(options.minimum_conversion_efficiency)?;
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
    let mut caller =
        AlignmentCaller::new(indexed_reference, references, options.minimum_base_quality);
    let mut counts = BTreeMap::<(BisulfiteStrand, u64, u8), Counts>::new();
    let mut stats = MbiasStats::default();
    visit_alignment_records(input, &mut reader, &header, selection.as_ref(), |record| {
        stats.input_records = increment(stats.input_records, "input record")?;
        if !filter.passes(&record)? {
            stats.filtered_records = increment(stats.filtered_records, "filtered record")?;
            return Ok(());
        }
        if let Some(conversion) = conversion.as_mut()
            && !conversion.passes(&record)?
        {
            stats.filtered_records = increment(stats.filtered_records, "filtered record")?;
            return Ok(());
        }
        let sequence_length = u64::try_from(record.sequence_len())
            .map_err(|error| RsomicsError::InvalidInput(error.to_string()))?;
        caller.visit(&record, |call| {
            if selection.as_ref().is_some_and(|selection| {
                !selection
                    .range
                    .contains(call.reference_id, call.reference_position)
            }) || bed.as_ref().is_some_and(|selection| {
                !selection.contains(
                    call.reference_id,
                    call.reference_position,
                    call.strand.is_top(),
                )
            }) || !includes(&options, call.context)
                || !options.trimming.includes(
                    call.strand,
                    call.read,
                    sequence_length,
                    call.query_position,
                )?
            {
                return Ok(());
            }
            let entry = counts
                .entry((call.strand, call.query_position, call.read))
                .or_default();
            if call.methylated {
                entry.methylated = increment(entry.methylated, "M-bias methylated count")?;
            } else {
                entry.unmethylated = increment(entry.unmethylated, "M-bias unmethylated count")?;
            }
            stats.calls = increment(stats.calls, "M-bias call")?;
            Ok(())
        })?;
        Ok(())
    })?;
    let metrics = counts
        .into_iter()
        .map(|((strand, query_position, read), counts)| {
            Ok(MbiasMetric {
                strand,
                read,
                position: query_position.checked_add(1).ok_or_else(|| {
                    RsomicsError::InvalidInput("M-bias position overflows".into())
                })?,
                methylated: counts.methylated,
                unmethylated: counts.unmethylated,
            })
        })
        .collect::<Result<Vec<_>>>()?;
    let suggestions = suggestions(&metrics)?;
    Ok(MbiasResult {
        metrics,
        suggestions,
        stats,
    })
}

fn includes(options: &MbiasOptions, context: SequenceContext) -> bool {
    match context {
        SequenceContext::Cpg => options.cpg,
        SequenceContext::Chg => options.chg,
        SequenceContext::Chh => options.chh,
    }
}

fn suggestions(metrics: &[MbiasMetric]) -> Result<Vec<MbiasSuggestion>> {
    let mut result = Vec::new();
    for strand in [
        BisulfiteStrand::Ot,
        BisulfiteStrand::Ob,
        BisulfiteStrand::Ctot,
        BisulfiteStrand::Ctob,
    ] {
        let length = metrics
            .iter()
            .filter(|metric| metric.strand == strand)
            .map(|metric| metric.position)
            .max()
            .unwrap_or(0);
        if length == 0 {
            continue;
        }
        let read_1 = threshold(metrics, strand, 1, length)?;
        let read_2 = threshold(metrics, strand, 2, length)?;
        result.push(MbiasSuggestion {
            strand,
            bounds: [read_1.0, read_1.1, read_2.0, read_2.1],
        });
    }
    Ok(result)
}

fn threshold(
    metrics: &[MbiasMetric],
    strand: BisulfiteStrand,
    read: u8,
    length: u64,
) -> Result<(u64, u64)> {
    let length =
        usize::try_from(length).map_err(|error| RsomicsError::InvalidInput(error.to_string()))?;
    let mut values = vec![None; length];
    for metric in metrics
        .iter()
        .filter(|metric| metric.strand == strand && metric.read == read)
    {
        let index = usize::try_from(metric.position - 1)
            .map_err(|error| RsomicsError::InvalidInput(error.to_string()))?;
        values[index] = Some(Counts {
            methylated: metric.methylated,
            unmethylated: metric.unmethylated,
        });
    }
    let middle = length / 2;
    let calibration = &values[length * 2 / 10..=length * 8 / 10];
    let mut total = 0u64;
    let mut average = 0.0;
    let mut minimum_ci = 1.0f64;
    let mut maximum_ci = 0.0f64;
    for counts in calibration.iter().flatten() {
        total = increment(total, "M-bias calibration position")?;
        average += fraction(*counts);
        minimum_ci = minimum_ci.min(confidence_interval(*counts, false));
        maximum_ci = maximum_ci.max(confidence_interval(*counts, true));
    }
    if total == 0 {
        return Ok((0, 0));
    }
    average /= total as f64;
    let left = (0..=middle)
        .rev()
        .find(|&index| is_outlier(values[index], average, minimum_ci, maximum_ci))
        .map_or(0, |index| index as u64 + 2);
    let right = (middle + 1..length)
        .find(|&index| is_outlier(values[index], average, minimum_ci, maximum_ci))
        .map_or(0, |index| index as u64);
    Ok((left, right))
}

fn is_outlier(counts: Option<Counts>, average: f64, minimum_ci: f64, maximum_ci: f64) -> bool {
    let Some(counts) = counts else {
        return false;
    };
    let value = fraction(counts);
    let separated = (confidence_interval(counts, true) < average && value < minimum_ci)
        || (confidence_interval(counts, false) > average && value > maximum_ci);
    separated && (value - average).abs() > 0.05
}

fn fraction(counts: Counts) -> f64 {
    counts.methylated as f64 / (counts.methylated + counts.unmethylated) as f64
}

fn confidence_interval(counts: Counts, upper: bool) -> f64 {
    const Z_SQUARED: f64 = 10.827_566_170_7;
    const Z: f64 = 3.290_526_731_5;
    let observations = (counts.methylated + counts.unmethylated) as f64;
    let adjusted_observations = observations + Z_SQUARED;
    let adjusted_fraction = (counts.methylated as f64 + 0.5 * Z_SQUARED) / adjusted_observations;
    let radius =
        Z * ((adjusted_fraction / adjusted_observations) * (1.0 - adjusted_fraction)).sqrt();
    if upper {
        (adjusted_fraction + radius).min(1.0)
    } else {
        (adjusted_fraction - radius).max(0.0)
    }
}

fn increment(value: u64, field: &str) -> Result<u64> {
    value
        .checked_add(1)
        .ok_or_else(|| RsomicsError::InvalidInput(format!("{field} overflows")))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture(name: &str) -> std::path::PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/extract")
            .join(name)
    }

    fn table(result: &MbiasResult) -> String {
        let mut output = "Strand\tRead\tPosition\tnMethylated\tnUnmethylated\n".to_owned();
        for metric in &result.metrics {
            output.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\n",
                metric.strand.label(),
                metric.read,
                metric.position,
                metric.methylated,
                metric.unmethylated
            ));
        }
        output
    }

    #[test]
    fn metrics_match_the_live_methyldackel_golden() {
        let result = mbias(
            &fixture("synthetic.bam"),
            &fixture("synthetic.fa"),
            MbiasOptions::default(),
        )
        .unwrap();
        assert_eq!(
            table(&result),
            std::fs::read_to_string(fixture("expected.mbias.tsv")).unwrap()
        );
        assert_eq!(result.stats.input_records, 8);
        assert_eq!(result.stats.filtered_records, 3);
        assert_eq!(result.stats.calls, 50);
        assert_eq!(
            result
                .suggestions
                .iter()
                .map(|suggestion| (suggestion.strand.label(), suggestion.bounds))
                .collect::<Vec<_>>(),
            [("OT", [0; 4]), ("OB", [0; 4])]
        );
    }

    #[test]
    fn region_limits_mapped_columns_without_rebasing_read_positions() {
        let result = mbias(
            &fixture("synthetic.bam"),
            &fixture("synthetic.fa"),
            MbiasOptions {
                region: Some("chrSynthetic:5-10".parse().unwrap()),
                ..MbiasOptions::default()
            },
        )
        .unwrap();
        let expected = std::fs::read_to_string(fixture("expected.mbias.region.tsv")).unwrap();
        assert_eq!(table(&result), expected);
    }

    #[test]
    fn suggestions_find_biased_read_ends() {
        let metrics = (1..=10)
            .map(|position| MbiasMetric {
                strand: BisulfiteStrand::Ot,
                read: 1,
                position,
                methylated: if position == 1 {
                    0
                } else if position == 10 {
                    2_000
                } else {
                    1_000
                },
                unmethylated: if position == 1 {
                    2_000
                } else if position == 10 {
                    0
                } else {
                    1_000
                },
            })
            .collect::<Vec<_>>();
        assert_eq!(
            threshold(&metrics, BisulfiteStrand::Ot, 1, 10).unwrap(),
            (2, 9)
        );
    }
}
