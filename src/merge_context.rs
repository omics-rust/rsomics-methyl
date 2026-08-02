use std::collections::{BTreeMap, HashSet};
use std::io::{BufRead, Write};
use std::ops::Range;
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

use crate::reference::IndexedReference;

const MAX_LINE_LENGTH: usize = 1024 * 1024;

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MergeContextStats {
    pub input_records: u64,
    pub output_records: u64,
    pub merged_records: u64,
}

#[derive(Clone, Debug)]
struct Metric {
    chromosome: String,
    start: u64,
    end: u64,
    methylated: u64,
    unmethylated: u64,
}

impl Metric {
    fn merge(&mut self, other: &Self, line_number: u64) -> Result<()> {
        self.methylated = self
            .methylated
            .checked_add(other.methylated)
            .ok_or_else(|| line_error(line_number, "methylated count overflows"))?;
        self.unmethylated = self
            .unmethylated
            .checked_add(other.unmethylated)
            .ok_or_else(|| line_error(line_number, "unmethylated count overflows"))?;
        Ok(())
    }

    fn write(&self, output: &mut dyn Write) -> Result<()> {
        let total = self
            .methylated
            .checked_add(self.unmethylated)
            .ok_or_else(|| RsomicsError::InvalidInput("methylation depth overflows".into()))?;
        if total == 0 {
            return Err(RsomicsError::InvalidInput(
                "methylation depth must be positive".into(),
            ));
        }
        let percentage = u128::from(self.methylated) * 100 / u128::from(total);
        writeln!(
            output,
            "{}\t{}\t{}\t{}\t{}\t{}",
            self.chromosome, self.start, self.end, percentage, self.methylated, self.unmethylated
        )
        .map_err(RsomicsError::Io)
    }
}

fn context_span(
    reference: &mut IndexedReference,
    chromosome: &str,
    position: u64,
) -> Result<Range<u64>> {
    let length = reference.length(chromosome)?;
    if position >= length {
        return Err(reference.error(format!(
            "{chromosome}:{position} is outside reference length {length}"
        )));
    }
    let position = usize::try_from(position).map_err(|error| reference.error(error))?;
    let length = usize::try_from(length).map_err(|error| reference.error(error))?;
    let start = position.saturating_sub(2);
    let end = position.saturating_add(3).min(length);
    let offset = position - start;
    let sequence = reference.sequence(chromosome, start..end)?;
    match sequence[offset].to_ascii_uppercase() {
        b'C' if sequence
            .get(offset + 1)
            .is_some_and(|base| base.eq_ignore_ascii_case(&b'G')) =>
        {
            Ok(position as u64..position as u64 + 2)
        }
        b'C' if sequence
            .get(offset + 2)
            .is_some_and(|base| base.eq_ignore_ascii_case(&b'G')) =>
        {
            Ok(position as u64..position as u64 + 3)
        }
        b'G' if offset >= 1 && sequence[offset - 1].eq_ignore_ascii_case(&b'C') => {
            Ok(position as u64 - 1..position as u64 + 1)
        }
        b'G' if offset >= 2 && sequence[offset - 2].eq_ignore_ascii_case(&b'C') => {
            Ok(position as u64 - 2..position as u64 + 1)
        }
        b'C' | b'G' => Ok(position as u64..position as u64 + 1),
        base => Err(reference.error(format!(
            "metric at {chromosome}:{position} targets reference base {} instead of C or G",
            char::from(base)
        ))),
    }
}

pub fn merge_context(
    mut input: impl BufRead,
    output: &mut dyn Write,
    reference: &Path,
) -> Result<MergeContextStats> {
    let mut reference = IndexedReference::open(reference)?;
    let mut stats = MergeContextStats::default();
    let mut pending = BTreeMap::<(u64, u64), Metric>::new();
    let mut seen_chromosomes = HashSet::new();
    let mut chromosome = None::<String>;
    let mut last_start = None::<u64>;
    let mut line = String::new();
    let mut line_number = 0u64;

    output
        .write_all(b"track type=\"bedGraph\" description=\"merged Methylation metrics\"\n")
        .map_err(RsomicsError::Io)?;

    loop {
        line.clear();
        let bytes = input.read_line(&mut line).map_err(RsomicsError::Io)?;
        if bytes == 0 {
            break;
        }
        line_number = line_number
            .checked_add(1)
            .ok_or_else(|| RsomicsError::InvalidInput("input line count overflows".into()))?;
        if line.len() > MAX_LINE_LENGTH {
            return Err(line_error(line_number, "line exceeds 1 MiB"));
        }
        let line = line.trim_end_matches(['\n', '\r']);
        if line.is_empty() || line.starts_with("track") {
            continue;
        }
        let mut metric = parse_metric(line, line_number)?;
        if chromosome.as_deref() != Some(metric.chromosome.as_str()) {
            flush_all(&mut pending, output, &mut stats)?;
            if !seen_chromosomes.insert(metric.chromosome.clone()) {
                return Err(line_error(
                    line_number,
                    format!(
                        "chromosome {} reappears after it was finalized",
                        metric.chromosome
                    ),
                ));
            }
            chromosome = Some(metric.chromosome.clone());
            last_start = None;
        }
        if last_start.is_some_and(|start| metric.start <= start) {
            return Err(line_error(
                line_number,
                "coordinates must be strictly increasing within each chromosome",
            ));
        }
        last_start = Some(metric.start);
        let span = context_span(&mut reference, &metric.chromosome, metric.start)?;
        metric.start = span.start;
        metric.end = span.end;
        stats.input_records = stats
            .input_records
            .checked_add(1)
            .ok_or_else(|| RsomicsError::InvalidInput("input record count overflows".into()))?;
        let key = (metric.start, metric.end);
        if let Some(existing) = pending.get_mut(&key) {
            existing.merge(&metric, line_number)?;
            stats.merged_records = stats.merged_records.checked_add(1).ok_or_else(|| {
                RsomicsError::InvalidInput("merged record count overflows".into())
            })?;
        } else {
            pending.insert(key, metric);
        }
        let settled = last_start
            .and_then(|start| start.checked_add(1))
            .ok_or_else(|| line_error(line_number, "coordinate overflows"))?;
        flush_settled(&mut pending, settled, output, &mut stats)?;
    }
    flush_all(&mut pending, output, &mut stats)?;
    Ok(stats)
}

fn parse_metric(line: &str, line_number: u64) -> Result<Metric> {
    let fields = line.split('\t').collect::<Vec<_>>();
    if fields.len() != 6 {
        return Err(line_error(
            line_number,
            format!("expected 6 tab-separated fields, found {}", fields.len()),
        ));
    }
    if fields[0].is_empty() {
        return Err(line_error(line_number, "chromosome name is empty"));
    }
    let start = parse_u64(fields[1], line_number, "start")?;
    let end = parse_u64(fields[2], line_number, "end")?;
    if end
        != start
            .checked_add(1)
            .ok_or_else(|| line_error(line_number, "start overflows"))?
    {
        return Err(line_error(
            line_number,
            "input metric must cover exactly one cytosine coordinate",
        ));
    }
    let percentage = parse_u64(fields[3], line_number, "percentage")?;
    if percentage > 100 {
        return Err(line_error(line_number, "percentage exceeds 100"));
    }
    let methylated = parse_u64(fields[4], line_number, "methylated count")?;
    let unmethylated = parse_u64(fields[5], line_number, "unmethylated count")?;
    let total = methylated
        .checked_add(unmethylated)
        .ok_or_else(|| line_error(line_number, "methylation depth overflows"))?;
    if total == 0 {
        return Err(line_error(
            line_number,
            "methylation depth must be positive",
        ));
    }
    Ok(Metric {
        chromosome: fields[0].to_owned(),
        start,
        end,
        methylated,
        unmethylated,
    })
}

fn parse_u64(value: &str, line_number: u64, field: &str) -> Result<u64> {
    value
        .parse()
        .map_err(|error| line_error(line_number, format!("invalid {field}: {error}")))
}

fn flush_settled(
    pending: &mut BTreeMap<(u64, u64), Metric>,
    settled: u64,
    output: &mut dyn Write,
    stats: &mut MergeContextStats,
) -> Result<()> {
    while pending
        .first_key_value()
        .is_some_and(|(_, metric)| metric.end <= settled)
    {
        let (_, metric) = pending.pop_first().unwrap();
        write_metric(metric, output, stats)?;
    }
    Ok(())
}

fn flush_all(
    pending: &mut BTreeMap<(u64, u64), Metric>,
    output: &mut dyn Write,
    stats: &mut MergeContextStats,
) -> Result<()> {
    while let Some((_, metric)) = pending.pop_first() {
        write_metric(metric, output, stats)?;
    }
    Ok(())
}

fn write_metric(
    metric: Metric,
    output: &mut dyn Write,
    stats: &mut MergeContextStats,
) -> Result<()> {
    metric.write(output)?;
    stats.output_records = stats
        .output_records
        .checked_add(1)
        .ok_or_else(|| RsomicsError::InvalidInput("output record count overflows".into()))?;
    Ok(())
}

fn line_error(line_number: u64, message: impl std::fmt::Display) -> RsomicsError {
    RsomicsError::InvalidInput(format!("line {line_number}: {message}"))
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::*;

    fn reference() -> tempfile::TempDir {
        let directory = tempfile::tempdir().unwrap();
        std::fs::write(
            directory.path().join("reference.fa"),
            b">chr1\nACGCCAGGCA\n>chr2\nCG\n",
        )
        .unwrap();
        std::fs::write(
            directory.path().join("reference.fa.fai"),
            b"chr1\t10\t6\t10\t11\nchr2\t2\t23\t2\t3\n",
        )
        .unwrap();
        directory
    }

    #[test]
    fn merges_cpg_and_chg_without_reordering_intervening_rows() {
        let directory = reference();
        let input = b"track type=\"bedGraph\"\n\
chr1\t1\t2\t50\t1\t1\n\
chr1\t2\t3\t100\t2\t0\n\
chr1\t3\t4\t0\t0\t1\n\
chr1\t4\t5\t50\t1\t1\n\
chr1\t6\t7\t0\t0\t2\n";
        let mut output = Vec::new();
        let stats = merge_context(
            Cursor::new(input),
            &mut output,
            &directory.path().join("reference.fa"),
        )
        .unwrap();
        assert_eq!(
            String::from_utf8(output).unwrap(),
            "track type=\"bedGraph\" description=\"merged Methylation metrics\"\n\
chr1\t1\t3\t75\t3\t1\n\
chr1\t3\t4\t0\t0\t1\n\
chr1\t4\t7\t25\t1\t3\n"
        );
        assert_eq!(
            stats,
            MergeContextStats {
                input_records: 5,
                output_records: 3,
                merged_records: 2,
            }
        );
    }

    #[test]
    fn rejects_unsorted_and_reappearing_chromosomes() {
        let directory = reference();
        let reference = directory.path().join("reference.fa");
        let mut output = Vec::new();
        let error = merge_context(
            Cursor::new(b"chr1\t2\t3\t50\t1\t1\nchr1\t1\t2\t50\t1\t1\n"),
            &mut output,
            &reference,
        )
        .unwrap_err();
        assert!(error.to_string().contains("strictly increasing"));

        let mut output = Vec::new();
        let error = merge_context(
            Cursor::new(b"chr1\t1\t2\t50\t1\t1\nchr2\t0\t1\t50\t1\t1\nchr1\t2\t3\t50\t1\t1\n"),
            &mut output,
            &reference,
        )
        .unwrap_err();
        assert!(error.to_string().contains("reappears"));
    }

    #[test]
    fn rejects_unknown_reference_and_invalid_rows() {
        let directory = reference();
        let reference = directory.path().join("reference.fa");
        let mut output = Vec::new();
        let error = merge_context(
            Cursor::new(b"missing\t0\t1\t50\t1\t1\n"),
            &mut output,
            &reference,
        )
        .unwrap_err();
        assert!(error.to_string().contains("unknown chromosome"));

        let mut output = Vec::new();
        let error = merge_context(
            Cursor::new(b"chr1\t0\t2\t50\t1\t1\n"),
            &mut output,
            &reference,
        )
        .unwrap_err();
        assert!(error.to_string().contains("exactly one"));
    }
}
