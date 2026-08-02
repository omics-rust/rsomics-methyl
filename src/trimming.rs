use std::str::FromStr;

use rsomics_common::{Result, RsomicsError};

use crate::strand::BisulfiteStrand;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ReadBounds([u64; 4]);

impl ReadBounds {
    pub fn values(self) -> [u64; 4] {
        self.0
    }

    fn pair(self, read: u8) -> (u64, u64) {
        if read == 2 {
            (self.0[2], self.0[3])
        } else {
            (self.0[0], self.0[1])
        }
    }
}

impl FromStr for ReadBounds {
    type Err = String;

    fn from_str(value: &str) -> std::result::Result<Self, Self::Err> {
        let values = value
            .split(',')
            .map(|field| {
                field
                    .parse::<u64>()
                    .map_err(|error| format!("invalid read bounds {value}: {error}"))
            })
            .collect::<std::result::Result<Vec<_>, _>>()?;
        let values = <[u64; 4]>::try_from(values)
            .map_err(|_| format!("read bounds require four comma-separated integers: {value}"))?;
        Ok(Self(values))
    }
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct TrimmingOptions {
    inclusion: [Option<ReadBounds>; 4],
    fixed_ends: [Option<ReadBounds>; 4],
}

impl TrimmingOptions {
    pub fn set_inclusion(&mut self, strand: BisulfiteStrand, bounds: ReadBounds) -> Result<()> {
        for (start, end) in [bounds.pair(1), bounds.pair(2)] {
            if start > 0 && end > 0 && start > end {
                return Err(RsomicsError::ConfigError(format!(
                    "{} inclusion start {start} exceeds end {end}",
                    strand.label()
                )));
            }
        }
        self.inclusion[index(strand)] = Some(bounds);
        Ok(())
    }

    pub fn set_fixed_ends(&mut self, strand: BisulfiteStrand, bounds: ReadBounds) {
        self.fixed_ends[index(strand)] = Some(bounds);
    }

    pub(crate) fn includes(
        &self,
        strand: BisulfiteStrand,
        read: u8,
        sequence_length: u64,
        query_position: u64,
    ) -> Result<bool> {
        let position = query_position
            .checked_add(1)
            .ok_or_else(|| RsomicsError::InvalidInput("read position overflows".into()))?;
        if let Some(bounds) = self.inclusion[index(strand)] {
            let (start, end) = bounds.pair(read);
            if (start > 0 && position < start) || (end > 0 && position > end) {
                return Ok(false);
            }
        }
        if let Some(bounds) = self.fixed_ends[index(strand)] {
            let (left, right) = bounds.pair(read);
            let right_start = sequence_length.saturating_sub(right);
            if query_position < left || query_position >= right_start {
                return Ok(false);
            }
        }
        Ok(true)
    }
}

fn index(strand: BisulfiteStrand) -> usize {
    match strand {
        BisulfiteStrand::Ot => 0,
        BisulfiteStrand::Ob => 1,
        BisulfiteStrand::Ctot => 2,
        BisulfiteStrand::Ctob => 3,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_zero_bounds_without_hidden_process_state() {
        assert_eq!(
            "5,0,0,0".parse::<ReadBounds>().unwrap().values(),
            [5, 0, 0, 0]
        );
        assert!("1,2,3".parse::<ReadBounds>().is_err());
        assert!("1,2,3,x".parse::<ReadBounds>().is_err());
    }

    #[test]
    fn inclusion_positions_are_one_based_and_inclusive() {
        let mut options = TrimmingOptions::default();
        options
            .set_inclusion(BisulfiteStrand::Ot, "5,8,0,0".parse().unwrap())
            .unwrap();
        assert!(!options.includes(BisulfiteStrand::Ot, 1, 10, 3).unwrap());
        assert!(options.includes(BisulfiteStrand::Ot, 1, 10, 4).unwrap());
        assert!(options.includes(BisulfiteStrand::Ot, 1, 10, 7).unwrap());
        assert!(!options.includes(BisulfiteStrand::Ot, 1, 10, 8).unwrap());
    }

    #[test]
    fn fixed_bounds_remove_counts_from_each_end() {
        let mut options = TrimmingOptions::default();
        options.set_fixed_ends(BisulfiteStrand::Ot, "5,2,0,0".parse().unwrap());
        assert!(!options.includes(BisulfiteStrand::Ot, 1, 10, 4).unwrap());
        assert!(options.includes(BisulfiteStrand::Ot, 1, 10, 5).unwrap());
        assert!(options.includes(BisulfiteStrand::Ot, 1, 10, 7).unwrap());
        assert!(!options.includes(BisulfiteStrand::Ot, 1, 10, 8).unwrap());
    }

    #[test]
    fn bounds_select_the_matching_strand_and_read() {
        let mut options = TrimmingOptions::default();
        options
            .set_inclusion(BisulfiteStrand::Ob, "0,0,3,7".parse().unwrap())
            .unwrap();
        assert!(!options.includes(BisulfiteStrand::Ob, 2, 10, 1).unwrap());
        assert!(options.includes(BisulfiteStrand::Ob, 2, 10, 2).unwrap());
        assert!(options.includes(BisulfiteStrand::Ob, 2, 10, 6).unwrap());
        assert!(!options.includes(BisulfiteStrand::Ob, 2, 10, 7).unwrap());
        assert!(options.includes(BisulfiteStrand::Ot, 2, 10, 1).unwrap());
    }

    #[test]
    fn invalid_inclusion_order_is_rejected() {
        let mut options = TrimmingOptions::default();
        assert!(
            options
                .set_inclusion(BisulfiteStrand::Ctob, "8,2,0,0".parse().unwrap())
                .is_err()
        );
    }
}
