use rsomics_bamio::raw::RawRecord;
use rsomics_common::{Result, RsomicsError};

const PAIRED: u16 = 0x1;
const PROPER_PAIR: u16 = 0x2;
const UNMAPPED: u16 = 0x4;
const MATE_UNMAPPED: u16 = 0x8;
pub(crate) const DUPLICATE: u16 = 0x400;

#[derive(Clone, Debug)]
pub(crate) struct AlignmentFilter {
    pub(crate) minimum_mapping_quality: u8,
    pub(crate) ignore_flags: u16,
    pub(crate) require_flags: u16,
    pub(crate) reject_duplicates: bool,
    pub(crate) reject_singletons: bool,
    pub(crate) reject_discordant: bool,
    pub(crate) reject_multimappers: bool,
}

impl AlignmentFilter {
    pub(crate) fn passes(&self, record: &RawRecord) -> Result<bool> {
        let flags = record.flags();
        if flags & UNMAPPED != 0
            || record.mapping_quality() < self.minimum_mapping_quality
            || flags & self.ignore_flags != 0
            || (self.require_flags != 0 && flags & self.require_flags != self.require_flags)
            || (self.reject_duplicates && flags & DUPLICATE != 0)
            || (self.reject_singletons
                && flags & (PAIRED | MATE_UNMAPPED) == (PAIRED | MATE_UNMAPPED))
            || (self.reject_discordant && flags & (PAIRED | PROPER_PAIR) == PAIRED)
        {
            return Ok(false);
        }
        if self.reject_multimappers
            && let Some(value) = aux_integer(record, *b"NH")?
            && value > 1
        {
            return Ok(false);
        }
        Ok(true)
    }
}

fn aux_integer(record: &RawRecord, tag: [u8; 2]) -> Result<Option<i64>> {
    let Some(kind) = record.aux_type(tag) else {
        return Ok(None);
    };
    let invalid = || {
        invalid_record(
            record,
            format!(
                "{}{} has invalid integer encoding",
                char::from(tag[0]),
                char::from(tag[1])
            ),
        )
    };
    let value = record.aux_value(tag).ok_or_else(invalid)?;
    let integer = match kind {
        b'c' => value
            .first()
            .map(|value| i64::from(i8::from_le_bytes([*value]))),
        b'C' => value.first().map(|value| i64::from(*value)),
        b's' => value
            .get(..2)
            .and_then(|value| <[u8; 2]>::try_from(value).ok())
            .map(i16::from_le_bytes)
            .map(i64::from),
        b'S' => value
            .get(..2)
            .and_then(|value| <[u8; 2]>::try_from(value).ok())
            .map(u16::from_le_bytes)
            .map(i64::from),
        b'i' => value
            .get(..4)
            .and_then(|value| <[u8; 4]>::try_from(value).ok())
            .map(i32::from_le_bytes)
            .map(i64::from),
        b'I' => value
            .get(..4)
            .and_then(|value| <[u8; 4]>::try_from(value).ok())
            .map(u32::from_le_bytes)
            .map(i64::from),
        _ => return Err(invalid()),
    };
    let integer = integer.ok_or_else(invalid)?;
    if integer < 0 {
        return Err(invalid_record(
            record,
            format!("{}{} must be nonnegative", tag[0] as char, tag[1] as char),
        ));
    }
    Ok(Some(integer))
}

pub(crate) fn invalid_record(record: &RawRecord, error: impl std::fmt::Display) -> RsomicsError {
    RsomicsError::InvalidInput(format!(
        "read {}: {error}",
        String::from_utf8_lossy(record.name())
    ))
}
