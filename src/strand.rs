use rsomics_bamio::raw::RawRecord;
use rsomics_common::{Result, RsomicsError};

const PAIRED: u16 = 0x1;
const REVERSE: u16 = 0x10;
const READ_1: u16 = 0x40;
const READ_2: u16 = 0x80;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum BisulfiteStrand {
    Ot,
    Ob,
    Ctot,
    Ctob,
}

impl BisulfiteStrand {
    pub(crate) fn is_top(self) -> bool {
        matches!(self, Self::Ot | Self::Ctot)
    }
}

pub(crate) fn bisulfite_strand(record: &RawRecord) -> Result<BisulfiteStrand> {
    let flags = record.flags();
    let converted = record
        .aux_value(*b"XG")
        .and_then(|value| value.first())
        .copied()
        .filter(|value| matches!(value, b'C' | b'G'));
    if converted == Some(b'C') {
        return Ok(if flags & 0x51 == 0x41 {
            BisulfiteStrand::Ot
        } else if flags & 0x51 == 0x51 || flags & 0x91 == 0x81 {
            BisulfiteStrand::Ctot
        } else if flags & 0x91 == 0x91 {
            BisulfiteStrand::Ot
        } else if flags & REVERSE != 0 {
            BisulfiteStrand::Ctot
        } else {
            BisulfiteStrand::Ot
        });
    }
    if converted == Some(b'G') {
        return Ok(if flags & 0x51 == 0x41 {
            BisulfiteStrand::Ctob
        } else if flags & 0x51 == 0x51 || flags & 0x91 == 0x81 {
            BisulfiteStrand::Ob
        } else if flags & 0x91 == 0x91 {
            BisulfiteStrand::Ctob
        } else if flags & REVERSE != 0 {
            BisulfiteStrand::Ob
        } else {
            BisulfiteStrand::Ctob
        });
    }
    if flags & PAIRED == 0 {
        return Ok(if flags & REVERSE == 0 {
            BisulfiteStrand::Ot
        } else {
            BisulfiteStrand::Ob
        });
    }
    if flags & 0x50 == 0x50 {
        Ok(BisulfiteStrand::Ob)
    } else if flags & READ_1 != 0 || flags & 0x90 == 0x90 {
        Ok(BisulfiteStrand::Ot)
    } else if flags & READ_2 != 0 {
        Ok(BisulfiteStrand::Ob)
    } else {
        Err(RsomicsError::InvalidInput(format!(
            "cannot determine bisulfite strand for read {}",
            String::from_utf8_lossy(record.name())
        )))
    }
}

pub(crate) fn aux_integer(record: &RawRecord, tag: [u8; 2]) -> Result<Option<i64>> {
    let Some(kind) = record.aux_type(tag) else {
        return Ok(None);
    };
    let invalid = || {
        invalid_read(
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
        return Err(invalid_read(record, "NH must be nonnegative"));
    }
    Ok(Some(integer))
}

pub(crate) fn invalid_read(record: &RawRecord, error: impl std::fmt::Display) -> RsomicsError {
    RsomicsError::InvalidInput(format!(
        "read {}: {error}",
        String::from_utf8_lossy(record.name())
    ))
}
