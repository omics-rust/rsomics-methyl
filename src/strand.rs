use rsomics_bamio::raw::RawRecord;
use rsomics_common::Result;

use crate::alignment::invalid_record;

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
        Err(invalid_record(record, "cannot determine bisulfite strand"))
    }
}
