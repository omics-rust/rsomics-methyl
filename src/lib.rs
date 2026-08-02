#![forbid(unsafe_code)]
#![allow(clippy::missing_errors_doc)]

mod alignment;
mod bed;
mod calling;
mod context;
pub mod extract;
pub mod mbias;
pub mod merge_context;
pub mod per_read;
mod reference;
mod selection;
mod strand;
mod trimming;

pub use context::{ReferenceStrand, SequenceContext};
pub use strand::BisulfiteStrand;
pub use trimming::{ReadBounds, TrimmingOptions};
