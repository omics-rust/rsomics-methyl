#![forbid(unsafe_code)]
#![allow(clippy::missing_errors_doc)]

mod context;
pub mod extract;
pub mod merge_context;
mod reference;
mod strand;

pub use context::{ReferenceStrand, SequenceContext};
