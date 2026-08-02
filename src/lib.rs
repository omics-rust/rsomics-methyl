#![forbid(unsafe_code)]
#![allow(clippy::missing_errors_doc)]

mod alignment;
mod context;
pub mod extract;
pub mod merge_context;
pub mod per_read;
mod reference;
mod strand;

pub use context::{ReferenceStrand, SequenceContext};
