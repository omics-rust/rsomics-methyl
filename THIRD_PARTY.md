# Third-party provenance

MethylDackel is Copyright (c) 2015-2023 Devon Ryan and contributors and is
available under the MIT License. Compatibility behavior is reviewed against
release 0.6.1 (`b6db120e96ec8cf9ab44e1b1074d2aa7af876932`) and corrected
upstream revision `3c77bda12141e99d80234d416e668a90ec70b3f7`.
The positional count contract, Agresti-Coull interval, inclusion-bound
suggestions, and exhaustive cytosine-context orientation are adapted from its
`MBias.c`, `svg.c`, and `extract.c` modules. BED strand-to-bisulfite-strand
mapping follows its `bed.c`, `extract.c`, and `MBias.c` behavior.
Non-CpG conversion-efficiency filtering follows its `common.c` contract.

The historical Rust implementation and fixtures are team-owned. Adapted
algorithms retain this provenance; the rsomics implementation is licensed
under MIT OR Apache-2.0.
