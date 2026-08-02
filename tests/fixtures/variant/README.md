# Opposite-strand variant fixture

The CpG and CHG references each have four usable opposite-strand observations
at each cytosine. Each forward cytosine has two non-reference bases plus one
`N`; each reverse has one non-reference base. `input.bam`, its index, and the
FASTA index were generated from the tracked source files with samtools 1.24.

With minimum opposite depth 4, current MethylDackel and rsomics agree at
maximum fractions 0.35 and 0.6. The other cases freeze three deliberate
corrections to the documented contract: equality is allowed, `N` contributes
to neither usable depth nor the non-reference numerator, and merged output is
excluded symmetrically when either complementary cytosine is likely variant.
The exhaustive report also omits an excluded cytosine instead of recreating it
as a zero-coverage gap.
