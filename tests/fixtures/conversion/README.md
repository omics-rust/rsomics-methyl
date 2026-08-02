# Conversion fixture

The reference contains forward and reverse CpG, CHG, and CHH contexts. Reads
have non-CpG conversion efficiencies of zero, one half, one, or no informative
non-CpG bases. `input.bam`, its index, and the FASTA index were generated from
the tracked source files with samtools 1.24. The 0.5 and 0.75 extraction and
M-bias goldens were produced by MethylDackel revision
`3c77bda12141e99d80234d416e668a90ec70b3f7`.
