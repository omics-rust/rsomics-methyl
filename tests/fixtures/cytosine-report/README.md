# Cytosine-report fixture

The project-owned reference covers forward and reverse CpG, CHG, and CHH
contexts, contig boundaries, and an ambiguous base. `empty.sam` makes every
reported cytosine a zero-coverage call. BAM, BAI, and FAI files were generated
with samtools 1.24. The expected report was produced by MethylDackel revision
`3c77bda12141e99d80234d416e668a90ec70b3f7` with all contexts enabled.
