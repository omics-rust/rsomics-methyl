# rsomics-methyl

`rsomics-methyl` is the rsomics product for bisulfite-sequencing methylation
extraction and bias QC. It consolidates the team-owned historical
`rsomics-methyldackel` implementation into one product family with
MethylDackel compatibility gates.

The source is under active reconstruction. The command tree currently exposes
checked CpG/CHG/CHH extraction with context-specific transactional bedGraph
or methylKit outputs, positional M-bias TSV/SVG reports, standalone context
merging, and per-alignment CpG metrics. Standard counts, fractions, total
depth, logit, and methylKit representations are available. Complementary CpG
and CHG strand calls can be merged before minimum-depth filtering. Indexed
1-based inclusive region selection is shared by extraction, M-bias, and
per-read reporting. Extraction and M-bias share strand- and read-specific
inclusion bounds and fixed-end trimming. Cytosine-report, BED selection,
conversion, and variant surfaces remain absent until their full behavior and
compatibility gates are implemented.

```console
rsomics-methyl extract reference.fa alignments.bam --output-prefix sample --region chr1:1-1000000
rsomics-methyl mbias reference.fa alignments.bam --output-prefix sample
rsomics-methyl extract reference.fa alignments.bam --output-prefix trimmed --OT 5,100,1,96
rsomics-methyl merge-context reference.fa sample_CpG.bedGraph --output merged.bedGraph
rsomics-methyl per-read reference.fa alignments.bam --output reads.tsv
```

## License

Licensed under either Apache License 2.0 or MIT, at your option. MethylDackel
behavior provenance and adapted-code attribution are recorded in
`THIRD_PARTY.md`.
