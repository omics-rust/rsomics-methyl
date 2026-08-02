# rsomics-methyl

`rsomics-methyl` is the rsomics product for bisulfite-sequencing methylation
extraction and bias QC. It consolidates the team-owned historical
`rsomics-methyldackel` implementation into one product family with
MethylDackel compatibility gates.

The source is under active reconstruction. The command tree currently exposes
checked CpG/CHG/CHH extraction with context-specific transactional bedGraph
outputs and standalone context merging. Bias, per-read, alternative-output,
region, trimming, conversion, and variant surfaces remain absent until their
full behavior and compatibility gates are implemented.

```console
rsomics-methyl extract reference.fa alignments.bam --output-prefix sample
rsomics-methyl merge-context reference.fa sample_CpG.bedGraph --output merged.bedGraph
```

## License

Licensed under either Apache License 2.0 or MIT, at your option. MethylDackel
behavior provenance and adapted-code attribution are recorded in
`THIRD_PARTY.md`.
