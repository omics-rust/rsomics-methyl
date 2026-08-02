# Extraction fixture

`generate_fixture.rs` creates a deterministic WGBS-like reference on the path
given as its first argument and writes coordinate-sorted SAM to standard
output. The record count is the number of alignments in `single` and
`targeted` modes and the number of fragments in `paired` mode. `targeted`
creates 512-fold, 100-base islands separated by one kilobase and includes CG,
CHG, and CHH reference contexts.

Keep generated references, alignments, outputs, and build artifacts on the
configured external volumes.

```console
rustc -O benchmarks/generate_fixture.rs -o "$TMPDIR/rsomics-methyl-fixture"
"$TMPDIR/rsomics-methyl-fixture" reference.fa 4000000 single \
  | samtools view -@ 1 -b -o alignments.bam -
samtools faidx reference.fa
samtools index alignments.bam
```

Benchmark direct release binaries at one thread, write both products to the
same output volume, and compare data-row hashes after every run. Alternate
which binary runs first within each pair so cache, temperature, and storage
drift do not consistently favor one implementation.
