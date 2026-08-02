# Extraction fixture

`synthetic.sam` and `synthetic.fa` define this project-owned fixture. The BAM,
BAI, and FAI files were generated with samtools 1.24. The expected standard,
fraction, count, logit, and methylKit outputs were produced by MethylDackel
revision `3c77bda12141e99d80234d416e668a90ec70b3f7` and exclude their
path-dependent headers. Forward and reverse-strand calls also exercise merged
context output and post-merge minimum-depth filtering. The fixed-end trimming
goldens use `--nOT 5,1,1,1` against the same revision. Inclusion-bound tests
follow its documented 1-based inclusive contract; the current upstream code
incorrectly excludes the first requested position. Exhaustive cytosine-report
goldens enable all three contexts; the zero-coverage variant uses `-q 61` to
filter every alignment.

The BED selection goldens use `[4,10)` in unstranded, top-only, and bottom-only
forms. Extraction and M-bias outputs match the same live MethylDackel revision
byte for byte after removing its path-dependent bedGraph header. Its exhaustive
cytosine report leaks zero-coverage sites outside the requested BED, and its
per-read path only checks processing chunks rather than individual alignments.
The corresponding rsomics goldens instead apply BED selection to every reported
cytosine and require each per-read alignment span and bisulfite strand to match.
