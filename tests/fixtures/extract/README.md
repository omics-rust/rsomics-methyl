# Extraction fixture

`synthetic.sam` and `synthetic.fa` define this project-owned fixture. The BAM,
BAI, and FAI files were generated with samtools 1.24. The expected standard,
fraction, count, logit, and methylKit outputs were produced by MethylDackel
revision `3c77bda12141e99d80234d416e668a90ec70b3f7` and exclude their
path-dependent headers. Forward and reverse-strand calls also exercise merged
context output and post-merge minimum-depth filtering.
