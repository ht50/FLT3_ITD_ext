# FLT3_ITD_ext
FLT3_ITD_ext is a perl script for detecting, annotating, and quantifying FLT3-ITDs in paired-end NGS data.

## Version 0.01
First commit.
Type "perl FLT3_ITD_ext.pl" for usage and options 
Input: bamfile or fastq files (R1 and R2)
Assumptions:
1. Adapters should be pre-trimmed from reads, and umis already collapsed and removed
2. Bamfiles should be relative to hg19 if extracting reads from FLT3 locus only
3. This version most applicable to hybrid capture data, e.g. need to add primer-specific logic and corrections for AMP or TSCA data.

## Version 1.0
Second commit.
Changes:
1. Platforms supported (flag -n): HC (hybrid-capture), amplicon, Archer*, NEB* (*assumes fasta and bwa index files of the primers) 
2. Three read populations may be processed from a bam (flag -t): "targeted" (default: alignments to FLT3 target locus), "loose" (also includes unmapped reads), or "all". 
3. Trims illumina adapters by default (AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT). May turn off trimming (flag: -a 0), but does not support other adapters.
4. UMI-aware bams from fgbio supported (flag -u "RX").
5. Bams relative to hg38 supported (flag -g hg38).
6. Following extraneous software tools are used and assumed in the path (otherwise script should be modified): bwa, samtools, sumaclust, fgbio, bbduk, bedtools.

## License
[MIT](https://choosealicense.com/licenses/mit/)
