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
1. Platforms supported (flag -n): HC (hybrid-capture), amplicon, Archer*, NEB* (*assumes fasta and bwa index files of the primers).
2. Three read populations may be processed from a bam (flag -t): "targeted" (default: alignments to FLT3 target locus), "loose" (also includes unmapped reads), or "all". 
3. Trims illumina adapters by default (AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT). May turn off trimming (flag: -a 0), but does not support other adapters.
4. UMI-aware bams from fgbio supported (flag -u "RX").
5. Bams relative to hg38 supported (flag -g hg38).

## Version 1.1
Changes:
1. Argument for fgbio umi strategy (flag -s), e.g. to support duplex umi bams
2. Minreads argument (flag -mr)
3. Conservative filters hardcoded to remove spurious ITD candidates from accidental extension (reported instead in "other_summary" output file)
4. Remove filtering of non-exonic ITDs but report in "other summary"
5. Handle more types of ITDs in VCF and fix adjustment check
6. Use R2 for readlengths in Archer
7. Improve handling of complex ITDs
8. Fix for missing MC tags in WT and Mut BAMs (thanks to jnktsj)

Before running perl script for first time:
1. Create a bwa index for the FLT3 target locus from the provided fasta file ("bwa index -p FLT3_dna_e1415 FLT3_dna_e14e15.fa") and modify script to provide path of the index if necessary (variable $refindex).
2. Download following software tools as needed and add to path (or modify script with paths): bwa, samtools, sumaclust, fgbio, bbduk, bedtools, java.

## License
[MIT](https://choosealicense.com/licenses/mit/)
