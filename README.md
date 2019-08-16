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

## License
[MIT](https://choosealicense.com/licenses/mit/)
