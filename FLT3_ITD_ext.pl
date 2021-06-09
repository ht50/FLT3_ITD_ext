#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(sum);
use List::Util qw(max);
use List::Util qw(min);
use Data::Dumper qw(Dumper);
use POSIX qw(ceil);
use POSIX qw(floor);
use Getopt::Long qw(GetOptions);
use Getopt::Long qw(HelpMessage);
use Pod::Usage;

=pod

=head1 NAME

FLT3_ITD_ext  Process bamfiles or fastqs for FLT3-ITDs

=head1 SYNOPSIS

  --bam, -b	Input bamfile (either this or fastq1+2 required)
  --typeb, -t	Reads to extract from input bam (defaults to "targeted" [FLT3-aligned]; or can be "loose" or "all")
  --fastq1, -f1	Input fastq1 (either fastq1+2 or bam required)
  --fastq2, -f2	Input fastq2 (either fastq1+2 or bam required)
  --output, -o  Output path (required)
  --ngstype, -n NGS platform type (defaults to "HC" [hybrid capture]; or can be "amplicon", "NEB", or "Archer")
  --genome, -g  Genome build (defaults to "hg19"; or can be "hg38")
  --adapter, -a	Trim adapters (defaults to true; assumes illumina)
  --web, -w	Create html webpages for each ITD call (defaults to false)
  --umitag, -u  BAM tag holding UMIs in the input bamfile for fgbio (defaults to ""; standard is "RX")
  --strat, -s   Strategy for UMI assignment used in fgbio GroupReadsByUmi (defaults to "adjacency" )
  --probes, -p  Probes/baits file basename (defaults to ""); assumes fasta file, bwa indexfiles
  --minreads, -mr  Minimum number of supporting reads to be included in VCF (umi-based if umitag set)
  --debug, -d	Save all intermediate files (defaults to false)
  --help, -h	Print this help

=head1 VERSION

1.1

=cut

GetOptions(
  "bam|b=s" => \(my $inbam = "" ),
  "typeb|t=s" => \(my $typeb = "targeted"),
  "fastq1|f1=s" => \(my $fastq1 = ""),
  "fastq2|f2=s" => \(my $fastq2 = ""),
  "output|o=s" => \ my $outpath,
  "ngstype|n=s" => \(my $ngstype = "HC"),
  "genome|g=s" => \(my $genome = "hg19"),
  "adapter|a=s" => \(my $adapter = 1),
  "web|w=s" => \(my $web = 0),
  "umitag|u=s" => \(my $umitag = ""),
  "strat|s=s" => \(my $strat = "adjacency"),
  "probes|p=s" => \(my $probes = ""),
  "minreads|mr=s" => \(my $minreads = 0),
  "debug|d" => \(my $debug = 0),
  "help|h" => sub { HelpMessage(0) }, 
) or HelpMessage(1);

HelpMessage(1) unless( $outpath && (length($inbam)>0 || (length($fastq1)>0 && length($fastq2)>0)) );
if( $typeb ne "targeted" && $typeb ne "loose" && $typeb ne "all" ) {
    print "Unrecognized typeb (must be targeted, loose, or all)\n";
    HelpMessage(1);
}
if( $genome ne "hg19" && $genome ne "hg38" ) {
    print "Unrecognized genome build (must be hg19 or hg38)\n";
    HelpMessage(1);
}
if( $ngstype ne "HC" && $ngstype ne "amplicon" && $ngstype ne "NEB" && $ngstype ne "Archer" ) {
    print "Unrecognized ngstype (must be HC, amplicon, NEB, or Archer)\n";
    HelpMessage(1);
}
if( $ngstype eq "HC" && $probes ne "" ) { print "Ignoring probes for HC...\n"; $probes = ""; }
if( $ngstype eq "NEB" && $probes eq "" ) { print "Need probes for NEB...\n"; HelpMessage(1); }
if( $ngstype eq "Archer" && $probes eq "" ) { print "Need probes for Archer...\n"; HelpMessage(1); }
if( $inbam eq "" && $umitag ne "" ) { print "Need bamfile with umi tags...\n"; HelpMessage(1); }

# CHANGE THIS TO LOCATION OF FLT3 INDEX (FLT3_dna_e14e15.fa, etc)
my $refindex = "FLT3_dna_e14e15"; 

# THESE COMMANDS OR FILES NEED TO BE IN $PATH; OTHERWISE CHANGE TO THEIR FULL PATHNAME 
my $samcmd = "samtools";
my $bedcmd = "bamToFastq";
my $bwacmd = "bwa";
my $clustercmd = "sumaclust";
my $javacmd = "java"; # the java cmd needs to be in $PATH (e.g. "/apps/java/jdk1.8.0_191/bin/java";)
my $fgbiojar = "fgbio.jar";
my $trimcmd = "bbduk.sh";
my $picardjar = "picard.jar";

my $minToCluster = 1;
my $maxEditDist = 5;  
my $maxInsAllowed = 2; # cannot have more than this in total ins to count as wt or itd
my $maxDelAllowed = 2; # cannot have more than this in total del to count as wt or itd
my $maxClipAllowed = 3; # cannot have more than this in max softclip to count as wt or itd
my $minClipToAlign = 9; 
my $manualItdCheckLength = 12; 
my $buffer = 10;  # buffer size past breakpoint needed for alignment to be included in counts
my $clipMatchRatioThreshold = 0.5; # ratio of clip matches to clip length needs to be at least this value in order to extend the clip
my $maxITDsize = 500;

# Conservative thresholds for filtering out ITD noise (e.g. from accidental extension)
my $noisyMinItdSize = 200;
my $noisyMinNocovRatio = 0.5; # ratio of the full duplicated sequence that is uncovered by mutant reads

# May be able to merge these, however staggered approach generates less false positives
my $bwaSeedLength = 15;
my $bwaThreshold = 12;
my $bwaClipsSeedLength = 6;
my $bwaClipsThreshold = 9;

#### reference FLT3 sequence
#### 1500bp centered around exons 14-15 (with 586bp buffer): chr13_28607438_28608937_minus_strand under hg19
my $refstart = 28607438; my $refend = 28608937; my $targetregion = "";  
if( $genome eq "hg38" ) { $refstart = 28033301; $refend = 28034800; }
my $reffasta = $refindex . ".fa"; my $refseq = ""; 
open (FI, $reffasta ) or die $!;
while(<FI>) { if( !/^>/ ) { chomp; $refseq .= $_; } }
close FI;

my $alignerLocal  = "bwa";  # "bwa" or "novo" or "bowtie2"
my $alignerClips  = "bwa";  # used for aligning softclips to reference
my $alignerExts   = "bwa";  # used for aligning extended reads to reference (bowtie2 is better at global)
my $alignerMiddle = "bwa";  # used for aligning middle portion of extended reads to reference (local alignment)
my $alignerToITDs = "bwa";  # used for aligning original reads to ITD candidates (bowtie2 is better at enforcing no-mixed)

my $mdcell = 12;  # for bwa ( $mdcell = 17 for bowtie2 )

my $proc_manual_clipalign = 0;
my $proc_index = 0;
my $proc_one_itd_per_size = 0;
my $proc_singleread_itds = 0;
my $proc_anylength_itds = 0;

my @line_cells; my @sub_cells; my @s_cells; my @a;
my $key; my $readname; my $cigar; my $seq; my $qa; my $i; my $j; my $k;
my %readlengths = ();

# protein sequence: p.534-647 in e13-e15 (full codons only)
my $aaseq = "PFPFIQDNISFYATIGVCLLFIVVLTLLICHKYKKQFRYESQLQMVQVTGSSDNEYFYVDFREYEYDLKWEFPRENLEFGKVLGSGAFGKVMNATAYGISKTGVSIQVAVKMLK";

# rna transcript c.1600-1941 in e13-e15 (full codons only)
my $rnaseq = "CCCTTCCCTTTCATCCAAGACAACATCTCATTCTATGCAACAATTGGTGTTTGTCTCCTCTTCATTGTCGTTTTAACCCTGCTAATTTGTCACAAGTACAAAAAGCAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGGAAGGTACTAGGATCAGGTGCTTTTGGAAAAGTGATGAACGCAACAGCTTATGGAATTAGCAAAACAGGAGTCTCAATCCAGGTTGCCGTCAAAATGCTGAAA";

# hash of intronic sequences by coordinates of 3'end of preceding exons
my %intronhash = (
    "1704" => "GTAAAAGCAAAGGTAAAAATTCATTATTCTTTCCTCTATCTGCAGAACTGCCTATTCCTAACTGACTCATCATTTCATCTCTGAAG",  # i13
    "1837" => "GTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTTGTCTTTGCAG",  # i14
    "1942" => "GTACAGTATAGTGGAAGGACAGCAACAAAGATGCACAAAAATGGGAGGCACAGTTTCCCACCCATGCCTTCTTCTCTTTTCCATCCTTTTAATGGTTACTGTTTGCCATGTTTCAAGGCTAAAAATGGAGTGGATTGGGGTGTCAACCCAGTCATGAATACAAAATTCAAGTCATAAATACAAAACTCCACCTCTTCTTCCACTTCTCTCTTTGTTTATTTTTATTTTCTTTTATTTTTTTGAGACAGGGTCTCACTCTGTCACCTAGGCTGGAGTACAGTAGTGTGATCATGGCTCACTGCAGCCTCCACCCCCTGGGCTCAGGCTCCTCAGCAGCTGGGGCTACAGGCGTGCACCACCACGCTCAGCTAAGTTTTTTATTTTTGGTAGAGACGGGGTTTCACCATGTTGCCTAGGCTGATCTGGGACTCCTGAGCTCAAACTATCTGCCCGCCTCGGCCTCCCAAAGTGTTGGGATTACAGGCGTGAGCCACTGTGCCAGGCCTATTTAGTTTTTATAGGCTTTGATATTTGGGGATACAGTTGCCACAGGCAGAAGTTCTTTCTTACATTTGTGTTGTTTTAT"  # i15
    );

my $offrna = 1600; my $offaa = 534;

my $vcfheader = "##fileformat=VCFv4.2\n" .
    "##INFO=<ID=GENE,Number=1,Type=String,Description=\"Gene name\">\n" .
    "##INFO=<ID=STRAND,Number=1,Type=String,Description=\"Gene strand\">\n" .
    "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n" .
    "##INFO=<ID=CDS,Number=1,Type=String,Description=\"CDS annotation (modified)\">\n" .
    "##INFO=<ID=AA,Number=1,Type=String,Description=\"Peptide annotation (modified)\">\n" .
    "##INFO=<ID=AR,Number=1,Type=Float,Description=\"Allelic Ratio\">\n" .
    "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">\n" .
    "##INFO=<ID=DP,Number=1,Type=Float,Description=\"Total Depth\">\n" .
    "##INFO=<ID=VD,Number=1,Type=Float,Description=\"Variant Depth\">\n" .
    "##INFO=<ID=AAR,Number=1,Type=Float,Description=\"Adjusted Allelic Ratio (Max Estimate)\">\n" .
    "##INFO=<ID=AAF,Number=1,Type=Float,Description=\"Adjusted Allele Frequency (Max Estimate)\">\n" .
    "##INFO=<ID=RAR,Number=1,Type=Float,Description=\"Raw Allelic Ratio\">\n" .
    "##INFO=<ID=RAF,Number=1,Type=Float,Description=\"Raw Allele Frequency\">\n" .
    "##INFO=<ID=RDP,Number=1,Type=Float,Description=\"Raw Total Depth\">\n" .
    "##INFO=<ID=RVD,Number=1,Type=Float,Description=\"Raw Variant Depth\">\n" .
    "##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Sample name\">\n" .
    "##FORMAT=<ID=CR,Number=6,Type=Integer,Description=\"Coverage Radius\">\n" .
    "##FORMAT=<ID=DB,Number=.,Type=String,Description=\"Duplicated Baits (among fwdP1-P3, revP1-P3, other)\">\n" .
    "##FORMAT=<ID=NB,Number=.,Type=String,Description=\"Non-binding Baits\">\n" .
    "##FORMAT=<ID=RB,Number=.,Type=String,Description=\"Reaching Baits capturing mut + wt junctions in R2 (unambiguous)\">\n" .
    "##FORMAT=<ID=AFB,Number=7,Type=Float,Description=\"Allele frequency by Baits\">\n" .
    "##FORMAT=<ID=Fwd1,Number=1,Type=String,Description=\"AF details for bait Fwd1\">\n" .
    "##FORMAT=<ID=Fwd2,Number=1,Type=String,Description=\"AF details for bait Fwd2\">\n" .
    "##FORMAT=<ID=Fwd3,Number=1,Type=String,Description=\"AF details for bait Fwd3\">\n" .
    "##FORMAT=<ID=Rev1,Number=1,Type=String,Description=\"AF details for bait Rev1\">\n" .
    "##FORMAT=<ID=Rev2,Number=1,Type=String,Description=\"AF details for bait Rev2\">\n" .
    "##FORMAT=<ID=Rev3,Number=1,Type=String,Description=\"AF details for bait Rev3\">\n" .
    "##FORMAT=<ID=Oth,Number=1,Type=String,Description=\"AF details for reads not assigned to a bait\">\n" .
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";


###################################
# REFERENCE DETAILS VERSUS hg19
#
# WT_ref
# start: 28608937
# end: 28607438
# ref coords: 1-1500
#
# Exon13: (length 107)
# start: 28608544
# end: 28608438
# cdot: c.1598-1704
# pdot: p.533-568 (codon2-codon3)
# ref coords: 394-500
#
# Intron13: (length 86)
# start: 28608437
# end: 28608352
# ref coords: 501-586
#
# Exon14: (length 133)
# start: 28608351
# end: 28608219
# cdot: c.1705-1837
# pdot: p.569-613 (codon1-codon1)
# ref coords: 587-719
#
# Intron14: (length 90)
# start: 28608218
# end:28608129
# ref coords: 720-809
#
# Exon15: (length 105)
# start:28608128
# end:28608024
# cdot: c.1838-1942
# pdot: p.613-648 (codon2-codon1)
# ref coords: 810-914
#
##################################

sub increment_cdot_coord {
    # argument #1: cdot coordinate (could be exonic "1837" or intronic "1203+20") 
    my $rpos = $_[0];

    my %inchash = (
	"1597+1087" => "1598", # end of i12
	"1704" => "1704+1",  # end of e13
	"1704+86" => "1705", # end of i13
	"1837" => "1837+1",  # end of e14
	"1837+90" => "1838", # end of i14
	"1942" => "1942+1",  # end of e15
	);

    if( exists( $inchash{$rpos} ) ) {
	return( $inchash{$rpos} );
    } else {
	my @cells = split(/\+/, $rpos);
	if( scalar(@cells) == 1 ) {
	    return( $rpos+1 );
	} else {
	    return( $cells[0] . "+" . ($cells[1]+1) );
	}
    }
}

sub decrement_cdot_coord {
    # argument #1: cdot coordinate (could be exonic "1837" or intronic "1203+20") 
    my $rpos = $_[0];

    my %dechash = (
	"1598" => "1597+1087", # start of e13
	"1704+1" => "1704",  # start of i13
	"1705" => "1704+86", # start of e14
	"1837+1" => "1837",  # start of i14
	"1838" => "1837+90", # start of e15
	"1942+1" => "1942",  # start of i15
	);

    if( exists( $dechash{$rpos} ) ) {
	return( $dechash{$rpos} );
    } else {
	my @cells = split(/\+/, $rpos);
	if( scalar(@cells) == 1 ) {
	    return( $rpos-1 );
	} else {
	    return( $cells[0] . "+" . ($cells[1]-1) );
	}
    }
}

sub proc_to_cdot_coords {
    # argument #1: reference coordinate
    my $rpos = $_[0];

    my @iestarts = (1, 394, 501, 587, 720, 810, 915, 1500);
    my @ielabs = ("i12", "e13", "i13", "e14", "i14", "e15", "i15" );
    my %eoffs = ( "e13" => 1204, "e14" => 1118, "e15" => 1028 );
    my %ioffs = ( "i12" => 694, "i13" => -500, "i14" => -719, "i15" => -914 );
    my %ibases = ( "i12" => 1597, "i13" => 1704, "i14" => 1837, "i15" => 1942 );

    my $ie = -1;
    for( my $i=0; $i<scalar(@iestarts)-1; $i++ ) {
	if( $rpos >= $iestarts[$i] && $rpos < $iestarts[$i+1] ) { $ie = $i; }
    }
    if( $rpos == 1500 ) { $ie = scalar(@ielabs)-1; }

    if( $ie >= 0 ) {
	if( substr($ielabs[$ie],0,1) eq "e" ) {
	    return( ( $rpos + $eoffs{$ielabs[$ie]}, $ielabs[$ie] ) );
	} else {
	    return( ( $ibases{$ielabs[$ie]} . "+" . ($rpos+$ioffs{$ielabs[$ie]}), $ielabs[$ie] ) );
	}
    } else {
	return( ( "", "" ) );
    }	    
}

sub proc_to_aa {
    # argument #1: DNA sequence
    my $dna = $_[0];
    if( length($dna) % 3 != 0 ) { return(""); }  # needs to be multiple of 3

    my %aahash = (
	TTT => "F", TTC => "F", TTA => "L", TTG => "L",
	CTT => "L", CTC => "L", CTA => "L", CTG => "L",
	ATT => "I", ATC => "I", ATA => "I", ATG => "M",
	GTT => "V", GTC => "V", GTA => "V", GTG => "V",
	TCT => "S", TCC => "S", TCA => "S", TCG => "S",
	CCT => "P", CCC => "P", CCA => "P", CCG => "P", 
	ACT => "T", ACC => "T", ACA => "T", ACG => "T", 
	GCT => "A", GCC => "A", GCA => "A", GCG => "A",
	TAT => "Y", TAC => "Y", TAA => "*", TAG => "*", 
	CAT => "H", CAC => "H", CAA => "Q", CAG => "Q",
	AAT => "N", AAC => "N", AAA => "K", AAG => "K",
	GAT => "D", GAC => "D", GAA => "E", GAG => "E",
	TGT => "C", TGC => "C", TGA => "*", TGG => "W", 
	CGT => "R", CGC => "R", CGA => "R", CGG => "R", 
	AGT => "S", AGC => "S", AGA => "R", AGG => "R", 
	GGT => "G", GGC => "G", GGA => "G", GGG => "G" 
	);

    my $aa = "";
    while( length($dna) > 0 ) {
	if( !exists($aahash{substr($dna,0,3)}) ) { return(""); }
	$aa .= $aahash{substr($dna,0,3)};
	if( $aahash{substr($dna,0,3)} eq "*" ) { return($aa); }
	$dna = substr($dna,3);
    }
    return($aa);
}

sub proc_subkey_to_altseq {
    # argument #1: substitutions key
    my $subkey = $_[0];
    $subkey =~ s/ref\_//g;
    $subkey =~ s/c\.//g;
    $subkey =~ s/c1\.//g;
    $subkey =~ s/c2\.//g;
    $subkey =~ s/\(//g;
    $subkey =~ s/\)//g;

    # argument #2: alt start coord
    my $altstart = $offrna;
    if( scalar(@_) > 1 ) { $altstart = $_[1]; }

    # argument #3: alt end coord
    my $altend = length($rnaseq)+$offrna-1;
    if(scalar(@_) > 2 ) { $altend = $_[2]; }

    my $altseq = $rnaseq;
    my %althash = (); foreach( keys %intronhash ) { $althash{$_} = $intronhash{$_}; }

    if( length($subkey) > 0 ) {
	my @subcells = split(/;/, $subkey);
	for(my $i=0; $i<scalar(@subcells); $i++ ) {
	    my @subsubcells = split(/:/, $subcells[$i]);
	    my @subsubsubcells = split(/>/, $subsubcells[1]);
	    if( !($subsubcells[0] =~ /\+/) ) {
		substr($altseq, $subsubcells[0]-$offrna, 1 ) = $subsubsubcells[1];
	    } else {
		my @subsubsubsubcells = split( /\+/, $subsubcells[0] );
		if( exists( $althash{$subsubsubsubcells[0]} ) ) {
		    substr( $althash{$subsubsubsubcells[0]}, $subsubsubsubcells[1]-1, 1 ) = $subsubsubcells[1]; 
		}
	    }
	}
    }

    my @subcells0 = split(/\+/, $altstart);
    my @subcells1 = split(/\+/, $altend);
    $altseq = substr($altseq, $subcells0[0]-$offrna, $subcells1[0]-$subcells0[0]+1);

    if( scalar(@subcells0)>1 && exists($althash{$subcells0[0]} )) {
	$altseq = substr($althash{$subcells0[0]}, $subcells0[1]-1) . $altseq;
    }

    if( scalar(@subcells1)>1 && exists($althash{$subcells1[0]} )) {
	$altseq .= substr($althash{$subcells1[0]}, 0, $subcells1[1]);
    }

    return( $altseq );
}

sub proc_subkey_to_subhash {
    # argument #1: substitutions key
    my $subkey = $_[0];
    $subkey =~ s/ref\_//g;
    $subkey =~ s/c\.//g;
    $subkey =~ s/c1\.//g;
    $subkey =~ s/c2\.//g;
    $subkey =~ s/\(//g;
    $subkey =~ s/\)//g;

    my %subhash = ();
    my @subcells = split(/;/, $subkey);
    for(my $i=0; $i<scalar(@subcells); $i++ ) {
	my @subsubcells = split(/:/, $subcells[$i]);
	my @subsubsubcells = split(/>/, $subsubcells[1]);
	$subhash{$subsubcells[0]} = $subsubsubcells[1];
    }

    return %subhash;
}

sub proc_cdot_to_pdot {
    # argument #1: cdot key
    my $cdot = $_[0];
    my @subcells; my @subsubcells; my @subsubsubcells; my @subsubsubsubcells;
    my $pdot = ""; my $aa = ""; my $subkey = "";

    my @cells = split( /\|/, substr($cdot,2) ); # strip leading c.
    my %subhash = (); 
    if( scalar(@cells)>1 ) { %subhash = proc_subkey_to_subhash($cells[1]); }

    # turn dup into ins nomenclature if not a perfect dup (underlying mutations)
    @subcells = split( /\(/, $cells[0] );
    if( $cells[0] =~ /dup/ && scalar(@subcells) > 1 ) {
	my %subhash1 = (); my %subhash2 = (); my $subkey1 =""; my $subkey2 = "";
	for( my $i=1; $i<scalar(@subcells); $i++ ) {
	    if( $subcells[$i] =~ /c1/ ) {
		%subhash1 = proc_subkey_to_subhash($subcells[$i]);
		$subkey1 = $subcells[$i];
		$subkey1 =~ s/c1/c/g;
		$subkey1 =~ s/\(//g;
	    } elsif( $subcells[$i] =~ /c2/ ) {
		%subhash2 = proc_subkey_to_subhash($subcells[$i]);
		$subkey2 = $subcells[$i];
		$subkey2 =~ s/c2/c/g;
		$subkey2 =~ s/\(//g;
	    }
	}

	$cells[0] = "";
	if( scalar( keys %subhash1 ) > 0 || scalar( keys %subhash2 ) > 0 ) {
	    $subcells[0] =~ s/dup//g;
	    $subcells[0] =~ s/\[\+i14\]//g;
	    @subsubcells = split( /\_/, $subcells[0] );
	    if( ( $subsubcells[0] =~ /\+/ || scalar( keys %subhash1 ) <= scalar( keys %subhash2 ) ) &&
		!($subsubcells[1] =~ /\+/) ) {
		foreach( keys %subhash1 ) { $subhash{$_} = $subhash1{$_}; }
		$cells[0] = ($subsubcells[1]) . "_" . increment_cdot_coord($subsubcells[1]) . 
		    "ins" . $subsubcells[0] . "-" . $subsubcells[1];
		if( length($subkey2)>0 ) { $cells[0] .= "(" . $subkey2 . ")"; }
	    } elsif( ( $subsubcells[1] =~ /\+/ || scalar( keys %subhash1 ) > scalar( keys %subhash2 ) ) &&
		!($subsubcells[0] =~ /\+/) ) {
		foreach( keys %subhash2 ) { $subhash{$_} = $subhash2{$_}; }
		$cells[0] = decrement_cdot_coord($subsubcells[0]) . "_" . ($subsubcells[0]) . 
		    "ins" . $subsubcells[0] . "-" . $subsubcells[1];
		if( length($subkey1)>0 ) { $cells[0] .= "(" . $subkey1 . ")"; }
	    }
	}
    }

    # process dup nomenclatures
    if( $cells[0] =~ /dup/ ) {  
	@subcells = split(/dup/, $cells[0]);
	@subsubcells = split(/\_/, $subcells[0] );

	# duplication of region within exons
	if( !($subsubcells[0] =~ /\+/) && !($subsubcells[1] =~ /\+/) ) {
	    if( ($subsubcells[0] % 3) == 1 && ($subsubcells[1] % 3) == 0 ) {
		$pdot = "p." . ($subsubcells[0]+2)/3 . "_" . ($subsubcells[1]/3) . "dup";
	    } elsif( ($subsubcells[0] % 3) == 2 && ($subsubcells[1] % 3) == 1 ) {
		$aa = proc_to_aa( substr($rnaseq, $subsubcells[1]-$offrna, 1) . substr($rnaseq, $subsubcells[0]-$offrna, 2) );
		if( substr($aaseq, ($subsubcells[1]+2)/3-$offaa, 1) eq $aa ) {
		    $pdot = "p." . ($subsubcells[0]+4)/3 . "_" . ($subsubcells[1]+2)/3 . "dup";
		} elsif( substr($aaseq, ($subsubcells[0]+1)/3-$offaa, 1) eq $aa ) {
		    $pdot = "p." . ($subsubcells[0]+1)/3 . "_" . ($subsubcells[1]-1)/3 . "dup";
		} elsif( $aa eq "*" ) {
		    $pdot = "p." . ($subsubcells[1]-1)/3 . "_" . ($subsubcells[1]+2)/3 . "ins" . $aa;
		} else {
		    $pdot = "p." . ($subsubcells[1]-1)/3 . "_" . ($subsubcells[1]+2)/3 . "ins" . $aa . 
			"/" . ($subsubcells[0]+4)/3 . "-" . ($subsubcells[1]-1)/3;
		}
	    } elsif( ($subsubcells[0] % 3) == 0 && ($subsubcells[1] % 3) == 2 ) {
		$aa = proc_to_aa( substr($rnaseq, $subsubcells[1]-1-$offrna, 2) . substr($rnaseq, $subsubcells[0]-$offrna, 1) );
		if( substr($aaseq, ($subsubcells[1]+1)/3-$offaa, 1) eq $aa ) {
		    $pdot = "p." . ($subsubcells[0]+3)/3 . "_" . ($subsubcells[1]+1)/3 . "dup";
		} elsif( substr($aaseq, $subsubcells[0]/3-$offaa, 1) eq $aa ) {
		    $pdot = "p." . $subsubcells[0]/3 . "_" . ($subsubcells[1]-2)/3 . "dup";
		} elsif( $aa eq "*" ) {
		    $pdot = "p." . ($subsubcells[1]-2)/3 . "_" . ($subsubcells[1]+1)/3 . "ins" . $aa;
		} else {
		    $pdot = "p." . ($subsubcells[1]-2)/3 . "_" . ($subsubcells[1]+1)/3 . "ins" . $aa .
			"/" . ($subsubcells[0]+3)/3 . "-" . ($subsubcells[1]-2)/3;
		}
	    }
	# duplication of region within an exon to an intron
	} elsif( !($subsubcells[0] =~ /\+/ ) && $subsubcells[1] =~ /\+/ ) {
	    my @subsubsubcells = split(/\+/, $subsubcells[1]);
	    if( exists($intronhash{$subsubsubcells[0]}) ) {
		my $dna = substr($rnaseq, $subsubcells[0]-$offrna, $subsubsubcells[0]-$subsubcells[0]+1) .
		    substr( $intronhash{$subsubsubcells[0]}, 0, $subsubsubcells[1] );
		if( (length($dna) % 3) == 0 && ($subsubcells[0] % 3 ) == 1 ) {
		    $pdot= "p." . ($subsubcells[0]-1)/3 . "_" . ($subsubcells[0]+2)/3 . "ins" . proc_to_aa( $dna );
		} elsif( (length($dna) % 3) == 0 && ($subsubcells[0] % 3 ) == 2 ) {
		    # may be off when $subsubcells[0] is at exon/intron junction (but only an issue for e13i13)
		    $pdot= "p." . ($subsubcells[0]+1)/3 . "_" . ($subsubcells[0]+4)/3 . "ins" . 
			proc_to_aa( substr($dna, 2) . substr($dna, 0, 2) );
		} elsif( (length($dna) % 3) == 0 && ($subsubcells[0] % 3 ) == 0 ) {
		    $pdot= "p." . $subsubcells[0]/3 . "_" . ($subsubcells[0]+3)/3 . "ins" . 
			proc_to_aa( substr($dna, 1) . substr($dna, 0, 1) );
		}
	    }
	# duplication of region within an intron to an exon
	} elsif( $subsubcells[0] =~ /\+/ && !($subsubcells[1] =~ /\+/) ) {
	    my @subsubsubcells = split(/\+/, $subsubcells[0]);
	    my $dna = substr( $intronhash{$subsubsubcells[0]}, $subsubsubcells[1]-1 ) .
		substr($rnaseq, $subsubsubcells[0]+1-$offrna, $subsubcells[1]-$subsubsubcells[0]);
	    if( exists($intronhash{$subsubsubcells[0]}) ) {
		if( (length($dna) % 3) == 0 && ($subsubcells[1] % 3 ) == 0 ) {
		    $pdot= "p." . $subsubcells[1]/3 . "_" . ($subsubcells[1]+3)/3 . "ins" . proc_to_aa( $dna );
		} elsif( (length($dna) % 3) == 0 && ($subsubcells[0] % 3 ) == 1 ) {
		    $pdot= "p." . ($subsubcells[1]-1)/3 . "_" . ($subsubcells[1]+2)/3 . "ins" . 
			proc_to_aa( substr($dna, length($dna)-1, 1) . substr($dna, 0, length($dna)-1) );
		} elsif( (length($dna) % 3) == 0 && ($subsubcells[0] % 3 ) == 2 ) {
		    # may be off when $subsubcells[1] is at intron/exon junction
		    $pdot= "p." . ($subsubcells[1]-2)/3 . "_" . ($subsubcells[1]+1)/3 . "ins" . 
			proc_to_aa( substr($dna, length($dna)-2, 2) . substr($dna, 0, length($dna)-2) );
		}		    
	    }
	} elsif( $subsubcells[0] =~ /\+/ && $subsubcells[1] =~ /\+/ ) {
	    return( "not_exonic" );
	}
    # process the ins/delins nomenclature
    } elsif( $cells[0] =~ /ins/ ) {
	@subcells = split(/ins/, $cells[0]);
	@subsubcells = split(/\//, $subcells[1] );
	$subcells[0] =~ s/del//g;
	my $iseq = "";
	for( my $i=0; $i<scalar(@subsubcells); $i++ ) {
	    if( !($subsubcells[$i] =~ /\-/) ) {
		$iseq .= $subsubcells[$i];
	    } else {
		@subsubsubcells = split(/\(/, $subsubcells[$i]);
		@subsubsubsubcells = split( /\-/, $subsubsubcells[0] );
		$subkey = "";
		if( scalar(@subsubsubcells)>1 ) { $subkey = $subsubsubcells[1]; }
		$subsubsubsubcells[1] =~ s/\[\+i14\]//g;
		$iseq .= proc_subkey_to_altseq( $subkey, $subsubsubsubcells[0], $subsubsubsubcells[1] );
	    }
	}

	my $delins = ""; my $aa; my $indelaa;
	@subsubcells = split(/\_/, $subcells[0]);

	# Allow insertions at exon/intron or intron/exon boundaries by translating to exonic coords
	if( $subsubcells[1] eq ($subsubcells[0] . "+1") ) { $subsubcells[1] = $subsubcells[0]+1; }
	if( !($subsubcells[1] =~ /\+/) && exists( $intronhash{$subsubcells[1]-1} ) &&
	    $subsubcells[0] eq ($subsubcells[1]-1) . "+" . length($intronhash{$subsubcells[1]-1}) ) {
	    $subsubcells[0] = $subsubcells[1]-1;
	}

	# Only report exonic insertions
	if( $subsubcells[0] =~ /\+/ || $subsubcells[1] =~ /\+/ ) { return("not_exonic"); }

	# Boundaries of insertion (left and right coordinates)
	my $leftC = $subsubcells[0]; my $rightC = $subsubcells[1];
	if( $cells[0] =~ /delins/ ) {
	    $leftC = $subsubcells[0]-1;
	    $rightC = $subsubcells[1]+1;
	}

	# Only report in-frame delins
	if( (($rightC - $leftC + length($iseq) - 1) % 3) != 0 ) { return("frameshift"); }

	# Take the nucleotides completing the left codon
	# Mod 3 = 0 : no nucleotides needed
	# Mod 3 = 1 : 1 nucleotide needed (position 1)
	# Mod 3 = 2 : 2 nucleotides needed (positions 1 & 2)
	my $leftnt = $leftC % 3;
	my $leftaa = ceil($leftC/3);
	for( my $i = $leftnt-1; $i >= 0; $i-- ) {
	    if( exists($subhash{$leftC-$i}) ) {
		$delins .= $subhash{$leftC-$i};
		delete $subhash{$leftC-$i};
	    } else {
		$delins .= substr($rnaseq, $leftC-$i-$offrna, 1);
	    }
	}

	# Add in the insertion seq
	$delins .= $iseq;

	# Take the nucleotides completing the right codon
	# Mod 3 = 1 : no nucleotides needed
	# Mod 3 = 2 : 2 nucleotides needed (positions 2 & 3)
	# Mod 3 = 3 : 1 nucleotides needed (position 3)
	my $rightnt = 0; 
	my $rightaa = ceil($rightC/3);
	if( ($rightC % 3) == 0 ) { $rightnt = 1; } elsif( ($rightC % 3) == 2 ) { $rightnt = 2; }
	for( my $i = 0; $i < $rightnt; $i++ ) {
	    if( exists($subhash{$rightC+$i}) ) {
		$delins .= $subhash{$rightC+$i};
		delete $subhash{$rightC+$i};
	    } else {
		$delins .= substr($rnaseq, $rightC+$i-$offrna, 1);
	    }
	}

	my $aaL = substr($aaseq, $leftaa-$offaa, 1);
	my $aaR = substr($aaseq, $rightaa-$offaa, 1);
	my $aaLinc = substr($aaseq, $leftaa+1-$offaa, 1);
	my $aaRdec = substr($aaseq, $rightaa-1-$offaa, 1);
	$indelaa = proc_to_aa($delins);

	if( ($leftC % 3)==0 && ($rightC % 3)==1 ) {
	    if( $rightaa == $leftaa+1 ) {
		$pdot = "p." . $aaL . $leftaa . "_" . $aaR . $rightaa . "ins" . $indelaa;
	    } elsif( $rightaa == $leftaa+2 )  {
		$pdot = "p." . $aaLinc . ($leftaa+1) . "delins" . $indelaa;
	    } elsif( $rightaa > $leftaa+2 )  {
		$pdot = "p." . $aaLinc . ($leftaa+1) . "_" . $aaRdec . ($rightaa-1) . "delins" . $indelaa;
	    } else {
		return( "" );
	    }
	} elsif( ($leftC % 3)==0 ) {
	    if( $rightaa == $leftaa+1 ) {
		if( substr($indelaa,0,1) eq $aaR ) {
		    $pdot = "p." . $aaL . $leftaa . "_" . $aaR . $rightaa . "ins" . substr($indelaa,1);
		} elsif( substr($indelaa,length($indelaa)-1,1) eq $aaR ) {
		    $pdot = "p." . $aaL . $leftaa . "_" . $aaR . $rightaa . "ins" . substr($indelaa,0,length($indelaa-1));
		} else {
		    $pdot = "p." . $aaR . $rightaa . "delins" . $indelaa;
		}
	    } elsif( $rightaa == $leftaa+2 )  {
		if( substr($indelaa,0,1) eq $aaR ) {
		    $pdot = "p." . $aaLinc . ($leftaa+1). "delins" . substr($indelaa,1);
		} elsif( substr($indelaa,length($indelaa)-1,1) eq $aaR ) {
		    $pdot = "p." . $aaLinc . ($leftaa+1) . "delins" . substr($indelaa,0,length($indelaa-1));
		} else {
		    $pdot = "p." . $aaLinc . ($leftaa+1) . "_" . $aaR . $rightaa . "delins" . $indelaa;
		}
	    } elsif( $rightaa > $leftaa+2 )  {
		if( substr($indelaa,0,1) eq $aaR ) {
		    $pdot = "p." . $aaLinc . ($leftaa+1) . "_" . $aaRdec . ($rightaa-1) . "delins" . substr($indelaa,1);
		} elsif( substr($indelaa,length($indelaa)-1,1) eq $aaR ) {
		    $pdot = "p." . $aaLinc . ($leftaa+1) . "_" . $aaRdec . ($rightaa-1) . "delins" . substr($indelaa,0,length($indelaa-1));
		} else {
		    $pdot = "p." . $aaLinc . ($leftaa+1) . "_" . $aaR . $rightaa . "delins" . $indelaa;
		}
	    } else {
		return( "" );
	    }
	} elsif( ($rightC % 3)==1 ) {
	    if( $rightaa == $leftaa+1 ) {
		if( substr($indelaa,0,1) eq $aaL ) {
		    $pdot = "p." . $aaL . $leftaa . "_" . $aaR . $rightaa . "ins" . substr($indelaa,1);
		} elsif( substr($indelaa,length($indelaa)-1,1) eq $aaL ) {
		    $pdot = "p." . $aaL . $leftaa . "_" . $aaR . $rightaa . "ins" . substr($indelaa,0,length($indelaa-1));
		} else {
		    $pdot = "p." . $aaL . $leftaa . "delins" . $indelaa;
		}
	    } elsif( $rightaa == $leftaa+2 )  {
		if( substr($indelaa,0,1) eq $aaL ) {
		    $pdot = "p." . $aaLinc . ($leftaa+1). "delins" . substr($indelaa,1);
		} elsif( substr($indelaa,length($indelaa)-1,1) eq $aaL ) {
		    $pdot = "p." . $aaLinc . ($leftaa+1) . "delins" . substr($indelaa,0,length($indelaa-1));
		} else {
		    $pdot = "p." . $aaL . $leftaa . "_" . $aaRdec . ($rightaa-1) . "delins" . $indelaa;
		}
	    } elsif( $rightaa > $leftaa+2 )  {
		if( substr($indelaa,0,1) eq $aaL ) {
		    $pdot = "p." . $aaLinc . ($leftaa+1) . "_" . $aaRdec . ($rightaa-1) . "delins" . substr($indelaa,1);
		} elsif( substr($indelaa,length($indelaa)-1,1) eq $aaL ) {
		    $pdot = "p." . $aaLinc . ($leftaa+1) . "_" . $aaRdec . ($rightaa-1) . "delins" . substr($indelaa,0,length($indelaa-1));
		} else {
		    $pdot = "p." . $aaL . $leftaa . "_" . $aaRdec . ($rightaa-1) . "delins" . $indelaa;
		}
	    } else {
		return( "" );
	    }
	} else {
	    if( $leftaa == $rightaa ) {
		if( substr($indelaa,0,1) eq $aaL ) {
		    $pdot = "p." . $aaL . $leftaa . "_" . $aaLinc . ($leftaa+1) . "ins" . substr($indelaa,1);
		} elsif( substr($indelaa,length($indelaa)-1,1) eq $aaR ) {
		    $pdot = "p." . $aaRdec . ($rightaa-1) . "_" . $aaR . $rightaa . "ins" . substr($indelaa,0,length($indelaa)-1);
		} else {
		    $pdot = "p." . $aaL . $leftaa . "delins" . $indelaa;
		}
	    } else {
		if( substr($indelaa,0,1) eq $aaL && substr($indelaa,length($indelaa)-1,1) eq $aaR ) {
		    if( length($indelaa)==1 ) {
			if( $rightaa == $leftaa+1 ) {
			    $pdot = "p." . $aaR . $rightaa . "del";
			} elsif( $rightaa > $leftaa+1 ) {
			    $pdot = "p." . $aaLinc . ($leftaa+1) .  "_" . $aaR . $rightaa . "del";
			} else {
			    return( "" );
			}
		    } elsif( length($indelaa)==2 )  {
			if( $rightaa == $leftaa+1 ) {
			    $pdot = "p.synonymous";
			} elsif( $rightaa == $leftaa+2 ) {
			    $pdot = "p." . $aaLinc . ($leftaa+1) . "del";
			} else {
			    $pdot = "p." . $aaLinc . ($leftaa+1) . "_" . $aaRdec . ($rightaa-1) . "_" . "del";
			}
		    } else {
			if( $rightaa == $leftaa+1 ) {
			    $pdot = "p." . $aaL . $leftaa . "_" . $aaR . $rightaa . "ins" .
				substr($indelaa, 1, length($indelaa)-2);
			    $pdot = "p.synonymous";
			} elsif( $rightaa == $leftaa+2 ) {
			    $pdot = "p." . $aaLinc . ($leftaa+1) . "delins" .
				substr($indelaa, 1, length($indelaa)-2);
			} else {
			    $pdot = "p." . $aaLinc . ($leftaa+1) . "_" . $aaRdec . ($rightaa-1) . "_" . "delins" .
				substr($indelaa, 1, length($indelaa)-2);
			}
		    }
		} elsif( substr($indelaa,0,1) eq $aaL ) {
		    if( length($indelaa)==1 ) {
			if( $rightaa == $leftaa+1 ) {
			    $pdot = "p." . $aaR . $rightaa . "del";
			} else {
			    $pdot = "p." . $aaLinc . ($leftaa+1) . "_" . $aaR . $rightaa . "del";
			}
		    } else {
			if( $rightaa == $leftaa+1 ) {
			    $pdot = "p." . $aaR . $rightaa . "delins" . substr($indelaa, 1);
			} else {
			    $pdot = "p." . $aaLinc . ($leftaa+1) . "_" . $aaR . $rightaa . "delins" . substr($indelaa, 1);;
			}
		    }
		} elsif( substr($indelaa,length($indelaa)-1,1) eq $aaR ) {
		    if( length($indelaa)==1 ) {
			if( $rightaa == $leftaa+1 ) {
			    $pdot = "p." . $aaL . $leftaa . "del";
			} else {
			    $pdot = "p." . $aaL . $leftaa . "_" . $aaRdec . ($rightaa-1) . "del";
			}
		    } else {
			if( $rightaa == $leftaa+1 ) {
			    $pdot = "p." . $aaL . $leftaa . "delins" . substr($indelaa, 0, length($indelaa)-1);
			} else {
			    $pdot = "p." . $aaL . $leftaa . "_" . $aaRdec . ($rightaa-1) . "delins" .
				substr($indelaa, 0, length($indelaa)-1);
			}
		    }
		} else {
		    $pdot = "p." . $aaL . $leftaa . "_" . $rightaa . "delins" . $indelaa;
		}
	    }
	}
    # process the del nomenclature
    } elsif( $cells[0] =~ /del/ ) {
	my $delins = ""; my $aa; my $indelaa;
	@subcells = split(/del/, $cells[0]);
	@subsubcells = split(/\_/, $subcells[0]);

	# Only report exonic in-frame deletions
	if( $subsubcells[0] =~ /\+/ || $subsubcells[1] =~ /\+/ ) { return( "not_exonic" ); }
	if( (($subsubcells[1] - $subsubcells[0] + 1) % 3) != 0 ) { return( "frameshift" ); }

	# Boundaries of insertion (left and right coordinates)
	my $leftC = $subsubcells[0]-1; my $rightC = $subsubcells[1]+1;

	# Take the nucleotides completing the left codon
	# Mod 3 = 0 : no nucleotides needed
	# Mod 3 = 1 : 1 nucleotide needed (position 1)
	# Mod 3 = 2 : 2 nucleotides needed (positions 1 & 2)
	my $leftnt = $leftC % 3;
	my $leftaa = ceil($leftC/3);
	for( my $i = $leftnt-1; $i >= 0; $i-- ) {
	    if( exists($subhash{$leftC-$i}) ) {
		$delins .= $subhash{$leftC-$i};
		delete $subhash{$leftC-$i};
	    } else {
		$delins .= substr($rnaseq, $leftC-$i-$offrna, 1);
	    }
	}

	# Take the nucleotides completing the right codon
	# Mod 3 = 1 : no nucleotides needed
	# Mod 3 = 2 : 2 nucleotides needed (positions 2 & 3)
	# Mod 3 = 3 : 1 nucleotides needed (position 3)
	my $rightnt = 0; 
	my $rightaa = ceil($rightC/3);
	if( ($rightC % 3) == 0 ) { $rightnt = 1; } elsif( ($rightC % 3) == 2 ) { $rightnt = 2; }
	for( my $i = 0; $i < $rightnt; $i++ ) {
	    if( exists($subhash{$rightC+$i}) ) {
		$delins .= $subhash{$rightC+$i};
		delete $subhash{$rightC+$i};
	    } else {
		$delins .= substr($rnaseq, $rightC+$i-$offrna, 1);
	    }
	}

	$indelaa = proc_to_aa($delins);

	if( $leftnt + $rightnt == 0 ) {
	    if( $rightaa == $leftaa+1 ) {
		$pdot = "p.synonymous";
	    } elsif( $rightaa == $leftaa+2 ) {
		$pdot = "p." . substr($aaseq, $leftaa+1-$offaa, 1) . ($leftaa+1) . "del";
	    } else {
		$pdot = "p." . substr($aaseq, $leftaa+1-$offaa, 1) . ($leftaa+1) . "_" .
		    substr($aaseq, $rightaa-1-$offaa, 1) . ($rightaa-1) . "del";
	    }
	} else {
	    my $aaL = substr($aaseq, $leftaa-$offaa, 1);
	    my $aaR = substr($aaseq, $rightaa-$offaa, 1);
	    if( $indelaa eq $aaL ) {
		if( $rightaa == $leftaa+1 ) {
		    $pdot = "p." . $aaR . $rightaa . "del";
		} else {
		    $pdot = "p." . substr($aaseq, $leftaa+1-$offaa, 1) . ($leftaa+1) . "_" . $aaR . $rightaa . "del";
		}
	    } elsif( $indelaa eq $aaR ) {
		if( $rightaa == $leftaa+1 ) {
		    $pdot = "p." . $aaL . $leftaa . "del";
		} else {
		    $pdot = "p." . $aaL . $leftaa . "_" . substr($aaseq, $rightaa-1-$offaa, 1) . ($rightaa-1) . "del";
		}
	    } else {
		$pdot = "p." . $aaL . $leftaa . "_" . $aaR . $rightaa . "delins" . $indelaa;
	    }
	}
    }

    $subkey = "";
    foreach(keys %subhash) { if( $_ =~ /\+/ ) { delete $subhash{$_}; } }
    foreach( sort {$a<=>$b} keys %subhash ) {
	if( length($subkey)>0 ) { $subkey .= ";" }
	$subkey .= $_ . ":" . substr($rnaseq,$_-$offrna,1) . ">" .$subhash{$_};
    }

    if( length($subkey)>0 ) {
	my $altseq = proc_subkey_to_altseq($subkey);
	my $altaa = proc_to_aa($altseq);
	my $altkey = "";
	if( $altaa ne $aaseq ) {
	    for( my $i=0; $i<length($altaa); $i++ ) {
		if( substr($altaa,$i,1) ne substr($aaseq,$i,1) ) {
		    if( length($altkey)>0 ) { $altkey .= ";"; }
		    $altkey .= ($i+$offaa) . ":" . substr($aaseq,$i,1) . ">" . substr($altaa,$i,1);
		}
	    }
	}
	if( length($altkey)>0 ) { $pdot .= "|p." . $altkey; }
    }

    return($pdot);
}

sub proc_alignment_dets_hash {
    # argument #1: reference to hash containing alignment details by seq coordinate
    my %adet = %{$_[0]};

    # argument #2: extended sequence of candidate itd
    my $candidateseq = $_[1];

    # merge consecutive hash entries when approrpriate:
    #  -- start positions in aligned seq & ref seq of hash entry equal to end positions of prev hash entry plus the
    #     same offset C allowed to be between 1 and 3
    my %adetNew = (); my %hits  = (); my $icurr;
    my @starts = sort {$a<=>$b} keys %adet;
    for( my $i=1; $i<scalar(@starts); $i++ ) {
	if( ($adet{$starts[$i]}{s0}-$adet{$starts[$i-1]}{s1}) > 0 &&
	    ($adet{$starts[$i]}{s0}-$adet{$starts[$i-1]}{s1}) <= 4 &&
	    ($adet{$starts[$i]}{s0}-$adet{$starts[$i-1]}{s1}) == ($adet{$starts[$i]}{r0}-$adet{$starts[$i-1]}{r1}) ) {
	    $hits{$i}=1;
	}
    }    
    for( my $i=0; $i<scalar(@starts); $i++ ) {
	if( !exists($hits{$i}) ) {
	    # make a new copy to avoid overwriting original hash passed by reference
	    $icurr = $i;
	    $adetNew{$starts[$icurr]}{s0} = $adet{$starts[$icurr]}{s0}; 
	    $adetNew{$starts[$icurr]}{s1} = $adet{$starts[$icurr]}{s1}; 
	    $adetNew{$starts[$icurr]}{r0} = $adet{$starts[$icurr]}{r0}; 
	    $adetNew{$starts[$icurr]}{r1} = $adet{$starts[$icurr]}{r1}; 
	    $adetNew{$starts[$icurr]}{mm} = $adet{$starts[$icurr]}{mm};
	} else {
	    $adetNew{$starts[$icurr]}{s1} = $adet{$starts[$i]}{s1};
	    $adetNew{$starts[$icurr]}{r1} = $adet{$starts[$i]}{r1};
	    my $offset = $adet{$starts[$i]}{s0} - $adet{$starts[$i-1]}{s1} - 1;
	    if( $offset > 0 ) {
		for( my $j=1; $j<=$offset; $j++ ) {
		    if( substr( $candidateseq, $adet{$starts[$i-1]}{s1}+$j-1, 1 ) ne 
			substr( $refseq, $adet{$starts[$i-1]}{r1}+$j-1, 1 ) ) {
			if( length($adetNew{$starts[$icurr]}{mm}) > 0 ) {
			    $adetNew{$starts[$icurr]}{mm} .= "+";
			}
			$adetNew{$starts[$icurr]}{mm} .= ($adet{$starts[$i-1]}{r1}+$j) . ":" . 
			    substr( $refseq, $adet{$starts[$i-1]}{r1}+$j-1, 1 ) . ">" .
			    substr( $candidateseq, $adet{$starts[$i-1]}{s1}+$j-1, 1 );
		    }
		}
	    }
	    if( length($adet{$starts[$i]}{mm}) > 0 ) {
		if( length($adetNew{$starts[$icurr]}{mm}) > 0 ) {
		    $adetNew{$starts[$icurr]}{mm} .= "+";
		}
		$adetNew{$starts[$icurr]}{mm} .= $adet{$starts[$i]}{mm};
	    }
	}
    }

    # processed alignment details will be returned in new hash %cdet (candidate itd details)
    my %cdet = ( 
	keyfull => "",
	keyrdot => "",
	keycdot => "", 
	keypdot => "", 
	keycdotalt => "",
	keypdotalt => "",
	dupseq => "",
	insseq => "",
	ie => "", # introns and exons overlapping itd
	duplen => 0,
	inslen => 0,
	totlen => length($candidateseq),
	itdsize => length($candidateseq)-1500,
	dup_r0 => -1,
	dup_r1 => -1,
	mutbk_s0 => -1,
	mutbk_s1 => -1,
	mismatches => 0
	);

    my $diff0 = 0;
    my $mmtot = 0;
    my $rstartNew = 0;
    my @startsNew = sort {$a<=>$b} keys %adetNew;

    # process the full key
    my $fullkey = "ref" . $adetNew{$startsNew[0]}{r0} . "-" . $adetNew{$startsNew[0]}{r1};
    for( my $i=0; $i<scalar(@startsNew); $i++ ) {
	if( $i>0 ) {
	    $diff0 = $adetNew{$startsNew[$i]}{s0} - $adetNew{$startsNew[$i-1]}{s1} - 1;
	    $rstartNew = $adetNew{$startsNew[$i]}{r0} - min($diff0, 0);

	    if( $diff0 > 0 ) { 
		$fullkey .= "_ins" . substr( $candidateseq, $adetNew{$startsNew[$i-1]}{s1}, $diff0 );
	    }
	    $fullkey .= "_ref" . $rstartNew . "-" . $adetNew{$startsNew[$i]}{r1};	    
	}

	if( length($adetNew{$startsNew[$i]}{mm}) > 0 ) {
	    $fullkey .= "(" . $adetNew{$startsNew[$i]}{mm} . ")"; 
	    my @tmparr = split( /\+/, $adetNew{$startsNew[$i]}{mm} );
	    $mmtot += scalar( @tmparr );
	}
	# Need logic for when $diff0 < 0 and mismatch is in this overlap region
    }
    $cdet{keyfull} = $fullkey;
    $cdet{mismatches} = $mmtot;

    # create rdot, cdot, and pdot keys
    my $leftr0 = ""; my $leftr1 = ""; my $leftmm = "";
    my $rightr0 = ""; my $rightr1 = ""; my $rightmm = "";
    my @cells; my @subcells; my @subsubcells;
    my $mmkey;

    @cells = split(/\_/, $fullkey);
    @subcells = split(/\(/, $cells[0]);
    @subsubcells = split(/\-/, $subcells[0]);
    if( substr($subsubcells[0],0,3) eq "ref" ) { $leftr0 = substr($subsubcells[0],3); $leftr1 = $subsubcells[1]; }
    if( scalar(@subcells) > 1 ) { $leftmm = substr($subcells[1],0,length($subcells[1])-1); }
    @subcells = split(/\(/, $cells[scalar(@cells)-1]);
    @subsubcells = split(/\-/, $subcells[0]);
    if( substr($subsubcells[0],0,3) eq "ref" ) { $rightr0 = substr($subsubcells[0],3); $rightr1 = $subsubcells[1]; }
    if( scalar(@subcells) > 1 ) { $rightmm = substr($subcells[1],0,length($subcells[1])-1); }

    my @cdotL1 = proc_to_cdot_coords( $leftr1 );
    my @cdotR0 = proc_to_cdot_coords( $rightr0 );
    $cdet{dup_r0} = $rightr0;
    $cdet{dup_r1} = $leftr1;
    $cdet{mutbk_s0} = $leftr1;

    if( $rightr0 > $leftr1 ) {
	$cdet{dup_r0} = $leftr1;
	$cdet{dup_r1} = $rightr0;
	$cdet{mutbk_s0} = $leftr1;
	$cdet{mutbk_s1} = $cdet{totlen} - (1500-$rightr0);

	@subcells = (); my $sub0 = -1; my $sub1 = -1; my $mutflag = 0;

	# should improve this to make more general
	if( scalar(@cells) > 3 && substr($cells[1],0,3) eq "ins" && substr($cells[2],0,3) eq "ref" ) {
	    @subcells = split(/\-/, substr($cells[2],3) ); 
	} elsif( scalar(@cells) > 2 && substr($cells[1],0,3) eq "ref" ) { 
	    @subcells = split(/\-/, substr($cells[1],3) );
	}
	if( scalar(@subcells) > 0 ) { 
	    $sub0 = $subcells[0];
	    @subsubcells = split(/\(/, $subcells[1]);
	    $sub1 = $subsubcells[0];

	    if( $sub0 < $leftr1 && $leftr1-$sub0+1 <= $cdet{totlen}-1500 ) { 
		$cdet{dup_r0} = $sub0;
		$cdet{dup_r1} = $leftr1;
		$cdet{duplen} = $cdet{dup_r1}-$cdet{dup_r0}+1;
		$cdet{dupseq} = substr( $refseq, $cdet{dup_r0}-1, $cdet{duplen} );
		if( substr($cells[1],0,3) eq "ins" ) {
		    $cdet{insseq} = substr($cells[1],3);
		    $cdet{inslen} = length($cdet{insseq});
		}
		$cdet{mutbk_s1} = $cdet{mutbk_s0} + $cdet{inslen} + 1;
		$mutflag = 1;
	    }
	}

	@subcells = (); $sub0 = -1; $sub1 = -1;	
	if( scalar(@cells) > 3 && substr($cells[scalar(@cells)-2],0,3) eq "ins" && 
	    substr($cells[scalar(@cells)-3],0,3) eq "ref" ) {
	    @subcells = split(/\-/, substr($cells[scalar(@cells)-3],3) ); 
	} elsif( scalar(@cells) > 2 && substr($cells[scalar(@cells)-2],0,3) eq "ref" ) {
	    @subcells = split(/\-/, substr($cells[scalar(@cells)-2],3) );
	}
	if( scalar(@subcells) > 0 && !$mutflag ) { 
	    $sub0 = $subcells[0];
	    @subsubcells = split(/\(/, $subcells[1]);
	    $sub1 = $subsubcells[0];
	    if( $sub1 > $rightr0 && $sub1-$rightr0 <= $cdet{totlen}-1500 && substr($cells[scalar(@cells)-2],0,3) eq "ins" ) {
		$cdet{dup_r0} = $rightr0;
		$cdet{dup_r1} = $sub1;
		$cdet{duplen} = $cdet{dup_r1}-$cdet{dup_r0}+1;
		$cdet{dupseq} = substr( $refseq, $cdet{dup_r0}-1, $cdet{duplen} );
		if( substr($cells[scalar(@cells)-2],0,3) eq "ins" ) {
		    $cdet{insseq} = substr($cells[scalar(@cells)-2],3);
		    $cdet{inslen} = length($cdet{insseq});
		}
		$cdet{mutbk_s0} = $cdet{mutbk_s1} - $cdet{inslen} - 1;
	    }
	}

	if( $rightr0 == $leftr1+1 ) {
	    $cdet{keyrdot} .= "r." . $leftr1 . "_" . $rightr0;
	    if( scalar(@cells) > 2 ) { $cdet{keyrdot} .=  "ins"; }
	} else {
	    $cdet{keyrdot} .= "r." . ($leftr1+1) . "_" . ($rightr0-1) . "del";
	    if( scalar(@cells) > 2 ) { $cdet{keyrdot} .=  "ins"; }
	}
	for( my $i=1; $i<scalar(@cells)-1; $i++ ) {
	    if($i>1) { $cdet{keyrdot} .= "/"; }
	    if( substr($cells[$i],0,3) eq "ins" || substr($cells[$i],0,3) eq "ref" ) { 
		$cdet{keyrdot} .= substr($cells[$i],3); 
	    }
	}
	if( length($leftmm)>0 ) { $cdet{keyrdot} .= "|left(r." . $leftmm . ")"; }
	if( length($rightmm)>0 ) { $cdet{keyrdot} .= "|right(r." . $rightmm . ")";  }

	if( length($cdotL1[0])>0 && length($cdotR0[0])>0 ) {
	    if( $rightr0 == $leftr1+1 ) {
		$cdet{keycdot} .= "c." . $cdotL1[0] . "_" . $cdotR0[0];
		if( scalar(@cells) > 2 ) { $cdet{keycdot} .=  "ins"; }
	    } else {
		$cdet{keycdot} .= "c." . increment_cdot_coord($cdotL1[0]) . 
		    "_" . decrement_cdot_coord($cdotR0[0]) . "del";
		if( scalar(@cells) > 2 ) { $cdet{keycdot} .=  "ins"; }
	    }
	    $cdet{ie} = $cdotL1[1];
	    if( scalar(@cells) == 2 && $cdotR0[1] ne $cdotL1[1] ) { $cdet{ie} .= "-" . $cdotR0[1]; }

	    for( my $i=1; $i<scalar(@cells)-1; $i++ ) {
		if($i>1) { $cdet{keycdot} .= "/"; }
		if( substr($cells[$i],0,3) eq "ins" ) { $cdet{keycdot} .= substr($cells[$i],3); }
		if( substr($cells[$i],0,3) eq "ref" ) { 
		    @subcells = split(/\(/, substr($cells[$i],3) );
		    @subsubcells = split(/\-/, $subcells[0]);
		    my @cdot0 = proc_to_cdot_coords( $subsubcells[0] );
		    my @cdot1 = proc_to_cdot_coords( $subsubcells[1] );
		    if( length($cdot0[0])>0 && length($cdot1[0])>0 ) {
			$cdet{keycdot} .= $cdot0[0] . "-" . $cdot1[0];
			if( $cdot0[1] eq $cdet{ie} || $cdot1[1] eq $cdet{ie} ) {
			    $cdet{ie} = $cdot0[1];
			} else {
			    $cdet{ie} .= "_" . $cdot0[1];
			}
			if( $cdot0[1] ne $cdot1[1] ) { $cdet{ie} .= "-" . $cdot1[1]; }
			if( scalar(@subcells) > 1 ) {
			    $mmkey = "";
			    @subsubcells = split(/\+/, substr($subcells[1],0,length($subcells[1])-1));
			    for( my $j=0; $j<scalar(@subsubcells); $j++ ) {
				my @subsubsubcells = split( /:/, $subsubcells[$j] );
				my @cdot2 = proc_to_cdot_coords( $subsubsubcells[0] );
				if( length($cdot2[0])>0 ) {
				    if( length($mmkey)>0 ) { $mmkey .= ";"; }
				    $mmkey .= $cdot2[0] . ":" . $subsubsubcells[1];
				}
			    }
			    if( length($mmkey)>0 ) { $cdet{keycdot} .= "(c." . $mmkey .  ")"; }
			} 
		    } else {
			$cdet{keycdot} .= "region(needs_manual_review)"; 
		    }
		}
	    }

	    $mmkey = "";
	    if( length($leftmm)>0 ) {
		@subcells = split( /\+/, $leftmm );
		for( my $i=0; $i<scalar(@subcells); $i++ ) {
		    @subsubcells = split( /:/, $subcells[$i] );
		    my @cdotmm = proc_to_cdot_coords($subsubcells[0]);
		    if( length($cdotmm[0]) > 0 ) {
			if(length($mmkey)>0) { $mmkey .= ";"; }
			$mmkey .= $cdotmm[0] . ":" . $subsubcells[1];
		    }
		}
	    }	
	    if( length($rightmm)>0 ) {
		@subcells = split( /\+/, $rightmm );
		for( my $i=0; $i<scalar(@subcells); $i++ ) {
		    @subsubcells = split( /:/, $subcells[$i] );
		    my @cdotmm = proc_to_cdot_coords($subsubcells[0]);
		    if( length($cdotmm[0]) > 0 ) {
			if(length($mmkey)>0) { $mmkey .= ";"; }
			$mmkey .= $cdotmm[0] . ":" . $subsubcells[1];
		    }
		}
	    }	
	    if( length($mmkey)>0 ) { $cdet{keycdot} .= "|ref_c." . $mmkey; }
	}
    } elsif( scalar(@cells) == 2 || ( scalar(@cells) == 3 && substr($cells[1],0,3) eq "ins" ) ) {
	my $rightmm0 = ""; my $rightmm1 = "";
	@subcells = split( /\+/, $rightmm );
	for( my $i=0;$i<scalar(@subcells); $i++ ) {
	    @subsubcells = split( /:/, $subcells[$i] );
	    if( $subsubcells[0] <= $leftr1 ) {
		if( length($rightmm0)>0 ) { $rightmm0 .= "+"; }
		$rightmm0 .= $subcells[$i];
	    } else {
		if( length($rightmm1)>0 ) { $rightmm1 .= "+"; }
		$rightmm1 .= $subcells[$i];
	    }
	}
	my $leftmm0 = ""; my $leftmm1 = "";
	@subcells = split( /\+/, $leftmm );
	for( my $i=0;$i<scalar(@subcells); $i++ ) {
	    @subsubcells = split( /:/, $subcells[$i] );
	    if( $subsubcells[0] >= $rightr0 ) {
		if( length($leftmm1)>0 ) { $leftmm1 .= "+"; }
		$leftmm1 .= $subcells[$i];
	    } else {
		if( length($leftmm0)>0 ) { $leftmm0 .= "+"; }
		$leftmm0 .= $subcells[$i];
	    }
	}

	$cdet{keyrdot} .= "r." . $leftr1 . "_" . ($leftr1+1);
	if( scalar(@cells) == 3 ) { $cdet{keyrdot} .= $cells[1] . "/"; }
	$cdet{keyrdot} .= "itd" . $rightr0 . "-" . $leftr1;
	if( length($rightmm0)>0 ) { $cdet{keyrdot} .= "(" . $rightmm0 . ")"; }
	if( length($leftmm)>0 ) { $cdet{keyrdot} .= "|left(r." . $leftmm . ")"; }
	if( length($rightmm1)>0 ) { $cdet{keyrdot} .= "|right(r." . $rightmm1 . ")"; }
	
	if( scalar(@cells) == 3 ) {
	    $cdet{insseq} = substr($cells[1],3);
	    $cdet{inslen} = length($cdet{insseq});
	    $cdet{mutbk_s1} = $leftr1 + $cdet{inslen} + 1;
	} else {
	    $cdet{mutbk_s1} = $leftr1 + 1;
	}
	$cdet{duplen} = $cdet{dup_r1}-$cdet{dup_r0}+1;
	$cdet{dupseq} = substr( $refseq, $cdet{dup_r0}-1, $cdet{duplen} );
	
	if( length($cdotL1[0])>0 && length($cdotR0[0])>0 ) {
	    $cdet{ie} = $cdotR0[1];
	    if( $cdotL1[1] ne $cdotR0[1] ) { $cdet{ie} .= "-" . $cdotL1[1]; }
	    if( scalar(@cells)==2 ) {
		$cdet{keycdot} = "c." . $cdotR0[0] . "_" . $cdotL1[0] . "dup";
	    } else {
		$cdet{keycdot} = "c." . $cdotL1[0] . "_" . increment_cdot_coord($cdotL1[0]) . $cells[1] . 
		    "/" . $cdotR0[0] . "-" . $cdotL1[0];  
		$cdet{keycdotalt} = "c." . decrement_cdot_coord($cdotR0[0]) . "_" . $cdotR0[0] .
		    "ins" . $cdotR0[0] . "-" . $cdotL1[0];
	    }
	    if( $cdotR0[1] eq "e14" && $cdotL1[1] eq "e15" ) { $cdet{keycdot} .= "[+i14]"; }
	    
	    if( length($leftmm1)>0 ) {
		$mmkey = "";
		@subcells = split( /\+/, $leftmm1 );
		for( my $i=0; $i<scalar(@subcells); $i++ ) {
		    @subsubcells = split( /:/, $subcells[$i] );
		    my @cdotmm = proc_to_cdot_coords($subsubcells[0]);
		    if( length($cdotmm[0]) > 0 ) {
			if(length($mmkey)>0) { $mmkey .= ";"; }
			$mmkey .= $cdotmm[0] . ":" . $subsubcells[1];
		    }
		}
		if( length($mmkey)>0 ) { 
		    if( scalar(@cells)==2 ) {
			$cdet{keycdot} .= "(c1." . $mmkey . ")"; 
		    } else {
			$cdet{keycdotalt} .= "(c." . $mmkey . ")"; 
		    }
		}
	    }
	    
	    if( scalar(@cells) > 2 ) { $cdet{keycdotalt} .= "/" . substr($cells[1],3); }  
	    
	    if( length($rightmm0)>0 ) {
		$mmkey = "";
		@subcells = split( /\+/, $rightmm0 );
		for( my $i=0; $i<scalar(@subcells); $i++ ) {
		    @subsubcells = split( /:/, $subcells[$i] );
		    my @cdotmm = proc_to_cdot_coords($subsubcells[0]);
		    if( length($cdotmm[0]) > 0 ) {
			if(length($mmkey)>0) { $mmkey .= ";"; }
			$mmkey .= $cdotmm[0] . ":" . $subsubcells[1];
		    }
		}
		if( length($mmkey)>0 ) {
		    if( scalar(@cells)==2 ) {
			$cdet{keycdot} .= "(c2." . $mmkey . ")"; 
		    } else {
			$cdet{keycdot} .= "(c." . $mmkey . ")"; 
		    }
		}
	    }

	    my $mmproc = ""; $mmkey = "";
	    if( scalar(@cells)==2 ) { $mmproc = $leftmm0; } else { $mmproc = $leftmm; }
	    if( length($rightmm1)>0 ) {
		if( length($mmproc)>0 ) { $mmproc .= "+"; }
		$mmproc .= $rightmm1;
	    }
	    if( length($mmproc) > 0 ) {
		@subcells = split( /\+/, $mmproc );
		for( my $i=0; $i<scalar(@subcells); $i++ ) {
		    @subsubcells = split( /:/, $subcells[$i] );
		    my @cdotmm = proc_to_cdot_coords($subsubcells[0]);
		    if( length($cdotmm[0]) > 0 ) {
			if(length($mmkey)>0) { $mmkey .= ";"; }
			$mmkey .= $cdotmm[0] . ":" . $subsubcells[1];
		    }
		}
	    }
	    if( length($mmkey)>0 ) { $cdet{keycdot} .= "|ref_c." . $mmkey; }	

	    if( scalar(@cells)==3 ) { 
		$mmkey = "";
		$mmproc = $leftmm0;
		if( length($rightmm)>0 ) {
		    if( length($mmproc)>0 ) { $mmproc .= "+"; }
		    $mmproc .= $rightmm;
		}
		if( length($mmproc) > 0 ) {
		    @subcells = split( /\+/, $mmproc );
		    for( my $i=0; $i<scalar(@subcells); $i++ ) {
			@subsubcells = split( /:/, $subcells[$i] );
			my @cdotmm = proc_to_cdot_coords($subsubcells[0]);
			if( length($cdotmm[0]) > 0 ) {
			    if(length($mmkey)>0) { $mmkey .= ";"; }
			    $mmkey .= $cdotmm[0] . ":" . $subsubcells[1];
			}
		    }
		}
		if( length($mmkey)>0 ) { $cdet{keycdotalt} .= "|ref_c." . $mmkey; }	
	    }
	}
    } else {
	my $insstr = "";
	my $inslen = 0;
	for( my $i=1; $i<scalar(@cells)-1; $i++ ) {
	    if( length($insstr)>0 ) { $insstr .= "/"; }
	    if( substr($cells[$i],0,3) eq "ins" ) {
		$insstr .= substr($cells[$i],3);
		$inslen += length(substr($cells[$i],3));
	    } elsif( substr($cells[$i],0,3) eq "ref" ) {
		@subcells = split(/\(/, substr($cells[$i],3) );
		@subsubcells = split(/\-/, $subcells[0]);
		$inslen += $subsubcells[1]-$subsubcells[0]+1;
		my @cdot0 = proc_to_cdot_coords( $subsubcells[0] );
		my @cdot1 = proc_to_cdot_coords( $subsubcells[1] );
		if( length($cdot0[0])>0 && length($cdot1[0])>0 ) {
		    $insstr .= $cdot0[0] . "-" . $cdot1[0];
		    if( $cdot0[1] eq $cdet{ie} || $cdot1[1] eq $cdet{ie} ) {
			$cdet{ie} = $cdot0[1];
		    } else {
			$cdet{ie} .= "_" . $cdot0[1];
		    }
		    if( $cdot0[1] ne $cdot1[1] ) { $cdet{ie} .= "-" . $cdot1[1]; }
		    if( scalar(@subcells) > 1 ) {
			$mmkey = "";
			@subsubcells = split(/\+/, substr($subcells[1],0,length($subcells[1])-1));
			for( my $j=0; $j<scalar(@subsubcells); $j++ ) {
			    my @subsubsubcells = split( /:/, $subsubcells[$j] );
			    my @cdot2 = proc_to_cdot_coords( $subsubsubcells[0] );
			    if( length($cdot2[0])>0 ) {
				if( length($mmkey)>0 ) { $mmkey .= ";"; }
				$mmkey .= $cdot2[0] . ":" . $subsubsubcells[1];
			    }
			}
			if( length($mmkey)>0 ) { $insstr .= "(c." . $mmkey .  ")"; }
		    } 
		} else {
		    $insstr .= "region(needs_manual_review)"; 
		}
	    }

	    $cdet{insseq} = "needs_manual_review";
	    $cdet{inslen} = $inslen;
	    $cdet{mutbk_s1} = $cdet{mutbk_s0} + $inslen + 1;
	    $cdet{keycdot} = "c." . $cdotL1[0] . "_" . increment_cdot_coord($cdotL1[0]) . "ins" . 
		$insstr . "/" . $cdotR0[0] . "-" . $cdotL1[0];
	    $cdet{keycdotalt} = "c." . decrement_cdot_coord($cdotR0[0]) . "_" . $cdotR0[0] . "ins" .
		$cdotR0[0] . "-" . $cdotL1[0];

	    my $rightmm0 = ""; my $rightmm1 = "";
	    @subcells = split( /\+/, $rightmm );
	    for( my $i=0;$i<scalar(@subcells); $i++ ) {
		@subsubcells = split( /:/, $subcells[$i] );
		if( $subsubcells[0] <= $leftr1 ) {
		    if( length($rightmm0)>0 ) { $rightmm0 .= "+"; }
		    $rightmm0 .= $subcells[$i];
		} else {
		    if( length($rightmm1)>0 ) { $rightmm1 .= "+"; }
		    $rightmm1 .= $subcells[$i];
		}
	    }
	    my $leftmm0 = ""; my $leftmm1 = "";
	    @subcells = split( /\+/, $leftmm );
	    for( my $i=0;$i<scalar(@subcells); $i++ ) {
		@subsubcells = split( /:/, $subcells[$i] );
		if( $subsubcells[0] >= $rightr0 ) {
		    if( length($leftmm1)>0 ) { $leftmm1 .= "+"; }
		    $leftmm1 .= $subcells[$i];
		} else {
		    if( length($leftmm0)>0 ) { $leftmm0 .= "+"; }
		    $leftmm0 .= $subcells[$i];
		}
	    }

	    if( length($leftmm1)>0 ) {
		$mmkey = "";
		@subcells = split( /\+/, $leftmm1 );
		for( my $i=0; $i<scalar(@subcells); $i++ ) {
		    @subsubcells = split( /:/, $subcells[$i] );
		    my @cdotmm = proc_to_cdot_coords($subsubcells[0]);
		    if( length($cdotmm[0]) > 0 ) {
			if(length($mmkey)>0) { $mmkey .= ";"; }
			$mmkey .= $cdotmm[0] . ":" . $subsubcells[1];
		    }
		}
		if( length($mmkey)>0 ) { 
		    $cdet{keycdotalt} .= "(c." . $mmkey . ")"; 
		}
	    }
	    $cdet{keycdotalt} .= "/" . $insstr; 
	    
	    if( length($rightmm0)>0 ) {
		$mmkey = "";
		@subcells = split( /\+/, $rightmm0 );
		for( my $i=0; $i<scalar(@subcells); $i++ ) {
		    @subsubcells = split( /:/, $subcells[$i] );
		    my @cdotmm = proc_to_cdot_coords($subsubcells[0]);
		    if( length($cdotmm[0]) > 0 ) {
			if(length($mmkey)>0) { $mmkey .= ";"; }
			$mmkey .= $cdotmm[0] . ":" . $subsubcells[1];
		    }
		}
		if( length($mmkey)>0 ) {
		    $cdet{keycdot} .= "(c." . $mmkey . ")"; 
		}
	    }

	    my $mmproc = $leftmm; $mmkey = "";
	    if( length($rightmm1)>0 ) {
		if( length($mmproc)>0 ) { $mmproc .= "+"; }
		$mmproc .= $rightmm1;
	    }
	    if( length($mmproc) > 0 ) {
		@subcells = split( /\+/, $mmproc );
		for( my $i=0; $i<scalar(@subcells); $i++ ) {
		    @subsubcells = split( /:/, $subcells[$i] );
		    my @cdotmm = proc_to_cdot_coords($subsubcells[0]);
		    if( length($cdotmm[0]) > 0 ) {
			if(length($mmkey)>0) { $mmkey .= ";"; }
			$mmkey .= $cdotmm[0] . ":" . $subsubcells[1];
		    }
		}
	    }
	    if( length($mmkey)>0 ) { $cdet{keycdot} .= "|ref_c." . $mmkey; }	

	    $mmproc = $leftmm0; $mmkey = "";
	    if( length($rightmm)>0 ) {
		if( length($mmproc)>0 ) { $mmproc .= "+"; }
		$mmproc .= $rightmm;
	    }
	    if( length($mmproc) > 0 ) {
		@subcells = split( /\+/, $mmproc );
		for( my $i=0; $i<scalar(@subcells); $i++ ) {
		    @subsubcells = split( /:/, $subcells[$i] );
		    my @cdotmm = proc_to_cdot_coords($subsubcells[0]);
		    if( length($cdotmm[0]) > 0 ) {
			if(length($mmkey)>0) { $mmkey .= ";"; }
			$mmkey .= $cdotmm[0] . ":" . $subsubcells[1];
		    }
		}
	    }
	    if( length($mmkey)>0 ) { $cdet{keycdotalt} .= "|ref_c." . $mmkey; }	
	}
    } 

    if( $cdet{keycdot} ne "" ) { $cdet{keypdot} = proc_cdot_to_pdot( $cdet{keycdot} ); }
    if( $cdet{keycdotalt} ne "" ) { $cdet{keypdotalt} = proc_cdot_to_pdot( $cdet{keycdotalt} ); }

    return( %cdet );
}

sub proc_vcf_from_keyfull {
    # argument #1: keyfull
    my $keyf = $_[0];
    #argument #2: keycdot
    my $keyc = $_[1];

    my @regions = split(/\_/, $keyf);
    my $ref0 = $regions[0]; my $ref1 = $regions[scalar(@regions)-1];
    if( scalar(@regions) < 2 || $ref0 !~ /^ref/ || $ref1 !~ /^ref/ ) { return(""); }
    
    $ref0 =~ s/ref//; $ref1 =~ s/ref//;
    my $mm0 = ""; my $mm1 = ""; 
    my @cells = split( /\(/, $ref0 );
    my @coords0 = split( /\-/, $cells[0] );
    if( scalar(@cells) > 1 ) { $mm0 = $cells[1]; $mm0 =~ s/\)//; }

    @cells = split( /\(/, $ref1 );
    my @coords1 = split( /\-/, $cells[0] );
    if( scalar(@cells) > 1 ) { $mm1 = $cells[1]; $mm1 =~ s/\)//; }

    my $seq0 = $refseq; my $seq1 = $refseq; my @mmlocs0 = (); my @mmlocs1 = (); 
    if( length($mm0) > 0 ) {
	@cells = split(/\+/, $mm0);
	foreach my $mm ( @cells ) {
	    my @subcells = split( /:/, $mm );
	    push @mmlocs0, $subcells[0];
	    my @subsubcells = split (/>/, $subcells[1] );
	    substr( $seq0, $subcells[0]-1, 1 ) = $subsubcells[1];
	}
    }

    if( length($mm1) > 0 ) {
	@cells = split(/\+/, $mm1);
	foreach my $mm ( @cells ) {
	    my @subcells = split( /\:/, $mm );
	    push @mmlocs1, $subcells[0];
	    my @subsubcells = split (/>/, $subcells[1] );
	    substr( $seq1, $subcells[0]-1, 1 ) = $subsubcells[1];
	}
    }

    my $inspos = $coords0[1];
    if( length($mm0) > 0 ) { $inspos = min(@mmlocs0)-1; }
    my $dupstart = $coords1[0];
    if( length($mm1) > 0 ) { $dupstart = max(@mmlocs1)+1; }

    my $ins = "";
    if( scalar(@regions) > 2 ) {
	for( my $i=1; $i < scalar(@regions)-1; $i++ ) {
	    my $tmpseq = $regions[$i];
	    if( $tmpseq =~ /^ins/ ) {
		$tmpseq =~ s/ins//;
	    } elsif( $tmpseq =~ /^ref/ ) {
		$tmpseq =~ s/ref//;		
		@cells = split( /\(/, $tmpseq );
		my @coords = split( /\-/, $cells[0] ); my $seq2 = $refseq;
		if( scalar(@cells) > 1 ) {
		    my $mm2 = $cells[1]; $mm2 =~ s/\)//;
		    @cells = split(/\+/, $mm2);
		    foreach my $mm ( @cells ) {
			my @subcells = split( /\:/, $mm );
			my @subsubcells = split (/>/, $subcells[1] );
			substr( $seq2, $subcells[0]-1, 1 ) = $subsubcells[1];
		    }
		}
		$tmpseq = substr( $seq2, $coords[0]-1, $coords[1]-$coords[0]+1 );
	    } else {
		return( "" );
	    }
	    $ins .= $tmpseq;
	}
    }
    
    my $alt;
    my $ref;
    if( $inspos+1 >= $dupstart ) {
	$ref = substr( $refseq, $inspos, 1 );
	$alt = substr($seq0, $inspos, $coords0[1]-$inspos) . $ins . 
	    substr($seq1, $coords1[0]-1, $inspos-$coords1[0]+1) . $ref;
    } else {
	$ref = substr( $refseq, $inspos, $dupstart-$inspos-1 );
	$alt = substr($seq0, $inspos, $coords0[1]-$inspos) . $ins . 
	    substr($seq1, $coords1[0]-1, $dupstart-$coords1[0]);
    }
    $alt = reverse $alt; $alt =~ tr/ATGCatgc/TACGtacg/; 
    $ref = reverse $ref; $ref =~ tr/ATGCatgc/TACGtacg/;

    return("13\t" . ($refend-$inspos) . "\t" . $keyc . "\t" . $ref . "\t" . $alt . "\t.\t." );
}

####################################################################
####################################################################
# BEGIN MAIN SCRIPT
####################################################################
####################################################################

my $fbase; my $outname; my %umitags = (); my %umisWt = (); my %umisMut = ();
if( length($inbam) > 0 ) {
  @line_cells = split(/\//, $inbam);
  $outname = $line_cells[scalar(@line_cells)-1]; $outname =~ s/\.bam$//; $outname .= "_FLT3";
  $fbase = $outpath . "/" . $outname;
  $fastq1 = $fbase . ".R1.fastq";
  $fastq2 = $fbase . ".R2.fastq";

  if( $typeb eq "targeted" || $typeb eq "loose" ) {
    # get chr13 name from bam file
    if( system( sprintf( "%s view -H %s | grep 'SN:chr13\\|SN:13' | grep -v 'SN:chr13_\\|SN:13_' > %s_header_chr13.txt", $samcmd, $inbam, $fbase ) ) ) {
      die "Failed to extract header from original input bam. Exiting...\n";
    }
    open(FI, $fbase . "_header_chr13.txt" ) or die $!;
    while(<FI>) { if($_ =~ /SN:(\w+)/) { $targetregion = $1 . ":" . $refstart .  "-" . $refend; } }
    close FI;
    if( system( sprintf( "rm %s_header_chr13.txt", $fbase )) ) {
      die "Failed to remove header chr13 file. Exiting...\n";
    }
  }

  if( $typeb eq "targeted" ) {
    # Extract alignments to the FLT3 target locus from the input bamfile
    # This is assumed to include all mutant reads (which should be true if a local aligner was used)
    if( system( sprintf( "%s view -b -F 0x900 %s %s > %s.bam", $samcmd, $inbam, $targetregion, $fbase ) ) ||
        system( sprintf( "%s index %s.bam", $samcmd, $fbase ) ) ||
        system( sprintf( "%s sort -n %s.bam -o %s.sorted.bam", $samcmd, $fbase, $fbase ) ) ) {
      die "Failed to extract existing FLT3 targeted locus alignments from original input bam. Exiting...\n";
    }
    if( system( sprintf( "rm %s.bam", $fbase )) ||
 	system( sprintf( "rm %s.bam.bai", $fbase )) ) {
      die "Failed to remove presorted target bam files. Exiting...\n";
    }
  } elsif( $typeb eq "loose" ) {
    # Extract alignments to the FLT3 target locus from the input bamfile
    print "Extracting FLT3 target locus...\n";
    if( system( sprintf( "%s view -b -F 0x900 %s %s > %s.bam", $samcmd, $inbam, $targetregion, $fbase ) ) ||
        system( sprintf( "%s index %s.bam", $samcmd, $fbase ) ) ) {
      die "Failed to extract existing FLT3 targeted locus alignments from original input bam. Exiting...\n";
    }

    # Extract unmapped alignments from the input bamfile
    print "Extracting unmapped reads...\n";
    if( system( sprintf( "%s view -bf 4 -F 0x900 %s > %s_unmapped.bam", $samcmd, $inbam, $fbase ) ) ||
	system( sprintf( "%s index %s_unmapped.bam", $samcmd, $fbase ) ) ) {
      die "Failed to extract unmapped alignments from original input bam. Exiting...\n";
    }

    # Merge target and unmapped files
    print "Merging FLT3 target locus and unmapped bams...\n";
    if( system( sprintf( "%s merge %s_merged.bam %s.bam %s_unmapped.bam", $samcmd, $fbase, $fbase, $fbase ) ) ) {
	die "Failed to merge target and unmapped files. Exiting...\n";
    }

    # Index and sort merged file
    if( system( sprintf( "%s index %s_merged.bam", $samcmd, $fbase ) ) ||
	system( sprintf( "%s sort -n %s_merged.bam -o %s.sorted.bam", $samcmd, $fbase, $fbase ) ) ) {
      die "Failed to extract existing FLT3 targeted locus alignments from original input bam. Exiting...\n";
    }

    # Remove presorted files
    if( system( sprintf( "rm %s.bam", $fbase )) ||
 	system( sprintf( "rm %s.bam.bai", $fbase )) ||
	system( sprintf( "rm %s_unmapped.bam", $fbase )) ||
 	system( sprintf( "rm %s_unmapped.bam.bai", $fbase )) ||
	system( sprintf( "rm %s_merged.bam", $fbase )) ||
 	system( sprintf( "rm %s_merged.bam.bai", $fbase )) ) {
      die "Failed to remove presorted bam files. Exiting...\n";
    }
  } elsif( $typeb eq "all" ) {
    if( system( sprintf( "%s sort -n %s -o %s.sorted.bam", $samcmd, $inbam, $fbase ) ) ) {
      die "Failed to sort original input bam by readnames. Exiting...\n";
    }
  } else {
      print "ERROR: typeb must either equal \"targeted\", \"loose\", or \"all\"...";
      HelpMessage(1);
  }

  if( $umitag ne "" ) {
    if( system( sprintf( "%s view %s.sorted.bam > %s.sorted.sam", 
	$samcmd, $fbase, $fbase ) ) ) {
      die "Failed to create sam file for extracting umi tags";
    }

    open( FI, $fbase . ".sorted.sam" ) or die $!;
    while(<FI>) {
      chomp;
      @line_cells = split( /\t/, $_ );
      if( /${umitag}(\S+)/ && !exists( $umitags{$line_cells[0]} ) ) {
        $umitags{$line_cells[0]} = $umitag . $1;
      }
    }
    close FI;
    if( system( sprintf( "rm %s.sorted.sam", $fbase )) ) {
      die "Failed to remove sam file for extracting umi tags";
    }
  }

  if( system( sprintf( "%s -i %s.sorted.bam -fq %s -fq2 %s", 
		 $bedcmd, $fbase, $fastq1, $fastq2 ) ) ) {
    die "Failed to convert bamfile to fastq file. Exiting...\n";
  }

  if( system( sprintf( "rm %s.sorted.bam", $fbase )) ) {
    die "Failed to remove intermediate bam files. Exiting...\n";
  }
} else {
  @line_cells = split(/\//, $fastq1 );
  $outname = $line_cells[scalar(@line_cells)-1];
  $outname =~ s/\.gz$//; $outname =~ s/\.fastq$//;
  $outname =~ s/\.R1$//; $outname =~ s/\_R1$//;
  $outname .= "_FLT3";
  $fbase = $outpath . "/" . $outname;
}

my $shortkey = $outname; $shortkey =~ s/\.consensus//; $shortkey =~ s/\.raw//; 
$shortkey =~ s/\.filtered\_FLT3//; $shortkey =~ s/\.mapped\_FLT3//; $shortkey =~ s/\.aligned_FLT3//;

# clean up old files
if( glob($fbase . "*.html") ) { system( "rm " . $fbase . "*.html" ); }
if( glob($fbase . "*.vcf") ) { system( "rm " . $fbase . "*.vcf" ); }
if( glob($fbase . "*summary.txt") ) { system( "rm " . $fbase . "*summary.txt" ); }

if( $adapter ) {
  my $trimfq1 = $fastq1; $trimfq1 =~ s/\.R1\.fastq/\.trimmed\.R1\.fastq/;
  my $trimfq2 = $fastq2; $trimfq2 =~ s/\.R2\.fastq/\.trimmed\.R2\.fastq/;
  if(
      system( sprintf( "%s -Xmx2g in1=%s in2=%s out1=%s out2=%s literal='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT' ktrim=r k=23 hdist=1 mink=11 hdist=1 ordered=t ignorebadquality", $trimcmd, $fastq1, $fastq2, $trimfq1, $trimfq2 ) ) ||
      system( sprintf( "mv %s %s", $trimfq1, $fastq1 ) ) ||
      system( sprintf( "mv %s %s", $trimfq2, $fastq2 ) ) ) {
    die "Failed to trim adapters. Exiting...\n";
  }
}

# Local alignments to reference
if( $proc_index && system( sprintf( "%s index -p %s %s", $bwacmd, $refindex, $reffasta ) ) ) {
  die "Failed to create bwa index for reference.  Exiting...\n" ;
}

if( system( sprintf( "%s mem -k %s -M -O 6 -T %s %s %s %s | grep FLT3_dna_e14e15 > %s_%s.sam", 
	$bwacmd, $bwaSeedLength, $bwaThreshold, $refindex, $fastq1, $fastq2, $fbase, $alignerLocal ) ) ) {
  die "Failed to locally align initial reads to targeted FLT3 locus. Exiting...\n";
}

my $readsuffix; my $str; # strand
my %endsFwdR2 = (); my %endsRevR2 = ();
my %psam = (); my %ssam = (); my %pseq = (); my %pqa = (); my %pstr = (); my %spos = ();
my %pcig = (); my %scig = (); my %ploc = (); my %sloc = (); my %sstr = (); my %ppos = ();
my %resolvedWt = (); my %resolvedWtPaired = (); my %readsPassed = (); my %pbase = (); 
open( FI, $fbase . "_" . $alignerLocal . ".sam" ) or die $!;
while(<FI>){
    if( !/^@/ ) {
	@line_cells = split( /\t/, $_ );

	if( $line_cells[2] eq "*" ) { next; }

	if( ( $line_cells[1] & 64 ) == 64 ) {
	    $readsuffix = "_1";
	} elsif( ( $line_cells[1] & 128 ) == 128 ) {
	    $readsuffix = "_2";
	} else {
	    next;
	}
	$readname = $line_cells[0] . $readsuffix;

	$cigar = $line_cells[5];
	if( $cigar =~ /^(\d+)H(.*)H$/ ) { next; }
	# For Archer, assume GSP2s are on R2
	if( $ngstype ne "Archer" || $readsuffix eq "_2" ) { $readlengths{length( $line_cells[9] )} += 1; }
	my $rpos = $line_cells[3]-1; 
	while( $cigar =~ /(\d+)([SHDIMN])(.*)/ ) {
	    $cigar = $3;
	    if( $2 eq "M" || $2 eq "N" || $2 eq "D" ) { $rpos += $1; }
	}

	if( ( $line_cells[1] & 16 ) == 16 ) { $str = "-"; } else { $str = "+"; }

	if( ( $line_cells[1] & 256 ) == 256 ) {
	    $ssam{$readname} = $_;
	    $sloc{$readname} = $line_cells[3]-1;
	    $scig{$readname} = $line_cells[5];
	    $spos{$readname} = $rpos;
	    $sstr{$readname} = $str;
	} elsif( ( $line_cells[1] & 2048 ) != 2048 ) {
	    $pbase{$line_cells[0]} = 1;
	    $psam{$readname} = $_;
	    $ploc{$readname} = $line_cells[3]-1;
	    $pcig{$readname} = $line_cells[5];
	    $pseq{$readname} = $line_cells[9]; 
	    $pqa{$readname} = $line_cells[10]; 
	    $ppos{$readname} = $rpos;
	    $pstr{$readname} = $str;

	    if( $readsuffix eq "_2" ) {
		if( ( $line_cells[1] & 16 ) == 0 ) {
		  $endsFwdR2{$line_cells[3]} += 1;
		} else {
		  $endsRevR2{$rpos} += 1;
		}
	    }
	}
    }
}
close FI;

if( !$debug && system( "rm " . $fbase . "_" . $alignerLocal . ".sam" ) ) {
  die "Error removing local alignment file to FLT3 target locus. Exiting...";
}

my %extBySeq = (); my %extByLength = (); my %clipsfastq = (); my %lloc = (); my %rloc = ();
foreach my $read ( keys %psam ) {
    my $extSeq = ""; my $extLoc = "";
    if( exists( $ssam{$read} ) ) {
	if( $pcig{$read} =~ /^(\d+)(S)(.*)/ && $pcig{$read} =~ /(\d+)(M)$/ && $scig{$read} =~ /^(\d+)(M)(.*)/ ) {
	    $lloc{$read} = $sloc{$read};
	    $rloc{$read} = $ppos{$read};
	} elsif( $pcig{$read} =~ /(\d+)(S)$/ && $pcig{$read} =~ /^(\d+)(M)(.*)/ && $scig{$read} =~ /(\d+)(M)$/ ) {
	    $lloc{$read} = $ploc{$read};
	    $rloc{$read} = $spos{$read};
	} else {
	    next;
	}
    } else {
	$cigar = $pcig{$read}; 
	my $clipLen = 0;
	my $insLen = 0;
	my $delLen = 0;
	if( @a = $cigar =~ /(\d+)S/g ) { $clipLen = max(@a); }
	if( @a = $cigar =~ /(\d+)I/g ) { foreach my $x (@a) { $insLen += $x; } }
	if( @a = $cigar =~ /(\d+)D/g ) { foreach my $x (@a) { $delLen += $x; } }

	if( $insLen-$delLen >= 3 && $clipLen == 0 ) {  # net insertion of 3+ bp
	    $lloc{$read} = $ploc{$read};
	    $rloc{$read} = $ppos{$read};
	} elsif( $clipLen >= $minClipToAlign ) { # clips without seconday alignments
	    if( $pcig{$read} =~ /^(\d+)(S)/ && $1 == $clipLen && $pcig{$read} =~ /(\d+)(M)$/ ) { 
		$clipsfastq{$read} = "@" . $read . "\n" . substr( $pseq{$read}, 0, $clipLen ) . 
		    "\n+\n" . substr( $pqa{$read}, 0 , $clipLen ) . "\n";
	    } elsif( $pcig{$read} =~ /(\d+)(S)$/ && $1 == $clipLen && $pcig{$read} =~ /^(\d+)(M)/ ) { 
		$cigar = $pcig{$read}; my $spos = 0;
		while( $cigar =~ /(\d+)([SHDIMN])(.*)/ ) {
		    $cigar = $3;
		    if( $2 eq "M" || $2 eq "N" || $2 eq "I" || $2 eq "S" || $2 eq "H" ) { $spos += $1; }
		}
		$clipsfastq{$read} = "@" . $read . "\n" . substr( $pseq{$read}, $spos-$clipLen ) . 
		    "\n+\n" . substr( $pqa{$read}, $spos-$clipLen ) . "\n";
	    }
	} elsif( $insLen <= $maxInsAllowed && $delLen <= $maxDelAllowed && $clipLen <= $maxClipAllowed && 
		 ($psam{$read} =~ /NM:i:(\d+)/ && $1 <= $maxEditDist) ) {
	    $resolvedWt{$read} = 1;  # reads aligning "stringently" to reference 
	}
    }
}

foreach my $read (keys %resolvedWt) {
    $read =~ s/\_1$//g; $read =~ s/\_2$//g;
    if( exists($resolvedWt{$read."_1"}) && exists($resolvedWt{$read."_2"}) && $pstr{$read."_1"} ne $pstr{$read."_2"} ) {
	$resolvedWtPaired{$read} = 1;
	$readsPassed{$read} = 1;
    }
}

if( $umitag ne "" ) {
  open( FO, ">" . $fbase . "_wt.sam" ) or die $!;
  foreach( sort keys %resolvedWtPaired ) {
    my $sam = $psam{$_ . "_1"}; chomp( $sam );
    print FO $sam;
    if( exists( $umitags{$_} ) ) { print FO "\tMQ:i:60\t" . $umitags{$_}; }
    print FO  "\n";
    $sam = $psam{$_ . "_2"}; chomp( $sam );
    print FO $sam;
    if( exists( $umitags{$_}  ) ) { print FO  "\tMQ:i:60\t" . $umitags{$_}; }
    print FO "\n";
  }
  close FO;

  if( system( sprintf( "%s view -bT %s %s_wt.sam | " .
                       "%s -jar %s SortSam I=/dev/stdin SO=queryname O=/dev/stdout | " .
                       "%s -jar %s SetMateInformation -i /dev/stdin -o %s_wt.bam",
        $samcmd, $reffasta, $fbase, $javacmd, $picardjar, $javacmd, $fgbiojar, $fbase ) ) ) {
    die "Failed to convert wt sam to bam. Exiting...\n";
  }

  if( system( sprintf( "%s -jar %s GroupReadsByUmi -i %s_wt.bam -o %s_wt_umigroups.bam -s %s -f %s_wt_umisizes.txt", 
	$javacmd, $fgbiojar, $fbase, $fbase, $strat, $fbase ) ) ) {
    die "Failed to perform fgbio GroupReadsByUmi on wildtype reads...\n";
  }

  if( system( sprintf( "%s view %s_wt_umigroups.bam > %s_wt_umigroups.sam", 
	$samcmd, $fbase, $fbase ) ) ) {
    die "Failed to convert wt umigroups bam to sam...\n";
  }

  open( FI, $fbase . "_wt_umigroups.sam" ) or die $!;
  while(<FI>) {
    @line_cells = split(/\t/, $_);
    if( $_ =~ /MI:Z:(\d+)/ ) { $umisWt{$line_cells[0]} = $1; }
  }

  if( !$debug && (
    system( sprintf( "rm %s_wt_umisizes.txt", $fbase ) ) ||
    system( sprintf( "rm %s_wt_umigroups.sam", $fbase ) ) ||
    system( sprintf( "rm %s_wt_umigroups.bam", $fbase ) ) ) ) {
    die "Failed to remove wt umigroups bam...\n"; 
  }

  if( !$debug && system( sprintf( "rm %s_wt.*", $fbase ) ) ) { die "Failed to remove wt alignment files...\n"; }
}

my $entryFA = "";
open(FO1, ">" . $fbase . ".nowt.R1.fastq" ) or die $!;
open(FO2, ">" . $fbase . ".nowt.R2.fastq" ) or die $!;
foreach my $read (sort keys %pbase ) {
    if( exists($psam{$read . "_1"}) && exists($psam{$read . "_2"}) ) {
	@line_cells = split(/\t/, $psam{$read . "_1"});
	$seq = $line_cells[9]; $qa = $line_cells[10];
       	if( ($line_cells[1] & 16) == 16 ) { $seq = reverse $seq; $seq =~ tr/ATGCatgc/TACGtacg/; $qa = reverse $qa; }
	$entryFA = "@" . $read . "/1\n" . $seq . "\n" . "+" . "\n" . $qa . "\n";
        if( !exists($resolvedWtPaired{$read}) ) { print FO1 $entryFA; }
	@line_cells = split(/\t/, $psam{$read . "_2"});
	$seq = $line_cells[9]; $qa = $line_cells[10];
       	if( ($line_cells[1] & 16) == 16 ) { $seq = reverse $seq; $seq =~ tr/ATGCatgc/TACGtacg/; $qa = reverse $qa; }	
	$entryFA = "@" . $read . "/2\n" . $seq . "\n" . "+" . "\n" . $qa . "\n";
	if( !exists($resolvedWtPaired{$read}) ) { print FO2 $entryFA; }
    }
}
close FO1;
close FO2;

if( scalar( keys %clipsfastq ) > 0 ) {
  open (FO, ">" . $fbase . "_clips.fastq" ) or die $!;
  foreach( values %clipsfastq ) { print FO $_; }
  close FO;

  # Local alignment of clips to reference
  if( system( sprintf( "%s mem -k %s -M -O 6 -T %s %s %s_clips.fastq > %s_clips.%s.sam", 
	$bwacmd, $bwaClipsSeedLength, $bwaClipsThreshold, $refindex , $fbase, $fbase, $alignerClips ) ) ) {
    die "Error locally aligning clips to FLT3 target locus.  Exiting...";
  }

  # TO-DO: explore multiple hits case; these are in XA:Z fields
  my $clipsUnresolved = "";
  open( FI, $fbase . "_clips." . $alignerClips . ".sam" ) or die $!;
  while(<FI>) {
    if( !/^@/ ) {
      @s_cells = split(/\t/, $_ );

      # all clips taken relative to reference strand alignment
      if( $s_cells[5] eq "*" || ( $s_cells[1] & 16 ) == 16 ) { $clipsUnresolved .= $_; next; }
      if( ( $s_cells[1] & 256 ) == 256 ) { next; }

      # criteria to anchor a clip to the reference
      my $matchlen = 0;
      my $s_rpos = $s_cells[3] - 1; $cigar = $s_cells[5];
      while( $cigar =~ /(\d+)([SHDIMN])(.*)/ ) {
	$cigar = $3;
	if( $2 eq "M" || $2 eq "N" || $2 eq "D" ) { $s_rpos += $1; }
	if( $2 eq "M" ) { $matchlen += $1; }
      }
      if( $matchlen / length($s_cells[9]) < $clipMatchRatioThreshold ) { $clipsUnresolved .= $_; next; }
      my $read = $s_cells[0];
      if( $pcig{$read} =~ /^(\d+)(S)/ && $pcig{$read} =~ /(\d+)(M)$/ && $s_cells[5] =~ /^(\d+)(M)/ ) {
	$lloc{$read} = $s_cells[3]-1; $rloc{$read} = $ppos{$read};
      } elsif( $pcig{$read} =~ /(\d+)(S)$/ && $pcig{$read} =~ /^(\d+)(M)/ && $s_cells[5] =~ /(\d+)(M)$/ ) {
	$lloc{$read} = $ploc{$read}; $rloc{$read} = $s_rpos;
      } else {
	$clipsUnresolved .= $_;
      }
    }
  }

  open( FO, ">" . $fbase . "_clipsUnresolved." . $alignerClips . ".sam" ) or die $!;
  print FO $clipsUnresolved;
  close FO;

  if( $proc_manual_clipalign ) {
    open( FI, $fbase . "_clipsUnresolved." . $alignerClips . ".sam" ) or die $!;
    while(<FI>) {
      if( !/^@/ ) {
	@s_cells = split(/\t/, $_ );
	my $read = $s_cells[0];
	my $testseq = $s_cells[9];
	if( ($s_cells[1] & 16) == 16 ) { $testseq = reverse $testseq; $testseq =~ tr/ATGCatgc/TACGtacg/; }
	my $ind; my $offset;
	if( $pcig{$read} =~ /^(\d+)(S)/ && $pcig{$read} =~ /(\d+)(M)$/ ) {
	  OUTER:
	    while( length($testseq) >= $manualItdCheckLength ) {
	      $ind = index( $refseq, $testseq ); # TO-DO: logic for multiple hits?
	      if( $ind >= 0 ) { last OUTER; }
	      for( $i = 0; $i < length($testseq)-2; $i++ ) {
		foreach( ("A","T","C","G") ) {
		  my $mutseq = substr( $testseq, 0, $i ) . $_ . substr( $testseq, $i+1 );
		  $ind = index( $refseq, $mutseq );
		  if( $ind >= 0 ) { last OUTER; }
		}
	      }
	      $testseq = substr( $testseq, 0, length( $testseq )-1 );
	    }
	    if( length($testseq) < $manualItdCheckLength ) { next; }
	    $lloc{$read} = $ind; $rloc{$read} = $ppos{$read};
	    print $read ."\t" . $testseq . "\n";
	} elsif( $pcig{$read} =~ /(\d+)(S)$/ && $pcig{$read} =~ /^(\d+)(M)/ ) {
	  OUTER:
	    while( length($testseq) >= $manualItdCheckLength ) {
	      $ind = index( $refseq, $testseq ); # TO-DO: logic for multiple hits?
	      if( $ind >= 0 ) { last OUTER; }
	      for( $i = 2; $i < length($testseq); $i++ ) {
		foreach( ("A","T","C","G") ) {
		  my $mutseq = substr( $testseq, 0, $i ) . $_ . substr( $testseq, $i+1 );
		  $ind = index( $refseq, $mutseq );
		  if( $ind >= 0 ) { last OUTER; }
		}
	      }
	      $testseq = substr( $testseq, 1 );
	    }
	    if( length($testseq) < $manualItdCheckLength ) { next; }
	    $ind = $ind + length($testseq);
	    $lloc{$read} = $ploc{$read}; $rloc{$read} = $ind;
	}
      } 
    }
  }
  
  if( !$debug && (
      system( "rm " . $fbase . "_clips.fastq" ) ||
      system( "rm " . $fbase . "_clips." . $alignerClips . ".sam" ) ||
      system( "rm " . $fbase . "_clipsUnresolved." . $alignerClips . ".sam" ) ) ) {
    die "Error removing unresolved clips. Exiting...";
  }
}

foreach my $read ( keys %lloc ) {
  my $itdsize = length($pseq{$read}) - $rloc{$read} + $lloc{$read};
  if( $itdsize > 0 && $itdsize <= $maxITDsize ) {
    my $extSeq = substr($refseq,0,$lloc{$read}) . $pseq{$read} . substr($refseq,$rloc{$read});
    $extByLength{$itdsize} += 1;
    $extBySeq{$itdsize}{$extSeq} += 1;
  }
}

my %sumaHash = (); my $sumaclust = ""; my %extsizeHash = (); my %extpreHash = ();
$i = 0;
foreach $j ( sort {$extByLength{$b}<=>$extByLength{$a}} keys %extByLength ) {
    foreach ( sort{ $extBySeq{$j}{$b}<=>$extBySeq{$j}{$a} } keys %{$extBySeq{$j}} ) {
	if( $extBySeq{$j}{$_} >= $minToCluster ){
	    for( $k=0; $k<$extBySeq{$j}{$_}; $k++ ) {
		$i += 1;
		$key = $i . "_length" . $j . "_size=" . $extBySeq{$j}{$_};
		$sumaHash{$key} = $_;
		$extsizeHash{$key} = length($_);
		$extpreHash{$key} = $extBySeq{$j}{$_}; # num of initial reads generating the extended read
		$sumaclust .= ">" . $key . "\n" . $_ . "\n";
	    }
	}
    }
}

if( length( $sumaclust ) == 0 ) {
    print "NO ITD CANDIDATE CLUSTERS GENERATED. Exiting...\n"; 
    system( "touch " . $fbase . "_ITD_none.vcf" );
    if( !$debug && ( system( "rm " . $fbase . ".nowt.R*.fastq" ) || 
        ( $inbam ne "" && system( "rm " . $fbase . ".R*.fastq" ) ) ) ) {
      die "Error removing fastq files. Exiting...";
    }
    exit(0);
}

open (FO, ">" . $fbase . "_for_sumaclust.fasta" ) or die $!;
print FO $sumaclust;
close FO;

if( system( $clustercmd . " -d -r -t 5 " . $fbase . "_for_sumaclust.fasta -O " . $fbase . "_sumaclust.clusters" .  " -F " . $fbase . "_sumaclust.fasta" ) ) {
  print "Sumaclust command failed. Proceeding manually with all candidates.\n"; 
  open(FO, ">" . $fbase . "_sumaclust.clusters" ) or die $!;
  foreach my $atmp (sort keys %sumaHash) { print FO $atmp . "\n"; }
  close FO;
}

open (FI, $fbase . "_sumaclust.clusters" ) or die $!;
open (FO, ">" . $fbase . "_candidate_clusters.fa" ) or die $!;
while(<FI>) {
    chomp;
    @line_cells = split( /\t/ , $_ );
    print FO ">" . $line_cells[0] . "\n";
    print FO $sumaHash{$line_cells[0]} . "\n";
}
close FI;
close FO;

# Create candidateNames/Seqs hash for tview / html (later)
my %candidateNames; my %candidateSeqs;
my $candidate = ""; my $candidatename; my $candidateseq;
open (FI, $fbase . "_candidate_clusters.fa") or die $!;
while(<FI>) {
    chomp;
    if( $_ =~ /^>/ ) {
	if( $candidate ne "" ) {
	    $candidateSeqs{$candidate} = $candidateseq;
	    $candidateNames{$candidate} = $candidatename;
	}
	$candidate = substr $_, 1;
	@line_cells = split( /;/, $candidate );
	$candidatename = $line_cells[0];
	$candidateseq = "";
    } else {
	$candidateseq .= $_;
    }
}
if( $candidate ne "" ) {
    $candidateSeqs{$candidate} = $candidateseq;
    $candidateNames{$candidate} = $candidatename;
}
close FI;

# Alignment of candidate ITD clusters against reference
if( system( sprintf( "%s mem -M -O 6 -T 20 %s %s_candidate_clusters.fa > %s_candidate_clusters_%s.sam", 
	$bwacmd, $refindex , $fbase, $fbase, $alignerExts ) ) ) {
  die "Error aligning extended reads to target FLT3 locus. Exiting...";
}

my %adets = ();
open (FI, $fbase . "_candidate_clusters_" . $alignerExts . ".sam") or die $!;
while(<FI>) {
    if( !/^@/ ) {
	@line_cells = split( /\t/, $_ );

	my %mmlocs = (); my $mmloc = $line_cells[3]-1;
	my $mdstring = substr($line_cells[$mdcell],5); #strip off MD:Z:
	while( $mdstring =~ /^(\d+)([ACGT]|\^[ACGT]+)(.*)/ ) {
	    $mdstring = $3; 
	    if( substr( $2, 0, 1 ) eq "^" ) {
		$mmloc += $1 + length($2) - 1;
	    } else {
		$mmloc += $1 + 1; $mmlocs{$mmloc} = $2; 
	    }
	}

	my $spos = 1; my $rpos = $line_cells[3]; my $hpos = 0; my $i=0;
	$cigar = $line_cells[5];
	print $cigar . "\n";
	while( $cigar =~ /(\d+)([SHDIMN])(.*)/ ) {
	    my $lloc = -1; my $rloc = -1; my $mmkey = ""; my $akey = "";
	    $cigar = $3;
	    if( $2 eq "M" ) {
		if( $i==0 || $cigar eq "" || $1 >= 12 ) {
		    foreach $mmloc ( sort {$a<=>$b} keys %mmlocs ) {
			if( $mmloc >= $rpos && $mmloc < $rpos+$1 ) { 
			    if( length($mmkey) > 0 ) { $mmkey .= "+"; }
			    $mmkey .= $mmloc . ":" . $mmlocs{$mmloc} . ">" . 
				substr( $line_cells[9], $mmloc+$spos-$rpos-$hpos-1, 1);
			}
		    }
		    $adets{$line_cells[0]}{$spos} = {r0=>$rpos+0, r1=>$rpos+$1-1, s0=>$spos+0, s1=>$spos+$1-1, mm=>$mmkey };
		}
		$spos += $1; $rpos += $1;
	    } elsif( $2 eq "I" || $2 eq "S" || $2 eq "H" ) {
		$spos += $1;
		if( $2 eq "H" ) { $hpos += $1; }
	    } elsif( $2 eq "D" ) {
		$rpos += $1;
	    }
	    $i++;
	}
    }
}
close FI;

if( !$debug && (
    system( "rm " . $fbase . "_for_sumaclust.fasta" ) ||
    system( "rm " . $fbase . "_sumaclust.clusters" ) ||
    system( "rm " . $fbase . "_candidate_clusters.fa" ) ||
    system( "rm " . $fbase . "_candidate_clusters_" . $alignerExts . ".sam" ) ) ) {
    die "Error removing candidate_cluster and sumaclust files. Exiting...";
}

if( !$debug ) { system( "rm " . $fbase . "_sumaclust.fasta" ); } # File may not exist if sumaclust fails

open (FO, ">" . $fbase . "_candidate_clusters_middle.fa") or die $!;
foreach( keys %adets ) {
    my @starts = sort {$a<=>$b} keys %{$adets{$_}};
    for( my $i=1; $i<scalar(@starts); $i++ ) {
	if( ($adets{$_}{$starts[$i]}{s0}-$adets{$_}{$starts[$i-1]}{s1}-1) >= 10 ) {
	    print FO ">" . $_ . "|start" . ($i-1) . "\n" . 
		substr( $sumaHash{$_}, $adets{$_}{$starts[$i-1]}{s1}, 
			$adets{$_}{$starts[$i]}{s0}-$adets{$_}{$starts[$i-1]}{s1}-1 ) . "\n";
	}
    }    
}
close FO;

# Align middles (inserts) against FLT3_reference to characterize
if( system( sprintf( "%s mem -k 10 -M -O 6 -T 10 %s %s_candidate_clusters_middle.fa > %s_candidate_clusters_middle_%s.sam", 
	$bwacmd, $refindex , $fbase, $fbase, $alignerMiddle ) ) ) {
  die "Error aligning middles (inserts) against FLT3 target locus. Exiting...";
}

open (FI, $fbase . "_candidate_clusters_middle_" . $alignerMiddle . ".sam") or die $!;
while(<FI>) {
    if( !/^@/ ) {
	@line_cells = split( /\t/, $_ );
	if( ($line_cells[1] & 16) == 16 ) { next; }

	my %mmlocs = (); my $mmloc = $line_cells[3]-1;
	my $mdstring = substr($line_cells[$mdcell],5); #strip off MD:Z:
	while( $mdstring =~ /^(\d+)([ACGT]|\^[ACGT]+)(.*)/ ) {
	    $mdstring = $3; 
	    if( substr( $2, 0, 1 ) eq "^" ) {
		$mmloc += $1 + length($2) - 1;
	    } else {
		$mmloc += $1 + 1; $mmlocs{$mmloc} = $2; 
	    }
	}

	@sub_cells = split( /\|start/, $line_cells[0] );
	my @starts = sort {$a<=>$b} keys %{$adets{$sub_cells[0]}};

	my $spos = $adets{$sub_cells[0]}{$starts[$sub_cells[1]]}{s1}+1; my $rpos = $line_cells[3]; my $hpos = 0; my $i=0;
	$cigar = $line_cells[5];

	while( $cigar =~ /(\d+)([SHDIMN])(.*)/ ) {
	    my $lloc = -1; my $rloc = -1; my $mmkey = ""; my $akey = "";
	    $cigar = $3;
	    if( $2 eq "M" ) {
		foreach $mmloc ( sort {$a<=>$b} keys %mmlocs ) {
		    if( $mmloc >= $rpos && $mmloc < $rpos+$1 ) { 
			if( length($mmkey) > 0 ) { $mmkey .= "+"; }
			$mmkey .= $mmloc . ":" . $mmlocs{$mmloc} . ">" . 
			    substr( $line_cells[9], $mmloc+$spos-$adets{$sub_cells[0]}{$starts[$sub_cells[1]]}{s1}-$rpos-$hpos-1, 1);
		    }
		}
		$adets{$sub_cells[0]}{$spos} = {r0=>$rpos+0, r1=>$rpos+$1-1, s0=>$spos+0, s1=>$spos+$1-1, mm=>$mmkey };
		$spos += $1; $rpos += $1;
	    } elsif( $2 eq "I" || $2 eq "S" || $2 eq "H" ) {
		$spos += $1;
		if( $2 eq "H" ) { $hpos += $1; }
	    } elsif( $2 eq "D" ) {
		$rpos += $1;
	    }
	    $i++;
	}
    }
}
if( !$debug && (
    system( "rm " . $fbase . "_candidate_clusters_middle.fa" ) ||
    system( "rm " . $fbase . "_candidate_clusters_middle_" . $alignerMiddle . ".sam" ) ) ) {
  die "Error removing candidate_clusters_middle (insert) files. Exiting...";
}

my %candidatedets = (); my %itdkeys = (); my %rdotkeys = (); my %cdotkeys = (); my %caltkeys = ();
foreach( keys %adets ) {
    my %cdet0 = proc_alignment_dets_hash( \%{$adets{$_}}, $sumaHash{$_} );
    $candidatedets{$_} = \%cdet0;
    $itdkeys{$_} = $candidatedets{$_}{keyfull};
    $rdotkeys{$_} = $candidatedets{$_}{keyrdot};
    $cdotkeys{$_} = $candidatedets{$_}{keycdot};
    $caltkeys{$_} = $candidatedets{$_}{keycdotalt};
}
if( scalar( keys %itdkeys ) == 0 ) {
    print "NO ITD CANDIDATES. Exiting...\n";
    system( "touch " . $fbase . "_ITD_none.vcf" );
    if( !$debug && ( system( "rm " . $fbase . ".nowt.R*.fastq" ) || 
        ( $inbam ne "" && system( "rm " . $fbase . ".R*.fastq" ) ) ) ) {
      die "Error removing fastq files. Exiting...";
    }
    exit(0);
}

my %nonitdkeys = ();
foreach( keys %itdkeys ) {
    if( $itdkeys{$_} =~ /dup/ || $itdkeys{$_} =~ /itd/ ) { next; }
    $nonitdkeys{$_} = $itdkeys{$_};
}


#### FILTER OUT RELATED ITDS (E.G. PREDICTED SEQUENCING ERRORS)
# Process each ITD size separately
# Then sort by underlying structure with greatest number of extended reads
# Then sort by representative candidate ITD with greatest number of extended reads and fewest mismatches

my %svkeysByBin = ();  # sort into nearest in-frame sizes and remove predicted sequencing errors
my %svkeysByBinTot = ();
my %cnamesByBinMax = ();
foreach( keys %candidatedets ) {
    my $insize = $candidatedets{$_}{totlen} - ($candidatedets{$_}{totlen} % 3);
    if( ($candidatedets{$_}{totlen} % 3) == 2 ) { $insize += 3; }
    @line_cells = split( /\_/, $candidatedets{$_}{keyfull} );
    my $svkey = "";
    for( my $i=0; $i<scalar(@line_cells); $i++ ) {
	@sub_cells = split(/\(/, $line_cells[$i] );
	if( length($svkey)>0 ) { $svkey .= "_"; }
	$svkey .= $sub_cells[0];
    }
    @line_cells = split( /size=/, $_ );
    $svkeysByBin{$insize}{$svkey}{$_} = ($line_cells[1]+0) + 1/($candidatedets{$_}{mismatches}+1);
    $svkeysByBinTot{$insize}{$svkey} += $line_cells[1];

    if( !exists( $cnamesByBinMax{$insize} ) || $line_cells[1] > $cnamesByBinMax{$insize}{max} ||
	( $line_cells[1] == $cnamesByBinMax{$insize}{max} && 
	  $candidatedets{$_}{mismatches} < $cnamesByBinMax{$insize}{mismatches} ) ) {
	$cnamesByBinMax{$insize}{cname} = $_;
	$cnamesByBinMax{$insize}{max} = $line_cells[1];
	$cnamesByBinMax{$insize}{mismatches} = $candidatedets{$_}{mismatches};
    } 
}

my %itdkeysfiltered = ();
foreach my $size ( keys %svkeysByBinTot ) {
    if( $proc_one_itd_per_size ) {
	$itdkeysfiltered{ $cnamesByBinMax{$size}{cname} } = 1;
    } elsif( scalar( keys %{$svkeysByBin{$size}} ) > 1 ) {
	my $sumaclust = ""; my $i=0; my $j=0;
	foreach my $sv ( sort{$svkeysByBinTot{$size}{$b} <=> $svkeysByBinTot{$size}{$a}} keys %{$svkeysByBinTot{$size}} ) {
	    my @sorted = sort{$svkeysByBin{$size}{$sv}{$b} <=> $svkeysByBin{$size}{$sv}{$a}} keys %{$svkeysByBin{$size}{$sv}};
	    $seq = "";
	    @line_cells = split(/\_/, $sv );
	    for( my $k=0; $k<scalar(@line_cells); $k++ ) {
		if( substr($line_cells[$k], 0, 3) eq "ins" ) {
		    $seq .= substr($line_cells[$k], 3); 
		} elsif( substr($line_cells[$k], 0, 3) eq "ref" ) {
		    @sub_cells = split( /\-/, substr($line_cells[$k], 3) );
		    $seq .= substr($refseq, $sub_cells[0]-1, $sub_cells[1]-$sub_cells[0]+1);
		}
	    }
	    for( my $k=0; $k<sum( values %{$svkeysByBin{$size}{$sv}} ); $k++ ) {
		if( ( $j==0 || $extpreHash{$sorted[0]} > 1 || $proc_singleread_itds ) && 
		    ( ($extsizeHash{$sorted[0]} % 3)==0 || $proc_anylength_itds ) ) {
		    $i++;
		    $sumaclust .= ">" . $sorted[0] . ":" . $i . "\n" . $seq . "\n";
		}
		$j++;
	    }
	}
	if( length($sumaclust)>0 ) {
	    open (FO, ">" . $fbase . "_for_sumaclust_len" . $size  . ".fasta" ) or die $!;
	    print FO $sumaclust;
	    close FO;
	    system( $clustercmd . " -d -r -t 5 " . $fbase . "_for_sumaclust_len". $size . ".fasta -O " . 
		    $fbase . "_sumaclust_len" . $size . ".clusters" .  " -F " . 
		    $fbase . "_sumaclust_len" . $size . ".fasta" );

	    open (FI, $fbase . "_sumaclust_len" . $size . ".clusters" ) or die $!;
	    while(<FI>) {
		chomp;
		@line_cells = split( /\t/ , $_ );
		@sub_cells = split( /\:/, $line_cells[0] ); 
		$itdkeysfiltered{ $sub_cells[0] } = 1;
	    }
	    close FI;

	    if( !$debug && (
		system( "rm " . $fbase . "_sumaclust_len" . $size . ".*" ) ||
		system( "rm " . $fbase . "_for_sumaclust_len" . $size . ".fasta" ) ) ) {
		die "Error removing binned sumaclust files. Exiting...";
	    } 
	}
    } else {
	@line_cells = keys %{$svkeysByBin{$size}};
	my $sv = $line_cells[0];
	my @sorted = sort{$svkeysByBin{$size}{$sv}{$b} <=> $svkeysByBin{$size}{$sv}{$a}} keys %{$svkeysByBin{$size}{$sv}};
	if( ($extsizeHash{$sorted[0]} % 3)==0 || $proc_anylength_itds ) { $itdkeysfiltered{ $sorted[0] } = 1; }
    }
}

if( scalar( keys %itdkeysfiltered ) == 0 ) {
    print "NO *FILTERED* ITD CANDIDATES. Exiting...\n";
    system( "touch " . $fbase . "_ITD_none.vcf" );
    if( !$debug && ( system( "rm " . $fbase . ".nowt.R*.fastq" ) || 
        ( $inbam ne "" && system( "rm " . $fbase . ".R*.fastq" ) ) ) ) {
      die "Error removing fastq files. Exiting...";
    }
}


# ITERATE THROUGH ITD CANDIDATES TO PRECLUDE MIXED PAIRS
# E.G. COULD ALTERNATIVELY USE BOWTIE2, WHICH HAS NO-DISCORDANT FLAG
my $itdN = 0; my $headerMut = "";
my %resolvedMutPaired = (); my %psamMut = (); my %rposMut = ();
foreach my $itd ( sort {$extpreHash{$b} <=> $extpreHash{$a} } keys %itdkeysfiltered ) {
  $itdN++;
  open(FO, ">" . $fbase . "_candidate_ITD_" . $itdN . ".fa" ) or die $!;
  print FO ">" . $itd . "\n" . $sumaHash{$itd} . "\n"; 
  close FO;

  # Build bwa candidates index and perform bwa-mem against candidates
  if ( system( sprintf( "%s index -p %s_candidate_ITD_%s %s_candidate_ITD_%s.fa", $bwacmd, $fbase, $itdN, $fbase, $itdN ) ) ||
       system( sprintf( "%s mem -M -O 6 -T 20 %s_candidate_ITD_%s %s.nowt.R1.fastq %s.nowt.R2.fastq > %s_to_candidate_ITD_%s_%s.sam", 
		 $bwacmd, $fbase, $itdN, $fbase, $fbase, $fbase, $itdN, $alignerToITDs ) ) ) {
    die "Error aligning reads against ITD candidates. Exiting...";
  }

  my %resolvedMut = (); my %resolvedMutMJ = (); my %strMut = ();
  open (FI, $fbase . "_to_candidate_ITD_" . $itdN . "_" . $alignerToITDs . ".sam" ) or die $!;
  while(<FI>) {
    if( /^@/ ) {
	if( /^\@SQ/ ) { $headerMut .= $_ };
    } else {
	chomp;
	@line_cells = split( /\t/, $_ );
	if( exists( $resolvedMutPaired{$line_cells[0]} ) || $line_cells[2] eq "*" ||
	    ( $line_cells[1] & 256 ) == 256 || ( $line_cells[1] & 2048 ) == 2048 ) {
	  next; 
	}
	if( ( $line_cells[1] & 64 ) == 64 ) {
	    $readsuffix = "_1";
	} elsif( ( $line_cells[1] & 128 ) == 128 ) {
	    $readsuffix = "_2";
	} else {
	    next;
	}

	$readname = $line_cells[0] . $readsuffix;
        $cigar = $line_cells[5];
	my $clipLen = 0;
	my $insLen = 0;
	my $delLen = 0;
	if( @a = $cigar =~ /(\d+)S/g ) { $clipLen = max(@a); }
	if( @a = $cigar =~ /(\d+)I/g ) { foreach my $x (@a) { $insLen += $x; } }
	if( @a = $cigar =~ /(\d+)D/g ) { foreach my $x (@a) { $delLen += $x; } }

	if( $insLen <= $maxInsAllowed && $delLen <= $maxDelAllowed && $clipLen <= $maxClipAllowed && 
		 ($_ =~ /NM:i:(\d+)/ && $1 <= $maxEditDist) ) {
	  my $rpos = $line_cells[3]-1;
	  while( $cigar =~ /(\d+)([SHDIMN])(.*)/ ) {
	    $cigar = $3;
	    if( $2 eq "M" || $2 eq "N" || $2 eq "D" ) { $rpos += $1; }
	  }
	  $resolvedMut{$readname} = 1;  # reads aligning "stringently" to itd genome 
	  $psamMut{$readname} = $_;
	  $rposMut{$readname} = $rpos;
	  if( ( $line_cells[1] & 16 ) == 16 ) {
	    $strMut{$readname} = "-";
	  } else {
	    $strMut{$readname} = "+";
	  }
	  if( ( $candidatedets{$itd}{mutbk_s0} >= $line_cells[3]+$buffer &&
		$candidatedets{$itd}{mutbk_s0} < $rpos-$buffer ) ||
    	      ( $candidatedets{$itd}{mutbk_s1} > $line_cells[3]+$buffer && 
		$candidatedets{$itd}{mutbk_s1} <= $rpos-$buffer ) ) {
	    $resolvedMutMJ{$readname} = 1;
	  }
        }
    }
  }
  close FI;

  foreach my $read (keys %resolvedMut) {
    $read =~ s/\_1$//g; $read =~ s/\_2$//g;
    if( exists($resolvedMut{$read."_1"}) && exists($resolvedMut{$read."_2"}) &&
	$strMut{$read."_1"} ne $strMut{$read."_2"} && 
        ($resolvedMutMJ{$read."_1"} || $resolvedMutMJ{$read."_2"} ) ) {
	$resolvedMutPaired{$read} = 1;
	$readsPassed{$read} = 1;
    }
  }

  if( !$debug && (
    system( "rm " . $fbase . "_candidate_ITD_" . $itdN . ".fa" ) ||
    system( "rm " . $fbase . "_candidate_ITD_" . $itdN . ".amb" ) ||
    system( "rm " . $fbase . "_candidate_ITD_" . $itdN . ".ann" ) ||
    system( "rm " . $fbase . "_candidate_ITD_" . $itdN . ".bwt" ) ||
    system( "rm " . $fbase . "_candidate_ITD_" . $itdN . ".pac" ) ||
    system( "rm " . $fbase . "_candidate_ITD_" . $itdN . ".sa" ) ||
    system( "rm " . $fbase . "_to_candidate_ITD_" . $itdN . "_" . $alignerToITDs . ".sam" ) ) ) {
    die "Error removing candidate_ITDs fasta, index, and alignment files. Exiting...";
  }
}

if( !$debug && ( system( "rm " . $fbase . ".nowt.R*.fastq" ) || 
  ( $inbam ne "" && system( "rm " . $fbase . ".R*.fastq" ) ) ) ) {
  die "Error removing fastq files. Exiting...";
}

open (FO, ">" . $fbase . "_to_candidate_ITDs.filtered.sam" ) or die $!;
print FO $headerMut;
foreach( sort keys %resolvedMutPaired ) {
    print FO $psamMut{$_."_1"};
    if( exists( $umitags{$_}) ) { print FO "\tMQ:i:60\t" . $umitags{$_}; } 
    print FO "\n";
    print FO $psamMut{$_."_2"};
    if( exists( $umitags{$_}) ) { print FO "\tMQ:i:60\t" . $umitags{$_}; }
    print FO "\n";
}
close FO;

if( $umitag ne "" ) {
  if( system( sprintf( "%s view -bS %s_to_candidate_ITDs.filtered.sam | " .
                       "%s -jar %s SortSam I=/dev/stdin SO=queryname O=/dev/stdout | " .
                       "%s -jar %s SetMateInformation -i /dev/stdin -o %s_mut.bam",
              $samcmd, $fbase, $javacmd, $picardjar, $javacmd, $fgbiojar, $fbase ) ) ) {
    die "Failed to convert mutant sam to bam. Exiting...\n";
  }

  if( system( sprintf( "%s -jar %s GroupReadsByUmi -i %s_mut.bam -o %s_mut_umigroups.bam -s %s -f %s_mut_umisizes.txt -m 0",
		       $javacmd, $fgbiojar, $fbase, $fbase, $strat, $fbase ) ) ) {
    die "Failed to perform fgbio GroupReadsByUmi on mutant reads...\n";
  }

  if( system( sprintf ( "%s view %s_mut_umigroups.bam > %s_mut_umigroups.sam",
    $samcmd, $fbase, $fbase ) ) ) {
    die "Failed to convert mut umigroups bam to sam...\n";
  }

  open( FI, $fbase . "_mut_umigroups.sam" ) or die $!;
  while(<FI>) {
    @line_cells = split(/\t/, $_);
    if( $_ =~ /MI:Z:(\d+)/ ) { $umisMut{$line_cells[0]} = $1; }
  }

  if( !$debug &&
      ( system( sprintf( "rm %s_mut.bam", $fbase ) ) || 
        system( sprintf( "rm %s_mut_umisizes.txt", $fbase ) ) ||
        system( sprintf( "rm %s_mut_umigroups.sam", $fbase ) ) ||
        system( sprintf( "rm %s_mut_umigroups.bam", $fbase ) ) ) ) {
    die "Failed to remove mutant bam file...\n";
  }
}

my %readsMut = (); my %readsMutLeft = (); my %readsMutRight = ();
my %readsMutWt = (); my %readsMutWtLeft = (); my %readsMutWtRight = ();
my %readsMutMut = (); my %readsMutMutLeft = (); my %readsMutMutRight = ();
my %readsWt = (); my %readsWtLeft = (); my %readsWtRight = ();
my %coverageMut = (); my %coverageMutLeft = (); my %coverageMutRight = (); 
my %coverageMutWt = (); my %coverageMutWtLeft = (); my %coverageMutWtRight = (); 
my %coverageWt = ();  my %coverageWtLeft = ();  my %coverageWtRight = (); 
my %readsMutByR2 = (); my %readsWtByR2 = (); 

my %umisMutLeft = (); my %umisMutRight = (); my %umisMutAll = ();
my %umisMutWtLeft = (); my %umisMutWtRight = ();
my %umisMutMutLeft = (); my %umisMutMutRight = (); my %umisMutMutAll= ();
my %umisWtLeft = (); my %umisWtRight = (); my %umisWtAll = ();
my %covUmisWtLeft = (); my %covUmisWtRight = (); my %covUmisWt = ();
my %covUmisMutLeft = (); my %covUmisMutRight = (); my %covUmisMut = ();
my %covUmisMutWtLeft = (); my %covUmisMutWtRight = ();
my $lposAdj; my $rposAdj;

open( FI, $fbase . "_to_candidate_ITDs.filtered.sam" ) or die $!;
while(<FI>) {
    if( /^@/ ) { next; }
    @line_cells = split( /\t/, $_ );

    $readname = $line_cells[0];
    my $cname = $line_cells[2];
    my $rpos=0;
    if( ( $line_cells[1] & 64 ) == 64 ) {
      $rpos = $rposMut{$readname."_1"};
    } elsif( ( $line_cells[1] & 128 ) == 128 ) {
      $rpos = $rposMut{$readname."_2"};
    } else {
      next;
    }
 
    # ADD LOGIC FOR NOBIND CASES (ALLOW RELAXED BUFFER?)

    if( $candidatedets{$cname}{mutbk_s0} >= $line_cells[3]+$buffer && 
	$candidatedets{$cname}{mutbk_s0} < $rpos-$buffer ) {
      if( ( $line_cells[1] & 128 ) == 128 ) { $readsMutByR2{$cname}{$readname} = 1; }
      if( !exists($readsMutLeft{$cname}{$readname}) ) {
	$readsMutLeft{$cname}{$readname} = 1;
	$readsMut{$cname}{$readname} = 1;
	if( exists($umisMut{$readname}) ) { 
	  $umisMutLeft{$cname}{$umisMut{$readname}} = 1;
	  $umisMutAll{$cname}{$umisMut{$readname}} = 1;
	} else {
	  $umisMutLeft{$cname}{other} = 1;
	  $umisMutAll{$cname}{other} = 1;
	}
      }
    }
    if( $candidatedets{$cname}{mutbk_s1} > $line_cells[3]+$buffer && 
	$candidatedets{$cname}{mutbk_s1} <= $rpos-$buffer ) {
      if( ( $line_cells[1] & 128 ) == 128 ) { $readsMutByR2{$cname}{$readname} = 1; }
      if( !exists($readsMutRight{$cname}{$readname}) ) {
	$readsMutRight{$cname}{$readname} = 1;
	$readsMut{$cname}{$readname} = 1;
	if( exists($umisMut{$readname}) ) {
	  $umisMutRight{$cname}{$umisMut{$readname}} = 1;
	  $umisMutAll{$cname}{$umisMut{$readname}} = 1;
	} else {
	  $umisMutRight{$cname}{other} = 1;
	  $umisMutAll{$cname}{other} = 1;
	}
      }
    }
    if( !exists($readsMutWtLeft{$cname}{$readname}) && 
	$candidatedets{$cname}{dup_r0} > $line_cells[3]+$buffer && 
	$candidatedets{$cname}{dup_r0} <= $rpos-$buffer ) {
	$readsMutWtLeft{$cname}{$readname} = 1;
	$readsMutWt{$cname}{$readname} = 1;
	if( exists($umisMut{$readname}) ) {
	  $umisMutWtLeft{$cname}{$umisMut{$readname}} = 1;
	} else {
	  $umisMutWtLeft{$cname}{other} = 1;
	}
    }
    if( !exists($readsMutWtRight{$cname}{$readname}) && 
	$candidatedets{$cname}{dup_r1}+$candidatedets{$cname}{itdsize} >= $line_cells[3]+$buffer && 
	$candidatedets{$cname}{dup_r1}+$candidatedets{$cname}{itdsize} < $rpos-$buffer ) {
	$readsMutWtRight{$cname}{$readname} = 1;
	$readsMutWt{$cname}{$readname} = 1;
	if( exists($umisMut{$readname}) ) {
	  $umisMutWtRight{$cname}{$umisMut{$readname}} = 1;
	} else {
	  $umisMutWtRight{$cname}{other} = 1;
	}
    }

    if( $rpos <= $candidatedets{$cname}{dup_r1} ) {
      $rposAdj = $rpos;
      $lposAdj = $line_cells[3];
    } elsif( $line_cells[3] >= $candidatedets{$cname}{dup_r1} ) {
      $rposAdj = $rpos - $candidatedets{$cname}{itdsize};
      $lposAdj = $line_cells[3] - $candidatedets{$cname}{itdsize};
    } else {
      $rposAdj = $candidatedets{$cname}{dup_r1};
      $lposAdj = $candidatedets{$cname}{dup_r0};
    }

    # Determine if mutant read also covers other mutant boundaries for depth purposes
    if( exists( $readsMut{$cname}{$readname} ) ) {
      foreach my $itd ( keys %itdkeysfiltered ) {
	if( $itd ne $cname ) {
          if( !exists($readsMutMutLeft{$cname}{$itd}{$readname}) && 
	    $candidatedets{$itd}{dup_r0} > $lposAdj+$buffer && 
	    $candidatedets{$itd}{dup_r0} <= $rposAdj-$buffer ) {
	    $readsMutMutLeft{$cname}{$itd}{$readname} = 1;
	    $readsMutMut{$cname}{$itd}{$readname} = 1;
	    if( exists($umisMut{$readname}) ) {
	      $umisMutMutLeft{$cname}{$itd}{$umisMut{$readname}} = 1;
	      $umisMutMutAll{$cname}{$itd}{$umisMut{$readname}} = 1;
	    } else {
	      $umisMutMutLeft{$cname}{$itd}{other} = 1;
	      $umisMutMutAll{$cname}{$itd}{other} = 1;
	    }
	  }
	  if( !exists($readsMutMutRight{$cname}{$itd}{$readname}) &&
	    $candidatedets{$cname}{dup_r1} >= $lposAdj+$buffer && 
	    $candidatedets{$cname}{dup_r1} < $rposAdj-$buffer ) {
	    $readsMutMutRight{$cname}{$itd}{$readname} = 1;
	    $readsMutMut{$cname}{$itd}{$readname} = 1;
	    if( exists($umisMut{$readname}) ) {
	      $umisMutMutRight{$cname}{$itd}{$umisMut{$readname}} = 1;
	      $umisMutMutAll{$cname}{$itd}{$umisMut{$readname}} = 1;
	    } else {
	      $umisMutMutRight{$cname}{$itd}{other} = 1;
	      $umisMutMutAll{$cname}{$itd}{other} = 1;
	    }
	  }
	}
      }
    }
}
close FI;

my %baitHash = (); my $bait = ""; my @baitsAll = keys %baitHash;
my %baitIssues = (); my %baitDists = (); my %baitsByRead = ();
my %countsByBait = (); my %umisByBait = (); my %umisFlag = ();
my %countsByBaitAdj = (); my %umisByBaitAdj =();
my @lens = sort {$readlengths{$b} <=> $readlengths{$a}} keys %readlengths;
my $readlength = $lens[0];
if( $probes ne "" ) {
  open( FI, $probes . ".fa" ) or die $!;
  while(<FI>) {
    chomp;
    if (/\>/) {
      $bait = substr($_,1); 
    } else {
      $baitHash{$bait}{seq} = $_;
      $baitHash{$bait}{len} = length($_);
      $i = index( $refseq, $_ );
      if( $i >= 0 ) {
	$baitHash{$bait}{dir} = "fwd";
	$baitHash{$bait}{start} = $i + 1;
	$baitHash{$bait}{end} = $i + length($_);
      } else {
	$seq = reverse $_; $seq =~ tr/ATGCatgc/TACGtacg/;
	$i = index( $refseq, $seq );
	if( $i >= 0 ) {
	  $baitHash{$bait}{dir} = "rev";
	  $baitHash{$bait}{start} = $i + 1;
	  $baitHash{$bait}{end} = $i + length($_);
	}
      }
    }
  }
  close FI;
  @baitsAll = sort keys %baitHash; push @baitsAll, "other";
    
  foreach( keys %itdkeysfiltered ) {
    foreach my $bait ( keys %baitHash ) {
      if( $baitHash{$bait}{dir} eq "fwd" ) {
  	  $baitDists{$_}{$bait} = ($baitHash{$bait}{len}) .
	      "_" . ($candidatedets{$_}{dup_r0}-$baitHash{$bait}{start}) . 
	      "_" . ($candidatedets{$_}{dup_r1}-$baitHash{$bait}{start}+1) . "bp"; 
      } elsif( $baitHash{$bait}{dir} eq "rev" ) {
	  $baitDists{$_}{$bait} = ($baitHash{$bait}{len}) .
	      "_" . ($baitHash{$bait}{end}-$candidatedets{$_}{dup_r1}) . 
	      "_" . ($baitHash{$bait}{end}-$candidatedets{$_}{dup_r0}+1) . "bp"; 
      }

      if( $candidatedets{$_}{dup_r0} <= $baitHash{$bait}{start} && 
	  $candidatedets{$_}{dup_r1} >= $baitHash{$bait}{end} ) {
	  $baitIssues{$_}{dupbait}{$bait} = 1;
	  $baitDists{$_}{$bait} .= "_DB";
      } elsif( ( $baitHash{$bait}{dir} eq "fwd" &&
		 $candidatedets{$_}{dup_r1} < $baitHash{$bait}{end} && 
		 $candidatedets{$_}{dup_r1} >= $baitHash{$bait}{start} ) ||
	       ( $baitHash{$bait}{dir} eq "rev" &&
		 $candidatedets{$_}{dup_r0} > $baitHash{$bait}{start} && 
		 $candidatedets{$_}{dup_r0} <= $baitHash{$bait}{end} ) ) {
	  $baitIssues{$_}{nobind}{$bait} = 1;
	  $baitDists{$_}{$bait} .= "_NB";
      }

      if( ( $baitHash{$bait}{dir} eq "fwd" &&
	    $candidatedets{$_}{dup_r1} >= $baitHash{$bait}{end} &&
	    $candidatedets{$_}{dup_r1} <  $baitHash{$bait}{start}+$readlength-$buffer ) ||
	  ( $baitHash{$bait}{dir} eq "rev" &&
	    $candidatedets{$_}{dup_r0} <= $baitHash{$bait}{start} &&
	    $candidatedets{$_}{dup_r0} >  $baitHash{$bait}{end}-$readlength+$buffer ) ) {
	  $baitIssues{$_}{unambig}{$bait} = 1;  # these baits sequence through entire duplicated sequence (assming not a duplicated bait)
	  $baitDists{$_}{$bait} .= "_RB";
      }
    }
  }
    
  open( FO, ">" . $fbase . "_to_baits.fa" ) or die $!;
  foreach $readname ( sort keys %readsPassed ) {
    my $seq = $pseq{$readname . "_2"}; # uses info from alignments to reference FLT3 but should be okay
    if( $pstr{$readname . "_2"} eq "-" ) {
	$seq = reverse $seq; $seq =~ tr/ATGCatgc/TACGtacg/;
    }
    print FO ">" . $readname . "\n" . substr( $seq, 0, 49 ) . "\n";
  }
  close FO;

  if( system( sprintf( "%s mem -k %s -M -O 6 -T %s %s %s_to_baits.fa > %s_to_baits.sam",
		       $bwacmd, $bwaSeedLength, $bwaThreshold, $probes, $fbase, $fbase ) ) ) {
    die "Failed to locally align to baits. Exiting...\n";
  }

  open( FI, $fbase . "_to_baits.sam" ) or die $!;
  while(<FI>) {
    if( /AS:i:(\d+)/ && $1 >= 40 ) {
      @line_cells = split( /\t/, $_ );
      $baitsByRead{$line_cells[0]} = $line_cells[2];
    }
  }

  if( !$debug && system( "rm " . $fbase . "_to_baits.*" ) ) {
    die "Error removing baits by read files. Exiting...";
  }
}

foreach $readname ( sort keys %resolvedWtPaired ) {
  foreach $readsuffix ( ( "_1", "_2" ) ) {
    @line_cells = split( /\t/, $psam{$readname . $readsuffix} );
    $cigar = $line_cells[5];
    my $rpos = $line_cells[3] - 1;
    while( $cigar =~ /(\d+)([SHDIMN])(.*)/ ) {
      $cigar = $3;
      if( $2 eq "M" || $2 eq "N" || $2 eq "D" ) { $rpos += $1; }
    }

    foreach my $cname ( keys %itdkeysfiltered ) {
      if( !exists( $readsWtLeft{$cname}{$readname} ) &&
	  $candidatedets{$cname}{dup_r0} > $line_cells[3]+$buffer && 
	  $candidatedets{$cname}{dup_r0} <= $rpos - $buffer ) {
	$readsWtLeft{$cname}{$readname} = 1;
	$readsWt{$cname}{$readname} = 1;
	if( exists($umisWt{$readname}) ) {
	  $umisWtLeft{$cname}{$umisWt{$readname}} = 1;
	  $umisWtAll{$cname}{$umisWt{$readname}} = 1;
	} else {
	  $umisWtLeft{$cname}{other} = 1;
	  $umisWtAll{$cname}{other} = 1;
	}
      }
      if( !exists( $readsWtRight{$cname}{$readname} ) &&
	  $candidatedets{$cname}{dup_r1} >= $line_cells[3]+$buffer && 
	  $candidatedets{$cname}{dup_r1} < $rpos - $buffer ) {
	$readsWtRight{$cname}{$readname} = 1;
	$readsWt{$cname}{$readname} = 1;
	if( exists($umisWt{$readname}) ) {
	  $umisWtRight{$cname}{$umisWt{$readname}} = 1;
	  $umisWtAll{$cname}{$umisWt{$readname}} = 1;
	} else {
	  $umisWtRight{$cname}{other} = 1;
	  $umisWtAll{$cname}{other} = 1;
	}
      }

      if( $probes ne "" && ($line_cells[1] & 128) == 128 && exists($baitsByRead{$readname}) ) {
	my $bait = $baitsByRead{$readname};
	if( exists($baitIssues{$cname}{nobind}{$bait}) ||
	    ( exists($baitIssues{$cname}{unambig}{$bait}) &&
	      ( length($line_cells[9]) == $readlength ||
		( $baitHash{$bait}{dir} eq "fwd" &&
		  #$candidatedets{$cname}{dup_r1} >= $baitHash{$bait}{end} &&
		  $candidatedets{$cname}{dup_r1} < 
		  $baitHash{$bait}{start}+length($line_cells[9])-$buffer ) ||
		( $baitHash{$bait}{dir} eq "rev" &&
		  #$candidatedets{$cname}{dup_r0} <= $baitHash{$bait}{start} &&
		  $candidatedets{$cname}{dup_r0} >
		  $baitHash{$bait}{end}-length($line_cells[9])+$buffer ) ) ) ) {
	  $readsWtByR2{$cname}{$readname} = 1; 
	} 
      }
    }
  }
}

if( $probes ne "" ) {
foreach ( keys %itdkeysfiltered ) {
  foreach my $bait ( @baitsAll ) {
    $countsByBait{$_}{WtLeft}{$bait} = 0; $countsByBait{$_}{WtRight}{$bait} = 0;
    $countsByBait{$_}{MutLeft}{$bait} = 0; $countsByBait{$_}{MutRight}{$bait} = 0;
    $countsByBait{$_}{MutWtLeft}{$bait} = 0; $countsByBait{$_}{MutWtRight}{$bait} = 0;
    $countsByBait{$_}{MutMutLeft}{$bait} = 0; $countsByBait{$_}{MutMutRight}{$bait} = 0;
    $countsByBait{$_}{WtR2}{$bait} = 0; $countsByBait{$_}{MutR2}{$bait} = 0; 
    $countsByBait{$_}{Mut}{$bait} = 0; $countsByBait{$_}{All}{$bait} = 0; 
    $umisByBait{$_}{WtLeft}{$bait} = 0; $umisByBait{$_}{WtRight}{$bait} = 0;
    $umisByBait{$_}{MutLeft}{$bait} = 0; $umisByBait{$_}{MutRight}{$bait} = 0;
    $umisByBait{$_}{MutWtLeft}{$bait} = 0; $umisByBait{$_}{MutWtRight}{$bait} = 0;
    $umisByBait{$_}{MutMutLeft}{$bait} = 0; $umisByBait{$_}{MutMutRight}{$bait} = 0;
    $umisByBait{$_}{WtR2}{$bait} = 0; $umisByBait{$_}{MutR2}{$bait} = 0; 
    $umisByBait{$_}{Mut}{$bait} = 0; $umisByBait{$_}{All}{$bait} = 0; 
  }
  foreach my $read ( keys %{$readsWtLeft{$_}} ) {
    if( exists( $baitsByRead{$read} ) ) {
      $countsByBait{$_}{WtLeft}{$baitsByRead{$read}} += 1;
      if( exists($umisWt{$read}) ) {
	if( !exists($umisFlag{$_}{WtLeft}{$umisWt{$read}}) ) {
	  $umisFlag{$_}{WtLeft}{$umisWt{$read}} = 1;
	  $umisByBait{$_}{WtLeft}{$baitsByRead{$read}} += 1;
	}
      }
    } else {
      $countsByBait{$_}{WtLeft}{other} += 1;
      if( exists($umisWt{$read}) ) {
	if( !exists($umisFlag{$_}{WtLeft}{$umisWt{$read}}) ) {
	  $umisFlag{$_}{WtLeft}{$umisWt{$read}} = 1;
	  $umisByBait{$_}{WtLeft}{other} += 1;
	}
      }
    }
  }
  foreach my $read ( keys %{$readsWtRight{$_}} ) {
    if( exists( $baitsByRead{$read} ) ) {
      $countsByBait{$_}{WtRight}{$baitsByRead{$read}} += 1;
      if( exists($umisWt{$read}) ) {
	if( !exists($umisFlag{$_}{WtRight}{$umisWt{$read}}) ) {
	  $umisFlag{$_}{WtRight}{$umisWt{$read}} = 1;
	  $umisByBait{$_}{WtRight}{$baitsByRead{$read}} += 1;
	}
      }
    } else {
      $countsByBait{$_}{WtRight}{other} += 1;
      if( exists($umisWt{$read}) ) {
	if( !exists($umisFlag{$_}{WtRight}{$umisWt{$read}}) ) {
	  $umisFlag{$_}{WtRight}{$umisWt{$read}} = 1;
	  $umisByBait{$_}{WtRight}{other} += 1;
	}
      }
    }
  }
  foreach my $read ( keys %{$readsMutWtLeft{$_}} ) {
    if( exists( $baitsByRead{$read} ) ) {
      $countsByBait{$_}{MutWtLeft}{$baitsByRead{$read}} += 1;
      if( exists($umisMut{$read}) ) {
	if( !exists($umisFlag{$_}{MutWtLeft}{$umisMut{$read}}) ) {
	  $umisFlag{$_}{MutWtLeft}{$umisMut{$read}} = 1;
	  $umisByBait{$_}{MutWtLeft}{$baitsByRead{$read}} += 1;
	}
      }
    } else {
      $countsByBait{$_}{MutWtLeft}{other} += 1;
      if( exists($umisMut{$read}) ) {
	if( !exists($umisFlag{$_}{MutWtLeft}{$umisMut{$read}}) ) {
	  $umisFlag{$_}{MutWtLeft}{$umisMut{$read}} = 1;
	  $umisByBait{$_}{MutWtLeft}{other} += 1;
	}
      }
    }
  }
  foreach my $read ( keys %{$readsMutWtRight{$_}} ) {
    if( exists( $baitsByRead{$read} ) ) {
      $countsByBait{$_}{MutWtRight}{$baitsByRead{$read}} += 1;
      if( exists($umisMut{$read}) ) {
	if( !exists($umisFlag{$_}{MutWtRight}{$umisMut{$read}}) ) {
	  $umisFlag{$_}{MutWtRight}{$umisMut{$read}} = 1;
	  $umisByBait{$_}{MutWtRight}{$baitsByRead{$read}} += 1;
	}
      }
    } else {
      $countsByBait{$_}{MutWtRight}{other} += 1;
      if( exists($umisMut{$read}) ) {
	if( !exists($umisFlag{$_}{MutWtRight}{$umisMut{$read}}) ) {
	  $umisFlag{$_}{MutWtRight}{$umisMut{$read}} = 1;
	  $umisByBait{$_}{MutWtRight}{other} += 1;
	}
      }
    }
  }
  foreach my $read ( keys %{$readsMutLeft{$_}} ) {
    if( exists( $baitsByRead{$read} ) ) {
      $countsByBait{$_}{MutLeft}{$baitsByRead{$read}} += 1;
      if( exists($umisMut{$read}) ) {
	if( !exists($umisFlag{$_}{MutLeft}{$umisMut{$read}}) ) {
	  $umisFlag{$_}{MutLeft}{$umisMut{$read}} = 1;
	  $umisByBait{$_}{MutLeft}{$baitsByRead{$read}} += 1;
	}
      }
    } else {
      $countsByBait{$_}{MutLeft}{other} += 1;
      if( exists($umisMut{$read}) ) {
	if( !exists($umisFlag{$_}{MutLeft}{$umisMut{$read}}) ) {
	  $umisFlag{$_}{MutLeft}{$umisMut{$read}} = 1;
	  $umisByBait{$_}{MutLeft}{other} += 1;
	}
      }
    }
  }
  foreach my $read ( keys %{$readsMutRight{$_}} ) {
    if( exists( $baitsByRead{$read} ) ) {
      $countsByBait{$_}{MutRight}{$baitsByRead{$read}} += 1;
      if( exists($umisMut{$read}) ) {
	if( !exists($umisFlag{$_}{MutRight}{$umisMut{$read}}) ) {
	  $umisFlag{$_}{MutRight}{$umisMut{$read}} = 1;
	  $umisByBait{$_}{MutRight}{$baitsByRead{$read}} += 1;
	}
      }
    } else {
      $countsByBait{$_}{MutRight}{other} += 1;
      if( exists($umisMut{$read}) ) {
	if( !exists($umisFlag{$_}{MutRight}{$umisMut{$read}}) ) {
	  $umisFlag{$_}{MutRight}{$umisMut{$read}} = 1;
	  $umisByBait{$_}{MutRight}{other} += 1;
	}
      }
    }
  }
  foreach my $read ( sort keys %{$readsWtByR2{$_}} ) {
    if( exists( $baitsByRead{$read} ) ) {
      $countsByBait{$_}{WtR2}{$baitsByRead{$read}} += 1;
      if( exists($umisWt{$read}) ) {
	if( !exists($umisFlag{$_}{WtR2}{$umisWt{$read}}) ) {
	  $umisFlag{$_}{WtR2}{$umisWt{$read}} = 1;
	  $umisByBait{$_}{WtR2}{$baitsByRead{$read}} += 1;
	}
      }
    } else {
      $countsByBait{$_}{WtR2}{other} += 1;
      if( exists($umisWt{$read}) ) {
	if( !exists($umisFlag{$_}{WtR2}{$umisWt{$read}}) ) {
	  $umisFlag{$_}{WtR2}{$umisWt{$read}} = 1;
	  $umisByBait{$_}{WtR2}{other} += 1;
	}
      }
    }
  }
  foreach my $read ( keys %{$readsMutByR2{$_}} ) {
    if( exists( $baitsByRead{$read} ) ) {
      $countsByBait{$_}{MutR2}{$baitsByRead{$read}} += 1;
      if( exists($umisMut{$read}) ) {
	if( !exists($umisFlag{$_}{MutR2}{$umisMut{$read}}) ) {
	  $umisFlag{$_}{MutR2}{$umisMut{$read}} = 1;
	  $umisByBait{$_}{MutR2}{$baitsByRead{$read}} += 1;
	}
      }
    } else {
      $countsByBait{$_}{MutR2}{other} += 1;
      if( exists($umisMut{$read}) ) {
	if( !exists($umisFlag{$_}{MutR2}{$umisMut{$read}}) ) {
	  $umisFlag{$_}{MutR2}{$umisMut{$read}} = 1;
	  $umisByBait{$_}{MutR2}{other} += 1;
	}
      }
    }
  }

  foreach my $cname ( keys %itdkeysfiltered ) {
    if( $cname ne $_ ) { 
      if( exists( $readsMutMutLeft{$cname}{$_} ) ) {
	foreach my $read ( keys %{$readsMutMutLeft{$cname}{$_}} ) {
	  if( exists( $baitsByRead{$read} ) ) {
	    $countsByBait{$_}{MutMutLeft}{$baitsByRead{$read}} += 1;
	    if( exists($umisMut{$read}) ) {
	      if( !exists($umisFlag{$_}{MutMutLeft}{$umisMut{$read}}) ) {
	        $umisFlag{$_}{MutMutLeft}{$umisMut{$read}} = 1;
	        $umisByBait{$_}{MutMutLeft}{$baitsByRead{$read}} += 1;
	      }
	    }
	  } else {
	    $countsByBait{$_}{MutMutLeft}{other} += 1;
	    if( exists($umisMut{$read}) ) {
	      if( !exists($umisFlag{$_}{MutMutLeft}{$umisMut{$read}}) ) {
	        $umisFlag{$_}{MutMutLeft}{$umisMut{$read}} = 1;
	        $umisByBait{$_}{MutMutLeft}{other} += 1;
	      }
	    }
	  }
        }
      }
      if( exists( $readsMutMutRight{$cname}{$_} ) ) {
	foreach my $read ( keys %{$readsMutMutRight{$cname}{$_}} ) {
	  if( exists( $baitsByRead{$read} ) ) {
	    $countsByBait{$_}{MutMutRight}{$baitsByRead{$read}} += 1;
	    if( exists($umisMut{$read}) ) {
	      if( !exists($umisFlag{$_}{MutMutRight}{$umisMut{$read}}) ) {
	        $umisFlag{$_}{MutMutRight}{$umisMut{$read}} = 1;
	        $umisByBait{$_}{MutMutRight}{$baitsByRead{$read}} += 1;
	      }
	    }
	  } else {
	    $countsByBait{$_}{MutMutRight}{other} += 1;
	    if( exists($umisMut{$read}) ) {
	      if( !exists($umisFlag{$_}{MutMutRight}{$umisMut{$read}}) ) {
	        $umisFlag{$_}{MutMutRight}{$umisMut{$read}} = 1;
	        $umisByBait{$_}{MutMutRight}{other} += 1;
	      }
	    }
	  }
        }
      }
    }
  }

  foreach my $bait ( @baitsAll ) {
    $countsByBaitAdj{$_}{Mut}{$bait} = 0; $umisByBaitAdj{$_}{Mut}{$bait} = 0;
    $countsByBaitAdj{$_}{All}{$bait} = 0; $umisByBaitAdj{$_}{All}{$bait} = 0;
    if( exists( $baitIssues{$_}{unambig}{$bait} ) ||
        exists( $baitIssues{$_}{nobind}{$bait} ) ) {
      $countsByBait{$_}{Mut}{$bait} = $countsByBait{$_}{MutR2}{$bait} + 0;
      $countsByBait{$_}{All}{$bait} = $countsByBait{$_}{MutR2}{$bait} + $countsByBait{$_}{WtR2}{$bait};
      $umisByBait{$_}{Mut}{$bait} = $umisByBait{$_}{MutR2}{$bait} + 0;
      $umisByBait{$_}{All}{$bait} = $umisByBait{$_}{MutR2}{$bait} + $umisByBait{$_}{WtR2}{$bait};
      # Simple adjustment for baits with multiple potential binding positions
      if( exists( $baitIssues{$_}{dupbait}{$bait} ) ) {
	$countsByBaitAdj{$_}{Mut}{$bait} = $countsByBait{$_}{Mut}{$bait} + 0;
	$umisByBaitAdj{$_}{Mut}{$bait} = $umisByBait{$_}{Mut}{$bait} + 0;
      } elsif( exists( $baitIssues{$_}{nobind}{$bait} ) ) {
	$countsByBaitAdj{$_}{Mut}{$bait} = -$countsByBait{$_}{Mut}{$bait} + 0;
	$countsByBaitAdj{$_}{All}{$bait} = -$countsByBait{$_}{All}{$bait} + 0;
	$umisByBaitAdj{$_}{Mut}{$bait} = -$umisByBait{$_}{Mut}{$bait} + 0;
	$umisByBaitAdj{$_}{All}{$bait} = -$umisByBait{$_}{All}{$bait} + 0;
      }
    } else {
      $countsByBait{$_}{Mut}{$bait} = ( $countsByBait{$_}{MutLeft}{$bait} +
                                        $countsByBait{$_}{MutRight}{$bait} ) / 2;
      $umisByBait{$_}{Mut}{$bait} = ( $umisByBait{$_}{MutLeft}{$bait} +
                                      $umisByBait{$_}{MutRight}{$bait} ) / 2;

      if( exists( $baitIssues{$_}{dupbait}{$bait} ) ) {
        # for dup baits, one of WtLeft or WtRight will be 0
        $countsByBait{$_}{All}{$bait} = ( $countsByBait{$_}{WtLeft}{$bait} +
					  $countsByBait{$_}{WtRight}{$bait} +
					  $countsByBait{$_}{MutLeft}{$bait} / 2 +
                                          $countsByBait{$_}{MutRight}{$bait} / 2 ); 
        $umisByBait{$_}{All}{$bait} = ( $umisByBait{$_}{WtLeft}{$bait} +
				        $umisByBait{$_}{WtRight}{$bait} +
				        $umisByBait{$_}{MutLeft}{$bait} / 2 +
                                        $umisByBait{$_}{MutRight}{$bait} / 2 ); 
        $countsByBaitAdj{$_}{Mut}{$bait} = $countsByBait{$_}{Mut}{$bait} + 0;
        $umisByBaitAdj{$_}{Mut}{$bait} = $umisByBait{$_}{Mut}{$bait} + 0;
      } else {
        $countsByBait{$_}{All}{$bait} = ( $countsByBait{$_}{WtLeft}{$bait} +
                                          $countsByBait{$_}{WtRight}{$bait} +
                                          $countsByBait{$_}{MutWtLeft}{$bait} +
                                          $countsByBait{$_}{MutWtRight}{$bait} ) / 2;
        $umisByBait{$_}{All}{$bait} = ( $umisByBait{$_}{WtLeft}{$bait} +
                                        $umisByBait{$_}{WtRight}{$bait} +
                                        $umisByBait{$_}{MutWtLeft}{$bait} +
                                        $umisByBait{$_}{MutWtRight}{$bait} ) / 2;
      }
    }
  }
}

# Adjust bait depths by other mutant reads for unambig and nobind reads
# Would be good to do this in future for properly assessed mutA_mutB_wt reads
# Also would be good to properly incorporate adjusted counts
foreach ( keys %itdkeysfiltered ) {
  foreach my $bait ( @baitsAll ) { 
    if( exists( $baitIssues{$_}{unambig}{$bait} ) ||
        exists( $baitIssues{$_}{nobind}{$bait} ) ) {
      foreach my $itd ( keys %itdkeysfiltered ) {
        if( $itd ne $_ ) { 
	  $countsByBait{$_}{All}{$bait} += $countsByBait{$itd}{Mut}{$bait};
	  $umisByBait{$_}{All}{$bait} += $umisByBait{$itd}{Mut}{$bait};
        }  
      }
    } else {
      foreach my $itd ( keys %itdkeysfiltered ) {
        if( $itd ne $_ ) { 
	  $countsByBait{$_}{All}{$bait} += ( $countsByBait{$_}{MutMutLeft}{$bait} +
					     $countsByBait{$_}{MutMutRight}{$bait} ) / 2; 
	  $umisByBait{$_}{All}{$bait} += ( $umisByBait{$_}{MutMutLeft}{$bait} +
					   $umisByBait{$_}{MutMutRight}{$bait} ) / 2; 
        }
      }
    }
  }
}

# Make sure countsByBait of Mut is not greater than countsByBait of All
foreach ( keys %itdkeysfiltered ) {
  foreach my $bait ( @baitsAll ) { 
    if( $countsByBait{$_}{All}{$bait} < $countsByBait{$_}{Mut}{$bait} ) {
      $countsByBait{$_}{All}{$bait} = $countsByBait{$_}{Mut}{$bait};
    }
    if( $umisByBait{$_}{All}{$bait} < $umisByBait{$_}{Mut}{$bait} ) {
      $umisByBait{$_}{All}{$bait} = $umisByBait{$_}{Mut}{$bait};
    }
    if( $countsByBaitAdj{$_}{All}{$bait} + $countsByBait{$_}{All}{$bait} < 
	$countsByBaitAdj{$_}{Mut}{$bait} + $countsByBait{$_}{Mut}{$bait} ) {
      $countsByBaitAdj{$_}{All}{$bait} = $countsByBaitAdj{$_}{Mut}{$bait} +
	  $countsByBait{$_}{Mut}{$bait} - $countsByBait{$_}{All}{$bait};
    }
    if( $umisByBaitAdj{$_}{All}{$bait} + $umisByBait{$_}{All}{$bait} < 
	$umisByBaitAdj{$_}{Mut}{$bait} + $umisByBait{$_}{Mut}{$bait} ) {
      $umisByBaitAdj{$_}{All}{$bait} = $umisByBaitAdj{$_}{Mut}{$bait} +
	  $umisByBait{$_}{Mut}{$bait} - $umisByBait{$_}{All}{$bait};
    }
  }
}
}

foreach( keys %itdkeysfiltered ) {
  $coverageMut{$_} = 0; $coverageMutLeft{$_} = 0; $coverageMutRight{$_} = 0;
  $coverageMutWt{$_} = 0; $coverageMutWtLeft{$_} = 0; $coverageMutWtRight{$_} = 0;
  $coverageWt{$_} = 0; $coverageWtLeft{$_} = 0; $coverageWtRight{$_} = 0;
  $covUmisMut{$_} = 0; $covUmisMutLeft{$_} = 0; $covUmisMutRight{$_} = 0;
  $covUmisMutWtLeft{$_} = 0; $covUmisMutWtRight{$_} = 0;
  $covUmisWt{$_} = 0; $covUmisWtLeft{$_} = 0; $covUmisWtRight{$_} = 0;
  
  if(exists($readsMut{$_})) { $coverageMut{$_} = scalar( keys %{$readsMut{$_}} ); }
  if(exists($readsMutLeft{$_})) { $coverageMutLeft{$_} = scalar( keys %{$readsMutLeft{$_}} ); }
  if(exists($readsMutRight{$_})) { $coverageMutRight{$_} = scalar( keys %{$readsMutRight{$_}} ); }
  if(exists($readsMutWt{$_})) { $coverageMutWt{$_} = scalar( keys %{$readsMutWt{$_}} ); }
  if(exists($readsMutWtLeft{$_})) { $coverageMutWtLeft{$_} = scalar( keys %{$readsMutWtLeft{$_}} ); }
  if(exists($readsMutWtRight{$_})) { $coverageMutWtRight{$_} = scalar( keys %{$readsMutWtRight{$_}} ); }
  if(exists($readsWt{$_})) { $coverageWt{$_} = scalar( keys %{$readsWt{$_}} ); }
  if(exists($readsWtLeft{$_})) { $coverageWtLeft{$_} = scalar( keys %{$readsWtLeft{$_}} ); }
  if(exists($readsWtRight{$_})) { $coverageWtRight{$_} = scalar( keys %{$readsWtRight{$_}} ); }
  if(exists($umisMutAll{$_})) { $covUmisMut{$_} = scalar( keys %{$umisMutAll{$_}} ); }
  if(exists($umisMutLeft{$_})) { $covUmisMutLeft{$_} = scalar( keys %{$umisMutLeft{$_}} ); }
  if(exists($umisMutRight{$_})) { $covUmisMutRight{$_} = scalar( keys %{$umisMutRight{$_}} ); }
  if(exists($umisMutWtLeft{$_})) { $covUmisMutWtLeft{$_} = scalar( keys %{$umisMutWtLeft{$_}} ); }
  if(exists($umisMutWtRight{$_})) { $covUmisMutWtRight{$_} = scalar( keys %{$umisMutWtRight{$_}} ); }
  if(exists($umisWtAll{$_})) { $covUmisWt{$_} = scalar( keys %{$umisWtAll{$_}} ); }
  if(exists($umisWtLeft{$_})) { $covUmisWtLeft{$_} = scalar( keys %{$umisWtLeft{$_}} ); }
  if(exists($umisWtRight{$_})) { $covUmisWtRight{$_} = scalar( keys %{$umisWtRight{$_}} ); }
}

my %coverageRadius = ();
foreach my $cname (keys %itdkeysfiltered) {
    @line_cells = split( /\_/, $cname );
    $coverageRadius{$cname} = [(0) x (substr($line_cells[1], 6)+length($refseq))];
}

open( FI, $fbase . "_to_candidate_ITDs.filtered.sam" ) or die $!;
open( FO, ">" . $fbase . "_to_candidate_ITDs.mutreads.sam" ) or die $!;
while(<FI>) {
    if( /^@/ ) {
	print FO $_; 
    } else {
	@line_cells = split( /\t/, $_ );
	if( exists($readsMut{$line_cells[2]}{$line_cells[0]}) ) {
	    print FO $_;
	    for( my $i=$line_cells[3]; $i < $line_cells[3]+length($line_cells[9]); $i++ ) {
		$coverageRadius{$line_cells[2]}[$i-1] = max($coverageRadius{$line_cells[2]}[$i-1], 
							  min($i-$line_cells[3]+1, $line_cells[3]+length($line_cells[9])-$i) );
	    }
	}
    }
}
close FO;
close FI;

if( !$debug && system( "rm " . $fbase . "_to_candidate_ITDs.filtered.sam" ) ) {
  die "Error removing to_candidate_ITDs_filtered alignment file. Exiting...";
}

my $rAFsum = 0; my $uAFsum = 0; my $auAFsum = 0;
my %coverageDets = (); my %coverageTots = ();
foreach my $cname ( keys %itdkeysfiltered ) {
    $coverageDets{$cname} = { "0"=>0, "1-5"=>0, "6-10"=>0, "11-25"=>0, "26-50"=>0, ">50"=>0 };
    for( my $i = $candidatedets{$cname}{dup_r0};
	 $i <= 2 * $candidatedets{$cname}{dup_r1} - $candidatedets{$cname}{dup_r0} + $candidatedets{$cname}{inslen} + 1; $i++ ) {
	if( $coverageRadius{$cname}[$i] == 0 ) {
	    $coverageDets{$cname}{"0"} += 1; 
	} elsif( $coverageRadius{$cname}[$i] <= 5 ) {
	    $coverageDets{$cname}{"1-5"} += 1; 
	} elsif( $coverageRadius{$cname}[$i] <= 10 ) {
	    $coverageDets{$cname}{"6-10"} += 1; 
	} elsif( $coverageRadius{$cname}[$i] <= 25 ) {
	    $coverageDets{$cname}{"11-25"} += 1; 
	} elsif( $coverageRadius{$cname}[$i] <= 50 ) {
	    $coverageDets{$cname}{"26-50"} += 1; 
	} else {
	    $coverageDets{$cname}{">50"} += 1; 
	}
    }
    my $mutCount = 0; my $mwCount = 0; my $mutCountUmi = 0; my $mwCountUmi = 0;
    my $mutCountAdj = 0; my $mutCountUmiAdj = 0; my $mwCountAdj = 0; my $mwCountUmiAdj = 0;
    if( $probes ne "" ) {
	foreach my $bait ( @baitsAll ) {
	    $mutCount += $countsByBait{$cname}{Mut}{$bait};
	    $mwCount += $countsByBait{$cname}{All}{$bait};
	    $mutCountUmi += $umisByBait{$cname}{Mut}{$bait};
	    $mwCountUmi += $umisByBait{$cname}{All}{$bait};
	    $mutCountAdj += $countsByBaitAdj{$cname}{Mut}{$bait};
	    $mwCountAdj += $countsByBaitAdj{$cname}{All}{$bait};
	    $mutCountUmiAdj += $umisByBaitAdj{$cname}{Mut}{$bait};
	    $mwCountUmiAdj += $umisByBaitAdj{$cname}{All}{$bait};
	}
	$mutCountAdj += $mutCount;
	$mwCountAdj += $mwCount;
	$mutCountUmiAdj += $mutCountUmi;
	$mwCountUmiAdj += $mwCountUmi;
	$coverageTots{$cname}{rAdj} = $mutCountAdj;
	$coverageTots{$cname}{uAdj} = $mutCountUmiAdj;
	$coverageTots{$cname}{arAdj} = $mutCountAdj/$mwCountAdj;
    } elsif( $ngstype eq "HC" ) {
	$mutCount = ($coverageMutLeft{$cname}+$coverageMutRight{$cname})/2;
	$mwCount = ($coverageWtLeft{$cname}+$coverageWtRight{$cname}+
		    $coverageMutWtLeft{$cname}+$coverageMutWtRight{$cname})/2;
	$mutCountUmi = ($covUmisMutLeft{$cname}+$covUmisMutRight{$cname})/2;
	$mwCountUmi = ($covUmisWtLeft{$cname}+$covUmisWtRight{$cname}+
		       $covUmisMutWtLeft{$cname}+$covUmisMutWtRight{$cname})/2;
    } elsif( $ngstype eq "amplicon" ) {
	$mutCount = $coverageMut{$cname};
	$mwCount = $coverageWt{$cname}+$coverageMut{$cname};
	$mutCountUmi = $covUmisMut{$cname};
	$mwCountUmi = $covUmisWt{$cname}+$covUmisMut{$cname};
    }

    $mwCount = max( $mutCount, $mwCount );
    $coverageTots{$cname}{rMut} = $mutCount;
    $coverageTots{$cname}{rAll}= $mwCount;
    if( $mwCount == 0 ) { $coverageTots{$cname}{rAF} = 0; } else { $coverageTots{$cname}{rAF} = $mutCount/$mwCount; }
    
    if( $umitag ne "" ) {
	$mwCountUmi = max( $mutCountUmi, $mwCountUmi );
	$mwCountUmiAdj = max( $mutCountUmiAdj, $mwCountUmiAdj );
	$coverageTots{$cname}{uMut} = $mutCountUmi;
	$coverageTots{$cname}{uAll} = $mwCountUmi;
	if( $mwCountUmi == 0 ) { $coverageTots{$cname}{uAF} = 0; } else { $coverageTots{$cname}{uAF} = $mutCountUmi/$mwCountUmi; }
	if( $mwCountUmiAdj == 0 ) { $coverageTots{$cname}{auAF} = 0; } else { $coverageTots{$cname}{auAF} = $mutCountUmiAdj/$mwCountUmiAdj; }
    } else {
	$coverageTots{$cname}{uMut} = 0;
	$coverageTots{$cname}{uAll} = 0;
	$coverageTots{$cname}{uAF} = 0;
	$coverageTots{$cname}{auAF} = 0;
    }
    $rAFsum += $coverageTots{$cname}{rAF};
    $uAFsum += $coverageTots{$cname}{uAF};
    $auAFsum += $coverageTots{$cname}{auAF};
}

# May need to re-consider edge cases where $rAFsum, $uAFsum, or $auAFsum is greater than 1 (e.g. alternatively scale down)
foreach my $cname (keys %itdkeysfiltered) {
  if( $rAFsum < 1 ) {
      $coverageTots{$cname}{rAR} = $coverageTots{$cname}{rAF}/(1-$rAFsum);
  } elsif( $coverageTots{$cname}{rAF} < 1 ) {
      $coverageTots{$cname}{rAR} = $coverageTots{$cname}{rAF}/(1-$coverageTots{$cname}{rAF});
  } else {
      $coverageTots{$cname}{rAR} = "Inf";
  }
  if( $uAFsum < 1 ) {
      $coverageTots{$cname}{uAR} = $coverageTots{$cname}{uAF}/(1-$uAFsum);
  } elsif( $coverageTots{$cname}{uAF} < 1 ) {
      $coverageTots{$cname}{uAR} = $coverageTots{$cname}{uAF}/(1-$coverageTots{$cname}{uAF});
  } else {
      $coverageTots{$cname}{uAR} = "Inf";
  }
  if( $auAFsum < 1 ) {
      $coverageTots{$cname}{auAR} = $coverageTots{$cname}{auAF}/(1-$auAFsum);
  } elsif( $coverageTots{$cname}{auAF} < 1 ) {
      $coverageTots{$cname}{auAR} = $coverageTots{$cname}{auAF}/(1-$coverageTots{$cname}{auAF}); 
  } else {
      $coverageTots{$cname}{auAR} = "Inf";
  }
}

# Remove itds with no mutant reads and flag ITDs where neither endpoints is within an exon + buffer
my %itdother = ();
foreach( keys %itdkeysfiltered ) {
    if( $coverageMut{$_} == 0 || $coverageTots{$_}{rMut} == 0 ) {
	delete $itdkeysfiltered{$_};
    } elsif( !($candidatedets{$_}{dup_r0} > 585 && $candidatedets{$_}{dup_r0} < 721 ) &&
	     !($candidatedets{$_}{dup_r1} > 585 && $candidatedets{$_}{dup_r1} < 721 ) &&
	     !($candidatedets{$_}{dup_r0} > 808 && $candidatedets{$_}{dup_r0} < 916 ) &&
	     !($candidatedets{$_}{dup_r1} > 808 && $candidatedets{$_}{dup_r1} < 916 ) ) {
	$itdother{$_} = 1;
    }
}

my $outvcf = ""; my $outsummary = ""; my $outother = ""; my $afb; my $uafb; my $tmpaf;
foreach( keys %itdkeysfiltered ) {
    my $nobind = "None"; my $dupbait = "None"; my $unambig = "None"; my @baits = ();
    my $ucovb; my $rcovb; my $udetb;
    
    my $cdot = $candidatedets{$_}{keycdot}; my $pdot = $candidatedets{$_}{keypdot};
    if( $pdot eq "not_exonic" && $candidatedets{$_}{keypdotalt} ne "not_exonic" ) {
	$pdot = "not_exonic(alt_ " . $candidatedets{$_}{keypdotalt} . ")";
    }
    
    if( $probes ne "" ) {
      if( exists( $baitIssues{$_}{nobind} ) ) {
        @baits = keys %{$baitIssues{$_}{nobind}};
	if( scalar(@baits) > 0 ) { $nobind = $baits[0]; }
	if( scalar(@baits) > 1 ) {
          for( $i=1; $i<scalar(@baits); $i++ ) { $nobind .= "|" . $baits[$i]; }
	} 
      }
      if( exists( $baitIssues{$_}{dupbait} ) ) {
        @baits = keys %{$baitIssues{$_}{dupbait}};
	if( scalar(@baits) > 0 ) { $dupbait = $baits[0]; }
	if( scalar(@baits) > 1 ) {
          for( $i=1; $i<scalar(@baits); $i++ ) { $dupbait .= "|" . $baits[$i]; }
	} 
      }
      if( exists( $baitIssues{$_}{unambig} ) ) {
        @baits = keys %{$baitIssues{$_}{unambig}};
	if( scalar(@baits) > 0 ) { $unambig = $baits[0]; }
	if( scalar(@baits) > 1 ) {
          for( $i=1; $i<scalar(@baits); $i++ ) { $unambig .= "|" . $baits[$i]; }
	} 
      }
	
      my $baitdist = "";
      if(exists($baitDists{$_}{$baitsAll[0]})) { $baitdist=$baitDists{$_}{$baitsAll[0]}; }

      my $tb = $countsByBait{$_}{All}{$baitsAll[0]};
      my $vb = $countsByBait{$_}{Mut}{$baitsAll[0]};
      my $utb = $umisByBait{$_}{All}{$baitsAll[0]};
      my $uvb = $umisByBait{$_}{Mut}{$baitsAll[0]};
      if( $tb > 0 ) { $afb = sprintf( "%.3g", $vb / $tb ); } else { $afb = "NA"; }
      if( $utb > 0 ) { $uafb = sprintf( "%.3g", $uvb / $utb ); } else { $uafb = "NA"; }
      $ucovb = $uafb . "(" . $uvb . "/" . $utb . ")[" . $baitdist . "]";
      $rcovb = $afb . "(" . $vb . "/" . $tb . ")";
      $udetb = $uvb . "/" . $utb;

      for( $i=1; $i<scalar(@baitsAll); $i++ ) {
        $baitdist = "";
	if(exists($baitDists{$_}{$baitsAll[$i]})) { $baitdist=$baitDists{$_}{$baitsAll[$i]}; }

	my $tbi = $countsByBait{$_}{All}{$baitsAll[$i]};
	my $vbi = $countsByBait{$_}{Mut}{$baitsAll[$i]};
	my $utbi = $umisByBait{$_}{All}{$baitsAll[$i]};
	my $uvbi = $umisByBait{$_}{Mut}{$baitsAll[$i]};

	if( $tbi > 0 ) { $tmpaf = sprintf( "%.3g", $vbi/$tbi ); } else { $tmpaf = "NA"; }
	$afb .= "," . $tmpaf;
	$rcovb .= ";" . $tmpaf . "(" . $vbi . "/" . $tbi . ")";

	if( $utbi > 0 ) { $tmpaf = sprintf( "%.3g", $uvbi/$utbi ); } else { $tmpaf = "NA"; }
	$uafb .= "," . $tmpaf;
	$ucovb .= ";" . $tmpaf . "(" . $uvbi . "/" . $utbi . ")[" . $baitdist . "]";
	$udetb .= "," . $uvbi . "/" . $utbi;
      }
    }
    
    my $vcfstring = proc_vcf_from_keyfull($itdkeys{$_}, $cdot ) . "\t" .
	"GENE=FLT3;STRAND=-;SVLEN=" . $candidatedets{$_}{itdsize} . 
	";CDS=" . $cdot . ";AA=" . $pdot .
        ";AR=" . sprintf( "%.4g", $coverageTots{$_}{uAR} ) . 
	";AF=" . sprintf( "%.4g", $coverageTots{$_}{uAF} ) . 
        ";DP=" . $coverageTots{$_}{uAll} . ";VD=" . $coverageTots{$_}{uMut} . 
 	";AAR=" . sprintf( "%.4g", $coverageTots{$_}{auAR} ) . 
 	";AAF=" . sprintf( "%.4g", $coverageTots{$_}{auAF} ) . 
        ";RAR=" . sprintf( "%.4g", $coverageTots{$_}{rAR} ) . 
	";RAF=" . sprintf( "%.4g", $coverageTots{$_}{rAF} ) . 
        ";RDP=" . $coverageTots{$_}{rAll} . ";RVD=" . $coverageTots{$_}{rMut} .
	";SAMPLE=" . $shortkey .
	"\tAR:AF";
    if( $probes ne "" ) { $vcfstring .= ":AFB"; }
    $vcfstring .= ":DP:VD:AAR:AAF:RAR:RAF:RDP:RVD:CR";
    if( $probes ne "" ) { $vcfstring .= ":DB:NB:RB:AFB1:Fwd1:Fwd2:Fwd3:Rev1:Rev2:Rev3:Oth"; }
    $vcfstring .= "\t" . 
	sprintf( "%.4g", $coverageTots{$_}{uAR} ) . ":" . 
	sprintf( "%.4g", $coverageTots{$_}{uAF} );
    if( $probes ne "" ) { $vcfstring .= ":" . $uafb; }
    $vcfstring .=  ":" .
	$coverageTots{$_}{uAll} . ":" . $coverageTots{$_}{uMut} . ":" . 
	sprintf( "%.4g", $coverageTots{$_}{auAR} ) . ":" . 
	sprintf( "%.4g", $coverageTots{$_}{auAF} ) . ":" . 
	sprintf( "%.4g", $coverageTots{$_}{rAR} ) . ":" .
	sprintf( "%.4g", $coverageTots{$_}{rAF} ) . ":" . 
	$coverageTots{$_}{rAll} . ":" . $coverageTots{$_}{rMut} . ":" . 
	$coverageDets{$_}{"0"} . "," . $coverageDets{$_}{"1-5"} . "," .
	$coverageDets{$_}{"6-10"} . "," . $coverageDets{$_}{"11-25"} . "," .
	$coverageDets{$_}{"26-50"} . "," . $coverageDets{$_}{">50"};
    if( $probes ne "" ) { $vcfstring .= ":" . $dupbait . ":" . $nobind . ":". $unambig . ":" . $ucovb; }
    $vcfstring .= "\n"; 

    my $tabsummary = $shortkey . "\t" . $candidatedets{$_}{itdsize} . "\t" . $cdot . "\t" . $pdot . "\t";
    if( $probes ne "" ) { $tabsummary .= $dupbait . "\t" . $nobind . "\t". $unambig . "\t"; }
    $tabsummary .= $coverageDets{$_}{"0"} . "," . $coverageDets{$_}{"1-5"} . "," .
	$coverageDets{$_}{"6-10"} . "," . $coverageDets{$_}{"11-25"} . "," .
	$coverageDets{$_}{"26-50"} . "," . $coverageDets{$_}{">50"} . "\t" .
 	sprintf( "%.4g", $coverageTots{$_}{auAR} ) . "\t" . 
 	sprintf( "%.4g", $coverageTots{$_}{auAF} ) . "\t" .
        sprintf( "%.4g", $coverageTots{$_}{uAR} ) . "\t" .
	sprintf( "%.4g", $coverageTots{$_}{uAF} ) . "\t";
    if( $probes ne "" ) { $tabsummary .= $uafb . "\t"; }
    $tabsummary .= $coverageTots{$_}{uMut} . "\t" . $coverageTots{$_}{uAll} . "\t" .
        sprintf( "%.4g", $coverageTots{$_}{rAR} ) . "\t" .
	sprintf( "%.4g", $coverageTots{$_}{rAF} ) . "\t" .
        $coverageTots{$_}{rMut} . "\t" . $coverageTots{$_}{rAll};
    if( $probes ne "" ) { $tabsummary .= "\t" . $ucovb . "\t". $rcovb; }
    $tabsummary .= "\n";

    my $covLength = sum values %{$coverageDets{$_}};
    if( exists( $itdother{$_} ) ||
	( $umitag ne "" && $coverageTots{$_}{uMut} < $minreads ) ||
	( $umitag eq "" && $coverageTots{$_}{rMut} < $minreads ) ||
	( $candidatedets{$_}{itdsize} >= $noisyMinItdSize && $covLength > 0 && 
	  $coverageDets{$_}{"0"}/$covLength >= $noisyMinNocovRatio ) ) {
	$outother .= $tabsummary;
    } else {
	$outvcf .= $vcfstring;
	$outsummary .= $tabsummary;
    }
}

foreach( keys %itdkeysfiltered ) {
    if( $candidatedets{$_}{dup_r0} == -1 && $candidatedets{$_}{dup_r1} == -1 &&
        $candidatedets{$_}{mutbk_s0} == -1 && $candidatedets{$_}{mutbk_s1} == -1 ) {
	$outother .= $shortkey . "\t" . $candidatedets{$_}{itdsize} . "\t" .
	    $candidatedets{$_}{keycdot} . "\t" . $candidatedets{$_}{keypdot} . "\t";
	if( $probes ne "" ) { $outother .= "?\t?\t?\t?\t?\t?\t"; }
	$outother .= "?\t?\t?\t?\t?\t?\t?\t?\t?\t?\t?\n";
    }
}

if( $outvcf ne "" ) {
  open( FO, ">" . $fbase . "_ITD.vcf" ) or die $!;
  print FO $vcfheader . $shortkey . "\n";
  print FO $outvcf;
  close FO;
} else {
  system( "touch " . $fbase . "_ITD_none.vcf" );
}

if( $outsummary ne "" ) {
  open( FO, ">" . $fbase . "_ITD_summary.txt" ) or die $!;
  print FO $outsummary;
  close FO;
}
  
if( $outother ne "" ) {
  open( FO, ">" . $fbase . "_other_summary.txt" ) or die $!;
  print FO $outother;
  close FO;
}

#### GENERATE TVIEW HTML FILES
if( $web ) {
  foreach( keys %itdkeysfiltered ) {
    $candidate = $_;
    $candidateseq = $candidateSeqs{$_}; 
    $candidatename = $candidateNames{$_}; 
    open (FO, ">" . $fbase . "_candidate" . $candidatename . ".fa" ) or die $!;
    print FO ">" . $candidate . "\n" . $candidateseq . "\n";
    close FO;
  }

  if( system( sprintf( "%s view -bS %s > %s", $samcmd,
		 $fbase . "_to_candidate_ITDs.mutreads.sam",
		 $fbase . "_to_candidate_ITDs.mutreads.bam" ) ) ||
    system( sprintf( "%s sort %s -o %s", $samcmd,
		 $fbase . "_to_candidate_ITDs.mutreads.bam",
		 $fbase . "_to_candidate_ITDs.mutreads.sorted.bam" ) ) ||
    system( sprintf( "%s index %s", $samcmd, $fbase . "_to_candidate_ITDs.mutreads.sorted.bam" ) ) ) {
    die "Error making sorted candidate_ITDs mutreads alignment file. Exiting...";
  }

  foreach( keys %itdkeysfiltered ) {
    my $fi = $fbase . "_other_candidate" . $_ . "_reads=" . $coverageMut{$_} . ".html";
    @line_cells = split( /\_/, $_ );
    if( system( sprintf( "COLUMNS=%s %s tview -d H %s %s > %s",
		     $extsizeHash{$_},
		     $samcmd,
		     $fbase . "_to_candidate_ITDs.mutreads.sorted.bam",
		     $fbase . "_candidate" . $_ . ".fa",
		     $fi ) ) ) {
	die "Error executing samtools tview. Exiting...";
    }

    if( !$debug && (
	system( "rm " . $fbase . "_candidate" . $_ . ".fa" ) ||
	system( "rm " . $fbase . "_candidate" . $_ . ".fa.fai" ) ) ) {
	die "Error removing candidate fasta files. Exiting...";
    }

    my $replaceseq = ""; my $replacekey = ""; my $ckey = ""; my $fo;
    if( exists( $candidatedets{$_} ) ) {
	$replaceseq = $candidatedets{$_}{dupseq};
	$replacekey = $_;
	#$ckey = $candidatedets{$_}{keycdot}; $ckey =~ s/c\./c/; $ckey =~ s/\//ins/g;
	#$fo = $fbase . "_" . $ckey . "_length=" . $candidatedets{$_}{itdsize} . "_reads=" . $coverageMut{$_} . ".html";
	$ckey = $candidatedets{$_}{keycdot} . "_length=" . $candidatedets{$_}{itdsize} . "_rawcov=" . $coverageTots{$_}{rMut}; 
	$fo = $fbase . "_candidate_" . $_ . ".html";
	print $fo . "\n";
    } 

    if( length( $replaceseq ) > 0 ) {
	my $htmlold = $replaceseq;
	my $htmlnew = "<font color=\"fuchsia\">" . $replaceseq . "</font>";
	open (FI, $fi) or die $!;
	open (FO, ">$fo" ) or die $!;
	while(<FI>) {
	    $_ =~ s/$htmlold/$htmlnew/g;
	    $_ =~ s/"fuchsia"/"blue"/;
	    $_ =~ s/$replacekey/$ckey/g;
	    print FO $_;
	}
	close FI;
	close FO;

	if( system( sprintf( "rm %s", $fi ) ) ) { die "Error removing htmlold. Exiting..."; }
    }
  }

  if( !$debug && (
    system( "rm " . $fbase . "_to_candidate_ITDs.mutreads.bam" ) ||
    system( "rm " . $fbase . "_to_candidate_ITDs.mutreads.sorted.bam" ) ||
    system( "rm " . $fbase . "_to_candidate_ITDs.mutreads.sorted.bam.bai" ) ) ) {
    die "Error removing to_candidte_ITDs mutreads alignment files. Exiting..." ;
  }
}

if( !$debug && system( "rm " . $fbase . "_to_candidate_ITDs.mutreads.sam" ) ) {
    die "Error removing to_candidte_ITDs mutreads alignment files. Exiting..." ;
}
