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

FLT3_ITD_ext.pl  Process bamfiles or paired fastq files for FLT3-ITDs

=head1 SYNOPSIS

  --bam, -b	Input bamfile (either this or fastq1+2 required)
  --typeb, -t	Reads to extract from input bam (defaults to "targeted" for FLT3 alignments; or can be "all")
  --fastq1, -f1	Input fastq1 (either fastq1+2 or bam required)
  --fastq2, -f2	Input fastq2 (either fastq1+2 or bam required)
  --output, -o  Output path (required)
  --arch, -a	Archive small fastq files (defaults to false)
  --debug, -d	Debugging files (defaults to false)
  --web, -w	Create html webpages for each ITD call (defaults to true)
  --help, -h	Print this help

=head1 VERSION

0.01

=cut

GetOptions(
  "bam|b=s" => \(my $inbam = "" ),
  "typeb|t=s" => \ (my $typeb = "targeted"),
  "fastq1|f1=s" => \(my $fastq1 = ""),
  "fastq2|f2=s" => \(my $fastq2 = ""),
  "output|o=s" => \ my $outpath,
  "arch|a" => \(my $archive = 0 ),
  "debug|d" => \(my $debug = 0),
  "web|w" => \(my $web = 1),
  "help|h" => sub { HelpMessage(0) }, 
) or HelpMessage(1);

HelpMessage(1) unless( $outpath && (length($inbam)>0 || (length($fastq1)>0 && length($fastq2)>0)) );

my $samcmd = "/apps/samtools/1.9/bin/samtools";
my $bedcmd = "/apps/bedtools2/2.23.0/bamToFastq";
my $bwacmd = "/apps/bwa/0.7.17/bwa";
my $clustercmd = "/home/ht50/tools/sumaclust_v1.0.31/sumaclust";

my $minToCluster = 1;
my $maxEditDist = 5;  #$maxEditDist = 4;
my $maxInsAllowed = 2; # cannot have more than this in total ins to count as wt or itd
my $maxDelAllowed = 2; # cannot have more than this in total del to count as wt or itd
my $maxClipAllowed = 2; # cannot have more than this in max softclip to count as wt or itd
my $clipFilterRatio = 0.6;  # ratio of total clipping length to sequence length needs to be greater than this value to pass filter
my $mutFilterRatio = 0.05;  # ratio of total number of mismatches + deletions to matched sequence length must be greater than this value to pass filter
my $minClipToAlign = 9;
my $manualItdCheckLength = 12;
my $buffer = 10;  # buffer size past breakpoint needed for alignment to be included in counts

my $clipMatchRatioThreshold = 0.5; # ratio of clip matches to clip length needs to be at least this value in order to extend the clip

# May be able to merge these, however staggered approach generates less false positives
my $bwaSeedLength = 15;
my $bwaThreshold = 12;
my $bwaClipsSeedLength = 6;
my $bwaClipsThreshold = 9;

#### reference FLT3 sequence (1500bp region centered around e14e15): chr13_28607438_28608937_minus_strand
my $basedir = "/home/ht50/FLT3/";
my $gene = "FLT3"; my $transcript = "NM_004119"; 
my $refkey = "FLT3_dna_e14e15";
my $refindex = $basedir . $refkey;
my $refstart = 28607438; my $refend = 28608937;  # 1500bp centered around exons 14-15 (with 586bp buffer)
my $target_region = "13:" . $refstart .  "-" . $refend;
my $reffasta = $basedir . $refkey . ".fa";
my $refseq = ""; 
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
my $proc_alignerLocal = 1;
my $proc_output_rhp = 0;

my $maxITDsize = 500;
my $maxInsFrac = 0.9; #0.5;  # do not allow ectopic inserts more than this fractional length of the ITD

my @line_cells; my @sub_cells; my @p_cells; my @s_cells; my @a;
my $key; my $readname; my $cigar; my $seq; my $qa; my $i; my $j; my $k;

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
    "##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">\n" .
    "##FORMAT=<ID=DP,Number=1,Type=Float,Description=\"Total Depth\">\n" .
    "##FORMAT=<ID=VD,Number=1,Type=Float,Description=\"Variant Depth\">\n" .
    "##INFO=<ID=GENE,Number=1,Type=String,Description=\"Gene name\">\n" .
    "##INFO=<ID=STRAND,Number=1,Type=String,Description=\"Gene strand\">\n" .
    "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n" .
    "##INFO=<ID=CDS,Number=1,Type=String,Description=\"CDS annotation (modified)\">\n" .
    "##INFO=<ID=AA,Number=1,Type=String,Description=\"Peptide annotation (modified)\">\n" .
#    "##INFO=<ID=ALTCDS,Number=1,Type=String,Description=\"Alternative CDS annotation (modified)\">\n" .
#    "##INFO=<ID=ALTAA,Number=1,Type=String,Description=\"Alternative peptide annotation (modified)\">\n" .
    "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">\n" .
    "##INFO=<ID=DP,Number=1,Type=Float,Description=\"Total Depth\">\n" .
    "##INFO=<ID=VD,Number=1,Type=Float,Description=\"Variant Depth\">\n" .
    "##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Sample name\">\n" .
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
	
	$cdet{dup_r0} = $rightr0;
	$cdet{dup_r1} = $leftr1;
	$cdet{mutbk_s0} = $leftr1;
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
	for( my $i=1; $i<scalar(@cells)-1; $i++ ) {
	    if( length($insstr)>0 ) { $insstr .= "/"; }
	    if( substr($cells[$i],0,3) eq "ins" ) {
		$insstr .= substr($cells[$i],3);
	    } elsif( substr($cells[$i],0,3) eq "ref" ) {
		@subcells = split(/\(/, substr($cells[$i],3) );
		@subsubcells = split(/\-/, $subcells[0]);
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

	    #"ins" . substr( $candidateseq, $leftr1-$leftr0+1, $cdet{itdsize} );
	    #"ins" . substr( $candidateseq, $leftr1-$leftr0, $cdet{itdsize} );
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

    my @cells = split(/\_/, $keyf);

    my $ins = ""; my $ref0 = ""; my $ref1 = ""; my $mm0 = ""; my $mm1 = "";

    $ref0 = $cells[0];
    if( scalar(@cells) == 2 ) {
	$ref1 = $cells[1];
    } elsif( scalar(@cells) == 3 ) {
	$ins = $cells[1];
	$ref1 = $cells[2];
    }

    if( $ref0 !~ /^ref/ || $ref1 !~ /^ref/ ||
	(length($ins)>0 && $ins !~ /^ins/ ) ) {
	return( "" );
    }

    $ref0 =~ s/ref//; $ref1 =~ s/ref//; $ins =~ s/ins//;

    @cells = split( /\(/, $ref0 );
    my @coords0 = split( /\-/, $cells[0] );
    if( scalar(@cells) > 1 ) { $mm0 = $cells[1]; $mm0 =~ s/\)//; }

    @cells = split( /\(/, $ref1 );
    my @coords1 = split( /\-/, $cells[0] );
    if( scalar(@cells) > 1 ) { $mm1 = $cells[1]; $mm1 =~ s/\)//; }

    my @mmlocs0 = (); my @mmlocs1 = (); my $seq0 = $refseq; my $seq1 = $refseq;
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

    my $insseq = "";
    if( $inspos+1 >= $dupstart ) {
	$insseq = substr($seq0, $inspos, $coords0[1]-$inspos) . $ins . substr($seq1, $coords1[0]-1, $inspos-$coords1[0]+1);
	#my $refbp = substr( $refseq, $inspos-1, 1 );
	#my $alt = $refbp . $insseq;
	my $refbp = substr( $refseq, $inspos, 1 );
	my $alt = $insseq . $refbp;
	$alt = reverse $alt; $alt =~ tr/ATGCatgc/TACGtacg/; $refbp =~ tr/ATGCatgc/TACGtacg/;
	return("13\t" . ($refend-$inspos) . "\t" . $keyc . "\t" . $refbp . "\t" . $alt . "\t.\t." );
	#return("13\t" . ($refend-$inspos+1) . "\t.\t" . $refbp . "\t" . $alt . "\t.\t." );
    }

    #$insseq = $ins . substr($refseq, $coords1[0]-1, $coords0[1]-$coords1[0]+1);
    return( ".\t.\t.\t.\t.\t.\t.\t" );
}

####################################################################
####################################################################
# BEGIN MAIN SCRIPT
####################################################################
####################################################################

my $fbase; my $outname;
if( length($inbam) > 0 ) {
  @line_cells = split(/\//, $inbam);
  $outname = $line_cells[scalar(@line_cells)-1]; $outname =~ s/\.bam$//; $outname .= "_FLT3";
  $fbase = $outpath . "/" . $outname;
  $fastq1 = $fbase . ".R1.fastq";
  $fastq2 = $fbase . ".R2.fastq";

  if( $typeb eq "targeted" ) {
    # Extract alignments to the FLT3 target locus from the input bamfile
    # This is assumed to include all mutant reads (which should be true if a local aligner was used)
    if( system( sprintf( "%s view -b -F 0x900 %s 13:28607438-28608937 > %s.bam", $samcmd, $inbam, $fbase ) ) ||
        system( sprintf( "%s index %s.bam", $samcmd, $fbase ) ) ||
        system( sprintf( "%s sort -n %s.bam -o %s.sorted.bam", $samcmd, $fbase, $fbase ) ) ) {
      die "Failed to extract existing FLT3 targeted locus alignments from original input bam. Exiting...\n";
    }
  } elsif( $typeb eq "all" ) {
    if( system( sprintf( "%s sort -n %s -o %s.sorted.bam", $samcmd, $inbam, $fbase ) ) ) {
      die "Failed to sort original input bam by readnames. Exiting...\n";
    }
  } else {
      print "ERROR: typeb must either equal \"targeted\" or \"all\"...";
      HelpMessage(1);
  }

  if( system( sprintf( "%s -i %s.sorted.bam -fq %s -fq2 %s", 
		 $bedcmd, $fbase, $fastq1, $fastq2 ) ) ) {
    die "Failed to convert bamfile to fastq file. Exiting...\n";
  }

  if( system( sprintf( "rm %s.bam", $fbase )) ||
      system( sprintf( "rm %s.bam.bai", $fbase )) ||
      system( sprintf( "rm %s.sorted.bam", $fbase )) ) {
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

my $shortkey = $outname; $shortkey =~ s/\.consensus//; $shortkey =~ s/\.filtered\_FLT3//; $shortkey =~ s/\.mapped\_FLT3//;

# clean up old files
if( glob($fbase . "*.html") ) { system( "rm " . $fbase . "*.html" ); }
if( glob($fbase . "*.vcf") ) { system( "rm " . $fbase . "*.vcf" ); }
if( glob($fbase . "*.fastq.gz") ) { system( "rm " . $fbase . "*.fastq.gz" ); }

# Local alignments to reference
if( $proc_index && system( sprintf( "%s index -p %s %s", $bwacmd, $refindex, $reffasta ) ) ) {
  die "Failed to create bwa index for reference.  Exiting...\n" ;
}

if( system( sprintf( "%s mem -k %s -M -O 6 -T %s %s %s %s > %s_%s.sam", 
	$bwacmd, $bwaSeedLength, $bwaThreshold, $refindex, $fastq1, $fastq2, $fbase, $alignerLocal ) ) ) {
  die "Failed to locally align initial reads to targeted FLT3 locus. Exiting...\n";
}

my %endsFwdR2 = (); my %endsRevR2 = ();
my $readsuffix; my $str; # strand
my %psam = (); my %ssam = (); my %pseq = (); my %pqa = (); my %pstr = (); my %spos = ();
my %pcig = (); my %scig = (); my %ploc = (); my %sloc = (); my %sstr = (); my %ppos = ();
my %resolvedWt = (); my %resolvedWtPaired = (); my %pbase = ();
open( FI, $fbase . "_" . $alignerLocal . ".sam" ) or die $!;
while(<FI>){
    if( !/^@/ ) {
	@line_cells = split( /\t/, $_ );

	if( ( $line_cells[1] & 64 ) == 64 ) {
	    $readsuffix = "_1";
	} elsif( ( $line_cells[1] & 128 ) == 128 ) {
	    $readsuffix = "_2";
	} else {
	    next;
	}

	$readname = $line_cells[0] . $readsuffix;
	my $rpos = $line_cells[3]-1;
	$cigar = $line_cells[5];
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
    }
}

my $entryFA = "";
open(FO1, ">" . $fbase . ".nowt.R1.fastq" ) or die $!;
open(FO2, ">" . $fbase . ".nowt.R2.fastq" ) or die $!;
if( $archive ) {
  open(FO3, ">" . $fbase . "_" . $alignerLocal . ".R1.fastq" ) or die $!; 
  open(FO4, ">" . $fbase . "_" . $alignerLocal . ".R2.fastq" ) or die $!; 
}
foreach my $read (keys %pbase ) {
    if( exists($psam{$read . "_1"}) && exists($psam{$read . "_2"}) ) {
	@line_cells = split(/\t/, $psam{$read . "_1"});
	$seq = $line_cells[9]; $qa = $line_cells[10];
       	if( ($line_cells[1] & 16) == 16 ) { $seq = reverse $seq; $seq =~ tr/ATGCatgc/TACGtacg/; $qa = reverse $qa; }
	$entryFA = "@" . $read . "/1\n" . $seq . "\n" . "+" . "\n" . $qa . "\n";
        if( !exists($resolvedWtPaired{$read}) ) { print FO1 $entryFA; }
	if( $archive ) { print FO3 $entryFA; }
	@line_cells = split(/\t/, $psam{$read . "_2"});
	$seq = $line_cells[9]; $qa = $line_cells[10];
       	if( ($line_cells[1] & 16) == 16 ) { $seq = reverse $seq; $seq =~ tr/ATGCatgc/TACGtacg/; $qa = reverse $qa; }	
	$entryFA = "@" . $read . "/2\n" . $seq . "\n" . "+" . "\n" . $qa . "\n";
	if( !exists($resolvedWtPaired{$read}) ) { print FO2 $entryFA; }
	if( $archive ) { print FO4 $entryFA; }
    }
}
close FO1;
close FO2;
if( $archive ) { close FO3; close FO4; }
if( $archive && system( "gzip -f " . $fbase . "_" . $alignerLocal . ".*.fastq" ) ) {
  die "Error gzipping archive fastqs. Exiting...";
}

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

open( FO, ">" . $fbase . "_dets_extseqs.txt" ) or die $! ;
foreach my $read ( keys %lloc ) {
  my $itdsize = length($pseq{$read}) - $rloc{$read} + $lloc{$read};
  if( $itdsize > 0 && $itdsize <= $maxITDsize ) {
    my $extSeq = substr($refseq,0,$lloc{$read}) . $pseq{$read} . substr($refseq,$rloc{$read});
    $extByLength{$itdsize} += 1;
    $extBySeq{$itdsize}{$extSeq} += 1;
    print FO $itdsize . "_" . $lloc{$read} . "_" . $rloc{$read} . "\t" . $read . "\t" . $pcig{$read} . "\t" . $pseq{$read} . "\t" . $extSeq . "\n";
  }
}
close FO;
if( !$debug && system( "rm " . $fbase . "_dets_extseqs.txt" ) ) {
  die "Error removing dets_extseqs file (details on extended reads). Exiting...";
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
    system( "rm " . $fbase . "_sumaclust.fasta" ) ||
    system( "rm " . $fbase . "_for_sumaclust.fasta" ) ||
    system( "rm " . $fbase . "_sumaclust.clusters" ) ||
    system( "rm " . $fbase . "_candidate_clusters.fa" ) ||
    system( "rm " . $fbase . "_candidate_clusters_" . $alignerExts . ".sam" ) ) ) {
    die "Error removing candidate_cluster and sumaclust files. Exiting...";
}

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
    exit(0);
}

open(FO, ">" . $fbase . "_candidate_ITDs.fa" ) or die $!;
foreach( sort keys %itdkeysfiltered ) { print FO ">" . $_ . "\n" . $sumaHash{$_} . "\n"; }
close FO;

# Build bwa candidates index and perform bwa-mem against candidates
if ( system( sprintf( "%s index -p %s_candidate_ITDs %s_candidate_ITDs.fa", $bwacmd, $fbase, $fbase ) ) ||
     system( sprintf( "%s mem -M -O 6 -T 20 %s_candidate_ITDs %s.nowt.R1.fastq %s.nowt.R2.fastq > %s_to_candidate_ITDs_%s.sam", 
		 $bwacmd, $fbase, $fbase, $fbase, $fbase, $alignerToITDs ) ) ) {
  die "Error aligning reads against ITD candidates. Exiting...";
}

if( !$debug && ( system( "rm " . $fbase . ".nowt.R*.fastq" ) || 
    ( $inbam ne "" && system( "rm " . $fbase . ".R*.fastq" ) ) ) ) {
  die "Error removing fastq files. Exiting...";
}

# Filter out suboptimal alignments to candidate clusters
my $header = "";
my %filterFlag;
my %alignmentsFwd = (); my %alignmentsRev = ();
my %regionsFwd = (); my %regionsRev = ();
my %strandsFwd = (); my %strandsRev = ();
open (FI, $fbase . "_to_candidate_ITDs_" . $alignerToITDs . ".sam" ) or die $!;
while(<FI>) {
    if( /^@/ ) {
	$header .= $_;
    } else {
	@line_cells = split( /\t/, $_ );
	$readname = $line_cells[0];
	$filterFlag{$readname} += 0;
	my $str = "+";
	if( ($line_cells[1] & 16) == 16 ) { $str = "-"; }
	if( ($line_cells[1] & 64) == 64 ) {
	    $alignmentsFwd{$readname} = $_;
	    $regionsFwd{$readname} = $line_cells[2];
	    $strandsFwd{$readname} = $str;
	} elsif( ($line_cells[1] & 128) == 128 ) {
	    $alignmentsRev{$readname} = $_;
	    $regionsRev{$readname} = $line_cells[2];
	    $strandsRev{$readname} = $str;
	} else {
	    print( "******** UNEXPECTED ALIGNMENT: NOT A PAIRED READ? ********** \n" );
	}

	$cigar = $line_cells[5];
	my $clipLen = 0;
	my $insLen = 0;
	my $delLen = 0;
	if( @a = $cigar =~ /(\d+)S/g ) { $clipLen = max(@a); }
	if( @a = $cigar =~ /(\d+)I/g ) { foreach my $x (@a) { $insLen += $x; } }
	if( @a = $cigar =~ /(\d+)D/g ) { foreach my $x (@a) { $delLen += $x; } }

	if( $insLen > $maxInsAllowed || $delLen > $maxDelAllowed || $clipLen > $maxClipAllowed ||
	    ($_ =~ /NM:i:(\d+)/ && $1 > $maxEditDist) ) {
	    $filterFlag{$readname} += 1;
	}
    }
}
close FI;

if( !$debug && (
    system( "rm " . $fbase . "_candidate_ITDs.*" ) ||
    system( "rm " . $fbase . "_to_candidate_ITDs_" . $alignerToITDs . ".sam" ) ) ) {
  die "Error removing candidate_ITDs fasta and alignment files. Exiting...";
}

my %candidateReads;
my %candidateReadsAll;
open (FO, ">" . $fbase . "_to_candidate_ITDs.filtered.sam" ) or die $!;
print FO $header;
foreach( keys %filterFlag) {
    if( $filterFlag{$_} == 0 && $regionsFwd{$_} eq $regionsRev{$_} && $strandsFwd{$_} ne $strandsRev{$_} ) {
	print FO $alignmentsFwd{$_};
	print FO $alignmentsRev{$_};
	@line_cells = split(/\t/, $alignmentsFwd{$_});
	$candidateReads{$line_cells[2]}+=1;
	$candidateReadsAll{$line_cells[2]}{$_} = 1;
    }
}
close FO;

my %readsMut = (); my %readsMutLeft = (); my %readsMutRight = ();
my %readsMutWt = (); my %readsMutWtLeft = (); my %readsMutWtRight = ();
my %readsWt = (); my %readsWtLeft = (); my %readsWtRight = ();
my %coverageMut = (); my %coverageMutLeft = (); my %coverageMutRight = (); 
my %coverageMutWt = (); my %coverageMutWtLeft = (); my %coverageMutWtRight = (); 
my %coverageWt = ();  my %coverageWtLeft = ();  my %coverageWtRight = (); 
my %readsPassed = ();

# FUTURE USE
my %readsMutByPrimerR2 = (); my %readsWtByPrimerR2 = ();
my %readsMutByR2 = (); my %readsWtByR2 = (); my %allreadsWtByR2 = (); 
my %locFwd0 = (); my %locFwd1 = (); my %locRev0 = (); my %locRev1 = ();
# END FUTURE USE

open( FI, $fbase . "_to_candidate_ITDs.filtered.sam" ) or die $!;
while(<FI>) {
    if( /^@/ ) { next; }
    @line_cells = split( /\t/, $_ );
    if( $line_cells[2] eq "*" || ( $line_cells[1] & 256 ) == 256 || ( $line_cells[1] & 2048 ) == 2048 ) { next; }

    my $cname = $line_cells[2];
    my $rpos = $line_cells[3]-1;
    $readname = $line_cells[0];
    $cigar = $line_cells[5];
    while( $cigar =~ /(\d+)([SHDIMN])(.*)/ ) {
	$cigar = $3;
	if( $2 eq "M" || $2 eq "N" || $2 eq "D" ) { $rpos += $1; }
    }
    if( !exists($readsMutLeft{$cname}{$readname}) &&
	$candidatedets{$cname}{mutbk_s0} >= $line_cells[3]+$buffer && 
	$candidatedets{$cname}{mutbk_s0} < $rpos-$buffer ) {
	$readsMutLeft{$cname}{$readname} = 1;
	$readsMut{$cname}{$readname} = 1;
	$readsPassed{$readname} = 1;
    }
    if( !exists($readsMutRight{$cname}{$readname}) &&
	$candidatedets{$cname}{mutbk_s1} > $line_cells[3]+$buffer && 
	$candidatedets{$cname}{mutbk_s1} <= $rpos-$buffer ) {
	$readsMutRight{$cname}{$readname} = 1;
	$readsMut{$cname}{$readname} = 1;
	$readsPassed{$readname} = 1;
    }
    if( !exists($readsMutWtLeft{$cname}{$readname}) && 
	$candidatedets{$cname}{dup_r0} > $line_cells[3]+$buffer && 
	$candidatedets{$cname}{dup_r0} <= $rpos-$buffer ) {
	$readsMutWtLeft{$cname}{$readname} = 1;
	$readsMutWt{$cname}{$readname} = 1;
	$readsPassed{$readname} = 1;
    }
    if( !exists($readsMutWtRight{$cname}{$readname}) && 
	$candidatedets{$cname}{dup_r1}+$candidatedets{$cname}{itdsize} >= $line_cells[3]+$buffer && 
	$candidatedets{$cname}{dup_r1}+$candidatedets{$cname}{itdsize} < $rpos-$buffer ) {
	$readsMutWtRight{$cname}{$readname} = 1;
	$readsMutWt{$cname}{$readname} = 1;
	$readsPassed{$readname} = 1;
    }

    # FUTURE USE
    if( 0 && ( $line_cells[1] & 64 ) == 64 ) {
	if( ($line_cells[1] & 16) == 16 ) { 
	    $locFwd0{$readname}=$rpos; $locFwd1{$readname}=$line_cells[3]-1;
	} else {
	    $locFwd0{$readname}=$line_cells[3]-1; $locFwd1{$readname}=$rpos;
	}
    } elsif( 0 && ( $line_cells[1] & 128 ) == 128 ) {
	if( ($line_cells[1] & 16) == 16 ) {
	    $locRev0{$readname}=$rpos; $locRev1{$readname}=$line_cells[3]-1;
	    if( $candidatedets{$cname}{mutbk_s1} > $line_cells[3]+$buffer && 
		$candidatedets{$cname}{mutbk_s1} <= $rpos-$buffer ) {
		$readsMutByPrimerR2{$cname}{R}{$rpos}{$readname} = 1;
		$readsMutByR2{$cname}{$readname} = 1;
	    }
	} else {
	    $locRev0{$readname}=$line_cells[3]-1; $locRev1{$readname}=$rpos;
	    if( $candidatedets{$cname}{mutbk_s0} >= $line_cells[3]+$buffer && 
		$candidatedets{$cname}{mutbk_s0} < $rpos-$buffer ) {
		$readsMutByPrimerR2{$cname}{L}{$line_cells[3]}{$readname} = 1;
		$readsMutByR2{$cname}{$readname} = 1;
	    }
	}
    }
    # END FUTURE USE
}
close FI;

foreach $readname ( keys %resolvedWtPaired ) {
    $readsPassed{$readname} = 1;
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
	    }
	    if( !exists( $readsWtRight{$cname}{$readname} ) &&
		$candidatedets{$cname}{dup_r1} >= $line_cells[3]+$buffer && 
		$candidatedets{$cname}{dup_r1} < $rpos - $buffer ) {
		$readsWtRight{$cname}{$readname} = 1;
		$readsWt{$cname}{$readname} = 1;
	    }

	    # FUTURE USE
	    if( 0 && ($line_cells[1] & 128) == 128 ) {
		if( $candidatedets{$cname}{dup_r0} > $line_cells[3]+$buffer && 
		    $candidatedets{$cname}{dup_r0} <= $rpos - $buffer ) {
		    if( ($line_cells[1] & 16 ) == 16 ) {
			$readsWtByPrimerR2{$cname}{R}{$rpos}{$readname} = 1;
			$readsWtByR2{$cname}{$readname} = 1;
		    }
		}
		if( $candidatedets{$cname}{dup_r1} >= $line_cells[3]+$buffer && 
			 $candidatedets{$cname}{dup_r1} < $rpos - $buffer ) {
		    if( ($line_cells[1] & 16 ) == 0 ) {
			$readsWtByPrimerR2{$cname}{L}{$line_cells[3]}{$readname} = 1;
			$readsWtByR2{$cname}{$readname} = 1;
		    }
		}
	    }
	    # END FUTURE USE
	}
        # FUTURE USE
	if( 0 && ($line_cells[1] & 128) == 128 ) {
	    if( ($line_cells[1] & 16 ) == 16 ) {
		$allreadsWtByR2{R}{$rpos}{$readname} = 1;
	    } else {
		$allreadsWtByR2{L}{$line_cells[3]}{$readname} = 1;
	    }
	}
	# END FUTURE USE
    }
}

foreach( keys %readsMut ) { $coverageMut{$_} = scalar( keys %{$readsMut{$_}} ); }
foreach( keys %readsMutLeft ) { $coverageMutLeft{$_} = scalar( keys %{$readsMutLeft{$_}} ); }
foreach( keys %readsMutRight ) { $coverageMutRight{$_} = scalar( keys %{$readsMutRight{$_}} ); }
foreach( keys %readsMutWt ) { $coverageMutWt{$_} = scalar( keys %{$readsMutWt{$_}} ); }
foreach( keys %readsMutWtLeft ) { $coverageMutWtLeft{$_} = scalar( keys %{$readsMutWtLeft{$_}} ); }
foreach( keys %readsMutWtRight ) { $coverageMutWtRight{$_} = scalar( keys %{$readsMutWtRight{$_}} ); }
foreach( keys %readsWt ) { $coverageWt{$_} = scalar( keys %{$readsWt{$_}} ); }
foreach( keys %readsWtLeft ) { $coverageWtLeft{$_} = scalar( keys %{$readsWtLeft{$_}} ); }
foreach( keys %readsWtRight ) { $coverageWtRight{$_} = scalar( keys %{$readsWtRight{$_}} ); }

my %coverageRadius = ();
foreach my $cname (keys %coverageMut) {
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

my %coverageDets = ();
foreach my $cname ( keys %coverageRadius ) {
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
}

my $vcfstring = $vcfheader . $shortkey . "\n";
foreach( keys %coverageDets ) {
    my $mutCount = ($coverageMutLeft{$_}+$coverageMutRight{$_})/2;
    my $mwCount = ($coverageWtLeft{$_}+$coverageWtRight{$_}+$coverageMutWtLeft{$_}+$coverageMutWtRight{$_})/2;
    $vcfstring .= proc_vcf_from_keyfull($itdkeys{$_}, $candidatedets{$_}{keycdot} ) . "\t" .
	"GENE=FLT3;STRAND=-;SVLEN=" . $candidatedets{$_}{itdsize} . 
	";CDS=" . $candidatedets{$_}{keycdot} . ";AA=" . $candidatedets{$_}{keypdot} .
#	";ALTCDS=" . $candidatedets{$_}{keycdotalt} . ";ALTAA=" . $candidatedets{$_}{keypdotalt} .
	";AF=" . sprintf( "%.4g", $mutCount/$mwCount ) . ";DP=" . $mwCount . ";VD=" . $mutCount . ";SAMPLE=" . $shortkey .
	"\tAF:DP:VD\t" . sprintf( "%.4g", $mutCount/$mwCount ) . ":" . $mwCount . ":" . $mutCount . "\n";
}
open( FO, ">" . $fbase . "_ITD.vcf" ) or die $!;
print FO $vcfstring;
close FO;

if( $debug ) {
  open(FO, ">" . $fbase . "_itddets.txt" ) or die $!;
  foreach( sort {$coverageMut{$b}<=>$coverageMut{$a}} keys %coverageMut ) {
    print FO $_ . "\t" . $itdkeysfiltered{$_} . "\t" . ($coverageMutLeft{$_} + $coverageMutRight{$_})/2 . "\t";
    if( exists( $coverageWt{$_} ) ) { print FO ($coverageWtLeft{$_} + $coverageWtRight{$_})/2; }  else { print FO "NA"; }
    print FO "\t" . $candidatedets{$_}{dup_r0} . "\t" . $candidatedets{$_}{dup_r1} . "\t" .
	$candidatedets{$_}{inslen} . "\n";
  }
  close FO;

  my $coverageOut = ""; my $tmpOut;
  foreach( keys %coverageDets ) {
    $tmpOut = $_ . "\t" . $candidatedets{$_}{itdsize} . "\t" . $candidatedets{$_}{ie};
    if( $proc_output_rhp ) {
	$tmpOut .= "\t" . $coverageMutLeft{$_} . "\t" . $coverageMutRight{$_} . "\t" . 
	    $coverageWtLeft{$_} . "\t" . $coverageWtRight{$_};
    }
    $tmpOut .= "\t" . ($coverageMutLeft{$_}+$coverageMutRight{$_})/2 . "\t" . 
	($coverageWtLeft{$_}+$coverageWtRight{$_}+$coverageMutWtLeft{$_}+$coverageMutWtRight{$_})/2 . "\t" . 
	$coverageDets{$_}{"0"} ."\t" . $coverageDets{$_}{"1-5"} ."\t" . $coverageDets{$_}{"6-10"} ."\t" . 
	$coverageDets{$_}{"11-25"} ."\t" . $coverageDets{$_}{"26-50"} ."\t" . $coverageDets{$_}{">50"} . "\t" .
	$candidatedets{$_}{dup_r0} . "\t" . $candidatedets{$_}{dup_r1} . "\t" .
	$candidatedets{$_}{keycdot} . "\t" . $candidatedets{$_}{keypdot} . "\t" .
	$candidatedets{$_}{keycdotalt} . "\t" . $candidatedets{$_}{keypdotalt} . "\n";
    $coverageOut .= $tmpOut;
  }
  print $coverageOut; print $vcfstring . "\n";
}

#### GENERATE TVIEW HTML FILES
if( $web ) {
  foreach( keys %coverageMut ) {
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

  foreach( keys %coverageMut ) {
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

    my $replaceseq = ""; my $fo;
    if( exists( $candidatedets{$_} ) ) {
	$replaceseq = $candidatedets{$_}{dupseq};
	my $ckey = $candidatedets{$_}{keycdot}; $ckey =~ s/c\./c/; $ckey =~ s/\//ins/g;
	$fo = $fbase . "_" . $ckey . "_length=" . $candidatedets{$_}{itdsize} . "_reads=" . $coverageMut{$_} . ".html";
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
	    print FO $_;
	}
	close FI;
	close FO;

	if( system( sprintf( "rm %s", $fi ) ) ) { die "Error removing htmlold. Exiting..."; }
    }
  }

  if( !$debug && (
    system( "rm " . $fbase . "_to_candidate_ITDs.mutreads.sam" ) ||
    system( "rm " . $fbase . "_to_candidate_ITDs.mutreads.bam" ) ||
    system( "rm " . $fbase . "_to_candidate_ITDs.mutreads.sorted.bam" ) ||
    system( "rm " . $fbase . "_to_candidate_ITDs.mutreads.sorted.bam.bai" ) ) ) {
    die "Error removing to_candidte_ITDs mutreads alignment files. Exiting..." ;
  }
}

