#! /usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;
pod2usage(-verbose => 2) if ((@ARGV==0)&&(-t STDIN));

=head1 SYNOPSIS

perl calc_shared_SNV_fraction.pl firm4_1_filt.freq

=head1 DESCRIPTION

This script will calculate jaccard indices (i.e. the shared fraction) of SNVs between pairs of samples. 

Only presence of alleles is considered (not absences, which are far more frequent). Thus, if both samples contain a variant, its a share, if one of them contain the variant is a non-share, if none of them have the variant its not counted. The following output-file is generated:

[SDP]"_filt_shared_fraction.txt" 

Note, this analysis can be time-consuming, depending on the number of pairwise sample combinations (and the number of SNVs). Computational progress will be printed continously to stdout.

=cut

GetOptions (
         'help' => sub { pod2usage( -exitstatus => 0, -verbose => 2 ) },
    ) or pod2usage(2);

#Sub-routine for calculating the extent of shared SNPs (jaccard index) for all pairwise combinations of samples, **per snp line**.

sub calc_share {
    my @sample_list = @{$_[0]};
    my @snp_freq = split("\t",$_[1]);
    shift @snp_freq;
    my $nb_samples = @sample_list;
    my %nb_comparison;
    my %nb_share;
    my $i = 0;
    while($i < $nb_samples) {
	my $j=0;
	my $snp_data_sample1 = $snp_freq[$i];
	my $sample1 = $sample_list[$i];
	while($j < $nb_samples) {
	    my $snp_data_sample2 = $snp_freq[$j];
	    my $sample2 = $sample_list[$j];
	    if ($snp_data_sample1 > 0 || $snp_data_sample2 > 0) { #presence in at least one of the samples
		$nb_comparison{$sample1}{$sample2}=1;
		if ($snp_data_sample1 > 0 && $snp_data_sample2 > 0) {
		    $nb_share{$sample1}{$sample2}=1; #Shared allele
		}
		else {
		    $nb_share{$sample1}{$sample2}=0; #Allele present in one, but not both
		}
	    }
	    else { #Allele absent in both
		$nb_share{$sample1}{$sample2}=0;
		$nb_comparison{$sample1}{$sample2}=0;
	    }
	    ++$j;
	}
    ++$i;
    }
    return(\%nb_share, \%nb_comparison);
}



#Open the file with SNP data. Using a sub-routine on each line: for all pairwise sample-comparisons, check whether the SNP was detected in both samples. 
open(FILE,$ARGV[0]) or
    die "Cant open $ARGV[0]!";
my $line_count=0;
my @samples;
my %share_sum;
my %nb_comparison_sum;
   my $verbose_count = 5000;
while(<FILE>) {
    chomp;
    if ($line_count>$verbose_count) {
	print "Finished with dist-calc for: ",$verbose_count," snp-lines\n";
	$verbose_count += 5000;
    }
    if ($line_count==0) {
	my @split_header = split("\t",$_);
	shift @split_header;
	@samples = @split_header;
    }
    else {
	my ($share, $nb_comparison) = calc_share(\@samples,$_);
	my %share = %{$share};
	my %nb_comparison = %{$nb_comparison};
	foreach my $sample1(@samples) {
	    foreach my $sample2(@samples) {
		$share_sum{$sample1}{$sample2}+= $share{$sample1}{$sample2};
		$nb_comparison_sum{$sample1}{$sample2} += $nb_comparison{$sample1}{$sample2};
	    }
	}
    }
    ++$line_count;
}

#Print an output file
my @split_name = split('\.',$ARGV[0]);
my $outfile = $split_name[0]."_shared_fraction.txt";
open OUTFILE, '>'.$outfile or
    die "Cant open $outfile!";

print OUTFILE "Sample1\tSample2\tNb_shared_alleles\tNb_scored_alleles\tShared_fraction\n";
foreach my $sample1(@samples) {
    foreach my $sample2(@samples) {
	if ($nb_comparison_sum{$sample1}{$sample2} == 0) {
	    print "Warning: too few comparisons to calculate jaccard-index for samples: ",$sample1,"\t",$sample2,"\n";
	    print OUTFILE $sample1,"\t",$sample2,"\t",$share_sum{$sample1}{$sample2},"\t",$nb_comparison_sum{$sample1}{$sample2},"\t0\n";
	}
	else {
	    unless ($sample1 eq $sample2) {
		my $shared_fraction = sprintf("%.4f",($share_sum{$sample1}{$sample2}/$nb_comparison_sum{$sample1}{$sample2}));
		print OUTFILE $sample1,"\t",$sample2,"\t",$share_sum{$sample1}{$sample2},"\t",$nb_comparison_sum{$sample1}{$sample2},"\t",$shared_fraction,"\n";
	    }
	}
    }
}
close OUTFILE;
