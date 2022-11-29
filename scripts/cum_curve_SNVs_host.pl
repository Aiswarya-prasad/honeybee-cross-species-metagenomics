#! /usr/bin/perl -w
use strict;
use List::Util qw(max min sum shuffle);
use Getopt::Long;
use Pod::Usage;
pod2usage(-verbose => 2) if ((@ARGV==0)&&(-t STDIN));

=head1 SYNOPSIS

perl cum_curve_SNVs_host.pl table_core_length.txt firm4_filt.freq

=head1 DESCRIPTION

This script will calculate the cumulative fraction of SNVs relative to the number of samples. By default, 10 random sample orders are used (using the perl-module "shuffle"). The following output-file is generated:

[SDP]"_cum_curve.txt"

The columns represent:
1. Host affiliation
2. curve-id (from 1 to 10)
3. cumulative number of samples
4. cumulative fraction SNVs
5. SDP name

The script will run regardless of the number of samples profiled, but the output file is obviously not very informative if the number of samples is low.

=cut

GetOptions(
    'help' => sub { pod2usage( -exitstatus => 0, -verbose => 2 ) },
) or pod2usage(2);

#Read the file containing the core-lengths for all taxa
open(FILE, $ARGV[0]) or
    die "Cant open $ARGV[0]!";
my @file = <FILE>;
shift @file;
my %core_length;
foreach my $line(@file) {
    chomp($line);
    my @split = split("\t",$line);
    $core_length{$split[0]} = $split[2];
}
close FILE;

#Read the snp-file and gather data for analysis

open(FILE, $ARGV[1]) or
    die "Cant open $ARGV[1]!";
my $strain=substr $ARGV[1],0,-10;
my @split_strain = split("_",$strain);
my $phylo = $split_strain[0];
my @samples;
my $nb_samples;
my %var_per_sample; #contains all snps occurring with 0.1 to < 1 relative abundance, per sample
my %tot_var;
my %host_samples;
my $line_count=0;
while(<FILE>) {
    chomp;
    my @split = split("\t",$_);
    if ($line_count==0) {
	shift @split;
	@samples = @split;
	$nb_samples = @samples;
	foreach my $sample(@samples) {
	    my $host = substr $sample, 0, 2;
	    $host_samples{$host} = [] unless exists $host_samples{$host};
	    push @{$host_samples{$host}}, $sample;
	}
	++$line_count;
	}
    else {
	my $snp_desc=shift @split;
	my @split_snp = split(":",$snp_desc);
	my $snp_pos = $split_snp[1];
	++$tot_var{$snp_pos};
	my $i = 0;
	while($i < $nb_samples) {
	    my $freq = sprintf("%.2f",$split[$i]);
            $var_per_sample{$samples[$i]}{$snp_pos} = [] unless exists $var_per_sample{$samples[$i]}{$snp_pos};
            push @{$var_per_sample{$samples[$i]}{$snp_pos}},$freq;
	    ++$i;
	}
    }
}
close FILE;

my @sorted_pos = sort {$a <=> $b} keys %tot_var;
my $nb_curves = 10;
my %host_var_curves;
my $cum_curve_outfile = $strain."_cum_curve.txt";
open OUTFILE, '>'.$cum_curve_outfile or die "Cant open $cum_curve_outfile!";
foreach my $host(keys  %host_samples) {
    my @host_samples = @{$host_samples{$host}};
    my $nb_samples = @host_samples;
    my @index_range = (0..($nb_samples-1));
    my $i = 0;
    my %curves = ();
    while($i < $nb_curves) {
	$curves{$i} = [];
	my @shuffle_index = shuffle(@index_range);
	my %var_freq=();
	my %var_pos=();
	my @samples = ();
	foreach my $index(@shuffle_index) {
	    my $sample = $host_samples[$index];
	    push @samples,$sample;
	    foreach my $pos(@sorted_pos) {
		$var_freq{$pos} = [] unless exists $var_freq{$pos};
		my @sample_var_freq_pos = @{$var_per_sample{$sample}{$pos}};
		foreach my $freq(@sample_var_freq_pos) {
		    push @{$var_freq{$pos}}, $freq;
		}
	    }
	    foreach my $pos(@sorted_pos) {
		my @freqs = @{$var_freq{$pos}};
		my $nb_freqs = @freqs;
		my $sum_freq = sum @freqs;
		unless ($sum_freq == $nb_freqs || $sum_freq == 0) {
		    ++$var_pos{$pos};
		}
	    }
	    my $nb_var_pos = keys %var_pos;
	    my $fraction_var = sprintf("%.3f",($nb_var_pos/$core_length{$strain})*100);
	    my $curve_nb = $i +1;
	    my $curve_id = "Curve_".$curve_nb;
	    my $nb_samples = @samples;
	    print OUTFILE $host,"\t",$curve_id,"\t",$nb_samples,"\t",$fraction_var,"\t",$strain,"\n";
	}
	++$i;
    }
}
close OUTFILE;
