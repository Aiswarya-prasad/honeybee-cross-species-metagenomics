#! /usr/bin/perl -w
use strict;
use List::Util qw(max min sum);
use Getopt::Long;
use Pod::Usage;
pod2usage(-verbose => 2) if ((@ARGV==0)&&(-t STDIN));

=head1 SYNOPSIS

perl summarize_snps_host.pl table_core_length.txt firm4_1.filtered.freq firm4_1_corecov_coord.txt

=head1 DESCRIPTION

This script will quantify the number and fraction of SNVs per sample, and across all profiled samples. The data is separated by host-affiliation, in order to compare diversity across hosts. After splitting the data by host, the script checks whether the SNV is still polymorphic across samples.

Two output files are generated:

[SDP]"_tot_var.txt": total fraction of SNVs per host
[SDP]"_sample_var_host": fraction of SNVs per sample

For the total fraction file, the number of samples profiled and the cumulative coverage is given in the last two columns.

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
my %var_per_sample;
my %host_var;
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
	    if ($host ne "Ac") {
		$host = "Am";
	    }
	    $host_samples{$host} = [] unless exists $host_samples{$host};
	    push @{$host_samples{$host}}, $sample;
	}
	++$line_count;
	}
    else {
	my $snp_desc=shift @split;
	my @split_snp = split(":",$snp_desc);
	my $snp_pos = $split_snp[1];
	my $i = 0;
	my %host_line = ();
	while($i < $nb_samples) {
	    my $host = substr $samples[$i], 0,2;
	    if ($host ne "Ac") {
		$host = "Am";
	    }
	    $host_line{$host} = [] unless exists $host_line{$host}; #gather all values associated with host in array
	    push @{$host_line{$host}},$split[$i];
	    if ($split[$i] > 0 && $split[$i] <= 0.95) {
		++$var_per_sample{$samples[$i]}{$snp_pos};
	    }
	    ++$i;
	}
	foreach my $host (keys %host_line) {
	    my @var = @{$host_line{$host}};
	    my $nb_var = @var;
	    my $sum_var = sum @var;
	    unless ($nb_var == $sum_var || $sum_var == 0) { #determine whether a snp-pos varies across a host, using above generated host-array
		++$host_var{$host}{$snp_pos};
	    }
	}
    }
}
close FILE;

#Read the file containing ter-cov for all samples of the given taxa

open(FILE, $ARGV[2]) or die
    "Cant open $ARGV[2]!";
my @coord_file = <FILE>;
close FILE;
shift @coord_file;
my %cum_cov_host;
foreach (@coord_file) {
    chomp;
    my @split = split("\t",$_);
    my $cov_ter = $split[2];
    my $sample = $split[1];
    my $host = substr $sample, 0, 2;
    if ($cov_ter > 20 ) {
	$cum_cov_host{$host} += $cov_ter;
    }
}

#Print outfile with fraction var per host

my $tot_var_outfile = $strain."_tot_var.txt";
open OUTFILE, '>'.$tot_var_outfile or die "Cant open $tot_var_outfile!";
foreach my $host(keys %host_var) {
    my %host_var_pos = %{$host_var{$host}};
    my $nb_host_var_pos = keys %host_var_pos;
    my @host_samples = @{$host_samples{$host}};
    my $nb_host_samples = @host_samples;
    my $host_cum_cov = sprintf("%.2f", $cum_cov_host{$host});
    if ($nb_host_var_pos > 0) {
	my $fraction_var = sprintf("%.3f",($nb_host_var_pos/$core_length{$strain})*100);
	print OUTFILE $phylo,"\t",$strain,"\t",$host,"\t",$nb_host_var_pos,"\t",$fraction_var,"\t",$nb_host_samples,"\t",$host_cum_cov,"\n";
    }
    else {
	print OUTFILE $phylo,"\t",$strain,"\t",$host,"\t0\t0\t",$nb_host_samples,"\t",$host_cum_cov,"\n";
    }
}
close OUTFILE;

#Print outfile with fraction var per host, per sample

my $sample_var_outfile = $strain."_sample_var_host.txt";
open OUTFILE, '>'.$sample_var_outfile or die "Cant open $sample_var_outfile!";
foreach my $host( keys %host_samples) {
    my @host_samples = @{$host_samples{$host}};
    foreach my $sample(@host_samples) {
	my $colony = substr $sample,0,2;
	if (exists  $var_per_sample{$sample}) {
	    my $nb_var_pos = keys %{$var_per_sample{$sample}};
	    my $fraction = sprintf("%.2f",($nb_var_pos/$core_length{$strain})*100);
	    print OUTFILE $phylo,"\t",$strain,"\t",$host,"\t",$colony,"\t",$sample,"\t",$nb_var_pos,"\t",$fraction,"\n";
	}
	else {
	    print OUTFILE $phylo,"\t",$strain,"\t",$host,"\t",$colony,"\t",$sample,"\t0\t0\n";
	}
    }
}
close OUTFILE;
