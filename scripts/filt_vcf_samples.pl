#! /usr/bin/perl -w
use Getopt::Long;
use Pod::Usage;
pod2usage(-verbose => 2) if ((@ARGV==0)&&(-t STDIN));

=head1 SYNOPSIS

perl filt_vcf_samples.pl api_1_corecov_coord.txt freebayes_candSNVs_core.vcf

=head1 DESCRIPTION

This script is used as part of the pipeline for generating vcf-files for each SDP, subset for samples with at least 20x terminus coverage. The script works as follows:

Step 1: Read the file containing the terminus coverage for all samples ([SDP]"_corecov_coord.txt", and get the names of samples meeting the threshold

Step 2: Read the vcf-file, and generate a temporary vcf-file ([SDP]"_temp.vcf"), subset for lines corresponding to the SDP indicated by the first input file

Step 3: Print the names of the filtered samples to stdout

Thus, the script will output a temporary vcf-file, and provide the input for vcflib.

=cut

GetOptions(
    'help' => sub { pod2usage( -exitstatus => 0, -verbose => 2 ) },
) or pod2usage(2);

my %ref_ids = (
'firm5_1' => 'Ga0133563',
'firm5_2' => 'Ga0133561',
'firm5_3' => 'Ga0133562',
'firm5_4' => 'Ga0133564',
'firm5_5' => 'Ga0226852',
'firm5_6' => 'Ga0226847',
'firm5_7' => 'Ga0312303',
'gilli_1' => 'GAPWK',
'gilli_2' => 'Ga0133555',
'gilli_3' => 'Ga0133557',
'gilli_4' => 'Ga0307802',
'gilli_5' => 'Ga0307800',
'gilli_6' => 'Ga0307803',
'gilli_7' => 'A9G09', #No complete ref
'gilli_8' => 'A9G13', #No complete ref
'gilli_9' => 'Ga0227302', #No complete ref
'snod_1' => 'SALWKB2',
'snod_2' => 'BHC53', #No complete ref, contigs have been re-ordered
'snod_3' => 'SALWKB12',#No complete ref
'bifido_1' => 'BAST',
'bifido_2' => 'BINDI',
'bifido_3' => 'Ga0098206',
'firm4_1' => 'Ga0326456',
'firm4_2' => 'Ga0072399',
'fper' => 'Ga0077910',
'bapis' => 'BBC0122',
'api_1' => 'Ga0312307',
'api_2' => 'C4S76', #No complete ref
'api_3' =>  'Ga0061079', #No complete ref
'com_1' => 'Ga0216351',
'com_2' => 'CIN', #No complete ref
'com_3' => 'Ga0248239', #No complete ref
'bom_2' => 'Ga0216357', #No complete ref
'bom_3' => 'Ga0308518', #No complete ref
'bom_1' => 'Ga0372762',
'lkun' => 'Ga0098617',
    );

my $sdp = substr $ARGV[0],0,-18;
my $genome_id = $ref_ids{$sdp};

#Read the *coord.txt file, and get the names of samples with sufficient coverage.
open(FILE, $ARGV[0]) or die
    "Cant open $ARGV[0]!";
my @high_cov_samples;
while(<FILE>) {
    chomp;
    next if ($_ =~ /Cluster/);
    my @split = split("\t",$_);
    my $sample = $split[1];
    my $ter_cov = $split[2];
    if ($ter_cov >= 20) {
	push @high_cov_samples,$sample;
    }
}
close FILE;

#Read the vcf-file, and create a temporary subset vcf-file, containing data for the SDP of interest

open OUTFILE, '>$sdp_temp.vcf' or die $!;
open(FILE, $ARGV[1]) or die
    "Cant open $ARGV[1]!";
while(<FILE>) {
    chomp;
    if ($_ =~ /#/) {
	print OUTFILE $_,"\n";
    }
    else {
	my @split = split("\t",$_);
	if ($split[0] eq $genome_id) {
	    print OUTFILE $_,"\n";
	}
    }
}
close FILE;

#Print the list of samples with high coverage, as input for bash-script
print join(" ",@high_cov_samples);
