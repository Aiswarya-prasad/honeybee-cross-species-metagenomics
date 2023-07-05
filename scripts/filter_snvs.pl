#! /usr/bin/perl -w
use List::Util qw(sum);
use Getopt::Long;
use Pod::Usage;
pod2usage(-verbose => 2) if ((@ARGV==0)&&(-t STDIN));

=head1 SYNOPSIS

perl filter_snvs.pl fper.vcf

=head1 DESCRIPTION

This script is used to filter a vcf-file, working as follows:

Step 1: Using the sub-routine "extract_dp_data", the relative coverage of a candidate SNV is calculated for each sample. The coverage is set to zero if it less than 0.1 (10%). Additionally, the fraction of samples with missing data (due to low coverage) is calculated. If the missing fraction is larger than 0.1 (10%), the SNV is skipped

Step 2: Using the sub-routine "is_polymorphic", the script will check whether the candidate SNV is still polymorphic after rare SNVs are removed. If all samples have the values 0 or 1 (i.e. no polymorphism, or only polymorphic relative to the reference genome), the SNV is skipped

Step 3: If an SNV passes the first two steps, it is printed to the output file

The output file is named [SDP]"_filt.freq". In the first column, the reference genome, SNV positiion and base change relative to the reference is indicated. The remaining columns display the relative proportion of the SNV within each sample. In cases where there is more than one two alleles for the same position, multiple lines are generated.

Finally, the script will print the total number of candidate SNVs found to be polymorphic, and the number of positions that were skipped due to missing data.

=cut

GetOptions(
    'help' => sub { pod2usage( -exitstatus => 0, -verbose => 2 ) },
) or pod2usage(2);

sub extract_dp_data {
    my @data = @{$_[0]};
    my $nb_tabs = @data;
    my $missing_data_count=0;
    my @data_AO;
    foreach my $tab(@data) {
	my @split_tab = split(":",$tab);
	my $AO = $split_tab[5];
	if ($AO =~ /\./) {
	    ++$missing_data_count;
	    push @data_AO, -1;
	}
	else {
	    my $DP = $split_tab[1];
	    my $allele_relcov = $AO/$DP;
	    if ($allele_relcov <= 0.1) { #set rare SNVs to zero
		push @data_AO, 0;
	    }
	    else {
		push @data_AO,$allele_relcov;
	    }
	}	    
    }
    my $missing_fraction = $missing_data_count/$nb_tabs;
    return($missing_fraction, \@data_AO);
}

sub is_polymorphic {
    my @all_data = @{$_[0]};
    my @data;
    foreach my $tab(@all_data) {
	unless ($tab == -1) {
	    push @data,$tab;
	}
    }
    my $nb_data = @data;
    my $sum_data = sum @data;
    my $sub_result=1; #1 equals "yes", 0 equals "no"
    if ($sum_data == $nb_data || $sum_data == 0) {
	$sub_result=0;
    }
    return($sub_result);
}

my @split_filename = split('\.',$ARGV[0]);
my $outfile = $split_filename[0]."_filt.freq";
open OUTFILE, '>'.$outfile or die $!;

open(FILE, $ARGV[0]) or die
    "Cant open $ARGV[0]!";
my @samples;
my $nb_ND=0;
my $nb_var=0;
while(<FILE>) {
    chomp;
    next if ($_ =~ /##/);
    my @split = split("\t",$_);
    if ($split[0] eq "#CHROM") {
	my $nb_tabs = @split;
        @samples = splice @split,9;
	print OUTFILE "\t",join("\t",@samples),"\n";
    }
    else {
	
	my @data = splice @split,9;
	my ($missing_fraction,$data_AO) = extract_dp_data(\@data);
	my @data_AO = @{$data_AO};
	++$nb_ND if ($missing_fraction >= 0.1);
	next if ($missing_fraction >= 0.1);
	my $polymorphic_state = is_polymorphic(\@data_AO);
	next if ($polymorphic_state == 0); #skip position if not polymorphic after removing rare SNVs
	++$nb_var;
	my $genome_id = $split[0];
	my $pos = $split[1];
	my $ref_base = $split[3];
	my $alt_base = $split[4];
	my $snv_info_line = $genome_id.":".$pos.":".$ref_base.">".$alt_base;
	print OUTFILE $snv_info_line,"\t",join("\t",@data_AO),"\n";
    }
}

print "###\n";
print "Processing ",$ARGV[0],"\n";
print "A total of ",$nb_var," positions were found to be polymorphic (exclusing rare SNVs)\n";
print "A total of ",$nb_ND," positions were skipped, because of missing coverage data in more than 10% of the samples\n###\n\n";
