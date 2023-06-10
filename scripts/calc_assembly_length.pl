#! /usr/bin/perl -w
use strict;
use Bio::Seq;
use Bio::SeqIO;
use Getopt::Long;
use Pod::Usage;
pod2usage(-verbose => 2) if ((@ARGV==0)&&(-t STDIN));

=head1 SYNOPSIS

perl calc_assembly_length.pl AmAi01_1_parsed.fasta

=head1 DESCRIPTION

This script is used to calculate the total sequence length of a multi-fasta file, it will print the length (bp) to stdout. The script was used to calculate assembly lengths for rarefied filtered metagenomic assemblies.

=cut

GetOptions(
    'help' => sub { pod2usage( -exitstatus => 0, -verbose => 2 ) },
);

my $input_file = shift;
my $seq_io = Bio::SeqIO->new(
                      -file => $input_file, 
                       -format => 'fasta',
                       );
my $tot_length=0;
while (my $seq_obj = $seq_io->next_seq()) {
    my $seq_length = $seq_obj->length;
    $tot_length += $seq_length;
}

print $tot_length,"\n";
