#! /usr/bin/perl -w
use List::Util qw( min max );
use Getopt::Long;
use Pod::Usage;

=head1 SYNOPSIS

perl filter_sam_aln_length_unmapped.pl file.sam > file_unmapped.sam

=head1 DESCRIPTION

This script is used to extract unmapped reads from a sam-file. Here, unmapped reads are defined as the reads which do not meet the alignment length criterium (alignment length <= 50), as extracted from the CIGAR string in the sam-file.

=cut

GetOptions(
    'help' => sub { pod2usage( -exitstatus => 0, -verbose => 2 ) },
);

open (FILE, $ARGV[0]) or
  die "Can't open $ARGV[0]: $!";

my @match_numbers;
while(<FILE>) {
    chomp;
    # print the header line as it is
    my @split = split(" ",$_);
    if (/^\@/) {
	print $_,"\n";
    }
    # if it is not a header
    else {
	my $cigar = $split[5];
	my @matches = ($cigar =~ m/([0-9]+M)/g);
	# count the number of matches
	my $nb_matches = @matches;
	if ($nb_matches > 0) { #mapped
	    foreach my $match(@matches) {
		push @match_numbers, substr $match,0,-1;
	    }
	    my $max = max @match_numbers;
	    # if the maximum match length is less 50
	    if ($max <= 50) {
		print $_,"\n";
	    }
	}
	else { #not mapped
	    print $_,"\n";
	}
	@match_numbers = ();
    }
}
