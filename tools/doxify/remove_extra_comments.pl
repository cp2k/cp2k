#!/usr/bin/env perl

=head1 NAME

remove_extra_comments.pl

=head1 SYNOPSIS

remove_extra_comments.pl [options] [infile]

=head1 OPTIONS

=over 4

=item B<infile>

Input file to process (uses stdin if not specified)

=item B<--help|-h>

A little help

=item B<--verbose|-v>

Output debug messages (to stderr), repeat for even more output

=back

=head1 DESCRIPTION

Remove any double sets of lines starting ! ******
E.g. removes the lines
! ****
! ****
These extra lines get created by the fixcomments.pl script because it's not
always guaranteed that a header will have them at the start and end and this
it will always add in the header lines.

=cut

use strict;
use warnings;
use open qw(:std :utf8);
use Pod::Usage qw(pod2usage);
use Getopt::Long;

my $verbose = 0;
my $help = 0;

GetOptions(
    'verbose+' => \$verbose,
    'help|?' => \$help) or pod2usage(2);
pod2usage(1) if $help;

sub print_debug {
    print(STDERR "DEBUG: @_\n") if ($verbose > 0);
}


my $count = 0;

while (my $currline = <>) { # While there are still lines to read
    # We want to strip out any double sets (i.e. lines next to each other of lines starting ! ****
    if ($currline =~ m/^!\s+\*+/i) {
        my $nextline = <>;  # If we get a match read in the next line
        while ($nextline eq $currline) {  # Keep reading until we find a different line
            $count += 1;  # Count up how many lines get stripped out
            $nextline = <>;
        }
        print $currline; # Dump out the ! ** and the next line that is not a duplicate
        print $nextline;
    } else {
        print $currline;  # Print out all the remainder of the file
    }
}

print_debug("remove_extra_comments.pl:: removed a total of ", $count, " duplicate comment lines from file", "\n");
