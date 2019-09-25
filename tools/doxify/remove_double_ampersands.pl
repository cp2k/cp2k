#!/usr/bin/env perl

=head1 NAME

remove_double_ampersands.pl

=head1 SYNOPSIS

remove_double_ampersands.pl [options] [infile]

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

Want to locate any lines where we have an ampersand at the end of a line followed
by an ampersand at the start of the next line as these causes the fixcomments.pl
script to break as it uses the ampersand to decide whether there is more data to
read.
An example would be (from cp_lbfgs.F):
     SUBROUTINE setulb(n, m, x, l, u, nbd, f, g, factr, pgtol, wa, iwa,&
    &                 task, iprint, csave, lsave, isave, dsave)

=cut

use strict;
use warnings;
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
my $countampend = 0;
my $countampstart = 0;
my $currline = <>; # First time we enter the loop we need to read in two lines, read the first now
my $linesread = 1;

while (my $nextline=<>) { # While there are still lines to read, either from STDIN or the files specified
    $linesread += 1;

    if ( ($currline =~ m/&\s*$/i) && ($nextline =~ m/^\s*&\s*\w+/i) ) {
        $nextline =~ s/&/ /;
        print_debug("remove_double_ampersands.pl:: CURRLINE is ", $currline);
        print_debug("remove_double_ampersands.pl:: NEXTLINE is ", $nextline);
        $count += 1;
    }

    print $currline; # Write all lines except the last one - initial & replaced with space
    $currline = $nextline
}

print $currline;  # Write the last line (unaltered)

print_debug("remove_double_ampersand.pl:: ", $linesread, " lines read in, found a total of ", $count, " lines with the double ampersand issue in file", " \n");
