#!/usr/bin/perl

# Want to locate any lines where we have an ampersand at the end of a line followed
# by an ampersand at the start of the next line as these causes the fixcomments.pl 
# script to break as it uses the ampersand to decide whether there is more data to 
# read. 
# An example would be (from cp_lbfgs.F): 
#      SUBROUTINE setulb(n, m, x, l, u, nbd, f, g, factr, pgtol, wa, iwa,&
#     &                 task, iprint, csave, lsave, isave, dsave)

$DEBUG=$ARGV[2];
if ($DEBUG eq 1) {
    print "DEBUG = YES \n";
    print "remove_double_ampersands.pl running with INPUT & OUTPUT = ";
    print $ARGV[0];
    print " ";
    print $ARGV[1];
    print "\n";
}
open ($INPUT , "<" , $ARGV[0]) or die "Cant open $ARGV[0] $!";
open ($OUTPUT , ">" , $ARGV[1]) or die "Cant create $ARGV[1] $!";

$count = 0;
$linesread = 0;
$countampend = 0;
$countampstart = 0;
$firsttime = 1;
while (<$INPUT>) # While there are still lines to read
{

# First time we enter the loop we need to read in two lines
    if ($firsttime eq 1) {
	$currline = $_; # currline = the line we've just read in  
	$linesread = $linesread + 1;
	$nextline = <$INPUT>; # read in the next line too 
	$linesread = $linesread + 1;
	$firsttime = 0;
    } else {
# If not first iteration need to copy nextline into currline and read currline from $_
	$currline = $nextline;
	$nextline = $_;
	$linesread = $linesread + 1;
    }

    if ( ($currline =~ m/&\s*$/i) && ($nextline =~ m/^\s*&\s*\w+/i) ) {
	$nextline =~ s/&/ /;
	if ($DEBUG eq 1) { # Only print when debugging the script
	    print "remove_double_ampersands.pl:: CURRLINE is ",$currline;
	    print "remove_double_ampersands.pl:: NEXTLINE is ",$nextline;
	}
	    $count = $count + 1;
    }
    print $OUTPUT $currline; # Write all lines except the last one to file - initial & replaced with space

}

print $OUTPUT $nextline;  # Write the last line out to file


if ( $DEBUG eq 1 ) { # Only print when debugging the script 
    print "remove_double_ampersand.pl:: ", $linesread, " lines read in, found a total of ", $count, " lines with the double ampersand issue in file = ", $ARGV[0], " \n";
}


# Close the files
close $OUTPUT;
close $INPUT;

