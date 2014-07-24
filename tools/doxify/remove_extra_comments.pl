#!/usr/bin/perl

# Remove any double sets of lines starting ! ****** 
# E.g. removes the lines
# ! ****      
# ! ****
# These extra lines get created by the fixcomments.pl script because it's not 
# always guaranteed that a header will have them at the start and end and this 
# it will always add in the header lines. 

$DEBUG=$ARGV[2];
if ($DEBUG eq 1) {
    print "remove_extra_comments.pl running with INPUT & OUTPUT = ";
    print $ARGV[0];
    print " ";
    print $ARGV[1];
    print "\n";
}

open ($INPUT , "<" , $ARGV[0]) or die "Cant open $ARGV[0] $!";
open ($OUTPUT , ">" , $ARGV[1]) or die "Cant create $ARGV[1] $!";

$count = 0;
while (<$INPUT>) # While there are still lines to read
{
    $currline = $_; # currline = the line we've just read in  

# We want to strip out any double sets (i.e. lines next to each other of lines starting ! **** 
    if ($currline =~ m/^!\s+\*+/i) {
	$nextline = <$INPUT>;  	# If we get a match read in the next line
	while ($nextline eq $currline) {  # Keep reading until we find a different line
	    $count = $count + 1;  # Count up how many lines get stripped out
	    $nextline = <$INPUT>;
 	}
        print $OUTPUT $currline; # Dump out the ! ** and the next line that is not a duplicate
        print $OUTPUT $nextline; 
    }	else {	
	print $OUTPUT $currline;  # Print out all the remainder of the file
    }
}

if ($DEBUG eq 1) {
    print "remove_extra_comments.pl:: removed a total of ", $count, " duplicate comment lines from file = ", $ARGV[0], "\n";
}

close $OUTPUT;
close $INPUT;

