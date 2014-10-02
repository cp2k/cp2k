#!/usr/bin/env perl

# Script to add the missing comments to CP2K
# We want to be able to add comment blocks to FUNCTION, SUBROUTINE and possibly later on 
# MODULE and TYPE. 

$DEBUG=$ARGV[2];
if ($DEBUG eq 1) {
    print "fixcomments.pl running with INPUT & OUTPUT = "; 
    print $ARGV[0];
    print " ";
    print $ARGV[1];
    print "\n";
}
open ($INPUT , "<" , $ARGV[0]) or die "Cant open $ARGV[0] $!";
open ($OUTPUT , ">" , $ARGV[1]) or die "Cant create $ARGV[1] $!";

# Initialise toggles to 0 and ensure that the header data is defined globally. If you initialise inside { } then the
# variables are only in scope inside the { }
# param = toggle for whether header contains any data in the \param fields (function or subroutine arguments)
# brief =    toggle  "       "          "           "        "      "     "  "   \brief field (description of routine)
# par =      toggle  "       "          "           "        "      "     "  "   \par field (history entry)
# version = toggle "       "          "           "        "      "     "  "   \version field (code version)
# date =    toggle  "       "          "           "        "      "     "  "   \date field (date)
# author = toggle  "       "          "           "        "      "     "  "   \author field (author)
# note = toggle  "       "          "           "        "      "     "  "   \note field (notes)
# retval = toggle  "       "          "           "        "      "     "  "   \retval field (returned value of function)
# random = toggle  "       "          "           "        "      "     "  "   random field - text that has a preceeding \ but isn't in our standard header e.g. \version
# remainder = toggle  "       "          "           "        "      "     "  "   remainder field - anything else in the comment block that doesn't have a preceeding \
# ampersand = Variable used to specifiy if the SUBROUTINE/FUNCTION definition contains an & 
# insideinterface - determines whether we are in an interface block or not
# hasretvalasarg - whether the procedure 
$param=0;  
$brief=0;
$par=0;
$date=0;
$version=0;
$author=0;
$note=0;
$retval=0;  
$random=0;
$remainder=0;
$ampersand = 0;
$insideinterface=0;
$hasretvalasarg=0;

# Variables with s at the end contains the actual text read in from header
%$params = ();
%$matched = ();
$briefs="";
$dates="";
$versions="";
$pars="";
$authors="";
$notes="";
$retvals="";
$randoms="";
$remainders="";
$buffer="";        # Used to buffer data
$oldheader="";  # Variable to keep a copy of the old header e.g. for MODULE and TYPE where we just read it in and dump it back out again.  
$leftoverlines="";

while (<$INPUT>) # While there are still lines to read in our INPUT file
{
    $currline = $_; # currline = the line we've just read in  

    # We've found a procedure header so split it up - starts with !> (^ makes it match from the start, i.e. first column)
    if ($currline =~ m/^!>/i) {
	# Each time we find a header we need to reset the toggles to 0 and reset their data values
	$param=0;
	$brief=0;
	$par=0;
	$date=0;
	$version=0;
	$author=0;
	$note=0;
	$retval=0;
        $random=0;
	$remainder=0;
	undef %params;
	undef %matched;
	$briefs="";
	$dates="";
	$versions="";
	$pars="";
	$authors="";
	$notes="";
	$retvals="";
	$buffer="";
        $randoms="";
        $remainders="";
        $oldheader="";	
	do {
            $oldheader = $oldheader . $currline;  # Keep the headers safe, may need them! We use the oldheader variable to 
                                                                 # keep the complete headers for MODULE and TYPE intact such that we
                                                                 # can just dump them straight back out without making any changes. 
	    # Pick up parameters - matches the word param separated by spaces
	    if ( ($currline =~ m/!>\s\\param\s+(\S+)\s*/) || ($currline =~ m/!>\s\\param\[.*\]\s+(\S+)\s+/) ) {
		$paramtype = $1; # Contains the param name
		if (exists $params{$paramtype}){ # Covers the case where two arguments have the same name
		    $paramtype=$paramtype . "new";
		}
		$params{$paramtype}=$currline; # Stores the param into a hash indexed by its name
		$matched{$paramtype}=0;          # Set matched for this paramtype to zero as we don't have a match in the argument list - yet
		$param=1;
		$brief=0;
		$par=0;
		$date=0;
		$version=0;
		$author=0;
		$note=0;
		$retval=0;
                $random=0;
		$remainder=0;
	    }
	    elsif ($currline =~ m/!>\s\\brief\s*/) {
		$briefs=$briefs . $currline;
		$param=0;
		$brief=1;
		$par=0;
		$date=0;
		$version=0;
		$author=0;
		$note=0;
		$retval=0;
                $random=0;
		$remainder=0;
	    }
	    elsif ($currline =~ m/!>\s\\date*/) {
		$dates=$dates . $currline;
		$param=0;
		$brief=0;
		$par=0;
		$date=1;
		$version=0;
		$author=0;
		$note=0;
		$retval=0;
                $random=0;
		$remainder=0;
	    }
	    elsif ($currline =~ m/!>\s\\version\s*/) {
		$versions=$versions . $currline;
		$param=0;
		$brief=0;
		$par=0;
		$date=0;
		$version=1;
		$author=0;
		$note=0;
		$retval=0;
                $random=0;
		$remainder=0;
	    }
	    elsif ($currline =~ m/!>\s\\par\s*/) {
		$pars=$pars . $currline;
		$param=0;
		$brief=0;
		$par=1;
		$date=0;
		$version=0;
		$author=0;
		$note=0;
		$retval=0;
                $random=0;
		$remainder=0;
	    }
	    elsif ($currline =~ m/!>\s\\author\s*/) {
		$authors=$authors . $currline;
		$param=0;
		$brief=0;
		$par=0;
		$date=0;
		$version=0;
		$author=1;
		$note=0;
		$retval=0;
                $random=0;
		$remainder=0;
            }
	    elsif ($currline =~ m/!>\s\\note\s*/) {
		$notes=$notes . $currline;
		$param=0;
		$brief=0;
		$par=0;
		$date=0;
		$version=0;
		$author=0;
		$note=1;
                $random=0;
		$remainder=0;
            }
	    elsif ($currline =~ m/!>\s\\retval\s*/) {
		$retvals=$retvals . $currline;
		$param=0;
		$brief=0;
		$par=0;
		$date=0;
		$version=0;
		$author=0;
		$note=0;
		$retval=1;
                $random=0;
		$remainder=0;
            }
            elsif ( ($currline =~ m/!>\s\\\S+/) ) {  # randoms contains anything else that looks like a DOXYGEN header. with a \whatever
		if (($currline !~ m/UNKNOWN_DOXYGEN_COMMENT/) && ($currline !~ m/UNKNOWN_COMMENT/) && ($currline !~ m/!>\s*\\param\s*\n/)  ) { # Must check to see if line has already been commented. We also avoid commenting a blank \param line
		    $randoms=$randoms . $currline;
		    chomp($randoms);
		    $randoms=$randoms . " UNKNOWN_DOXYGEN_COMMENT\n";  # Add on text for output
		} else { 
		    $randoms=$randoms . $currline;  # Otherwise just add to randoms
		}
                $param=0;
                $brief=0;
                $par=0;
                $date=0;
                $version=0;
                $author=0;
                $note=0;
                $retval=0;
                $random=1;
		$remainder=0;
	    }
            # Handle multi line entries. Append what you find onto the previous one read in until you get to 
            # another \
	    # The \brief \param, \par, \author random and remainder entries can all be multi-line. 
	    # The \version and \date entries should be single line and thus we don't have an elsif for them.
	    elsif ($currline =~ m/^!>\s*/){  # Added a * to the \s in case of missing space in the header entry
		if ($param==1){
		    $params{$paramtype}=$params{$paramtype} . $currline; 
		}
		elsif ($brief==1){
		    $briefs=$briefs . $currline;
		}
		elsif ($par==1){
		    $pars=$pars . $currline;
		}
		elsif ($author==1){
		    $authors=$authors . $currline;
		}
		elsif ($note==1){
		    $notes=$notes . $currline;
		}
		elsif ($retval==1){
		    $retvals=$retvals . $currline;
		}
                elsif ($random==1){
		    if (($currline !~ m/UNKNOWN_DOXYGEN_COMMENT/) && ($currline !~ m/UNKNOWN_COMMENT/) ) { # Must check to see if line has already been commented
			$randoms=$randoms . $currline;
			chomp($randoms);
			$randoms=$randoms . " UNKNOWN_DOXYGEN_COMMENT\n"; # Add on text for output
		    } else {
			$randoms=$randoms . $currline; 
		    }
                } else {  # Get any header lines beginning with "!> some text but no \ thus not in DOXYGEN format 
		    $remainder=1;
		    if ( ($currline !~ m/UNKNOWN_COMMENT/) && ($currline !~ m/!>\s*\n/) ) {  # Must check to see if the line has already been commented and also that it's not empty
			$currline =~ s/!>/!> \\note UNKNOWN_COMMENT /g;  # Add on text for output
			$remainders=$remainders . $currline;
		    } else {
			$remainders=$remainders . $currline;
		    }
		}
	    }
            elsif (($currline !~ m/^!\s*\*/) && ($brief==1)) { # If there is more text in the \brief entry that begins with a ! and not !> then make sure this is added on too. 
		$currline =~ s/!/!>/; # Make sure the line starts with !> 
		$briefs=$briefs . $currline;
	    }
            elsif ($currline !~ m/^!\s*\*/) { # Finally pick up anything else in header thats not a *** line e.g. a line beginning "! some comment or other"
                $currline =~ s/^\s+//;
		if ($currline !~ m/UNKNOWN_COMMENT/) {  # Must check to see if line has already been commented
		    $currline =~ s/!/!> \\note UNKNOWN_COMMENT /g; # Add on text for output
		    $remainders=$remainders . $currline;
		} else {
		    $remainders=$remainders . $currline;
		}
            }
	    $currline = <$INPUT>;  # Get the next line in the comment block
	} while ($currline =~ m/^\s*!/i); # Rule as to when you have finished a header. Currently Anything beginning with !

    } # End of if ($currline =~ m/^!>/i) block
# End of header processing

    # Are we inside an interface block - if so we don't want to add any comments here - also allows for blank spaces.
    if ( $currline =~ m/^\s*INTERFACE\s*\n/ ){
        $insideinterface=1;
    } 
    if ( $currline =~ m/^\s*END\s+INTERFACE\s*\n/ ){
        $insideinterface=0;
    }

    if ( ((($currline =~ m/(SUBROUTINE|FUNCTION)\s+(\w|\[|\])+\s*\(/ && $currline !~ m/^\!/ && $currline !~ m/!\s*(SUBROUTINE|FUNCTION)/ ) || $ampersand eq 1) && ($insideinterface eq 0)  ) || ($currline =~ m/^\s*SUBROUTINE\s*(\w|\[|\])+\s*.*/ && $currline !~ m/SUBROUTINE\s+(\w|\[|\])+\s*\(/ ) ) { # We only work on lines that contain subroutines/functions. 
# Assuming that lines contain SUBROUTINE or FUNCTION followed by space and then the name of the procedure which may have a space before the (
# We don't add comments to code inside an interface block
# We also protect against adding comments to commented out SUBROUTINE/FUNCTION calls
# Need to allow for SUBROUTINE without any arguments or brackets too
	if ($ampersand ne 1){
	    $currline =~ /((\bFUNCTION\b|\bSUBROUTINE\b).+)/;  # Remove anything preceeding the SUBROUTINE or FUNCTION definition e.g. RECURSIVE, REAL, INTEGER, (KIND=dp), ELEMENTAL etc etc. output ends up in $1
	    $tmp=$1; # $1 contains whatever remains after removing anything before the word SUBROUTINE or FUNCTION 
	    $tmp =~ s/\bBIND\b\(\w+\W*\w*\W*\w*\W*\)//; # Remove BIND(*) from the definition as we don't want to 
	    $functionline = $tmp; # functionline contains the procedure definition with any unrequired text stripped off, $currline still contains the actual code
	} else { # We also need to check any other lines in subroutine definition for BIND() too, e.g. if the definition extends over more than one line
	    $tmp =$currline;
	    $tmp =~ s/\bBIND\b\(\w+\W*\w*\W*\w*\W*\)//;  # Remove the BIND(*) text as above
	    $functionline = $tmp; # Again functionline is the line we'll work with to find out the procedure arguments etc. 
	}
        chomp($functionline); # Strip the newline char from the end
# Check to see if functionline contains RESULT( 
	if ($functionline =~ m/RESULT\s*\(\w+\)/) {
	    $hasretvalasarg = 1;
	} else {
	    $hasretvalasarg = 0;
	}

        $ampersand = 0;
        $lelement = "";
	
        # Split the subroutine or function definition by space, comma, brackets
        my @string = split /[(),\s+]+/, $functionline;
        foreach(@string) {
	    $p = $_;
	    $matched{$p}=1; # Update the matched hash table so that all arguments in subroutine/function header are set as true
	    # If we encounter an & then it spans multiple lines
	    if ($p eq "&"){
		$ampersand = 1;
		$buffer = "$buffer$currline"; # Buffer the line details as we cant print out yet
	    } elsif ($p !~ m/((\bSUBROUTINE\b|\bFUNCTION\b)|^\z)/i) { 

		if ($lelement eq "FUNCTION" || $lelement eq "SUBROUTINE") { # If previous element is FUNCTION then this must be the name
		    if (($briefs ne "") && ($briefs ne "!> \\brief\n")) { # If brief exists and contains text
			print $OUTPUT $briefs;
		    } elsif ($briefs eq "!> \\brief\n") {  # If \brief is present but empty still need to add text
			print $OUTPUT "! *****************************************************************************\n";
			print $OUTPUT "!> \\brief ...\n";
		    } else { # If brief doesn't exist
			print $OUTPUT "! *****************************************************************************\n";
			print $OUTPUT "!> \\brief ...\n";
		    }
		} else {
		    if ( (exists $params{$p}) && ($params{$p} !~ m/\\param\s*(\w+|\[.*\]\s+\w+)\s*\n/ ) ) { # Entry must exist and contain some text after the parameter
			print $OUTPUT $params{$p};
		    } else {
			if ( ($p ne "RESULT") && ($lelement ne "RESULT") ) {  # Print out all regular procedure arguments, RESULT gets treated as a parameter but we don't want to print it out so protect against this happening. 
                                                                                                        # We don't print the RESULT value under \param as we handle this separately. 
			    # If the entry for this parameter is missing we use the standard text for a missing entry
			    if ($params{$p} eq ""){
				print $OUTPUT "!> \\param $p ...\n";
			    } else { 
                                # We need to guard against \param entries which have no text but have a blank line "!>" line appended on to them this block of code handles this 
				if ( ($params{$p} =~ m/\\param\s*(\w+|\[.*\]\s+\w+)\s*\n/) && ( $params{$p} =~ m/!>\s*\n$/) ) {
				    my @tmpstring = split /\n/, $params{$p};  # Need to split this parameter into it's individual lines
				    $tmplen = scalar @tmpstring;  # Find out how many lines we have in this entry
				    for ($i = 0; $i < $tmplen; $i++){  # Loop over lines
					if ($i eq 0) { # Append on the ... to the first line but leave all the other lines alone
					    $tmpstring[$i] = $tmpstring[$i] . " ..."; 
					}
					print $OUTPUT "$tmpstring[$i]\n";  # Remember to re-add the carriage return to each line
				    }
				} else {
				    chomp($params{$p}); 
				    print $OUTPUT $params{$p}," ...\n";
				}
			    }
			}
			if ( ($lelement eq "RESULT") && ($retvals eq "") && ($hasretvalasarg eq 1) ) {    # If the previous element was RESULT then this element must be whatever gets returned by the procedure - if no retval data is available print out the ... to header
                            #retval is stored for now so it is always printed in the same place (after the \params)
                            $retvals = "!> \\retval $p ...\n";
			}
		    }

		} # End of if $lelement conditional
	    } # End of if $p eq "&" conditional
	    $lelement = $p; # Take a note of the array element for comparison
        } # End of for loop over @string

# Code to loop over all the elements of matched hash table. If any remain that were not in the function/subroutine header 
# we want to make sure we keep this data in the header. It could be a comment or it could refer to old arguments 
# that no longer exist. However, we don't want to throw the text away.  We only print out if the argument doesn't match one 
# of the arguments *and* there is no more of the procedure definition left to read, i.e. ampersand is 0 
	foreach $paramtype ( sort keys %params )   # We need to sort the keys otherwise the output order of the hash is given in internal order 
	{
	    if (($matched{$paramtype} eq 0) && ($ampersand eq 0)) {
		if ($params{$paramtype} eq '!> \param '.$paramtype." ...\n") {
		  # there was no comment, just drop parameter
		} elsif ($params{$paramtype} !~ m/UNMATCHED_PROCEDURE_ARGUMENT/) { # Must protect against updating an existing comment
		    chomp($params{$paramtype}); # Get rid of \n so UNMATCHED* text can be appended on. 
		    print $OUTPUT $params{$paramtype} . " UNMATCHED_PROCEDURE_ARGUMENT: please check \n";
		} else {
		    print $OUTPUT $params{$paramtype};
		}
	    }
	}

        # If after looping through the elements there is no ampersand
        # we close off the comment and write out the procedure definition
        if ($ampersand ne 1) {
            if ($retvals ne ""){   # Print RESULT value first so that it should come straight after the \param definitions
                print $OUTPUT $retvals;
            }           
            if (($dates eq "") || ($dates eq "!> \\date\n")){  # dates entry empty or exists and contains no text
###                print $OUTPUT "!> \\date MISSING_COMMENT: Unknown\n";  # Use this line if you want to add text to the entry
                print $OUTPUT $dates;
            } else {
                print $OUTPUT $dates;
            }
            if (($pars eq "") || ($pars eq "!> \\par History\n")){   # pars entry empty or exists but contains no text 
###                print $OUTPUT "!> \\par History\n"; # Use this line if you want to add text to the entry
###                print $OUTPUT "!>     MISSING_COMMENT: Unknown\n"; # Use this line if you want to add text to the entry
                print $OUTPUT $pars;
            } else {
                print $OUTPUT $pars;
            }
            if (($authors eq "") || ($authors eq "!> \\author\n") ){ # authors empty or exists but contains no text
###                print $OUTPUT "!> \\author MISSING_COMMENT: Unknown\n"; # Use this line if you want to add text to the entry
                print $OUTPUT $authors;
            } else {
                print $OUTPUT $authors;
            }
            if ($versions ne ""){
                print $OUTPUT $versions;
            }           
            if ($randoms ne ""){
                print $OUTPUT $randoms;
            }                
            if ($notes ne ""){
                print $OUTPUT $notes;
            }           
	    if ($remainders ne "") {
		print $OUTPUT $remainders; # Dumps out whatever else remainded in the header (e.g. stuff begining !> without a \ or stuff beginning with just a !) for the SUBROUTINE/FUNCTION at the end
	    }
            print $OUTPUT "! *****************************************************************************\n";
            print $OUTPUT "$leftoverlines$buffer$currline";

# Reset all the variables after writing out the header 
            $buffer = "";
	    $param=0;
	    $brief=0;
	    $par=0;
	    $date=0;
	    $version=0;
	    $author=0;
	    $note=0;
	    $retval=0;
            $random=0;
	    $remainder=0;
	    undef %params;
	    undef %matched;
	    $briefs="";
	    $dates="";
	    $versions="";
	    $pars="";
	    $authors="";
	    $notes="";
	    $retvals="";
            $oldheader="";
            $randoms="";
            $remainders="";
            $leftoverlines="";
        }
    } elsif ( ($oldheader ne "") && ($currline !~ m/^\s*$/) && ($currline !~ m/#if/) ){ # If it was a MODULE or TYPE header we still need to write it back out to file. We also guard against blank lines and #if lines here. 
        if ($currline !~ m/^\s*!\s*/){   # Added ^\s* to prevent lines with inline comments from being moved around
	    print $OUTPUT $oldheader;
	    print $OUTPUT $leftoverlines;
	    print $OUTPUT $currline;
# Again reset variables after output 
	    $param=0;
	    $brief=0;
	    $par=0;
	    $date=0;
	    $version=0;
	    $author=0;
	    $note=0;
	    $retval=0;
            $random=0;
	    $remainder=0;
	    undef %params;
	    undef %matched;
	    $briefs="";
	    $dates="";
	    $versions="";
	    $pars="";
	    $authors="";
	    $notes="";
	    $retvals="";
            $randoms="";
            $remainders="";
            $oldheader="";
            $leftoverlines="";
         } else {
            $oldheader = $oldheader . $currline;
         }
    } else {
        if ($oldheader ne ""){
	    $leftoverlines = $leftoverlines . $currline;  # If we have a procedure header and there's anything else between it and the FUNCTION/SUBROUTINE block we need to store it for printing
        } else {
            print $OUTPUT $currline;   # For all the other lines just write out whatever was read in from file
        }
    } # End of the if (($currline is SUBROUTINE or FUNCTION ) block
    $lastline=$currline; # Keep a copy of last line 
} # End of main while loop

close $OUTPUT;
close $INPUT;

# Perhaps do a second pass to remove any double occurrences of !> *** type lines
