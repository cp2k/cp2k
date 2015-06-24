#!/usr/bin/env perl

=head1 NAME

fixcomments.pl

=head1 SYNOPSIS

fixcomments.pl infile outfile [verbosity]

=head1 OPTIONS

=over 4

=item B<infile>

Input file to process.

=item B<outfile>

Output file to write after processing.

=item B<verbosity>

Optional verbosity level: 0=quiet, 1=some info, 2=lots of info

=back

=head1 DESCRIPTION

This script will read a Fortran file and attempt to add doxify comments
to both FUNCTIONs and SUBROUTINEs.

=cut

use strict;
use warnings;
use FindBin;
use lib "$FindBin::RealBin/lib";
use Pod::Usage qw(pod2usage);

# Check options
if ((@ARGV != 2) && (@ARGV != 3)) {
    pod2usage(-verbose => 1, -message => "$0: Wrong number of arguments given.\n");
}

my $DBG=0;
if (defined($ARGV[2])) {
    $DBG = $ARGV[2];
}

sub print_debug{
    if($DBG>0){ print "DEBUG: @_\n";}
}

print_debug("fixcomments.pl running with");
print_debug("Input:  $ARGV[0]");
print_debug("Output: $ARGV[1]");

# Regular expressions for matching
my $DOXYGEN_HEADER = "^!>";

# The empty string
my ($EMPTY);
$EMPTY=q{};

# Toggle variables to keep track of which doxygen item is being processed
# in the current header.
my ($hasParam,
    $hasBrief,
    $hasPar,
    $hasDate,
    $hasVersion,
    $hasAuthor,
    $hasNote,
    $hasRetVal,
    $hasRandom,
    $hasRemainder);

# Variables with s at the end contains the actual text read in from
# existing doxygen header
my (%params,
    $briefs,
    $dates,
    $versions,
    $pars,
    $authors,
    $notes,
    $retVals,
    $randoms,
    $remainders);

# If currently in INTERFACE block
my $insideInterface = 0;

# Whether the SUBROUTINE/FUNCTION definition contains an & line continuation
my $hasAmpersand = 0;

# Whether when processing subroutine definition we are in brackets
my $inParentheses = 0;

# Whether the procedure definition contains the RETURN value
my $hasRetValAsArg = 0;

# Keeps track of which parameters in the current procedure have
# a matching doxygen header
my %matched;

# Buffer for lines
my $buffer = $EMPTY;

# The old doxygen header as read from file
my $oldheader = $EMPTY;

# Any unprocessed lines
my $leftoverlines = $EMPTY;

# The name of the procedure being processed
my $procedureName = $EMPTY;

# Whether the procedure is a function
my $isFunction;

initVariables();

open (my $INPUT , "<" , $ARGV[0]) or die("Cant open $ARGV[0] $!");
open (my $OUTPUT , ">" , $ARGV[1]) or die("Cant create $ARGV[1] $!");

# While there are still lines to read in our INPUT file
while (<$INPUT>) {
    # Get the line we've just read in
    my $currline = $_;

    # If an existing doxygen header is encountered then process
    if (matchDoxygenHeaderDefinition($currline)) {
        print_debug("Matched doxygen header");
        $currline = processDoxygenHeader($currline);
    }

    # Are we inside an interface block?
    # If so we don't want to add any comments here
    if ( $currline =~ m/^\s*INTERFACE\s*\n/xms ) {
        print_debug("Start of INTERFACE encountered");
        $insideInterface = 1;
    }
    if ( $currline =~ m/^\s*END\s+INTERFACE\s*\n/xms ) {
        print_debug("End of INTERFACE encountered");
        $insideInterface = 0;
    }

    # Look for procedure (SUBROUTINE/FUNCTION definitions)
    # These can be the initial definiton line, or a continuation line
    # We don't add comments to code inside an interface block
    if ((matchSubroutineDefinition($currline) || $hasAmpersand)
        && !($insideInterface)) {
        print_debug("Matched subroutine definition");
        processSubroutineDefinition($currline);
    } else {
        print_debug("Subroutine definition not matched");
        if ($oldheader eq $EMPTY) {
            # No header remaining so just print out line as read
            print_debug("Empty old header, writing out line:");
            print_debug($currline);
            print $OUTPUT $currline;
        } else {
            # Header has been processed, need to do something with remaining lines
            print_debug("Non-empty old header for line:");
            print_debug($currline);
            if (($currline =~ m/^\s*$/xms) ||
                ($currline =~ m/\#if/xms)) {
                # Blank line or #if line, so store for printing later
                print_debug("Empty line or #if, storing for later");
                $leftoverlines = $leftoverlines . $currline;
            } elsif ($currline =~ m/^\s*!\s*/xms) {
                # Have a comment so add to header
                print_debug("Inline comment, adding to header");
                $oldheader = $oldheader . $currline;
            } else {
                # If it was a MODULE or TYPE header we still need to write
                # it back out to file.
                print_debug("Writing existing non FUNCTION/SUBROUTINE header to file");
                print $OUTPUT $oldheader;
                print $OUTPUT $leftoverlines;
                print $OUTPUT $currline;
                # Reset variables after output
                initVariables();
            }
        }
    } # End of the if (($currline is SUBROUTINE or FUNCTION ) block
} # End of main while loop

close $OUTPUT;
close $INPUT;

# Perhaps do a second pass to remove any double occurrences of !> *** type lines

sub initToggles {
    $hasParam = 0;
    $hasBrief = 0;
    $hasPar = 0;
    $hasDate = 0;
    $hasVersion = 0;
    $hasAuthor = 0;
    $hasNote = 0;
    $hasRetVal = 0;
    $hasRandom = 0;
    $hasRemainder = 0;
    return;
}

sub initVariables {
    initToggles();
    undef %params;
    undef %matched;
    $briefs = $EMPTY;
    $dates = $EMPTY;
    $versions = $EMPTY;
    $pars = $EMPTY;
    $authors = $EMPTY;
    $notes = $EMPTY;
    $retVals = $EMPTY;
    $randoms = $EMPTY;
    $remainders = $EMPTY;
    $buffer = $EMPTY;
    $oldheader = $EMPTY;
    $leftoverlines = $EMPTY;
    $procedureName = $EMPTY;
    $isFunction = 0;
    $inParentheses = 0;
    return;
}

# Look for a doxygen header definition in an input line.
sub matchDoxygenHeaderDefinition {
    my ($lineToProcess) = @_;

    print_debug("Trying to match doxygen header $DOXYGEN_HEADER against line");
    print_debug($lineToProcess);

    my $match = 0;

    if ($lineToProcess =~ m/$DOXYGEN_HEADER/xms) {
        $match = 1;
    }

    print_debug("Doxygen header match value: $match");
    return $match;
}

# Look for SUBROUTINE or FUNCTION definitions.
sub matchSubroutineDefinition {
    my ($lineToProcess) = @_;

    print_debug("Trying to match subroutine definition against line:");
    print_debug($lineToProcess);

    # Immediately discount lines with comments at the start
    my $commentPattern = '^\!';
    if ($lineToProcess =~ m/$commentPattern/xms) {
        print_debug("Matched comment");
        return 0;
    }
    # Assume that lines contain SUBROUTINE or FUNCTION followed by space
    # and then the name of the procedure which may have a space before the (
    my $patternProcWithBrackets = '(SUBROUTINE|FUNCTION)'    # Subroutine or function
                                 .'\s+'                      # followed by one or more whitespace characters
                                 .'(\w|\[|\])+'              # followed by one or more of (word or [ or ])
                                 .'\s*'                      # followed by zero or more whitespace characters
                                 .'\(';                      # followed by open bracket.

    # Need to allow for SUBROUTINE without any arguments or brackets too
    my $patternSubNoBracketsAtStartOfLine = '^\s*'           # Start with zero or more whitespace
                                           .'SUBROUTINE'     # then subroutine
                                           .'\s*'            # then zero or more whitespace
                                           .'(\w|\[|\])+'    # then one or more of (word or [ or ])
                                           .'\s*.*';         # then zero or more whitespace

    # We also protect against adding comments to commented out SUBROUTINE/FUNCTION calls
    my $patternCommentedOut = '!\s*(SUBROUTINE|FUNCTION)';

    # Protect against definitions inside quotes
    my $patternInQuotes = '"\s*(SUBROUTINE|FUNCTION)';

    my $match1 = 0;
    my $match2 = 0;
    my $match3 = 0;
    my $match4 = 0;

    if ($lineToProcess =~ m/$patternProcWithBrackets/xms) {
        $match1 = 1;
    }
    if ($lineToProcess =~ m/$patternSubNoBracketsAtStartOfLine/xms) {
        $match2 = 1;
    }
    if ($lineToProcess =~ m/$patternCommentedOut/xms) {
        $match3 = 1;
    }
    if ($lineToProcess =~ m/$patternInQuotes/xms) {
        $match4 = 1;
    }

    my $match = (($match1 || $match2) && !($match3) && !($match4));
    print_debug("Subroutine definition match value: $match");
    return $match
}

sub processDoxygenHeader {

    my ($currline) = @_;

    print_debug("Processing Doxygen Header");

    # Each time we find a header we need to reset toggles and their data
    initVariables();
    my $paramName = $EMPTY;

    # Start of do-while over match on ($currline =~ m/^\s*!/i)
    do {
        print_debug("Processing header line:");
        print_debug($currline);
        # Keep the headers safe, may need them! We use the oldheader variable to
        # keep the complete headers for MODULE and TYPE intact such that we
        # can just dump them straight back out without making any changes.
        $oldheader = $oldheader . $currline;
        # Pick up parameters - matches the word param separated by spaces
        if (($currline =~ m/!>\s\\param\s+(\S+)\s*/xms) ||
            ($currline =~ m/!>\s\\param\[.*\]\s+(\S+)\s*/xms)) {

            # Get the param name
            $paramName = $1;

            print_debug("Got header for parameter $paramName");

            # Cover the case where two arguments have the same name
            if (exists $params{$paramName}) {
                $paramName = $paramName . "new";
            }
            # Store the param into a hash indexed by its name
            $params{$paramName} = $currline;
            # Set matched for this param to zero as we don't
            # yet have a match in the argument list
            $matched{$paramName} = 0;
            initToggles();
            $hasParam = 1;
        } elsif ($currline =~ m/!>\s\\brief\s*/xms) {
            print_debug("Got briefs header");
            $briefs = $briefs . $currline;
            initToggles();
            $hasBrief = 1;
        } elsif ($currline =~ m/!>\s\\date*/xms) {
            print_debug("Got dates header");
            $dates = $dates . $currline;
            initToggles();
            $hasDate = 1;
        } elsif ($currline =~ m/!>\s\\version\s*/xms) {
            print_debug("Got version header");
            $versions = $versions . $currline;
            initToggles();
            $hasVersion = 1;
        } elsif ($currline =~ m/!>\s\\par\s*/xms) {
            print_debug("Got par header");
            $pars = $pars . $currline;
            initToggles();
            $hasPar = 1;
        } elsif ($currline =~ m/!>\s\\author\s*/xms) {
            print_debug("Got author header");
            $authors = $authors . $currline;
            initToggles();
            $hasAuthor= 1;
        } elsif ($currline =~ m/!>\s\\note\s*/xms) {
            print_debug("Got note header");
            $notes = $notes . $currline;
            initToggles();
            $hasNote = 1;
        } elsif ($currline =~ m/!>\s\\retval\s*/xms) {
            print_debug("Got retval header");
            $retVals = $retVals . $currline;
            initToggles();
            $hasRetVal = 1;
        } elsif ($currline =~ m/!>\s\\\S+/xms) {
            # Randoms contains anything else that looks
            # like a DOXYGEN header. with a \whatever
            # Must check to see if line has already been commented.
            # We also avoid commenting a blank \param line
            print_debug("Got random header");
            if (($currline !~ m/UNKNOWN_DOXYGEN_COMMENT/xms) &&
                ($currline !~ m/UNKNOWN_COMMENT/xms) &&
                ($currline !~ m/!>\s*\\param\s*\n/xms)) {
                $randoms = $randoms . $currline;
                chomp($randoms);
                # Add on text for output
                $randoms = $randoms . " UNKNOWN_DOXYGEN_COMMENT\n";
            } else {
                # Otherwise just add to randoms
                $randoms = $randoms . $currline;
            }
            initToggles();
            $hasRandom = 1;
        } elsif ($currline =~ m/^!>\s*/xms) {
            # Handle multi line entries. Append what you find onto the
            # previous one read in until you get to another \
            # The \brief, \param, \par, \author, random and remainder
            # entries can all be multi-line.
            # The \version and \date entries should be single line and
            # thus we don't have an elsif for them.
            print_debug("Got line continuation");
            if ($hasParam) {
                $params{$paramName} = $params{$paramName} . $currline;
            } elsif ($hasBrief) {
                $briefs = $briefs . $currline;
            } elsif ($hasPar) {
                $pars = $pars . $currline;
            } elsif ($hasAuthor) {
                $authors = $authors . $currline;
            } elsif ($hasNote){
                $notes = $notes . $currline;
            } elsif ($hasRetVal) {
                $retVals = $retVals . $currline;
            } elsif ($hasRandom) {
                # Must check to see if line has already been commented
                if (($currline !~ m/UNKNOWN_DOXYGEN_COMMENT/xms) &&
                    ($currline !~ m/UNKNOWN_COMMENT/xms)) {
                    $randoms = $randoms . $currline;
                    chomp($randoms);
                    # Add on text for output
                    $randoms = $randoms . " UNKNOWN_DOXYGEN_COMMENT\n";
                } else {
                    $randoms = $randoms . $currline;
                }
            } else {
                # Get any header lines beginning with "!> some text but
                # no \ thus not in DOXYGEN format
                $hasRemainder = 1;
                # Must check to see if the line has already been
                # commented and also that it's not empty
                if (($currline !~ m/UNKNOWN_COMMENT/xms) &&
                    ($currline !~ m/!>\s*\n/xms)) {
                    # Add on text for output
                    $currline =~ s/!>/!> \\note UNKNOWN_COMMENT /gx;
                    $remainders = $remainders . $currline;
                } else {
                    $remainders = $remainders . $currline;
                }
            }
        } elsif ($currline !~ m/^!\s*\*/xms) {
            # Any other header line that's not "***..." or a Doxygen
            # comment, ie "! some comment or other"
            print_debug("Got non-doxygen line");
            if ($hasBrief) {
                # Put comment in brief, replacing ! with !>
                $currline =~ s/!/!>/xms;
                $briefs = $briefs . $currline;
            } else {
                $currline =~ s/^\s+//;
                # Must check to see if line has already been commented
                if ($currline !~ m/UNKNOWN_COMMENT/xms) {
                    # Add on text for output, changing ! for !>
                    $currline =~ s/!/!> \\note UNKNOWN_COMMENT /gx;
                    $remainders = $remainders . $currline;
                } else {
                    $remainders = $remainders . $currline;
                }
            }
        }
        # Get the next line in the header block
        $currline = <$INPUT>;
    } while ($currline =~ m/^\s*!/xms); # Rule as to when you have finished a header. Currently Anything beginning with !
    return $currline;
} # End of header processing subroutine

sub processSubroutineDefinition {

    my ($currline) = @_;

    print_debug("Processing Subroutine line:");
    print_debug($currline);
    # functionLine will contain the procedure definition with any
    # unrequired text stripped off, $currline will still contain the actual code
    my $functionLine = $EMPTY;
    if (!($hasAmpersand)) {
        # Remove anything preceeding the SUBROUTINE or FUNCTION definition
        # e.g. RECURSIVE, REAL, INTEGER, (KIND=dp), ELEMENTAL etc etc.
        $currline =~ /((\bFUNCTION\b|\bSUBROUTINE\b).+)/x;
        # $1 contains whatever remains after removing anything before the
        # word SUBROUTINE or FUNCTION
        $functionLine = $1;
    } else {
        $functionLine = $currline;
    }
    # Remove BIND(*) from the definition
    my $tmp = $functionLine;
    $tmp =~ s/\bBIND\b\(\w+\W*\w*\W*\w*\W*\)//x;
    $functionLine = $tmp;

    # Strip the newline char from the end
    chomp($functionLine);

    $hasAmpersand = 0;
    $hasRetValAsArg = 0;
    my $lelement = $EMPTY;

    # Split the subroutine or function definition by space or comma
    my @string = split(/([,\(\)\s+]+)/x, $functionLine);
    foreach (@string) {
        if (defined) {
        my $p = $_;
        $p =~ s/^\s+|\s+$//gx;
        if (($p ne $EMPTY) && ($p ne ",")) {
        print_debug("Processing item $p");
        if ($p eq "&") {
            # If we encounter an & then it spans multiple lines
            print_debug("Got & line continuation");
            $hasAmpersand = 1;
            # Buffer the line details as we can't print out yet
            $buffer = "$buffer$currline";
        } elsif ($p eq "(") {
            print_debug("Entered parentheses");
            $inParentheses = 1;
        } elsif ($p eq ")") {
            print_debug("Left parentheses");
            $inParentheses = 0;
        } elsif ($p eq ",") {
            print_debug("Comma");
        } else {
            if ($inParentheses) {
                print_debug("In parentheses, parameter: $p, previous parameter: $lelement");
                # Must have either a parameter of a function return value
                if (($lelement =~ m/RESULT/ixms) && ($hasRetValAsArg)) {
                    print_debug("Got return parameter $p");
                    if ($retVals eq $EMPTY) {
                        # If the previous element was RESULT outside of
                        # parantheses then this element must be whatever
                        # gets returned by the procedure.
                        # if no retval data is available print out the ...
                        # to header
                        # Only stored for now so it is always printed
                        # in the same place (after the \params)
                        $retVals = "!> \\retval $p ...\n";
                    }
                } else {
                    # Must be parameter
                    # Update the matched hash table so that all arguments
                    # in subroutine/function header are set as true
                    print_debug("Processing parameter $p");
                    $matched{$p} = 1;
                    if (!(exists($params{$p})) || !(defined($params{$p})) || ($params{$p} eq $EMPTY)) {
                        # If the entry for this parameter is missing we use
                        # the standard text for a missing entry
                        print_debug("Missing entry for parameter $p, creating blank");
                        print $OUTPUT "!> \\param $p ...\n";
                    } else {
                        print_debug("Using existing entry for parameter $p");
                        if ($params{$p} !~ m/\\param\s*(\w+|\[.*\]\s+\w+)\s*\n/xms ) {
                            # Entry must contain some text after the parameter name
                            print_debug("Using entry unchanged");
                            print $OUTPUT $params{$p};
                        } else {
                            if ($params{$p} =~ m/!>\s*\n$/x) {
                                # We need to guard against \param entries which have
                                # no text but have a blank line "!>" line appended
                                # Need to split this parameter into it's individual lines
                                my @tmpString = split(/\n/x, $params{$p});
                                # Find out how many lines we have in this entry
                                my $tmpLen = scalar(@tmpString);
                                # Append on the ... to the first line
                                $tmpString[0] = $tmpString[0] . " ...";
                                for (my $i = 0; $i < $tmpLen; $i++) {
                                    # Re-add the carriage return to each line
                                    print $OUTPUT "$tmpString[$i]\n";
                                }
                            } else {
                                chomp($params{$p});
                                print $OUTPUT $params{$p}," ...\n";
                            }
                        }
                    }
                }
            } else {
                print_debug("Processing element not in parentheses: $p");
                if ($lelement =~ m/(SUBROUTINE|FUNCTION)/ixms) {
                    # Previous element is FUNCTION or SUBROUTINE so this must be name
                    $procedureName = $p;
                    print_debug("Got procedure name $procedureName");
                    if ($lelement =~ m/FUNCTION/ixms) {
                        $isFunction = 1;
                    }
                    if (($briefs eq $EMPTY) or ($briefs eq "!> \\brief\n")) {
                        # \brief does not exist or is present but empty - add text
                        print $OUTPUT "! *****************************************************************************\n";
                        print $OUTPUT "!> \\brief ...\n";
                    } else {
                        # \brief exists and contains text
                        print $OUTPUT $briefs;
                    }
                } elsif ($p =~ m/RESULT/ixms) {
                    # Check to see if parameter is RESULT for a FUNCTION
                    print_debug("Need to get function result parameter");
                    $hasRetValAsArg = 1;
                }
            }
            $lelement = $p; # Take a note of the array element for comparison, ignoring & and ()
        } # End of if $p eq "&" conditional
    }
        }
    } # End of for loop over @string

    # Code to loop over all the elements of matched hash table. If any
    # remain that were not in the function/subroutine header we want to
    # make sure we keep this data in the header. It could be a comment
    # or it could refer to old arguments that no longer exist.
    # However, we don't want to throw the text away. We only print out
    # if the argument doesn't match one of the arguments *and* there is
    # no more of the procedure definition left to read, i.e. hasAmpersand is false
    foreach my $paramName (sort keys %params) {
        # We need to sort the keys otherwise the output order of the hash
        # is given in internal order
        if (!($matched{$paramName}) && !($hasAmpersand)) {
            if ($params{$paramName} eq '!> \param '.$paramName." ...\n") {
                # there was no comment, just drop parameter
            } elsif ($params{$paramName} !~ m/UNMATCHED_PROCEDURE_ARGUMENT/) {
                # Must protect against updating an existing comment
                # Get rid of \n so UNMATCHED* text can be appended on.
                chomp($params{$paramName});
                print $OUTPUT $params{$paramName} . " UNMATCHED_PROCEDURE_ARGUMENT: please check \n";
            } else {
                print $OUTPUT $params{$paramName};
            }
        }
    }

    # If after looping through the elements there is no ampersand
    # we close off the comment and write out the procedure definition
    if (!($hasAmpersand)) {
        if ($retVals ne $EMPTY) {
            # Print RESULT value first so that it should come straight after the \param definitions
            print $OUTPUT $retVals;
        } else {
            # Get return value from function name
            if ($isFunction) {
                print $OUTPUT "!> \\retval $procedureName ...\n";
            }
        }
        if (($dates eq $EMPTY) || ($dates eq "!> \\date\n")) {
            # dates entry empty or exists and contains no text
###            print $OUTPUT "!> \\date MISSING_COMMENT: Unknown\n";  # Use this line if you want to add text to the entry
            print $OUTPUT $dates;
        } else {
            print $OUTPUT $dates;
        }
        if (($pars eq $EMPTY) || ($pars eq "!> \\par History\n")) {
            # pars entry empty or exists but contains no text
###            print $OUTPUT "!> \\par History\n"; # Use this line if you want to add text to the entry
###            print $OUTPUT "!>     MISSING_COMMENT: Unknown\n"; # Use this line if you want to add text to the entry
            print $OUTPUT $pars;
        } else {
            print $OUTPUT $pars;
        }
        if (($authors eq $EMPTY) || ($authors eq "!> \\author\n")) {
            # authors empty or exists but contains no text
###            print $OUTPUT "!> \\author MISSING_COMMENT: Unknown\n"; # Use this line if you want to add text to the entry
            print $OUTPUT $authors;
        } else {
            print $OUTPUT $authors;
        }
        if ($versions ne $EMPTY) {
            print $OUTPUT $versions;
        }
        if ($randoms ne $EMPTY) {
            print $OUTPUT $randoms;
        }
        if ($notes ne $EMPTY){
            print $OUTPUT $notes;
        }
        if ($remainders ne $EMPTY) {
            # Dumps out whatever else remainded in the header (e.g. stuff begining !> without a \ or stuff beginning with just a !) for the SUBROUTINE/FUNCTION at the end
            print $OUTPUT $remainders;
        }
        print $OUTPUT "! *****************************************************************************\n";
        print $OUTPUT "$leftoverlines$buffer$currline";
        # Reset all the variables after writing out the header
        initVariables();
        return;
    }
}
