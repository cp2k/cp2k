#!/bin/bash

export DEBUG=0 # Set to 0 for no debugging, 1 for debugging
SCRIPTDIR=$(cd $(dirname "$0"); pwd) # Pick up full path to scripts from wherever doxyify.sh lives

for file in "$@"
do 

# First apply the pre-processing script to get rid of any double & type lines
    if [ $DEBUG == 1 ] 
    then
	$SCRIPTDIR/remove_double_ampersands.pl $file $file.preprocessed 1 
    else 
	$SCRIPTDIR/remove_double_ampersands.pl $file $file.preprocessed 
    fi

# Run the fixcomments.pl script. This adds comment blocks to any subroutine/function 
# definitions that don't have any and checks that existing comments are complete, 
# fixing these if required. 
    if [ $DEBUG == 1 ] 
    then
	$SCRIPTDIR/fixcomments.pl $file.preprocessed $file.tmp 1
    else
	$SCRIPTDIR/fixcomments.pl $file.preprocessed $file.tmp
    fi

# After adding comments, remove any double comment header lines
    if [ $DEBUG == 1 ] 
    then
	$SCRIPTDIR/remove_extra_comments.pl $file.tmp $file.new 1
    else
	$SCRIPTDIR/remove_extra_comments.pl $file.tmp $file.new
    fi

# DEBUG: dump difference of original and final commented file into diffs.txt so we 
# can look at what was changed. Used for debugging and can probably be omitted for 
# the final version
# Also runs tkdiff. You only want to run tkdiff on a small number of files
# Change the path to location of tkdiff on your system. 
    if [ $DEBUG == 1 ] 
    then 
	echo "Diff for file $file" >> diffs.txt
	diff $file.preprocessed $file.new >> diffs.txt
	echo " " >> $diffs.txt
	cp $file $file.orig  # copy original input to *.orig as and tkdiff this
	/Applications/TkDiff.app/Contents/MacOS/tkdiff $file.orig $file &
    fi

# Copy the final modified source file on top of the original file
    cp $file.new $file

# Remove any uneeded output files 
    rm -f $file.tmp $file.preprocessed $file.new 
    if [ $DEBUG == 0 ] # We only remove .orig file if Tkdiff is not used
    then 
	rm -f $file.orig
    fi

done 
