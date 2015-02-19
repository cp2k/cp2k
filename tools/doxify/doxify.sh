#!/bin/bash

DEBUG=0 # Set to 0 for no debugging, 1 for debugging
SCRIPTDIR=$(cd $(dirname "$0"); pwd) # Pick up full path to scripts from wherever doxify.sh lives

for file in "$@"; do

   # generate temp-file names
   tmp_file1=`mktemp`
   tmp_file2=`mktemp`
   tmp_file3=`mktemp`

   # First apply the pre-processing script to get rid of any double & type lines
   $SCRIPTDIR/remove_double_ampersands.pl $file $tmp_file1 $DEBUG

   # Run the fixcomments.pl script. This adds comment blocks to any subroutine/function
   # definitions that don't have any and checks that existing comments are complete,
   # fixing these if required.
   $SCRIPTDIR/fixcomments.pl $tmp_file1 $tmp_file2 $DEBUG

   # After adding comments, remove any double comment header lines
   $SCRIPTDIR/remove_extra_comments.pl $tmp_file2 $tmp_file3 $DEBUG

   # Copy the final modified source file on top of the original file
   if (! cmp -s $file $tmp_file3) ; then
       cp $tmp_file3 $file
   fi

   # Remove temp-files
   rm -f $tmp_file1 $tmp_file2 $tmp_file3

done

#EOF
