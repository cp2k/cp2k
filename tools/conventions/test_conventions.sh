#!/bin/bash

# Test if coding conventions are fulfilled
# author: Ole Schuett

echo -n "Date: "
date --utc --rfc-3339=seconds

(
  set -e # abort if error is encountered
  cd ../../makefiles
  make -j ARCH=Linux-x86-64-gfortran VERSION="dumpast" > make_conventions1.out
  make -j ARCH=local_warn            VERSION="psmp"    > make_conventions2.out
)
MAKE_EXIT_CODE=$?

if (( $MAKE_EXIT_CODE )); then
   cd ../../obj/local_warn/psmp/
   echo ""
   grep -B 2 Error *.warn
   echo ""
   echo "Summary: Compilation failed."
   echo "Status: FAILED"
else
   rm -f *.issues
   ./analyze_gfortran_ast.py      ../../obj/Linux-x86-64-gfortran/dumpast/*.ast > ast.issues
   ./analyze_gfortran_warnings.py ../../obj/local_warn/psmp/*.warn              > warn.issues
   ./analyze_src.py               ../../                                        > src.issues
   ./summarize_issues.py --suppressions=conventions.supp *.issues
fi
#EOF
