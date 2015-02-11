#!/bin/bash -e

# Test if coding conventions are fulfilled
# author: Ole Schuett

echo -n "Date: "
date --utc --rfc-3339=seconds


#svn up
#svn info

cd ../../makefiles
make -j ARCH=Linux-x86-64-gfortran VERSION="dumpast warn" > make_conventions.out

cd ../tools/conventions

rm -f *.issues
./analyze_gfortran_ast.py      ../../obj/Linux-x86-64-gfortran/dumpast/*.ast > ast.issues
./analyze_gfortran_warnings.py ../../obj/Linux-x86-64-gfortran/warn/*.warn   > warn.issues
./analyze_src.py              ../../                                         > src.issues

./summarize_issues.py --suppressions=conventions.supp *.issues

#EOF
