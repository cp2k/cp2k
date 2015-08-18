#!/bin/bash -e

# Test if coding conventions are fulfilled
# author: Ole Schuett

echo -n "Date: "
date --utc --rfc-3339=seconds


#svn up
#svn info

cd ../../makefiles
make -j ARCH=Linux-x86-64-gfortran VERSION="dumpast" > make_conventions1.out
make -j ARCH=local_cuda_warn       VERSION="psmp"    > make_conventions2.out

cd ../tools/conventions

rm -f *.issues
./analyze_gfortran_ast.py      ../../obj/Linux-x86-64-gfortran/dumpast/*.ast > ast.issues
./analyze_gfortran_warnings.py ../../obj/local_cuda_warn/psmp/*.warn         > warn.issues
./analyze_src.py               ../../                                        > src.issues

./summarize_issues.py --suppressions=conventions.supp *.issues

#EOF
