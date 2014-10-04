#!/bin/bash -e

# Test if coding conventions are fulfilled
# author: Ole Schuett

echo -n "Date: "
date --utc --rfc-3339=seconds

#svn revert -R .
svn up
svn info

cd makefiles
make -j ARCH=Linux-x86-64-gfortran VERSION=dumpast

cd ../obj/Linux-x86-64-gfortran/dumpast
echo "Checking conventions..."
../../../tools/analyze_gfortran_ast.py --suppressions=../../../tools/dashboard/conventions.supp *.ast

#EOF
