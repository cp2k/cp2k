#!/bin/bash

set -o nounset
set -o pipefail

# Test if coding conventions are fulfilled
# author: Ole Schuett

echo -n "Date: "
date --utc --rfc-3339=seconds

(
  set -e # abort if error is encountered
  make -j ARCH=Linux-x86-64-gfortran VERSION="dumpast" > make_conventions1.out
  make -j ARCH=local_warn VERSION="psmp" > make_conventions2.out
)
MAKE_EXIT_CODE=$?

if ((MAKE_EXIT_CODE)); then
  echo ""
  grep -B 2 Error ./obj/local_warn/psmp/*.warn
  echo ""
  echo "Summary: Compilation failed."
  echo "Status: FAILED"
else
  rm -f ./*.issues

  set -o errexit

  ./tools/conventions/analyze_gfortran_ast.py ./obj/Linux-x86-64-gfortran/dumpast/*.ast > ast.issues
  ./tools/conventions/analyze_gfortran_warnings.py ./obj/local_warn/psmp/*.warn > warn.issues
  ./tools/conventions/summarize_issues.py --suppressions=./tools/conventions/conventions.supp ./*.issues
fi
