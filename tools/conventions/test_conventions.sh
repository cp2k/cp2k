#!/bin/bash -e

# author: Ole Schuett

./tools/conventions/analyze_gfortran_ast.py ./build/*.ast > ast.issues
./tools/conventions/analyze_gfortran_warnings.py ./build/*.warn > warn.issues
./tools/conventions/summarize_issues.py --suppressions=./tools/conventions/conventions.supp ./*.issues

# EOF
