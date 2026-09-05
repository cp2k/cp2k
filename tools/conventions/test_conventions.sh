#!/bin/bash -e

# author: Ole Schuett

JOBS=0
while [ "$#" -gt 0 ]; do
  case "$1" in
    -j | --jobs)
      JOBS="$2"
      shift 2
      ;;
    -j[0-9]*)
      JOBS="${1#-j}"
      shift
      ;;
    *)
      echo "Unknown option: $1" >&2
      exit 2
      ;;
  esac
done

if [[ $# -ne 0 || ! "${JOBS}" =~ ^[0-9]+$ ]]; then
  echo "Usage: $0 [-j jobs]" >&2
  exit 2
fi

./tools/conventions/analyze_gfortran_ast.py --jobs "${JOBS}" ./build/*.ast > ast.issues
./tools/conventions/analyze_gfortran_warnings.py --jobs "${JOBS}" ./build/*.warn > warn.issues
./tools/conventions/summarize_issues.py --suppressions=./tools/conventions/conventions.supp ./*.issues

# EOF
