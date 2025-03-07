#!/bin/bash

# author: Ole Schuett

if (($# != 1)); then
  echo "Usage: clang_format_wrapper.sh <file>"
  exit 1
fi

# For some files the clang-format requires too much memory. As a work-around we
# use the "clang-format off" marker to disable formatting. To actually reduce
# the memory usage this script hides the offending sections from clang-format.

# Look for "clang-format off" markers to split the file into sections.
if csplit --quiet --prefix="section" "$1" "/clang-format off/"; then
  # Format only the first section, ie. the one before "clang-format off".
  clang-format --style=llvm -i section00
  # Join all sections back together.
  cat section* > "$1"
  # Remove intermediary files.
  rm section*
else
  # No markers found, just call clang-format directly.
  clang-format --style=llvm -i "$1"
fi

#EOF
