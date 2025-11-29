#!/bin/bash

# author: Ole Schuett

if (($# != 2)); then
  echo "Usage: test_build.sh <PROFILE> <VERSION>"
  exit 1
fi

VERSION=$2

if [ -f build/bin/cp2k."${VERSION}" ]; then
  echo -e "\nSummary: Compilation works fine."
  echo -e "Status: OK\n"
else
  echo -e "\nSummary: Compilation failed."
  echo -e "Status: FAILED\n"
fi

exit 0 # Prevent CI from overwriting our summary message.

#EOF
