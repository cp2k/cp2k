#!/bin/bash -e

# author: Ole Schuett

echo -e "\n========== Generating Doxygen =========="
cd /opt/cp2k
if ! ./tools/doxify/generate_doxygen.sh &> doxygen.out; then
  echo -e "failed.\n\n"
  tail -n 100 doxygen.out
  mkdir -p /workspace/artifacts/
  cp make.out /workspace/artifacts/
  echo -e "\nSummary: Doxygen generation failed."
  echo -e "Status: FAILED\n"
  exit 0
fi

mkdir -p /workspace/artifacts
mv ./doxygen/html /workspace/artifacts/doxygen

echo "Summary: Doxygen generation works fine."
echo "Status: OK"

#EOF
