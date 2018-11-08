#!/bin/bash -e

# author: Ole Schuett

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

echo -e "\n========== Compiling CP2K =========="
cd /workspace/cp2k
make -j VERSION=pdbg cp2k_shell

echo -e "\n========== Installing ASE =========="
cd /opt/ase/
git pull
pip3 install .

echo -e "\n========== Running ASE Tests =========="
cd /opt/ase/
export ASE_CP2K_COMMAND="mpiexec -np 2 /workspace/cp2k/exe/local/cp2k_shell.pdbg"

set +e # disable error trapping for remainder of script
(
set -e # abort if error is encountered
for i in ./ase/test/cp2k/cp2k_*.py
do
  echo "Running $i ..."
  python3 "$i"
done
)

EXIT_CODE=$?

echo ""

ASE_REVISION=$(git rev-parse --short HEAD)
if (( EXIT_CODE )); then
    echo "Summary: Something is wrong with ASE commit ${ASE_REVISION}."
    echo "Status: FAILED"
else
    echo "Summary: ASE commit ${ASE_REVISION} works fine."
    echo "Status: OK"
fi

#EOF
