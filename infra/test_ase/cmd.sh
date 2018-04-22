#!/bin/bash

# author: Ole Schuett

set -e

echo -e "\n========== Rsyncing =========="
rsync --exclude="*~"          \
      --exclude=".*/"         \
      --exclude="*.pyc"       \
      --exclude=/cp2k/obj/    \
      --exclude=/cp2k/lib/    \
      --exclude=/cp2k/exe/    \
      --ignore-times          \
      --update                \
      --verbose               \
      --recursive             \
      --checksum              \
      /opt/cp2k-local/  /opt/cp2k-master/

echo -e "\n========== Compiling CP2K =========="
source /opt/cp2k-toolchain/install/setup
cd /opt/cp2k-master/cp2k/makefiles
make -j VERSION=pdbg cp2k_shell
ln -s /opt/cp2k-master/cp2k/exe/local/cp2k.pdbg /usr/bin/cp2k

echo -e "\n========== Installing ASE =========="
cd /opt/
git clone https://gitlab.com/ase/ase.git
pip3 install ./ase/

echo -e "\n========== Running ASE tests =========="
cd /opt/ase/
export ASE_CP2K_COMMAND="mpiexec -np 2 /opt/cp2k-master/cp2k/exe/local/cp2k_shell.pdbg"

(
set -e # abort if error is encountered
for i in ./ase/test/cp2k/cp2k_*.py
do
  echo "Running $i ..."
  python3 $i
done
)

EXIT_CODE=$?

echo ""

ASE_REVISION=`git rev-parse --short HEAD`
if (( $EXIT_CODE )); then
    echo "Summary: Something is wrong with ASE commit ${ASE_REVISION}."
    echo "Status: FAILED"
else
    echo "Summary: ASE commit ${ASE_REVISION} works fine."
    echo "Status: OK"
fi

#EOF
