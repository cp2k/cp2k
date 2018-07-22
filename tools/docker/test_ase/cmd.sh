#!/bin/bash

# author: Ole Schuett

set -e

echo -e "\n========== Copying Changed Files =========="
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

rsync --exclude="*~"          \
      --exclude=".*/"         \
      --exclude="*.pyc"       \
      --ignore-times          \
      --update                \
      --verbose               \
      --recursive             \
      --checksum              \
      /opt/cp2k-local/cp2k/tools/toolchain/  /opt/cp2k-toolchain/

echo -e "\n========== Updating Toolchain =========="
cd /opt/cp2k-toolchain/
./install_cp2k_toolchain.sh --install-all --with-make=no

echo -e "\n========== Compiling CP2K =========="
source /opt/cp2k-toolchain/install/setup
cd /opt/cp2k-master/cp2k/makefiles
make -j VERSION=pdbg cp2k_shell
ln -s /opt/cp2k-master/cp2k/exe/local/cp2k.pdbg /usr/bin/cp2k

echo -e "\n========== Installing ASE =========="
cd /opt/ase/
git pull
pip3 install .

echo -e "\n========== Running ASE Tests =========="
cd /opt/ase/
export ASE_CP2K_COMMAND="mpiexec -np 2 /opt/cp2k-master/cp2k/exe/local/cp2k_shell.pdbg"

set +e # disable error trapping for remainder of script
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
