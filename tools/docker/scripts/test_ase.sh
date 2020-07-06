#!/bin/bash -e

# author: Ole Schuett

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

cd /workspace/cp2k
echo -n "Compiling cp2k... "
if make -j VERSION=pdbg &> /dev/null ; then
   echo "done."
else
   echo -e "failed.\n\n"
   echo "Summary: Compilation failed."
   echo "Status: FAILED"
   exit
fi

echo -e "\n========== Installing CP2K =========="
cat > /usr/bin/cp2k_shell << EndOfMessage
#!/bin/bash -e
source /opt/cp2k-toolchain/install/setup
export OMP_NUM_THREADS=1
mpiexec -np 2 /workspace/cp2k/exe/local/cp2k_shell.pdbg "\$@"
EndOfMessage
chmod +x /usr/bin/cp2k_shell

# The cp2k main binary is used by ase/test/cp2k/cp2k_dcd.py.
# https://gitlab.com/ase/ase/merge_requests/1109
cat > /usr/bin/cp2k << EndOfMessage
#!/bin/bash -e
source /opt/cp2k-toolchain/install/setup
export OMP_NUM_THREADS=1
mpiexec -np 2 /workspace/cp2k/exe/local/cp2k.pdbg "\$@"
EndOfMessage
chmod +x /usr/bin/cp2k

mkdir -p ~/.ase
echo '{"cp2k": "/usr/bin/cp2k_shell"}' > ~/.ase/executables.json

echo -e "\n========== Installing ASE =========="
cd /opt/ase/
git pull
pip3 install ".[test]"

echo -e "\n========== Running ASE Tests =========="
cd /opt/ase/
ASE_REVISION=$(git rev-parse --short HEAD)

if ase test -j 0 -c cp2k cp2k ; then
    echo "Summary: ASE commit ${ASE_REVISION} works fine."
    echo "Status: OK"
else
    echo "Summary: Something is wrong with ASE commit ${ASE_REVISION}."
    echo "Status: FAILED"
fi

#EOF
