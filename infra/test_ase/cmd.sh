#!/bin/bash -e

echo "========== running ASE tests =========="
export OMP_NUM_THREADS=1
export ASE_CP2K_COMMAND="mpiexec -np 2 cp2k_shell.psmp"

(
for i in ./ase/test/cp2k/cp2k_*.py
do
  echo "Running $i ..."
  python $i
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
