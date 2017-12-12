#!/bin/bash -e

date --utc

echo "========== running ASE tests =========="
export ASE_CP2K_COMMAND="mpiexec -np 2 /cp2k/exe/local/cp2k_shell.popt"

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
    echo "Status: FAILED" >> $REPORT
else
    echo "Summary: ASE commit ${ASE_REVISION} works fine."
    echo "Status: OK"
fi

date --utc

#EOF
