#!/bin/bash -e

echo "========== running AiiDA-CP2K tests =========="

(
./run_tests.sh
)

EXIT_CODE=$?

echo ""

AIIDA_COMMIT=`git rev-parse --short HEAD`
if (( $EXIT_CODE )); then
    echo "Summary: Something is wrong with aiida-cp2k commit ${AIIDA_COMMIT}."
    echo "Status: FAILED"
else
    echo "Summary: aiida-cp2k commit ${AIIDA_COMMIT} works fine."
    echo "Status: OK"
fi

#EOF
