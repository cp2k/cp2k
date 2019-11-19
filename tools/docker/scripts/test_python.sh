#!/bin/bash -e

# author: Ole Schuett

function run_tests {
    PYTHON=$1
    SCRIPTS=$2
    BASEDIR=$(pwd)
    VERSION=$(${PYTHON} -V 2>&1 | tr -d '\n')
    echo ""
    echo "========== Running ${VERSION} Tests =========="
    #find ./src/ ./tools/ -name "*.pyc" -exec rm {} \;
    for i in $SCRIPTS ; do
        set +e #disable error trapping
        if [[ $(head -n1 "$i") =~ "python3" && "$PYTHON" != "python3" ]]; then
          echo "Skipping $i - it's Python3 only."
          continue
        fi
        echo "Running $i"
        cd "$(dirname "$i")"
        SCRIPT=$(basename "$i")
        if ! $PYTHON -B "./${SCRIPT}" ; then
            echo "Test $i failed"
            ERRORS=$((ERRORS+1))
        fi
        set -e #re-enable error trapping
        cd "$BASEDIR"
    done
}

function run_pre_commit {
    set +e #disable error trapping
    if ! pre-commit run --all-files --hook-stage manual check-ast ; then
        echo "Syntax check failed for either Python 2 or 3"
        ERRORS=$((ERRORS+1))
    fi
    set -e #re-enable error trapping
}

#===============================================================================
ERRORS=0

cd /workspace/cp2k

# find executable python scripts
ALL_TEST_SCRIPTS=$(find ./src/ ./tools/ -name "*_test.py"  -executable)
ESSENTIAL_TEST_SCRIPTS=$(find ./tools/build_utils -name "*_test.py"  -executable)

# python 2.7
run_tests python2.7 "${ALL_TEST_SCRIPTS}"

# python 3.x
run_tests python3 "${ALL_TEST_SCRIPTS}"

# run manual pre-commit checks
run_pre_commit

# print report
NUM_TEST_SCRIPTS=$(( $(wc -w <<< "${ALL_TEST_SCRIPTS}") + 1 ))
echo ""
echo "Summary: Ran ${NUM_TEST_SCRIPTS} Python tests, found ${ERRORS} errors."
if [ $ERRORS == 0 ]; then
    echo "Status: OK"
else
    echo "Status: FAILED"
fi

#EOF
