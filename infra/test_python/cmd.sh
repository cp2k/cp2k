#!/bin/bash

# author: Ole Schuett

set -e

run_selftests() {
    PYTHON=$1
    shift
    BASEDIR=`pwd`
    VERSION=`$PYTHON -V 2>&1 | tr -d '\n'`
    echo ""
    echo "========== Running ${VERSION} Tests =========="
    #find ./src/ ./tools/ -name "*.pyc" -exec rm {} \;
    for i in $@ ; do
        set +e #disable error trapping
        if [[ `head -n1 $i` =~ "python3" && "$PYTHON" != "python3" ]]; then
          echo "Skipping $i - it's Python3 only."
          continue
        fi
        echo "Selftesting" $i
        cd $(dirname $i)
        $PYTHON -B ./$(basename $i) --selftest
        if [ $? -ne 0 ]; then
            echo "Selftest of $i failed"
            ERRORS=$((ERRORS+1))
        fi
        set -e #re-enable error trapping
        cd $BASEDIR
    done
}


#===============================================================================
echo -e "\n========== Copying Changed Files =========="
mkdir /opt/cp2k-local-rw/
cp -r /opt/cp2k-local/cp2k/tools /opt/cp2k-local-rw/tools
cp -r /opt/cp2k-local/cp2k/src /opt/cp2k-local-rw/src
cd /opt/cp2k-local-rw/


ERRORS=0

# find executable python scripts
ALL_SCRIPTS=`find ./src/ ./tools/ -name "*.py"  -executable`
ESSENTIAL_SCRIPTS=`find ./tools/build_utils -name "*.py"  -executable`
ESSENTIAL_SCRIPTS="$ESSENTIAL_SCRIPTS ./src/dbcsr/libsmm_acc/libcusmm/generate.py"

# python 2.6
run_selftests python2.6 $ESSENTIAL_SCRIPTS

# python 2.7
run_selftests python2.7 $ALL_SCRIPTS

# python 3.x
run_selftests python3 $ALL_SCRIPTS

# print report
NUM_SCRIPTS=`echo $ALL_SCRIPTS | wc -w`
echo ""
echo "Summary: Checked ${NUM_SCRIPTS} Python scripts, found ${ERRORS} errors."
if [ $ERRORS == 0 ]; then
    echo "Status: OK"
else
    echo "Status: FAILED"
fi

#EOF
