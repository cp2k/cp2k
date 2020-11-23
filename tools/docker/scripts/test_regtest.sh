#!/bin/bash

# author: Ole Schuett

if (( $# != 2 )) ; then
    echo "Usage: test_regtest.sh <ARCH> <VERSION>"
    exit 1
fi

ARCH=$1
VERSION=$2

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

# Make OpenMPI happy.
if command -v ompi_info &> /dev/null ; then
    TESTOPTS="-mpiexec 'mpiexec --bind-to none --allow-run-as-root' ${TESTOPTS}"
    export OMPI_MCA_plm_rsh_agent=/bin/false
fi

# Switch to stable DBCSR version if requested.
if [ -n "${USE_STABLE_DBCSR}" ] ; then
    echo "Switching to stable DBCSR version..."
    if ! git -C cp2k/exts/dbcsr checkout -b v2.1.0-rc16 6c52382 ; then
        echo -e "\nSummary: Could not checkout stable DBCSR version."
        echo -e "Status: FAILED\n"
        exit 0
    fi
fi

# Compile cp2k.
echo -en "\nCompiling cp2k... "
cd /workspace/cp2k || exit 1
if make -j ARCH="${ARCH}" VERSION="${VERSION}" &> make.out ; then
    echo "done."
else
    echo -e "failed.\n\n"
    tail -n 100 make.out
    echo -e "\nSummary: Compilation failed."
    echo -e "Status: FAILED\n"
    exit 0
fi

# Check benchmark files for input errors.
if [[ "${ARCH}" == "local" ]] ; then
    echo -en "\nChecking benchmarks... "
    if ! ./tools/regtesting/check_inputs.py "./exe/${ARCH}/cp2k.${VERSION}" "./benchmarks/" ; then
        echo -e "\nSummary: Some benchmark inputs could not be parsed."
        echo -e "Status: FAILED\n"
        exit 0
    fi
fi

# Run regtests.
echo -e "\n========== Running Regtests =========="
make ARCH="${ARCH}" VERSION="${VERSION}" TESTOPTS="${TESTOPTS}" test

exit 0 # Prevent CI from overwriting do_regtest's summary message.

#EOF
