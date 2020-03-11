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
if which ompi_info &> /dev/null ; then
    TESTOPTS="-mpiexec 'mpiexec --bind-to none --allow-run-as-root' ${TESTOPTS}"
    export OMPI_MCA_plm_rsh_agent=/bin/false
fi

echo -en "\nCompiling cp2k... "
cd /workspace/cp2k
if make -j ARCH="${ARCH}" VERSION="${VERSION}" &> make.out ; then
    echo "done."
else
    echo "failed."
    cat make.out
    echo -e "\nSummary: Compilation failed."
    echo -e "Status: FAILED\n"
    exit 0
fi


if [[ "${ARCH}" == "local" ]] ; then
    echo -en "\nChecking benchmarks... "
    if ! ./tools/regtesting/check_inputs.py "./exe/${ARCH}/cp2k.${VERSION}" "./benchmarks/" ; then
        echo -e "\nSummary: Some benchmark inputs could not be parsed."
        echo -e "Status: FAILED\n"
        exit 0
    fi
fi


echo -e "\n========== Running Regtests =========="
make ARCH="${ARCH}" VERSION="${VERSION}" TESTOPTS="${TESTOPTS}" test

exit 0 # Prevent CI from overwriting do_regtest's summary message.

#EOF
