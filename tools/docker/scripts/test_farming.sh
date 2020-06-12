#!/bin/bash -e

# author: Ole Schuett

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

echo -e "\n========== Running Regtests =========="
cd /workspace/cp2k
make VERSION=pdbg test TESTOPTS="-farming -skipunittest -skipdir TMC/regtest_ana_on_the_fly -skipdir TMC/regtest_ana_post_proc -skipdir TMC/regtest ${TESTOPTS}"

#EOF
