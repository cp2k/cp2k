#!/bin/bash
# vim: dict+=/usr/share/beakerlib/dictionary.vim cpt=.,w,b,u,t,i,k
# shellcheck disable=all
source /usr/share/beakerlib/beakerlib.sh || exit 1

rlJournalStart
  rlPhaseStartTest
    [[ -z "${CP2K_VARIANT}" ]] && rlDie "No CP2K_VARIANT environment variable set"

    rlRun "pushd tests" 0 "cd to cp2k tests folder"

    rlRun "args=\"--maxtasks $(nproc) --ompthreads 2\"" 0 "Set base arguments"
    [[ "${CP2K_SKIP_UNITTEST,,}" == "true" ]] && rlRun "args=\"\$args --skip_unittests\"" 0 "Skip unit tests"
    [[ "${CP2K_SKIP_REGTEST,,}" == "true" ]] && rlRun "args=\"\$args --skip_regtests\"" 0 "Skip regression tests"
    [[ "${CP2K_SMOKE_ONLY,,}" == "true" ]] && rlRun "args=\"\$args --smoketest\"" 0 "Set smoketest flag"

    if [[ "${CP2K_VARIANT}" != "serial" ]]; then
      rlRun "module avail" 0 "Show available modules"
      rlRun "module load mpi/${CP2K_VARIANT}-$(arch)" 0 "Load MPI module: ${CP2K_VARIANT}"
      rlRun "args=\"\$args --mpiranks 2\"" 0 "Set MPI arguments"
      rlRun "args=\"\$args $MPI_BIN psmp\"" 0 "Set run specific arguments"
    else
      rlRun "args=\"\$args /usr/bin ssmp\"" 0 "Set run specific arguments"
    fi
    rlRun "./do_regtest.py $args" 0 "Run regression tests"

    rlRun "popd"
  rlPhaseEnd
rlJournalEnd
