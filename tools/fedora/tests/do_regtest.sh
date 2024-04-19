#!/bin/bash
# vim: dict+=/usr/share/beakerlib/dictionary.vim cpt=.,w,b,u,t,i,k
# shellcheck disable=all
source /usr/share/beakerlib/beakerlib.sh || exit 1

rlJournalStart
  rlPhaseStartTest
    [[ -z "${CP2K_VARIANT}" ]] && rlDie "No CP2K_VARIANT environment variable set"

    rlRun "cd ${TMT_TREE}/tests" 0 "cd to cp2k tests folder"

    rlRun "args=\"--maxtasks $(nproc) --ompthreads 2\"" 0 "Set base arguments"
    [[ "${CP2K_SKIP_UNITTEST,,}" == "true" ]] && rlRun "args=\"\$args --skip_unittests\"" 0 "Skip unit tests"
    [[ "${CP2K_SKIP_REGTEST,,}" == "true" ]] && rlRun "args=\"\$args --skip_regtests\"" 0 "Skip regression tests"
    [[ "${CP2K_SMOKE_ONLY,,}" == "true" ]] && rlRun "args=\"\$args --smoketest\"" 0 "Set smoketest flag"

    if [[ "${CP2K_VARIANT}" != "serial" ]]; then
      rlRun "module avail" 0 "Show available modules"
      rlRun "module load mpi/${CP2K_VARIANT}" 0 "Load MPI module: ${CP2K_VARIANT}"
      rlRun "export CP2K_STEM=$MPI_BIN/cp2k" 0 "Export CP2K_STEM"
      rlRun "args=\"\$args --mpiranks 2\"" 0 "Set MPI arguments"
      rlRun "args=\"\$args local_${CP2K_VARIANT} psmp\"" 0 "Set run specific arguments"
    else
      rlRun "export CP2K_STEM=/usr/bin/cp2k" 0 "Export CP2K_STEM"
      rlRun "args=\"\$args local ssmp\"" 0 "Set run specific arguments"
    fi
    rlRun "./do_regtest.py $args" 0 "Run regression tests"
  rlPhaseEnd
rlJournalEnd
