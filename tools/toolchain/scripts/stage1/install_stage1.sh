#!/bin/bash -e

./scripts/stage1/install_mpich.sh
./scripts/stage1/install_openmpi.sh
./scripts/stage1/install_intelmpi.sh
./scripts/stage1/install_valgrind.sh

#EOF
