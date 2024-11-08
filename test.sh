#!/usr/bin/env bash

set -x

export DIR=$(cd $(dirname $BASH_SOURCE) && pwd)
export CP2K_DATA_DIR="$DIR/data"
build_folder="$(pwd)"

# pushd "$DIR/tests/QS/regtest-ot-1/"
"$build_folder/bin/cp2k.psmp" -i "$DIR/tests/QS/regtest-ot/H2O-geo-ot-lumo-all.inp" #-o "$build_folder/H2O-OT-1.out"

# python tests/do_regtest.py \
#     --restrictdir QS/regtest-ot-1 \
#     --workbasedir TEST \
#     build/bin psmp

