#!/bin/bash -e

CP2K_ROOT=$(dirname "$0")

"${CP2K_ROOT}/tools/precommit/precommit.py" --allow-modifications

#EOF
