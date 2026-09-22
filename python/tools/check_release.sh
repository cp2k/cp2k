#!/bin/bash
# SPDX-License-Identifier: GPL-2.0-or-later

# Build and test an installed distribution. Never upload or create release tags.
set -euo pipefail

if (($# != 1)); then
  echo "Usage: bash python/tools/check_release.sh NEW_OUTPUT_DIRECTORY" >&2
  exit 2
fi
PACKAGE_DIR=$(cd "$(dirname "$0")/.." && pwd)
PYTHON=${PYTHON:-python3}
mkdir "$1" # Refuse to overwrite any previous release evidence.
OUTPUT_DIR=$(cd "$1" && pwd)
WORK_DIR=$(mktemp -d "${TMPDIR:-/tmp}/cp2k-python-release.XXXXXX")
trap 'rm -rf "$WORK_DIR"' EXIT
unset PYTHONPATH
export PYTHONNOUSERSITE=1

"$PYTHON" -m venv "$WORK_DIR/build"
"$WORK_DIR/build/bin/python" -m pip install build twine
"$WORK_DIR/build/bin/python" -m build "$PACKAGE_DIR" --outdir "$OUTPUT_DIR/dist"
"$WORK_DIR/build/bin/python" -m twine check "$OUTPUT_DIR"/dist/*
"$WORK_DIR/build/bin/python" -m pip freeze > "$OUTPUT_DIR/build-requirements.txt"

WHEELS=("$OUTPUT_DIR"/dist/*.whl)
if ((${#WHEELS[@]} != 1)) || [[ ! -f ${WHEELS[0]} ]]; then
  echo "Expected exactly one built wheel." >&2
  exit 1
fi
"$PYTHON" -m venv "$WORK_DIR/test"
"$WORK_DIR/test/bin/python" -m pip install "${WHEELS[0]}[test]"
"$WORK_DIR/test/bin/python" -m pip check
"$WORK_DIR/test/bin/python" -m pip freeze > "$OUTPUT_DIR/test-requirements.txt"
cd "$OUTPUT_DIR"
"$WORK_DIR/test/bin/python" -I - << 'PY'
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import sys

import cp2k
import cp2k.library

assert cp2k.__version__ == importlib.metadata.version("cp2k-python")
assert cp2k.library._runtime is None
assert "mpi4py.MPI" not in sys.modules
manifest = {
    "version": cp2k.__version__,
    "python": sys.version,
    "platform": platform.platform(),
    "installed_module": cp2k.__file__,
    "sha256": {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(Path("dist").iterdir())
    },
}
Path("manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
PY
"$WORK_DIR/test/bin/python" -I -m pytest "$PACKAGE_DIR/tests" \
  --import-mode=importlib -m 'not integration' -q -ra \
  --basetemp="$OUTPUT_DIR/tests" --junitxml="$OUTPUT_DIR/unit-tests.xml"
echo "Release artifacts and checks: $OUTPUT_DIR"
echo "Nothing was uploaded. Native integration and release approval are still required."
