#!/usr/bin/env python3
"""
Minimal Python wrapper for testing the dftd4 command line interface.

The wrapper will assume a specific order in the arguments rather than
providing a generic command line interface by itself since it is
supposed to be used by meson for testing purposes only.
"""

try:
    import subprocess, sys, json, os, pytest
except ImportError:
    exit(77)

if len(sys.argv) < 4:
    raise RuntimeError("Requires at least four arguments")

thr = 1.0e-9
prog = sys.argv[1]
outp = sys.argv[2]
args = sys.argv[3:]

stat = subprocess.call(
    [prog, "--json", os.path.basename(outp)] + args,
    shell=False,
    stdin=None,
    stderr=subprocess.STDOUT,
)
if stat != 0:
    raise RuntimeError("Calculation failed")

with open(outp) as f:
    ref = json.load(f)
    del ref["version"]

with open(os.path.basename(outp)) as f:
    res = json.load(f)

for key in ref:
    if key not in res:
        raise RuntimeError("Missing '" + key + "' entry in results")
    assert pytest.approx(res[key], abs=thr) == ref[key]
