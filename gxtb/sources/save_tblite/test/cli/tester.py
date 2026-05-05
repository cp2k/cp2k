#!/usr/bin/env python3
"""
Minimal Python wrapper for testing the tblite command line interface.

The wrapper will assume a specific order in the arguments rather than
providing a generic command line interface by itself since it is
supposed to be used by meson for testing purposes only.
"""

try:
    import subprocess, sys, json, os
except ImportError:
    exit(77)

if len(sys.argv) < 4:
    raise RuntimeError("Requires at least four arguments")


class approx:
    def __init__(self, value, abs):
        self.value = value
        self.abs = abs
        self.log = []

    def __eq__(self, other):
        def compare(a, b, ctx):
            if isinstance(a, list) and isinstance(b, list):
                return all(compare(x, y, f"{ctx}[{idx}]") for idx, (x, y) in enumerate(zip(a, b)))

            if isinstance(a, dict) and isinstance(b, dict):
                try:
                    return all(compare(a[k], b[k], f"{ctx}[{k}]") for k in a.keys())
                except KeyError as e:
                    self.log.append(f"{ctx}: Missing key {e} in dictionary")
                    return False

            if isinstance(a, (int, float)) and isinstance(b, (int, float)):
                stat = abs(a - b) < self.abs
                if not stat:
                    self.log.append(f"{ctx}: {a} != {b} (diff: {abs(a - b)})")
                return stat

            stat = a == b
            if not stat:
                self.log.append(f"{ctx}: {a} != {b}")
            return stat

        stat = compare(self.value, other, "")
        if not stat:
            print("\n".join(self.log))

        return stat


thr = 1.0e-6
prog = sys.argv[1]
outp = sys.argv[2]
with open(sys.argv[3]) as fd:
    wdir = os.path.dirname(fd.name)
    args = [arg.replace('$ORIGIN', wdir) for arg in fd.read().strip().split("\n")]

stat = subprocess.call(
    [prog, *args, "--json", os.path.basename(outp)],
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

assert approx(ref, abs=thr) == res
