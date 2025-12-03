#!/usr/bin/env python3

# author: Ole Schuett

# Redirects gfortran's stdout/stderr to files e.g. when run with -fdump-fortran-original.
# To use it pass the following argument to CMake:
#   -DCMAKE_Fortran_COMPILER_LAUNCHER="redirect_gfortran_output.py"

import sys
import os.path
import subprocess


# ======================================================================================
def main() -> None:
    assert sys.argv[-2] == "-o" and sys.argv[-1].endswith(".o")
    stem = os.path.basename(sys.argv[-1][:-2])
    ast_fn, warn_fn = f"{stem}.ast", f"{stem}.warn"
    with open(ast_fn, "w") as stdout, open(warn_fn, "w") as stderr:
        p = subprocess.run(sys.argv[1:], stdout=stdout, stderr=stderr)

    if p.returncode != 0:
        sys.stderr.write(open(warn_fn, "r").read())

    sys.exit(p.returncode)


# ======================================================================================
main()

# EOF
