#!/usr/bin/env python3

# author: Ole Schuett

# Runs "cp2k.x --check" on all input files below given path.

import os
import sys
import subprocess
from subprocess import PIPE, DEVNULL


def main():
    if len(sys.argv) != 3:
        print("Usage: check_inputs.py <cp2k-binary> <path-to-inputs>")
        sys.exit(1)

    orig_cwd = os.getcwd()
    cp2k_binary, path_to_check = sys.argv[1:]
    n_files, n_errors = 0, 0

    for root, dirs, files in os.walk(path_to_check):
        # Unpack bzip2 files.
        bz2_files = [fn for fn in files if fn.endswith(".bz2")]
        for fn in bz2_files:
            cmd = ["bzip2", "--decompress", "--keep", "--force", fn]
            subprocess.run(cmd, cwd=root, check=True)

        # Check input files within root directory.
        inp_files = [fn for fn in files if fn.endswith(".inp")]
        for fn in inp_files:
            cmd = [os.path.join(orig_cwd, cp2k_binary), "--check", fn]
            p = subprocess.run(cmd, cwd=root, stdout=PIPE, stderr=DEVNULL)
            n_files += 1
            if p.returncode != 0:
                n_errors += 1
                fn_full = os.path.join(root, fn)
                print(f"\n\nError in {fn_full}:")
                print(p.stdout.decode("utf8"))

    print(f"Found {n_files} input files and {n_errors} errors.")
    sys.exit(n_errors)


main()

# EOF
