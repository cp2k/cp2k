#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check that a static archive contains only the objects specified in the PACKAGES files.
"""

# author: Ole Schuett


import subprocess
import os
from os import path
import argparse

KNOWN_EXTENSIONS = ("F", "c", "cu", "cpp", "cxx", "cc")


def main(ar_exe, src_dir, lib_dir):
    # Search all files belonging to a given archive
    archives_files = {}

    for root, _, files in os.walk(src_dir):
        if "PACKAGE" in files:
            with open(path.join(root, "PACKAGE")) as fhandle:
                content = fhandle.read()

            package = eval(content)

            archive = "libdbcsr" + path.basename(root)
            if "archive" in package.keys():
                archive = package["archive"]

            file_parts = [fn.rsplit(".", 1) for fn in files]
            src_basenames = [
                parts[0] for parts in file_parts if parts[-1] in KNOWN_EXTENSIONS
            ]

            if archive in archives_files:
                archives_files[archive] |= set(src_basenames)
            else:
                archives_files[archive] = set(src_basenames)

    # Check if the symbols in each archive have a corresponding source file
    for archive in archives_files:
        archive_fn = path.join(lib_dir, archive + ".a")

        if not path.exists(archive_fn):
            continue

        output = subprocess.check_output([ar_exe, "t", archive_fn])
        for line in output.decode("utf8").strip().splitlines():
            if line == "__.SYMDEF SORTED":
                continue  # needed for MacOS

            assert line.endswith(
                ".o"
            ), "discovered a non-object file inside a static archive"

            if line[:-2] not in archives_files[archive]:
                print(
                    "Could not find source for object '{}' in archive '{}', removing archive.".format(
                        line, archive_fn
                    )
                )
                os.remove(archive_fn)
                break


# ===============================================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
        Parse files and package manifests in the source tree to create rules for objects and executables

        This script is part of the build utility scripts for DBCSR.
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("ar_executable", metavar="ar_executable", type=str)
    parser.add_argument("src_dir", metavar="src_dir", type=str)
    parser.add_argument("lib_dir", metavar="lib_dir", type=str)

    args = parser.parse_args()
    main(ar_exe=args.ar_executable, src_dir=args.src_dir, lib_dir=args.lib_dir)
