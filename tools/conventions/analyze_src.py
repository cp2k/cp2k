#!/usr/bin/env python3

# author: Ole Schuett

import argparse
import re
import sys
import os
from os import path
from datetime import datetime

flag_exceptions_re = re.compile(
    r"__COMPILE_REVISION|__COMPILE_DATE|__COMPILE_ARCH|__COMPILE_HOST|"
    r"__INTEL_COMPILER|__cplusplus|\$\{.*\}\$|__.*__"
)

portable_filename_re = re.compile(r"^[a-zA-Z0-9._/#~=+-]*$")

BANNER_F = """\
!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - {:4}  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!
"""

BANNER_Fypp = """\
#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright (C) 2000 - {:4}  CP2K developers group                                              !
#!-------------------------------------------------------------------------------------------------!
"""

BANNER_C = """\
/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - {:d}  CP2K developers group                         *
 *****************************************************************************/
"""

DEFAULT_EXCLUDED_DIRS = (
    ".git",
    "obj",
    "lib",
    "exe",
    "regtesting",
    "tools/toolchain/build",
    "tools/toolchain/install",
)


def validate(cp2k_dir, filelist=None, excluded_dirs=DEFAULT_EXCLUDED_DIRS):
    """
    Check the source files in the given directory/filelist for convention violations, like:

    - correct copyright headers
    - undocumented preprocessor flags
    - stray unicode characters

    :param cp2k_dir: base directory to look for files
    :param filelist: specific list of files to check in the given directory
    :param excluded_dirs: List of directories to exclude in the checks (usually relative to `cp2k_dir`)
    """
    # check flags and banners
    flags = set()
    year = datetime.utcnow().year
    warnings = []

    # directories to exclude, if given as relative dirs, they're relative to cp2k_dir
    excluded_dirs = [path.join(cp2k_dir, p) for p in DEFAULT_EXCLUDED_DIRS]

    # also exclude any Git submodules
    if path.exists(path.join(cp2k_dir, ".gitmodules")):
        with open(path.join(cp2k_dir, ".gitmodules"), "r") as fhandle:
            excluded_dirs += [
                path.join(cp2k_dir, l.split()[-1]) for l in fhandle if "path =" in l
            ]

    if filelist:
        # limited emulation of the os.walk generator:
        # triples of (dir, list-of-subdirs, list-of-files)
        fileiter = [
            (path.join(cp2k_dir, path.dirname(f)), [], [path.basename(f)])
            for f in filelist
        ]
    else:
        fileiter = os.walk(path.join(cp2k_dir, "src"))

    for root, dirs, files in fileiter:
        dirs[:] = [d for d in dirs if path.join(root, d) not in excluded_dirs]

        for fn in files:
            fn_ext = fn.rsplit(".", 1)[-1]
            if fn_ext in ("template", "instantiate"):
                continue

            with open(path.join(root, fn)) as fhandle:
                content = fhandle.read()

            # check banner
            if (
                (fn_ext in ("F",) and not content.startswith(BANNER_F.format(year)))
                or (
                    fn_ext in ("fypp",)
                    and not content.startswith(BANNER_Fypp.format(year))
                )
                or (
                    fn_ext in ("c", "cu", "cpp", "h", "hpp")
                    and not content.startswith(BANNER_C.format(year))
                )
            ):
                warnings += ["%s: Copyright banner malformed" % fn]

            # find all flags
            for line in content.split("\n"):
                if len(line) == 0:
                    continue
                if line[0] != "#":
                    continue
                if line.split()[0] not in ("#if", "#ifdef", "#ifndef", "#elif"):
                    continue
                line = line.split("//", 1)[0]
                line = re.sub("[|()!&><=*/+-]", " ", line)
                line = line.replace("defined", " ")
                for m in line.split()[1:]:
                    if re.match("[0-9]+[ulUL]*", m):
                        continue  # skip numbers
                    if fn_ext in ("h", "hpp") and fn.upper().replace(".", "_") == m:
                        continue  # ignore aptly named inclusion guards
                    flags.add(m)

    flags = [f for f in flags if not flag_exceptions_re.match(f)]

    with open(path.join(cp2k_dir, "INSTALL.md")) as fhandle:
        install_txt = fhandle.read()
    with open(path.join(cp2k_dir, "src/cp2k_info.F")) as fhandle:
        cp2k_info = fhandle.read()

    flags_src = re.search(
        r"FUNCTION cp2k_flags\(\)(.*)END FUNCTION cp2k_flags", cp2k_info, re.DOTALL
    ).group(1)

    for f in sorted(flags):
        if f not in install_txt:
            warnings += ["Flag %s not mentioned in INSTALL.md" % f]
        if f not in flags_src:
            warnings += ["Flag %s not mentioned in cp2k_flags()" % f]

    if not filelist:
        # check for copies of data files
        data_files = set()
        for _, _, files in os.walk(path.join(cp2k_dir, "data")):
            data_files.update(files)
        data_files.remove("README")
        for root, _, files in os.walk(path.join(cp2k_dir, "tests")):
            d = path.relpath(root, cp2k_dir)
            for c in data_files.intersection(files):
                warnings += ["Data file %s copied to %s" % (c, d)]

    if filelist:
        fileiter = [(cp2k_dir, [], filelist)]
    else:
        fileiter = os.walk(cp2k_dir)

    # check linebreaks and encoding
    for root, dirs, files in fileiter:
        dirs[:] = [d for d in dirs if path.join(root, d) not in excluded_dirs]

        for fn in files:
            absfn = path.join(root, fn)
            shortfn = path.relpath(absfn, cp2k_dir)

            if not portable_filename_re.match(fn):
                warnings += ["Filename %s not portable" % shortfn]

            if not path.exists(absfn):
                continue  # skip broken symlinks

            with open(absfn, "rb") as fhandle:
                content = fhandle.read()

            if b"\0" in content:
                continue  # skip binary files
            if b"\r\n" in content:
                warnings += ["Text file %s contains DOS linebreaks" % shortfn]

            # check for non-ascii chars
            if re.search(b"[\x80-\xFF]", content):
                if absfn.endswith(".py"):
                    continue  # python files are utf8 encoded by default

                if b"# -*- coding: utf-8 -*-" in content:
                    continue

                for lineno, line in enumerate(content.splitlines()):
                    m = re.search(b"[\x80-\xFF]", line)
                    if m:
                        warnings += [
                            "Found non-ascii char in %s line %d at position %d"
                            % (shortfn, lineno + 1, m.start(0) + 1)
                        ]

    return warnings


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Check source code for coding convention violations"
    )
    parser.add_argument(
        "-b", "--base-dir", type=str, default=".", help="CP2K base directory to check"
    )
    parser.add_argument(
        "file",
        type=str,
        nargs="*",
        help="Limit the test to given files (given as relative paths to the base dir, otherwise all relevant files will be scanned)",
    )
    parser.add_argument(
        "--fail",
        action="store_true",
        help="return non-0 exit code if warnings have been found (useful for pre-commit scripts)",
    )
    parser.add_argument(
        "--suppressions",
        type=argparse.FileType("r"),
        help="Specify a suppression file with messages to not count towards the number of warnings",
    )
    args = parser.parse_args()

    warnings = validate(args.base_dir, args.file)

    if args.suppressions:
        suppress = [
            l.rstrip() for l in args.suppressions if l and not l.startswith("#")
        ]
        warnings = [w for w in warnings if not any(s in w for s in suppress)]

    for warning in warnings:
        print(warning)

    if args.fail and warnings:
        sys.exit(1)
