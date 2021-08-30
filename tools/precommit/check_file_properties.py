#!/usr/bin/env python3

# author: Ole Schuett & Tiziano Müller

import argparse
import pathlib
import re
import sys
from datetime import datetime
from functools import lru_cache
import typing


# We assume this script is in tools/precommit/
CP2K_DIR = pathlib.Path(__file__).resolve().parents[2]

FLAG_EXCEPTIONS = (
    r"\$\{.*\}\$",
    r"__.*__",
    r"CUDA_VERSION",
    r"FD_DEBUG",
    r"GRID_DO_COLLOCATE",
    r"INTEL_MKL_VERSION",
    r"LIBINT2_MAX_AM_eri",
    r"LIBINT_CONTRACTED_INTS",
    r"XC_MAJOR_VERSION",
    r"XC_MINOR_VERSION",
    r"_OPENMP",
    r"__COMPILE_ARCH",
    r"__COMPILE_DATE",
    r"__COMPILE_HOST",
    r"__COMPILE_REVISION",
    r"__CRAY_PM_FAKE_ENERGY",
    r"__DATA_DIR",
    r"__FFTW3_UNALIGNED",
    r"__FORCE_USE_FAST_MATH",
    r"__HAS_smm_cnt",
    r"__HAS_smm_ctn",
    r"__HAS_smm_ctt",
    r"__HAS_smm_dnt",
    r"__HAS_smm_dtn",
    r"__HAS_smm_dtt",
    r"__HAS_smm_snt",
    r"__HAS_smm_stn",
    r"__HAS_smm_stt",
    r"__HAS_smm_znt",
    r"__HAS_smm_ztn",
    r"__HAS_smm_ztt",
    r"__INTEL_COMPILER",
    r"__OFFLOAD_CUDA",
    r"__OFFLOAD_HIP",
    r"__PILAENV_BLOCKSIZE",
    r"__PW_CUDA_NO_HOSTALLOC",
    r"__PW_CUDA_NO_HOSTALLOC",
    r"__RELEASE_VERSION",
    r"__T_C_G0",
    r"__YUKAWA",
    r"__cplusplus",
)

FLAG_EXCEPTIONS_RE = re.compile(r"|".join(FLAG_EXCEPTIONS))
PORTABLE_FILENAME_RE = re.compile(r"^[a-zA-Z0-9._/#~=+-]*$")
OP_RE = re.compile(r"[\\|()!&><=*/+-]")
NUM_RE = re.compile("[0-9]+[ulUL]*")
CP2K_FLAGS_RE = re.compile(
    r"FUNCTION cp2k_flags\(\)(.*)END FUNCTION cp2k_flags", re.DOTALL
)


BANNER_F = """\
!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-{:d} CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: {:s}                                                      !
!--------------------------------------------------------------------------------------------------!
"""

BANNER_FYPP = """\
#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-{:d} CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: {:s}                                                     !
#!-------------------------------------------------------------------------------------------------!
"""

BANNER_C = """\
/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-{:d} CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: {:s}                                 */
/*----------------------------------------------------------------------------*/
"""

C_EXTENSIONS = (".c", ".cu", ".cpp", ".cc", ".h", ".hpp")

BSD_DIRECTORIES = ("src/offload/", "src/grid/")


@lru_cache(maxsize=None)
def get_install_txt() -> str:
    return CP2K_DIR.joinpath("INSTALL.md").read_text()


@lru_cache(maxsize=None)
def get_flags_src() -> str:
    cp2k_info = CP2K_DIR.joinpath("src/cp2k_info.F").read_text()
    match = CP2K_FLAGS_RE.search(cp2k_info)
    assert match
    return match.group(1)


def check_file(path: pathlib.Path) -> typing.List[str]:
    """
    Check the given source file for convention violations, like:

    - correct copyright headers
    - undocumented preprocessor flags
    - stray unicode characters
    """
    warnings = []

    fn_ext = path.suffix
    abspath = path.resolve()
    basefn = path.name

    if not PORTABLE_FILENAME_RE.match(str(path)):
        warnings += [f"Filename '{path}' not portable"]

    if not abspath.exists():
        return warnings  # skip broken symlinks

    raw_content = abspath.read_bytes()

    if b"\0" in raw_content:
        return warnings  # skip binary files

    content = raw_content.decode("utf8")
    if "\r\n" in content:
        warnings += [f"{path}: contains DOS linebreaks"]

    if fn_ext not in (".pot", ".patch") and basefn != "Makefile" and "\t" in content:
        warnings += [f"{path}: contains tab character"]

    # check banner
    year = datetime.utcnow().year
    bsd_licensed = any(str(path).startswith(d) for d in BSD_DIRECTORIES)
    spdx = "BSD-3-Clause    " if bsd_licensed else "GPL-2.0-or-later"
    if fn_ext == ".F" and not content.startswith(BANNER_F.format(year, spdx)):
        warnings += [f"{path}: Copyright banner malformed"]
    if fn_ext == ".fypp" and not content.startswith(BANNER_FYPP.format(year, spdx)):
        warnings += [f"{path}: Copyright banner malformed"]
    if fn_ext in C_EXTENSIONS and not content.startswith(BANNER_C.format(year, spdx)):
        warnings += [f"{path}: Copyright banner malformed"]

    # find all flags
    flags = set()
    line_continuation = False
    for line in content.splitlines():
        line = line.lstrip()

        if not line_continuation:
            if not line or line[0] != "#":
                continue
            if line.split()[0] not in ("#if", "#ifdef", "#ifndef", "#elif"):
                continue

        line = line.split("//", 1)[0]
        line_continuation = line.rstrip().endswith("\\")
        line = OP_RE.sub(" ", line)
        line = line.replace("defined", " ")

        for word in line.split()[1:]:
            if NUM_RE.match(word):
                continue  # skip numbers

            if fn_ext in (".h", ".hpp") and word == basefn.upper().replace(".", "_"):
                continue  # ignore aptly named inclusion guards
            flags.add(word)

    flags = {flag for flag in flags if not FLAG_EXCEPTIONS_RE.match(flag)}

    for flag in sorted(flags):
        if flag not in get_install_txt():
            warnings += [f"{path}: Flag '{flag}' not mentioned in INSTALL.md"]
        if flag not in get_flags_src():
            warnings += [f"{path}: Flag '{flag}' not mentioned in cp2k_flags()"]

    return warnings


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Check the given FILENAME for conventions"
    )
    parser.add_argument("files", metavar="FILENAME", type=pathlib.Path, nargs="+")

    args = parser.parse_args()

    all_warnings = []

    for fpath in args.files:
        all_warnings += check_file(fpath)

    for warning in all_warnings:
        print(warning)

    if all_warnings:
        sys.exit(1)
