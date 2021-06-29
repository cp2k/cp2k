#!/usr/bin/env python3

# author: Ole Schuett

import argparse
import re
import sys
import os
from os import path
from datetime import datetime

FLAG_EXCEPTIONS = (
    "\$\{.*\}\$",
    "__.*__",
    "CUDA_VERSION",
    "FD_DEBUG",
    "GRID_DO_COLLOCATE",
    "GRID_DO_COLLOCATE",
    "INTEL_MKL_VERSION",
    "LIBINT2_MAX_AM_eri",
    "LIBINT_CONTRACTED_INTS",
    "XC_MAJOR_VERSION",
    "XC_MINOR_VERSION",
    "_OPENMP",
    "__COMPILE_ARCH",
    "__COMPILE_DATE",
    "__COMPILE_HOST",
    "__COMPILE_REVISION",
    "__CRAY_PM_FAKE_ENERGY",
    "__DATA_DIR",
    "__FFTW3_UNALIGNED",
    "__FFTW3_UNALIGNED",
    "__FORCE_USE_FAST_MATH",
    "__FORCE_USE_FAST_MATH",
    "__HAS_smm_cnt",
    "__HAS_smm_ctn",
    "__HAS_smm_ctt",
    "__HAS_smm_dnt",
    "__HAS_smm_dtn",
    "__HAS_smm_dtt",
    "__HAS_smm_snt",
    "__HAS_smm_stn",
    "__HAS_smm_stt",
    "__HAS_smm_znt",
    "__HAS_smm_ztn",
    "__HAS_smm_ztt",
    "__INTEL_COMPILER",
    "__OFFLOAD_CUDA",
    "__OFFLOAD_HIP",
    "__PILAENV_BLOCKSIZE",
    "__PW_CUDA_NO_HOSTALLOC",
    "__PW_CUDA_NO_HOSTALLOC",
    "__RELEASE_VERSION",
    "__RELEASE_VERSION",
    "__T_C_G0",
    "__YUKAWA",
    "__cplusplus",
)
flag_exceptions_re = re.compile("|".join(FLAG_EXCEPTIONS))

portable_filename_re = re.compile(r"^[a-zA-Z0-9._/#~=+-]*$")

BANNER_F = """\
!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-{:d} CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
"""

BANNER_Fypp = """\
#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-{:d} CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!
"""

BANNER_C = """\
/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-{:d} CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/
"""

C_FAMILY_EXTENSIONS = ("c", "cu", "cpp", "h", "hpp")

DEFAULT_EXCLUDED_DIRS = (
    ".git",
    "obj",
    "lib",
    "exe",
    "regtesting",
    "tools/toolchain/build",
    "tools/toolchain/install",
    "tools/autotools",
)

# =======================================================================================
def check_file(fn):
    """
    Check the source files in the given directory/filelist for convention violations, like:

    - correct copyright headers
    - undocumented preprocessor flags
    - stray unicode characters
    """
    warnings = []

    fn_ext = fn.rsplit(".", 1)[-1]
    absfn = path.abspath(fn)
    basefn = path.basename(fn)

    if not portable_filename_re.match(fn):
        warnings += ["Filename %s not portable" % fn]

    if not path.exists(absfn):
        return warnings  # skip broken symlinks

    with open(absfn, "rb") as fhandle:
        raw_content = fhandle.read()

    if b"\0" in raw_content:
        return warnings  # skip binary files

    content = raw_content.decode("utf8")
    if "\r\n" in content:
        warnings += ["Text file %s contains DOS linebreaks" % fn]

    if fn_ext not in ("pot", "patch") and basefn != "Makefile" and "\t" in content:
        warnings += ["Text file %s contains tab character" % fn]

    # check banner
    year = datetime.utcnow().year
    if fn_ext == "F" and not content.startswith(BANNER_F.format(year)):
        warnings += ["%s: Copyright banner malformed" % fn]
    if fn_ext == "fypp" and not content.startswith(BANNER_Fypp.format(year)):
        warnings += ["%s: Copyright banner malformed" % fn]

    if fn_ext in C_FAMILY_EXTENSIONS and not content.startswith(BANNER_C.format(year)):
        warnings += ["%s: Copyright banner malformed" % fn]

    # find all flags
    flags = set()
    line_continuation = False
    for line in content.split("\n"):
        if not line_continuation:
            if len(line) == 0 or line[0] != "#":
                continue
            if line.split()[0] not in ("#if", "#ifdef", "#ifndef", "#elif"):
                continue
        line = line.split("//", 1)[0]
        line_continuation = line.strip().endswith("\\")
        line = re.sub(r"[\\|()!&><=*/+-]", " ", line)
        line = line.replace("defined", " ")
        for m in line.split()[1:]:
            if re.match("[0-9]+[ulUL]*", m):
                continue  # skip numbers
            if fn_ext in ("h", "hpp") and basefn.upper().replace(".", "_") == m:
                continue  # ignore aptly named inclusion guards
            flags.add(m)

    flags = [f for f in flags if not flag_exceptions_re.match(f)]

    cp2k_dir = path.realpath(path.join(path.dirname(__file__), "../../"))
    with open(path.join(cp2k_dir, "INSTALL.md"), encoding="utf8") as fhandle:
        install_txt = fhandle.read()
    with open(path.join(cp2k_dir, "src/cp2k_info.F"), encoding="utf8") as fhandle:
        cp2k_info = fhandle.read()

    cp2k_flags_pattern = r"FUNCTION cp2k_flags\(\)(.*)END FUNCTION cp2k_flags"
    flags_src = re.search(cp2k_flags_pattern, cp2k_info, re.DOTALL).group(1)

    for f in sorted(flags):
        if f not in install_txt:
            warnings += ["Flag %s not mentioned in INSTALL.md" % f]
        if f not in flags_src:
            warnings += ["Flag %s not mentioned in cp2k_flags()" % f]

    return warnings


# =======================================================================================
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: check_file_properties.py <filename>")
        sys.exit(1)

    fn = sys.argv[1]
    warnings = check_file(fn)

    for warning in warnings:
        print(warning)

    if warnings:
        sys.exit(1)


# EOF
