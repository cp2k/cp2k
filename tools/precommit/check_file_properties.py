#!/usr/bin/env python3

# author: Ole Schuett & Tiziano Müller

import argparse
import os
import pathlib
import re
import sys
from datetime import datetime, timezone
from functools import lru_cache
import itertools
from typing import Tuple, List, TypeVar, Iterable

T = TypeVar("T")

# We assume this script is in tools/precommit/
CP2K_DIR = pathlib.Path(__file__).resolve().parents[2]

FLAG_EXCEPTIONS = (
    r"\$\{..*\}\$",
    r"__..*__",
    r"_M_..*",
    r"__ARM_ARCH",
    r"__ARM_FEATURE_..*",
    r"CUDA_VERSION",
    r"DBM_LIBXSMM_PREFETCH",
    r"DBM_VALIDATE_AGAINST_DBCSR",
    r"OPENCL_DBM_SOURCE_MULTIPLY",
    r"FD_DEBUG",
    r"GRID_DO_COLLOCATE",
    r"INTEL_MKL_VERSION",
    r"LIBINT2_MAX_AM_eri",
    r"LIBINT_CONTRACTED_INTS",
    r"XC_MAJOR_VERSION",
    r"XC_MINOR_VERSION",
    r"NDEBUG",
    r"OMP_DEFAULT_NONE_WITH_OOP",
    r"_OPENMP",
    r"__COMPILE_ARCH",
    r"__COMPILE_DATE",
    r"__COMPILE_HOST",
    r"__COMPILE_REVISION",
    r"__CRAY_PM_FAKE_ENERGY",
    r"__DATA_DIR",
    r"__FFTW3_UNALIGNED",
    r"__FORCE_USE_FAST_MATH",
    r"__INTEL_LLVM_COMPILER",
    r"__INTEL_COMPILER",
    r"OFFLOAD_CHECK",
    r"__OFFLOAD_CUDA",
    r"__OFFLOAD_HIP",
    r"__PILAENV_BLOCKSIZE",
    r"__PW_CUDA_NO_HOSTALLOC",
    r"__T_C_G0",
    r"__YUKAWA",
    r"__cplusplus",
    r"HIP_VERSION",
    r"LIBXSMM_GEMM_PREFETCH_NONE",
    r"LIBXSMM_VERSION_NUMBER",
    r"LIBXSMM_VERSION_MAJOR",
    r"LIBXSMM_VERSION_MINOR",
    r"LIBXSMM_VERSION_PATCH",
    r"LIBXSMM_VERSION2",
    r"LIBXSMM_VERSION3",
    r"LIBXSMM_VERSION4",
    r"CPVERSION",
)

FLAG_EXCEPTIONS_RE = re.compile(r"|".join(FLAG_EXCEPTIONS))
PORTABLE_FILENAME_RE = re.compile(r"^[a-zA-Z0-9._/#~=+-]*$")
OP_RE = re.compile(r"[\\|()!&><=*/+-]")
NUM_RE = re.compile(r"[0-9]+[ulUL]*")
CP2K_FLAGS_RE = re.compile(
    r"FUNCTION cp2k_flags\(\)(.*)END FUNCTION cp2k_flags", re.DOTALL
)
STR_END_NOSPACE_RE = re.compile(r'[^ ]"\s*//\s*&')
STR_BEGIN_NOSPACE_RE = re.compile(r'^\s*"[^ ]')
STR_END_SPACE_RE = re.compile(r' "\s*//\s*&')
STR_BEGIN_SPACE_RE = re.compile(r'^\s*" ')


BANNER_F = """\
!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-{:d} CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: {:s}                                                      !
!--------------------------------------------------------------------------------------------------!
"""

BANNER_SHELL = """\
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

BSD_DIRECTORIES = ("src/offload/", "src/grid/", "src/dbm/")


@lru_cache(maxsize=None)
def get_install_txt() -> str:
    return CP2K_DIR.joinpath("INSTALL.md").read_text(encoding="utf8")


@lru_cache(maxsize=None)
def get_flags_src() -> str:
    cp2k_info = CP2K_DIR.joinpath("src/cp2k_info.F").read_text(encoding="utf8")
    match = CP2K_FLAGS_RE.search(cp2k_info)
    assert match
    return match.group(1)


@lru_cache(maxsize=None)
def get_bibliography_dois() -> List[str]:
    bib = CP2K_DIR.joinpath("src/common/bibliography.F").read_text(encoding="utf8")
    matches = re.findall(r'DOI="([^"]+)"', bib, flags=re.IGNORECASE)
    return [doi for doi in matches if "/" in doi]  # filter invalid DOIs.


def check_file(path: pathlib.Path) -> List[str]:
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
    is_executable = os.access(abspath, os.X_OK)

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

    if (
        fn_ext not in (".pot", ".patch")
        and basefn != "Makefile"
        and basefn != "generate_arch_files.sh"
        and "\t" in content
    ):
        warnings += [f"{path}: contains tab character"]

    if fn_ext == ".cu" and "#if defined(_OMP_H)\n#error" not in content:
        warnings += [f"{path}: misses check against OpenMP usage"]

    # Check spaces in Fortran multi-line strings.
    if fn_ext == ".F":
        for i, (a, b) in enumerate(pairwise(content.split("\n"))):
            if STR_END_NOSPACE_RE.search(a) and STR_BEGIN_NOSPACE_RE.search(b):
                warnings += [f"{path}:{i+1} Missing space in multi-line string"]
            if STR_END_SPACE_RE.search(a) and STR_BEGIN_SPACE_RE.search(b):
                warnings += [f"{path}:{i+1} Double space in multi-line string"]

    # check banner
    year = datetime.now(timezone.utc).year
    bsd_licensed = any(str(path).startswith(d) for d in BSD_DIRECTORIES)
    spdx = "BSD-3-Clause    " if bsd_licensed else "GPL-2.0-or-later"
    if fn_ext == ".F" and not content.startswith(BANNER_F.format(year, spdx)):
        warnings += [f"{path}: Copyright banner malformed"]
    if fn_ext == ".fypp" and not content.startswith(BANNER_SHELL.format(year, spdx)):
        warnings += [f"{path}: Copyright banner malformed"]
    if fn_ext == ".cmake" or path.name == "CMakeLists.txt":
        if not content.startswith(BANNER_SHELL.format(year, spdx)):
            warnings += [f"{path}: Copyright banner malformed"]
    if fn_ext in C_EXTENSIONS and not content.startswith(BANNER_C.format(year, spdx)):
        warnings += [f"{path}: Copyright banner malformed"]
    if path.name == "LICENSE" and bsd_licensed and f"2000-{year}" not in content:
        warnings += [f"{path}: Copyright banner malformed"]
    if path.name == "cp2k_info.F" and f'cp2k_year = "{year}"' not in content:
        warnings += [f"{path}: Wrong year."]

    # check shebang
    PY_SHEBANG = "#!/usr/bin/env python3"
    if fn_ext == ".py" and is_executable and not content.startswith(f"{PY_SHEBANG}\n"):
        warnings += [f"{path}: Wrong shebang, please use '{PY_SHEBANG}'"]

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

        line = line.split("/*", 1)[0]  # C comment
        line = line.split("//", 1)[0]  # C++ comment
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
        if fn_ext == ".cl":  # usually compiled at RT (no direct user-control)
            continue
        if flag == "_OMP_H" and fn_ext == ".cu":
            continue
        if flag not in get_install_txt():
            warnings += [f"{path}: Flag '{flag}' not mentioned in INSTALL.md"]
        if flag not in get_flags_src():
            warnings += [f"{path}: Flag '{flag}' not mentioned in cp2k_flags()"]

    # Check for DOIs that could be a bibliography reference.
    if re.match(r"docs/[^/]+/.*\.md", str(path)) and "docs/CP2K_INPUT" not in str(path):
        for line in content.splitlines():
            for doi in get_bibliography_dois():
                if doi.lower() in line:
                    warnings += [f"{path}: Please replace doi:{doi} with biblio ref."]

    return warnings


# ======================================================================================
def pairwise(iterable: Iterable[T]) -> Iterable[Tuple[T, T]]:
    """itertools.pairwise is not available before Python 3.10."""
    # pairwise('ABCDEFG') --> AB BC CD DE EF FG
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)


# ======================================================================================
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
