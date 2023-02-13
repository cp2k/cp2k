#!/usr/bin/env python3

# author: Ole Schuett

from difflib import unified_diff
from pathlib import Path
from typing import List
import argparse
import sys

# ======================================================================================
def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    cp2k_root = Path(__file__).parent.parent
    output_path = cp2k_root / "src" / "CMakeLists.txt"
    old_lines = output_path.read_text(encoding="utf8").split("\n")
    banner = old_lines[:6]

    new_lines = (
        banner
        + preemble()
        + filelist(cp2k_root, "F", "CP2K_SRCS_F")
        + filelist(cp2k_root, "c", "CP2K_SRCS_C")
        + filelist(cp2k_root, "cu", "CP2K_SRCS_GPU")
        + ["include(cp2k_main_targets)", ""]
        + ["# EOF", ""]
    )

    # If requested, check existing file...
    if args.check:
        diff = list(unified_diff(old_lines, new_lines, "got", "expected", lineterm=""))
        if diff:
            print(f"File {output_path} is not consistent with generator script:")
            print("\n".join(diff))
            sys.exit(1)
        else:
            print(f"File {output_path} is consistent with generator script.")
            sys.exit(0)

    # ...otherwise, write new output.
    output_path.write_text("\n".join(new_lines), encoding="utf8")
    print(f"Wrote {output_path}")


# ======================================================================================
def filelist(cp2k_root: Path, file_extension: str, cmake_list_name: str) -> List[str]:
    src_dir = cp2k_root / "src"
    files_list = list(src_dir.glob(f"**/*.{file_extension}"))
    files_rel_list = [x.relative_to(src_dir) for x in files_list]
    files_str = [f"  {x}" for x in sorted(files_rel_list)]
    files_str[-1] = files_str[-1] + ")"
    return ["list(", "  APPEND", f"  {cmake_list_name}"] + files_str + [""]


# ======================================================================================
def preemble() -> List[str]:
    return [
        "",
        "#",
        "# This file was created by ../cmake/generate_cmakelists.py",
        "#",
        "",
        "configure_file(base/base_uses.f90 base/base_uses.f90 @ONLY)",
        "",
        "# Begrudgingly tolerated because PyTorch has no C API.",
        "list(APPEND CP2K_SRCS_CPP torch_c_api.cpp)",
        "",
    ]


# ======================================================================================
main()

# EOF
