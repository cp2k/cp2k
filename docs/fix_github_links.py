#!/usr/bin/env python3

# author: Ole Schuett

import sys
from pathlib import Path
import re

LINK_PATTERN = re.compile('<a .*? class="fa fa-github"> Edit on GitHub</a>')


# ======================================================================================
def main() -> None:
    if len(sys.argv) != 2:
        print("fix_github_links.py <build_dir>")
        sys.exit(1)

    print("Fixing GitHub links...")
    build_dir = Path(sys.argv[1])

    replace(build_dir / "bibliography.html", newlink("src/common/bibliography.F"))
    replace(build_dir / "units.html", newlink("src/common/cp_units.F"))

    message = "This page was generated. Please use the smaller [Edit on GitHub] links to see the original location of each string."
    replace(build_dir / "CP2K_INPUT.html", alert(message))
    cp2k_input_dir = build_dir / "CP2K_INPUT"
    for fn in cp2k_input_dir.rglob("*.html"):
        replace(fn, alert(message))


# ======================================================================================
def replace(filename: Path, replacement: str) -> None:
    orig_content = filename.read_text()
    new_content, count = LINK_PATTERN.subn(replacement, orig_content)
    assert count == 1
    filename.write_text(new_content)


# ======================================================================================
def newlink(src_file: str) -> str:
    github_url = f"https://github.com/cp2k/cp2k/blob/master/{src_file}"
    return f'<a href="{github_url}" class="fa fa-github"> Edit on GitHub</a>'


# ======================================================================================
def alert(message: str) -> str:
    return f'<a onclick="alert(\'{message}\')" class="fa fa-github"> Edit on GitHub</a>'


# ======================================================================================
main()

# EOF
