#!/usr/bin/env python3

# author: Ole Schuett

import re
import sys
from pathlib import Path


def main() -> None:
    if len(sys.argv) != 2:
        print("Usage: format_makefile.py <file>")
        sys.exit(1)
    makefile = Path(sys.argv[1])

    lines_out = []
    continuation = False
    for line in makefile.read_text(encoding="utf8").split("\n"):

        # Remove trailing whitespaces.
        line = line.rstrip()

        # Detect continued lines.
        prev_continuation = continuation
        continuation = line.endswith("\\")

        # Continued lines are indented 8 spaces.
        if prev_continuation:
            lines_out.append(" " * 8 + line.strip())

        # Tabbed lines are indented with excatly one tab.
        elif line.startswith("\t"):
            lines_out.append("\t" + line.strip())

        # All other lines are not indented.
        else:
            lines_out.append(line.strip())

    makefile.write_text("\n".join(lines_out), encoding="utf8")


main()

# EOF
