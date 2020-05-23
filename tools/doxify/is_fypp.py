#!/usr/bin/env python3

import sys, re

FYPP_SYMBOLS = r"(#|\$|@)"
FYPP_LINE = r"^\s*" + FYPP_SYMBOLS + r":"
FYPP_INLINE = r"(" + FYPP_SYMBOLS + r"{|}" + FYPP_SYMBOLS + r")"
FYPP_RE = re.compile(r"(" + FYPP_LINE + r"|" + FYPP_INLINE + r")")

with open(sys.argv[1], "r") as infile:
    if any(FYPP_RE.search(l) for l in infile):
        sys.exit(0)

sys.exit(1)
