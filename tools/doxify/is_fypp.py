import sys, re

FYPP_SYMBOLS = r"(#|\$|@)"
FYPP_LINE = r"^\s*" + FYPP_SYMBOLS + r":"
FYPP_INLINE = r"(" + FYPP_SYMBOLS + r"{|}" + FYPP_SYMBOLS + r")"
FYPP_RE = re.compile(r"(" + FYPP_LINE + r"|" + FYPP_INLINE + r")")

infile = open(sys.argv[1], 'r')
for line in infile.readlines():
    if FYPP_RE.search(line):
        sys.exit(0)

sys.exit(1)
