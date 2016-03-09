#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import re
from os.path import basename

total_matches = 0

#===============================================================================
def main():
    if(len(sys.argv) < 2):
        print("Tool for finding common mistakes with OMP pragmas.")
        print("Hint: It expects DEFAULT(NONE) in the first line of a PARALLEL pragma.")
        print("Usage: find_openmp_mistakes.py <file1.F> ... <fileN.F>")
        sys.exit(1)

    files = sys.argv[1:]
    for fn in files:
        check(fn)

    print("Found %d spots in %d files."%(total_matches, len(files)))

#===============================================================================
def check(fn):
    global total_matches

    f = open(fn)
    all_lines = f.read().split("\n")

    for lineno, line in enumerate(all_lines):
        m = re.search("^\s*!(.*)OMP\s", line, re.IGNORECASE)
        if(m and m.group(1)!="$"):
            total_matches += 1
            print('Found strange OMP stuff "%s" in %s:%d'%(m.group(0), basename(fn), lineno+1))

        m = re.search("!\$OMP\s+CRITICAL\s*$", line, re.IGNORECASE)
        if(m):
            total_matches += 1
            print('Found unnamed OMP CRITICAL in %s:%d'%(basename(fn), lineno+1))


        m = re.search("!\$OMP\s+PARALLEL\s+(.*)$", line, re.IGNORECASE)
        if(m):
            m2 = re.search("default\s*\(none\)", m.group(1), re.IGNORECASE)
            if(not m2):
                total_matches += 1
                print('Found OMP PARALLEL without DEFAULT(NONE) %s:%d'%(basename(fn), lineno+1))
                #print line



#===============================================================================
if(len(sys.argv)==2 and sys.argv[-1]=="--selftest"):
    pass #TODO implement selftest
else:
    main()
#EOF
