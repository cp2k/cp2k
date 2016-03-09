#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
from os.path import basename

threshold = 5
total_matches = 0

#===============================================================================
def main():
    if(len(sys.argv) < 2):
        print("Tool for finding large commented regions of code. Uses very simply heuristics.")
        print("Usage: find_bitrot_code <file1.F> ... <fileN.F>")
        sys.exit(1)

    files = sys.argv[1:]
    for fn in files:
        check(fn)

    print("Found %d spots in %d files."%(total_matches, len(files)))
#===============================================================================
def check(fn):
    f = open(fn)
    all_lines = f.read().split("\n")

    counter = 0
    start = 0

    def report():
        global total_matches
        total_matches += 1
        print(" +"+("-"*87)+"+")
        msg = "%s: line %d ... line %d"%(basename(fn), start+1, lineno+1)
        print(" | "+msg.ljust(85) +" |")
        print(" +"+("-"*87)+"+")
        for i in range(start, lineno):
            print(" | %4d "%(i)+all_lines[i].ljust(80) +" |")
        print(" +"+("-"*87)+"+")
        print("\n")

    for lineno, line in enumerate(all_lines):
        s = line.strip().lower()

        if(s.startswith("!>") or s.startswith("!$omp")):
            if(counter > threshold): report()
            counter = 0
        elif(s.startswith("!")):
            if(counter==0):
                start=lineno
            counter += 1
        elif(len(s) > 0):
            if(counter > threshold): report()
            counter = 0


#===============================================================================
if(len(sys.argv)==2 and sys.argv[-1]=="--selftest"):
    pass #TODO implement selftest
else:
    main()
#EOF
