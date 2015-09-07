#!/usr/bin/python
# -*- coding: utf-8 -*-

# author: Ole Schuett

import sys
from os import path

#===============================================================================
def main():
    if(len(sys.argv) != 2):
        print("Usage prettify_TEST_FILES.py <tests-dir>")
        sys.exit(1)

    tests_dir = sys.argv[1]
    assert(tests_dir.endswith("/"))

    lines = open(tests_dir+"TEST_DIRS").readlines()
    test_subdirs = [l.split()[0] for l in lines if l[0]!="#"]
    for d in test_subdirs:
        fn = tests_dir+d+"/TEST_FILES"
        print "Working on: "+fn
        content = open(fn).read()
        output = ""
        for line in content.strip().split("\n"):
            line = line.strip()
            if(line == "#EOF"):
                continue
            if(line.startswith("#")):
                output += line + "\n"
                continue
            parts = line.split()
            assert(len(parts) < 5)
            assert(len(parts[0]) < 50)
            assert(len(parts[1]) < 5)
            output += "%-50s %5s" %(parts[0], parts[1])
            if(len(parts) > 2):
                assert(len(parts[2]) < 10)
                output += " %10s" %parts[2]
            if(len(parts) > 3):
                assert(len(parts[3]) < 30)
                output += " %30s" %parts[3]
            output += "\n"
        output += "#EOF\n"

        f = open(fn, "w")
        f.write(output)
        f.close()

main()
#EOF
