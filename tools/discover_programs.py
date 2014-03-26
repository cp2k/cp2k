#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re, sys, os
from os import path

# pre-compiled regular expressions
re_program = re.compile(r"\n\s*end\s*program")
re_main = re.compile(r"\sint\s+main\s*\(")

#============================================================================
def main():
    if(len(sys.argv) not in (2,3)):
        print("Usage: discover_programs.py <src-dir> [--ignore-cuda-files]")
        sys.exit(1)

    srcdir = sys.argv[1]

    ignore_cu_files = False
    if(len(sys.argv) > 2):
        assert(sys.argv[2] == "--ignore-cuda-files")
        ignore_cu_files = True

    programs = []
    for root, dirs, files in os.walk(srcdir):
        if(root.endswith("/preprettify")):
            continue
        if("/.svn" in root):
            continue

        for fn in files:
            abs_fn = path.join(root, fn)
            if(fn[-2:] == ".F"):
                if(is_fortran_program(abs_fn)):
                    programs.append(abs_fn)

            elif(fn[-2:] == ".c"):
                if(has_main_function(abs_fn)):
                    programs.append(abs_fn)

            elif(fn[-3:] == ".cu" and not ignore_cu_files):
                if(has_main_function(abs_fn)):
                    programs.append(abs_fn)


    print(" ".join(programs))


#============================================================================
def is_fortran_program(fn):
    f = open(fn)
    s = path.getsize(fn)
    f.seek(max(0, s - 100))
    tail = f.read()
    f.close()
    m = re_program.search(tail.lower())
    return(m != None)

#============================================================================
def has_main_function(fn):
    f = open(fn)
    content = f.read()
    f.close()
    m = re_main.search(content)
    return(m != None)


main()
#EOF
