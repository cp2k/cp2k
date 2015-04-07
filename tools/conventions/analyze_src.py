#!/usr/bin/python
# -*- coding: utf-8 -*-

# author: Ole Schuett

import re
import os
import sys
from os import path
from datetime import datetime

flag_exceptions_re = re.compile("__COMPILE_.*|__SHORT_FILE__|__INTEL_COMPILER|"
                          +"__cplusplus|_OPENMP|_GNU_SOURCE|__CUDA_ARCH__|"
                          +"cl_.*|CL_VERSION_.*|__OPENCL_VERSION__|__OPENCL")

year = datetime.utcnow().year

BANNER_F = "!-----------------------------------------------------------------------------!\n" \
          +"!   CP2K: A general program to perform molecular dynamics simulations         !\n" \
          +"!   Copyright (C) 2000 - %d  CP2K developers group                          !\n"%year \
          +"!-----------------------------------------------------------------------------!\n"

BANNER_C = "/*****************************************************************************\n" \
          +" *  CP2K: A general program to perform molecular dynamics simulations        *\n" \
          +" *  Copyright (C) 2000 - %d  CP2K developers group                         *\n"%year \
          +" *****************************************************************************/\n"

#===============================================================================
def main():
    if(len(sys.argv) < 2):
        print("Usage: analyse_src.py <cp2k-root-dir>")
        print("       This tool checks the source code for violations of the coding conventions.")
        sys.exit(1)

    cp2k_dir = sys.argv[1]

    flags = set()
    for root, dirs, files in os.walk(path.join(cp2k_dir, "src")):
        if(".svn" in root): continue
        #print("Scanning %s ..."%root)
        for fn in files:
            fn_ext = fn.rsplit(".",1)[-1]
            if(fn_ext in ("template", "instantiate")): continue
            content = open(path.join(root, fn)).read()

            # check banner
            if((fn_ext in ("F", ) and not content.startswith(BANNER_F)) or
               (fn_ext in ("c", "cu", "cpp", "h", "hpp") and not content.startswith(BANNER_C))):
                    print fn+": Copyright banner malformed"
                    #print '"'+ '"\n"'.join(content.split("\n")[:4]) + '"'

            # find all flags
            for line in content.split("\n"):
                if(len(line) == 0): continue
                if(line[0] != "#"): continue
                if(line.split()[0] not in ("#if","#ifdef","#ifndef","#elif")): continue
                line = line.split("//",1)[0]
                line = re.sub("[|()!&><=]", " ", line)
                line = line.replace("defined", " ")
                for m in line.split()[1:]:
                    if m.isdigit(): continue
                    if(fn_ext=="h" and fn.upper().replace(".", "_") == m): continue
                    flags.add(m)

    flags = [f for f in flags if not flag_exceptions_re.match(f)]

    #print("Found %d flags."%len(flags))
    #print flags

    install_txt = open(path.join(cp2k_dir, "INSTALL")).read()
    cp2k_info = open(path.join(cp2k_dir, "src/cp2k_info.F")).read()
    flags_src = re.search(r"FUNCTION cp2k_flags\(\)(.*)END FUNCTION cp2k_flags", cp2k_info, re.DOTALL).group(1)

    for f in sorted(flags):
        if(f not in install_txt):
            print("Flag %s not mentioned in INSTALL"%f)
        if(f not in flags_src):
            print("Flag %s not mentioned in cp2k_flags()"%f)

#===============================================================================

main()
#EOF
