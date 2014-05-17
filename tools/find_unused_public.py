#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re, sys
from os import path
from os.path import dirname, basename, normpath
import os

# pre-compiled regular expressions
re_module  = re.compile(r"(?:^|\n)\s*module\s+(\w+)\s.*\n\s*end\s*module",re.DOTALL)
re_useonly = re.compile(r"\n\s*use\s+(\w+)\s*,\s*only\s*:(.*)(?=\n)")
re_use     = re.compile(r"\n\s*use\s+(\w+)\s*(?=\n)")
re_pub     = re.compile(r"\n\s*public\s*::\s*(.*)(?=\n)")


#=============================================================================
def main():
    src_dir = "../src"
    parsed_files = {}
    for root, dir, files in os.walk(src_dir):
        if("preprettify" in root):
            continue
        if(".svn" in root):
            continue
        for fn in files:
            if(not fn.endswith(".F")):
                continue
            absfn = root+"/"+fn

            parsed_files[absfn] =  parse_file(absfn)

    all_used_symboles = set()
    for p in parsed_files.values():
        for m, syms in p['use']:
            for s in syms:
                all_used_symboles.add(m+"@"+s)

    n = 0
    for fn, p in parsed_files.items():
        if(len(p['mod']) != 1):
            continue
        m = p['mod'][0]
        if(m+"@*" in all_used_symboles):
            continue
        unused = []
        for s in p['pub']:
            if(m+"@"+s not in all_used_symboles):
                unused.append(s)
                n += 1
        if(len(unused) > 0):
            print("%s USElessly declares PUBLIC: "%fn+ ", ".join(unused)+"\n")


    print("Found %d unUSEd PUBLIC symbols."%n)


#=============================================================================
def parse_file(fn):
    # print("parsing "+fn)
    content = open(fn).read()
    # re.IGNORECASE is horribly expensive. Converting to lower-case upfront
    content = content.lower()
    content = re.sub("!.*\n", "\n", content)
    content = re.sub("&\s*\n", "", content)
    content = re.sub("&", " ", content)

    mods = re_module.findall(content)

    uses = []
    matches = re_use.findall(content)
    for m in matches:
        uses.append((m.strip(), ("*",)))

    matches = re_useonly.findall(content)
    for m in matches:
        syms = [p.strip() for p in m[1].split(",")]
        uses.append((m[0].strip(), syms))

    publics = []
    matches = re_pub.findall(content)
    for m in matches:
        publics += [p.strip() for p in m.split(",")]

    return({"mod":mods, "use":uses, "pub":publics})


#=============================================================================
main()

#EOF
