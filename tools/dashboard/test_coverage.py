#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# author: Ole Schuett

import os
import sys
from os import path
from pprint import pformat
from datetime import datetime

#===============================================================================
def main():
    if(len(sys.argv) != 3):
        print("Usage test_coverage.py <reference-file> <lcov-file>")
        sys.exit(1)

    print("Date: %s"%datetime.utcnow().replace(microsecond=0))
    os.system("svn info | grep Revision")

    ref_fn, lcov_fn = sys.argv[1:]
    content = open(lcov_fn).read()
    lines = content.split("\n")
    data_lines = lines[4:-3]
    assert(lines[2] == "="*80)

    coverage = {}
    for l in data_lines:
        parts = [p.strip() for p in l.split("|")]
        p = parts[1].split("%")[0]
        coverage[parts[0]] = float(p)

    if(not path.exists(ref_fn)):
        open(ref_fn, "w").write(pformat(coverage))
        print "Summary: Wrote new reference file"
        print "Status: UNKOWN"
        sys.exit(0)

    ref_coverage = eval(open(ref_fn).read())

    issues = 0
    new_ref_coverage = dict()
    for n, c in coverage.items():
        if(ref_coverage.has_key(n)):
            if(c < ref_coverage[n]):
                issues += 1
                print('Coverage of file "%s" decreased from %.1f%% to %.1f%%.'%(n, ref_coverage[n], c))
            new_ref_coverage[n] = max(c, ref_coverage[n])
        else:
            if(c > 90.0):
                new_ref_coverage[n] = c
            else:
                issues += 1
                print('New file "%s" has only %.1f%% coverage.'%(n, c))

    assert("Total:" in lines[-2])
    total_coverage = float(lines[-2].split("|")[1].split("%")[0])
    print "Summary: Found %d issues, total coverage is at %.f%%."%(issues, total_coverage)
    print "Status: "+ ("OK" if(issues==0) else "FAILED")

    open(ref_fn, "w").write(pformat(new_ref_coverage))

#===============================================================================
main()

#EOF
