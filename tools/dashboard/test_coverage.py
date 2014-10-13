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

    assert(lines[2] == "="*80)
    assert(lines[-3] == "="*80)
    assert(lines[3][0] == "[")

    coverage = {}
    for l in lines[4:-3]:
        assert(len(l.strip())>0)
        assert(l[0] != "[")
        parts = [p.strip() for p in l.split("|")]
        rate, nlines = parts[1].split()
        assert(rate[-1] == "%")
        coverage[parts[0]] = (float(rate[:-1]), int(nlines))

    if(not path.exists(ref_fn)):
        open(ref_fn, "w").write(pformat(coverage))
        print "Summary: Wrote new reference file"
        print "Status: UNKOWN"
        sys.exit(0)

    ref_coverage = eval(open(ref_fn).read())

    issues = 0
    new_ref_coverage = dict()
    for fn in coverage.keys():
        cov_rate, nlines = coverage[fn]
        uncov_lines = nlines*(100.0 - cov_rate)/100.0
        if(ref_coverage.has_key(fn)):
            cov_rate0, nlines0 = ref_coverage[fn]
            uncov_lines0 = nlines0*(100.0 - cov_rate0)/100.0
            tol = max(nlines, nlines0) * 0.001 # uncov_lines has limited precision
            if(uncov_lines-uncov_lines0>tol and cov_rate<cov_rate0):
                issues += 1
                print('Coverage of file "%s" decreased from %.1f%% to %.1f%%.'%(fn, cov_rate0, cov_rate))
                print('Number of untests lines of file "%s" increased from %d to %d.'%(fn, uncov_lines0, uncov_lines))
                new_ref_coverage[fn] = (cov_rate0, nlines0)
            else:
                new_ref_coverage[fn] = (cov_rate, nlines)
        else:
            if(cov_rate < 90.0):
                issues += 1
                print('New file "%s" has only %.1f%% coverage.'%(fn, cov_rate))
            else:
                new_ref_coverage[fn] = (cov_rate, nlines)

    assert("Total:" in lines[-2])
    total_coverage = float(lines[-2].split("|")[1].split("%")[0])
    print "Summary: Found %d issues, total coverage is at %.1f%%."%(issues, total_coverage)
    print "Status: "+ ("OK" if(issues==0) else "FAILED")

    open(ref_fn, "w").write(pformat(new_ref_coverage))

#===============================================================================
main()

#EOF
