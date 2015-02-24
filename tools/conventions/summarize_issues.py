#!/usr/bin/python
# -*- coding: utf-8 -*-

# author: Ole Schuett


import sys
from os import path

#===============================================================================
def main():
    if(len(sys.argv) < 2):
        print("Usage: summarize_issues.py [--suppressions=<supp-file>] <issue-file-1> ... <issue-file-N>")
        print("       This combines multiple files with issues into a dashboard-report.\n")
        sys.exit(1)

    suppress = []
    issue_files = sys.argv[1:]
    if sys.argv[1].startswith("--suppressions="):
        content = open(sys.argv[1].split("=")[1]).read()
        lines = [l.strip() for l in content.split("\n") if len(l.strip())>0 ]
        suppress = [l for l in lines if not l.startswith("#")]
        issue_files = sys.argv[2:]

    issues = []
    for fn in issue_files:
        issues += open(fn).read().split("\n")

    issues = [i for i in issues if len(i.strip())>0]
    issues = sorted(set(issues))
    issues_shown = [i for i in issues   if(i not in suppress)]
    issues_supp  = [i for i in issues   if(i     in suppress)]
    unused_supp  = [i for i in suppress if(i not in issues)]

    for i in issues_supp:
        print i+" (suppressed)"

    #for i in unused_supp:
    #    print i+" (unused suppression)"
    print "There are %d unused suppressions.\n"%len(unused_supp)

    for i in issues_shown:
        print i

    n = len(issues_shown)
    m = len(issues_supp)
    print "Summary: Found %d issues (%d suppressed)"%(n, m)
    print "Status: " + ("OK" if n==0 else "FAILED")


#===============================================================================
main()

#EOF
