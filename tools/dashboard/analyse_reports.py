#!/usr/bin/python
# -*- coding: utf-8 -*-

# author: Ole Schuett

import ConfigParser
import sys
import os
import re
import gzip
import numpy as np
from glob import glob

#===============================================================================
def main():
    config = ConfigParser.ConfigParser()
    config.read("dashboard.conf")

    valuesDB = dict()
    nreports = 0
    for s in config.sections():
        if(config.get(s,"report_type") != "regtest"):
            continue
        all_reports = glob("archive/%s/rev_*.txt.gz"%s)
        for fn in all_reports:
            rev = int(fn.rsplit("rev_",1)[1][:-7])
            if(rev < 15256):
                continue # there were no numeric results back then
            reported_values = parse_report(fn)
            if(reported_values):
                nreports += 1
                merge_values(valuesDB, reported_values)

    print("Parsed %d reports\n"%nreports)
    analyze(valuesDB)

#===============================================================================
def merge_values(valuesDB, new_values):
    for k, v in new_values.items():
        if(k not in valuesDB):
            valuesDB[k] = []
        valuesDB[k].append(v)

#===============================================================================
def analyze(valuesDB):
    print("%6s %-80s %-15s  %s"%("","Test Name", "Mean", "Tol."))
    print("-"*110)
    for k, samples in valuesDB.items():
        mean = np.mean(samples)
        max_diff = np.max(np.abs(samples-mean))
        max_diff_rel = max_diff / abs(mean)
        print("STATS: %-70s %20.10e  %10.1e"%(k, mean, max_diff_rel))

#===============================================================================
def parse_report(fn):
    print("Parsing: "+fn)

    values = dict()
    report_txt = gzip.open(fn, 'rb').read()
    m = re.search("\n-+regtesting cp2k-+\n(.*)\n-+ Summary -+\n", report_txt, re.DOTALL)
    if(not m):
        print("Report malformed, skipping.")
        return(None)

    main_part = m.group(1)
    curr_dir = None
    for line in main_part.split("\n"):
        if(line.startswith(">>>")):
            curr_dir = line.rsplit("/tests/")[1] + "/"
        elif(line.startswith("<<<")):
            curr_dir = None
        elif(curr_dir):
            parts = line.split()
            if(not parts[0].endswith(".inp")):
                print("Found strange line:\n"+line)
                continue
            if(parts[1] == "-"):
                continue
            values[curr_dir+parts[0]] = float(parts[1])
        else:
            pass # ignore line

    return values

#===============================================================================
main()
#EOF
