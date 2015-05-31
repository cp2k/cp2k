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
from os import path

#===============================================================================
def main():
    tolerances = parse_tolerances()

    config = ConfigParser.ConfigParser()
    config.read("dashboard.conf")

    ref_tester = "mkrack-pdbg"

    ref_report = sorted(glob("archive/%s/rev_*.txt.gz"%ref_tester))[-1]

    print("Using %s as reference."%ref_report)
    ref_values = parse_report(ref_report, set())
    ref_values_str = parse_report_str(ref_report, set())

    valuesDB = dict()
    testerDB = dict()
    nreports = 0
    outliers = dict()

    for s in config.sections():
        if(config.get(s,"report_type") != "regtest"):
            continue
        outliers[s] = []
        all_reports = glob("archive/%s/rev_*.txt.gz"%s)
        resets = set()
        for fn in sorted(all_reports, reverse=True):
            rev = int(fn.rsplit("rev_",1)[1][:-7])
            if(rev < 15256):
                continue # there were no numeric results back then
            reported_values = parse_report(fn, resets)
            if(reported_values):
                nreports += 1
                c = merge_values(valuesDB, ref_values, tolerances, reported_values)
                outliers[s].append(c)

    print("Parsed %d reports\n"%nreports)

    analyze(valuesDB, ref_values, ref_values_str, tolerances)
    print("\n")

    # seconds reports
    print("%7s %-30s %15s"%("","Tester Name", "Outlier Rate"))
    print("-"*55)
    def get_sortkey(x): return config.getint(x, "sortkey")
    for s in sorted(outliers.keys(), key=get_sortkey):
        print("STATS2: %-30s %15.3f"%(s, np.mean(outliers[s])))

#===============================================================================
def merge_values(valuesDB, ref_values, tolerances, new_values):
    outliers = 0
    for k, v in new_values.items():
        if(k not in valuesDB):
            valuesDB[k] = []
        valuesDB[k].append(v)
        if(k not in ref_values):
            continue
        diff = abs(v - ref_values[k])
        if(diff > tolerances[k]*abs(ref_values[k])):
            outliers += 1
            #print "outlier", diff,ref_values[k], v, tolerances[k]
    return(outliers)


#===============================================================================
def analyze(valuesDB, ref_values, ref_values_str, tolerance):
    n_ok=0; n_shaky=0
    print("%7s %-70s %10s %10s %10s %8s"%("","Test Name", "Max Diff", "Tol", "Out-Rate", "Grade"))
    print("-"*120)
    for k in sorted(ref_values.keys()):
        samples = np.array(valuesDB[k])
        tol = tolerance[k]
        diff = np.abs(samples - ref_values[k]) / abs(ref_values[k])
        max_diff = np.max(diff)
        mean_outsides = np.mean(diff > tol)
        if(max_diff < tol):
            stat = "ok"
            set_ref_value(k, ref_values_str[k])
            n_ok += 1
        else:
            stat ="shaky"
            n_shaky += 1
        print("STATS1: %-70s %10.1e %10.1e %10.4f %8s"%(k, max_diff, tol, mean_outsides, stat))

    print("Ok: %d  Shaky: %d"%(n_ok, n_shaky))

#===============================================================================
def set_ref_value(name, value):
    fn = "../../tests/"+path.dirname(name)+"/TEST_FILES"
    inp = path.basename(name)
    output = []
    content = open(fn).read()
    for line in content.split("\n"):
        if(line.startswith(inp)):
            parts = line.split()
            if(len(parts) == 2):
                parts += ["1.0E-14", value]
            elif(len(parts) == 3):
                parts += [value]
            else:
                assert(False)
            print "NEW:" +"    ".join(parts)
            output.append("    ".join(parts))
        else:
            output.append(line)

    f = open(fn, "w")
    f.write("\n".join(output))
    f.close()

#===============================================================================
def parse_tolerances():
    tolerances = dict()
    tests_root = "../../tests/"
    lines = open(tests_root+"TEST_DIRS").readlines()
    test_dirs = [l.split()[0] for l in lines if l[0]!="#"]
    for d in test_dirs:
        fn = tests_root+d+"/TEST_FILES"
        content = open(fn).read()
        for line in content.strip().split("\n"):
            if(line[0] == "#"):
                continue
            parts = line.split()
            name = d+"/"+parts[0]
            if(len(parts)==2):
                tolerances[name] = 1e-14
            elif(len(parts)==3):
                tolerances[name] = float(parts[2])
            else:
                raise(Exception("Found strange line in: "+fn))

    return(tolerances)

#===============================================================================
def parse_report_str(fn, resets):
    print("Parsing: "+fn)

    values = dict()
    report_txt = gzip.open(fn, 'rb').read()
    m = re.search("\n-+regtesting cp2k-+\n(.*)\n-+ Summary -+\n", report_txt, re.DOTALL)
    if(not m):
        print("Regtests not finished, skipping.")
        return(None)

    main_part = m.group(1)
    curr_dir = None
    for line in main_part.split("\n"):
        if("/UNIT/" in line):
            curr_dir = None # ignore unit-tests
        elif(line.startswith(">>>")):
            curr_dir = line.rsplit("/tests/")[1] + "/"
        elif(line.startswith("<<<")):
            curr_dir = None
        elif(curr_dir):
            parts = line.split()
            if(not parts[0].endswith(".inp")):
                print("Found strange line:\n"+line)
                continue
            if(parts[1]== "RUNTIME" and parts[2]=="FAIL"):
                continue  # ignore crashed tests
            if(parts[1] == "-"):
                continue  # test without numeric check

            test_name = curr_dir+parts[0]
            if(test_name in resets):
                #print("Dropping older values: "+test_name)
                continue # test was reseted, don't add older results
            if(parts[2] == "NEW"):
                resets.add(test_name)
                #print("Found NEW: "+test_name)
            values[test_name] = parts[1]
        else:
            pass # ignore line

    return values

#===============================================================================
def parse_report(fn, resets):
    print("Parsing: "+fn)

    values = dict()
    report_txt = gzip.open(fn, 'rb').read()
    m = re.search("\n-+regtesting cp2k-+\n(.*)\n-+ Summary -+\n", report_txt, re.DOTALL)
    if(not m):
        print("Regtests not finished, skipping.")
        return(None)

    main_part = m.group(1)
    curr_dir = None
    for line in main_part.split("\n"):
        if("/UNIT/" in line):
            curr_dir = None # ignore unit-tests
        elif(line.startswith(">>>")):
            curr_dir = line.rsplit("/tests/")[1] + "/"
        elif(line.startswith("<<<")):
            curr_dir = None
        elif(curr_dir):
            parts = line.split()
            if(not parts[0].endswith(".inp")):
                print("Found strange line:\n"+line)
                continue
            if(parts[1]== "RUNTIME" and parts[2]=="FAIL"):
                continue  # ignore crashed tests
            if(parts[1] == "-"):
                continue  # test without numeric check

            test_name = curr_dir+parts[0]
            if(test_name in resets):
                #print("Dropping older values: "+test_name)
                continue # test was reseted, don't add older results
            if(parts[2] == "NEW"):
                resets.add(test_name)
                #print("Found NEW: "+test_name)
            values[test_name] = float(parts[1])
        else:
            pass # ignore line

    return values

#===============================================================================
main()
#EOF
