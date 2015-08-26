#!/usr/bin/python
# -*- coding: utf-8 -*-

# author: Ole Schuett

from datetime import datetime
from glob import glob
import ConfigParser
import gzip
import re

#===============================================================================
def main():
    test_defs = parse_test_files()

    config = ConfigParser.ConfigParser()
    config.read("dashboard.conf")


    # parse latest reports
    latest_values = dict()
    revisions = dict()
    tester_names = list()
    inp_names = set()
    def get_sortkey(s): return config.getint(s, "sortkey")
    for s in sorted(config.sections(), key=get_sortkey):
        if(config.get(s,"report_type") != "regtest"): continue
        tester_names.append(s)
        latest_report_fn = sorted(glob("archive/%s/rev_*.txt.gz"%s))[-1]
        revisions[s] = int(latest_report_fn.rsplit("/rev_")[1][:5])
        report = parse_report(latest_report_fn)
        inp_names.update(report.keys())
        latest_values[s] = report

    # html-header
    output  = '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">\n'
    output += '<html><head>\n'
    output += '<meta http-equiv="Content-Type" content="text/html; charset=utf-8">\n'
    output += '<style type="text/css">\n'
    output += 'tr:hover { background-color: #ffff99; }\n'
    output += '</style>\n'
    output += '<title>CP2K Regtest Survey</title>\n'
    output += '</head><body>\n'
    output += '<center><h1>CP2K REGTEST SURVEY</h1></center>\n'

    # table-header
    output += '<table border="1" cellspacing="3" cellpadding="5">\n'
    output += '<tr align="center"><th>Name</th><th>Type</th><th>Tolerance</th><th>Reference</th>\n'
    for tname in tester_names:
        output += '<th>%s<br>rev %d</th>\n'%(tname, revisions[tname])
    output += '</tr>\n'

    # table-body
    for inp in sorted(inp_names):
        output += '<tr align="right"><th align="left">%s</th>\n'%inp
        output += '<td>%s</td>\n'%test_defs[inp]["type"]
        output += '<td>%s</td>\n'%test_defs[inp]["tolerance"]
        output += '<td>%s</td>\n'%test_defs[inp].get("ref-value", "")
        for tname in tester_names:
            value = latest_values[tname].get(inp, "")
            output += '<td>%s</td>\n'%value
        output += '</tr>\n'
    output += '</table>\n'

    #html-footer
    now = datetime.utcnow().replace(microsecond=0)
    output += '<p><small>Page last updated: %s</small></p>\n'%now.isoformat()
    output += '</body></html>'

    # write output file
    fn = "regtest_survey.html"
    f = open(fn, "w")
    f.write(output)
    f.close()
    print("Wrote: "+fn)

#===============================================================================
def parse_test_files():
    test_defs = dict()

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
            entry = {'type': parts[1]}
            if(len(parts)==2):
                entry['tolerance'] = 1e-14
            elif(len(parts) == 3):
                entry['tolerance'] = float(parts[2])
            elif(len(parts) == 4):
                entry['tolerance'] = float(parts[2])
                entry['ref-value'] = float(parts[3])
            else:
                raise(Exception("Found strange line in: "+fn))

            test_defs[name] = entry

    return(test_defs)

#===============================================================================
def parse_report(fn):
    print("Parsing: "+fn)

    values = dict()
    report_txt = gzip.open(fn, 'rb').read()
    m = re.search("\n-+ ?regtesting cp2k ?-+\n(.*)\n-+ Summary -+\n", report_txt, re.DOTALL)
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
            if(parts[1] == "-"):
                continue  # test without numeric check
            test_name = curr_dir+parts[0]
            values[test_name] = parts[1]
        else:
            pass # ignore line

    return values

#===============================================================================
main()
#EOF
