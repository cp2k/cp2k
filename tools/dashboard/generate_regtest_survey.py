#!/usr/bin/env python3

# author: Ole Schuett

from urllib.request import urlopen
from datetime import datetime
from glob import glob
import numpy as np
import gzip
import sys
import re

try:
    import configparser
except ImportError:
    import ConfigParser as configparser

# ===============================================================================
def main():
    if len(sys.argv) != 3:
        print("Usage generate_regtest_survey.py <dashboard.conf> <output-dir>")
        sys.exit(1)

    config_fn, outdir = sys.argv[1:]
    assert outdir.endswith("/")

    # parse ../../tests/*/*/TEST_FILES
    test_defs = parse_test_files()

    # parse ../../tests/TEST_TYPES
    test_types = parse_test_types()

    # find eligible testers by parsing latest reports
    tester_names = list()
    tester_values = dict()
    inp_names = set()
    config = configparser.ConfigParser()
    config.read(config_fn)

    def get_sortkey(s):
        return config.getint(s, "sortkey")

    for tname in sorted(config.sections(), key=get_sortkey):
        list_recent = open(outdir + "archive/%s/list_recent.txt" % tname).readlines()
        if not list_recent:
            continue  # list_recent is empty
        latest_report_url = list_recent[0].strip()
        report = parse_report(latest_report_url)
        if report:
            inp_names.update(report.keys())
            tester_values[tname] = report
            tester_names.append(tname)

    # remove outdated inp-names
    inp_names = inp_names.intersection(test_defs.keys())

    # html-header
    output = '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">\n'
    output += "<html><head>\n"
    output += '<meta http-equiv="Content-Type" content="text/html; charset=utf-8">\n'
    output += '<script type="text/javascript" src="https://code.jquery.com/jquery-2.1.4.min.js"></script>\n'
    output += '<script type="text/javascript" src="https://cdnjs.cloudflare.com/ajax/libs/jquery.tablesorter/2.23.2/js/jquery.tablesorter.min.js"></script>\n'
    output += '<script type="text/javascript">\n'
    output += "$(document).ready(function(){\n"
    output += '    $("table").tablesorter();\n'
    output += '    $("table").bind("sortStart",function(){ $(".waitmsg").show(); });\n'
    output += '    $("table").bind("sortEnd",  function(){ $(".waitmsg").hide(); });\n'
    output += "  }\n"
    output += ");\n"
    output += "</script>\n"
    output += '<style type="text/css">\n'
    output += ".nowrap { white-space: nowrap; }\n"
    output += "tr:hover { background-color: #ffff99; }\n"
    output += ".waitmsg {\n"
    output += "  color: red;\n"
    output += "  font: bold 20pt sans-serif;\n"
    output += "  display: none;\n"
    output += "}\n"
    output += "#factbox {\n"
    output += "  display: inline-block;\n"
    output += "  border-radius: 1em;\n"
    output += "  box-shadow: .2em .2em .7em 0 #777;\n"
    output += "  background: #f7f7f0;\n"
    output += "  padding: 1em;\n"
    output += "  margin: 20px;\n"
    output += "}\n"
    output += "#factbox h2 { margin: 0 0 0.5em 0; }\n"
    output += "</style>\n"
    output += "<title>CP2K Regtest Survey</title>\n"
    output += "</head><body>\n"
    output += "<center><h1>CP2K REGTEST SURVEY</h1></center>\n"

    # fun facts
    output += '<div id="factbox"><table>\n'
    ntests = len(test_defs)
    output += "<tr><td>Total number of test-cases</td>"
    output += '<td align="right">%d</td><td align="right">100.0%%</tr>\n' % ntests
    n = len([t for t in test_defs.values() if len(t["flags"]) != 0])
    output += "<tr><td>Tests which require flags</td>"
    output += '<td align="right">%d</td><td align="right">%.1f%%</td></tr>\n' % (
        n,
        n / (0.01 * ntests),
    )
    n = len([t for t in test_defs.values() if int(t["type"]) != 0])
    output += "<tr><td>Numeric tests, ie. type &ne; 0</td>"
    output += '<td align="right">%d</td><td align="right">%.1f%%</td></tr>\n' % (
        n,
        n / (0.01 * ntests),
    )
    n = len([t for t in test_defs.values() if "ref-value" in t])
    output += "<tr><td>Numeric tests with fixed reference</td>"
    output += '<td align="right">%d</td><td align="right">%.1f%%</td></tr>\n' % (
        n,
        n / (0.01 * ntests),
    )
    for i in range(14, 9, -1):
        tol = float("1.e-%d" % i)
        n = len(
            [
                t
                for t in test_defs.values()
                if int(t["type"]) != 0 and float(t["tolerance"]) <= tol
            ]
        )
        output += "<tr><td>Numeric tests with tolerance &le; 10<sup>-%d</sup></td>" % i
        output += '<td align="right">%d</td><td align="right">%.1f%%</td></tr>\n' % (
            n,
            n / (0.01 * ntests),
        )
    output += "</table></div>\n"

    # table-body
    tester_nskipped = dict([(n, 0) for n in tester_names])
    tester_nfailed = dict([(n, 0) for n in tester_names])
    tbody = ""
    for inp in sorted(inp_names):
        # calculate median and MAD
        values = list()
        for tname in tester_names:
            val = tester_values[tname].get(inp, None)
            if val:
                values.append(float(val))
        values = np.array(values)
        median_iterp = np.median(values)  # takes midpoint if len(values)%2==0
        median_idx = (np.abs(values - median_iterp)).argmin()  # find closest point
        median = values[median_idx]
        norm = median + 1.0 * (median == 0.0)
        rel_diff = abs((values - median) / norm)
        mad = np.amax(rel_diff)  # Maximum Absolute Deviation
        outlier = list(rel_diff > test_defs[inp]["tolerance"])

        # output table-row
        tbody += '<tr align="right">\n'
        style = 'bgcolor="#FF9900"' if any(outlier) else ""
        tbody += '<th align="left" %s>%s</th>\n' % (style, inp)
        ttype_num = test_defs[inp]["type"]
        ttype_re = test_types[int(test_defs[inp]["type"])].split("!")[0].strip()
        tbody += '<td title="%s" >%s</td>\n' % (ttype_re, ttype_num)
        tbody += "<td>%.1e</td>\n" % test_defs[inp]["tolerance"]
        tbody += "<td>%.1e</td>\n" % mad
        tol_mad = test_defs[inp]["tolerance"] / max(1e-14, mad)
        if tol_mad < 1.0:
            tbody += '<td bgcolor="#FF9900">%.1e</td>\n' % tol_mad
        elif tol_mad > 10.0:
            tbody += '<td bgcolor="#99FF00">%.1e</td>\n' % tol_mad
        else:
            tbody += "<td>%.1e</td>\n" % tol_mad
        tbody += "<td>%s</td>\n" % test_defs[inp].get("ref-value", "")
        tbody += "<td>%.17g</td>\n" % median
        tbody += "<td>%i</td>\n" % np.sum(outlier)
        for tname in tester_names:
            val = tester_values[tname].get(inp, None)
            if not val:
                tbody += "<td></td>\n"
                tester_nskipped[tname] += 1
            elif outlier.pop(0):
                tbody += '<td bgcolor="#FF9900">%s</td>' % val
                tester_nfailed[tname] += 1
            else:
                tbody += "<td>%s</td>" % val
        tbody += '<th align="left" %s>%s</th>\n' % (style, inp)
        tbody += "</tr>\n"

    # table-header
    theader = '<tr align="center"><th>Name</th><th>Type</th><th>Tolerance</th>'
    theader += '<th><abbr title="Maximum Absolute Deviation">MAD</abbr></th>'
    theader += "<th>Tol. / MAD</th><th>Reference</th><th>Median</th>"
    theader += (
        '<th><abbr title="Failures if reference were at median.">#failed</abbr></th>'
    )
    for tname in tester_names:
        theader += '<th><span class="nowrap">%s</span>' % config.get(tname, "name")
        theader += "<br>#failed: %d" % tester_nfailed[tname]
        theader += "<br>#skipped: %d" % tester_nskipped[tname]
        theader += "</th>\n"
    theader += "<th>Name</th>"
    theader += "</tr>\n"

    # assemble table
    output += '<div class="waitmsg">Sorting, please wait...</div>\n'
    output += "<p>Click on table header to sort by column.</p>\n"
    output += '<table border="1" cellpadding="5">\n'
    output += "<thead>" + theader + "</thead>"
    output += "<tfoot>" + theader + "</tfoot>"
    output += "<tbody>" + tbody + "</tbody>"
    output += "</table>\n"

    # html-footer
    now = datetime.utcnow().replace(microsecond=0)
    output += "<p><small>Page last updated: %s</small></p>\n" % now.isoformat()
    output += "</body></html>"

    # write output file
    fn = outdir + "regtest_survey.html"
    f = open(fn, "w")
    f.write(output)
    f.close()
    print("Wrote: " + fn)


# ===============================================================================
def parse_test_files():
    test_defs = dict()

    tests_root = "../../tests/"
    test_dir_lines = open(tests_root + "TEST_DIRS").readlines()
    for dline in test_dir_lines:
        dline = dline.strip()
        if len(dline) == 0:
            continue
        if dline.startswith("#"):
            continue
        d = dline.split()[0]
        flags = dline.split()[1:]  # flags requiremented by this test_dir
        fn = tests_root + d + "/TEST_FILES"
        content = open(fn).read()
        for line in content.strip().split("\n"):
            if line[0] == "#":
                continue
            parts = line.split()
            name = d + "/" + parts[0]
            entry = {"type": parts[1], "flags": flags}
            if len(parts) == 2:
                entry["tolerance"] = 1.0e-14  # default
            elif len(parts) == 3:
                entry["tolerance"] = float(parts[2])
            elif len(parts) == 4:
                entry["tolerance"] = float(parts[2])
                entry["ref-value"] = parts[3]  # do not parse float
            else:
                raise (Exception("Found strange line in: " + fn))

            test_defs[name] = entry

    return test_defs


# ===============================================================================
def parse_test_types():
    test_types = [None]
    lines = open("../../tests/TEST_TYPES").readlines()
    ntypes = int(lines[0])
    for i in range(1, ntypes + 1):
        test_types.append(lines[i])
    return test_types


# ===============================================================================
def parse_report(report_url):
    print("Parsing: " + report_url)
    data = urlopen(report_url, timeout=5).read()
    report_txt = gzip.decompress(data).decode("utf-8", errors="replace")

    m = re.search(
        "\n-+ ?regtesting cp2k ?-+\n(.*)\n-+ Summary -+\n", report_txt, re.DOTALL
    )
    if not m:
        print("Not a complete regtests report - skipping.")
        return None

    main_part = m.group(1)
    curr_dir = None
    values = dict()
    for line in main_part.split("\n"):
        if "/UNIT/" in line:
            curr_dir = None  # ignore unit-tests
        elif line.startswith(">>>"):
            curr_dir = line.rsplit("/tests/")[1] + "/"
        elif line.startswith("<<<"):
            curr_dir = None
        elif curr_dir:
            if "RUNTIME FAIL" in line:
                continue  # ignore crashed tests
            if "KILLED" in line:
                continue  # ignore timeouted tests
            if "FAILED START" in line:
                continue  # ignore crashed farming run
            parts = line.split()
            if parts[0].rsplit(".", 1)[1] not in ("inp", "restart"):
                print("Found strange line:\n" + line)
                continue
            test_name = curr_dir + parts[0]
            if parts[1] == "-":
                continue  # test without numeric check
            try:
                float(parts[1])  # try parsing float...
                values[test_name] = parts[1]  # ... but pass on the original string
            except:
                pass  # ignore values which can not be parsed
        else:
            pass  # ignore line

    return values


# ===============================================================================
main()
