#!/usr/bin/python
# -*- coding: utf-8 -*-

# Generates the CP2K Dashboard html page
# Inspired by Iain's cp2k_page_update.sh
#
# author: Ole Schuett

import ConfigParser
import sys
import urllib2
import re
from datetime import datetime, timedelta
import subprocess
import traceback
from urllib2 import urlopen

#===============================================================================
def main():
    if(len(sys.argv) != 2):
        print("Usage update_dashboard.py <config-file>")
        sys.exit(1)

    config_file = sys.argv[1]
    config = ConfigParser.ConfigParser()
    config.read(config_file)

    now = datetime.now().replace(microsecond=0)

    svn_info = check_output("svn info svn://svn.code.sf.net/p/cp2k/code/trunk".split())
    trunk_revision = int(re.search("Revision: (\d+)\n", svn_info).group(1))

    output  = '<html>\n'
    output += '<head><meta http-equiv="refresh" content="200"></head>\n'
    output += '<body>\n'
    output += '<center><h1>CP2K DASHBOARD</h1>\n'
    output += '<table border="1">\n'
    output += '<tr><th>Name</th><th>Host</th><th>Status</th>'
    output += '<th>Revision</th><th>Date</th><th>Summary</th></tr>\n\n'

    for s in sorted(config.sections()):
        print "Working on: "+s
        name        = config.get(s,"name")
        host        = config.get(s,"host")
        report_type = config.get(s,"report_type")
        report_url  = config.get(s,"report_url")
        link_url    = config.get(s,"link_url") if(config.has_option(s,"link_url")) else None

        try:
            report_txt = urlopen(report_url).read()

            if(report_type == "regtest"):
                report = parse_regtest_report(report_txt)
            elif(report_type == "generic"):
                report = parse_generic_report(report_txt)
            else:
                raise(Exception("Unkown report_type"))

            report_age = now - report['date']
            report_fresh = report['revision']==trunk_revision or report_age<timedelta(hours=25)
            if(not report_fresh):
                report['status'] = "OUTDATED"
        except:
            print traceback.print_exc()
            report = dict()
            report['status'] = "UNKOWN"  
            report['summary'] = "Could not retrieve and parse report"

        output += '<tr align="center">'
        if(link_url):
            output += '<td align="left"><a href="%s">%s</a></td>'%(link_url, name)
        else:
            output += '<td align="left">%s</td>'%name
        output += '<td align="left">%s</a></td>'%host
        if(report['status'] == "OK"):
            bgcolor = "#00FF00"
        elif(report['status'] == "FAILED"):
            bgcolor = "#FF0000"
        else:
            bgcolor = "#d3d3d3"
        output += '<td bgcolor=%s><a href="%s">%s</a></td>'%(bgcolor, report_url, report['status'])
        if(report.has_key('date')):
            rev_url = "http://sourceforge.net/p/cp2k/code/%d/"%report['revision']
            rev_delta = "(%d)"%(report['revision'] - trunk_revision)
            output += '<td align="left"><a href="%s">%s</a> %s</td>'%(rev_url, report['revision'], rev_delta)
            output += '<td>%s</td>'%report['date']
        else:
            output += '<td>N/A</td><td>N/A</td>'
        output += '<td align="left">%s</td>'%report['summary']
        output += '</tr>\n\n'

    output += '</table></center>\n'
    output += '<p><small>Page last updated: %s</small></p>\n'%now.isoformat()
    output += '</body></html>'

    open("index.html", "w").write(output)
    print("Wrote index.html")

#===============================================================================
def parse_regtest_report(report_txt):
    report = dict()
    report_date = re.search("\nDate: (.*)\n", report_txt).group(1)
    report_time = re.search("\nTime: (.*)\n", report_txt).group(1)
    report['revision']  = int(re.search("revision (\d+)\.\n", report_txt).group(1))
    report['date'] = datetime.strptime(report_date+' '+report_time, '%Y-%m-%d %H:%M:%S')

    m = re.search("GREPME (\d+) (\d+) (\d+) (\d+) (\d+) (.+)\n", report_txt)
    runtime_errors = int(m.group(1))
    wrong_results  = int(m.group(2))
    correct_tests  = int(m.group(3))
    new_inputs     = int(m.group(4))
    num_tests      = int(m.group(5))
    memory_leaks   = int(m.group(6).replace("X", "0"))

    report['summary'] = "correct: %d / %d"%(correct_tests, num_tests)
    if(new_inputs > 0):
        report['summary'] += "; new: %d"%new_inputs
    if(wrong_results > 0):
        report['summary'] += "; wrong: %d"%wrong_results
    if(runtime_errors > 0):
        report['summary'] += "; crashed: %d"%runtime_errors
    if(memory_leaks > 0):
        report['summary'] += "; memleaks: %d"%memory_leaks

    if(wrong_results>0 or runtime_errors>0 or memory_leaks>0):
        report['status'] = "FAILED"
    else:
        report['status'] = "OK"

    return(report)


#===============================================================================
def parse_generic_report(report_txt):
    report = dict()
    date_str = re.search("(^|\n)Date: (.*)\n", report_txt).group(2)
    report['revision']  = int(re.search("\nRevision: (\d+)\n", report_txt).group(1))
    report['date'] = datetime.strptime(date_str, '%Y-%m-%d %H:%M:%S+00:00')
    report['summary'] = re.search("\nSummary: (.+)\n", report_txt).group(1)
    report['status'] = re.search("\nStatus: (.+)\n", report_txt).group(1)

    return(report)

#===============================================================================
def check_output(command):
    p = subprocess.Popen(command, stdout=subprocess.PIPE)
    output = p.communicate()[0]
    assert(p.wait() == 0)
    return(output)

main()
#EOF
