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
from os import path
from pprint import pformat
from xml.dom import minidom

#===============================================================================
def main():
    if(len(sys.argv) != 4):
        print("Usage update_dashboard.py <config-file> <status-file> <output-file>")
        sys.exit(1)

    config_fn, status_fn, output_fn = sys.argv[1:]

    config = ConfigParser.ConfigParser()
    config.read(config_fn)

    if(path.exists(status_fn)):
        last_change = eval(open(status_fn).read())
    else:
        last_change = dict()

    svn_info = check_output("svn info svn://svn.code.sf.net/p/cp2k/code/trunk".split())
    trunk_revision = int(re.search('Last Changed Rev: (\d+)\n', svn_info).group(1))

    # find latest revision that should have been tested by now
    now = datetime.utcnow().replace(microsecond=0)
    freshness_threshold = now - timedelta(hours=25)
    log = svn_log()
    revs_beyond_threshold = [r for r, d in log if d < freshness_threshold]
    threshold_rev = revs_beyond_threshold[0]
    print "threshold_rev: ", threshold_rev

    output  = '<html>\n'
    output += '<head><meta http-equiv="refresh" content="200"></head>\n'
    output += '<body>\n'
    output += '<center><h1>CP2K DASHBOARD</h1>\n'
    output += '<table border="1">\n'
    output += '<tr><th>Name</th><th>Host</th><th>Status</th>'
    output += '<th>Revision</th><th>Summary</th><th>Since</th></tr>\n\n'

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

            if(report['revision'] < threshold_rev):
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

        #Status
        if(report['status'] == "OK"):
            bgcolor = "#00FF00"
        elif(report['status'] == "FAILED"):
            bgcolor = "#FF0000"
        else:
            bgcolor = "#d3d3d3"
        output += '<td bgcolor=%s><a href="%s">%s</a></td>'%(bgcolor, report_url, report['status'])

        #Revision
        if(report.has_key('revision')):
            output += rev_link(report['revision'], trunk_revision)
            if(not last_change.has_key(s) or last_change[s][0]!=report['status']):
                last_change[s] = (report['status'], report['revision'])
        else:
            output += '<td>N/A</td>'

        #Summary
        output += '<td align="left">%s</td>'%report['summary']

        #Since
        if(last_change.has_key(s)):
            output += rev_link(last_change[s][1], trunk_revision)
        else:
            output += '<td>N/A</td>'

        output += '</tr>\n\n'

    output += '</table></center>\n'
    output += '<p><small>Page last updated: %s</small></p>\n'%now.isoformat()
    output += '</body></html>'

    open(status_fn, "w").write(pformat(last_change))
    print("Wrote "+status_fn)

    open(output_fn, "w").write(output)
    print("Wrote: "+output_fn)

#===============================================================================
def svn_log(limit=100):
    # xml version contains nice UTC timestamp
    cmd = "svn log --limit %d svn://svn.code.sf.net/p/cp2k/code/trunk --xml"%limit
    log_xml = check_output(cmd.split())
    dom = minidom.parseString(log_xml)
    revisions = []
    for entry in dom.getElementsByTagName("logentry"):
        rev_num = int(entry.attributes['revision'].value)
        rev_date_str = entry.getElementsByTagName('date')[0].firstChild.nodeValue
        rev_date = datetime.strptime(rev_date_str[:19], '%Y-%m-%dT%H:%M:%S')
        revisions.append( (rev_num, rev_date) )
    return(revisions)

#===============================================================================
def rev_link(rev, trunk_rev):
    rev_url = "http://sourceforge.net/p/cp2k/code/%d/"%rev
    rev_delta = "(%d)"%(rev - trunk_rev)
    output = '<td align="left"><a href="%s">%s</a> %s</td>'%(rev_url, rev, rev_delta)
    return(output)

#===============================================================================
def parse_regtest_report(report_txt):
    report = dict()
    report['revision']  = int(re.search("revision (\d+)\.\n", report_txt).group(1))

    m = re.search("\nGREPME (\d+) (\d+) (\d+) (\d+) (\d+) (.+)\n", report_txt)
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
    report['revision']  = int(re.search("(^|\n)Revision: (\d+)\n", report_txt).group(2))
    report['summary'] = re.search("(^|\n)Summary: (.+)\n", report_txt).group(2)
    report['status'] = re.search("(^|\n)Status: (.+)\n", report_txt).group(2)

    return(report)

#===============================================================================
def check_output(command):
    p = subprocess.Popen(command, stdout=subprocess.PIPE)
    output = p.communicate()[0]
    assert(p.wait() == 0)
    return(output)

#===============================================================================
main()
#EOF
