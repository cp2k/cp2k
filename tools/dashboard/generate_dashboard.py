#!/usr/bin/python
# -*- coding: utf-8 -*-

# Generates the CP2K Dashboard html page
# Inspired by Iain's cp2k_page_update.sh
#
# author: Ole Schuett

import ConfigParser
import sys
import os
import urllib2
import re
import gzip
from datetime import datetime, timedelta
import subprocess
import traceback
from urllib2 import urlopen
from os import path
from pprint import pformat
from xml.dom import minidom
from glob import glob

#===============================================================================
def main():
    if(len(sys.argv) != 4):
        print("Usage update_dashboard.py <config-file> <status-file> <output-dir>")
        sys.exit(1)

    config_fn, status_fn, outdir = sys.argv[1:]
    assert(outdir.endswith("/"))

    config = ConfigParser.ConfigParser()
    config.read(config_fn)

    if(path.exists(status_fn)):
        last_ok = eval(open(status_fn).read())
    else:
        last_ok = dict()

    svn_info = check_output("svn info svn://svn.code.sf.net/p/cp2k/code/trunk".split())
    trunk_revision = int(re.search('Last Changed Rev: (\d+)\n', svn_info).group(1))

    # find latest revision that should have been tested by now
    now = datetime.utcnow().replace(microsecond=0)
    freshness_threshold = now - timedelta(hours=25)
    log = svn_log()
    revs_beyond_threshold = [r for r, d in log if d < freshness_threshold]
    threshold_rev = revs_beyond_threshold[0]
    print "threshold_rev: ", threshold_rev

    output  = html_header(title="CP2K Dashboard")
    output += '<center><table border="1" cellspacing="3" cellpadding="5">\n'
    output += '<tr><th>Name</th><th>Host</th><th>Status</th>'
    output += '<th>Revision</th><th>Summary</th><th>Last OK</th></tr>\n\n'

    def get_sortkey(s):
        return config.getint(s, "sortkey")

    for s in sorted(config.sections(), key=get_sortkey):
        print "Working on: "+s
        name        = config.get(s,"name")
        host        = config.get(s,"host")
        report_type = config.get(s,"report_type")
        report_url  = config.get(s,"report_url")
        info_url    = config.get(s,"info_url") if(config.has_option(s,"info_url")) else None

        report_txt = retrieve_report(report_url)
        report = parse_report(report_txt, report_type)
        if(report['revision']):
            if(report['revision']<threshold_rev):
                report['status'] = "OUTDATED"
            else:
                # store only fresh reports, prevents overwritting archive
                fn = outdir+"archive/%s/rev_%d.txt.gz"%(s,report['revision'])
                write_file(fn, report_txt, gz=True)

        output += '<tr align="center">'
        output += '<td align="left"><a href="archive/%s/index.html">%s</a></td>'%(s, name)
        output += '<td align="left">%s</td>'%host
        output += status_cell(report['status'], report_url)

        #Revision
        if(report.has_key('revision')):
            output += revision_cell(report['revision'], trunk_revision)
            if(report['status'] == "OK"):
                last_ok[s] = report['revision']
        else:
            output += '<td>N/A</td>'

        #Summary
        output += '<td align="left">%s</td>'%report['summary']

        #Last OK
        if(report['status'] != "OK"):
            if(last_ok.has_key(s)):
                output += revision_cell(last_ok[s], trunk_revision)
            else:
                output += '<td>N/A</td>'
        else:
            output += '<td></td>'

        output += '</tr>\n\n'

        # generate archive index
        archive_output = html_header(title=name)
        archive_output += '<p>Go back to <a href="../../index.html">main page</a></p>'
        if(info_url):
            archive_output += '<p>Get <a href="%s">more information</a></p>'%info_url
        archive_output += '<table border="1" cellspacing="3" cellpadding="5">\n'
        archive_output += '<tr><th>Revision</th><th>Status</th><th>Summary</th></tr>\n\n'
        for fn in sorted(glob(outdir+"archive/%s/rev_*.txt.gz"%s), reverse=True):
            archive_output += '<tr align="center">'
            report_txt = gzip.open(fn, 'rb').read()
            report = parse_report(report_txt, report_type)
            archive_output += revision_cell(report['revision'], trunk_revision)
            archive_output += status_cell(report['status'], path.basename(fn)[:-3])
            archive_output += '<td align="left">%s</td>'%report['summary']
            archive_output += '</tr>\n\n'
        archive_output += '</table>\n' + html_footer(now)
        write_file(outdir+"archive/%s/index.html"%s, archive_output)

    output += '</table></center>\n' + html_footer(now)
    write_file(status_fn, pformat(last_ok))
    write_file(outdir+"index.html", output)

#===============================================================================
def retrieve_report(report_url):
    try:
        return urlopen(report_url, timeout=5).read()
    except:
        print traceback.print_exc()
        return None

#===============================================================================
def parse_report(report_txt, report_type):
    if(report_txt==None):
        return( {'status':'UNKOWN', 'summary':'Error while retrieving report.'} )
    try:
        if(report_type == "regtest"):
            return parse_regtest_report(report_txt)
        elif(report_type == "generic"):
            return parse_generic_report(report_txt)
        else:
            raise(Exception("Unkown report_type"))
    except:
        print traceback.print_exc()
        return( {'status':'UNKOWN', 'summary':'Error while parsing report.'} )

#===============================================================================
def html_header(title):
    output  = '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">\n'
    output += '<html><head>\n'
    output += '<meta http-equiv="Content-Type" content="text/html; charset=utf-8">\n'
    output += '<meta http-equiv="refresh" content="200">\n'
    output += '<title>%s</title>\n'%title
    output += '</head><body>\n'
    output += '<center><h1>%s</h1></center>\n'%title.upper()
    return(output)

#===============================================================================
def html_footer(now):
    output  = '<p><small>Page last updated: %s</small></p>\n'%now.isoformat()
    output += '</body></html>'
    return(output)

#===============================================================================
def write_file(fn, content, gz=False):
    d = path.dirname(fn)
    if(len(d) > 0 and not path.exists(d)):
        os.makedirs(d)
        print("Created dir: "+d)
    f = gzip.open(fn, 'wb') if(gz) else open(fn, "w")
    f.write(content)
    f.close()
    print("Wrote: "+fn)

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
def status_cell(status, report_url):
    if(status == "OK"):
        bgcolor = "#00FF00"
    elif(status == "FAILED"):
        bgcolor = "#FF0000"
    else:
        bgcolor = "#d3d3d3"
    return('<td bgcolor="%s"><a href="%s">%s</a></td>'%(bgcolor, report_url, status))

#===============================================================================
def revision_cell(rev, trunk_rev):
    rev_url = "http://sourceforge.net/p/cp2k/code/%d/"%rev
    rev_delta = "(%d)"%(rev - trunk_rev)
    output = '<td align="left"><a href="%s">%s</a> %s</td>'%(rev_url, rev, rev_delta)
    return(output)

#===============================================================================
def parse_regtest_report(report_txt):
    report = dict()
    report['revision']  = int(re.search("(revision|Revision:) (\d+)\.?\n", report_txt).group(2))

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
