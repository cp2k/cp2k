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
import smtplib
from email.mime.text import MIMEText
import re
import gzip
from datetime import datetime, timedelta
import subprocess
import traceback
from urllib import urlencode
from urllib2 import urlopen
from os import path
from pprint import pformat
from xml.dom import minidom
from glob import glob


#===============================================================================
def main():
    if(len(sys.argv) != 5):
        print("Usage update_dashboard.py <config-file> <addressbook> <status-file> <output-dir>")
        sys.exit(1)

    config_fn, abook_fn, status_fn, outdir = sys.argv[1:]
    assert(outdir.endswith("/"))

    assert(path.exists(config_fn))
    config = ConfigParser.ConfigParser()
    config.read(config_fn)

    addressbook = dict([line.split() for line in open(abook_fn).readlines()])

    if(path.exists(status_fn)):
        status = eval(open(status_fn).read())
    else:
        status = dict()

    log = svn_log()
    trunk_revision = log[0]['num']
    log_index = dict([(r['num'], r) for r in log])

    # find latest revision that should have been tested by now
    now = datetime.utcnow().replace(microsecond=0)
    freshness_threshold = now - timedelta(hours=25)
    revs_beyond_threshold = [ r['num'] for r in log if r['date'] < freshness_threshold ]
    threshold_rev = revs_beyond_threshold[0]
    print "threshold_rev: ", threshold_rev

    output  = html_header(title="CP2K Dashboard")
    output += '<center><table border="1" cellspacing="3" cellpadding="5">\n'
    output += '<tr><th>Name</th><th>Host</th><th>Status</th>'
    output += '<th>Revision</th><th>Summary</th><th>Last OK</th><th>Tickets</th></tr>\n\n'

    def get_sortkey(s):
        return config.getint(s, "sortkey")

    for s in sorted(config.sections(), key=get_sortkey):
        print "Working on: "+s
        name        = config.get(s,"name")
        host        = config.get(s,"host")
        report_type = config.get(s,"report_type")
        report_url  = config.get(s,"report_url")
        info_url    = config.get(s,"info_url") if(config.has_option(s,"info_url")) else None
        do_notify   = config.getboolean(s,"notify")

        report_txt = retrieve_report(report_url)
        report = parse_report(report_txt, report_type)

        if(s not in status.keys()):
            status[s] = {'last_ok': None, 'notified': False}

        if(report['status'] == "OK"):
            status[s]['last_ok'] = report['revision']
            status[s]['notified'] = False
        elif(do_notify and not status[s]['notified']):
            send_notification(report, addressbook, status[s]['last_ok'], log_index, name, s)
            status[s]['notified'] = True

        uptodate = report['revision']==trunk_revision # report from latest commit?
        if(report['revision'] != None):
            if(report['revision']<threshold_rev):
                report['status'] = "OUTDATED"
            elif(report['status'] in ("OK", "FAILED")):
                # store only useful and fresh reports, prevents overwriting archive
                fn = outdir+"archive/%s/rev_%d.txt.gz"%(s,report['revision'])
                write_file(fn, report_txt, gz=True)
                check_output("cd %s/archive; tar -cf %s/%s_reports.tar %s/*.txt.gz"%(outdir,s,s,s), shell=True)

        output += '<tr align="center">'
        output += '<td align="left"><a href="archive/%s/index.html">%s</a></td>'%(s, name)
        output += '<td align="left">%s</td>'%host
        output += status_cell(report['status'], report_url, uptodate)

        #Revision
        output += revision_cell(report['revision'], trunk_revision)

        #Summary
        output += '<td align="left">%s</td>'%report['summary']

        #Last OK
        if(report['status'] != "OK"):
            output += revision_cell(status[s]['last_ok'], trunk_revision)
        else:
            output += '<td></td>'

        output += ticket_cell(label=s)

        output += '</tr>\n\n'

        # generate archive index
        archive_output = html_header(title=name)
        archive_output += '<p>Go back to <a href="../../index.html">main page</a></p>'
        if(info_url):
            archive_output += '<p>Get <a href="%s">more information</a></p>'%info_url
        archive_output += '<p>Download <a href="%s_reports.tar">all reports</a></p>'%s
        archive_output += '<table border="1" cellspacing="3" cellpadding="5">\n'
        archive_output += '<tr><th>Revision</th><th>Status</th><th>Summary</th><th>Author</th><th>Commit Message</th></tr>\n\n'

        # read all archived reports
        archive_reports = dict()
        for fn in sorted(glob(outdir+"archive/%s/rev_*.txt.gz"%s), reverse=True):
            report_txt = gzip.open(fn, 'rb').read()
            report = parse_report(report_txt, report_type)
            report['url'] = path.basename(fn)[:-3]
            archive_reports[report['revision']] = report

        # loop over all relevant revisions
        rev_start = max(min(archive_reports.keys()), min(log_index.keys()))
        rev_end = max(log_index.keys())
        for r in range(rev_end, rev_start-1, -1):
            archive_output += '<tr>'
            archive_output += revision_cell(r, trunk_revision)
            if(archive_reports.has_key(r)):
                report = archive_reports[r]
                archive_output += status_cell(report['status'], report['url'])
                archive_output += '<td align="left">%s</td>'%report['summary']
            else:
                archive_output += 2*'<td></td>'
            svn_rev = log_index[r]
            archive_output += '<td align="left">%s</td>'%svn_rev['author']
            archive_output += '<td align="left">%s</td>'%svn_rev['msg'].split("\n")[0]
            archive_output += '</tr>\n\n'
        archive_output += '</table>\n' + html_footer(now)
        write_file(outdir+"archive/%s/index.html"%s, archive_output)

    output += '</table></center>\n' + html_footer(now)
    write_file(status_fn, pformat(status))
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
        return( {'status':'UNKNOWN', 'summary':'Error while retrieving report.', 'revision':None} )
    try:
        if(report_type == "regtest"):
            return parse_regtest_report(report_txt)
        elif(report_type == "generic"):
            return parse_generic_report(report_txt)
        else:
            raise(Exception("Unknown report_type"))
    except:
        print traceback.print_exc()
        return( {'status':'UNKNOWN', 'summary':'Error while parsing report.', 'revision':None} )

#===============================================================================
def send_notification(report, addressbook, last_ok, svn_log, name, s):
    rev_end = report['revision'] if(report['revision']) else max(svn_log.keys())
    authors = set([svn_log[r]['author'] for r in range(last_ok+1, rev_end+1)])
    emails = [addressbook[a] for a in authors]
    print("Sending email to: "+", ".join(emails))

    msg_txt  = "Dear CP2K developer,\n\n"
    msg_txt += "the dashboard has detected a problem that one of your recent commits might have introduced.\n\n"
    msg_txt += "   test name:      %s\n"%name
    msg_txt += "   report state:   %s\n"%report['status']
    msg_txt += "   report summary: %s\n"%report['summary']
    msg_txt += "   last OK rev:    %d\n\n"%last_ok
    msg_txt += "For more information visit:\n"
    msg_txt += "   http://www.cp2k.org/static/dashboard/archive/%s/index.html \n\n"%s
    msg_txt += "Sincerely,\n"
    msg_txt += "  your CP2K Dashboard ;-)\n"

    msg = MIMEText(msg_txt)
    msg['Subject'] = "Problem with "+name
    msg['From']    = "CP2K Dashboard <dashboard@cp2k.org>"
    msg['To']      = ", ".join(emails)

    smtp_conn = smtplib.SMTP('localhost')
    smtp_conn.sendmail(msg['From'], emails, msg.as_string())
    smtp_conn.quit()

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
def svn_log(limit=3000):
    sys.stdout.write("Fetching svn log... ")
    sys.stdout.flush()
    # xml version contains nice UTC timestamp
    cmd = "svn log --limit %d svn://svn.code.sf.net/p/cp2k/code --xml"%limit
    log_xml = check_output(cmd.split())
    dom = minidom.parseString(log_xml)
    revisions = []
    for entry in dom.getElementsByTagName("logentry"):
        rev = dict()
        rev['num'] = int(entry.attributes['revision'].value)
        rev_date_str = entry.getElementsByTagName('date')[0].firstChild.nodeValue
        rev['date'] = datetime.strptime(rev_date_str[:19], '%Y-%m-%dT%H:%M:%S')
        rev['author'] = entry.getElementsByTagName('author')[0].firstChild.nodeValue
        rev['msg'] = entry.getElementsByTagName('msg')[0].firstChild.nodeValue
        revisions.append(rev)
    print("done.")
    return(revisions)

#===============================================================================
def status_cell(status, report_url, uptodate=True):
    if(status == "OK"):
        bgcolor = "#00FF00" if(uptodate) else "#8CE18C"
    elif(status == "FAILED"):
        bgcolor = "#FF0000" if(uptodate) else "#E18C8C"
    else:
        bgcolor = "#d3d3d3"
    return('<td bgcolor="%s"><a href="%s">%s</a></td>'%(bgcolor, report_url, status))

#===============================================================================
def revision_cell(rev, trunk_rev):
    if(rev == None):
        return('<td>N/A</td>')
    rev_url = "http://sourceforge.net/p/cp2k/code/%d/"%rev
    rev_delta = "(%d)"%(rev - trunk_rev)
    output = '<td align="left"><a href="%s">%s</a> %s</td>'%(rev_url, rev, rev_delta)
    return(output)

#===============================================================================
def ticket_cell(label):
    base_url = "https://sourceforge.net/p/cp2k/bugs"
    new_url = base_url+"/new/?" + urlencode({'labels':label})
    query = urlencode({'q':'!status:wont-fix && !status:closed && labels:"%s"'%label})
    feed_url = base_url+"/search_feed/?limit=25&sort=ticket_num_i+asc&" + query
    output = '<td  align="right">'
    try:
        # sometime the http-request to sourceforge times out
        tickets_xml = urlopen(feed_url, timeout=5).read()
        dom = minidom.parseString(tickets_xml)
        for entry in dom.getElementsByTagName("item"):
            title = entry.getElementsByTagName('title')[0].firstChild.nodeValue
            link = entry.getElementsByTagName('link')[0].firstChild.nodeValue
            tid = int(link.strip("/ ").split("/")[-1])
            output += '<a href="%s" title="%s">#%d</a>, '%(link, title, tid)
    except:
        print traceback.print_exc()
        output += "??? "
    output += '<a href="%s"'%new_url
    output += ' style="text-decoration:none;font-weight:bold;font-size:larger;"'
    output += ' title="Create a new Ticket">+</a></td>'
    return(output)

#===============================================================================
def parse_regtest_report(report_txt):
    report = dict()

    m = re.search("svn: E000111: (Can't connect .*):", report_txt)
    if(m):
        report['revision'] = None
        report['status'] = "UNKNOWN"
        report['summary'] = m.group(1)
        return(report)

    report['revision'] = int(re.search("(revision|Revision:) (\d+)\.?\n", report_txt).group(2))

    if("LOCKFILE" in report_txt):
        report['status'] = "UNKNOWN"
        report['summary'] = "Test directory is locked."
        return(report)

    m = re.search("make: .* Error .*", report_txt)
    if(m):
        report['status'] = "FAILED"
        report['summary'] = "Compilation failed."
        return(report)

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
        report['summary'] += "; failed: %d"%runtime_errors
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
def check_output(command, **kwargs):
    p = subprocess.Popen(command, stdout=subprocess.PIPE, **kwargs)
    output = p.communicate()[0]
    assert(p.wait() == 0)
    return(output)

#===============================================================================
main()
#EOF
