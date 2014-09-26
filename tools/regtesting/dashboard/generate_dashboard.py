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
    trunk_rev = int(re.search("Revision: (\d+)\n", svn_info).group(1))

    output  = '<html>\n'
    output += '<head><meta http-equiv="refresh" content="200"></head>\n'
    output += '<body>\n'
    output += '<center><h1>CP2K DASHBOARD</h1>\n'
    output += '<table border="1">\n'
    output += '<tr><th>Automated regtester name</th><th>Configuration</th>'
    output += '<th>Number of tests</th><th>Ok</th><th>Wrong</th><th>New</th>'
    output += '<th>Runtime fail</th><th>Memory leaks</th><th>Revision</th><th>Timestamp</th></tr>\n\n'

    for s in config.sections():
        print "Working on: "+s
        org = config.get(s, "org")
        base_url = config.get(s, "base_url")
        report_url = base_url + "/regtest-0"
        test_ok = False
        test_fresh = False
        try:
            report = urlopen(report_url).read()
            test_date = re.search("\nDate: (.*)\n", report).group(1)
            test_time = re.search("\nTime: (.*)\n", report).group(1)
            test_rev  = int(re.search("revision (\d+)\.\n", report).group(1))
            test_datetime = datetime.strptime(test_date+' '+test_time, '%Y-%m-%d %H:%M:%S')
            test_age = now - test_datetime
            test_fresh = test_rev==trunk_rev or test_age<timedelta(hours=25)

            m = re.search("GREPME (\d+) (\d+) (\d+) (\d+) (\d+) (.+)\n", report)
            if(m):
              runtime_errors = int(m.group(1))
              wrong_results  = int(m.group(2))
              correct_tests  = int(m.group(3))
              new_inputs     = int(m.group(4))
              num_tests      = int(m.group(5))
              memory_leaks   = m.group(6).replace("X", "N/A")
              test_ok        = True
        except:
            traceback.print_exc()
            pass

        arch_url = config.get(s, "arch_url")
        arch_label = config.get(s, "arch_label")
        output += '<tr align="center">'
        output += '<td align="left"><a href="%s">%s</a></td>'%(base_url, org)
        output += '<td align="left"><a href="%s">%s</a></td>'%(arch_url, arch_label)

        if(test_ok and test_fresh):
            output += '<td>%s</td>'%num_tests
            bgcolor = "#00FF00" if(correct_tests==num_tests) else "#FF0000"
            output += '<td bgcolor=%s>%s</td>'%(bgcolor, correct_tests)
            output += '<td>%s</td>'%wrong_results
            output += '<td>%s</td>'%new_inputs
            output += '<td>%s</td>'%runtime_errors
            output += '<td>' if(memory_leaks in ("0", "N/A")) else '<td bgcolor=#FF0000>'
            output +=     '%s</td>'%memory_leaks
            output += '<td>%d</td>'%test_rev
            output += '<td>%s %s</td>'%(test_date, test_time)
        elif(test_ok):
            output += '<td colspan=8 bgcolor=#d3d3d3>Out of date</td>'
        else:
            output += '<td colspan=8 bgcolor=#d3d3d3>N/A</td>'

        output += '</tr>\n\n'

    output += '</table></center>\n'
    output += '<p><small>Page last updated: %s</small></p>\n'%now.isoformat()
    output += '</body></html>'

    open("index.html", "w").write(output)
    print("Wrote index.html")

#===============================================================================
def check_output(command):
    p = subprocess.Popen(command, stdout=subprocess.PIPE)
    output = p.communicate()[0]
    assert(p.wait() == 0)
    return(output)

main()
#EOF
