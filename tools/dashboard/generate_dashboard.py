#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Generates the CP2K Dashboard html page
# Inspired by Iain's cp2k_page_update.sh
#
# author: Ole Schuett

import sys
import os
import smtplib
from email.mime.text import MIMEText
import re
import gzip
from datetime import datetime, timedelta
from subprocess import check_output
import traceback
from os import path
from pprint import pformat
from glob import glob
import itertools
from urllib.parse import urlencode
from urllib.request import urlopen
import configparser
import pickle

import matplotlib as mpl
mpl.use('Agg')  # change backend, to run without X11
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

#===============================================================================
class GitLog(list):
    def __init__(self):
        cmd = ["git", "log", "--pretty=format:%H%n%ct%n%an%n%ae%n%s%n%b<--seperator-->"]
        outbytes = check_output(cmd)
        output = outbytes.decode("utf-8", errors='replace')
        for entry in output.split("<--seperator-->")[:-1]:
            lines = entry.strip().split("\n")
            commit = dict()
            commit['git-sha'] = lines[0]
            commit['date'] = datetime.fromtimestamp(float(lines[1]))
            commit['author-name'] = lines[2]
            commit['author-email'] = lines[3]
            commit['msg'] = lines[4]
            m = re.match("git-svn-id: svn://svn.code.sf.net/p/cp2k/code/trunk@(\d+) .+", lines[-1])
            if(m):
                commit['svn-rev'] = int(m.group(1))
            self.append(commit)

        # git-log outputs entries from new to old.
        self.index = { c['git-sha']: i for i, c in enumerate(self) }
        self.svn2git = { c['svn-rev']: c['git-sha'] for c in self if 'svn-rev' in c }
        print("done.")

#===============================================================================
def main():
    if(len(sys.argv) != 4):
        print("Usage update_dashboard.py <config-file> <status-file> <output-dir>")
        sys.exit(1)

    config_fn, status_fn, outdir = sys.argv[1:]
    assert(outdir.endswith("/"))
    assert(path.exists(config_fn))

    config = configparser.ConfigParser()
    config.read(config_fn)
    log = GitLog()  # Reads history from local git repo.

    gen_frontpage(config, log, status_fn, outdir)
    gen_archive(config, log, outdir)
    gen_url_list(config, outdir)

#===============================================================================
def gen_frontpage(config, log, status_fn, outdir):
    if(path.exists(status_fn)):
        status = eval(open(status_fn).read())
    else:
        status = dict()

    output  = html_header(title="CP2K Dashboard")
    output += '<div id="flex-container"><div>\n'
    output += html_gitbox(log)
    output += html_linkbox()
    output += '</div>\n'
    output += '<table border="1" cellspacing="3" cellpadding="5">\n'
    output += '<tr><th>Name</th><th>Host</th><th>Status</th>'
    output += '<th>Commit</th><th>Summary</th><th>Last OK</th></tr>\n\n'

    def get_sortkey(s):
        return config.getint(s, "sortkey")

    now = datetime.utcnow().replace(microsecond=0)

    for s in sorted(config.sections(), key=get_sortkey):
        print("Working on summary entry of: "+s)
        name        = config.get(s,"name")
        host        = config.get(s,"host")
        report_type = config.get(s,"report_type")
        report_url  = config.get(s,"report_url")
        do_notify   = config.getboolean(s,"notify") if(config.has_option(s,"notify")) else False
        timeout     = config.getint(s,"timeout") if(config.has_option(s,"timeout")) else 24

        # find latest commit that should have been tested by now
        freshness_threshold = now - timedelta(hours=timeout)
        ages_beyond_threshold = [ i for i, c in enumerate(log) if c['date'] < freshness_threshold ]
        threshold_age = ages_beyond_threshold[0]

        # get and parse report
        report_txt = retrieve_report(report_url)
        report = parse_report(report_txt, report_type, log)

        if(s not in status):
            status[s] = {'last_ok': None, 'notified': False}

        if(report['status'] == "OK"):
            status[s]['last_ok'] = report['git-sha']
            status[s]['notified'] = False
        elif(do_notify and not status[s]['notified']):
            send_notification(report, status[s]['last_ok'], log, name, s)
            status[s]['notified'] = True

        if(report['git-sha']):
            age = log.index[report['git-sha']]
            if(age > threshold_age):
                report['status'] = "OUTDATED"
            elif(report['status'] in ("OK", "FAILED")):
                # store only useful and fresh reports, prevents overwriting archive
                store_report(report, report_txt, s, outdir)

        uptodate = report['git-sha'] == log[0]['git-sha'] # report from latest commit?
        output += '<tr align="center">'
        output += '<td align="left"><a href="archive/%s/index.html">%s</a></td>'%(s, name)
        output += '<td align="left">%s</td>'%host
        output += status_cell(report['status'], report_url, uptodate)

        #Commit
        output += commit_cell(report['git-sha'], log)

        #Summary
        output += '<td align="left">%s</td>'%report['summary']

        #Last OK
        if(report['status'] != "OK"):
            output += commit_cell(status[s]['last_ok'], log)
        else:
            output += '<td></td>'

        output += '</tr>\n\n'

    output += '</table>\n'
    output += '<div id="dummybox"></div></div>\n' # complete flex-container
    output += html_footer()
    write_file(outdir+"index.html", output)
    write_file(status_fn, pformat(status))

#===============================================================================
def gen_archive(config, log, outdir):

    for s in config.sections():
        print("Working on archive page of: "+s)
        name          = config.get(s,"name")
        report_type   = config.get(s,"report_type")
        info_url      = config.get(s,"info_url") if(config.has_option(s,"info_url")) else None
        archive_files = glob(outdir+"archive/%s/rev_*.txt.gz"%s) + \
                        glob(outdir+"archive/%s/commit_*.txt.gz"%s)

        # read cache
        cache_fn = outdir+"archive/%s/reports.cache"%s
        if not path.exists(cache_fn):
            reports_cache = dict()
        else:
            reports_cache = pickle.load(open(cache_fn, "rb"))
            cache_age = path.getmtime(cache_fn)
            # remove outdated cache entries
            reports_cache = {k:v for k,v in reports_cache.items() if path.getmtime(k) < cache_age }

        # read all archived reports
        archive_reports = dict()
        for fn in archive_files:
            if fn in reports_cache:
                report = reports_cache[fn]
            else:
                report_txt = gzip.open(fn, 'rb').read().decode("utf-8", errors='replace')
                report = parse_report(report_txt, report_type, log)
                report['url'] = path.basename(fn)[:-3]
                reports_cache[fn] = report
            sha = report['git-sha']
            if sha is None:
                continue # Skipping report, it's usually a svn commit from a different branch
            assert sha not in archive_reports
            archive_reports[sha] = report

        # write cache
        pickle.dump(reports_cache, open(cache_fn, "wb"))

        # loop over all relevant commits
        all_url_rows = []
        all_html_rows = []
        max_age = 1 + max([log.index[sha] for sha in archive_reports.keys()])
        for commit in log[:max_age]:
            sha = commit['git-sha']
            html_row  = '<tr>'
            html_row += commit_cell(sha, log)
            if(sha in archive_reports):
                report = archive_reports[sha]
                html_row += status_cell(report['status'], report['url'])
                html_row += '<td align="left">%s</td>'%report['summary']
                url_row = "https://dashboard.cp2k.org/archive/%s/%s.gz\n"%(s, report['url'])
            else:
                html_row += 2*'<td></td>'
                url_row = ""
            html_row += '<td align="left">%s</td>'%commit['author-name']
            html_row += '<td align="left">%s</td>'%commit['msg']
            html_row += '</tr>\n\n'
            all_html_rows.append(html_row)
            all_url_rows.append(url_row)

        # generate html pages
        for full_archive in (False, True):
            if(full_archive):
                html_out_postfix = "index_full.html"
                urls_out_postfix = "list_full.txt"
                other_index_link = '<p>View <a href="index.html">recent archive</a></p>'
                max_age = None  # output all
            else:
                html_out_postfix = "index.html"
                urls_out_postfix = "list_recent.txt"
                other_index_link = '<p>View <a href="index_full.html">full archive</a></p>'
                max_age = 100

            # generate archive index
            output = html_header(title=name)
            output += '<p>Go back to <a href="../../index.html">main page</a></p>\n'
            if(info_url):
                output += '<p>Get <a href="%s">more information</a></p>\n'%info_url
            output += gen_plots(archive_reports, log, outdir+"archive/"+s+"/", full_archive)
            output += other_index_link
            output += '<table border="1" cellspacing="3" cellpadding="5">\n'
            output += '<tr><th>Commit</th><th>Status</th><th>Summary</th><th>Author</th><th>Commit Message</th></tr>\n\n'
            output += "".join(all_html_rows[:max_age])
            output += '</table>\n'
            output += other_index_link
            output += html_footer()
            html_out_fn = outdir+"archive/%s/%s"%(s,html_out_postfix)
            write_file(html_out_fn, output)

            url_list = "".join(all_url_rows[:max_age])
            urls_out_fn = outdir+"archive/%s/%s"%(s,urls_out_postfix)
            write_file(urls_out_fn, url_list)

#===============================================================================
def gen_url_list(config, outdir):
    print("Working on url lists.")
    for postfix in ("list_full.txt", "list_recent.txt"):
        url_list = ""
        for s in config.sections():
            fn = outdir+"archive/%s/%s"%(s, postfix)
            if(not path.exists(fn)):
                continue
            url_list += open(fn).read()
        write_file(outdir+"archive/" + postfix, url_list)

#===============================================================================
def gen_plots(archive_reports, log, outdir, full_archive):

    ordered_shas = [commit['git-sha'] for commit in log]
    ordered_reports = [archive_reports[sha] for sha in ordered_shas if sha in archive_reports]

    # collect plot data
    plots = {}
    for report in ordered_reports:
        for p in report['plots']:
            if(p['name'] not in plots.keys()):
                plots[p['name']] = {'curves':{}}
            plots[p['name']]['title'] = p['title'] # update title
            plots[p['name']]['ylabel'] = p['ylabel'] # update label
        for pp in report['plotpoints']:
            p = plots[pp['plot']]
            if(pp['name'] not in p['curves'].keys()):
                p['curves'][pp['name']] = {'x':[], 'y':[], 'yerr':[]}
            c = p['curves'][pp['name']]
            age = log.index[report['git-sha']]
            c['x'].append(-age)
            c['y'].append(pp['y'])
            c['yerr'].append(pp['yerr'])
            c['label'] = pp['label'] # update label

    # write raw data
    tags = sorted([(pname, cname) for pname, p in plots.items() for cname in p['curves'].keys()])
    if(tags):
        raw_output = "# %6s    %40s" % ("age", "commit")
        for pname, cname in tags:
            raw_output += "   %18s   %22s"%(pname+"/"+cname,pname+"/"+cname+"_err")
        raw_output += "\n"
        for report in reversed(ordered_reports):
            age = log.index[report['git-sha']]
            raw_output +=  "%8d    %40s"%(-age, report['git-sha'])
            for pname, cname in tags:
                pp = [pp for pp in report['plotpoints'] if(pp['plot']==pname and pp['name']==cname)]
                assert(len(pp)<=1)
                if(pp):
                    raw_output += "   %18f   %22f"%(pp[0]['y'],pp[0]['yerr'])
                else:
                    raw_output += "   %18s   %22s"%("?","?")
            raw_output += "\n"
        write_file(outdir+"plot_data.txt", raw_output)

    # create png images
    if (full_archive):
        fig_ext = "_full.png"
        max_age = max([log.index[sha] for sha in archive_reports.keys()])
    else:
        fig_ext = ".png"
        max_age = 100

    for pname, p in plots.items():
        print("Working on plot: "+pname)
        fig = plt.figure(figsize=(12,4))
        fig.subplots_adjust(bottom=0.18, left=0.06, right=0.70)
        fig.suptitle(p['title'], fontsize=14, fontweight='bold', x=0.4)
        ax = fig.add_subplot(111)
        ax.set_xlabel('Commit Age')
        ax.set_ylabel(p['ylabel'])
        markers = itertools.cycle('os>^*')
        for cname in sorted(p['curves'].keys()):
            c = p['curves'][cname]
            if(full_archive):
                ax.plot(c['x'], c['y'], label=c['label'], linewidth=2) # less crowded
            else:
                ax.errorbar(c['x'], c['y'], yerr=c['yerr'], label=c['label'],
                            marker=next(markers), linewidth=2, markersize=6)
        ax.set_xlim(-max_age-1, 0)
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left',
                  numpoints=1, fancybox=True, shadow=True, borderaxespad=0.0)
        visibles = [[y for x,y in zip(c['x'],c['y']) if x>=-max_age] for c in p['curves'].values()] # visible y-values
        ymin = min([min(ys) for ys in visibles if ys]) # lowest point from lowest curve
        ymax = max([max(ys) for ys in visibles if ys]) # highest point from highest curve
        if(full_archive):
            ax.set_ylim(0.98*ymin, 1.02*ymax)
        else:
            ymax2 = max([min(ys) for ys in visibles if ys]) # lowest point from highest curve
            ax.set_ylim(0.98*ymin, min(1.02*ymax, 1.3*ymax2))  # protect against outlayers
        fig.savefig(outdir+pname+fig_ext)
        plt.close(fig)

    # write html output
    html_output = ""
    for pname in sorted(plots.keys()):
        html_output += '<a href="plot_data.txt"><img src="%s" alt="%s"></a>\n'%(pname+fig_ext, plots[pname]['title'])
    return(html_output)

#===============================================================================
def send_notification(report, last_ok, log, name, s):
    idx_end = log.index[report['git-sha']] if(report['git-sha']) else 0
    idx_last_ok = log.index[last_ok]
    if(idx_end == idx_last_ok): return # probably a flapping tester
    emails = set([log[i]['author-email'] for i in range(idx_end, idx_last_ok)])
    print("Sending email to: "+", ".join(emails))

    msg_txt  = "Dear CP2K developer,\n\n"
    msg_txt += "the dashboard has detected a problem that one of your recent commits might have introduced.\n\n"
    msg_txt += "   test name:      %s\n"%name
    msg_txt += "   report state:   %s\n"%report['status']
    msg_txt += "   report summary: %s\n"%report['summary']
    msg_txt += "   last OK commit: %s\n\n"%last_ok[:7]
    msg_txt += "For more information visit:\n"
    msg_txt += "   https://dashboard.cp2k.org/archive/%s/index.html \n\n"%s
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
    output += '<link rel="icon" type="image/x-icon" href="data:image/x-icon;base64,AAABAAEAEBAQAAAAAAAoAQAAFgAAACgAAAAQAAAAIAAAAAEABAAAAAAAgAAAAAAAAAAAAAAAEAAAAAAAAAAAAAAAAAD/AAmRCQAAb/8AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAiIgAAAAAAACIiAAAAAAAAIiIAAAAAAAAAAAAAAAAAADMzAAAAAAAAMzMAAAAAAAAzMwAAAAAAAAAAAAAAAAAAEREAAAAAAAAREQAAAAAAABERAAAAAAAAAAAAAAD+fwAA/n8AAPw/AAD4HwAA+B8AAPgfAAD4HwAA+B8AAPgfAAD4HwAA+B8AAPgfAAD4HwAA+B8AAPgfAAD4HwAA">\n'
    output += '<style type="text/css">\n'
    output += '.ribbon {\n'
    output += '  overflow: hidden;\n'
    output += '  position: absolute;\n'
    output += '  right:0px;\n'
    output += '  top: 0px;\n'
    output += '  width: 200px;\n'
    output += '  height: 200px;\n'
    output += '}\n'
    output += '.ribbon a {\n'
    output += '  position: relative;\n'
    output += '  white-space: nowrap;\n'
    output += '  background-color: #a00;\n'
    output += '  border: 1px solid #faa;\n'
    output += '  color: #fff;\n'
    output += '  display: block;\n'
    output += '  font: bold 11pt sans-serif;\n'
    output += '  padding: 7px;\n'
    output += '  top: 35px;\n'
    output += '  right: 10px;\n'
    output += '  width: 300px;\n'
    output += '  text-align: center;\n'
    output += '  text-decoration: none;\n'
    output += '  transform: rotate(45deg);\n'
    output += '  box-shadow: 0 0 10px #888;\n'
    output += '}\n'
    output += '#flex-container {\n'
    output += '  display: -webkit-flex; /* Safari */\n'
    output += '  display: flex;\n'
    output += '  -webkit-flex-flow: row wrap-reverse; /* Safari */\n'
    output += '  flex-flow:         row wrap-reverse;\n'
    output += '  -webkit-justify-content: space-around; /* Safari */\n'
    output += '  justify-content:         space-around;\n'
    output += '  -webkit-align-items: flex-end; /* Safari */\n'
    output += '  align-items:         flex-end;\n'
    output += '}\n'
    output += '.sidebox {\n'
    output += '  width: 15em;\n'
    output += '  border-radius: 1em;\n'
    output += '  box-shadow: .2em .2em .7em 0 #777;\n'
    output += '  background: #f7f7f0;\n'
    output += '  padding: 1em;\n'
    output += '  margin: 40px 20px;\n'
    output += '}\n'
    output += '.sidebox h2 {\n'
    output += '  margin: 0 0 0.5em 0;\n'
    output += '}\n'
    output += '.sidebox p {\n'
    output += '  margin: 0.5em;\n'
    output += '}\n'
    output += '#dummybox {\n'
    output += '  width: 15em;\n'
    output += '}\n'
    output += '</style>\n'
    output += '<title>%s</title>\n'%title
    output += '</head><body>\n'
    output += '<div class="ribbon"><a href="https://cp2k.org/dev:dashboard">Need Help?</a></div>\n'
    output += '<center><h1>%s</h1></center>\n'%title.upper()
    return(output)

#===============================================================================
def html_linkbox():
    output  = '<div class="sidebox">\n'
    output += '<h2>More...</h2>\n'
    output += '<a href="regtest_survey.html">Regtest Survey</a><br>\n'
    output += '<a href="https://www.cp2k.org/static/coverage/">Test Coverage</a><br>\n'
    output += '<a href="discontinued_tests.html">Discontinued Tests</a><br>\n'
    output += '</div>\n'
    return(output)

#===============================================================================
def html_gitbox(log):
    now = datetime.utcnow()
    output  = '<div class="sidebox">\n'
    output += '<h2>Recent Commits</h2>\n'
    for commit in log[0:10]:
        url = "https://github.com/cp2k/cp2k/commit/" + commit['git-sha']
        msg = commit['msg']
        if(len(msg) > 27):
            msg = msg[:26] + "..."
        output += '<p><a title="%s" href="%s">%s</a><br>\n'%(commit['msg'], url, msg)
        delta = now - commit['date']
        age = delta.days*24.0 + delta.seconds/3600.0
        output += '<small>git:' + commit['git-sha'][:7]
        if 'svn-rev' in commit:
            output += ' / svn:%d'%commit['svn-rev']
        output += '<br>\n%s %.1fh ago.</small></p>\n'%(commit['author-name'], age)
    output += '</div>\n'
    return(output)

#===============================================================================
def html_footer():
    now = datetime.utcnow().replace(microsecond=0)
    output  = '<p><small>Page last updated: %s</small></p>\n'%now.isoformat()
    output += '</body></html>'
    return(output)

#===============================================================================
def write_file(fn, content, gz=False):
    d = path.dirname(fn)
    if(len(d) > 0 and not path.exists(d)):
        os.makedirs(d)
        print("Created dir: "+d)
    if(path.exists(fn)):
        old_bytes = gzip.open(fn, 'rb').read() if(gz) else open(fn, 'rb').read()
        old_content = old_bytes.decode('utf-8', errors='replace')
        if(old_content == content):
            print("File did not change: "+fn)
            return
    f = gzip.open(fn, 'wb') if(gz) else open(fn, "wb")
    f.write(content.encode("utf-8"))
    f.close()
    print("Wrote: "+fn)

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
def commit_cell(git_sha, log):
    if(git_sha is None):
        return('<td>N/A</td>')
    idx = log.index[git_sha]
    commit = log[idx]
    git_url = "https://github.com/cp2k/cp2k/commit/" + git_sha
    output = '<td align="left"><a href="%s">%s</a>'%(git_url, git_sha[:7])
    if 'svn-rev' in commit:
        svn_rev = commit['svn-rev']
        svn_url = "https://sourceforge.net/p/cp2k/code/%d/"%svn_rev
        output += ' / <a href="%s">%d</a>'%(svn_url, svn_rev)
    output += ' (%d)</td>'%(-idx)
    return(output)

#===============================================================================
def retrieve_report(report_url):
    try:
        return urlopen(report_url, timeout=5).read().decode("utf-8", errors='replace')
    except:
        print(traceback.print_exc())
        return None

#===============================================================================
def store_report(report, report_txt, section, outdir):
    if ('svn-rev' in report):
        fn = outdir+"archive/%s/rev_%d.txt.gz"%(section, report['svn-rev'])
    else:
        fn = outdir+"archive/%s/commit_%s.txt.gz"%(section, report['git-sha'])
    write_file(fn, report_txt, gz=True)

#===============================================================================
def parse_report(report_txt, report_type, log):
    if(report_txt is None):
        return( {'status':'UNKNOWN', 'summary':'Error while retrieving report.', 'git-sha':None} )
    try:
        if(report_type == "regtest"):
            report = parse_regtest_report(report_txt)
        elif(report_type == "generic"):
            report = parse_generic_report(report_txt)
        else:
            raise(Exception("Unknown report_type"))

        if ('svn-rev' in report):
            assert 'git-sha' not in report
            rev = report['svn-rev']
            if rev not in log.svn2git:
                return( {'status':'UNKNOWN', 'summary':'Could not convert svn revision to git commit.', 'git-sha':None} )
            report['git-sha'] = log.svn2git[rev]

        report.update(parse_plots(report_txt))
        return(report)
    except:
        print(traceback.print_exc())
        return( {'status':'UNKNOWN', 'summary':'Error while parsing report.', 'git-sha':None} )

#===============================================================================
def parse_regtest_report(report_txt):
    m = re.search("svn: .*(Can't connect .*):", report_txt)
    if(m):
        return({'status':'UNKNOWN', 'summary':m.group(1)})

    report = dict()
    m = re.search("(^|\n)CommitSHA: (\w{40})\n", report_txt)
    if (m):
        report['git-sha'] = m.group(2)
    else:
        m = re.search("(revision|Revision:) (\d+)\.?\n", report_txt)
        report['svn-rev'] = int(m.group(2))

    if("LOCKFILE" in report_txt):
        report['status'] = "UNKNOWN"
        report['summary'] = "Test directory is locked."
        return(report)

    m = re.findall("\nGREPME (\d+) (\d+) (\d+) (\d+) (\d+) (.+)\n", report_txt)
    if(not m and re.search("make: .* Error .*", report_txt)):
        report['status'] = "FAILED"
        report['summary'] = "Compilation failed."
        return(report)

    runtime_errors = int(m[-1][0])
    wrong_results  = int(m[-1][1])
    correct_tests  = int(m[-1][2])
    new_inputs     = int(m[-1][3])
    num_tests      = int(m[-1][4])
    memory_leaks   = int(m[-1][5].replace("X", "0"))

    report['summary'] = "correct: %d / %d"%(correct_tests, num_tests)
    if(new_inputs > 0):
        report['summary'] += "; new: %d"%new_inputs
    if(wrong_results > 0):
        report['summary'] += "; wrong: %d"%wrong_results
    if(runtime_errors > 0):
        report['summary'] += "; failed: %d"%runtime_errors
    if(memory_leaks > 0):
        report['summary'] += "; memleaks: %d"%memory_leaks

    runtimes = [float(m) for m in re.findall("\nRegtest took (.+) seconds.\n", report_txt)]
    report['summary'] += "; %.0fmin"%(sum(runtimes)/60.0)

    if(wrong_results>0 or runtime_errors>0 or memory_leaks>0):
        report['status'] = "FAILED"
    else:
        report['status'] = "OK"

    return(report)

#===============================================================================
def parse_generic_report(report_txt):
    m = re.search("svn: .*(Can't connect .*):", report_txt)
    if(m):
        return({'status':'UNKNOWN', 'summary':m.group(1)})
    report = dict()

    m = re.search("(^|\n)CommitSHA: (\w{40})\n", report_txt)
    if (m):
        report['git-sha'] = m.group(2)
    else:
        m = re.search("(^|\n)Revision: (\d+)\n", report_txt)
        report['svn-rev'] = int(m.group(2))

    report['summary'] = re.findall("(^|\n)Summary: (.+)\n", report_txt)[-1][1]
    report['status'] = re.findall("(^|\n)Status: (.+)\n", report_txt)[-1][1]
    return(report)

#===============================================================================
def parse_plots(report_txt):
    plots = [eval("dict(%s)"%m[1]) for m in re.findall("(^|\n)Plot: (.+)(?=\n)", report_txt)]
    plotpoints = [eval("dict(%s)"%m[1]) for m in re.findall("(^|\n)PlotPoint: (.+)(?=\n)", report_txt)]
    return({'plots': plots, 'plotpoints': plotpoints})

#===============================================================================
if(len(sys.argv)==2 and sys.argv[-1]=="--selftest"):
    pass #TODO implement selftest
else:
    main()
#EOF
