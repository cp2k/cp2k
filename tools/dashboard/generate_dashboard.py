#!/usr/bin/env python3

# Generates the CP2K Dashboard html page
# Inspired by Iain's cp2k_page_update.sh
#
# author: Ole Schuett

import sys
import os
import smtplib
from email.mime.text import MIMEText
import html
import re
import gzip
from datetime import datetime, timedelta
from subprocess import check_output
import traceback
from os import path
from pprint import pformat
from glob import glob
import itertools
import configparser
import pickle
from pathlib import Path
import requests
import hashlib
from collections import OrderedDict

import matplotlib as mpl

mpl.use("Agg")  # change backend, to run without X11
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

# ===============================================================================
class GitLog(list):
    def __init__(self):
        cmd = ["git", "log", "--pretty=format:%H%n%ct%n%an%n%ae%n%s%n%b<--separator-->"]
        outbytes = check_output(cmd)
        output = outbytes.decode("utf-8", errors="replace")
        for entry in output.split("<--separator-->")[:-1]:
            lines = entry.strip().split("\n")
            commit = dict()
            commit["git-sha"] = lines[0]
            commit["date"] = datetime.fromtimestamp(float(lines[1]))
            commit["author-name"] = lines[2]
            commit["author-email"] = lines[3]
            commit["msg"] = lines[4]
            self.append(commit)

        # git-log outputs entries from new to old.
        self.index = {c["git-sha"]: i for i, c in enumerate(self)}


# ===============================================================================
def main():
    if len(sys.argv) < 4:
        print(
            "Usage update_dashboard.py <config-file> <status-file> <output-dir> [--send-emails]"
        )
        sys.exit(1)

    config_fn, status_fn, outdir = sys.argv[1:4]
    assert outdir.endswith("/")
    assert path.exists(config_fn)

    global send_emails
    if len(sys.argv) == 5:
        assert sys.argv[4] == "--send-emails"
        send_emails = True
    else:
        send_emails = False

    config = configparser.ConfigParser()
    config.read(config_fn)
    log = GitLog()  # Reads history from local git repo.

    gen_frontpage(config, log, status_fn, outdir)
    gen_archive(config, log, outdir)
    gen_url_list(config, outdir)


# ===============================================================================
def gen_frontpage(config, log, status_fn, outdir):
    if path.exists(status_fn):
        status = eval(open(status_fn).read())
    else:
        status = dict()

    output = html_header(title="CP2K Dashboard")
    output += '<div id="flex-container"><div>\n'
    output += html_gitbox(log)
    output += html_linkbox()
    output += "</div>\n"
    output += '<table border="1" cellspacing="3" cellpadding="5">\n'
    output += "<tr><th>Name</th><th>Host</th><th>Status</th>"
    output += "<th>Commit</th><th>Summary</th><th>Last OK</th></tr>\n\n"

    def get_sortkey(s):
        return config.getint(s, "sortkey")

    now = datetime.utcnow().replace(microsecond=0)

    for s in sorted(config.sections(), key=get_sortkey):
        print("Working on summary entry of: " + s)
        name = config.get(s, "name")
        host = config.get(s, "host")
        report_url = config.get(s, "report_url")
        do_notify = (
            config.getboolean(s, "notify") if (config.has_option(s, "notify")) else True
        )
        timeout = (
            config.getint(s, "timeout") if (config.has_option(s, "timeout")) else 24
        )

        # find latest commit that should have been tested by now
        freshness_threshold = now - timedelta(hours=timeout)
        ages_beyond_threshold = [
            i for i, c in enumerate(log) if c["date"] < freshness_threshold
        ]
        threshold_age = ages_beyond_threshold[0]

        # get and parse report
        report_txt = retrieve_report(report_url)
        report = parse_report(report_txt, log)

        if s not in status:
            status[s] = {"last_ok": None, "notified": False}

        if report["status"] == "OK":
            status[s]["last_ok"] = report["git-sha"]
            status[s]["notified"] = False
        elif do_notify and not status[s]["notified"]:
            send_notification(report, status[s]["last_ok"], log, name, s)
            status[s]["notified"] = True

        if report["git-sha"]:
            age = log.index[report["git-sha"]]
            if age > threshold_age:
                report["status"] = "OUTDATED"
            elif report["status"] in ("OK", "FAILED"):
                # store only useful and fresh reports, prevents overwriting archive
                store_report(report, report_txt, s, outdir)

        uptodate = report["git-sha"] == log[0]["git-sha"]  # report from latest commit?
        output += '<tr align="center">'
        output += '<td align="left"><a href="archive/%s/index.html">%s</a></td>' % (
            s,
            name,
        )
        output += '<td align="left">%s</td>' % host
        output += status_cell(report["status"], report_url, uptodate)

        # Commit
        output += commit_cell(report["git-sha"], log)

        # Summary
        output += '<td align="left">%s</td>' % report["summary"]

        # Last OK
        if report["status"] != "OK":
            output += commit_cell(status[s]["last_ok"], log)
        else:
            output += "<td></td>"

        output += "</tr>\n\n"

    output += "</table>\n"
    output += '<div id="dummybox"></div></div>\n'  # complete flex-container
    output += html_footer()
    write_file(outdir + "index.html", output)
    write_file(status_fn, pformat(status))


# ===============================================================================
def gen_archive(config, log, outdir):

    for s in config.sections():
        print("Working on archive page of: " + s)
        name = config.get(s, "name")
        info_url = (
            config.get(s, "info_url") if (config.has_option(s, "info_url")) else None
        )
        archive_files = glob(outdir + "archive/%s/rev_*.txt.gz" % s) + glob(
            outdir + "archive/%s/commit_*.txt.gz" % s
        )

        # read cache
        cache_fn = outdir + "archive/%s/reports.cache" % s
        if not path.exists(cache_fn):
            reports_cache = dict()
        else:
            reports_cache = pickle.load(open(cache_fn, "rb"))
            cache_age = path.getmtime(cache_fn)
            # remove outdated cache entries
            reports_cache = {
                k: v for k, v in reports_cache.items() if path.getmtime(k) < cache_age
            }

        # read all archived reports
        archive_reports = dict()
        for fn in archive_files:
            if fn in reports_cache:
                report = reports_cache[fn]
            else:
                report_txt = (
                    gzip.open(fn, "rb").read().decode("utf-8", errors="replace")
                )
                report = parse_report(report_txt, log)
                report["url"] = path.basename(fn)[:-3]
                reports_cache[fn] = report
            sha = report["git-sha"]
            assert sha not in archive_reports
            archive_reports[sha] = report

        # write cache
        if reports_cache:
            pickle.dump(reports_cache, open(cache_fn, "wb"))

        # loop over all relevant commits
        all_url_rows = []
        all_html_rows = []
        max_age = 1 + max([-1] + [log.index[sha] for sha in archive_reports.keys()])
        for commit in log[:max_age]:
            sha = commit["git-sha"]
            html_row = "<tr>"
            html_row += commit_cell(sha, log)
            if sha in archive_reports:
                report = archive_reports[sha]
                html_row += status_cell(report["status"], report["url"])
                html_row += '<td align="left">%s</td>' % report["summary"]
                url_row = "https://dashboard.cp2k.org/archive/%s/%s.gz\n" % (
                    s,
                    report["url"],
                )
            else:
                html_row += 2 * "<td></td>"
                url_row = ""
            html_row += '<td align="left">%s</td>' % html.escape(commit["author-name"])
            html_row += '<td align="left">%s</td>' % html.escape(commit["msg"])
            html_row += "</tr>\n\n"
            all_html_rows.append(html_row)
            all_url_rows.append(url_row)

        # generate html pages
        for full_archive in (False, True):
            if full_archive:
                html_out_postfix = "index_full.html"
                urls_out_postfix = "list_full.txt"
                other_index_link = '<p>View <a href="index.html">recent archive</a></p>'
                max_age = None  # output all
            else:
                html_out_postfix = "index.html"
                urls_out_postfix = "list_recent.txt"
                other_index_link = (
                    '<p>View <a href="index_full.html">full archive</a></p>'
                )
                max_age = 100

            # generate archive index
            output = html_header(title=name)
            output += '<p>Go back to <a href="../../index.html">main page</a></p>\n'
            if info_url:
                output += '<p>Get <a href="%s">more information</a></p>\n' % info_url
            output += gen_plots(
                archive_reports, log, outdir + "archive/" + s + "/", full_archive
            )
            output += other_index_link
            output += '<table border="1" cellspacing="3" cellpadding="5">\n'
            output += "<tr><th>Commit</th><th>Status</th><th>Summary</th><th>Author</th><th>Commit Message</th></tr>\n\n"
            output += "".join(all_html_rows[:max_age])
            output += "</table>\n"
            output += other_index_link
            output += html_footer()
            html_out_fn = outdir + "archive/%s/%s" % (s, html_out_postfix)
            write_file(html_out_fn, output)

            url_list = "".join(all_url_rows[:max_age])
            urls_out_fn = outdir + "archive/%s/%s" % (s, urls_out_postfix)
            write_file(urls_out_fn, url_list)


# ===============================================================================
def gen_url_list(config, outdir):
    print("Working on url lists.")
    for postfix in ("list_full.txt", "list_recent.txt"):
        url_list = ""
        for s in config.sections():
            fn = outdir + "archive/%s/%s" % (s, postfix)
            if not path.exists(fn):
                continue
            url_list += open(fn).read()
        write_file(outdir + "archive/" + postfix, url_list)


# ===============================================================================
def gen_plots(archive_reports, log, outdir, full_archive):

    ordered_shas = [commit["git-sha"] for commit in log]
    ordered_reports = [
        archive_reports[sha] for sha in ordered_shas if sha in archive_reports
    ]

    # collect plot data
    plots = OrderedDict()
    for report in ordered_reports:
        for p in report["plots"]:
            if p["name"] not in plots.keys():
                plots[p["name"]] = {
                    "curves": OrderedDict(),
                    "title": p["title"],
                    "ylabel": p["ylabel"],
                }
        for pp in report["plotpoints"]:
            p = plots[pp["plot"]]
            if pp["name"] not in p["curves"].keys():
                p["curves"][pp["name"]] = {
                    "x": [],
                    "y": [],
                    "yerr": [],
                    "label": pp["label"],
                }
            c = p["curves"][pp["name"]]
            age = log.index[report["git-sha"]]
            c["x"].append(-age)
            c["y"].append(pp["y"])
            c["yerr"].append(pp["yerr"])

    # write raw data
    tags = sorted(
        [(pname, cname) for pname, p in plots.items() for cname in p["curves"].keys()]
    )
    if tags:
        raw_output = "# %6s    %40s" % ("age", "commit")
        for pname, cname in tags:
            raw_output += "   %18s   %22s" % (
                pname + "/" + cname,
                pname + "/" + cname + "_err",
            )
        raw_output += "\n"
        for report in reversed(ordered_reports):
            age = log.index[report["git-sha"]]
            raw_output += "%8d    %40s" % (-age, report["git-sha"])
            for pname, cname in tags:
                pp = [
                    pp
                    for pp in report["plotpoints"]
                    if (pp["plot"] == pname and pp["name"] == cname)
                ]
                if len(pp) > 1:
                    print("Warning: Found redundant plot points.")
                if pp:
                    raw_output += "   %18f   %22f" % (pp[-1]["y"], pp[-1]["yerr"])
                else:
                    raw_output += "   %18s   %22s" % ("?", "?")
            raw_output += "\n"
        write_file(outdir + "plot_data.txt", raw_output)

    # create png images
    if full_archive:
        fig_ext = "_full.png"
        max_age = max([-1] + [log.index[sha] for sha in archive_reports.keys()])
    else:
        fig_ext = ".png"
        max_age = 100

    for pname, p in plots.items():
        print("Working on plot: " + pname)
        fig = plt.figure(figsize=(12, 4))
        fig.subplots_adjust(bottom=0.18, left=0.06, right=0.70)
        fig.suptitle(p["title"], fontsize=14, fontweight="bold", x=0.4)
        ax = fig.add_subplot(111)
        ax.set_xlabel("Commit Age")
        ax.set_ylabel(p["ylabel"])
        markers = itertools.cycle("os>^*")
        for cname in p["curves"].keys():
            c = p["curves"][cname]
            if full_archive:
                ax.plot(c["x"], c["y"], label=c["label"], linewidth=2)  # less crowded
            else:
                ax.errorbar(
                    c["x"],
                    c["y"],
                    yerr=c["yerr"],
                    label=c["label"],
                    marker=next(markers),
                    linewidth=2,
                    markersize=6,
                )
        ax.set_xlim(-max_age - 1, 0)
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.legend(
            bbox_to_anchor=(1.01, 1),
            loc="upper left",
            numpoints=1,
            fancybox=True,
            shadow=True,
            borderaxespad=0.0,
        )
        visibles = [
            [y for x, y in zip(c["x"], c["y"]) if x >= -max_age]
            for c in p["curves"].values()
        ]  # visible y-values
        visibles = [ys for ys in visibles if ys]  # remove completely invisible curves
        if not visibles:
            print("Warning: Found no visible plot curve.")
        else:
            ymin = min([min(ys) for ys in visibles])  # lowest point from lowest curve
            ymax = max([max(ys) for ys in visibles])  # highest point from highest curve
            if full_archive:
                ax.set_ylim(0.98 * ymin, 1.02 * ymax)
            else:
                ymax2 = max(
                    [min(ys) for ys in visibles]
                )  # lowest point from highest curve
                ax.set_ylim(
                    0.98 * ymin, min(1.02 * ymax, 1.3 * ymax2)
                )  # protect against outlayers
        fig.savefig(outdir + pname + fig_ext)
        plt.close(fig)

    # write html output
    html_output = ""
    for pname in sorted(plots.keys()):
        html_output += '<a href="plot_data.txt"><img src="%s" alt="%s"></a>\n' % (
            pname + fig_ext,
            plots[pname]["title"],
        )
    return html_output


# ===============================================================================
def send_notification(report, last_ok, log, name, s):
    idx_end = log.index[report["git-sha"]] if (report["git-sha"]) else 0
    if not last_ok:
        return  # we don't know when this started
    idx_last_ok = log.index[last_ok]
    if idx_end == idx_last_ok:
        return  # probably a flapping tester
    emails = set([log[i]["author-email"] for i in range(idx_end, idx_last_ok)])
    emails = [e for e in emails if "noreply" not in e]
    emails_str = ", ".join(emails)
    if not emails:
        return  # no author emails found
    if len(emails) > 3:
        print("Spam protection, found more than three authors: " + emails_str)
        return
    if not send_emails:
        print("Email sending disabled, would otherwise send to: " + emails_str)
        return

    print("Sending email to: " + emails_str)

    msg_txt = "Dear CP2K developer,\n\n"
    msg_txt += "the dashboard has detected a problem that one of your recent commits might have introduced.\n\n"
    msg_txt += "   test name:      %s\n" % name
    msg_txt += "   report state:   %s\n" % report["status"]
    msg_txt += "   report summary: %s\n" % report["summary"]
    msg_txt += "   last OK commit: %s\n\n" % last_ok[:7]
    msg_txt += "For more information visit:\n"
    msg_txt += "   https://dashboard.cp2k.org/archive/%s/index.html \n\n" % s
    msg_txt += "Sincerely,\n"
    msg_txt += "  your CP2K Dashboard ;-)\n"

    msg = MIMEText(msg_txt)
    msg["Subject"] = "Problem with " + name
    msg["From"] = "CP2K Dashboard <dashboard@cp2k.org>"
    msg["To"] = ", ".join(emails)

    smtp_conn = smtplib.SMTP("localhost")
    smtp_conn.sendmail(msg["From"], emails, msg.as_string())
    smtp_conn.quit()


# ===============================================================================
def html_header(title):
    output = '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">\n'
    output += "<html><head>\n"
    output += '<meta http-equiv="Content-Type" content="text/html; charset=utf-8">\n'
    output += '<meta http-equiv="refresh" content="200">\n'
    output += '<link rel="icon" type="image/x-icon" href="data:image/x-icon;base64,AAABAAEAEBAQAAAAAAAoAQAAFgAAACgAAAAQAAAAIAAAAAEABAAAAAAAgAAAAAAAAAAAAAAAEAAAAAAAAAAAAAAAAAD/AAmRCQAAb/8AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAiIgAAAAAAACIiAAAAAAAAIiIAAAAAAAAAAAAAAAAAADMzAAAAAAAAMzMAAAAAAAAzMwAAAAAAAAAAAAAAAAAAEREAAAAAAAAREQAAAAAAABERAAAAAAAAAAAAAAD+fwAA/n8AAPw/AAD4HwAA+B8AAPgfAAD4HwAA+B8AAPgfAAD4HwAA+B8AAPgfAAD4HwAA+B8AAPgfAAD4HwAA">\n'
    output += '<style type="text/css">\n'
    output += ".ribbon {\n"
    output += "  overflow: hidden;\n"
    output += "  position: absolute;\n"
    output += "  right:0px;\n"
    output += "  top: 0px;\n"
    output += "  width: 200px;\n"
    output += "  height: 200px;\n"
    output += "}\n"
    output += ".ribbon a {\n"
    output += "  position: relative;\n"
    output += "  white-space: nowrap;\n"
    output += "  background-color: #a00;\n"
    output += "  border: 1px solid #faa;\n"
    output += "  color: #fff;\n"
    output += "  display: block;\n"
    output += "  font: bold 11pt sans-serif;\n"
    output += "  padding: 7px;\n"
    output += "  top: 35px;\n"
    output += "  right: 10px;\n"
    output += "  width: 300px;\n"
    output += "  text-align: center;\n"
    output += "  text-decoration: none;\n"
    output += "  transform: rotate(45deg);\n"
    output += "  box-shadow: 0 0 10px #888;\n"
    output += "}\n"
    output += "#flex-container {\n"
    output += "  display: -webkit-flex; /* Safari */\n"
    output += "  display: flex;\n"
    output += "  -webkit-flex-flow: row wrap-reverse; /* Safari */\n"
    output += "  flex-flow:         row wrap-reverse;\n"
    output += "  -webkit-justify-content: space-around; /* Safari */\n"
    output += "  justify-content:         space-around;\n"
    output += "  -webkit-align-items: flex-end; /* Safari */\n"
    output += "  align-items:         flex-end;\n"
    output += "}\n"
    output += ".sidebox {\n"
    output += "  width: 15em;\n"
    output += "  border-radius: 1em;\n"
    output += "  box-shadow: .2em .2em .7em 0 #777;\n"
    output += "  background: #f7f7f0;\n"
    output += "  padding: 1em;\n"
    output += "  margin: 40px 20px;\n"
    output += "}\n"
    output += ".sidebox h2 {\n"
    output += "  margin: 0 0 0.5em 0;\n"
    output += "}\n"
    output += ".sidebox p {\n"
    output += "  margin: 0.5em;\n"
    output += "}\n"
    output += "#dummybox {\n"
    output += "  width: 15em;\n"
    output += "}\n"
    output += "</style>\n"
    output += "<title>%s</title>\n" % title
    output += "</head><body>\n"
    output += '<div class="ribbon"><a href="https://cp2k.org/dev:dashboard">Need Help?</a></div>\n'
    output += "<center><h1>%s</h1></center>\n" % title.upper()
    return output


# ===============================================================================
def html_linkbox():
    output = '<div class="sidebox">\n'
    output += "<h2>More...</h2>\n"
    output += '<a href="regtest_survey.html">Regtest Survey</a><br>\n'
    output += '<a href="https://www.cp2k.org/static/coverage/">Test Coverage</a><br>\n'
    output += '<a href="discontinued_tests.html">Discontinued Tests</a><br>\n'
    output += "</div>\n"
    return output


# ===============================================================================
def html_gitbox(log):
    now = datetime.utcnow()
    output = '<div class="sidebox">\n'
    output += "<h2>Recent Commits</h2>\n"
    for commit in log[0:10]:
        url = "https://github.com/cp2k/cp2k/commit/" + commit["git-sha"]
        msg = commit["msg"]
        if len(msg) > 27:
            msg = msg[:26] + "..."
        output += '<p><a title="%s" href="%s">%s</a><br>\n' % (
            html.escape(commit["msg"]),
            url,
            html.escape(msg),
        )
        delta = now - commit["date"]
        age = delta.days * 24.0 + delta.seconds / 3600.0
        output += "<small>git:" + commit["git-sha"][:7]
        output += "<br>\n%s %.1fh ago.</small></p>\n" % (
            html.escape(commit["author-name"]),
            age,
        )
    output += "</div>\n"
    return output


# ===============================================================================
def html_footer():
    now = datetime.utcnow().replace(microsecond=0)
    output = "<p><small>Page last updated: %s</small></p>\n" % now.isoformat()
    output += "</body></html>"
    return output


# ===============================================================================
def write_file(fn, content, gz=False):
    d = path.dirname(fn)
    if len(d) > 0 and not path.exists(d):
        os.makedirs(d)
        print("Created dir: " + d)
    if path.exists(fn):
        old_bytes = gzip.open(fn, "rb").read() if (gz) else open(fn, "rb").read()
        old_content = old_bytes.decode("utf-8", errors="replace")
        if old_content == content:
            print("File did not change: " + fn)
            return
    f = gzip.open(fn, "wb") if (gz) else open(fn, "wb")
    f.write(content.encode("utf-8"))
    f.close()
    print("Wrote: " + fn)


# ===============================================================================
def status_cell(status, report_url, uptodate=True):
    if status == "OK":
        bgcolor = "#00FF00" if (uptodate) else "#8CE18C"
    elif status == "FAILED":
        bgcolor = "#FF0000" if (uptodate) else "#E18C8C"
    else:
        bgcolor = "#d3d3d3"
    return '<td bgcolor="%s"><a href="%s">%s</a></td>' % (bgcolor, report_url, status)


# ===============================================================================
def commit_cell(git_sha, log):
    if git_sha is None:
        return "<td>N/A</td>"
    idx = log.index[git_sha]
    commit = log[idx]
    git_url = "https://github.com/cp2k/cp2k/commit/" + git_sha
    output = '<td align="left"><a href="%s">%s</a>' % (git_url, git_sha[:7])
    output += " (%d)</td>" % (-idx)
    return output


# ===============================================================================
def retrieve_report(url):
    try:
        # see if we have a cached entry
        h = hashlib.md5(url.encode("utf8")).hexdigest()
        etag_file = Path("/tmp/dashboard_retrieval_cache_" + h + ".etag")
        data_file = Path("/tmp/dashboard_retrieval_cache_" + h + ".data")
        etag = etag_file.read_text() if etag_file.exists() else ""

        # make conditional http request
        r = requests.get(url, headers={"If-None-Match": etag}, timeout=5)
        r.raise_for_status()
        if r.status_code == 304:  # Not Modified - cache hit
            return data_file.read_text()

        # check report size
        report_size = int(r.headers["Content-Length"])
        assert report_size < 3 * 1024 * 1024  # 3 MB

        # cache miss - store response
        if "ETag" in r.headers:
            data_file.write_text(r.text)
            etag_file.write_text(r.headers["ETag"])
        return r.text

    except:
        print(traceback.print_exc())
        return None


# ===============================================================================
def store_report(report, report_txt, section, outdir):
    fn = outdir + "archive/%s/commit_%s.txt.gz" % (section, report["git-sha"])
    write_file(fn, report_txt, gz=True)


# ===============================================================================
def parse_report(report_txt, log):
    if report_txt is None:
        return {
            "status": "UNKNOWN",
            "summary": "Error while retrieving report.",
            "git-sha": None,
        }
    try:
        report = dict()
        report["git-sha"] = re.search("(^|\n)CommitSHA: (\w{40})\n", report_txt).group(
            2
        )
        report["summary"] = re.findall("(^|\n)Summary: (.+)\n", report_txt)[-1][1]
        report["status"] = re.findall("(^|\n)Status: (.+)\n", report_txt)[-1][1]
        report["plots"] = [
            eval("dict(%s)" % m[1])
            for m in re.findall("(^|\n)Plot: (.+)(?=\n)", report_txt)
        ]
        report["plotpoints"] = [
            eval("dict(%s)" % m[1])
            for m in re.findall("(^|\n)PlotPoint: (.+)(?=\n)", report_txt)
        ]

        # Check that every plot has at least one PlotPoint
        for plot in report["plots"]:
            points = [pp for pp in report["plotpoints"] if pp["plot"] == plot["name"]]
            if not points:
                report["status"] = "FAILED"
                report["summary"] = 'Plot "%s" has no PlotPoints.' % plot["name"]

        # Check that CommitSHA belongs to the master branch.
        if report["git-sha"] not in log.index:
            report["git-sha"] = None
            report["status"] = "FAILED"
            report["summary"] = "Unknown CommitSHA."

        return report
    except:
        print(traceback.print_exc())
        return {
            "status": "UNKNOWN",
            "summary": "Error while parsing report.",
            "git-sha": None,
        }


# ===============================================================================
main()
