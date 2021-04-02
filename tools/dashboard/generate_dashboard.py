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

import matplotlib as mpl

mpl.use("Agg")  # change backend, to run without X11
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

# ======================================================================================
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


# ======================================================================================
def main():
    if len(sys.argv) < 4:
        usage = "<config-file> <status-file> <output-dir> [--send-emails]"
        print(f"Usage: generate_dashboard.py {usage}")
        sys.exit(1)

    config_fn, status_fn, outdir = sys.argv[1:4]
    outdir = outdir.strip("/")
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


# ======================================================================================
def gen_frontpage(config, log, status_fn, outdir):
    if path.exists(status_fn):
        status = eval(open(status_fn, encoding="utf8").read())
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
        do_notify = config.getboolean(s, "notify", fallback=True)
        timeout = config.getint(s, "timeout", fallback=24)

        # find latest commit that should have been tested by now
        threshold = now - timedelta(hours=timeout)
        ages_beyond_threshold = [i for i, c in enumerate(log) if c["date"] < threshold]
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
        output += f'<td align="left"><a href="archive/{s}/index.html">{name}</a></td>'
        output += f'<td align="left">{host}</td>'
        output += status_cell(report["status"], report_url, uptodate)

        # Commit
        output += commit_cell(report["git-sha"], log)

        # Summary
        output += f'<td align="left">{report["summary"]}</td>'

        # Last OK
        if report["status"] != "OK":
            output += commit_cell(status[s]["last_ok"], log)
        else:
            output += "<td></td>"

        output += "</tr>\n\n"

    output += "</table>\n"
    output += '<div id="dummybox"></div></div>\n'  # complete flex-container
    output += html_footer()
    write_file(f"{outdir}/index.html", output)
    write_file(status_fn, pformat(status))


# ======================================================================================
def gen_archive(config, log, outdir):

    for s in config.sections():
        print(f"Working on archive page of: {s}")
        name = config.get(s, "name")
        info_url = config.get(s, "info_url", fallback=None)
        archive_files = glob(f"{outdir}/archive/{s}/commit_*.txt.gz")

        # read cache
        cache_fn = f"{outdir}/archive/{s}/reports.cache"
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
                with gzip.open(fn, "rb") as f:
                    report_txt = f.read().decode("utf-8", errors="replace")
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
                html_row += f'<td align="left">{report["summary"]}</td>'
                url_row = f'https://dashboard.cp2k.org/archive/{s}/{report["url"]}.gz\n'
            else:
                html_row += 2 * "<td></td>"
                url_row = ""
            html_row += f'<td align="left">{html.escape(commit["author-name"])}</td>'
            html_row += f'<td align="left">{html.escape(commit["msg"])}</td>'
            html_row += "</tr>\n\n"
            all_html_rows.append(html_row)
            all_url_rows.append(url_row)

        # generate html pages
        for full_archive in (False, True):
            if full_archive:
                html_out_postfix = "index_full.html"
                urls_out_postfix = "list_full.txt"
                toggle_link = '<p>View <a href="index.html">recent archive</a></p>'
                max_age = None  # output all
            else:
                html_out_postfix = "index.html"
                urls_out_postfix = "list_recent.txt"
                toggle_link = '<p>View <a href="index_full.html">full archive</a></p>'
                max_age = 100

            # generate archive index
            output = html_header(title=name)
            output += '<p>Go back to <a href="../../index.html">main page</a></p>\n'
            if info_url:
                output += f'<p>Get <a href="{info_url}">more information</a></p>\n'
            output += gen_plots(
                archive_reports, log, f"{outdir}/archive/{s}", full_archive
            )
            output += toggle_link
            output += '<table border="1" cellspacing="3" cellpadding="5">\n'
            output += "<tr><th>Commit</th><th>Status</th><th>Summary</th>"
            output += "<th>Author</th><th>Commit Message</th></tr>\n\n"
            output += "".join(all_html_rows[:max_age])
            output += "</table>\n"
            output += toggle_link
            output += html_footer()
            write_file(f"{outdir}/archive/{s}/{html_out_postfix}", output)

            url_list = "".join(all_url_rows[:max_age])
            write_file(f"{outdir}/archive/{s}/{urls_out_postfix}", url_list)


# ======================================================================================
def gen_url_list(config, outdir):
    print("Working on url lists.")
    for postfix in ("list_full.txt", "list_recent.txt"):
        url_list = ""
        for s in config.sections():
            fn = f"{outdir}/archive/{s}/{postfix}"
            if not path.exists(fn):
                continue
            url_list += open(fn, encoding="utf8").read()
        write_file(f"{outdir}/archive/{postfix}", url_list)


# ======================================================================================
def gen_plots(archive_reports, log, outdir, full_archive):

    ordered_shas = [commit["git-sha"] for commit in log]
    ordered_reports = [archive_reports[s] for s in ordered_shas if s in archive_reports]

    # collect plot data
    plots = {}
    for report in ordered_reports:
        for p in report["plots"]:
            if p["name"] not in plots.keys():
                new_plot = {"curves": {}, "title": p["title"], "ylabel": p["ylabel"]}
                plots[p["name"]] = new_plot
        for pp in report["plotpoints"]:
            p = plots[pp["plot"]]
            if pp["name"] not in p["curves"].keys():
                new_curve = {"x": [], "y": [], "yerr": [], "label": pp["label"]}
                p["curves"][pp["name"]] = new_curve
            c = p["curves"][pp["name"]]
            age = log.index[report["git-sha"]]
            c["x"].append(-age)
            c["y"].append(pp["y"])
            c["yerr"].append(pp["yerr"])

    max_age_full = max([-1] + [log.index[sha] for sha in archive_reports.keys()])
    max_age = max_age_full if full_archive else 100

    output = ""
    for pname, p in plots.items():
        print(f"Working on plot: {pname}")
        fn = f"{pname}_full" if full_archive else pname
        title = p["title"]
        make_plot_data(f"{outdir}/{fn}", title, p["curves"], max_age, log)
        make_plot_image(
            f"{outdir}/{fn}", title, p["ylabel"], p["curves"], max_age, full_archive
        )
        output += f'<a href="{fn}.txt"><img src="{fn}.png" alt="{title}"></a>\n'

    return output


# ======================================================================================
def make_plot_data(fn, title, curves, max_age, log):
    output = f"# {title}\n"
    header = ["age", "commit"] + list(curves.keys())
    output += "# " + "\t".join(header) + "\n"
    for age in range(max_age + 1):
        values = []
        for label, c in curves.items():
            if -age in c["x"]:
                i = c["x"].index(-age)
                value = "{:.1f}".format(c["y"][i])
            else:
                value = "?"
            values.append(value.rjust(len(label)))
        if any([v.strip() != "?" for v in values]):
            index_colums = [f"{-age:5d}", log[age]["git-sha"][:6]]
            output += "\t".join(index_colums + values) + "\n"

    write_file(f"{fn}.txt", output)


# ======================================================================================
def make_plot_image(fn, title, ylabel, curves, max_age, full_archive):

    # Setup figure.
    fig = plt.figure(figsize=(12, 4))
    fig.subplots_adjust(bottom=0.18, left=0.06, right=0.70)
    fig.suptitle(title, fontsize=14, fontweight="bold", x=0.4)
    ax = fig.add_subplot(111)
    ax.set_xlabel("Commit Age")
    ax.set_ylabel(ylabel)
    markers = itertools.cycle("os>^*")

    # Draw curves.
    for cname, c in curves.items():
        if full_archive:
            ax.plot(c["x"], c["y"], label=c["label"], linewidth=2)  # less crowded
        else:
            style = dict(marker=next(markers), linewidth=2, markersize=6)
            ax.errorbar(c["x"], c["y"], yerr=c["yerr"], label=c["label"], **style)
    ax.set_xlim(-max_age - 1, 0)
    ax.xaxis.set_minor_locator(AutoMinorLocator())

    # Draw legend.
    legend_style = dict(numpoints=1, fancybox=True, shadow=True, borderaxespad=0.0)
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", **legend_style)

    # Determine y-range such that all curves are visible while ignoring outlayers.
    visibles = []
    for c in curves.values():
        ys = [y for x, y in zip(c["x"], c["y"]) if x >= -max_age]  # visible y-values
        visibles += [ys] if ys else []  # ignore completely invisible curves
    if not visibles:
        print("Warning: Found no visible plot curve.")
    else:
        ymin = min([min(ys) for ys in visibles])  # lowest point from lowest curve
        ymax = max([max(ys) for ys in visibles])  # highest point from highest curve
        if full_archive:
            ax.set_ylim(0.98 * ymin, 1.02 * ymax)
        else:
            ymax2 = max([min(ys) for ys in visibles])  # lowest point from highest curve
            ax.set_ylim(0.98 * ymin, min(1.02 * ymax, 1.3 * ymax2))  # ignore outlayers

    # Save figure to file.
    fig.savefig(fn + ".png")
    plt.close(fig)


# ======================================================================================
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
    msg_txt += "the dashboard has detected a problem"
    msg_txt += " that one of your recent commits might have introduced.\n\n"
    msg_txt += f"   test name:      {name}\n"
    msg_txt += f"   report state:   {report['status']}\n"
    msg_txt += f"   report summary: {report['summary']}\n"
    msg_txt += f"   last OK commit: {last_ok[:7]}\n\n"
    msg_txt += "For more information visit:\n"
    msg_txt += f"   https://dashboard.cp2k.org/archive/{s}/index.html \n\n"
    msg_txt += "Sincerely,\n"
    msg_txt += "  your CP2K Dashboard ;-)\n"

    msg = MIMEText(msg_txt)
    msg["Subject"] = f"Problem with {name}"
    msg["From"] = "CP2K Dashboard <dashboard@cp2k.org>"
    msg["To"] = ", ".join(emails)

    smtp_conn = smtplib.SMTP("localhost")
    smtp_conn.sendmail(msg["From"], emails, msg.as_string())
    smtp_conn.quit()


# ======================================================================================
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
    output += f"<title>{title}</title>\n"
    output += "</head><body>\n"
    output += '<div class="ribbon"><a href="https://cp2k.org/dev:dashboard">Need Help?</a></div>\n'
    output += f"<center><h1>{title.upper()}</h1></center>\n"
    return output


# ======================================================================================
def html_linkbox():
    output = '<div class="sidebox">\n'
    output += "<h2>More...</h2>\n"
    output += '<a href="regtest_survey.html">Regtest Survey</a><br>\n'
    output += '<a href="https://www.cp2k.org/static/coverage/">Test Coverage</a><br>\n'
    output += '<a href="discontinued_tests.html">Discontinued Tests</a><br>\n'
    output += "</div>\n"
    return output


# ======================================================================================
def html_gitbox(log):
    now = datetime.utcnow()
    output = '<div class="sidebox">\n'
    output += "<h2>Recent Commits</h2>\n"
    for commit in log[0:10]:
        url = "https://github.com/cp2k/cp2k/commit/" + commit["git-sha"]
        msg = commit["msg"]
        autor_name = html.escape(commit["author-name"])
        msg = commit["msg"] if len(commit["msg"]) < 28 else commit["msg"][:26] + "..."
        title = html.escape(commit["msg"])
        output += f'<p><a title="{title}" href="{url}">{html.escape(msg)}</a><br>\n'
        delta = now - commit["date"]
        age = delta.days * 24.0 + delta.seconds / 3600.0
        output += "<small>git:" + commit["git-sha"][:7]
        output += f"<br>\n{autor_name} {age:.1f}h ago.</small></p>\n"
    output += "</div>\n"
    return output


# ======================================================================================
def html_footer():
    now = datetime.utcnow().replace(microsecond=0)
    output = f"<p><small>Page last updated: {now.isoformat()}</small></p>\n"
    output += "</body></html>"
    return output


# ======================================================================================
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


# ======================================================================================
def status_cell(status, report_url, uptodate=True):
    if status == "OK":
        bgcolor = "#00FF00" if (uptodate) else "#8CE18C"
    elif status == "FAILED":
        bgcolor = "#FF0000" if (uptodate) else "#E18C8C"
    else:
        bgcolor = "#d3d3d3"
    return f'<td bgcolor="{bgcolor}"><a href="{report_url}">{status}</a></td>'


# ======================================================================================
def commit_cell(git_sha, log):
    if git_sha is None:
        return "<td>N/A</td>"
    idx = log.index[git_sha]
    commit = log[idx]
    git_url = "https://github.com/cp2k/cp2k/commit/" + git_sha
    return f'<td align="left"><a href="{git_url}">{git_sha[:7]}</a> ({-idx})</td>'


# ======================================================================================
def retrieve_report(url):
    try:
        # see if we have a cached entry
        h = hashlib.md5(url.encode("utf8")).hexdigest()
        etag_file = Path(f"/tmp/dashboard_retrieval_cache_{h}.etag")
        data_file = Path(f"/tmp/dashboard_retrieval_cache_{h}.data")
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


# ======================================================================================
def store_report(report, report_txt, section, outdir):
    fn = f'{outdir}/archive/{section}/commit_{report["git-sha"]}.txt.gz'
    write_file(fn, report_txt, gz=True)


# ======================================================================================
def parse_report(report_txt, log):
    if report_txt is None:
        return {
            "status": "UNKNOWN",
            "summary": "Error while retrieving report.",
            "git-sha": None,
        }
    try:
        report = dict()
        git_sha = re.search("(^|\n)CommitSHA: (\w{40})\n", report_txt).group(2)
        report["git-sha"] = git_sha
        report["summary"] = re.findall("(^|[\n\r])Summary: (.+)\n", report_txt)[-1][1]
        report["status"] = re.findall("(^|\n)Status: (.+)\n", report_txt)[-1][1]
        all_plots = re.findall("(^|\n)Plot: (.+)(?=\n)", report_txt)
        report["plots"] = [eval(f"dict({m[1]})") for m in all_plots]
        all_plotpoints = re.findall("(^|\n)PlotPoint: (.+)(?=\n)", report_txt)
        report["plotpoints"] = [eval(f"dict({m[1]})") for m in all_plotpoints]

        # Check that every plot has at least one PlotPoint
        for plot in report["plots"]:
            points = [pp for pp in report["plotpoints"] if pp["plot"] == plot["name"]]
            if not points:
                report["status"] = "FAILED"
                report["summary"] = f"Plot {plot['name']} has no PlotPoints."

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


# ======================================================================================
main()
