#!/usr/bin/env python3

# Generates the CP2K Dashboard html page
# Inspired by Iain's cp2k_page_update.sh
#
# author: Ole Schuett

from configparser import ConfigParser
from datetime import datetime, timedelta
from email.mime.text import MIMEText
from os import path
from pathlib import Path
from pprint import pformat
from subprocess import check_output
import argparse
import dataclasses
import gzip
import hashlib
import html
import itertools
import os
import pickle
import re
import requests
import smtplib
import sys
import traceback
from typing import Any, List, Dict, cast, Optional, ValuesView, NewType

try:
    from typing import Literal
except:
    from typing_extensions import Literal  # type: ignore

import matplotlib as mpl  # type: ignore

mpl.use("Agg")  # change backend, to run without X11
import matplotlib.pyplot as plt  # type: ignore
from matplotlib.ticker import AutoMinorLocator  # type: ignore

# ======================================================================================
GitSha = NewType("GitSha", str)
ReportStatus = Literal["OK", "FAILED", "OUTDATED", "UNKNOWN"]

# ======================================================================================
class GitLog:
    def __init__(self) -> None:
        self.commits: List[Commit] = []
        cmd = ["git", "log", "--pretty=format:%H%n%ct%n%an%n%ae%n%s%n%b<--separator-->"]
        outbytes = check_output(cmd)
        output = outbytes.decode("utf-8", errors="replace")
        for entry in output.split("<--separator-->")[:-1]:
            lines = entry.strip().split("\n")
            commit = Commit(
                sha=GitSha(lines[0]),
                date=datetime.fromtimestamp(float(lines[1])),
                author_name=lines[2],
                author_email=lines[3],
                message=lines[4],
            )
            self.commits.append(commit)

        # git-log outputs entries from new to old.
        self.index_by_sha = {c.sha: i for i, c in enumerate(self.commits)}


# ======================================================================================
@dataclasses.dataclass
class Commit:
    sha: GitSha
    author_name: str
    author_email: str
    date: datetime
    message: str


# ======================================================================================
@dataclasses.dataclass
class PlotCurve:
    label: str
    x: List[float] = dataclasses.field(default_factory=list)
    y: List[float] = dataclasses.field(default_factory=list)
    yerr: List[float] = dataclasses.field(default_factory=list)


# ======================================================================================
class Plot:
    def __init__(self, name: str, title: str, ylabel: str, **_: Any):
        self.name = name
        self.title = title
        self.ylabel = ylabel


# ======================================================================================
class AggregatedPlot:
    def __init__(self, plot: Plot):
        self.name = plot.name
        self.title = plot.title
        self.ylabel = plot.ylabel
        self.curves_by_name: Dict[str, PlotCurve] = {}

    @property
    def curves(self) -> ValuesView[PlotCurve]:
        return self.curves_by_name.values()


# ======================================================================================
class PlotPoint:
    def __init__(
        self, plot: str, name: str, label: str, y: float, yerr: float, **_: Any
    ):
        self.plot_name = plot
        self.curve_name = name
        self.curve_label = label
        self.y = y
        self.yerr = yerr


# ======================================================================================
@dataclasses.dataclass
class Report:
    status: ReportStatus
    summary: str
    sha: Optional[GitSha] = None
    plots: List[Plot] = dataclasses.field(default_factory=list)
    plot_points: List[PlotPoint] = dataclasses.field(default_factory=list)
    archive_path: Optional[str] = None


# ======================================================================================
@dataclasses.dataclass
class Status:
    notified: bool = False
    last_ok: Optional[GitSha] = None


# ======================================================================================
def main() -> None:
    parser = argparse.ArgumentParser(description="Generates the CP2K Dashboard.")
    parser.add_argument("config_file", type=Path)
    parser.add_argument("status_file", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--send-emails", action="store_true")
    args = parser.parse_args()

    assert args.config_file.exists()
    config = ConfigParser()
    config.read(args.config_file)
    log = GitLog()  # Reads history from local git repo.

    gen_frontpage(config, log, args.status_file, args.output_dir, args.send_emails)
    gen_archive(config, log, args.output_dir)
    gen_url_list(config, args.output_dir)


# ======================================================================================
def gen_frontpage(
    config: ConfigParser, log: GitLog, status_fn: Path, outdir: Path, send_emails: bool
) -> None:

    status: Dict[str, Status] = {}
    if path.exists(status_fn):
        with open(status_fn, "rb") as f:
            status = pickle.load(f)

    output = html_header(title="CP2K Dashboard")
    output += '<div id="flex-container"><div>\n'
    output += html_gitbox(log)
    output += html_linkbox()
    output += "</div>\n"
    output += '<table border="1" cellspacing="3" cellpadding="5">\n'
    output += "<tr><th>Name</th><th>Host</th><th>Status</th>"
    output += "<th>Commit</th><th>Summary</th><th>Last OK</th></tr>\n\n"

    now = datetime.utcnow().replace(microsecond=0)

    for s in sorted(config.sections(), key=lambda s: config.getint(s, "sortkey")):
        print("Working on summary entry of: " + s)
        name = config.get(s, "name")
        host = config.get(s, "host")
        report_url = config.get(s, "report_url")
        do_notify = config.getboolean(s, "notify", fallback=True)
        timeout = config.getint(s, "timeout", fallback=24)

        # find latest commit that should have been tested by now
        threshold = now - timedelta(hours=timeout)
        threshold_age = next(i for i, c in enumerate(log.commits) if c.date < threshold)

        # get and parse report
        report_txt = retrieve_report(report_url)
        report = parse_report(report_txt, log)

        if s not in status:
            status[s] = Status()

        if report.status == "OK":
            status[s].last_ok = report.sha
            status[s].notified = False
        elif do_notify and not status[s].notified:
            send_notification(report, status[s].last_ok, log, name, s, send_emails)
            status[s].notified = True

        if report.sha:
            age = log.index_by_sha[report.sha]
            if age > threshold_age:
                report.status = "OUTDATED"
            elif report.status in ("OK", "FAILED"):
                # store only useful and fresh reports, prevents overwriting archive
                assert report_txt
                store_report(report, report_txt, s, outdir)

        up_to_date = report.sha == log.commits[0].sha  # report from latest commit?
        output += '<tr align="center">'
        output += f'<td align="left"><a href="archive/{s}/index.html">{name}</a></td>'
        output += f'<td align="left">{host}</td>'
        output += status_cell(report.status, report_url, up_to_date)

        # Commit
        output += commit_cell(report.sha, log)

        # Summary
        output += f'<td align="left">{report.summary}</td>'

        # Last OK
        if report.status != "OK":
            output += commit_cell(status[s].last_ok, log)
        else:
            output += "<td></td>"

        output += "</tr>\n\n"

    output += "</table>\n"
    output += '<div id="dummybox"></div></div>\n'  # complete flex-container
    output += html_footer()
    write_file(outdir / "index.html", output)
    with open(status_fn, "wb") as f:
        pickle.dump(status, f)


# ======================================================================================
def gen_archive(config: ConfigParser, log: GitLog, outdir: Path) -> None:

    for s in config.sections():
        print(f"Working on archive page of: {s}")
        name = config.get(s, "name")
        info_url = config.get(s, "info_url", fallback=None)
        archive_dir = outdir / f"archive/{s}"
        archive_files = archive_dir.glob("commit_*.txt.gz")

        # read cache
        reports_cache: Dict[Path, Report] = {}
        cache_fn = f"{archive_dir}/reports.cache"
        if path.exists(cache_fn):
            with open(cache_fn, "rb") as f:
                reports_cache = pickle.load(f)
            cache_age = path.getmtime(cache_fn)
            # remove outdated cache entries
            reports_cache = {
                k: v for k, v in reports_cache.items() if path.getmtime(k) < cache_age
            }

        # read all archived reports
        archive_reports: Dict[GitSha, Report] = {}
        for fn in archive_files:
            report: Report
            if fn in reports_cache:
                report = reports_cache[fn]
            else:
                with open(fn, "rb") as f:
                    content = gzip.decompress(f.read())
                report = parse_report(content.decode("utf-8", errors="replace"), log)
                report.archive_path = path.basename(fn)[:-3]
                reports_cache[fn] = report
            assert report.sha and report.sha not in archive_reports
            archive_reports[report.sha] = report

        # write cache
        if reports_cache:
            with open(cache_fn, "wb") as f:
                pickle.dump(reports_cache, f)

        # loop over all relevant commits
        all_url_rows = []
        all_html_rows = []
        max_age_full = max([-1] + [log.index_by_sha[sha] for sha in archive_reports])
        for commit in log.commits[: max_age_full + 1]:
            sha = commit.sha
            html_row = "<tr>"
            html_row += commit_cell(sha, log)
            if sha in archive_reports:
                report = archive_reports[sha]
                assert report.archive_path
                html_row += status_cell(report.status, report.archive_path)
                html_row += f'<td align="left">{report.summary}</td>'
                archive_base_url = "https://dashboard.cp2k.org/archive"
                url_row = f"{archive_base_url}/{s}/{report.archive_path}.gz\n"
            else:
                html_row += 2 * "<td></td>"
                url_row = ""
            html_row += f'<td align="left">{html.escape(commit.author_name)}</td>'
            html_row += f'<td align="left">{html.escape(commit.message)}</td>'
            html_row += "</tr>\n\n"
            all_html_rows.append(html_row)
            all_url_rows.append(url_row)

        # generate html pages
        for full_archive in (False, True):
            if full_archive:
                html_out_postfix = "index_full.html"
                urls_out_postfix = "list_full.txt"
                toggle_link = '<p>View <a href="index.html">recent archive</a></p>'
                max_age = max_age_full
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
            output += gen_plots(archive_reports, log, archive_dir, full_archive)
            output += toggle_link
            output += '<table border="1" cellspacing="3" cellpadding="5">\n'
            output += "<tr><th>Commit</th><th>Status</th><th>Summary</th>"
            output += "<th>Author</th><th>Commit Message</th></tr>\n\n"
            output += "".join(all_html_rows[: max_age + 1])
            output += "</table>\n"
            output += toggle_link
            output += html_footer()
            write_file(archive_dir / html_out_postfix, output)

            url_list = "".join(all_url_rows[: max_age + 1])
            write_file(archive_dir / urls_out_postfix, url_list)


# ======================================================================================
def gen_url_list(config: ConfigParser, outdir: Path) -> None:
    print("Working on url lists.")
    for postfix in ("list_full.txt", "list_recent.txt"):
        url_list = ""
        for s in config.sections():
            fn = outdir / f"archive/{s}/{postfix}"
            if not path.exists(fn):
                continue
            with open(fn, encoding="utf8") as f:
                url_list += f.read()
        write_file(outdir / f"archive/{postfix}", url_list)


# ======================================================================================
def gen_plots(
    archive_reports: Dict[GitSha, Report], log: GitLog, outdir: Path, full_archive: bool
) -> str:

    ordered_shas = [c.sha for c in log.commits]
    ordered_reports = [archive_reports[s] for s in ordered_shas if s in archive_reports]

    # aggregate plot data
    aggregated_plots: Dict[str, AggregatedPlot] = {}
    for report in ordered_reports:
        for p in report.plots:
            if p.name not in aggregated_plots:
                aggregated_plots[p.name] = AggregatedPlot(p)
        for pp in report.plot_points:
            plot = aggregated_plots[pp.plot_name]
            cname = pp.curve_name
            if cname not in plot.curves_by_name:
                plot.curves_by_name[cname] = PlotCurve(label=pp.curve_label)
            assert report.sha
            age = log.index_by_sha[report.sha]
            plot.curves_by_name[cname].x.append(-age)
            plot.curves_by_name[cname].y.append(pp.y)
            plot.curves_by_name[cname].yerr.append(pp.yerr)

    max_age_full = max([-1] + [log.index_by_sha[sha] for sha in archive_reports])
    max_age = max_age_full if full_archive else 100

    output = ""
    for pname, plot in aggregated_plots.items():
        print(f"Working on plot: {pname}")
        fn = f"{pname}_full" if full_archive else pname
        make_plot_data(outdir / Path(f"{fn}.txt"), plot, max_age, log)
        make_plot_image(outdir / Path(f"{fn}.png"), plot, max_age, full_archive)
        output += f'<a href="{fn}.txt"><img src="{fn}.png" alt="{plot.title}"></a>\n'

    return output


# ======================================================================================
def make_plot_data(fn: Path, plot: AggregatedPlot, max_age: int, log: GitLog) -> None:
    output = f"# {plot.title}\n"
    headers = [f"{c.label:>6}" for c in plot.curves]
    output += "# " + "\t".join(["age", "commit"] + headers) + "\n"
    for age in range(max_age + 1):
        values = []
        for curve, header in zip(plot.curves, headers):
            if -age in curve.x:
                i = curve.x.index(-age)
                value = "{:.1f}".format(curve.y[i])
            else:
                value = "?"
            values.append(value.rjust(len(header)))
        if any([v.strip() != "?" for v in values]):
            index_colums = [f"{-age:5d}", log.commits[age].sha[:6]]
            output += "\t".join(index_colums + values) + "\n"

    write_file(fn, output)


# ======================================================================================
def make_plot_image(
    fn: Path, plot: AggregatedPlot, max_age: int, full_archive: bool
) -> None:

    # Setup figure.
    fig = plt.figure(figsize=(12, 4))
    fig.subplots_adjust(bottom=0.18, left=0.06, right=0.70)
    fig.suptitle(plot.title, fontsize=14, fontweight="bold", x=0.4)
    ax = fig.add_subplot(111)
    ax.set_xlabel("Commit Age")
    ax.set_ylabel(plot.ylabel)
    markers = itertools.cycle("os>^*")

    # Draw curves.
    for curve in plot.curves:
        if full_archive:
            ax.plot(curve.x, curve.y, label=curve.label, linewidth=2)  # less crowded
        else:
            style = dict(marker=next(markers), linewidth=2, markersize=6)
            ax.errorbar(curve.x, curve.y, yerr=curve.yerr, label=curve.label, **style)
    ax.set_xlim(-max_age - 1, 0)
    ax.xaxis.set_minor_locator(AutoMinorLocator())

    # Draw legend.
    legend_style = dict(numpoints=1, fancybox=True, shadow=True, borderaxespad=0.0)
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", **legend_style)

    # Determine y-range such that all curves are visible while ignoring outlayers.
    visibles = []
    for curve in plot.curves:
        ys = [y for x, y in zip(curve.x, curve.y) if x >= -max_age]  # visible y-values
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
    fig.savefig(fn)
    plt.close(fig)


# ======================================================================================
def send_notification(
    report: Report,
    last_ok: Optional[GitSha],
    log: GitLog,
    name: str,
    s: str,
    send_emails: bool,
) -> None:

    idx_end = log.index_by_sha[report.sha] if report.sha else 0
    if not last_ok:
        return  # we don't know when this started
    idx_last_ok = log.index_by_sha[last_ok]
    if idx_end == idx_last_ok:
        return  # probably a flapping tester
    emails_raw = set(c.author_email for c in log.commits[idx_end:idx_last_ok])
    emails = [e for e in emails_raw if "noreply" not in e]
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
    msg_txt += f"   report state:   {report.status}\n"
    msg_txt += f"   report summary: {report.summary}\n"
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
def html_header(title: str) -> str:
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
def html_linkbox() -> str:
    output = '<div class="sidebox">\n'
    output += "<h2>More...</h2>\n"
    output += '<a href="regtest_survey.html">Regtest Survey</a><br>\n'
    output += '<a href="https://www.cp2k.org/static/coverage/">Test Coverage</a><br>\n'
    output += '<a href="discontinued_tests.html">Discontinued Tests</a><br>\n'
    output += "</div>\n"
    return output


# ======================================================================================
def html_gitbox(log: GitLog) -> str:
    now = datetime.utcnow()
    output = '<div class="sidebox">\n'
    output += "<h2>Recent Commits</h2>\n"
    for commit in log.commits[0:10]:
        url = "https://github.com/cp2k/cp2k/commit/" + commit.sha
        autor_name = html.escape(commit.author_name)
        title = html.escape(commit.message)
        cm = commit.message if len(commit.message) < 28 else commit.message[:26] + "..."
        output += f'<p><a title="{title}" href="{url}">{html.escape(cm)}</a><br>\n'
        delta = now - commit.date
        age = delta.days * 24.0 + delta.seconds / 3600.0
        output += "<small>git:" + commit.sha[:7]
        output += f"<br>\n{autor_name} {age:.1f}h ago.</small></p>\n"
    output += "</div>\n"
    return output


# ======================================================================================
def html_footer() -> str:
    now = datetime.utcnow().replace(microsecond=0)
    output = f"<p><small>Page last updated: {now.isoformat()}</small></p>\n"
    output += "</body></html>"
    return output


# ======================================================================================
def write_file(fn: Path, content: str, gz: bool = False) -> None:
    d = fn.parent
    if not d.exists():
        os.makedirs(d)
        print(f"Created dir: {d}")
    if fn.exists():
        with open(fn, "rb") as f:
            old_bytes = gzip.decompress(f.read()) if gz else f.read()
        old_content = old_bytes.decode("utf-8", errors="replace")
        if old_content == content:
            print(f"File did not change: {fn}")
            return

    encoded_content = content.encode("utf-8")
    with open(fn, "wb") as f:
        f.write(gzip.compress(encoded_content) if gz else encoded_content)
    print(f"Wrote: {fn}")


# ======================================================================================
def status_cell(status: ReportStatus, report_url: str, up_to_date: bool = True) -> str:
    if status == "OK":
        bgcolor = "#00FF00" if (up_to_date) else "#8CE18C"
    elif status == "FAILED":
        bgcolor = "#FF0000" if (up_to_date) else "#E18C8C"
    else:
        bgcolor = "#d3d3d3"
    return f'<td bgcolor="{bgcolor}"><a href="{report_url}">{status}</a></td>'


# ======================================================================================
def commit_cell(sha: Optional[GitSha], log: GitLog) -> str:
    if sha is None:
        return "<td>N/A</td>"
    idx = log.index_by_sha[sha]
    commit = log.commits[idx]
    github_url = "https://github.com/cp2k/cp2k/commit/" + sha
    return f'<td align="left"><a href="{github_url}">{sha[:7]}</a> ({-idx})</td>'


# ======================================================================================
def retrieve_report(url: str) -> Optional[str]:
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
        traceback.print_exc()
        return None


# ======================================================================================
def store_report(report: Report, report_txt: str, section: str, outdir: Path) -> None:
    fn = outdir / f"archive/{section}/commit_{report.sha}.txt.gz"
    write_file(fn, report_txt, gz=True)


# ======================================================================================
def parse_report(report_txt: Optional[str], log: GitLog) -> Report:
    if report_txt is None:
        return Report(status="UNKNOWN", summary="Error while retrieving report.")
    try:
        m = re.search("(^|\n)CommitSHA: (\w{40})\n", report_txt)
        sha = GitSha(m.group(2)) if m else None
        summary = re.findall("(^|[\n\r])Summary: (.+)\n", report_txt)[-1][1]
        status = re.findall("(^|\n)Status: (.+)\n", report_txt)[-1][1]
        plot_matches = re.findall("(^|\n)Plot: (.+)(?=\n)", report_txt)
        plots = [cast(Plot, eval(f"Plot({m[1]})")) for m in plot_matches]
        point_matches = re.findall("(^|\n)PlotPoint: (.+)(?=\n)", report_txt)
        points = [cast(PlotPoint, eval(f"PlotPoint({m[1]})")) for m in point_matches]

        # Check that every plot has at least one PlotPoint
        for plot in plots:
            if not any([p.plot_name == plot.name for p in points]):
                return Report("FAILED", f"Plot {plot.name} has no PlotPoints.")

        # Check that CommitSHA belongs to the master branch.
        if sha not in log.index_by_sha:
            return Report("FAILED", "Unknown CommitSHA.")

        return Report(status, summary, sha=sha, plots=plots, plot_points=points)
    except:
        traceback.print_exc()
        return Report("UNKNOWN", "Error while parsing report.")


# ======================================================================================
main()

# EOF
