#!/usr/bin/env python3

# author: Ole Schuett


import argparse
from collections import Counter


def summarize(issue_files, suppressions):
    suppress = []

    if suppressions:
        with open(suppressions, encoding="utf8") as fhandle:
            suppress = (line.rstrip() for line in fhandle)
            # ignore empty and commented out lines
            suppress = [line for line in suppress if line and not line.startswith("#")]

    issues = []

    for fn in issue_files:
        with open(fn, encoding="utf8") as fhandle:
            lines = (line.strip() for line in fhandle)
            # only add non-empty lines
            issues += [line for line in lines if line]

    issues = sorted(set(issues))
    issues_shown = [i for i in issues if i not in suppress]
    issues_supp = [i for i in issues if i in suppress]
    unused_supp = [i for i in suppress if i not in issues]
    n_shown = len(issues_shown)
    n_supp = len(issues_supp)
    n_unused = len(unused_supp)
    counts_by_conv: Counter[str] = Counter()
    for i in issues:
        conv = "other"
        if "https://cp2k.org/conv#" in i:
            conv = i.split("https://cp2k.org/conv#")[1]
        counts_by_conv[conv] += 1

    if n_supp:
        print("\nRegistered suppression entries:\n")
        for i in issues_supp:
            print(f"  {i}")

    if n_unused:
        print(f"\nThere are {n_unused} registered but unused suppressions:\n")
        for i in unused_supp:
            print(f"  {i}")

    if n_shown:
        print(f"\nFound {n_shown} unregistered unsuppressed issues:\n")
        for i in issues_shown:
            print(f"  {i}")

    # Plotting was formerly only the total number of suppressed issues;
    # now the number of issues per convention is also plotted, but the
    # plot title and the label for n_supp remain the same for backward
    # compatibility (which has a long history in the dashboard).
    print(
        '\nPlot: name="supps", title="Suppressed Convention Violations", ylabel="# suppressions"'
    )
    print(
        f'PlotPoint: name="coding", plot="supps", label="Coding Conventions", y={n_supp}, yerr=0'
    )
    for conv, count in counts_by_conv.most_common():
        print(
            f'PlotPoint: name="{conv}", plot="supps", label="{conv}", y={count}, yerr=0'
        )

    # Print summary and status
    print(f"Summary: Found {n_shown} issues ({n_supp} suppressed)")
    print("Status: {status}".format(status="FAILED" if issues_shown else "OK"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Combines multiple files with issues into a dashboard-report"
    )
    parser.add_argument(
        "files",
        metavar="<issue-file>",
        type=str,
        nargs="+",
        help="files containing the logs of detected issues",
    )
    parser.add_argument("--suppressions", type=str)
    args = parser.parse_args()

    summarize(args.files, args.suppressions)
