#!/usr/bin/env python3

# author: Ole Schuett


import argparse


def summarize(issue_files, suppressions):

    suppress = []

    if suppressions:
        with open(suppressions) as fhandle:
            suppress = (line.rstrip() for line in fhandle)
            # ignore empty and commented out lines
            suppress = [line for line in suppress if line and not line.startswith("#")]

    issues = []

    for fn in issue_files:
        with open(fn) as fhandle:
            lines = (line.strip() for line in fhandle)
            # only add non-empty lines
            issues += [line for line in lines if line]

    issues = sorted(set(issues))
    issues_shown = [i for i in issues if i not in suppress]
    issues_supp = [i for i in issues if i in suppress]
    unused_supp = [i for i in suppress if i not in issues]

    print("Suppressed entries:\n")
    for i in issues_supp:
        print("  {}".format(i))
    print()

    print("There are %d unused suppressions:\n" % len(unused_supp))
    for i in unused_supp:
        print("  {}".format(i))
    print()

    print("Found {} unsuppressed issues:\n".format(len(issues_shown)))
    for i in issues_shown:
        print("  {}".format(i))

    print(
        """
Plot: name="supps", title="Suppressed Convention Violations", ylabel="# suppressions"
PlotPoint: name="coding", plot="supps", label="Coding Conventions", y={suppressed}, yerr=0
Summary: Found {issues} issues ({suppressed} suppressed)
Status: {status}""".format(
            issues=len(issues_shown),
            suppressed=len(issues_supp),
            status="FAILED" if issues_shown else "OK",
        )
    )


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
