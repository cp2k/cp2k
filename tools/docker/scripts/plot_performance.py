#!/usr/bin/env python3

# author: Ole Schuett

import re
import sys
from collections import OrderedDict

# ======================================================================================
def main():
    if len(sys.argv) < 4 or (len(sys.argv) - 1) % 3 != 0:
        print(
            "Usage: plot_performance.py <title_1> <plot_1> <file_1> ... "
            "<title_N> <plot_N> <file_N>"
        )
        sys.exit(1)

    titles, plots, timings = [], [], []
    for i in range((len(sys.argv) - 1) // 3):
        titles.append(sys.argv[3 * i + 1])
        plots.append(sys.argv[3 * i + 2])
        timings.append(parse_timings(sys.argv[3 * i + 3]))

    routines = list(set([r for t in timings for r in list(t.keys())[:6]]))
    routines.remove("total")
    routines.sort(reverse=True, key=lambda r: timings[0].get(r, 0.0))

    for titel, plot, timing in zip(titles, plots, timings):
        timing["rest"] = timing["total"] - sum([timing.get(r, 0.0) for r in routines])

        full_title = f"Timings of {titel}"
        print(f'Plot: name="{plot}", title="{full_title}", ylabel="time [s]"')
        for r in ["rest"] + routines:
            t = timing.get(r, 0.0)
            print(f'PlotPoint: plot="{plot}", name="{r}", label="{r}", y={t}, yerr=0.0')
        print("")


# ======================================================================================
def parse_timings(out_fn):
    output = open(out_fn, encoding="utf8").read()

    pattern = r"\n( -+\n - +-\n - +T I M I N G +-\n([^\n]*\n){4}.*? -+)\n"
    match = re.search(pattern, output, re.DOTALL)
    report_lines = match.group(1).split("\n")[7:-1]
    print("\nFrom {}:\n{}\n".format(out_fn, match.group(0).strip()))

    # Extract average self time.
    timings = {}
    for line in report_lines:
        parts = line.split()
        timings[parts[0]] = float(parts[3])

    # Add total runtime.
    assert report_lines[0].split()[0] == "CP2K"
    timings["total"] = float(report_lines[0].split()[5])

    # Sort by time, longest first.
    return OrderedDict(sorted(timings.items(), reverse=True, key=lambda kv: kv[1]))


# ======================================================================================
main()

# EOF
