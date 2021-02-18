#!/usr/bin/env python3

# author: Ole Schuett

import re
import sys
from collections import OrderedDict

# ======================================================================================
def main():
    if len(sys.argv) != 5:
        print("Usage: plot_performance.py <title> <label> <output_omp> <output_mpi>")
        sys.exit(1)
    title, label, output_omp, output_mpi = sys.argv[1:]

    timings_omp = parse_timings(output_omp)
    timings_mpi = parse_timings(output_mpi)
    routines = list(timings_omp.keys())[:6]

    full_title = f"Timings of {title} with OpenMP"
    plot = f"{label}_timings"
    print(f'Plot: name="{plot}", title="{full_title}", ylabel="time [s]"')
    for r in routines:
        t = timings_omp[r]
        print(f'PlotPoint: plot="{plot}", name="{r}", label="{r}", y={t}, yerr=0.0')
    print("")

    full_title = f"Parallel Efficiency of {title}"
    plot = f"{label}_scaling"
    print(f'Plot: name="{plot}", title="{full_title}", ylabel="MPI / OMP time"')
    for r in routines:
        e = timings_mpi[r] / timings_omp[r]
        print(f'PlotPoint: plot="{plot}", name="{r}", label="{r}", y={e}, yerr=0.0')


# ======================================================================================
def parse_timings(out_fn):
    output = open(out_fn, encoding="utf8").read()

    pattern = r"\n( -+\n - +-\n - +T I M I N G +-\n([^\n]*\n){4}.*? -+)\n"
    report_lines = re.search(pattern, output, re.DOTALL).group(1).split("\n")[7:-1]

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
