#!/usr/bin/env python3

# author: Ole Schuett

import re
import sys
from collections import OrderedDict

# ======================================================================================
def main():
    if len(sys.argv) != 6:
        print(
            "Usage: plot_performance.py <nprocs> <title> <label> <output_omp> <output_mpi>"
        )
        sys.exit(1)
    nprocs, title, label, output_omp, output_mpi = sys.argv[1:]

    timings_omp = parse_timings(output_omp)
    timings_mpi = parse_timings(output_mpi)
    routines = list(set(list(timings_omp.keys())[:6] + list(timings_mpi.keys())[:6]))
    routines.remove("total")
    routines.sort(reverse=True, key=lambda r: timings_omp.get(r, 0.0))

    sum_routines_omp = sum([timings_omp.get(r, 0.0) for r in routines])
    sum_routines_mpi = sum([timings_mpi.get(r, 0.0) for r in routines])
    timings_omp["rest"] = timings_omp["total"] - sum_routines_omp
    timings_mpi["rest"] = timings_mpi["total"] - sum_routines_mpi

    full_title = f"Timings of {title} with {nprocs} OpenMP Threads"
    plot = f"{label}_timings_{nprocs}omp"
    print(f'Plot: name="{plot}", title="{full_title}", ylabel="time [s]"')
    for r in ["rest"] + routines:
        t = timings_omp[r]
        print(f'PlotPoint: plot="{plot}", name="{r}", label="{r}", y={t}, yerr=0.0')
    print("")

    full_title = f"Timings of {title} with {nprocs} MPI Ranks"
    plot = f"{label}_timings_{nprocs}mpi"
    print(f'Plot: name="{plot}", title="{full_title}", ylabel="time [s]"')
    for r in ["rest"] + routines:
        t = timings_mpi[r]
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
