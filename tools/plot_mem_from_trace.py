#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import re
import matplotlib.pyplot as plt

# ===============================================================================
def main():
    if len(sys.argv) != 2:
        print("Usage: plot_mem_from_trace.py <trace.out>")
        sys.exit(1)

    re_mem = re.compile("Hostmem: (\d+) MB GPUmem: (\d+) MB")
    trace_fn = sys.argv[1]
    f = open(trace_fn)
    hostmem_stats = []
    gpumem_stats = []
    for line in f.readlines():
        m = re_mem.search(line)
        if not m:
            continue
        hostmem_stats.append(m.group(1))
        gpumem_stats.append(m.group(2))

    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    ax1.plot(hostmem_stats, color="red", label="Host")
    ax1.set_ylabel("Host Memory Usage [MB]")
    ax1.set_xlabel("time")

    ax2 = ax1.twinx()
    ax2.plot(gpumem_stats, color="green", label="GPU")
    ax2.set_ylabel("GPU Memory Usage [MB]")

    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    plt.legend(handles1 + handles2, labels1 + labels2)
    plt.show()


# ===============================================================================
main()
# EOF
