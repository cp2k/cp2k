#!/usr/bin/env python3

# author: Ole Schuett

import sys
import math

# ===============================================================================
def main():
    if len(sys.argv) != 2:
        print("Usage: timings.py <timings.txt>")
        sys.exit(1)

    filename = sys.argv[1]
    with open(filename) as fhandle:
        timings = sorted(float(line.split()[0]) for line in fhandle.readlines())

    print('Plot: name="timings", title="Timing Distribution", ylabel="time [s]"')
    tmpl = 'PlotPoint: name="{}th_percentile", plot="timings", label="{}th %ile", y={}, yerr=0.0'
    for p in (100, 99, 98, 95, 90, 80):
        v = percentile(timings, p / 100.0)
        print(tmpl.format(p, p, v))


# ===============================================================================
def percentile(values, percent):
    k = (len(values) - 1) * percent
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return values[int(k)]
    d0 = values[int(f)] * (c - k)
    d1 = values[int(c)] * (k - f)
    return d0 + d1


# ===============================================================================
main()
