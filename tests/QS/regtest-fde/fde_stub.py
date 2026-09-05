#!/usr/bin/env python3
import math

GRID_FILE, RHOA_FILE, RES_FILE = "EMBGRID", "EMBRHOA", "EMBRES"

R0 = (3.77945, 3.77945, 4.77945)
ALPHA = 0.6
NELEC = 2.0
Z = 2.0
E_CASSCF = -2.85
E_EMB = 0.0
E_NUC = 0.0

with open(GRID_FILE) as fh:
    npoints = int(fh.readline().split()[0])
    coords = [tuple(float(v) for v in fh.readline().split()[:3]) for _ in range(npoints)]

with open(RHOA_FILE, "w") as fh:
    fh.write("%d\n" % npoints)
    for (x, y, z) in coords:
        d2 = (x - R0[0])**2 + (y - R0[1])**2 + (z - R0[2])**2
        fh.write("%.10E\n" % math.exp(-ALPHA * d2))

with open(RES_FILE, "w") as fh:
    fh.write("%.10f\n%.10f\n%.10f\n%.10f\n1\n" % (E_CASSCF, E_EMB, E_NUC, NELEC))
    fh.write("%.6f %.8f %.8f %.8f\n" % (Z, R0[0], R0[1], R0[2]))
