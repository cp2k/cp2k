#!/usr/bin/env python3


"""
Make a plot of CP2K benchmark data.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker
import argparse
from itertools import cycle

# Base directory
BASE = "/home/jeremy/cp2k-uk/results/"
YSHIFT = 1.1


def makePlot(benchmarkName, machineNames, shifts, showPlot):
    """
    Make a comparison plot of benchmark results.

    Arguments:
    benchmarkName -- the name of the benchmark to plot results for
    machineNames -- the names of the machines to get results for
    shifts -- the shifts of the data point label on the y-axis
    showPlot -- whether to show the plot in a window when done
    """
    # Turn interactive mode off (required in ipython)
    plt.ioff()

    markers = ["o", "v", "^", "s", "*"]
    markerCycler = cycle(markers)
    shiftCycler = cycle(shifts)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    title = "Performance of the " + benchmarkName + " benchmark"
    plt.title(title)
    plt.xlabel("Number of nodes used")
    plt.ylabel("Time (seconds)")

    for machineName in machineNames:
        fileName = (
            BASE
            + machineName.lower()
            + "_benchmarks/"
            + benchmarkName
            + "/"
            + benchmarkName
            + "_besttimes.txt"
        )
        print("Processing benchmark results from file:")
        print(fileName)
        data = np.loadtxt(
            fileName,
            dtype={"names": ("nodes", "time", "label"), "formats": ("i4", "f4", "a5")},
        )
        ax.loglog(
            data["nodes"], data["time"], marker=next(markerCycler), label=machineName
        )

        shift = next(shiftCycler)
        alignment = "left"
        if shift < 1.0:
            alignment = "right"
        for i, labelText in enumerate(data["label"]):
            ax.annotate(
                labelText,
                xy=(data["nodes"][i], data["time"][i]),
                xytext=(data["nodes"][i], shift * data["time"][i]),
                horizontalalignment=alignment,
            )

    ax.get_xaxis().set_major_formatter(matplotlib.ticker.FormatStrFormatter("%d"))
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.FormatStrFormatter("%d"))
    plt.legend()

    outFileName = benchmarkName + "-comparison"
    for machineName in machineNames:
        outFileName = outFileName + "-" + machineName.lower()
    outFileName = outFileName + ".png"

    plt.savefig(outFileName)

    if showPlot:
        plt.show()
    plt.clf()


if __name__ == "__main__":
    # Command line arguments
    parser = argparse.ArgumentParser(description="Plot benchmark results.")
    parser.add_argument(
        "--name", default="fayalite-FIST", help="name of the benchmark to process"
    )
    parser.add_argument(
        "--machines",
        default="Magnus",
        help="name of the machines used to obtain results",
        nargs="*",
    )
    parser.add_argument(
        "--shifts",
        default=[YSHIFT],
        help="Fractional shift of label in y-direction",
        nargs="*",
        type=float,
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="show the plot window as well as saving to file",
    )

    args = parser.parse_args()

    benchmarkName = args.name
    machineNames = args.machines
    shifts = args.shifts
    showPlot = args.show

    makePlot(benchmarkName, machineNames, shifts, showPlot)
