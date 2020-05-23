#!/usr/bin/env python3


"""
Make a plot of CP2K benchmark data.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker
import argparse

# Base directory
BASE = "/home/jeremy/cp2k-uk/results/"
YSHIFT = 1.1


def makePlot(benchmarkName, machineName, shift, showPlot):
    """
    Make a plot of benchmark results.

    Arguments:
    benchmarkName -- the name of the benchmark to plot results for
    machineName -- the name of the machine to get results for
    shift -- the shift of the data point label on the y-axis
    showPlot -- whether to show the plot in a window when done
    """

    # Turn interactive mode off (required in ipython)
    plt.ioff()

    # Construct file name based on machine name and benchmark name
    fileName = (
        BASE
        + machineName.lower()
        + "_benchmarks/"
        + benchmarkName
        + "/"
        + benchmarkName
        + "_besttimes.txt"
    )

    print("Processing benchmark results for %s from file:" % machineName)
    print(fileName)

    # Load data from file
    data = np.loadtxt(
        fileName,
        dtype={"names": ("nodes", "time", "label"), "formats": ("i4", "f4", "a5")},
    )

    # Create plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # Plot data
    ax.loglog(data["nodes"], data["time"], marker="o")
    # Change axes ticks to integers
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.FormatStrFormatter("%d"))
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.FormatStrFormatter("%d"))
    # Add title and axes labels
    title = "Performance of the " + benchmarkName + " benchmark on " + machineName
    plt.title(title)
    plt.xlabel("Number of nodes used")
    plt.ylabel("Time (seconds)")
    # Label the points with the best performance
    for i, labelText in enumerate(data["label"]):
        ax.annotate(
            labelText,
            xy=(data["nodes"][i], data["time"][i]),
            xytext=(data["nodes"][i], shift * data["time"][i]),
        )

    # Save plot to file
    outFileName = benchmarkName + "-" + machineName.lower() + ".png"

    plt.savefig(outFileName)
    if showPlot:
        plt.show()
    plt.clf()


if __name__ == "__main__":
    # Command line options
    parser = argparse.ArgumentParser(description="Plot benchmark results.")
    parser.add_argument(
        "--name", default="fayalite-FIST", help="name of the benchmark to process"
    )
    parser.add_argument(
        "--machine", default="Magnus", help="name of the machine used to obtain results"
    )
    parser.add_argument(
        "--shift",
        default=YSHIFT,
        help="Fractional shift of label in y-direction",
        type=float,
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="show the plot window as well as saving to file",
    )

    args = parser.parse_args()

    benchmarkName = args.name
    machineName = args.machine
    shift = args.shift
    showPlot = args.show

    makePlot(benchmarkName, machineName, shift, showPlot)
