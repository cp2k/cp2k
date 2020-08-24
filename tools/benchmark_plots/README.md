# Plotting benchmark results

Two python scripts are provided to help plot benchmark results:

- `plot_benchmark.py`  - Plots results for an individual benchmark
- `plot_comparison.py` - Plots a comparison of results for an individual benchmark

Library requirements:

- `numpy`
- `matplotlib`

## `plot_benchmark.py` usage

```shell
> python plot_benchmark.py --name <benchmarkname> --machine <machinename>  \
  [--shift <YSHIFT>] [--show]
```

## `plot_benchmark.py` arguments

- `--name`: Name of the benchmark to plot results for, e.g. fayalite-FIST
- `--machine`: Machine to plot results for, e.g. Magnus
- `--shift`: Fractional shift of label in y-direction, for fine tuning, e.g. 1.1
- `--show`: Show the plot window as well as saving to a file

## `plot_benchmark.py` assumptions

The results will be read from files called:
`<machine>_benchmarks/<name>/<name>_besttimes.txt`
(converted to lower case) in the directory defined as `BASE`.
This `BASE` directory should be edited directly in the script to set the
correct location for the benchmark results. See below for the format
of this file.

## `plot_comparison.py` usage

```shell
> python plot_comparison.py --name <benchmarkname> --machines <machinenames>  \
  [--shifts <YSHIFT>] [--show]
```

## `plot_comparison.py` arguments

- `--name`: Name of the benchmark to plot results for, e.g. fayalite-FIST
- `--machine`: List of machines to plot results for, e.g. Magnus ARCHER HECToR
- `--shift`: Fractional shifts of label in y-direction, for fine tuning,
  e.g. `1.1 1.1 0.9`.
  If less shifts are specified than machines then they will be used in rotation.
- `--show`: Show the plot window as well as saving to a file

## `plot_comparison.py` assumptions

The results will be read from files called:
`<machine>_benchmarks/<name>/<name>_besttimes.txt`
(converted to lower case) in the directory defined as `BASE`.
This `BASE` directory should be edited directly in the script to set the
correct location for the benchmark results. See below for the format
of this file.

## Result File Format

The result files should contain the result for a single run on each line in
the form:

```output
<nodes> <time> <configuration>
```

Where:

- `<nodes>`: is the number of nodes run on
- `<time>`: is the time taken in seconds
- `<configuration>`: is the configuration for the best time, eg `2_TH` if the
  best result is obtained for 2 OpenMP threads per MPI task.
