# CP2K Benchmarks

This directory contains input files for CP2K's benchmarks.

For measurements from different machines, please refer to [CP2K benchmark suite](https://www.cp2k.org/performance), and for documentation on CP2K's input files, please refer to the [Input Reference Manual](https://manual.cp2k.org/). Python scripts for generating the scaling graphs are provided in [tools/benchmark_plots/](../tools/benchmark_plots).

**Note:** the benchmark names make common use of acronyms. For explanations, please refer to the [Glossary of Acronyms and Abbreviations](https://www.cp2k.org/acronyms).

## Introduction

The purpose of the CP2K benchmark suite is to provide performance which can be used to guide users towards the best configuration (e.g. machine, number of MPI processors, number of OpenMP threads) for a particular problem, and give a good estimation for the parallel performance of the code for different types of methods.

The systems used to obtain the benchmark results are described on the [systems page](https://www.cp2k.org/performance:systems).

## Benchmarks

See the `README.md` inside each benchmark sub-directory for descriptions of each benchmark along with performance numbers.

Benchmarks currently available:

- [Fayalite-FIST](Fayalite-FIST)
- [QS](QS)
- [QS_DM_LS](QS_DM_LS)
- [QS_HFX](QS_HFX)
- [QS_diag](QS_diag)
- [QS_mp2_rpa](QS_mp2_rpa)
- [QS_ot_ls](QS_ot_ls)
- [QS_pao_ml_tio2](QS_pao_ml_tio2)
- [QS_rubri](QS_rubri)
- [QS_single_node](QS_single_node)
- [QS_stmv](QS_stmv)

### Run Benchmarks

Some benchmarks require a preliminary step to generate an input file, e.g. a wavefunction. When that is the case, it is specified in the `README.md` inside the benchmark's sub-directory.

The general way to run the benchmarks with the hybrid parallel executable is, e.g. for 2 threads per rank:

```
export OMP_NUM_THREADS=2
parallel_launcher launcher_options path_to_cp2k.psmp -i inputfile.inp -o logfile.log
```

where:

- The parallel_launcher is mpirun, mpiexec, or some variant such as aprun on Cray systems or srun when using Slurm.
- `launcher_options` specifies parallel placement in terms of total numbers of nodes, MPI ranks/tasks, tasks per node, and OpenMP threads per task (which should be equal to the value given to OMP_NUM_THREADS). This is not necessary if parallel runtime options are picked up by the launcher from the job environment.

The reported walltime for a given run can be obtained by querying the resulting `.log` file for CP2K's internal timing, as follows:
```
$ grep "CP2K     "  *.log
```

### Plotting

Python scripts for generating the scaling graphs are provided in [cp2k/tools/benchmark_plots/](../tools/benchmark_plots/).

## Contributing

We encourage you to contribute benchmark results from your own local cluster or HPC system - just run the inputs and add timings in the relevant sections below.
Please also update the [list of machines](https://www.cp2k.org/performance:systems) for which benchmark data is provided.

