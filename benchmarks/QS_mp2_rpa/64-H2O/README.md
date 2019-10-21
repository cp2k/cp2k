# Quickstep Moeller-Plesset perturbation theory to 2nd order Random-Phase Approximation - 64 H2O (Regression Tests)

Hybrid benchmark for RI-MP2 and RI-dRPA (64-H2O-TZ).

## Description of Input Files

- `H2O-64-PBE-TZ.inp`: needed to generate an initial wfn for the SCF runs
- `H2O-64-RI-MP2-TZ.inp`: actual RI-MP2 benchmark
- `H2O-64-RI-dRPA-TZ.inp`: actual RI-dRPA benchmark

## Additional files

- `BASIS_H2O`: contains the primary and auxiliary(RI) basis sets
- `POTENTIAL_H2O`: contains the GTH pseudo potentials
- `H2O-64.xyz`: geometry in xyz format

the additional files [`t_c_g.dat`](cp2k/data/t_c_g.dat) is needed for the RI-MP2 run, and can be found in the cp2k/data directory.

## Benchmark Requirements

To run these benchmarks, CP2K needs to be compiled with libint support (-D__LIBINT). Libint library has to be compiled such that higher angular momentum can be computed (see: tools/hfx_tools/libint_tools/README_LIBINT), use, for example, --with-libint-max-am=6 --with-libderiv-max-am1=5.

It is advantageouss to have a OMP/MPI hybrid code (cp2k.psmp).
In particular the RI-MP2 and RI-dRPA inputs are suitable for being used with 2 threads per task.

## How to Run the Benchmark

1) run H2O-64-PBE-TZ.inp, this will generate the file H2O-64-PBE-TZ-RESTART.wfn, necessary for the other two runs.
2) run H2O-64-RI-MP2-TZ.inp for RI-MP2
3) run H2O-64-RI-dRPA-TZ.inp for RI-dRPA

## Results Archive

| Input File            | Configuration             | Total Number of Cores| Runtime |
| --------------------- | -------------------------:| --------------------:| -------:|
| H2O-64-PBE-TZ.inp     |           128 MPI x 2 OMP |                  256 | 50 s    |
| H2O-64-RI-MP2-TZ.inp  |          1024 MPI x 2 OMP |                 2048 | ~18 min |
| H2O-64-RI-dRPA-TZ.inp |          1024 MPI x 2 OMP |                 2048 | ~12 min |

*) The timings have been obtained on CRAY-XE6 (MonteRosa@CSCS)
