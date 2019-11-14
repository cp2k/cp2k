# Quickstep Moeller-Plesset perturbation theory to 2nd order Random-Phase Approximation - 128 H2O

Hybrid benchmark for RI-MP2 and RI-dRPA.

## Description

This benchmark is a single-point energy calculation using 2nd order Møller-Plesset perturbation theory (MP2) with the Resolution-of-the-Identity approximation to calculate the exchange-correlation energy.

## Description of Input Files

- [`H2O-128-PBE-TZ.inp`](H2O-128-PBE-TZ.inp): needed to generate an initial wfn for the SCF runs
- [`H2O-128-RI-dRPA-TZ.inp`](H2O-128-RI-dRPA-TZ.inp): actual RI-dRPA benchmark

## Additional files

- [`BASIS_H2O`](BASIS_H2O): contains the primary and auxiliary(RI) basis sets
- [`POTENTIAL_H2O`](POTENTIAL_H2O): contains the GTH pseudo potentials
- [`H2O-128.xyz`](H2O-128.xyz): geometry in xyz format

## Benchmark Requirements

To run these benchmarks, CP2K needs to be compiled with libint support (-D__LIBINT). Libint library has to be compiled such that higher angular momentum can be computed (see: tools/hfx_tools/libint_tools/README_LIBINT), use, for example, --with-libint-max-am=6 --with-libderiv-max-am1=5.

It is advantageouss to have a OMP/MPI hybrid code (cp2k.psmp).
In particular the RI-MP2 and RI-dRPA inputs are suitable for being used with 2 threads per task.

## How to Run the Benchmark

1) run `H2O-128-PBE-TZ.inp`, this will generate the file `H2O-128-PBE-TZ-RESTART.wfn`, necessary for the other two runs.
3) run `H2O-128-RI-dRPA-TZ.inp` for RI-dRPA

## Results

### Results on Piz Daint, CSCS

| Input File             | Date       | CP2K Git SHA | Number of nodes | Node Configuration  | Runtime |
| ---------------------- | ---------- | ------------:| ---------------:| ------------------- | ------- |
| H2O-128-PBE-TZ.inp     | 2019-08-19 | 4519a8ad7    | 4 nodes         | 12 MPI x 1 OMP      | ~2 min  |
| H2O-128-RI-dRPA-TZ.inp | 2019-08-19 | 4519a8ad7    | 128 nodes       | 2 MPI x 6 OMP       | 80 min  |

*) The timings have been obtained on CRAY-XC50 (PizDaint@CSCS, GPU partition)

