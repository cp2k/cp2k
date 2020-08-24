# Quickstep Moeller-Plesset perturbation theory to 2nd order Random-Phase Approximation

Hybrid benchmark for RI-MP2 and RI-dRPA.

## Description

This benchmark is a single-point energy calculation using 2nd order Moeller-Plesset
perturbation theory (MP2) with the Resolution-of-the-Identity approximation to
calculate the exchange-correlation energy.

## Benchmark Requirements

To run these benchmarks, CP2K needs to be compiled with libint support (`-D__LIBINT`).
Libint library has to be compiled such that higher angular momentum can be computed
(see: [libint-cp2k/README](https://github.com/cp2k/libint-cp2k), use,
for example, `--with-libint-max-am=6 --with-libderiv-max-am1=5`.

It is advantageouss to have a OMP/MPI hybrid code (cp2k.psmp). In particular the
RI-MP2 and RI-dRPA inputs are suitable for being used with 2 threads per task.
