# Quickstep Moeller-Plesset perturbation theory to 2nd order Random-Phase Approximation - 32 H2O (Regression Tests)

Hybrid benchmark for RI-MP2 and RI-dRPA (32-H2O-TZ).

## Description of Input Files

- `RI-MP2.inp`: farming input for measuring MP2 time
- `RI-RPA.inp`: farming input for measuring RPA time

### Additional files

- `BASIS_H2O`: contains the primary and auxiliary(RI) basis sets
- `H2O-32.xyz`: geometry in xyz format
- `H2O-32-PBE-TZ.inp`: needed to generate an initial DFT wfn (RPA, MP2)
- `H2O-32-HF-TZ.inp`: needed to refine DFT wfn at HF level  (MP2)
- `H2O-32-RI-MP2-TZ.inp`: actual RI-MP2 benchmark               (MP2)
- `H2O-32-RI-RPA-TZ.inp`: actual RI-RPA benchmark               (RPA)

the additional files [`t_c_g.dat`](cp2k/data/t_c_g.dat) and [`POTENTIAL`](cp2k/data/POTENTIAL) are taken from cp2k/data directory.

## Benchmark Requirements

To run these benchmarks, CP2K needs to be compiled with libint support (-D__LIBINT). Libint library has to be compiled such that higher angular momentum can be computed (see: tools/hfx_tools/libint_tools/README_LIBINT), use, for example, --with-libint-max-am=6 --with-libderiv-max-am1=5.

## How to Run the Benchmark

1) run `RI-MP2.inp`
2) run `RI-RPA.inp`

## Results Archive

| Input File | Configuration             | Total Number of Cores| Runtime [s]  |
| ---------- | -------------------------:| --------------------:| ------------:|
| RI-MP2.inp | 16 nodes x 16 MPI x 1 OMP |                  256 |          392 |
| RI-RPA.inp | 16 nodes x 16 MPI x 1 OMP |                  256 |          221 |

*) The timings have been obtained on CRAY-XC40 (PizDora@CSCS)
