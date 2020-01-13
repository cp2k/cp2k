# Quickstep Moeller-Plesset perturbation theory to 2nd order Random-Phase Approximation - 32 H2O

## Description of Input Files

- [`RI-MP2.inp`](RI-MP2.inp): farming input for measuring MP2 time
- [`RI-RPA.inp`](RI-RPA.inp): farming input for measuring RPA time

### Additional files

- [`BASIS_H2O`](BASIS_H2O): contains the primary and auxiliary(RI) basis sets
- [`H2O-32.xyz`](H2O-32.xyz): geometry in xyz format
- [`H2O-32-PBE-TZ.inp`](H2O-32-PBE-TZ.inp): needed to generate an initial DFT wfn (RPA, MP2)
- [`H2O-32-HF-TZ.inp`](H2O-32-HF-TZ.inp): needed to refine DFT wfn at HF level (MP2)
- [`H2O-32-RI-MP2-TZ.inp`](H2O-32-RI-MP2-TZ.inp): actual RI-MP2 benchmark (MP2)
- [`H2O-32-RI-dRPA-TZ.inp`](H2O-32-RI-dRPA-TZ.inp): actual RI-RPA benchmark (RPA)

the additional files [`t_c_g.dat`](../../../data/t_c_g.dat) and [`POTENTIAL`](../../../data/POTENTIAL) are taken from [cp2k/data](../../../data) directory.

## How to Run the Benchmark

1) run `RI-MP2.inp`
2) run `RI-RPA.inp`

## Results

### Results on Piz Dora, CSCS

| Input File | Configuration             | Total Number of Cores| Runtime [s]  |
| ---------- | -------------------------:| --------------------:| ------------:|
| RI-MP2.inp | 16 nodes x 16 MPI x 1 OMP |                  256 |          392 |
| RI-RPA.inp | 16 nodes x 16 MPI x 1 OMP |                  256 |          221 |

*) The timings have been obtained on CRAY-XC40 (PizDora@CSCS)
