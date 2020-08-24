# Low-scaling GW performance test for a graphene nanoribbon

## Description

This benchmark is a low-scaling G0W0 calculation of a graphene nanoribbon.

## Benchmark Requirements

To run the benchmark, CP2K needs to be compiled with libint support (`-D__LIBINT`).

## Description of Input Files

- [`GW.inp`](GW.inp): input for measuring the low-scaling GW time

## Additional files

- [`nanoribbon.xyz`](nanoribbon.xyz): geometry in xyz format

## Results on Piz Dora, CSCS

| Input File | Configuration             | Total Number of Cores| Runtime [s]  |
| ---------- | -------------------------:| --------------------:| ------------:|
| GW.inp     | 16 nodes x 36 MPI x 1 OMP |                  576 |          305 |

The timings have been obtained on CRAY-XC40 (PizDora@CSCS)
