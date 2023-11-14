# Quickstep low-scaling RPA and SOS-MP2 gradients benchmark

## Description

This benchmark set showcases sparse tensor based low-scaling post-HF ENERGY_FORCE calculations. Both
the RPA and SOS-MP2 methods are represented. The implementation relies on the RI approximation in
the AO basis. For increased sparsity and performance, the overlap RI metric is used. All input files
are of production level. The methods make use of the DBT and DBM libraries, which are GPU
accelerated.

## Files

There are three folders for the simulation of periodic water systems of increasing size, from 32 to
128 molecules. Both the RPA and the SOS-MP2 calculations restart from a PBE wavefunction, which has
to be calculated first (H2O-#-PBE-TZ.inp). All calculations are based on the triple-zeta cc-TZ and
RI_TZ basis sets. The SOS-MP2 correlation energy is calculated on top of a HF SCF, and the RPA
correlation energy is calculated on top of a PBE SCF (without EXX).

## Reference timings

These benchmarks were run during the pilot phase of the
[LUMI-G supercomputer](https://docs.lumi-supercomputer.eu/hardware/compute/lumig/) (December 2022),
under
[git commit d8f3624](https://github.com/cp2k/cp2k/commit/d8f36242127b9a828f127550ecd9613eefb3f1cc).
All calculations were run with 4 GPUs (8 GCDs), 16 MPI ranks, and 3 OMP threads per node. Note that
this amounts to a total of 48 out of the 64 CPUs of a LUMI-G node, because full nodes where not
available at the time.

### 32-H2O

| Input file            | Wall time (s) on 2 nodes | Wall time (s) on 4 nodes |
| --------------------- | ------------------------ | ------------------------ |
| H2O-32-RPA-TZ.inp     | 1531                     | 788                      |
| H2O-32-SOS-MP2-TZ.inp | 1787                     | 918                      |

### 64-H2O

| Input file            | Wall time (s) on 8 nodes | Wall time (s) on 16 nodes |
| --------------------- | ------------------------ | ------------------------- |
| H2O-64-RPA-TZ.inp     | 2570                     | 1415                      |
| H2O-64-SOS-MP2-TZ.inp | 2785                     | 1505                      |

### 128-H2O

| Input file             | Wall time (s) on 32 nodes | Wall time (s) on 64 nodes |
| ---------------------- | ------------------------- | ------------------------- |
| H2O-128-RPA-TZ.inp     | 3710                      | 2578                      |
| H2O-128-SOS-MP2-TZ.inp | 3702                      | 2539                      |

## Notes on parameter choice and performance

Note that the choice of the overlap RI metric is the most performant because it leads to the highest
sparsity. Alternatively, the truncated Coulomb operator could be used, for example with a cutoff
radius of 1-2 Angstroms, which could lead to more accurate results, especially for system with
smaller basis sets (e.g. double-zeta). Using the exact same input files, but with a truncated
Coulomb metric of range 1 Ang, leads to calculations 20% to 30% more expensive.

Most input files use 4 as the value for MEMORY_CUT. This controls the batching procedure of sparse
tensor contractions such that the memory footprint remains limited. Note that using a high value for
this keyword leads to overheads. As an illustration, the H2O-128-RPA-TZ.inp input file with
MEMORY_CUT set to 3 does not run on 32 LUMI-G nodes because of memory. However, it runs on 64 nodes
with a wall time of 2149 seconds (compared to 2578 seconds with MEMORY_CUT set to 4).

Finally, the value of 6 was chosen for the MIN_BLOCK_SIZE keyword, while the default is 5. This is a
trade-off between better sparsity for smaller block size, and higher efficiency of block-wise
matrix-matrix mulitplications for larger sizes. This mostly matters for the force calculations of
large system, because the bottleneck of the calculation involves the contraction of fairly dense
tensors. The H2O-128-RPA-TZ.inp input file with MIN_BLOCK_SIZE 5 runs in 4105 seconds on 32 nodes,
versus 3710 seconds with a block size of 6.
