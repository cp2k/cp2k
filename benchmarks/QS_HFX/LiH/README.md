# Quickstep Hartree Fock eXchange - LiH (Regression Tests)

Hybrid benchmark to test CP2K scaling up to 10000s of cores

## Files description

- `input_bulk_HFX_3.inp`: is the actual benchmark
- `input_bulk_B88_3.inp`: needed to generate an initial wfn (wave function) for the HFX runs (not part of the benchmark)
- the additional files [`t_c_g.dat`](../../../data/t_c_g.dat) and [`POTENTIAL`](../../..//data/POTENTIAL) are needed, and can be found in the `cp2k/data/` directory.

## Benchmark Requirements

To run these this benchmark, CP2K needs to be compiled with libint support (-D__LIBINT), and it is advantageous to have a OMP/MPI hybrid code (cp2k.psmp). 

## How to Run the Benchmark

1) run `input_bulk_B88_3.inp` and copy the last `LiH_bulk_3-RESTART.wfn` to `B88.wfn` (5min. with 256x1 mpixomp tasks)
2) run `input_bulk_HFX_3.inp` on a number of nodes (about 30min on 2048 cores)

The obtained energies should be similar to (obtained on 2048 MPI x 8 OMP cores):

```
     1  OT CG          0.15E+00    139.386     0.0014482195     -870.3838179520
     2  OT LS          0.30E+00     19.767     1.0000000000     -870.7951031236
     3  OT CG          0.30E+00     20.724     0.0001614074     -870.9275954835
     4  OT LS          0.31E+00     19.665     1.0000000000     -870.9346994255
     5  OT CG          0.31E+00     20.641     0.0000168390     -870.9347195049
     6  OT LS          0.29E+00     19.760     1.0000000000     -870.9347906736
     7  OT CG          0.29E+00     21.192     0.0000009576     -870.9347911831
```

## Results Archive

### 2009-04-15

Running on ORNL's Cray XT5 (Jaguar) the following runtime has been obtained in a setup using 8 OMP threads per node (8 cores per node / 16 Gb per node).

| Cores | Full CP2K[s] | HFX[s]  | local HFX[s] | Mem/node[Mb] |
| -----:| ------------:| -------:| ------------:| ------------:|
|  1024 |      6569.07 | 6488.24 |      6406.99 |     14845.00 |
|  2048 |      1369.92 | 1314.86 |      1259.82 |     13468.00 |
|  4096 |       722.24 |  670.18 |       631.11 |      7293.00 |
|  8192 |       402.18 |  352.21 |       317.59 |      4306.00 |
| 16384 |       274.41 |  213.88 |       172.88 |      2801.00 |
| 32768 |       201.71 |  135.47 |        85.29 |      2046.00 |
| 65536 |       255.97 |  117.33 |        41.43 |      1673.00 |
