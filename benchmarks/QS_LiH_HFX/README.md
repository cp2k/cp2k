# Quickstep Hartree Fock eXchange - LiH (Regression Tests)

Hybrid benchmark to test CP2K scaling up to 10000s of cores

## Description

This is a single-point DFT energy calculation using Quickstep GAPW (Gaussian and Augmented Plane-Waves) with hybrid Hartree-Fock exchange. It consists of a 216 atom Lithium Hydride crystal with 432 electrons in a 12.3 cubic angstrom cell. These types of calculations are generally around one hundred times the computational cost of a standard local DFT calculation, although this can be reduced using the Auxiliary Density Matrix Method (ADMM). Using OpenMP is of particular benefit here as the HFX implementation requires a large amount of memory to store partial integrals. By using several threads, fewer MPI processes share the available memory on the node and thus enough memory is available to avoid recomputing any integrals on-the-fly, improving performance.

## Files description

- [`input_bulk_B88_3.inp`](input_bulk_B88_3.inp): needed to generate an initial wfn (wave function) file for the HFX runs (this should be run once before running the actual HFX benchmark, and is not a part of the benchmark)
- [`input_bulk_HFX_3.inp`](input_bulk_HFX_3.inp): the actual input file for the HFX benchmark
- the additional files [`t_c_g.dat`](../../data/t_c_g.dat) and [`POTENTIAL`](../../data/POTENTIAL) are needed, and can be found in the `cp2k/data/` directory.

## Benchmark Requirements

To run these this benchmark, CP2K needs to be compiled with libint support (-D__LIBINT), and it is advantageous to have a OMP/MPI hybrid code (cp2k.psmp). 

## How to Run the Benchmark

1) as a preliminary step, not part of the benchmark: run `input_bulk_B88_3.inp` (5min. with 256x1 mpixomp tasks) and rename the resulting wavefunction file `LiH_bulk_3-RESTART.wfn` to `B88.wfn`

```
cp LiH_bulk_3-RESTART.wfn B88.wfn
```

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

## Notes

The amount of memory available per MPI process must be altered according to the number of MPI processes being used. If this is not done the benchmark will crash with an out of memory (OOM) error. The input file keyword `MAX_MEMORY` in `input_bulk_HFX_3.inp` needs to be changed as follows:

```
MAX_MEMORY 14000
```

should be changed to

```
MAX_MEMORY new_value
```

The new value of `MAX_MEMORY` is chosen by dividing the total amount of memory available on a node by the number of MPI processes being used per node.
If a shorter runtime is desirable, the following line in `input_bulk_HFX_3.inp`:

```
MAX_SCF 20
```

may be changed to

```
MAX_SCF 1
```

in order to reduce the maximum number of SCF cycles and hence the execution time.
If the runtime or required memory needs to be reduced so the benchmark can run on a smaller number of nodes, the OPT1 basis set can be used instead of the default OPT2. To this end, the line

```
BASIS_SET OPT2
```

in `input_bulk_B88_3.inp` and in `input_bulk_HFX_3.inp` should be changed to

```
BASIS_SET OPT1
```

## Results

### Best Configurations

The best configurations are shown below. Click the links under "Detailed Results" to see more detail.

| Machine Name | Architecture | Date       | SVN Revision | Fastest time (s) | Number of cores | Number of threads                  | Detailed results |
| ------------:| ------------:| ----------:| ------------:| ----------------:| ---------------:| ----------------------------------:| ----------------:|
| HECToR       | Cray XE6     | 21/1/2014  | 13196(*)     | 121.362          | 65536           | 8 OMP threads per MPI task	        | [hector-lih-hfx](https://www.cp2k.org/performance:hector-lih-hfx) |
| ARCHER	   | Cray XC30	  | 9/1/2014   | 13473(*)	  | 51.172	         | 49152           | 6 OMP threads per MPI task	        | [archer-lih-hfx](https://www.cp2k.org/performance:archer-lih-hfx) |
| Magnus	   | Cray XC40	  | 10/11/2014 | 14377(*)	  | 62.075	         | 24576           | 4 OMP threads per MPI task	        | [magnus-lih-hfx](https://www.cp2k.org/performance:magnus-lih-hfx) |
| Piz Daint	   | Cray XC30	  | 12/05/2015 | 15268        | 66.051	         | 32768           | 4 OMP threads per MPI task, no GPU	| [piz-daint-lih-hfx](https://www.cp2k.org/performance:piz-daint-lih-hfx) |
| Cirrus	   | SGI ICE XA	  | 24/11/2016 | 17566	      | 483.676	         | 2016            | 6 OMP threads per MPI task	        | [cirrus-lih-hfx](https://www.cp2k.org/performance:cirrus-lih-hfx) |
| Noctua	   | Cray CS500	  | 25/09/2019 | 9f58d81      | 131.290	         | 10240           | 4 OMP thread per MPI task	        | [noctua-lih-hfx](https://www.cp2k.org/performance:noctua-lih-hfx) |

(*) Prior to r14945, a bug resulted in an underestimation of the number of ERIs which should be computed (by roughly 50% for this benchmark. Therefore these results cannot be compared directly with later ones.

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
