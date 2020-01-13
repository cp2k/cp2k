# Quickstep Moeller-Plesset perturbation theory to 2nd order Random-Phase Approximation - 64 H2O

## Description of Input Files

- [`H2O-64-PBE-TZ.inp`](H2O-64-PBE-TZ.inp): needed to generate an initial wfn for the SCF runs
- [`H2O-64-RI-MP2-TZ.inp`](H2O-64-RI-MP2-TZ.inp): actual RI-MP2 benchmark: the system consists of 64 water molecules in a 12.4 cubic angstrom cell. This is exactly the same system as used in the [Quickstep H2O-64](../../QS/H2O-64.inp) benchmark but using a much more accurate model, which is around 100 times more computationally demanding than standard DFT calculations.
- [`H2O-64-RI-dRPA-TZ.inp`](H2O-64-RI-dRPA-TZ.inp): actual RI-dRPA benchmark

## Additional files

- [`BASIS_H2O`](BASIS_H2O): contains the primary and auxiliary(RI) basis sets
- [`POTENTIAL_H2O`](POTENTIAL_H2O): contains the GTH pseudo potentials
- [`H2O-64.xyz`](H2O-64.xyz): geometry in xyz format

the additional files [`t_c_g.dat`](../../../data/t_c_g.dat) is needed for the RI-MP2 run, and can be found in the [cp2k/data](../../../data) directory.

## How to Run the Benchmark

1) run `H2O-64-PBE-TZ.inp`: this will generate the file `H2O-64-PBE-TZ-RESTART.wfn`, necessary for the two benchmark runs.
2) run `H2O-64-RI-MP2-TZ.inp` for the RI-MP2 benchmark.
3) and/or run `H2O-64-RI-dRPA-TZ.inp` for the RI-RPA benchmark.

## Results

### Best Configurations

The best configurations are shown below. Click the links under "Detailed Results"to see more detail.

| Machine Name | Architecture | Date       | SVN Revision | Fastest time (s) | Number of Cores | Number of Threads                  | Detailed Results |
| ------------ | ------------ | ---------- | ------------ | ---------------- | --------------- | ---------------------------------- | ---------------- |
| HECToR       | Cray XE6     | 13/1/2014  | 13196	      | 141.633          | 49152           | 8 OMP threads per MPI task	        | [hector-h2o-64-ri-mp2](https://www.cp2k.org/performance:hector-h2o-64-ri-mp2) |
| ARCHER	   | Cray XC30	  | 9/1/2014   | 13473	      | 83.945	         | 36864           | 4 OMP threads per MPI task	        | [archer-h2o-64-ri-mp2](https://www.cp2k.org/performance:archer-h2o-64-ri-mp2) |
| Magnus	   | Cray XC40	  | 4/11/2014  | 14377	      | 63.891	         | 24576           | 6 OMP threads per MPI task	        | [magnus-h2o-64-ri-mp2](https://www.cp2k.org/performance:magnus-h2o-64-ri-mp2) |
| Piz Daint	   | Cray XC30	  | 12/05/2015 | 15268	      | 48.15	         | 32768           | 8 OMP threads per MPI task, no GPU | [piz-daint-h2o-64-ri-mp2](https://www.cp2k.org/performance:piz-daint-h2o-64-ri-mp2) |
| Cirrus	   | SGI ICE XA	  | 24/11/2016 | 17566	      | 303.571	         | 2016            | 1 OMP thread per MPI task	        | [cirrus-h2o-64-ri-mp2](https://www.cp2k.org/performance:cirrus-h2o-64-ri-mp2) |
| Noctua	   | Cray CS500	  | 25/09/2019 | 9f58d81      | 82.571	         | 10240           | 2 OMP thread per MPI task	        | [noctua-h2o-64-ri-mp2](https://www.cp2k.org/performance:noctua-h2o-64-ri-mp2) |

### Results on Monte Rosa, CSCS

| Input File            | Configuration             | Total Number of Cores| Runtime |
| --------------------- | -------------------------:| --------------------:| -------:|
| H2O-64-PBE-TZ.inp     |           128 MPI x 2 OMP |                  256 | 50 s    |
| H2O-64-RI-MP2-TZ.inp  |          1024 MPI x 2 OMP |                 2048 | ~18 min |
| H2O-64-RI-dRPA-TZ.inp |          1024 MPI x 2 OMP |                 2048 | ~12 min |

*) The timings have been obtained on CRAY-XE6 (MonteRosa@CSCS)

### Results on Piz Daint, CSCS

| Input File            | Date       | CP2K Git SHA | Number of nodes | Node Configuration  | Runtime |
| --------------------- | ---------- | ------------:| ---------------:| ------------------- | ------- |
| H2O-64-PBE-TZ.inp     | 2019-08-07 | 4519a8ad7    | 4 nodes         | 12 MPI x 1 OMP      | 69 s    |
| H2O-64-RI-dRPA-TZ.inp | 2019-08-07 | 4519a8ad7    | 128 nodes       | 2 MPI x 6 OMP       | 455 s   |

*) The timings have been obtained on CRAY-XC50 (PizDaint@CSCS, GPU partition)

