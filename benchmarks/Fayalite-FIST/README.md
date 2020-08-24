# Fayalite-FIST

## Description

This is a short molecular dynamics run of 1'000 time steps in a NPT ensemble at
300K. It consists of 28'000 atoms - a 103 supercell with 28 atoms of iron silicate
(Fe2SiO4, also known as Fayalite) per unit cell. The simulation employs a classical
potential (Morse with a hard-core repulsive term and 5.5 angstrom cutoff) with
long-range electrostatics using Smoothed Particle Mesh Ewald (SPME) summation.
While CP2K does support classical potentials via the Frontiers In Simulation
Technology (FIST) module, this is not a typical calculation for CP2K but is
included to give an impression of the performance difference between machines
for the MM part of a QM/MM calculation.

## Benchmarks

- [`fayalite.inp`](fayalite.inp)

## Results

The best configurations are shown below.
Click the links under "Detailed Results" to see more detail.

<!-- markdownlint-disable MD013 -->
| Machine Name | Architecture | Date       | SVN Revision | Fastest time (s) | Number of Cores | Number of Threads                 | Detailed Results |
| ------------ | ------------ | ---------- | ------------ | ---------------- | --------------- | --------------------------------- | ---------------- |
| HECToR       | Cray XE6     | 21/1/2014  | 13196        | 403.928          | 512 cores       | 2 OMP threads per MPI task        | [hector-h2o-64](https://www.cp2k.org/performance:hector-h2o-64) |
| ARCHER       | Cray XC30    | 9/1/2014   | 13473        | 197.117          | 576 cores       | 1 OMP thread per MPI task         | [archer-h2o-64](https://www.cp2k.org/performance:archer-h2o-64) |
| Magnus       | Cray XC40    | 6/11/2014  | 14377        | 150.493          | 384 cores       | 1 OMP thread per MPI task         | [magnus-h2o-64](https://www.cp2k.org/performance:magnus-h2o-64) |
| Piz Daint    | Cray XC30    | 12/05/2015 | 15268        | 207.972          | 192 cores       | 1 OMP thread per MPI task, no GPU | [piz-daint-h2o-64](https://www.cp2k.org/performance:piz-daint-h2o-64) |
| Cirrus       | SGI ICE XA   | 24/11/2016 | 7566         | 166.192          | 1152 cores      | 9 OMP threads per MPI task        | [cirrus-h2o-64](https://www.cp2k.org/performance:cirrus-h2o-64) |
| Noctua       | Cray CS500   | 25/09/2019 | 9f58d81      | 119.820          | 640 cores       | 10 OMP thread per MPI task        | [noctua-h2o-64](https://www.cp2k.org/performance:noctua-h2o-64) |
<!-- markdownlint-enable MD013 -->
