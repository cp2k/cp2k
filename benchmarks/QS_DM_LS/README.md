# Quickstep Density Matrix Linear Scaling

## Description

This is a single-point energy calculation using linear-scaling DFT.

For large systems the linear-scaling approach for solving Self-Consistent-Field equations will be much cheaper computationally than using standard DFT and allows scaling up to 1 million atoms for simple systems. The linear scaling cost results from the fact that the algorithm is based on an iteration on the density matrix. The cubically-scaling orthogonalisation step of standard Quickstep DFT using OT is avoided and the key operation is sparse matrix-matrix multiplications, which have a number of non-zero entries that scale linearly with system size. These are implemented efficiently in the DBCSR library.

The problem size can be tuned by the parameter `NREP` in the input file, whereby the number of atoms scales cubically with `NREP`.

## Files Description

- [H2O-dft-ls.inp](H2O-dft-ls.inp) (NREP=6): H20 density functional theory linear scaling consisting of 20'736 atoms in a 59 Å³ box (6'912 water molecules in total). An LDA functional is used with a DZVP MOLOPT basis set and a 300 Ry cut-off. 
- [H2O-dft-ls.NREP4.inp](H2O-dft-ls.NREP4.inp): H20 density functional theory linear scaling consisting of 6'144 atoms in a 39 Å³ box (2'048 water molecules in total). An LDA functional is used with a DZVP MOLOPT basis set and a 300 Ry cut-off. 
- [H2O-dft-ls.NREP2.inp](H2O-dft-ls.NREP2.inp): H20 density functional theory linear scaling consisting of 6'144 atoms in a 39 Å³ box (2'048 water molecules in total). An LDA functional is used with a DZVP MOLOPT basis set and a 300 Ry cut-off (a smaller version of the H2O-dft-ls benchmark, with NREP=2, meant to run on 1 node).
- [TiO2.inp](TiO2.inp)
- [amorph.inp](amorph.inp)

## Results

### NREP=4

The best configurations are shown below. Click the links under "Detailed Results" to see more detail.

| Machine Name | Architecture | Date       | SVN Revision | Fastest time (s) | Number of Cores | Number of Threads                  | Detailed Results |
| ------------ | ------------ | ---------- | ------------ | ---------------- | --------------- | ---------------------------------- | ---------------- |
| HECToR       | Cray XE6     | 16/1/2014  | 13196        | 98.256           | 65536           | 8 OMP threads per MPI task	        | [hector-h2o-dft-ls](https://www.cp2k.org/performance:hector-h2o-dft-ls) |
| ARCHER	   | Cray XC30	  | 8/1/2014   | 13473	      | 28.476	         | 49152           | 4 OMP threads per MPI task	        | [archer-h2o-dft-ls](https://www.cp2k.org/performance:archer-h2o-dft-ls) |
| Magnus	   | Cray XC40	  | 3/12/2014  | 14377	      | 30.921	         | 24576           | 2 OMP threads per MPI task	        | [magnus-h2o-dft-ls](https://www.cp2k.org/performance:magnus-h2o-dft-ls) |
| Piz Daint	   | Cray XC30	  | 12/05/2015 | 15268	      | 27.900	         | 32768           | 2 OMP threads per MPI task, no GPU	| [piz-daint-h2o-dft-ls](https://www.cp2k.org/performance:piz-daint-h2o-dft-ls) |
| Cirrus	   | SGI ICE XA	  | 24/11/2016 | 17566	      | 543.032	         | 2016            | 2 OMP threads per MPI task	        | [cirrus-h2o-dft-ls](https://www.cp2k.org/performance:cirrus-h2o-dft-ls) |
| Noctua	   | Cray CS500	  | 25/09/2019 | 9f58d81      | 37.730	         | 10240           | 10 OMP thread per MPI task	        | [noctua-h2o-dft-ls](https://www.cp2k.org/performance:noctua-h2o-dft-ls) |

