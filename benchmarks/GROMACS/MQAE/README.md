## GROMACS/CP2K QM/MM benchmark MQAE

- N-(6-methoxyquinolyl) acetoethyl ester in solution
- 16,396 atoms
- 34 QM atoms

### About

This benchmark consists of a short QM/MM simulation using the GROMACS/CP2K interface. 5 MD steps are
performed with a time step of 1 fs. The following XC functional set ups are included:

- BLYP - using DVZP-MOLOPT-GTH

`mqae.top` - The Gromacs topology file.

`mqae.ndx` - The Gromacs index file.

`mqae.gro` - The Gromacs coordinates and velocities file.

`mqae.mdp` - The Gromacs MD parameter file.

`mqae.inp` - The CP2K input file. Contains QM parameters.

`mqae_cp2k.pdb` - The pdb coordinates for CP2K.

### How to run the benchmark (8 MPI ranks with 2 OpenMP threads each)

- Change to benchmark folder

  - `cd ${CP2K_ROOT}/benchmarks/GROMACS/MQAE`

- Run with container

  - `podman run -it --rm -v ${PWD}:/mnt spack_gromacs mpiexec -n 1 gmx_mpi grompp -f mqae.mdp -p mqae.top -c mqae.gro -n mqae.ndx -qmi mqae_cp2k.inp -o mqae.tpr -maxwarn 1`
  - `podman run -it --rm -v ${PWD}:/mnt spack_gromacs mpiexec -n 8 gmx_mpi mdrun -s mqae.tpr`

- Run with launch script

  - `${CP2K_ROOT}/install/bin/launch mpiexec -n 1 gmx_mpi grompp -f mqae.mdp -p mqae.top -c mqae.gro -n mqae.ndx -qmi mqae_cp2k.inp -o mqae.tpr -maxwarn 1`
  - `${CP2K_ROOT}/install/bin/launch mpiexec -n 8 gmx_mpi mdrun -s mqae.tpr`
