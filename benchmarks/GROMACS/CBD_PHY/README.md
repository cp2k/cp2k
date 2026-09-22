# GROMACS/CP2K QM/MM benchmark CBD-PHY

- PBD ID: 4O0P (adapted)
- 167,923 atoms
- 68 QM atoms

### About

This benchmark consists of a short QM/MM simulation using the GROMACS/CP2K interface. 5 MD steps are
performed with a time step of 1 fs. The following XC functional set ups are included:

- PBE - using DVZP-MOLOPT-GTH

- `CBD_PHY.top` - The GROMACS topology file

- `CBD_PHY.ndx` - The GROMACS index file

- `CBD_PHY.gro` - The GROMACS coordinates file

- `CBD_PHY.mdp` - The GROMACS MD parameter file

- `CBD_PHY_cp2k.inp` - The CP2K input file. Contains QM parameters

- `CBD_PHY_cp2k.pdb` - The PDB coordinates for CP2K

### How to run the benchmark (32 MPI ranks with 2 OpenMP threads each)

- Change to benchmark folder

  - `cd ${CP2K_ROOT}/benchmarks/GROMACS/CBD_PHY`

- Run with container

  - `podman run -it --rm -v ${PWD}:/mnt spack_gromacs mpiexec -n 1 gmx_mpi grompp -f CBD_PHY.mdp -p CBD_PHY.top -c CBD_PHY.gro -n CBD_PHY.ndx -qmi CBD_PHY_cp2k.inp -o CBD_PHY.tpr -maxwarn 1`
  - `podman run -it --rm -v ${PWD}:/mnt spack_gromacs mpiexec -n 32 gmx_mpi mdrun -s CBD_PHY.tpr`

- Run with launch script

  - `${CP2K_ROOT}/install/bin/launch mpiexec -n 1 gmx_mpi grompp -f CBD_PHY.mdp -p CBD_PHY.top -c CBD_PHY.gro -n CBD_PHY.ndx -qmi CBD_PHY_cp2k.inp -o CBD_PHY.tpr -maxwarn 1`
  - `${CP2K_ROOT}/install/bin/launch mpiexec -n 32 gmx_mpi mdrun -s CBD_PHY.tpr`
