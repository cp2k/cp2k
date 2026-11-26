# GROMACS/CP2K QM/MM benchmark CBD-PHY

- PBD ID: 4O0P (adapted)
- 167,923 atoms
- 68 QM atoms

### About

This benchmark consists of a short QM/MM simulation using the GROMACS/CP2K interface. 5 MD steps are
performed with a time step of 1 fs. The XC functional PBE with the DVZP-MOLOPT-GTH basis set is used
for the QM calculation.

This directory contains the following files:

- `CBD_PHY.top.bz2` - The compressed GROMACS topology file

- `CBD_PHY.ndx.bz2` - The compressed GROMACS index file

- `CBD_PHY.gro.bz2` - The compressed GROMACS coordinates file

- `CBD_PHY.mdp` - The GROMACS MD parameter file

- `CBD_PHY_cp2k.inp` - The CP2K input file. Contains QM parameters

- `CBD_PHY_cp2k.pdb.bz2` - The compressed file of atomic coordinates in PDB format for CP2K

### How to run the benchmark (32 MPI ranks with 2 OpenMP threads each)

- Change to benchmark folder

  - `cd ${CP2K_ROOT}/benchmarks/GROMACS/CBD_PHY`

- Use `bunzip2` to decompress the files

```shell
for i in ./*.bz2 ; do bunzip2 $i ; done
```

- Run with container

  - `podman run -it --rm -v ${PWD}:/mnt spack_gromacs mpiexec -n 1 gmx_mpi grompp -f CBD_PHY.mdp -p CBD_PHY.top -c CBD_PHY.gro -n CBD_PHY.ndx -qmi CBD_PHY_cp2k.inp -o CBD_PHY.tpr -maxwarn 1`
  - `podman run -it --rm -v ${PWD}:/mnt spack_gromacs mpiexec -n 32 gmx_mpi mdrun -s CBD_PHY.tpr`

- Run with launch script

  - `${CP2K_ROOT}/install/bin/launch mpiexec -n 1 gmx_mpi grompp -f CBD_PHY.mdp -p CBD_PHY.top -c CBD_PHY.gro -n CBD_PHY.ndx -qmi CBD_PHY_cp2k.inp -o CBD_PHY.tpr -maxwarn 1`
  - `${CP2K_ROOT}/install/bin/launch mpiexec -n 32 gmx_mpi mdrun -s CBD_PHY.tpr`
