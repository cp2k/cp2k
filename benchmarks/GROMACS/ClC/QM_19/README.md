# GROMACS/CP2K QM/MM benchmark ClC-19

- PBD ID: 1KPK
- 150,925 atoms
- 19 QM atoms

### About

This benchmark consists of a short QM/MM simulation using the GROMACS/CP2K interface. 5 MD steps are
performed with a time step of 1 fs. The XC functional BLYP with the DVZP-MOLOPT-GTH basis set is
used for the QM calculation.

This directory contains the following files:

- `ClC.top.bz2` - The compressed GROMACS topology file

- `ClC.ndx.bz2` - The compressed GROMACS index file

- `ClC.gro.bz2` - The compressed GROMACS coordinates file

- `ClC.mdp` - The GROMACS MD parameter file

- `ClC_cp2k.inp` - The CP2K input file. Contains QM parameters

- `ClC_cp2k.pdb.bz2` - The compressed file of atomic coordinates in PDB format for CP2K

### How to run the benchmark (32 MPI ranks with 2 OpenMP threads each)

- Change to benchmark folder

  - `cd ${CP2K_ROOT}/benchmarks/GROMACS/ClC/QM_19`

- Use `bunzip2` to decompress the files

```shell
for i in ./*.bz2 ; do bunzip2 $i ; done
```

- Run with container

  - `podman run -it --rm -v ${PWD}:/mnt spack_gromacs mpiexec -n 1 gmx_mpi grompp -f ClC.mdp -p ClC.top -c ClC.gro -n ClC.ndx -qmi ClC_cp2k.inp -o ClC.tpr -maxwarn 1`
  - `podman run -it --rm -v ${PWD}:/mnt spack_gromacs mpiexec -n 32 gmx_mpi mdrun -s ClC.tpr`

- Run with launch script

  - `${CP2K_ROOT}/install/bin/launch mpiexec -n 1 gmx_mpi grompp -f ClC.mdp -p ClC.top -c ClC.gro -n ClC.ndx -qmi ClC_cp2k.inp -o ClC.tpr -maxwarn 1`
  - `${CP2K_ROOT}/install/bin/launch mpiexec -n 32 gmx_mpi mdrun -s ClC.tpr`
