# CP2K & Apptainer (Singularity)

CP2K Apptainers are especially suited to be run on HPC systems.

## Prerequisites

Install the latest [Apptainer](https://apptainer.org/) version as described [here](https://apptainer.org/docs/admin/latest/installation.html#installation-on-linux). Most Linux distributions provide rpm packages for Apptainer.

## Download of pre-built CP2K containers

Apptainer `sif` files pre-built with the definition files in this folder can be downloaded for CP2K version 2023.2. There are `sif` files compiled using MPICH for a `generic` (NEHALEM)

```
wget https://gitlab.psi.ch/krack/cp2k/-/raw/master/apptainer/download/cp2k-2023.2_mpich_generic_psmp.sif
```

and for the `x86_64` Intel CPUs `haswell`

```
wget https://gitlab.psi.ch/krack/cp2k/-/raw/master/apptainer/download/cp2k-2023.2_mpich_haswell_psmp.sif
```

and `skylake-avx512`

```
wget https://gitlab.psi.ch/krack/cp2k/-/raw/master/apptainer/download/cp2k-2023.2_mpich_skylake-avx512_psmp.sif
```

are available. After the download, use the command `chmod a+x` to make the `sif` file executable, e.g.

```
chmod a+x cp2k-2023.2_mpich_generic_psmp.sif
```

## Run CP2K with Apptainer

Check if the pre-built container runs on your system with

```
apptainer run -B $PWD cp2k-2023.2_mpich_generic_psmp.sif cp2k -h -v
```

or simply

```
cp2k-2023.2_mpich_generic_psmp.sif cp2k -h
```

should work as well, if the `sif` file is executable. This test should not give any error messages.
A more extensive testing can be performed with the `test` command

```
apptainer test -B $PWD:/mnt cp2k-2023.2_mpich_generic_psmp.sif
```

which launches a full CP2K regression test run as it is defined in the `%test` section of the Apptainer definition file.

### Running MPI within the container

The MPI of the container can be employed to run CP2K within a compute node, e.g.

```
apptainer run -B $PWD cp2k-2023.2_mpich_generic_psmp.sif mpiexec -n 4 -genv OMP_NUM_THREADS=1 cp2k -i H2O-32.inp
```

### Running MPI outside the container

For multi-node runs on HPC cluster systems, it is required to use the MPI of the host system

```
export OMP_NUM_THREADS=1
mpiexec -n 4 apptainer run -B $PWD cp2k-2023.2_mpich_generic_psmp.sif cp2k -i H2O-32.inp
```

to achieve best performance, but incompabilities, e.g. because of proprietary drivers or installations, might disable runing the pre-built container with the host MPI. If the host system has installed SLURM as a scheduler, `srun` can (should) be used instead of `mpiexec` (or `mpirun`)

```
srun apptainer run -B $PWD cp2k-2023.2_mpich_generic_psmp.sif cp2k -i H2O-32.inp
```

With SLURM, `srun` is usually the proper way to launch a production run in batch mode using a CP2K `sif` file.

## Building your own CP2K container with apptainer

Each Apptainer definition file in this folder provides a usage description which is the command to build the corresponding `sif` file, e.g.

```
apptainer build -B $PWD:/mnt cp2k-2023.2_mpich_generic_psmp.sif 2023.2_mpich_generic_psmp.def | tee cp2k-2023.2_mpich_generic_psmp.log
```
