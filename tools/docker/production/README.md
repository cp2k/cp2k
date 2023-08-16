# CP2K docker containers for production

CP2K docker containers are suited for many host systems.

## Prerequisites

Install the latest [Docker](https://docs.docker.com/get-docker/) version on the host system. Most Linux distributions provide rpm packages for Docker.

## Download of pre-built CP2K containers

Pre-built CP2K Docker containers can be downloaded from [Dockerhub](https://hub.docker.com/r/mkrack/cp2k/tags/). There are containers for compiled using
MPICH and OpenMPI for various target CPUs like `generic` (NEHALEM) and for the `x86_64` Intel CPUs `haswell` and `skylake-avx512` are available.
After the download, use the command `chmod a+x` to make the `sif` file executable, e.g.

## Run a CP2K docker container

Check if the pre-built container runs on your system with

```
docker run -it --shm-size=1g -v $PWD:/mnt -u $(id -u $USER):$(id -g $USER) mkrack/cp2k:2023.2_mpich_generic_psmp
```

or

```
docker run -it --shm-size=1g -v $PWD:/mnt -u $(id -u $USER):$(id -g $USER) mkrack/cp2k:2023.2_mpich_generic_psmp cp2k -v
```

This test should not give any error messages. A more extensive testing can be performed, if you downloaded a docker container which includes the
CP2K regression tests. These docker containers are slightly larger, but allow for running a full CP2K regression test with `run_test`

```
docker run -it --shm-size=1g -v $PWD:/mnt -u $(id -u $USER):$(id -g $USER) mkrack/cp2k:2023.2_mpich_generic_psmp run_tests
```

### Running MPI within the container

The MPI of the container can be employed to run CP2K within a compute node, e.g.

```
docker run -it --shm-size=1g -v $PWD:/mnt -u $(id -u $USER):$(id -g $USER) mkrack/cp2k:2023.2_mpich_generic_psmp mpiexec -n 4 -genv OMP_NUM_THREADS=1 cp2k -i H2O-32.inp
```

### Running MPI outside the container

For multi-node runs on HPC cluster systems, it is required to use the MPI of the host system

```
export OMP_NUM_THREADS=1
mpiexec -n 4 docker run -it --shm-size=1g -v $PWD:/mnt -u $(id -u $USER):$(id -g $USER) mkrack/cp2k:2023.2_mpich_generic_psmp cp2k -i H2O-32.inp
```

to achieve best performance, but incompabilities, e.g. because of proprietary drivers or installations, might disable runing the pre-built container with the host MPI.
If the host system has installed SLURM as a scheduler, `srun` can (should) be used instead of `mpiexec` (or `mpirun`)

```
srun docker run -it --shm-size=1g -v $PWD:/mnt -u $(id -u $USER):$(id -g $USER) mkrack/cp2k:2023.2_mpich_generic_psmp cp2k -i H2O-32.inp
```

With SLURM, `srun` is usually the proper way to launch a production run in batch mode.

## Building your own CP2K docker container

You can build your own CP2K docker container, e.g. using the Docker file `Dockerfile.2023.2_mpich_generic_psmp`, by running

```
docker build -f ./Dockerfile.2023.2_mpich_generic_psmp -t $USER/cp2k:2023.2_mpich_generic_psmp .
```
