# CP2K docker containers for production

CP2K docker containers are suited for many host systems.

### Prerequisites

Install the latest [docker](https://docs.docker.com/get-docker/) version on the host system. Most Linux distributions provide rpm packages for docker.

If you do not want to install docker or if you cannot install it, because your computer center does not allow for it, alternatively, you can download the docker containers with **apptainer** (formerly singularity) using

```
apptainer pull docker://mkrack/cp2k:2023.2_mpich_generic_psmp
```

After the download make sure that `.sif` file is executable. Use the command `chmod a+x <tag name>.sif` to make the `.sif` file executable. Check the apptainer [README](https://github.com/mkrack/cp2k/tree/master/tools/apptainer/README.md) of CP2K for more details about apptainer, its installation and how to [run](https://github.com/mkrack/cp2k/tree/master/tools/apptainer#run-cp2k-with-apptainer) the `.sif` container files.

### Download of pre-built CP2K containers with docker

Pre-built CP2K docker containers can also be downloaded from [Dockerhub](https://hub.docker.com/r/mkrack/cp2k/tags/) using docker. There are containers compiled using MPICH and OpenMPI for various target CPUs like `generic` (NEHALEM) and for the `x86_64` Intel CPUs `haswell` and `skylake-avx512` are available.

### Run a CP2K docker container

Check if the pre-built container runs on your system using docker, e.g. with

```
docker run -it --shm-size=1g -v $PWD:/mnt -u $(id -u $USER):$(id -g $USER) mkrack/cp2k:2023.2_mpich_generic_psmp
```

or

```
docker run -it --shm-size=1g -v $PWD:/mnt -u $(id -u $USER):$(id -g $USER) mkrack/cp2k:2023.2_mpich_generic_psmp cp2k -v
```

This test should not give any error messages. A more extensive testing can be performed, if you downloaded a docker container which includes the
CP2K regression tests. These docker containers are slightly larger, but allow for running a full CP2K regression test with `run_tests`

```
docker run -it --shm-size=1g -v $PWD:/mnt -u $(id -u $USER):$(id -g $USER) mkrack/cp2k:2023.2_mpich_generic_psmp run_tests
```

### Running MPI within the container

The MPI of the container can be employed to run CP2K within a compute node. The containers built with MPICH can be run with

```
docker run -it --shm-size=1g -v $PWD:/mnt -u $(id -u $USER):$(id -g $USER) mkrack/cp2k:2023.2_mpich_generic_psmp mpirun -n 4 -genv OMP_NUM_THREADS=1 cp2k -i H2O-32.inp
```

whereas the containers built with OpenMPI can be run, e.g. using

```
docker run -it --shm-size=1g -v $PWD:/mnt -u $(id -u $USER):$(id -g $USER) mkrack/cp2k:2023.2_openmpi_generic_psmp mpirun -bind-to none -np 4 -x OMP_NUM_THREADS=1 cp2k -i H2O-32.inp
```

All containers provide as further binaries:

- `cp2k.popt` (enforces `OMP_NUM_THREADS=1` by contrast to `cp2k.psmp` or just `cp2k`)
- `cp2k_shell` (runs the CP2K shell tool)
- `dumpdcd` (processing and conversion of `.dcd` files)
- `xyz2dcd` (convert `.xyz` files to `.dcd` files)
- `graph` (FES post-processing tool for Metadynamics runs)

### Running MPI outside the container

For multi-node runs on HPC cluster systems, it is required to use the MPI of the host system

```
export OMP_NUM_THREADS=1
mpiexec -n 4 docker run -it --shm-size=1g -v $PWD:/mnt -u $(id -u $USER):$(id -g $USER) mkrack/cp2k:2023.2_mpich_generic_psmp cp2k -i H2O-32.inp
```

to achieve best performance, but incompabilities, e.g. because of proprietary drivers or installations, might disable runing the pre-built container with the host MPI.
If the host system has installed SLURM as a scheduler, `srun` can (should) be used instead of `mpirun` (or `mpiexec`)

```
srun docker run -it --shm-size=1g -v $PWD:/mnt -u $(id -u $USER):$(id -g $USER) mkrack/cp2k:2023.2_mpich_generic_psmp cp2k -i H2O-32.inp
```

With SLURM, `srun` is usually the proper way to launch a production run in batch mode.

### Building your own CP2K docker container

You can build your own CP2K docker container, e.g. using the Docker file `Dockerfile.2023.2_mpich_generic_psmp` for the CP2K version 2023.2, by running

```
docker build -f ./Dockerfile.2023.2_mpich_generic_psmp -t $USER/cp2k:2023.2_mpich_generic_psmp .
```

or for the current CP2K trunk version (see [master](https://github.com/cp2k/cp2k/tree/master) branch) with

```
docker build -f ./Dockerfile.master_mpich_generic_psmp -t cp2k/cp2k:master$(date +%Y%m%d)_mpich_generic_psmp .
```

Each Docker file includes a `Usage` line in its header with the `docker build` command to run that file.

If the docker files will include any CP2K regression test during its build, the log file size will be much larger. In that case, it is recommended either to use docker's legacy builder

```
DOCKER_BUILDKIT=0 docker build -f ./Dockerfile.2023.2_mpich_generic_psmp -t cp2k/cp2k:2023.2_mpich_generic_psmp .
```

or to add the `--progress` flag when using the new docker builder `buildx`, e.g.

```
docker buildx build --progress=plain -f ./Dockerfile.2023.2_mpich_generic_psmp -t cp2k/cp2k:2023.2_mpich_generic_psmp .
```

### (Re-)Generate the docker files

The Python script [`generate_docker_files`](https://github.com/cp2k/cp2k/tree/master/tools/docker/production/generate_docker_files.py) generates all docker files in this [folder](https://github.com/cp2k/cp2k/tree/master/tools/docker/production/) if run without any further option. It provides a few command line options allowing for an adaptation of the created docker files to the actual needs, e.g. with

```
./generate_docker_files.py -j 28 --test
```

it will create docker files using 28 instead of the default 8 CPU cores for building the docker containers and the container will run a full CP2K regression test after the container has been built. Run `generate_docker_files -h` for more details about the available options.
