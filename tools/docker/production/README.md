# CP2K docker containers for production

CP2K docker containers are suited for many host systems. A convenient way to use them is to run such a container with [Apptainer](https://apptainer.org/) which is described in the following sections. Further down, you will find the information how to run and build CP2K production containers with [docker](https://docs.docker.com/).

## Apptainer (Singularity)

It is possible to download and run CP2K docker containers with [Apptainer](https://apptainer.org/) as `.sif` files. Such Apptainer files are especially suited to be run on HPC systems or systems where [docker](https://docs.docker.com/) is not installed.

### Prerequisites for the use of Apptainer

Install the latest [Apptainer](https://apptainer.org/) version as described [here](https://apptainer.org/docs/admin/latest/installation.html#installation-on-linux). Most Linux distributions provide rpm packages for Apptainer.

### Download a pre-built CP2K docker container with Apptainer

The available pre-built CP2K production docker containers can be found [here](https://hub.docker.com/repository/docker/mkrack/cp2k/tags?page=1&ordering=last_updated). The name of a docker container indicates CP2K version, MPI implementation (MPICH or OpenMPI), target CPU (`generic`, `haswell`, or `skylake-avx512`), CUDA support, and cp2k binary version (e.g. `psmp` for MPI parallelized with OpenMP support). If you don't know your CPU model, then try the CPU target `generic` (NEHALEM) first:

```
apptainer pull docker://mkrack/cp2k:2023.2_mpich_generic_psmp
```

### Run CP2K with Apptainer

Check if the pre-built container runs on your system with

```
apptainer run cp2k_2023.2_mpich_generic_psmp.sif cp2k -h -v
```

or simply

```
cp2k_2023.2_mpich_generic_psmp.sif
```

should work as well, if the `.sif` file is executable. This test should not give any error messages. A warning like `setlocale: LC_ALL: cannot change locale` will disappear when `LC_ALL=C` is set in advance:

```
LC_ALL=C cp2k_2023.2_mpich_generic_psmp.sif cp2k
```

If the `.sif` is not executable, then you can use the command `chmod a+x` to make it executable:

```
chmod a+x cp2k_2023.2_mpich_generic_psmp.sif
```

A more extensive testing, a kind of self-test, can be performed with the `run_tests` command

```
apptainer run -B $PWD:/mnt cp2k_2023.2_mpich_generic_psmp.sif run_tests
```

which launches a full [CP2K regression test](https://www.cp2k.org/dev:regtesting/) run within the container. The test run will use 8 tasks (CPU cores) by default. If you have more or less CPU cores available, you can select any multiple of 4, e.g,

```
apptainer run -B $PWD:/mnt cp2k_2023.2_mpich_generic_psmp.sif run_tests --maxtasks 32
```

### Running CUDA enabled containers with Apptainer

CUDA enabled containers can be employed on host systems which provide also NVIDIA GPU resources. Pull a CUDA enabled CP2K docker container matching the GPU type of your system, e.g. for Pascal GPUs (P100)

```
apptainer pull docker://mkrack/cp2k:2023.2_mpich_generic_cuda_P100_psmp
```

Run the container with the `--nv` flag using the `nvidia-smi` command from the container

```
apptainer run --nv ./cp2k_2023.2_mpich_generic_cuda_P100_psmp.sif nvidia-smi
```

and compare its output with the ouput of the `nvidia-smi` command from the host system. These outputs should agree when the container is able to recognize the GPU resources of the host system, correctly. Check also the CUDA version of the host system shown (top right) in the output of `nvidia-smi` command. That CUDA version should match or be newer than the CUDA version installed in the container. Run the command

```
apptainer run --nv ./cp2k_2023.2_mpich_generic_cuda_P100_psmp.sif nvcc --version
```

to retrieve the CUDA version of the actual container.  With the `run_tests` command

```
apptainer run --nv -B $PWD:/mnt ./cp2k_2023.2_mpich_generic_cuda_P100_psmp.sif run_tests
```

one can check eventually if the container works correctly and the GPU usage can be monitored on the host system with the `nvidia-smi` command while running.

### Running MPI within the container

The MPI of the container can be employed to run CP2K within a compute node, e.g. using 4 MPI ranks with 2 OpenMP threads each

```
apptainer run -B $PWD cp2k_2023.2_mpich_generic_psmp.sif mpiexec -n 4 -genv OMP_NUM_THREADS=2 cp2k -i H2O-32.inp
```

with MPICH and similarly with OpenMPI

```
apptainer run -B $PWD cp2k_2023.2_openmpi_generic_psmp.sif mpiexec -n 4 -x OMP_NUM_THREADS=2 cp2k -i H2O-32.inp
```

### Running MPI outside the container

For multi-node runs on HPC cluster systems, it is required to use the MPI of the host system

```
export OMP_NUM_THREADS=2
mpiexec -n 4 apptainer run -B $PWD cp2k_2023.2_mpich_generic_psmp.sif cp2k -i H2O-32.inp
```

to achieve best performance, but incompabilities, e.g. because of proprietary drivers or installations, might disable runing the pre-built container with the host MPI. If the host system has installed SLURM as a scheduler, `srun` can (should) be used instead of `mpiexec` (or `mpirun`)

```
srun -n 4 apptainer run -B $PWD cp2k_2023.2_mpich_generic_psmp.sif cp2k -i H2O-32.inp
```

With SLURM, `srun` is usually the proper way to launch a production run in batch mode using a CP2K `sif` file.

### CP2K tools

The containers provide also various CP2K tools like `dumpdcd`, `graph`, and `xyz2dcd`. Run

```
cp2k_2023.2_mpich_generic_psmp.sif dumpdcd -h
```

for more information and [see below](#further-binaries).

## Docker

### Prerequisites for the use of docker

Install the latest [docker](https://docs.docker.com/get-docker/) version on the host system. Most Linux distributions provide rpm packages for docker.

### Download of pre-built CP2K containers with docker

Pre-built CP2K docker containers can also be downloaded from [Dockerhub](https://hub.docker.com/r/mkrack/cp2k/tags/) using docker. There are containers compiled using MPICH and OpenMPI for various target CPUs like `generic` (NEHALEM) and for the `x86_64` Intel CPUs `haswell` and `skylake-avx512` are available.

### Run a CP2K docker container

Check if the pre-built container runs on your system using docker, e.g. with

```
docker run -it --rm mkrack/cp2k:2023.2_mpich_generic_psmp
```

or

```
docker run -it --rm mkrack/cp2k:2023.2_mpich_generic_psmp cp2k -h -v
```

This test should not give any error messages. A more extensive testing can be performed with `run_tests`

```
docker run -it --rm --shm-size=1g -v $PWD:/mnt -u $(id -u $USER):$(id -g $USER) mkrack/cp2k:2023.2_mpich_generic_psmp run_tests
```

The flag `-it` keeps the interactive terminal attached which allows for an interrupt of the run with `Ctrl-C`. The `--rm` flag removes automatically the created container when it exits.

### Running MPI within the container

The MPI of the container can be employed to run CP2K within a compute node. The containers built with MPICH can be run with

```
docker run -it --rm --shm-size=1g -v $PWD:/mnt -u $(id -u $USER):$(id -g $USER) mkrack/cp2k:2023.2_mpich_generic_psmp mpirun -n 4 -genv OMP_NUM_THREADS=2 cp2k -i H2O-32.inp
```

whereas the containers built with OpenMPI can be run, e.g. using

```
docker run -it --rm --shm-size=1g -v $PWD:/mnt -u $(id -u $USER):$(id -g $USER) mkrack/cp2k:2023.2_openmpi_generic_psmp mpirun -bind-to none -np 4 -x OMP_NUM_THREADS=2 cp2k -i H2O-32.inp
```

### Running CUDA enabled docker containers

After pulling a CUDA enabled CP2K docker container, e.g. for Pascal GPUs (P100), similarly to the [section above](#running-cuda-enabled-containers-with-apptainer), the following docker commands will check

```
docker run -it --rm --gpus all mkrack/cp2k:2023.2_mpich_generic_cuda_P100_psmp nvidia-smi

docker run -it --rm --gpus all mrack/cp2k:2023.2_mpich_generic_cuda_P100_psmp nvcc --version

docker run -it --rm --gpus all --shm-size=1g -v $PWD:/mnt -u $(id -u $USER):$(id -g $USER) mkrack/cp2k:2023.2_mpich_generic_cuda_P100_psmp run_tests
```

if the container is working correctly.

### Further binaries

All containers provide as further binaries:

- `cp2k.popt` (enforces `OMP_NUM_THREADS=1` by contrast to `cp2k.psmp` or just `cp2k`)
- `cp2k_shell` (runs the CP2K shell tool)
- `dumpdcd` (processing and conversion of `.dcd` files)
- `xyz2dcd` (convert `.xyz` files to `.dcd` files and optionally include cell information from CP2K's `.cell` file)
- `graph` (FES post-processing tool for Metadynamics runs)

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
