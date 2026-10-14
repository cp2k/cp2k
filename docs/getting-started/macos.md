# Install on macOS

This page describes how to run CP2K on a Mac with Apple Silicon (M1, M2, M3, M4, and later) using a
pre-built CP2K container image and [Podman](https://podman.io/). The steps are:

1. [Install Homebrew](#install-homebrew), the package manager for macOS
1. [Install Podman](#install-podman) and create a Podman machine
1. [Install CP2K](#install-cp2k) by pulling the CP2K container image
1. [Run CP2K](#run-cp2k)

The CP2K container image `docker.io/cp2k/cp2k:2026.2_mpich_aarch64_psmp` contains a Linux `aarch64`
(ARM 64-bit) build of CP2K 2026.2 with MPICH and OpenMP support (`psmp`). On macOS, Podman runs this
container in a lightweight Linux virtual machine, the so-called Podman machine, which uses the
native virtualization framework of macOS. No emulation is involved, i.e. CP2K runs at the native
speed of the Apple Silicon CPU cores. The same container image runs on all Apple Silicon Macs.

```{note}
A native macOS build of CP2K is also available via [Homebrew](./distributions.md#homebrew)
(`brew install cp2k`). The container provides a Linux environment that is identical on every
machine, which is convenient for reproducible results and for testing input files before running
them on a Linux cluster.
```

## Install Homebrew

Skip this step if Homebrew is already installed, i.e. if `brew --version` prints a version number.

Open the Terminal app and run the installation script from [brew.sh](https://brew.sh/)

```shell
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

The script asks for the password of your macOS user account, which needs administrator rights, and
installs the Xcode Command Line Tools if they are missing. On Apple Silicon, Homebrew is installed
in `/opt/homebrew`. Follow the "Next steps" printed at the end of the installation to add Homebrew
to your `PATH`.

Check the installation with

```shell
brew --version
brew doctor
```

## Install Podman

Install Podman with Homebrew

```shell
brew install podman
podman --version
```

```{note}
Alternatively, Podman can be installed with the official macOS installer package available from
[podman.io](https://podman.io/docs/installation#macos). Podman Desktop provides an optional
graphical user interface, but it is not needed for the steps below.
```

### Create a Podman machine

Containers on macOS run inside a Linux virtual machine. Its resources, i.e. the number of CPU cores,
the memory, and the disk size, are set when the machine is created. Check the number of CPU cores
and the memory (in bytes) of your Mac with

```shell
sysctl -n hw.ncpu hw.perflevel0.physicalcpu hw.perflevel1.physicalcpu hw.memsize
```

which prints the total number of CPU cores, the number of performance cores, the number of
efficiency cores, and the memory size. Then create and start a Podman machine with, e.g., all CPU
cores and 12 GiB of memory (`--memory` is given in MiB and `--disk-size` in GiB)

```shell
podman machine init --cpus $(sysctl -n hw.ncpu) --memory 12288 --disk-size 100
podman machine start
```

Leave at least 4 GiB of memory to macOS, e.g. use up to 12 GiB on a Mac with 16 GiB of memory. The
memory of the Podman machine limits the size of the systems that can be calculated with CP2K.

Check that Podman machine is up and working with

```shell
podman machine list
podman run --rm docker.io/library/hello-world
```

The architecture should be reported as `arm64`.

The resources of an existing Podman machine can be changed after stopping it, e.g.

```shell
podman machine stop
podman machine set --cpus 8 --memory 16384
podman machine start
```

```{note}
The Podman machine does not start automatically after a reboot of the Mac. Run
`podman machine start` before using CP2K, and `podman machine stop` to release its memory when
CP2K is not needed.
```

## Install CP2K

Pull the CP2K container image from Docker Hub (download size about 1.6 GB)

```shell
podman pull docker.io/cp2k/cp2k:2026.2_mpich_aarch64_psmp
```

Run the container without arguments to check the installation

```shell
podman run --rm docker.io/cp2k/cp2k:2026.2_mpich_aarch64_psmp
```

which prints the CP2K help message and the version information, including `CP2K version 2026.2`.
Other CP2K versions are available from [hub.docker.com](https://hub.docker.com/r/cp2k/cp2k/tags).
Select only tags containing `aarch64`.

### Build the container image locally

Alternatively, a CP2K container image can be built locally from the CP2K source code, e.g. to run
the current development version of CP2K. Clone the CP2K repository from GitHub and change to the
`tools/docker` directory

```shell
git clone https://github.com/cp2k/cp2k.git
cd cp2k/tools/docker
```

The Dockerfile
[Dockerfile.test_spack_psmp-gcc15](https://github.com/cp2k/cp2k/blob/master/tools/docker/Dockerfile.test_spack_psmp-gcc15)
builds CP2K (`psmp`) with GCC 15 and MPICH, and all its dependencies with
[Spack](https://spack.readthedocs.io) using `make_cp2k.sh` (see
[](./build-from-source.md#spack-based-build-via-make_cp2ksh)). The Spack packages built are stored
in a local Spack cache, which is provided by a [MinIO](https://min.io/) server running in a separate
container. Start this Spack cache container first with

```shell
./spack_cache_start.sh
```

and then build the CP2K container image with

```shell
podman build --shm-size=1g -t spack_psmp-gcc15 -f ./Dockerfile.test_spack_psmp-gcc15 ../../
```

The last argument `../../` passes the root directory of the CP2K source code as build context.
Building all dependencies from scratch takes a while (up to a few hours, depending on the number of
CPU cores of the Podman machine). Subsequent builds, e.g. after a `git pull`, are much faster, since
the Spack packages are taken from the Spack cache, if they have not changed. As a final build step,
the CP2K regression tests are run. Their summary can be checked with

```shell
podman run --rm spack_psmp-gcc15 tail -n 20 /opt/cp2k/report.log
```

The resulting image `localhost/spack_psmp-gcc15` can be used like the image downloaded from Docker
Hub, i.e. just replace `docker.io/cp2k/cp2k:2026.2_mpich_aarch64_psmp` with `spack_psmp-gcc15` in
the `podman` commands below, e.g.

```shell
podman run --rm spack_psmp-gcc15
```

The Spack cache container can be stopped with `./spack_cache_stop.sh` and its content listed with
`./spack_cache_list.sh`. The script `./spack_cache_remove.sh` removes the stopped Spack cache
container including all cached Spack packages.

```{note}
Without a running Spack cache container, the build stops with an error message. The Spack cache can
be disabled by adding the option `--build-arg SPACK_CACHE=""` to the `podman build` command. In
that case, all dependencies are rebuilt from scratch for each new build.
```

## Run CP2K

Podman shares the `/Users` directory of the Mac with the Podman machine. Therefore, the input files
must be located in a directory below `/Users`, e.g. in your home directory. As an example, create a
working directory and save the input file `h2o.inp` from [](./first-calculation) in it

```shell
mkdir -p ~/cp2k-work/h2o
cd ~/cp2k-work/h2o
```

Run CP2K with 4 MPI processes and 2 OpenMP threads per process

```shell
podman run --rm --shm-size=1g -v "$PWD:/mnt" docker.io/cp2k/cp2k:2026.2_mpich_aarch64_psmp \
  mpiexec -n 4 -genv OMP_NUM_THREADS=2 cp2k -i h2o.inp -o h2o.out
```

with

- `--rm` removes the container after CP2K has finished
- `--shm-size=1g` increases the shared memory size needed for the MPI communication
- `-v "$PWD:/mnt"` mounts the current directory in the container at `/mnt`, which is the working
  directory of CP2K within the container
- `mpiexec -n 4 -genv OMP_NUM_THREADS=2` starts 4 MPI processes with 2 OpenMP threads each using the
  MPICH of the container

All output files, like `h2o.out`, are written to the current directory and are owned by your macOS
user. The final energy can be found with

```shell
grep "Total energy:" h2o.out
```

The CP2K data files, like basis sets and pseudopotentials, are included in the container image. Add
the option `-it` to be able to interrupt a running calculation with `Ctrl-C`.

### Choice of MPI processes and OpenMP threads

The number of MPI processes times the number of OpenMP threads should not exceed the number of CPU
cores of the Podman machine. Apple Silicon CPUs have fast performance cores and slower efficiency
cores. In parallel runs, all processes wait for the slowest one. Therefore, using only the
performance cores is often as fast as or faster than using all cores.

It is worth trying a few combinations for your system, e.g. with the
[benchmarks](https://github.com/cp2k/cp2k/tree/master/benchmarks), before starting production runs.

### Regression tests

The CP2K regression tests included in the container can be run with

```shell
podman run -it --rm --shm-size=1g docker.io/cp2k/cp2k:2026.2_mpich_aarch64_psmp run_tests
```

The full test suite takes a while to complete.

## Update and Uninstall

To switch to a newer CP2K version, pull the corresponding image and remove the old one, e.g.

```shell
podman images
podman rmi docker.io/cp2k/cp2k:2026.2_mpich_aarch64_psmp
```

To remove the Podman machine including all images and containers, and then Podman itself, run

```shell
podman machine rm
brew uninstall podman
```

## Troubleshooting

- **Cannot connect to Podman**: the Podman machine is not running. Start it with
  `podman machine start`, and check its state with `podman machine list`.
- **Input file not found**: the current directory is not below `/Users`, e.g. it is in `/Volumes` on
  an external disk, which is not shared with the Podman machine by default. Move the input files to
  a directory in your home directory.
- **CP2K is killed or exits with code 137**: the Podman machine ran out of memory. Increase its
  memory with `podman machine set --memory` as shown [above](#create-a-podman-machine), or reduce
  the number of MPI processes.
- **MPI errors about shared memory** (`/dev/shm`): the `--shm-size=1g` option is missing or too
  small for the number of MPI processes.
- **Calculations are very slow**: check that the image tag contains `aarch64` and that
  `podman image inspect --format '{{.Architecture}}' <image>` prints `arm64`.

See also [](./troubleshooting), [podman.io](https://podman.io/docs/installation#macos), and
[cp2k-containers](https://github.com/cp2k/cp2k-containers).
