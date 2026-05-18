# Install from Distribution

The easiest way to get started with CP2K is to install it from a distribution.

```{note}
The pre-compiled package from distributions may not performs as well as package
built from source, as it usually adopts a conservative compilation strategy to
ensure broad compatibility.
```

## Arch

```shell
$ pacman -S cp2k
```

See also [archlinux.org](https://aur.archlinux.org/packages/cp2k).

## Conda

[Download](https://conda-forge.org/download/) and install [conda](https://conda-forge.org/). The
available CP2K versions can be listed with

```shell
$ conda search conda-forge::cp2k
```

It is recommended to select the latest CP2K version and to install it into a new, dedicated conda
environment, e.g. `cp2k_env`, to avoid dependency conflicts.

```shell
$ conda create -n cp2k_env conda-forge::cp2k=2026.1
```

Activate the new conda environment

```shell
$ conda activate cp2k_env
$ cp2k -h -v
```

and run either OpenMP parallel (e.g. with 4 OpenMP threads)

```shell
$ export OMP_NUM_THREADS=4
$ cp2k <cp2k.inp>
```

or MPI/OpenMP parallel CP2K jobs

```shell
$ export OMP_NUM_THREADS=2
$ mpiexec -n 4 cp2k <cp2k.inp>
```

See also [conda-forge.org](https://conda-forge.org/packages/) and
[conda-metadata-app.streamlit.app](https://conda-metadata-app.streamlit.app/?q=conda-forge/cp2k).

## Debian / Ubuntu

```shell
$ apt-get install cp2k
```

See also [debian.org](https://packages.debian.org/search?keywords=cp2k) and
[ubuntu.com](https://packages.ubuntu.com/search?keywords=cp2k).

## Docker

```shell
$ docker pull cp2k/cp2k
```

See also [hub.docker.com](https://hub.docker.com/r/cp2k/cp2k) and
[cp2k-containers](https://github.com/cp2k/cp2k-containers).

## Easybuild

```shell
$ eb --software-name=CP2K --toolchain=foss,2023a
```

See also [easybuild.io](https://docs.easybuild.io/version-specific/supported-software/c/CP2K/).

## Fedora

```shell
$ dnf install cp2k
```

See also [fedoraproject.org](https://src.fedoraproject.org/rpms/cp2k)

## FreeBSD

```shell
$ pkg install cp2k
```

See also [freshports.org](https://www.freshports.org/science/cp2k/).

## Homebrew

```shell
$ brew install cp2k
```

See also [brew.sh](https://formulae.brew.sh/formula/cp2k).

## NixOS

```shell
$ nix-shell -p cp2k
```

See also [nixos.org](https://search.nixos.org/packages?query=cp2k).

## Nvidia GPU Cloud

```shell
$ docker pull nvcr.io/hpc/cp2k:v2023.2
```

See also [ngc.nvidia.com](https://catalog.ngc.nvidia.com/orgs/hpc/containers/cp2k).

## Spack

```shell
$ spack install cp2k
```

See also [spack.io](https://packages.spack.io/package.html?name=cp2k) and [](./build-with-spack).
