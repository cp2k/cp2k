# Installation

The first pivotal step with CP2K is a correct installation, comprising preparation of dependencies
and linking and building to CP2K itself, for which there are several workflows available. This page
gives an overview of considerations, while the step-by-step instructions are kept brief; further
explanations and alternative methods have been relegated to dedicated pages linked at the end.

For visibility, the soft link `./INSTALL.md` redirects to this page.

## Considerations

Here is a checklist for crucial premises and factors. The label "(assumed)" marks the conditions
implied for most of the documentation.

1. Which flavor of operating system is the target machine?

- **Linux** (assumed): the workhorse for high-performance computing, a majority of compiling and
  regression tests are hosted in some Linux distribution.
- **Windows**: for lightweight beginner experiments, the Windows Subsystem for Linux (WSL) or any
  other virtual machine platforms can provide a Linux environment on Windows.

1. Does the user have root privilege or sudo power on the machine?

- **Yes**: proceed with caution, especially if the privilege means that actions may disrupt the
  environment for other users and softwares.
- **No** (assumed): it is totally fine to install as a normal user in a convenient directory with
  read and write permission as well as sufficient disk space.

1. Does the machine have internet connection?

- **Yes** (assumed): okay, new packages will be downloaded automatically to prepare dependencies;
  however, once ready, running CP2K itself does not require internet connection.
- **No**: if access to external public websites is restricted, unstable or outright unavailable, an
  offline installation may start with downloading all packages elsewhere, transferring them to
  designated directory on the machine and validating file integrity (such as with `sha256sum`).

1. Does the machine have a module system that manages multiple existing compilers, MPI libraries,
   math libraries and other dependent packages? (Module systems like LMod and Environment Module
   enable loading or unloading modules to control active environment variables and paths in runtime
   without conflicts; one can check with `module avail`, `module show` or commands alike.)

- **Yes**: it is recommended to check for availability and compatibility beforehand and activate
  relevant modules for installation and execution; this not only avoids having multiple copies of
  the same package, but also get to leverage custom optimizations tailored to the machine, in
  particular hardware-dependent libfabric and MPI configurations for parallel execution.
- **No** (assumed): if no such modules exists, or the modules are outdated or incompatible, there is
  of course the option to prepare everything from scratch with default configuration.

1. Does the machine have distinct types of nodes? (For instance, "login node" where users log in to
   perform tasks with low workload including compiling programs, and "compute node" where the more
   resource-intensive computation jobs are carried out; the key is whether program environment and
   hardware specification are consistently the same.)

- **Yes**: heterogeneous computing with distinct CPU architecture and supported instruction sets
  require extra care for the optimization target, and effective environment variables and paths seen
  in one environment may need to be explicitly activated in another via submission scripts.
- **No** (assumed): optimization against native (host) CPU is appropriate, and environment setup is
  easier.

1. Does the machine have some workstation GPU? (A "workstation GPU" focuses on the double-precision
   float-point arithmetic and is suitable for professional computing, in contrast to a "gaming GPU"
   with less optimized performance in this regard.)

- **Yes**: some of the [eigensolvers](../technologies/eigensolvers/index),
  [accelerators](../technologies/eigensolvers/index) and other
  [libraries](../technologies/libraries) have GPU acceleration support ready or under development
  for matrix operations in CP2K's quantum chemical calculations; note that the requirement for CP2K
  is different from that for other molecular dynamics (molecular mechanics) codes.
- **No** (assumed): on the safe side, a CPU-only build is always possible and is not necessarily
  slower without GPU.

## Instructions

Installation takes as simple as 3 steps given the assumed points above. Prerequisite tools that are
commonly available from package managers (`apt-get`, `dnf` and the like) are:

- A compiler suite compliant with the Fortran 2008 standard (see also the github page on
  [compiler support](https://github.com/cp2k/cp2k/wiki/Compiler-support));
- Utilities like `python3`, `git`, `less`, `make`, `tar`, `wget`, `zlib`, `unzip`, `bzip2` and `xz`
  as the universal infrastructure. (This is not a complete list.)

```{note}
Starting from 2026, [CMake](https://cmake.org) is the only build system supported and maintained by
CP2K. Versions from 2025 or prior used a custom system based on GNU Make and ARCH files, which has
been deprecated and removed from documentation.
```

### Step 1: obtain source code from github

- For end users, the preferred method is to go to [Releases](https://github.com/cp2k/cp2k/releases/)
  and download and extract one of the versioned tarballs `cp2k-X.Y.tar.bz2`.
- For developers, the preferred method is to fork the [`cp2k/cp2k`](https://github.com/cp2k/cp2k)
  repository and download via `git clone`.

### Step 2: set up dependencies and install CP2K

Enter the CP2K directory. Two branches of workflow are available:

- One uses the `make_cp2k.sh` script under the main directory;
- Another uses the `install_cp2k_toolchain.sh` and `build_cp2k.sh` scripts under the directory
  `tools/toolchain/`.

#### Method 1: make_cp2k.sh

[`make_cp2k.sh`](../../make_cp2k.sh) is recently added and available since the 2026.2 release. It is
based on [Spack](https://spack.io/), and it manages dependencies and builds CP2K locally within the
current working directory unless the variable `CP2K_ROOT` is set.

To display the introduction and complete list of options, run

```shell
./make_cp2k.sh --help
```

To perform by default a full build of CP2K including (almost) all features, run

```shell
./make_cp2k.sh
```

For a minimal build from scratch, run

```shell
./make_cp2k.sh -bd -df all
```

Based on which desired features can be added explicitly for building a tailored CP2K binary.

```shell
./make_cp2k.sh -bd -df all -ef libint -ef libxc -ef spglib -ef tblite
```

#### Method 2: Toolchain

On the other hand, [`install_cp2k_toolchain.sh`](../../tools/toolchain/install_cp2k_toolchain.sh)
(with a dedicated [README](../../tools/toolchain/README.md) page) is an system in existence for
years that prepares most of the dependencies by building from source. Nowadays it is followed by
[`build_cp2k.sh`](../../tools/toolchain/build_cp2k.sh), available also since the 2026.2 release,
which configures links and CMake flags for compiling CP2K.

Enter the directory, and to read help and default options, run

```shell
cd ./tools/toolchain/
./install_cp2k_toolchain.sh --help
./build_cp2k.sh --help
```

To start preparing toolchain, run

```shell
./install_cp2k_toolchain.sh
```

After that, run

```shell
./build_cp2k.sh
```

It is also available to make installed program and dependencies outside the source tree.

```shell
./install_cp2k_toolchain.sh --install-dir=/opt/cp2k/toolchain
./build_cp2k.sh --prefix /opt/cp2k
```

### Step 3: do regtesting for validation

If compilation works fine, it is recommended to test the generated CP2K binary with a collection of
input files under the `./tests` directory. The CMake build procedure concludes with a message about
how to use the python script [`do_regtest.py`](../../tests/do_regtest.py), and further notes can be
found at [](../development/regtesting).

## Further explanations and alternative methods

```{toctree}
---
titlesonly:
maxdepth: 2
---
build-from-source
build-with-spack
distributions
```
