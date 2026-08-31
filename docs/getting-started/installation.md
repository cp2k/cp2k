# Installation

The first pivotal step with CP2K is a correct installation, comprising preparation of dependencies
and linking and building to CP2K itself, for which there are several workflows available. This page
gives an overview of considerations before linking the detailed instruction pages at the end.

For visibility, the soft link `./INSTALL.md` redirects to this page.

## Considerations

Here is a checklist for crucial premises and factors. The label "(assumed)" marks the conditions
implied for most of the documentation.

1. Which flavor of operating system is the target machine?
   - **Linux** (assumed): the workhorse for high-performance computing, a majority of compiling and
     regression tests are hosted in some Linux distribution.
   - **Darwin**: MacOS on Apple Silicon M1 is also regularly tested, and installation methods based
     on Spack or Homebrew are recommended.
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
   - **No**: if access to external public websites is restricted, unstable or outright unavailable,
     an offline installation may start with downloading all packages elsewhere, transferring them to
     designated directory on the machine and validating file integrity (such as with `sha256sum`).
1. Does the machine have a module system that manages multiple existing compilers, MPI libraries,
   math libraries and other dependent packages? (Module systems like LMod and Environment Module
   enable loading or unloading modules to control active environment variables and paths in runtime
   without conflicts; one can check with `module avail`, `module show` or commands alike.)
   - **Yes**: it is recommended to check for availability and compatibility beforehand and activate
     relevant modules for installation and execution; this not only avoids having multiple copies of
     the same package, but also get to leverage custom optimizations tailored to the machine, in
     particular hardware-dependent libfabric and MPI configurations for parallel execution.
   - **No** (assumed): if no such modules exists, or the modules are outdated or incompatible, there
     is of course the option to prepare everything from scratch with default configuration.
1. Does the machine have distinct types of nodes? (For instance, "login node" where users log in to
   perform tasks with low workload including compiling programs, and "compute node" where the more
   resource-intensive computation jobs are carried out; the key is whether program environment and
   hardware specification are consistently the same.)
   - **Yes**: heterogeneous computing with distinct CPU architecture and supported instruction sets
     require extra care for the optimization target, and effective environment variables and paths
     seen in one environment may need to be explicitly activated in another via submission scripts.
   - **No** (assumed): optimization against native (host) CPU is appropriate, and environment setup
     is easier.
1. Does the machine have some workstation GPU? (A "workstation GPU" or "professional GPU" focuses on
   the double-precision float-point arithmetic and is suitable for scientific computing, in contrast
   to a "gaming GPU" or "consumer GPU" with less optimized double-precision performance.)
   - **Yes**: some of the [eigensolvers](../technologies/eigensolvers/index),
     [accelerators](../technologies/accelerators/index) and other
     [libraries](../technologies/libraries) have GPU acceleration support ready or under development
     for matrix operations in CP2K's quantum chemical calculations; note that the requirement for
     CP2K is different from that for other molecular dynamics codes based on molecular mechanics.
   - **No** (assumed): on the safe side, a CPU-only build is always possible and is not necessarily
     slower without GPU.

## Instructions

```{toctree}
---
titlesonly:
maxdepth: 2
---
build-from-source
build-with-spack
distributions
```
