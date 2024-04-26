# CP2K

[![Release Status][release-badge]][release-link] [![Debian Status][debian-badge]][debian-link]
[![Fedora Status][fedora-badge]][fedora-link] [![Ubuntu Status][ubuntu-badge]][ubuntu-link]
[![Arch Status][arch-badge]][arch-link] [![Homebrew Status][homebrew-badge]][homebrew-link]
[![Docker Status][docker-badge]][docker-link] [![Spack Status][spack-badge]][spack-link]
[![Conda Status][conda-badge]][conda-link]

CP2K is a quantum chemistry and solid state physics software package that can perform atomistic
simulations of solid state, liquid, molecular, periodic, material, crystal, and biological systems.
CP2K provides a general framework for different modeling methods such as DFT using the mixed
Gaussian and plane waves approaches GPW and GAPW. Supported theory levels include DFT, MP2, RPA, GW,
tight-binding (xTB, DFTB), semi-empirical methods (AM1, PM3, PM6, RM1, MNDO, ...), and classical
force fields (AMBER, CHARMM, ...). CP2K can do simulations of molecular dynamics, metadynamics,
Monte Carlo, Ehrenfest dynamics, vibrational analysis, core level spectroscopy, energy minimization,
and transition state optimization using NEB or dimer method.

CP2K is written in Fortran 2008 and can be run efficiently in parallel using a combination of
multi-threading, MPI, and CUDA.

## Downloading CP2K source code

To clone the current master (development version):

```shell
git clone --recursive https://github.com/cp2k/cp2k.git cp2k
```

Note the `--recursive` flag that is needed because CP2K uses git submodules.

To clone a release version v*x.y*:

```shell
git clone -b support/vx.y --recursive https://github.com/cp2k/cp2k.git cp2k
```

For more information on downloading CP2K, see [Downloading CP2K](https://www.cp2k.org/download). For
help on git, see [Git Tips & Tricks](https://github.com/cp2k/cp2k/wiki/Git-Tips-&-Tricks).

## Install CP2K

The easiest way to build CP2K with all of its dependencies is as a
[Docker container](./tools/docker/README.md).

For building CP2K from scratch see the [installation instructions](./INSTALL.md).

## Links

- [CP2K.org](https://www.cp2k.org) for showcases of scientific work, tutorials, exercises,
  presentation slides, etc.
- [The manual](https://manual.cp2k.org/) with descriptions of all the keywords for the CP2K input
  file
- [The dashboard](https://dashboard.cp2k.org) to get an overview of the currently tested
  architectures
- [The Google group](https://groups.google.com/group/cp2k) to get help if you could not find an
  answer in one of the previous links
- [Acknowledgements](https://www.cp2k.org/funding) for list of institutions and grants that help to
  fund the development of CP2K

## Directory organization

- [`arch`](./arch): Collection of definitions for different architectures and compilers
- [`benchmarks`](./benchmarks): Inputs for benchmarks
- [`data`](./data): Simulation parameters e.g. basis sets and pseudopotentials
- [`exts`](./exts): Access to external libraries via GIT submodules
- [`src`](./src): The source code
- [`tests`](./tests): Inputs for tests and regression tests
- [`tools`](./tools): Mixed collection of useful scripts related to cp2k

Additional directories created during build process:

- `lib`: Libraries built during compilation
- `obj`: Objects and other intermediate compilation-time files
- `exe`: Where the executables will be located

[arch-badge]: https://img.shields.io/aur/version/cp2k
[arch-link]: https://aur.archlinux.org/packages/cp2k
[conda-badge]: https://img.shields.io/conda/vn/conda-forge/cp2k
[conda-link]: https://anaconda.org/conda-forge/cp2kcp2k
[debian-badge]: https://img.shields.io/debian/v/cp2k
[debian-link]: https://packages.debian.org/search?keywords=cp2k
[docker-badge]: https://img.shields.io/docker/v/cp2k/cp2k?label=docker
[docker-link]: https://hub.docker.com/r/cp2k/cp2k
[fedora-badge]: https://img.shields.io/fedora/v/cp2k
[fedora-link]: https://src.fedoraproject.org/rpms/cp2k
[homebrew-badge]: https://img.shields.io/homebrew/v/cp2k
[homebrew-link]: https://formulae.brew.sh/formula/cp2k
[release-badge]: https://img.shields.io/github/v/release/cp2k/cp2k
[release-link]: https://github.com/cp2k/cp2k/releases
[spack-badge]: https://img.shields.io/spack/v/cp2k
[spack-link]: https://packages.spack.io/package.html?name=cp2k
[ubuntu-badge]: https://img.shields.io/ubuntu/v/cp2k
[ubuntu-link]: https://packages.ubuntu.com/search?keywords=cp2k
