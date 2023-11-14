# Spack

CP2K can be built and installed with [Spack]. [Spack] is a package manager designed to support
multiple versions and configurations on a wide variety of platforms and environments, with focus on
HPC.

To install CP2K with [Spack], you need to [install Spack].

## Install CP2K with Spack

A barebone version of CP2K can be installed with [Spack] as follows:

```bash
spack install cp2k
```

This command will build CP2K, as well as all the necessary dependencies. If it is the first time you
run this command, building all the dependencies might take a while.

In order to use CP2K installed with [Spack], you can simply type

```bash
spack load cp2k
```

This command will add the appropriate directories to `PATH` and `MANPATH`. See
[using installed Spack packages] for more information.

### Customizing Installation

#### Variants

[Spack] allows to fully customize an installation. The [CP2K Spack package] has several options for
customization (called "variants", see [Spack package variants]). For example, to install CP2K with
`libint` and `libxc`, one can type

```bash
spack install cp2k +libint +libxc
```

`+VARIANT` will enable a boolean variant, while `~VARIANT` will disable a boolean variant. Non
boolean variants can be specified with the `VARIANT=VALUE` syntax (see [Spack package variants] for
more details). The previous installation command takes care of building CP2K with `libint` and
`libxc` support. More importantly, it takes care of building the appropriate versions of `libint`
and `libxc` to work with CP2K (Fortran support, ...).

#### Versions

Versions in [Spack] can be specified with `@` following the package name (see \[version
specifier\]). The following installs version `2023.2` of CP2K:

```bash
spack install cp2k@2023.2
```

A more complete installation of CP2K can be installed with the following:

```bash
spack install cp2k@2023.2 +libint +libxc +dlaf +sirius +cosma +spglib lmax=6 
```

The `cp2k@2023.2 +libint +libxc +dlaf +sirius +cosma +spglib lmax=6` string is called a [spec] in
[Spack] lingo.

#### CUDA and ROCm Support

The `cuda` and `rocm` variants are available for CP2K. Therefore, CUDA support can be enabled with
`+cuda` and ROCm support can be enabled with `+rocm`. `cuda_arch` and `amdgpu_target` allow
specifying the GPU architecture:

```bash
spack install cp2k +cuda cuda_arch=80
spack install cp2k +rocm amdgpu_target=gfx90a
```

[Spack] is designed to support the installation of different versions of the same software,
therefore there is no problem with running _both_ commands above. However, `spack load cp2k` will no
longer work, you will need to be a bit more specific:

```bash
spack load cp2k +cuda
```

### Managing Dependencies

Sometimes you need to control dependencies too. Dependencies are also [Spack] packages, and their
installation can be configured in the same way as for CP2K. A dependency spec is defined by `^`.

For example, if you want to install CP2K with CUDA support but [DBCSR] without CUDA support you can
do

```bash
spack install cp2k +cuda cuda_arch=80 ^dbcsr ~cuda
```

Another example, is the choice of vendor libraries. In order to use Intel oneAPI MKL, it can be
specified as a dependency:

```bash
spack install cp2k ^intel-oneapi-mkl +cluster
```

This command will build CP2K with Intel oneAPI MKL. `+cluster` is a variant of the
[Intel oneAPI Spack package] enabling cluster support (ScaLAPACK, BLACS, ...).

## Developer Workflow

[Spack] has support for developer workflows. One way to use [Spack] to develop CP2K, is to use
[Spack] to install all the necessary dependencies and then build CP2K manually.

To install only the dependencies of CP2K (and not CP2K itself) you can use

```bash
spack install --only=dependencies CP2K_SPEC
```

Once all dependencies are installed, you can get a shell with all the dependencies set up:

```bash
spack build-env CP2K_SPEC -- bash
```

From here, you can iterate CP2K development using the native build system as Spack would use it.

[cp2k spack package]: https://packages.spack.io/package.html?name=cp2k
[dbcsr]: https://cp2k.github.io/dbcsr/develop/
[install spack]: https://spack.readthedocs.io/en/latest/getting_started.html#installation
[intel oneapi spack package]: https://packages.spack.io/package.html?name=intel-oneapi-mkl
[spack]: https://spack.readthedocs.io/en/latest/
[spack package variants]: https://spack.readthedocs.io/en/latest/basic_usage.html#variants
[spec]: https://spack.readthedocs.io/en/latest/basic_usage.html#specs-dependencies
[using installed spack packages]: https://spack.readthedocs.io/en/latest/basic_usage.html#using-installed-packages
