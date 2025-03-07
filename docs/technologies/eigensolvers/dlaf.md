# DLA-Future

[DLA-Future] is a distributed linear algebra library implemented using the C++26
[`std::execution` library] as provided by [pika].

DLA-Future provides a ScaLAPACK-like Fortran interface in
[DLA-Future-Fortran](https://github.com/eth-cscs/DLA-Future-Fortran), which can be used as a drop-in
replacement for ScaLAPACK (with a subset of ScaLAPACK arguments, e.g. workspace arguments are not
present).

## Dependencies

[DLA-Future] has several dependencies. It is officially distributed via the [Spack] package manager,
which is the recommended installation option. The [DLA-Future Spack package] lists all the required
and optional dependencies.

## CMake

[DLA-Future] is enabled in CP2K with the following CMake option:

```bash
-DCP2K_USE_DLAF=ON
```

## CP2K Input File

[DLA-Future] can be selected in the CP2K input file with the
[PREFERRED_DIAG_LIBRARY](#CP2K_INPUT.GLOBAL.PREFERRED_DIAG_LIBRARY) keyword:

```
&GLOBAL
  PREFERRED_DIAG_LIBRARY DLAF
  DLAF_NEIGVEC_MIN 1024
  [...]
&END GLOBAL
```

The [DLAF_NEIGVEC_MIN](#CP2K_INPUT.GLOBAL.DLAF_NEIGVEC_MIN) indicates the minimum matrix size for
which DLA-Future is used instead of ScaLAPACK. DLA-Future is optimized for large matrices on GPUs,
and might be slower than ScaLAPACK for small matrices or CPU matrices.

### Block Size

The default block size of CP2K might be sub-optimal for DLA-Future. While DLA-Future benefits from
larger block sizes, the optimal block size depends on the CP2K calculation. The block size can be
adjusted with the following keywords:

```
&GLOBAL
  [...]
  &FM
    FORCE_BLOCK_SIZE .TRUE.
    NCOL_BLOCKS 1024
    NCOL_BLOCKS 1024
    [...]
  &END FM
&END GLOBAL
```

## Environment Variables

## pika

[DLA-Future] is built on top of [pika], a C++ library based on the C++26 [`std::execution` library]
(see [P2300] proposal), providing a CPU runtime with user-level threads, as well as integration with
CUDA/HIP, and MPI. [pika]'s behavior can be controlled by command line options, or environment
variables. Please refer to [pika]'s documentation for more information about
[controlling the number of threads and thread bindings](https://pikacpp.org/usage.html#controlling-the-number-of-threads-and-thread-bindings).

[dla-future]: https://github.com/eth-cscs/DLA-Future
[dla-future spack package]: https://packages.spack.io/package.html?name=dla-future
[p2300]: https://cplusplus.github.io/sender-receiver/execution.html
[pika]: https://pikacpp.org/
[spack]: https://spack.readthedocs.io/en/latest/
[`std::execution` library]: https://eel.is/c++draft/#exec
