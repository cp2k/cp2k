# DLA-Future

[DLA-Future] is a distributed linear algebra library implemented using C++ `std::execution` [P2300](https://cplusplus.github.io/sender-receiver/execution.html) as implemented in [pika]().

DLA-Future provides a ScaLAPACK-like Fortran interface in [DLA-Future-Fortran](https://github.com/eth-cscs/DLA-Future-Fortran), which can be used as a drop-in replacement for ScaLAPACK (with a subset of ScaLAPACK arguments, e.g. workspace arguments are not present).

## Dependencies

[DLA-Future] has several dependencies. It is officially distributed via the [Spack] package manager, which is the recommended installation option. The [DLA-Future Spack package] lists all the required and optional dependencies. 

## CMake

[DLA-Future] is enabled with the following CMake option:

```bash
-DCP2K_USE_DLAF=ON
```

## CP2K Input File

[DLA-Future] is selected in the CP2K input file with the following keyword:
```
&GLOBAL
  PREFERRED_DIAG_LIBRARY DLAF
  DLAF_NEIGVEC_MIN 1024
  [...]
&END GLOBAL
```

The `DLAF_NEIGVEC_MIN` indicates the minimum matrix size for which DLA-Future is used instead of ScaLAPACK. DLA-Future is efficient for large matrices, but it might be slower than ScaLAPACK for small matrices.

### Block Size

The default block size of CP2K might be sub-optimal for DLA-Future. While DLA-Future benefits from larger block sizes, the optimal block size depends on the CP2K calculation. The block size can be adjusted with the following keywords:
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

[DLA-Future]: https://github.com/eth-cscs/DLA-Future
