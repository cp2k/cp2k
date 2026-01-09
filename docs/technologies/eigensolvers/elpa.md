# ELPA

The [ELPA] library provides highly efficient and highly scalable direct eigensolvers for symmetric
(hermitian) matrices.

## CMake

[ELPA] is enabled with the following CMake option:

```bash
-DCP2K_USE_ELPA=ON
```

Pass `-DCP2K_ENABLE_ELPA_OPENMP_SUPPORT=ON` to CMake if ELPA was built with OpenMP support. Pass
`-DCP2K_ENABLE_ELPA_GPU=OFF` to CMake if ELPA was not built with GPU support.

[elpa]: https://elpa.mpcdf.mpg.de/
