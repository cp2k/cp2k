# cuSOLVERMp

NVIDIA [cuSOLVERMp] is a high-performance, distributed-memory, GPU-accelerated library that provides
tools for the solution of dense linear systems and eigenvalue problems.

## Dependencies

[cuSOLVERmp] < 0.7

- [CAL]: requires `libcal.\*` in the `$PATH`
- [UCC]: requires `libucc.\*` and `libucs.\*` in the `$PATH`

[cuSOLVERmp] >= 0.7

- [NCCL]: requires `libnccl.\*` in the `$PATH`
- [UCC]: requires `libucc.\*` and `libucs.\*` in the `$PATH`

## CMake

[cuSOLVERMp] is enabled with the following CMake option:

```bash
-DCP2K_USE_CUSOLVER_MP=ON
```

The `FindCuSolverMP.cmake` module tries to automatically deduce the [cuSOLVERmp] version from the
`cusolvermp.h` header file, and enables the `CP2K_CUSOLVERMP_USE_NCCL` CMake option in case it finds
cuSOLVERmp >= 0.7.

## Generalized diagonalization

The `DIRECT_GENERALIZED_DIAGONALIZATION` global input keyword enables direct generalized
diagonalization paths that avoid a CP2K-side Cholesky reduction where the selected eigensolver
supports them. With [cuSOLVERMp], this currently covers real symmetric and complex Hermitian
generalized eigenproblems through `cusolverMpSygvd`.

When [cuSOLVERMp] uses the NCCL backend, CP2K requires at most one local MPI rank per visible GPU.
Use one MPI rank per GPU or reduce the number of local MPI ranks when only one GPU is visible.

[cal]: https://developer.download.nvidia.com/compute/cublasmp/redist/libcal/
[cusolvermp]: https://docs.nvidia.com/cuda/cusolvermp/
[nccl]: https://developer.nvidia.com/nccl
[ucc]: https://github.com/openucx/ucc
