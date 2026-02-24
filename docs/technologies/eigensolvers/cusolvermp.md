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

[cal]: https://developer.download.nvidia.com/compute/cublasmp/redist/libcal/
[cusolvermp]: https://docs.nvidia.com/cuda/cusolvermp/
[nccl]: https://developer.nvidia.com/nccl
[ucc]: https://github.com/openucx/ucc
