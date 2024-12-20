# cuSOLVERMp

NVIDIA [cuSOLVERMp] is a high-performance, distributed-memory, GPU-accelerated library that provides
tools for the solution of dense linear systems and eigenvalue problems.

## Dependencies

- [CAL]: requires `libcal.\*` in the `$PATH`
- [UCC]: requires `libucc.\*` and `libucs.\*` in the `$PATH`

## CMake

[cuSOLVERMp] is enabled with the following CMake option:

```bash
-DCP2K_USE_CUSOLVER_MP=ON
```

[cal]: https://developer.download.nvidia.com/compute/cublasmp/redist/libcal/
[cusolvermp]: https://docs.nvidia.com/cuda/cusolvermp/
[ucc]: https://github.com/openucx/ucc
