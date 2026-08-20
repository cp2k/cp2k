# CUDA

- Use `-DCP2K_USE_ACCEL=CUDA` to generally enable support for NVIDIA GPUs
- Specify the CUDA compute capability with `-DCMAKE_CUDA_ARCHITECTURES`, for example `100` for B200.
  The deprecated `-DCP2K_WITH_GPU` selector remains available with the values K20X, K40, K80, P100,
  V100, A100, H100, B200, GB10, and A40.
- Use `-DCP2K_USE_NVHPC=ON` when building with the NVHPC kit.
- Use `-DCP2K_WITH_GPU_PROFILING` to turn on NVIDIA Tools Extensions. It requires to link
  `-lnvToolsExt`.
- Link to a BLAS/ScaLAPACK library that accelerates large DGEMMs (e.g., libsci_acc)
- Use `-DCP2K_ENABLE_GRID_GPU=OFF` to disable the GPU backend of the grid library.
- Use `-DCP2K_ENABLE_DBM_GPU=OFF` to disable the GPU backend of the sparse tensor library.
- Use `-DCP2K_ENABLE_PW_GPU=OFF` to disable the GPU backend of FFTs and associated gather/scatter
  operations.
- Use `-DCP2K_DBCSR_USE_CPU_ONLY=ON` to disable the GPU backend of DBCSR.
