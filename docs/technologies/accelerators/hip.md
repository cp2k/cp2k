# HIP / ROCm

The code for the HIP based grid backend was developed and tested on Mi100 but should work out of the
box on NVIDIA hardware as well.

- Use `-DCP2K_USE_ACCEL=HIP` to generally enable support for AMD GPUs
- Use `-DCP2K_ENABLE_GRID_GPU=OFF` to disable the GPU backend of the grid library.
- Use `-DCP2K_ENABLE_DBM_GPU=OFF` to disable the GPU backend of the sparse tensor library.
- Use `-DCP2K_ENABLE_PW_GPU=OFF` to disable the GPU backend of FFTs and associated gather/scatter
  operations.
- Use `-DCP2K_DBCSR_USE_CPU_ONLY=ON` to disable the GPU backend of DBCSR.
- Add `-DCP2K_USE_UNIFIED_MEMORY=ON` to enable unified memory support (experimental and only
  supports Mi250X and above)
- Add `-DCP2K_WITH_GPU==Mi50, Mi60, Mi100, Mi250, Mi300`. Architectures supported include Mi300(A,X)
  (gfx1103), Mi250 (gfx90a), Mi100 (gfx908), and Mi50 (gfx906). The HIP backend for the grid library
  supports NVIDIA hardware as well. It uses the same code and can be used to validate the backend in
  case only NVIDIA hardware is available.
- Use `-DCP2K_WITH_GPU_PROFILING` to turn on the AMD ROC TX and Tracer libray. It requires to link
  `-lroctx64 -lroctracer64`.

For comprehensive ROCm documentation, see: <https://rocm.docs.amd.com/en/latest/>.
