# OpenCL

OpenCL devices are currently supported for DBCSR and DBM/DBT, and can cover GPUs and other devices.
Kernels can be automatically tuned.

Note: the OpenCL backend uses some functionality from LIBXSMM (dependency). CP2K's offload-library
serving DBM/DBT and other libraries depends on DBCSR's OpenCL backend.

## Installing OpenCL and preparing the runtime environment

- Installing an OpenCL runtime depends on the operating system and the device vendor. Debian for
  instance brings two packages called `opencl-headers` and `ocl-icd-opencl-dev` which can be present
  in addition to a vendor-specific installation. The OpenCL header files are only necessary if
  CP2K/DBCSR is compiled from source. Please note, some implementations ship with outdated OpenCL
  headers which can prevent using latest features (if an application discovers such features only at
  compile-time). When building from source, for instance `libOpenCL.so` is sufficient at link-time
  (ICD loader). However, an Installable Client Driver (ICD) is finally necessary at runtime.
- NVIDIA CUDA, AMD HIP, and Intel OneAPI are fully equipped with an OpenCL runtime (if
  `opencl-headers` package is not installed, CPATH can be needed to point into the former
  installation, similarly `LIBRARY_PATH` for finding `libOpenCL.so` at link-time). Installing a
  minimal or stand-alone OpenCL is also possible, e.g., following the instructions for Debian (or
  Ubuntu) as given for every [release](https://github.com/intel/compute-runtime/releases) of the
  [Intel Compute Runtime](https://github.com/intel/compute-runtime).
- The environment variable `ACC_OPENCL_VERBOSE` prints information at runtime of CP2K about kernels
  generated (`ACC_OPENCL_VERBOSE=2`) or executed (`ACC_OPENCL_VERBOSE=3`) which can be used to check
  an installation.

## Building CP2K with OpenCL-based DBCSR

- CP2K's toolchain supports `--enable-opencl` to select DBCSR's OpenCL backend. This can be combined
  with `--enable-cuda` (`--gpu-ver` is then imposed) to use a GPU for CP2K's GRID and PW components
  (no OpenCL support yet) with DBM's CUDA implementation to be preferred.
- For manually writing an ARCH-file, add `-D__OPENCL` and `-D__DBCSR_ACC` to `CFLAGS` and add
  `-lOpenCL` to the `LIBS` variable, i.e., `OFFLOAD_CC` and `OFFLOAD_FLAGS` can duplicate `CC` and
  `CFLAGS` (no special offload compiler needed). Please also set `OFFLOAD_TARGET = opencl` to enable
  the OpenCL backend in DBCSR. For OpenCL, it is not necessary to specify a GPU version (e.g.,
  `GPUVER = V100` would map/limit to `exts/dbcsr/src/acc/opencl/smm/params/tune_multiply_V100.csv`).
  In fact, `GPUVER` limits tuned parameters to the specified GPU, whereas by default all tuned
  parameters are embedded (`exts/dbcsr/src/acc/opencl/smm/params/*.csv`) and applied at runtime. If
  auto-tuned parameters are not available for DBCSR, well-chosen defaults will be used to populate
  kernels at runtime.
- Auto-tuned parameters are embedded into the binary, i.e., CP2K does not rely on a hard-coded
  location. Setting `OPENCL_LIBSMM_SMM_PARAMS=/path/to/csv-file` environment variable can supply
  parameters for an already built application, or `OPENCL_LIBSMM_SMM_PARAMS=0` can disable using
  tuned parameters. Refer to <https://cp2k.github.io/dbcsr/> on how to tune kernels (parameters).

## Building CP2K with OpenCL-based DBM library

- Pass `-DCP2K_USE_ACCEL=OPENCL` to CMake in addition to following above instructions for "Building
  CP2K with OpenCL-based DBCSR". An additional Makefile rule can be necessary to transform OpenCL
  code into a ressource header file.
