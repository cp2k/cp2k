# OpenCL

CP2K can use OpenCL devices through the DBCSR and DBM/DBT offload paths. This supports GPUs and
other devices that provide a suitable OpenCL implementation. The OpenCL backend relies on
[LIBXS](../libraries.md#libxs-improved-performance-for-matrix-multiplication) and
[LIBXSTREAM](../libraries.md#libxstream-opencl-offload-runtime); LIBXSMM may additionally be used
through LIBXS.

OpenCL acceleration does not provide the CUDA/HIP GRID or PW GPU backends. It is therefore most
useful for the DBCSR and DBM workloads supported by the OpenCL backend.

## Installing OpenCL and preparing the runtime environment

An OpenCL build needs both development files at configuration time and an actual device driver at
runtime:

- CMake must be able to find OpenCL headers and the OpenCL loader library. On Debian and Ubuntu, the
  distribution packages are commonly named `opencl-headers` and `ocl-icd-opencl-dev`; package names
  vary on other systems.
- The loader library alone is not enough to run CP2K. The system also needs an Installable Client
  Driver (ICD) provided by the device vendor or another OpenCL implementation. Vendor SDKs and
  runtimes, such as NVIDIA CUDA, AMD ROCm, or Intel's compute runtime, may provide the necessary
  components, but their availability and installation layout are platform dependent.
- For a manually managed installation in a non-standard prefix, make OpenCL, LIBXS, and LIBXSTREAM
  discoverable by CMake. `CMAKE_PREFIX_PATH` is usually the most convenient mechanism; see
  [](../../getting-started/build-from-source.md) for the general CMake workflow.
- Set `ACC_OPENCL_VERBOSE=2` to print information about generated kernels, or `ACC_OPENCL_VERBOSE=3`
  to print information about executed kernels. These settings are useful for checking that CP2K
  reaches the intended OpenCL runtime and device.

## Building CP2K with OpenCL-based DBCSR

Configure CP2K with CMake as usual and add:

```bash
cmake -S . -B build -GNinja \
  -DCP2K_USE_ACCEL=OPENCL \
  -DCP2K_USE_LIBXS=ON \
  -DCMAKE_PREFIX_PATH=/path/to/dependencies
cmake --build build --parallel
```

`CP2K_USE_LIBXS=ON` is required for the OpenCL backend. CMake also requires OpenCL and LIBXSTREAM,
and the DBCSR installation selected by CMake must support the requested OpenCL configuration. The
CMake summary identifies the discovered OpenCL and LIBXSTREAM dependencies. No GPU architecture
needs to be specified for an OpenCL build.

CMake supplies the required compile definitions, include paths, and link dependencies. In
particular, do not duplicate these settings through manually added compiler or linker flags.

The DBCSR OpenCL kernels use tuned small-matrix-multiplication parameters where available. The
available parameter sets are embedded into the application, so the executable does not depend on a
fixed external parameter-file location. If no suitable tuned parameters are available, DBCSR uses
fallback defaults to construct kernels at runtime.

Set `OPENCL_LIBSMM_SMM_PARAMS=/path/to/csv-file` to override the built-in parameters for an existing
application, or set `OPENCL_LIBSMM_SMM_PARAMS=0` to disable them. See the
[DBCSR documentation](https://cp2k.github.io/dbcsr/) for details on kernel-parameter tuning.

## Building CP2K with the OpenCL-based DBM library

`-DCP2K_USE_ACCEL=OPENCL` also enables the OpenCL-capable DBM build path. The `CP2K_ENABLE_DBM_GPU`
option is enabled by default when an accelerator backend is selected; set
`-DCP2K_ENABLE_DBM_GPU=OFF` to build without DBM offload while retaining the rest of the selected
configuration.

During the build, CMake generates the header that embeds the DBM OpenCL kernel source and adds the
required target dependencies automatically. No additional build rule or generated source file has to
be maintained by the user.

The CUDA/HIP-specific GRID and PW acceleration options do not enable corresponding OpenCL backends.
To disable DBCSR acceleration explicitly, configure with `-DCP2K_DBCSR_USE_CPU_ONLY=ON`.
