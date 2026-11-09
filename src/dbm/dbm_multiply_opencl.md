# OpenCL Backend

The OpenCL backend of DBM is maintained in
[LIBXSTREAM](https://github.com/hfp/libxstream/tree/main/samples/dbm), which is a prerequisite of
CP2K's OpenCL support anyway. LIBXSTREAM enumerates its sources and kernels for CP2K's build
(`libxstream_add_dbm` in CMake, `samples/dbm/dbm.mk` for the Makefile), and `dbm_multiply_opencl.c`
merely adapts the backend to the interface shared with the CUDA and HIP backends. The backend and
its environment variables are documented alongside the code in LIBXSTREAM.
