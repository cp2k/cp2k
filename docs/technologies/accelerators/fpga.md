# FPGA

```{warning}
Support for FPGAs has not yet been ported to CMake.
```

- Use `-D__PW_FPGA` to enable FPGA support for PW (fft) calculations. Currently tested only for
  Intel Stratix 10 and Arria 10 GX1150 FPGAs.
- Supports single precision and double precision fft calculations with the use of dedicated APIs.
- Double precision is the default API chosen when set using the `-D__PW_FPGA` flag.
- Single precision can be set using an additional `-D__PW_FPGA_SP` flag along with the `-D__PW_FPGA`
  flag.
- Kernel code must be synthesized separately and copied to a specific location.
- See <https://github.com/pc2/fft3d-fpga/> for the kernel code and instructions for synthesis.
- Read `src/pw/fpga/README.md` for information on the specific location to copy the binaries to.
- Currently supported FFT3d sizes - 16^3, 32^3, 64^3.
- Include aocl compile flags and `-D__PW_FPGA -D__PW_FPGA_SP` to `CFLAGS`, aocl linker flags to
  `LDFLAGS` and aocl libs to `LIBS`.
- When building FPGA and OFFLOAD together then `-D__NO_OFFLOAD_PW` must be used.
