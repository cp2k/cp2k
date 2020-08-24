# FPGA

## Finding the kernel code

- See <https://github.com/pc2/fft3d-fpga> for the kernel code and instructions
  on synthesizing them for INTEL FPGAs.

### Copy the bitstream generated to a particular path

- The bitstream generated on synthesis can be found as the `.aocx` file.
- Copy this file to the following locations:
  - If single precision: `~/cp2k/fpgabitstream/fft3d/synthesis_sp/syn??/*.`
  - If double precision: `~/cp2k/fpgabitstream/fft3d/synthesis_dp/syn??/*.`
  - where `syn??` is dependent on the size of the FFT3d
  - Therefore, a 16^3 FFT3d file should be copied to `syn16` folder.
  - The sizes supported are 16^3, 32^3, 64^3.
- If the required FFT size is different from the default options:
  - modify the switch case in the `fft_fpga.c` to include the required size
    and the path to the location of the bitstream.
