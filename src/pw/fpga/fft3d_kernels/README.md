# FPGA

## Steps to generate bitstream for the given FFT3d kernel code 

### Set the precision for the FFT3d
 -  The `fft3d_config.h` has a macro called `TYPE_FLOAT` that toggles between
   single and double precision floating point variables.
 - Setting its value to 0 enable double precision floating point and 1 to single
   precision

### Set the size of the FFT3d 
 - The `fft3d_config.h` file contains the parameter `LOGN` to specify the log
    of the size of the FFT.
 - For example, LOGN of 6 represents a 64^3 FFT3d, considering 2^6 = 64
   signifies a dimension of the 3d FFT.
 - Not every kernel size may fit the FPGA device for the given kernel design.

### Synthesize kernels for a particular board
  Example:
	`aoc -fp-relaxed --board $(BOARDNAME) -g -v fft3d.cl -o fft3d.aocx`

### Synthesize kernels for emulation 
  Example:
	`aoc -march=emulator -board=$(BOARDNAME) -g -v fft3d.cl -o fft3d.aocx`


### Copy the bitstream generated to a particular path
 - The bitstream generated can be found as the `.aocx` file. 
 - Copy this file to the following locations:
    - If single precision: `~/cp2k/bitstream/fft3d/synthesis_sp/syn??/*. `
    - If double precision: `~/cp2k/bitstream/fft3d/synthesis_dp/syn??/*.` 
    - where `syn??` is dependent on the size of the FFT3d
    - Therefore, a 16^3 FFT3d file should be copied to `syn16` folder.
    - The sizes supported are 16^3, 32^3, 64^3.
 - If the required FFT size is different from the default options:
    - modify the switch case in the `fft_fpga.cpp` to include the required size
      and the path to the location of the bitstream.


## NOTE
 - Kernel code has been tested to work for uniform FFT3d sizes.

## Setting Up FPGA Fast Emulation Environment for OpenCL 

This can be used to check fpga emulation locally without the availability of
Intel FPGAs or Intel SDK licenses.

### Three steps to installation

 1. Install Intel CPU Runtime for OpenCL Applications
 2. Install Intel SDK for OpenCL Applications
 3. Install AOCL Pro from Intel FPGA SDK for OpenCL
 
### Intel CPU Runtime for OpenCL  Applications

This is the basic driver to use CPUs as OpenCL target device.

Download [Intel CPU Runtime for OpenCL Applications 18.1 for Linux* OS (64bit only)]( https://software.intel.com/en-us/articles/opencl-drivers)  


```
$ tar xzvf l_opencl_p_18.1.0.014.tgz
$ cd l_opencl_p_18.1.0.014
$ sudo ./install.sh

```

Pay attention to the missing packages.

*Note:* You may be required to install `alien` and `libtinfo5` in case `lsb-core` is missing.

### Intel SDK for OpenCL Applications

This installs the necessary tools to compile OpenCL kernel code.

Download [Intel SDK for OpenCL Applications](https://software.intel.com/en-us/intel-opencl) 
and follow the same steps as above. 

### AOCL Pro

This installs the requires device/libraries for `Intel(R) FPGA Emulation Platform for 
OpenCL(TM)` platform.

Download from [link](http://fpgasoftware.intel.com/opencl/) 

```
$ chmod +x  AOCLProSetup-19.1.0.240-linux.run
$ ./AOCLProSetup-19.1.0.240-linux.run
```

You do not need to be sudo to install it, however to write the file in `/opt`, sudo is required.

```

$ sudo sh -c 'echo '<inst dir>/intelFPGA_pro/19.1/hld/linux64/lib/libintelocl_emu.so' >  /etc/OpenCL/vendors/Intel_FPGA_SSG_Emulator.icd'

```

### CLINFO
Install `clinfo` and check if the FPGA emulation platform is recognized. 

If not, check: 

1. the permissions of the `icd` file 
2. if the `libintelocl_emu.so` file has all the shared libraries required ( `ldd libintelocl_emu.so` )
3. if the SDK is able to locate the `icd` file
4. Make sure the missing packages mentioned above are installed

### Actual fast emulation compilation

```
$ ioc64 -device=fpga_fast_emu -input=src/krnl_vadd.cl -ir=kernel.aocx
```

#### If the compilation runs out of resources, check `clinfo` for memory requirements related to the platform and set the following environment variables accordingly

```
CL_CONFIG_CPU_FORCE_GLOBAL_MEM_SIZE
CL_CONFIG_CPU_FORCE_LOCAL_MEM_SIZE 
CL_CONFIG_CPU_FORCE_MAX_MEM_ALLOC_SIZE 
CL_CONFIG_CPU_FORCE_PRIVATE_MEM_SIZE 
CPU_DEV_MAX_WG_PRIVATE_SIZE 
```
