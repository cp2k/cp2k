# fft: FFT wrapper for CP2K

This package hosts the FFT operations required by cp2k. The code is entirely written in C and can be
built stand-alone in order to be reusable by other projects.

This package offers the following main features:

- Local 1D, 2D and 3D FFTs
- out-of-place and in-place FFTs.
- Local and distributed FFTs using MPI
- Complex-to-complex, Real-to-complex and complex-to-real FFTs
- Fortran API

Currently, this package has its own tiny reference backend (not well-performing) and a FFTW-based
backends. This package currently supports FFTW-based backend by FFTW3 and MKL.

## Unit Test

The `fft_unittest.x` binary runs tests for the different components of this backend (local FFTs, MPI
backends if supported by the library, redistribution routines, compound FFTs, additional features).

```shell
$ cd cp2k/src/fft
$ make
$ ./fft_unittest.x

$ ./grid_unittest.x
Using reference FFT library.
The 1D FFT does work correctly (15 26)!
The 1D FFT does work correctly (18 22)!
The 1D FFT does work correctly (20 28)!
The 1D FFT does work correctly (14 13)!
The 1D R/C FFT does work correctly (15 26)!
...
The R2C-3D FFT with ray layout does work correctly (2 3 5)/(2 4 8)!
The R2C-3D FFT with ray layout does work correctly (8 4 2)/(8 4 2)!
The R2C-3D FFT with ray layout does work correctly (5 3 2)/(8 4 2)!

 The 3D FFT routines work correctly!
Time to test high-level FFTs: 2.032545
The addition between different grids works correctly (2 4 8)/(2 4 8)
The addition between different grids works correctly (8 4 2)/(7 3 2)
The addition between different grids works correctly (2 4 8)/(1 2 4)
The addition between different grids works correctly (11 7 5)/(5 3 2)
 -----------------------------------------------------------------------------------------------------
 -                                                                                                   -
 -                                         FFT TIMING REPORT                                         -
 -                                                                                                   -
 -----------------------------------------------------------------------------------------------------
 ROUTINE                                       CALLS    AVG TOTAL    MAX TOTAL    AVG SELF    MAX SELF
                                                           TIME         TIME         TIME        TIME
 coll_x_dist_yz_r_2_4_5_1                          1        0.000        0.000       0.000       0.000
...
 fft_3d_fw_c2c_local_8_8_8                       104        0.000        0.045       0.000       0.000
All tests have passed :-)
```

## FFT performance test

The `fft_perftest.x` binary runs performance tests for FFT sizes interesting for CP2K runs.
