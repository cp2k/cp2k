# fft: FFT wrapper for CP2K

This package hosts the FFT operations required by cp2k. The code is entirely written in C and can be
built stand-alone in order to be reusable by other projects.

This package offers the following main features:

- Local 1D, 2D and 3D FFTs
- out-of-place and in-place FFTs.
- Complex-to-complex, Real-to-complex and complex-to-real FFTs
- Fortran API

Currently, this package has a FFTW-based (including MKL) and a GPU (CUDA or HIP) backend.

## Unit Test

The `fft_unittest.x` binary runs tests for the different components of this backend (local FFTs with
1 to 3 dimensions, complex-to-complex/real-to-complex/complex-to-real, out-of-place/in-place)
