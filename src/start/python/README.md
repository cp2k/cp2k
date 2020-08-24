# CP2K Python Bindings

## Installation

There is a target `py-cython-bindings` in the global `Makefile` to build the
Python bindings. The shared object can be found in:
`<CP2K_SOURCE_DIR>/lib/<ARCH>/<VERSION>/python`

Only the Python headers and a NumPy installation are required.

## Development

To regenerate the C file from the `cp2k.pyx`, `Cython` is required and should
be called as follows:

```sh
cd <CP2K_SOURCE_DIR>/src/start/python
cython cp2k.pyx
```

Unittests can be found in the `test/` directory. They must be run in separate
Python interpreter instances due to side effects in the library.

## Known Issues

* If libcp2k is built with MPI support, you may get an MPI initialization error
  depending on your MPI implementation/configuration. In that case MPI must be
  initialized first by using Mpi4py and the Fortran MPI communicator handler
  must be passed down the CP2K via the respective `...comm` functions.
  The reason for this is documented here: <https://github.com/jhedev/mpi_python>
