---
project: Simple DFT-D3
summary: A simple reimplementation of the DFT-D3 dispersion model
project_github: https://github.com/dftd3/simple-dftd3
project_download: https://github.com/dftd3/simple-dftd3/releases
project_website: https://dftd3.readthedocs.io
author: Sebastian Ehlert
author_description: Quantum chemistry researcher developing semi-empirical quantum chemistry in Fortran and Python.
github: https://github.com/awvwgk
src_dir: ./src
         ./app
output_dir: ./_docs
exclude_dir: ./test
docmark: <
predocmark: >
source: true
graph: false
sort: alpha
preprocess: false
print_creation_date: true
extra_mods: iso_fortran_env:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html
            mctc_io:https://grimme-lab.github.io/mctc-lib/modules/mctc_io.html
            mctc_env:https://grimme-lab.github.io/mctc-lib/modules/mctc_env.html
creation_date: %Y-%m-%d %H:%M %z
md_extensions: markdown.extensions.toc
               markdown.extensions.smarty
---

A simple drop-in replacement for ``dftd3``.

This program provides a small and easy to use implementation of the DFT-D3 dispersion correction
(see [*JCP* **132**, 154104 (2010)](https://dx.doi.org/10.1063/1.3382344)
and [*JCC* **32**, 1456 (2011)](https://dx.doi.org/10.1002/jcc.21759) for details).

It is mostly based on the [`dftd4`](https://github.com/dftd4/dftd4) program and
borrows one or two ideas from the implementation in [`ased3`](https://github.com/ehermes/ased3).


[TOC]


## Usage

To use DFT-D3 in your application you can either use the Fortran API or the C API.

### Fortran API

To perform a D3(BJ)-ATM calculation in a Fortran program, the [[dftd3]] module should be imported.
Other data types to communicate with the library are obtained from the [MCTC library](https://grimme-lab.github.io/mctc-lib) using the environment module ([[mctc_env]]) and the IO module ([[mctc_io]]).

```fortran
subroutine calc_dftd3(mol, method, energy, gradient, sigma, error)
   use mctc_env
   use mctc_io
   use dftd3
   type(structure_type), intent(in) :: mol
   character(len=*), intent(in) :: method
   real(wp), intent(out) :: energy
   real(wp), intent(out) :: gradient(:, :)
   real(wp), intent(out) :: sigma(:, :)
   type(error_type), allocatable, intent(out) :: error
   type(d3_model) :: disp
   type(d3_param) :: inp
   type(rational_damping_param), allocatable :: param

   call get_rational_damping(inp, method, error, s9=1.0_wp)
   if (allocated(error)) return
   call new_rational_damping(param, inp)

   call new_d3_model(disp, mol)

   call get_dispersion(mol, disp, param, realspace_cutoff(), energy, &
      & gradient, sigma)

end subroutine calc_dftd3
```


### C API

An example wrapper for a DFT-D3(BJ)-ATM calculation is given here.

```c
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "dftd3.h"

static const buffersize = 512;

int
calc_dftd3(int natoms, int* numbers, double* positions,
           double* lattice, bool* periodic, char* method,
           double* energy, double* gradient, double* sigma)
{
   // Local API objects from the s-dftd3 library
   dftd3_error error = dftd3_new_error();
   dftd3_structure mol = NULL;
   dftd3_model disp = NULL;
   dftd3_param param = NULL;
   int stat = EXIT_SUCCESS;

   // Create a new geometry for the library to work with
   mol = dftd3_new_structure(error, natoms, numbers, positions, lattice, periodic);
   stat = dftd3_check_error(error);

   if (stat) {
      // Initialize the D3 dispersion model for the given structure,
      // this step depends on the atomic numbers, but not on the actual geometry
      disp = dftd3_new_d3_model(error, mol);
      stat = dftd3_check_error(error);
   }

   if (stat) {
      // Load D3(BJ)-ATM parameters for the given method from internal storage,
      // this step depends on the atomic numbers, but not on the actual geometry
      param = dftd3_load_rational_damping(error, mol, method, true);
      stat = dftd3_check_error(error);
   }

   if (stat) {
      // Evaluate the dispersion energy, gradient and virial,
      // the gradient and virial are optional and can be replaced by NULL
      dftd3_get_dispersion(error, mol, disp, param, &energy, gradient, sigma);
      stat = dftd3_check_error(error);
   }

   if (!stat) {
      char buffer[buffersize];
      dftd3_get_error(error, buffer, buffersize);
      printf("[Error] %s\n", buffer);
   }

   // Always free the used memory
   dftd3_delete_error(&error);
   dftd3_delete_structure(&mol);
   dftd3_delete_model(&disp);
   dftd3_delete_param(&param);

   return stat;
}
```

Overall, any DFT-D3 calculation requires the creation of four API objects.

The error handling is done with a ``dftd3_error`` handle, it can be checked using the ``dftd3_error_check`` function and returns zero values on success and non-zero values on failures.
On failures the error handle can be queried for the error message with ``dftd3_get_error``.
The library itself will not attempt to write or terminate on any encounter of an error.

To pass the geometry information a ``dftd3_structure`` object is used.
The object has immutable number of atoms, atomic numbers and boundary conditions, but allows updating the coordinates and lattice vectors after creation.
The structure object is required to initialize system specific data of other API objects.
Atomic numbers are allowed in a range of 1 to 118, but the used dispersion model might support a smaller range.

The actual D3 dispersion model is created as ``dftd3_model`` object using the ``new_d3_model`` constructor.
It contains the interpolation scheme for the dispersion coefficients, which is common to all DFT-D3 dispersion corrections independently of the damping function.
In principle this part can be replaced by a different scheme to evaluate the dispersion coefficients.
This step might require an internal initialization on the first invocation.
The internal initialization is done once and in an OpenMP thread-safe way in case the library is compiled OpenMP threading.
Otherwise guard this API call with a ``critical`` pragma if you require thread-safe execution.

Finally, to evaluate dispersion energies and gradients a damping function is required to connect the dispersion model to a method like a density functional.
The damping parameters are stored as ``dftd3_param`` object and can either be loaded from the internal storage of the library or created by supplying the parameters directly.

Make sure to delete the API objects after you are done with them.
The deconstructor will return without any action in case a null pointer is passed.
After deconstructing the API objects are overwritten with a null pointer.


## Getting Started

### Meson

Create a new meson project and include `s-dftd3` either as git-submodule in your subprojects directory or create a wrap file to fetch it from upstream:

```ini
[wrap-git]
directory = s-dftd3
url = https://github.com/dftd3/simple-dftd3
revision = head
```

To load the project the necessary boilerplate code for subprojects is just

<!--pygments doesn't know about meson, python highlighting looks okayish-->
```python
sdftd3_prj = subproject(
  's-dftd3',
  version: '>=0.1',
  default_options: [
    'default_library=static',
  ],
)
sdftd3_dep = sdftd3_prj.get_variable('sdftd3_dep')
```

Now you can add `sdftd3_dep` to your dependencies and access the public API by the `dftd3` module.

We recommend to set the default library type of `s-dftd3` to static when linking your applications or library against it.
Note for library type both and shared `s-dftd3` will install itself along with your project.

For more fine-tuned control you can access:

- the library target with `sdftd3_lib`
- the private include dir of this target, containing the Fortran module files, with `sdftd3_inc`
- the license files of `s-dftd3` with `sdftd3_lic`

If you are linking your application statically against `s-dftd3` and still want to distribute the license files of `s-dftd3` (thank you), just use

```python
install_data(
  sdftd3_prj.get_variable('sdftd3_lic'),
  install_dir: get_option('datadir')/'licenses'/meson.project_name()/'s-dftd3',
)
```


### Fortran Package Manager (fpm)

This project supports [fpm](https://github.com/fortran-lang/fpm) as build system as well.
Just add it to the dependencies in your `fpm.toml` file:

```toml
[dependencies]
s-dftd3.git = "https://github.com/dftd3/simple-dftd3"
```
