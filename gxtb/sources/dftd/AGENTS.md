# AGENTS.md - Development Guide for DFT-D4

This document provides guidance for AI agents and contributors working with the DFT-D4 codebase.

## Project Overview

**DFT-D4** (Generally Applicable Atomic-Charge Dependent London Dispersion Correction) is a Fortran-based library and CLI tool implementing the D4 dispersion model. The project delivers:
- Dispersion energies, gradients, and Hessian via numerical differentiation.
- Pairwise-resolved analysis and support for 3D periodic systems.
- Public APIs: Fortran modules, a C interface, and an optional Python CFFI extension with interfaces to ASE, PySCF and QCSchema.

## Repository Structure

```
dftd4/
├── app/               # CLI frontend (dftd4) and sample input files
├── assets/            # Parameter tables and example data
├── config/            # Build configuration scripts
├── doc/               # Sphinx documentation; ford.md for developer docs
├── include/           # Public C headers (dftd4.h)
├── man/               # Manual pages (asciidoc)
├── python/            # Optional Python extension (CFFI-based)
├── src/dftd4/         # Main library source code
│   ├── damping/       # Damping models and parameters
│   ├── data/          # Element and model data
│   ├── integrator/    # Numerical integration
│   ├── model/         # Dispersion model implementations
|   │   ├── reference/ # Reference polarizabilities
│   │   └── ...        # Definition of the dispersion models
│   ├── param/         # Functional specific damping parameters
│   ├── api.f90        # C API implementation (optional)
│   └── ...            # Core modules: cutoff, numdiff, output, utils, version
└── subprojects/       # Meson wrap dependencies (mctc-lib, mstore, multicharge, json-fortran)
├── test/              # Tests (unit/ for Fortran, api/ for C examples)
```

## Build Systems

DFT-D4 supports Meson (preferred), CMake, and fpm.

### Meson (Recommended)

```bash
# Configure
meson setup _build

# Build
meson compile -C _build

# Test
meson test -C _build --print-errorlogs

# Install
meson install -C _build
```

**Notes**:
- Requires Fortran 2008 compiler, Ninja backend, and a BLAS/LAPACK provider.
- GCC/Intel are tested; GCC 15.0.x–15.1.x have a known bug that can trigger interface mismatch errors.
- Meson version 1.8.0 has a known bug and is explicitly unsupported. Use any other version ≥ 0.55.
- Select compiler via `FC=... meson setup ...`.

### CMake

```bash
# Configure
cmake -B _build -G Ninja -DCMAKE_INSTALL_PREFIX=$HOME/.local

# Build
cmake --build _build

# Test
ctest --test-dir _build --output-on-failure

# Install
cmake --install _build
```

**Note**: The CMake build does not produce the Python extension module.

### fpm (Fortran Package Manager)

```bash
# Build
fpm build

# Test
fpm test

# Run application
fpm run -- --help
```

**Note**: The fpm build targets the standalone binary and Fortran library; it does not export the C API.

## Dependencies

- **BLAS/LAPACK** backend (MKL, OpenBLAS, Netlib, or custom).
- **mctc-lib** – structure handling and error management.
- **mstore** – data storage utilities.
- **multicharge** – atomic charge model.
- **json-fortran** – JSON parsing (used for JSON output and tests).
- **OpenMP** – enabled by default; disable with build options if needed.
- Optional: asciidoctor (man pages), FORD (developer docs), Python ≥ 3.6 with `cffi` (Python extension), C compiler (C API tests).

Dependencies are provided via Meson subprojects/wraps, CMake find modules, fpm, or external installs.

## APIs and Interaction Modes

### Fortran API (preferred)
- Consume as a Meson subproject (`dependency('dftd4')`) or link an installed library via pkg-config.
- With fpm, add to your manifest:
  ```toml
  [dependencies]
  dftd4.git = "https://github.com/dftd4/dftd4"
  ```
  then `use dftd4_*` modules (e.g., `dftd4_model`, `dftd4_damping`).

#### Core Data Types
- `structure_type` (`mctc_io`) - molecular structure container (numbers, positions, lattice, periodic flags).
- `dispersion_model` (`dftd4_model_type`) - abstract base for D4/D4S dispersion models.
- `d4_model` / `d4s_model` (`dftd4_model_d4`, `dftd4_model_d4s`) - concrete dispersion model instances.
- `damping_type` (`dftd4_damping_type`) - damping function for two- and three-body as well as parameters.
- `param_type` (`dftd4_param_type`) - damping parameters.
- `realspace_cutoff` (`dftd4_cutoff`) - cutoff settings for CN, two-body, and three-body terms.
- `error_type` (`mctc_env`) - error propagation across Fortran APIs.

#### Core Functions
- `read_structure` (`mctc_io`) - load a structure from disk into `structure_type`.
- `new_d4_model` / `new_d4s_model` (`dftd4_model_d4`, `dftd4_model_d4s`) - construct dispersion models.
- `get_damping_params` (`dftd4_param`) - fetch damping parameters by functional name or ID.
- `get_dispersion` (`dftd4_disp`) - compute dispersion energy (and optionally gradients/virial).
- `get_pairwise_dispersion` (`dftd4_disp`) - compute pairwise-resolved energies.
- `get_properties` (`dftd4_disp`) - compute CN, charges, C6, and polarizabilities.

### C API
- Enable with `-Dapi=true`; interfaces live in `src/dftd4/api.f90`, header in `include/dftd4.h`.
- Include the header and link with pkg-config:
  ```bash
  pkg-config --cflags --libs dftd4
  ```
  If pkg-config omits transitive deps, append `-lmulticharge -lmctc-lib -lmstore`.

### Python API
- Install from conda-forge (`conda install dftd4-python`), or build in-tree:
  ```bash
  meson setup _build -Dpython=true -Dpython_version=$(which python3)
  meson compile -C _build
  ```
  Out-of-tree build instructions live in `python/README.rst`.

### CLI (for scripting/verification)
- Build and run the installed or in-tree binary:
  ```bash
  meson compile -C _build
  _build/app/dftd4 --func pbe0 coord
  ```
  Add `--json` to emit machine-readable output for downstream tools.

## Testing

### Running Tests

Tests are organized by functionality in `test/`:
- `test/unit/test_model.f90` - dispersion model construction and core behavior.
- `test/unit/test_damping.f90` - damping function evaluation.
- `test/unit/test_dftd4.f90` - end-to-end D4 energy/gradient computations.
- `test/unit/test_pairwise.f90` - pairwise-resolved dispersion energies.
- `test/unit/test_param.f90` - parameter parsing and access from `assets/parameters.toml`.
- `test/unit/test_periodic.f90` - periodic boundary handling and cell-related paths.
- `test/api/example.c` - build/link sanity check for the C interface.

Commands:
```bash
meson test -C _build --print-errorlogs
ctest --test-dir _build --output-on-failure
fpm test
```

### Python API Tests

The Python tests live in `python/dftd4/test_*.py` and use `pytest`.

Recommended setup (editable install into the venv + test; no optional deps):
```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install meson meson-python ninja cffi numpy pkgconfig pytest

meson setup --wipe _build \
  -Dprefix="$VIRTUAL_ENV/_install" \
  -Dlibdir=lib \
  -Ddefault_library=shared \
  -Dapi=true \
  -Dpython=true \
  -Dpython_version="$(command -v python3)" \
  -Dbuildtype=debug
meson compile -C _build
meson install -C _build

PKG_CONFIG_PATH="$VIRTUAL_ENV/_install/lib/pkgconfig:$PKG_CONFIG_PATH" \
  pkg-config --modversion dftd4
PKG_CONFIG_PATH="$VIRTUAL_ENV/_install/lib/pkgconfig:$PKG_CONFIG_PATH" \
  python -m pip install -e python -Csetup-args="-Dbuildtype=debug"

LD_LIBRARY_PATH="$VIRTUAL_ENV/_install/lib:$LD_LIBRARY_PATH" \
  python -m pytest -vv --pyargs dftd4
```

Optional test modules require extra dependencies:
- `test_ase.py` → `ase`
- `test_qcschema.py` → `qcelemental`
- `test_pyscf.py` → `pyscf`

Out-of-tree build (Python-only):
```bash
cd python
meson setup _build -Dpython_version=$(which python3)
meson compile -C _build
```
Ensure `PKG_CONFIG_PATH` can find the installed `dftd4` library before building.

### Writing Tests

Fortran tests use `mctc_env_testing`:
```fortran
use mctc_env_testing, only : new_unittest, unittest_type, error_type, check

subroutine collect_my_tests(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)
   testsuite = [ &
      & new_unittest("test-name", test_procedure), &
      & new_unittest("expected-fail", test_fail, should_fail=.true.) &
   ]
end subroutine
```

Add new tests alongside feature changes and register them in the relevant suite.

## Code Style and Conventions

- Fortran 2008, free-form (`.f90`, `.F90` for preprocessed).
- Module naming: `dftd4_<subsystem>_<component>` (e.g., `dftd4_model_d4`, `dftd4_damping`).
- Default `private`; explicitly export public entities.
- Always `implicit none` in all program units.
- Error handling via `mctc_env:error_type` and `fatal_error` for propagation; avoid bare `error stop` except in test drivers.
- Maintain consistent formatting and clear intent; add brief comments where code is non-obvious.

### File Header

All source files must include the LGPL header and SPDX tag:

```fortran
! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! Lesser GNU General Public License for more details.
!
! You should have received a copy of the Lesser GNU General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.
```

### Error Handling

Use the `error_type` for error propagation:

```fortran
use mctc_env, only : error_type, fatal_error

subroutine my_routine(result, error)
   type(error_type), allocatable, intent(out) :: error
   
   if (some_error_condition) then
      call fatal_error(error, "Descriptive error message")
      return
   end if
end subroutine
```

### Documentation

- FORD docstrings for public interfaces (`!>` preceding, `!<` trailing).
- Sphinx docs in `doc/`; man pages in `man/` (asciidoc).
- Keep CLI help (`app/help.f90`), man pages, and README examples synchronized when changing behavior.

## CI/CD Workflow

GitHub Actions workflows are configured in `.github/workflows/`:
- Build systems: Meson, CMake, and fpm.
- Compilers: GCC and Intel.
- Documentation builds (Sphinx/FORD where applicable).
- Coverage upload to Codecov.
- DCO enforcement via `.github/dco.yml` (all commits must be signed off).

## Adding New Features

### Adding or Updating Parameters/Functionals
1. Update `assets/parameters.toml` with new entries and references.
2. Adjust parameter handling in `src/dftd4/param.f90` and related data modules if needed.
3. Add regression/unit tests in `test/unit/`.
4. Update documentation (`doc/`, man pages, README snippets) to reflect new options or defaults.

### Extending APIs
1. Fortran: add modules under `src/dftd4/` and export needed interfaces.
2. C API: extend `src/dftd4/api.f90` and `include/dftd4.h`.
3. Python: update CFFI bindings in `python/` and adjust `pyproject.toml`/`meson_options.txt` as needed.
4. Register new sources in `meson.build`/`CMakeLists.txt`.
5. Add tests (Fortran or C) and documentation.

### CLI Enhancements
1. Implement behavior in `app/` (e.g., `argument.f90`, `driver.f90`, `help.f90`).
2. Sync help text, man page (`man/dftd4.1.adoc`), and README/recipes.
3. Add examples to `app/` if helpful.
4. Cover with regression tests when feasible.

## Common Tasks

### Run a Dispersion Calculation
```bash
meson compile -C _build
_build/app/dftd4 --func pbe0 coord
```

### Dump Results to JSON (suppress .EDISP)
```bash
_build/app/dftd4 --func pbe0 --json --grad --noedisp struct.xyz
```

### Pairwise-Resolved Energies
```bash
_build/app/dftd4 --pair-resolved mol.xyz
```

### Fortran API Example (compute D4 energy)
```fortran
program d4_energy
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type, read_structure
   use dftd4, only : d4_model, new_d4_model, damping_type, &
      & get_damping_params, realspace_cutoff, get_dispersion, &
      & dftd_models, twobody_damping_function, threebody_damping_function, &
      & param_type
   implicit none

   type(structure_type) :: mol
   type(d4_model) :: model
   type(damping_type) :: damp
   type(param_type), allocatable :: param
   type(realspace_cutoff) :: cutoff
   type(error_type), allocatable :: error
   real(wp) :: energy

   call read_structure(mol, "coord", error)
   if (allocated(error)) then
      write (*, '(a)') trim(error%message)
      stop 1
   end if

   call new_d4_model(error, model, mol)
   if (allocated(error)) then
      write (*, '(a)') trim(error%message)
      stop 1
   end if

   call get_damping_params(error, "pbe", dftd_models%d4, d4%default_damping_2b, &
      & d4%default_damping_3b, param)
   call new_damping(error, damp, d4%default_damping_2b, d4%default_damping_3b)
   call get_dispersion(mol, model, damp, param, cutoff, energy)

   write (*, '(a,f18.12)') "D4 dispersion energy (Hartree): ", energy
end program d4_energy
```

### List Available Parameters
```bash
_build/app/dftd4 param --list
# or inspect assets/parameters.toml
```

### Linking into Other Codes (e.g., VASP)
- Use pkg-config: `pkg-config --cflags --libs dftd4`
- If transitive deps are missed, add `-lmulticharge -lmctc-lib -lmstore` explicitly.

## Build Options

### Meson Options (`meson_options.txt`)

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `lapack` | combo | `auto` | BLAS/LAPACK backend (`auto`, `mkl`, `mkl-rt`, `openblas`, `netlib`, `custom`) |
| `custom_libraries` | array | `[]` | Extra libraries for custom BLAS/LAPACK |
| `openmp` | boolean | `true` | Enable OpenMP parallelization |
| `api` | boolean | `true` | Build C API |
| `python` | boolean | `false` | Build Python extension module |
| `python_version` | string | `python3` | Python executable to link against |
| `ilp64` | boolean | `false` | Enable BLAS/LAPACK ILP64 |

### CMake Options (selected)

| Option | Description |
|--------|-------------|
| `WITH_OpenMP` | Enable OpenMP |
| `WITH_API` | Build C API |
| `WITH_PYTHON` | Build Python extension (may require out-of-tree steps) |

## Troubleshooting

1. **GCC 15.0.x–15.1.x interface mismatch bug**: use another GCC version.
2. **Meson 1.8.0 error**: upgrade or downgrade meson (`pip install meson!=1.8.0`).
3. **JSON tests fail** (mctc-lib JSON I/O): ensure the jonquil dependency is available, or disable JSON in the mctc-lib subproject with Meson (`-Djson=disabled`).
4. **Missing BLAS/LAPACK**: set `lapack` option (Meson) or pass provider paths; ensure `LD_LIBRARY_PATH`/`LIBRARY_PATH` includes the backend.
5. **Pkg-config misses deps**: link `-lmulticharge -lmctc-lib -lmstore` explicitly.
6. **OpenMP issues**: disable with `-Dopenmp=false` if encountering runtime problems.
7. **Python build failures**: ensure `python_version` points to a Python 3 with `cffi` installed.

## External Resources

- [GitHub Repository](https://github.com/dftd4/dftd4)
- [Documentation](https://dftd4.readthedocs.io/en/latest/)
- [API Docs](https://dftd4.github.io/dftd4/)
- [Issue Tracker](https://github.com/dftd4/dftd4/issues)

## Contributing

1. Fork the repository and create a feature branch.
2. Follow code style and licensing conventions (LGPL header, SPDX).
3. Add tests for new functionality and keep docs in sync.
4. Ensure all tests pass on the supported build systems.
5. Sign off commits (DCO requirement: `git commit -s -m "message"`).
6. Submit a pull request with a clear description and rationale.

Error messages are designed to be helpful; report unclear errors as bugs with reproduction steps.
