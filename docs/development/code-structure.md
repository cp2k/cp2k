# Code Structure

CP2K is a large, complex application which has many features, methods and algorithms implemented.
When looking at the code for the first time it can be very challenging to understand how it all
works, or even where to start looking! This page is intended for novice developers who have read and
understood the literature and wish to locate the relevant algorithms and data structures in the
code.

## File Names

Source files should have prefixes corresponding to their main functionalities. For example:

- `qs_*` for Quickstep related source codes (Hamiltonian construction, integration, collocation,
  energy minimisation and SCF cycle etc)
- `xc_*` for Exchange-Correlation functionals used by Quickstep
- `md_*` for Molecular Dynamics related source codes
- `mc_*` for Quantum Monte Carlo related source codes
- `fist_*` for FIST classical MD related source codes
- `input_*` for general input functions of CP2K
- `qmmm_*` for QM/MM related source codes
- `message_*` for MPI message passing related source codes
- `machine_*` for architecture dependent codes
- `admm_*` for auxilliary density matrix (ADMM) method related codes
- `ai_*` for integrals of the primitive cartesian Gaussians
- `atomic_*` for datatypes related to information on atoms in a simulation
- `atom_*` for atomic calculations

These prefixes are not exclusive, nor are they always logical. There are exceptions in code naming
conventions, for example:

- `realspace_grid_types.F` and `realspace_grid_cube.F` are both used in Quickstep, but do not have
  the corresponding `qs_` prefix

## Overall Structure

- Extensive use of Fortran modules, and there are **no** global variables
- Major parts of the CP2K code are compiled into separate libraries, for example:
  - `libcp2kmain`, `libcp2kbase`, `libcp2kdbcsrwrap`, `libcp2kfft` etc.

## Structure of Quickstep

Quickstep part of the CP2K code calculates the ab initio self-consistent Kohn-Sham energy and the
associated forces of a periodic system. The calculation involves

- Construction of the Kohn-Sham energy functional and Hamiltonian, which involves:
  - Mapping of operators represented as matrices in Gaussian basis onto the real space (RS)
    multi-grids (collocation). This is required for the computation of the Hartree potential, which
    is calculated in the planewave basis, and the exchange-correlation energy density functional
  - Mapping of functions defined on the RS grids into matrix elements represented in the Gaussian
    basis (integration)
  - Fast Fourier Transform that maps functions defined on each level of the RS multi-grid into the
    corresponding planewave coefficients; and its reverse operation
- Dense and Sparse linear algebra operations, e.g. matrix multiplications
- Minimisation of the Kohn-Sham energy with respect to the electronic density matrix (using matrix
  operations)
- Self-consistent cycle for the electronic charge density

Most of the computational time is spent on:

- Collocation
- Integration
- Linear algebraic operations
- Fast Fourier Transforms

## Data Structure of Key Variables

This subsection is on the modules containing the definition of the key data used in Quickstep
calculations

- Electronic density and its derivatives, in various representations: sparse matrix in Gaussian
  basis, function on RS multi-grid, and as planewave coefficients etc. All density data is contained
  in a single container derived type.
  - Module: `qs_rho_types`
  - File: `qs_rho_types.F`
  - Container type: `qs_rho_type`
