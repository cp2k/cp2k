# Run a First Calculation

This page walks through a small single-point energy calculation for a water molecule. It uses the
{term}`Quickstep` module, the Gaussian and plane wave ({term}`GPW`) method, a molecular Gaussian
basis set, and Goedecker-Teter-Hutter ({term}`GTH`) pseudopotentials. The example is intentionally
small enough to run in a few seconds while still showing the parts of a typical CP2K input file that
matter for larger calculations.

## Input File

Save the following input as `h2o.inp`. The same file is also available as [](h2o.inp).

```text
&GLOBAL
  PROJECT h2o
  RUN_TYPE ENERGY
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT
    POTENTIAL_FILE_NAME GTH_POTENTIALS
    &MGRID
      CUTOFF 400
      REL_CUTOFF 50
    &END MGRID
    &POISSON
      PERIODIC NONE
      PSOLVER MT
    &END POISSON
    &SCF
      EPS_SCF 1.0E-6
      MAX_SCF 50
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC 10.0 10.0 10.0
      PERIODIC NONE
    &END CELL
    &COORD
      O  5.0000  5.0000  5.0000
      H  5.7586  5.0000  5.5043
      H  4.2414  5.0000  5.5043
    &END COORD
    &KIND O
      BASIS_SET DZVP-MOLOPT-GTH
      POTENTIAL GTH-PBE-q6
    &END KIND
    &KIND H
      BASIS_SET DZVP-MOLOPT-GTH
      POTENTIAL GTH-PBE-q1
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
```

## Running CP2K

Run the calculation with one of the installed CP2K binaries:

```bash
OMP_NUM_THREADS=1 cp2k.psmp -i h2o.inp -o h2o.out
```

The executable name depends on how CP2K was built:

| executable  | meaning                     |
| ----------- | --------------------------- |
| `cp2k.psmp` | MPI + OpenMP parallel build |
| `cp2k.pdbg` | MPI + OpenMP debug build    |
| `cp2k.ssmp` | serial/OpenMP build         |
| `cp2k.sdbg` | serial/OpenMP debug build   |

For MPI-parallel runs, launch CP2K through the MPI launcher used on your system, for example:

```bash
mpirun -np 2 -x OMP_NUM_THREADS=1 cp2k.psmp -i h2o.inp -o h2o.out
```

`cp2k.psmp` supports both MPI and OpenMP. Setting `OMP_NUM_THREADS=1` keeps this first example in a
simple MPI-only layout. To see the output on screen while also saving it, replace `-o h2o.out` with
`| tee h2o.out`.

## What the Input Does

[RUN_TYPE](#CP2K_INPUT.GLOBAL.RUN_TYPE) is set to `ENERGY`, so CP2K evaluates the electronic ground
state energy without moving the atoms.

[METHOD](#CP2K_INPUT.FORCE_EVAL.METHOD) selects `Quickstep`, CP2K's electronic-structure module for
Gaussian-based density functional theory and related methods.

`BASIS_SET_FILE_NAME` and `POTENTIAL_FILE_NAME` tell CP2K where to find the basis-set and
pseudopotential libraries. The matching [KIND](#CP2K_INPUT.FORCE_EVAL.SUBSYS.KIND) sections then
choose one Gaussian basis set and one GTH pseudopotential for each element.

[CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.CUTOFF) and
[REL_CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.REL_CUTOFF) control the real-space integration grids
used by the GPW method. They are not a replacement for increasing the Gaussian basis quality; for
accurate work the basis set and grid parameters should be converged together.

The [POISSON](#CP2K_INPUT.FORCE_EVAL.DFT.POISSON) section and the
[CELL](#CP2K_INPUT.FORCE_EVAL.SUBSYS.CELL) section both use `PERIODIC NONE`, which is appropriate
for this isolated molecule in a large non-periodic box.

## Checking the Result

At the end of `h2o.out`, CP2K prints the total energy in Hartree:

```text
ENERGY| Total FORCE_EVAL ( QS ) energy [hartree]
```

You should also see a line stating that the self-consistent field ({term}`SCF`) cycle converged. If
the SCF cycle does not converge, increase `MAX_SCF`, improve the initial guess, or use a more robust
SCF setup.

The timing table printed at the end of every CP2K run is useful for a first performance check. For
larger calculations, compare timings between MPI/OpenMP layouts and watch whether most of the time
is spent in grid operations, sparse matrix operations, diagonalization, or communication.

## Next Steps

- Converge [CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.CUTOFF) and
  [REL_CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.REL_CUTOFF): [](../methods/dft/cutoff)
- Learn the idea behind GPW: [](../methods/dft/gpw)
- Learn about basis sets and pseudopotentials: [](../methods/dft/basis_sets),
  [](../methods/dft/pseudopotentials)
- Build or install CP2K: [](build-from-source), [](build-with-spack), [](distributions)
- Explore more complete examples: <https://github.com/cp2k/cp2k-examples>
- Read the practical CP2K overview paper:
  [The CP2K Program Package Made Simple](https://doi.org/10.1021/acs.jpcb.5c05851)

```{youtube} qMR-NAaUheg
---
url_parameters: ?start=45
align: center
privacy_mode:
---
```
