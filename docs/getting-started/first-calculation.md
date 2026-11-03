# Run a First Calculation

This page walks through a small single-point energy calculation for a water molecule. It uses the
{term}`Quickstep` module for the density functional theory ({term}`DFT`), the Gaussian and plane
wave ({term}`GPW`) method, a molecular Gaussian basis set, and Goedecker-Teter-Hutter ({term}`GTH`)
pseudopotentials. The example is intentionally small enough to run in a few seconds while still
showing the parts of a typical CP2K input file that matter for larger calculations.

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

The following assumes a successful installation of CP2K and a proper environment configuration. Run
the calculation with the installed CP2K binary as:

```shell
OMP_NUM_THREADS=1 cp2k.psmp -i h2o.inp -o h2o.out
```

For MPI-parallel runs, launch CP2K through the MPI launcher used on your system, for example:

```shell
mpirun -np 2 -x OMP_NUM_THREADS=1 cp2k.psmp -i h2o.inp -o h2o.out
```

The executable name depends on how CP2K was built, with the suffix one of the four possibilities of
`(p|s)(smp|dbg)` (`p` for parallel, `s` for serial, `smp` for symmetric multiprocessing and `dbg`
for debug):

| executable  |           meaning           |
| :---------: | :-------------------------: |
| `cp2k.psmp` | MPI + OpenMP parallel build |
| `cp2k.ssmp` |     serial/OpenMP build     |
| `cp2k.pdbg` |  MPI + OpenMP debug build   |
| `cp2k.sdbg` |  serial/OpenMP debug build  |

`cp2k.psmp` supports both MPI and OpenMP; in the `mpirun` command above, the `-np 2` argument starts
2 MPI processes and the `OMP_NUM_THREADS=1` declaration (either as an environment variable or as the
`-x` argument of `mpirun`) requests 1 OpenMP thread per process, thus requiring a total of 2x1=2
logical cores. (How it maps or binds to physical cores depends on the MPI library and CPU hardware.)

Setting `OMP_NUM_THREADS=1` keeps this first example in a simple MPI-only layout. In practice, the
number of MPI processes and OpenMP threads may be adjusted based on the CPU and RAM allocated for a
CP2K task against its computational cost. Note that not all functionalities have been developed with
support or optimization for OpenMP parallelization yet; if sufficient memory is available, the more
preferable way is to run a task mostly or fully with MPI processes for parallelization.

```{tip}
In addition to the suffixes above, a symbolic link is also available for the first two in the form
of `cp2k.popt -> cp2k.psmp` or `cp2k.sopt -> cp2k.ssmp`. If CP2K is launched with `cp2k.popt` or
`cp2k.sopt`, the OpenMP parallel execution will use 1 thread per process, equivalent to having set
`OMP_NUM_THREADS=1` for `cp2k.psmp` or `cp2k.ssmp`.
```

The `-i` option of `cp2k.psmp` provides the input file name and `-o` provides the output file name.
To see the output printed to screen while also saving it into a file, replace the `-o h2o.out` part
with `| tee h2o.out` (i.e. piping the output to the linux `tee` command).

## What the Input Does

The input of CP2K is well-organized using a hierarchical structure. The basic elements of the syntax
are **SECTION** and **KEYWORD** in layers. Each **SECTION** opens with `&SECTION_NAME` and ends with
`&END SECTION_NAME`, both starting with a `&` character to distinguish with **KEYWORD**. Indentation
is not required, though a unit of two whitespaces for each is commonly used for clarity. Each of the
sections and keywords accept certain type and number of input values. Defining a section explicitly
enables the corresponding functionality; if keywords under an explicit section are not explicitly
given, they will implicitly use default values.

```{note}
Most keyword use the syntax `KEYWORD_NAME VALUE`, but there are two special cases:

- The `SECTION_PARAMETERS` of a section needs a value to be specified on the line that the section
  begins, that is, right after `&SECTION_NAME`.
- The `DEFAULT_KEYWORD` in a section needs the value(s) to be directly specified in the indicated
  format, that is, not preceded by any `KEYWORD_NAME`.

For a logical-type keyword:

- Any of `0`, `F`, `.F.`, `FALSE`, `.FALSE.`, `N`, `NO`, `OFF` will be parsed as false.
- Any of `1`, `T`, `.T.`, `TRUE`, `.TRUE.`, `Y`, `YES`, `ON` will be parsed as true.
- The absence of any of these values, that is when a logical keyword is written as an explicit "lone
  keyword", is by default to be parsed as true.

For a real-type keyword in need of a float-point numerical value:

- Both the usual notation of decimal expansion 0.001 and the scientific notation 1E-3 can be used.
- Values will be treated in double precision.
- Simple products and fractions combined with certain functions of a single number can be used like
  `2/3`, `0.3*COS(60)` or `-SQRT(3)/2`. The Fortran intrinsic math functions `COS`, `EXP`, `LOG`,
  `LOG10`, `SIN`, `SQRT`, and `TAN` are available in uppercase.
- Some keywords allow for a unit specification in a square bracket between the keyword name and the
  numerical value; internal unit conversion to *atomic units* happens right during parsing before
  any subsequent calculation, and valid units are documented in the dedicated Units page.
```

The `h2o.inp` above starts with the [GLOBAL](#CP2K_INPUT.GLOBAL) section.

- [PROJECT](#CP2K_INPUT.GLOBAL.PROJECT_NAME) specifies the name of the project.
- [RUN_TYPE](#CP2K_INPUT.GLOBAL.RUN_TYPE) is set to `ENERGY`, so CP2K evaluates just a single-point
  energy; force and stress will not be evaluated and atoms will not undergo any motion.
- [PRINT_LEVEL](#CP2K_INPUT.GLOBAL.PRINT_LEVEL) is not specified explicitly, so the default `MEDIUM`
  will be used to control the amount of output: information of `SILENT`, `LOW` and `MEDIUM` level
  will be printed whereas that of `HIGH` and `DEBUG` level will not.

Subsequently, the [FORCE_EVAL](#CP2K_INPUT.FORCE_EVAL) section sets up the theoretical method and
system to evaluate the energy.

- The [METHOD](#CP2K_INPUT.FORCE_EVAL.METHOD) keyword selects `Quickstep`, the electronic-structure
  module of CP2K for density functional theory and related methods.
- The [DFT](#CP2K_INPUT.FORCE_EVAL.DFT) section contains specification for said method.
  - [BASIS_SET_FILE_NAME](#CP2K_INPUT.FORCE_EVAL.DFT.BASIS_SET_FILE_NAME) and
    [POTENTIAL_FILE_NAME](#CP2K_INPUT.FORCE_EVAL.DFT.POTENTIAL_FILE_NAME) tell CP2K where to find
    the basis-set and pseudopotential libraries; they are assumed to be in the directory that the
    CMake build option `CP2K_DATA_DIR` points to.
  - As the water molecule has a spin multiplicity of 1 and no net charge in its closed-shell ground
    state, the [CHARGE](#CP2K_INPUT.FORCE_EVAL.DFT.CHARGE) uses the default of 0 implicitly,
    [MULTIPLICITY](#CP2K_INPUT.FORCE_EVAL.DFT.MULTIPLICITY) uses the default of 1 (for an even total
    number of electrons) implicitly, and there is no need to request for an unrestricted
    (spin-polarized) calculation with [UKS](#CP2K_INPUT.FORCE_EVAL.DFT.UKS).
  - The [QS](#CP2K_INPUT.FORCE_EVAL.DFT.QS) section is implicitly used; the Gaussian and plane wave
    method is set because the default of [METHOD](#CP2K_INPUT.FORCE_EVAL.DFT.QS.METHOD) is `GPW`.
  - The [MGRID](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID) section controls the real-space integration grids
    used by the GPW method, where the [CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.CUTOFF) and
    [REL_CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.REL_CUTOFF) are explicitly given.
  - The [POISSON](#CP2K_INPUT.FORCE_EVAL.DFT.POISSON) section selects a Poisson solver for the
    electrostatic boundary condition; as an isolated water molecule is not periodic, the
    [PERIODIC](#CP2K_INPUT.FORCE_EVAL.DFT.POISSON.PERIODIC) is set as `NONE` and the
    [PSOLVER](#CP2K_INPUT.FORCE_EVAL.DFT.POISSON.POISSON_SOLVER) is set to `MT`.
  - The [SCF](#CP2K_INPUT.FORCE_EVAL.DFT.SCF) section controls the iterations of the self-consistent
    field ({term}`SCF`) method, with [EPS_SCF](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.EPS_SCF) for accuracy
    and [MAX_SCF](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.MAX_SCF) for maximum number of iterations.
  - The [XC](#CP2K_INPUT.FORCE_EVAL.DFT.XC) section selects the DFT exchange-correlation functional,
    in this case the GGA functional PBE with no Hartree-Fock exchange or dispersion correction.
- The [SUBSYS](#CP2K_INPUT.FORCE_EVAL.SUBSYS) section defines the system chemistry.
  - The [CELL](#CP2K_INPUT.FORCE_EVAL.SUBSYS.CELL) section sets up the cell vectors; although no
    periodic boundary condition ({term}`PBC`) is in use, a cell is necessary for the definition of
    integration grids. Here [PERIODIC](#CP2K_INPUT.FORCE_EVAL.SUBSYS.CELL.PERIODIC) is set to `NONE`
    which is consistent with the Poisson solver above.
  - The [COORD](#CP2K_INPUT.FORCE_EVAL.SUBSYS.COORD) gives the absolute Cartesian atomic coordinates
    with a default [UNIT](#CP2K_INPUT.FORCE_EVAL.SUBSYS.COORD.UNIT) of angstrom.
  - The [KIND](#CP2K_INPUT.FORCE_EVAL.SUBSYS.KIND) sections choose one Gaussian basis set and one
    GTH pseudopotential for each element.

The minimal example of `h2o.inp` does not demonstrate full possibilities of the syntax. In an input
file, inline or standalone comments can be marked with an initial `#` or `!` character, useful for
annotation purposes; a trailing `\` character at the end of a line marks line continuation, useful
for splitting a long list of values on multiple shorter lines for one keyword.

## Checking the Result

At the end of `h2o.out`, CP2K prints the total energy in Hartree:

```text
ENERGY| Total FORCE_EVAL ( QS ) energy [hartree]
```

You should also see a line stating that the self-consistent field ({term}`SCF`) cycle converged. If
the SCF cycle does not converge, increase `MAX_SCF`, improve the initial guess, or use a more robust
SCF setup.

```{important}
Unless explicitly requested in input or labelled otherwise, quantities appearing in the output file
should be interpreted as in atomic units (`a.u.`). Atomic units are defined for physical quantities;
a list of unit conversion factors and relevant constants used in CP2K can be printed by setting the
[PHYSCON](#CP2K_INPUT.GLOBAL.PRINT.PHYSCON) under the `&GLOBAL/&PRINT` section to true.

In this example, the unit of energy printed on the `ENERGY| Total FORCE_EVAL ( QS ) energy` line is
controlled by the [ENERGY_UNIT](#CP2K_INPUT.FORCE_EVAL.PRINT.PROGRAM_RUN_INFO.ENERGY_UNIT) keyword
under the `&FORCE_EVAL/&PRINT/&PROGRAM_RUN_INFO` section, which defaults to `hartree`, the atomic
unit of energy. Other energy units like `eV` and `kcalmol` may be used instead. This particular
keyword is independent of the method used to evaluate energy, and the label `QS` may be adapted to
methods other than Quickstep DFT.
```

The timing table printed at the end of every CP2K run is useful for a first performance check. For
larger calculations, compare timings between MPI/OpenMP layouts and watch whether most of the time
is spent in grid operations, sparse matrix operations, diagonalization, or communication.

A list of [REFERENCES](#CP2K_INPUT.GLOBAL.REFERENCES) will also be printed at the end for the large
number of relevant works of theoretical foundations and programming developments. Please acknowledge
them by giving proper and citations in the description of publications.

## Next Steps

- If this example calculation has been executed on the cloud, consider building or installing CP2K
  on your server: [](installation)
- Learn the idea behind GPW: [](../methods/dft/gpw)
- Learn about basis sets and pseudopotentials: [](../methods/dft/basis_sets),
  [](../methods/dft/pseudopotentials)
- Resolve or improve SCF convergence: [](../methods/dft/convergence)
- Test for convergence of [CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.CUTOFF) and
  [REL_CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.REL_CUTOFF): [](../methods/dft/cutoff)
- Learn about the periodicities and Poisson solvers: [](../methods/dft/electrostatics/index)
- Proceed to geometry optimization: [](../methods/optimization/geometry_and_cell_opt)
- Proceed to molecular dynamics simulation: [](../methods/sampling/molecular_dynamics)
- Explore more complete examples: <https://github.com/cp2k/cp2k-examples>
- Use heavier tasks to benchmark performance: see `benchmarks` directory

## Video

```{youtube} qMR-NAaUheg
---
url_parameters: ?start=45
align: center
privacy_mode:
---
```
