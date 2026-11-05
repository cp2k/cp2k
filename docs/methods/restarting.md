# Restarting CP2K Calculations

Restarting CP2K calculations is useful for continuing interrupted tasks, extending the trajectory of
molecular dynamics simulations, or using converged wavefunctions or geometry as starting points for
new calculations. This guide explains how to produce and utilize two major types of restart files:

- **Motion restart** file for tasks involving atomic motion (such as
  [optimization](./optimization/geometry_and_cell_opt) and
  [molecular dynamics](./sampling/molecular_dynamics)) with any force method;
- **Wavefunction restart** for tasks using the Quickstep module for the force method (such as
  [DFT](./dft/index), [post-HF](./post_hartree_fock/index) and [SE](./semiempiricals/index)) with or
  without atomic motion.

Other specialized functionalities may also define their custom restart files, which are not covered
by this guide. Refer to the dedicated documentation and input section for those cases.

## Prerequisites

### File forms

In Fortran terminology, there are two **forms** of file input/output:

- One is **formatted**, which is for ASCII plain-text, human-readable files;
- The other is **unformatted**, which is for binary machine-readable (but not human readable) files.

The formatted files are open to manual text editing and are portable between machines, while the
unformatted files are generally only intended to be consumed internally by the same program on the
same machine that has produced it. The data storage pattern in unformatted files may be subject to
changes during program development; there is no guarantee for backwards or forwards compatibility of
unformatted files across versions even with identical machine and compiler environments. Contact
developers if there is some advanced need to read and modify unformatted files.

The form of a file can be readily checked with the `file` command: returning `ASCII text` means
formatted, and returning something else like `SIMH tape data` means unformatted.

```shell
$ file ./project-1.restart
/path/to/project-1.restart: ASCII text
$ file ./project-RESTART.wfn
/path/to/project-RESTART.wfn: SIMH tape data
```

### Iteration levels in CP2K

Many methods use an iterative procedure to obtain the results, which may be nested into layers; for
example, when using Quickstep KS-DFT for a geometry optimization task, the optimizer refreshes the
geometry with the energy and gradient (force) information in a step-by-step manner, and each step of
evaluating energy and gradient with KS-DFT requires several SCF cycles of updating the wavefunction.

In CP2K, **iteration levels** are used to manage layers of iterative procedures and control the
frequency of printing output files. The names and descriptions of iteration levels can be found in
the input reference under any section named `EACH`, such as `QS_SCF` for the SCF cycles in Quickstep
module, `GEO_OPT` for the steps of fixed-cell geometry optimization, `CELL_OPT` for those of cell
optimization, `MD` for timesteps of molecular dynamics, etc. So for `GEO_OPT`, a value of 1 set in
the `EACH` section means printing once every 1 step of geometry optimization, 5 means once every 5
steps, and a non-positive integer means no printing at all. Different `EACH` sections can have
different default values for each iteration level, as can be checked in the input reference.

The presence of iteration levels depends primarily, but not solely, on the type of current task as
specified by [GLOBAL/RUN_TYPE](#CP2K_INPUT.GLOBAL.RUN_TYPE); for instance, a `GEO_OPT` task may
target a local minimum or a saddle point, the latter of which has the `ROT_OPT` iteration level as
defined by the dimer method in addition to the `GEO_OPT` level. To understand what iteration levels
are involved, adjust [ITER_LEVEL_VERBOSITY](#CP2K_INPUT.GLOBAL.ITER_LEVEL_VERBOSITY). The `EACH`
sections will ignore iteration levels not present in the current task even if specified. Also, while
the control of output files uses input sections with the `EACH` section, non-default settings of
`EACH` are only meaningful when the output is actually available *during* iterations; some like the
[FINAL_STRUCTURE](#CP2K_INPUT.MOTION.PRINT.FINAL_STRUCTURE) are only meant for a finished iteration.

## Motion restart

A **motion restart** file stores information of a snapshot, including atomic positions, velocities,
cell parameters, thermostat/barostat states, elastic band states, random number generator states,
etc. It is produced with the [MOTION/PRINT/RESTART](#CP2K_INPUT.MOTION.PRINT.RESTART) and
[MOTION/PRINT/RESTART_HISTORY](#CP2K_INPUT.MOTION.PRINT.RESTART_HISTORY) sections from a task with
particle motion. Take the following snippet of the input sections for a molecular dynamics (MD) task
as an example.

```text
&GLOBAL
  PRINT_LEVEL LOW
  PROJECT_NAME md_project
  RUN_TYPE MD
&END GLOBAL
&MOTION
  &MD
    ! Other important settings omitted for clarity
    STEPS 600
    STEP_START_VAL 0
  &END MD
  &PRINT
    &RESTART ON
      BACKUP_COPIES 0
      &EACH
        MD 1
      &END EACH
    &END RESTART
    &RESTART_HISTORY ON
      &EACH
        MD 200
      &END EACH
    &END RESTART_HISTORY
  &END PRINT
&END MOTION
```

The MD task will run for 600 steps, also called iterations, as specified by keyword `STEP`. The
`STEP_START_VAL` specifies that the counter of the steps starts at 0; even though 0 is just the
default value (if omitted from input), this keyword is written explicitly for later discussion. The
`STEP` value of 600 does *not* limit `STEP_START_VAL` to under 600.

After normal termination of such a task, check the presence of files using a `.restart` suffix with
the `find` command:

```shell
$ find . -type f -name "*.restart"
./md_project-1.restart
./md_project-1_200.restart
./md_project-1_400.restart
./md_project-1_600.restart
```

These are motion restart files that record information on some selected iterations. By default, they
are formatted as a regular CP2K input file and can be read and edited with common command-line tools
like `less`, `more`, `vim`, `nano`, etc. They contain not only the explicit settings provided in the
original input but also other implicit definitions of default values; moreover, the numerical values
are written in a scientific notation with many significant digits. Thus, it is expected that motion
restart file is generally longer than the original input.

The `STEP_START_VAL` value shows which iteration each restart file is at.

- The first one, `md_project-1.restart`, is produced due to setting `MOTION/PRINT/RESTART` to ON.
  Its record is from the *most recent* iteration; given `EACH/MD` set to 1, it is updated every 1
  iteration. If `EACH/MD` is set to 5, the record will be updated every 5 iterations instead.
- The other three, `md_project-1_{200,400,600}.restart`, are produced due to setting
  `MOTION/PRINT/RESTART_HISTORY` to ON. Their records are from the respective *historically visited*
  iterations; given `STEP` set to 600, `STEP_START_VAL` set to 0 and `EACH/MD` set to 200, they
  correspond to the iterations above 0 but not above (0+600) that is a multiple of 200, namely the
  200-th, 400-th, and 600-th iterations in this case.

The record in `md_project-1.restart` is *most recent* in a chronological sense, which depending on
the value of `EACH/MD` and how program is terminated may or may not be *the exact latest* iteration.
If the task has been run to completion as a normal termination, `md_project-1.restart` is the latest
record as can be verified with the `STEP_START_VAL` value, regardless of whether iterations counter
finishes at a multiple of the `EACH/MD` value. In case an out-of-memory, job scheduler interruption
or other disturbance triggers externally an erroneous termination, chances are that the unfortunate
timing does not permit rescuing the status by printing one more restart file; a smaller `EACH/MD`
value would leave a restart file more recent than a larger one would.

In practice, it is not always ideal to set the frequency to 1: a very long input file with a massive
atomistic model and complicated definition of the force method can lead to reduced performance and
significant burden on disk I/O, especially if the formatted motion restart file needs to be dumped
frequently. In such rare cases, see the hint at the end of this section.

Unlike `RESTART_HISTORY` that uses a filename containing the unique iteration counter, the `RESTART`
uses one single filename for printing on different iterations and overwrites in updates accordingly;
the filename is determined by the [FILENAME](#CP2K_INPUT.MOTION.PRINT.RESTART.FILENAME) keyword. The
[BACKUP_COPIES](#CP2K_INPUT.MOTION.PRINT.RESTART.BACKUP_COPIES) keyword accompanies the overwriting
behavior by creating backup files. In the snippet above it is explicitly set to 0 to disable backup.
If set to a positive number N instead, old files will be moved gradually out of the way and renamed
to add a suffix; in the end, apart from the most recent `md_project-1.restart`, a series of N backup
files ranging from the newer `md_project-1.restart.bak-1` to the older `md_project-1.restart.bak-N`
will also exist. While restart backup files do preserve several snapshots, they are not intended for
reconstruction of the entire trajectory from the very beginning.

A motion restart file can be utilized by directly passed to CP2K as an input file.

```shell
cp2k.psmp -i md_project-1.restart -o md_project-1.restart.out
```

A motion restart files can also be parsed in the [EXT_RESTART](#CP2K_INPUT.EXT_RESTART) section in
another freshly created input, where [RESTART_FILE_NAME](#CP2K_INPUT.EXT_RESTART.RESTART_FILE_NAME)
should point to the existing restart file.

The `EXT_RESTART` section offers fine-grained control over which components to restart. The default
is to restart all of them if available and meaningful for the current task; this can be inverted by
setting the keyword [RESTART_DEFAULT](#CP2K_INPUT.EXT_RESTART.RESTART_DEFAULT) to `.FALSE.` and then
selectively setting a few other items back to `.TRUE.`. In addition, there is the keyword
[CUSTOM_PATH](#CP2K_INPUT.EXT_RESTART.CUSTOM_PATH) that enables copying entire input sections from
the restart to the current input.

```{hint}
For situations where it is expected that formatted motion restart files will have extreme size, it
is possible to use an unformatted version produced by
[RESTART/SPLIT_RESTART_FILE](#CP2K_INPUT.MOTION.PRINT.RESTART.SPLIT_RESTART_FILE) and parsed by
[EXT_RESTART/BINARY_RESTART_FILE_NAME](#CP2K_INPUT.EXT_RESTART.BINARY_RESTART_FILE_NAME).
```

Generally, the best practices for using a motion restart is to check parameter consistency and keep
backup copies for long simulations. Values such as potential energies, temperatures, and conserved
quantities that vary continuously with steps should not show sudden jumps (discontinuity) before and
after the usage of a restart file. Motion restart file alone may not guarantee that the potential
energy is evaluated consistently if the force method is dependent on history and sensitive to the
initial guess; other restart files specific to the method are also needed, with the wavefunction
restart introduced in the next section for Quickstep DFT a primary example.

If a restarted run keeps the original `PROJECT_NAME` and is under the same directory with all of the
original output files in place, be careful that new outputs may either be appended to the end of the
existing files, or overwriting the files from the start; manual backups in other directories can
prevent accidental data loss or unintentional concatenation.

## Wavefunction restart

A **wavefunction restart** file is used by tasks with the Quickstep module for the force method,
which stores electronic structure information obtained from the SCF cycles. It is most commonly
produced with the [FORCE_EVAL/DFT/SCF/PRINT/RESTART](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.PRINT.RESTART)
and [FORCE_EVAL/DFT/SCF/PRINT/RESTART_HISTORY](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.PRINT.RESTART_HISTORY)
sections. Unlike formatted motion restart files, wavefunction restart files are **unformatted** and
are *not* intended for human reading, plaintext editing, porting to other binaries, or downstream
electronic-structure analyses. Other well-defined [](./electronic_structure/molecular_orbitals) is
available as the interface to post-processing programs for that last purpose.

The following snippet turns off the `RESTART_HISTORY` and keeps only the `RESTART` with no backups
for a task of single-point energy evaluation.

```text
&GLOBAL
  PRINT_LEVEL LOW
  PROJECT project_name
  RUN_TYPE ENERGY
&END GLOBAL
&FORCE_EVAL
  METHOD QS
  &DFT
    ! Other important settings omitted for clarity
    &SCF
      SCF_GUESS ATOMIC
      &PRINT
        &RESTART ON
          BACKUP_COPIES 0
        &END RESTART
        &RESTART_HISTORY OFF
        &END RESTART_HISTORY
      &END PRINT
    &END SCF
  &END DFT
&END FORCE_EVAL
```

The wavefunction restart file, say `project_name-RESTART.wfn`, can then be utilized by setting the
[FORCE_EVAL/DFT/WFN_RESTART_FILE_NAME](#CP2K_INPUT.FORCE_EVAL.DFT.WFN_RESTART_FILE_NAME) to it and
switching [FORCE_EVAL/DFT/SCF/SCF_GUESS](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.SCF_GUESS) to `RESTART`.

```text
&FORCE_EVAL
  &DFT
    WFN_RESTART_FILE_NAME project_name-RESTART.wfn
    &SCF
      SCF_GUESS RESTART
    &END SCF
  &END DFT
&END FORCE_EVAL
```

If `WFN_RESTART_FILE_NAME` is not explicitly given, the filename of the wavefunction restart will be
deduced from the `FORCE_EVAL/DFT/SCF/PRINT/RESTART/FILENAME` (regardless of whether the section is
explicitly given or not) and the `GLOBAL/PROJECT_NAME`. Using `SCF_GUESS RESTART` without a properly
named restart file will result in a warning and fallback to another method of initial guess.

```none
 *** WARNING in qs_initial_guess.F:312 :: User requested to restart the    ***
 *** wavefunction from the file named: project_name-RESTART.wfn. This file ***
 *** does not exist. Please check the existence of the file or change      ***
 *** properly the value of the keyword WFN_RESTART_FILE_NAME. Calculation  ***
 *** continues using ATOMIC GUESS.                                         ***
```

`SCF_GUESS RESTART` takes a single wavefunction restart file. Meanwhile, `SCF_GUESS HISTORY_RESTART`
takes a series of these files, starting with `project_name-RESTART.wfn` and continuing to backups
`project_name-RESTART.wfn.bak-1`, ..., `project_name-RESTART.wfn.bak-N`, as a positive number of
`BACKUP_COPIES` in `FORCE_EVAL/DFT/SCF/PRINT/RESTART` would generate. The number of files required
for such a history restart depends on the Quickstep extrapolation method and order, as specified by
[EXTRAPOLATION](#CP2K_INPUT.FORCE_EVAL.DFT.QS.EXTRAPOLATION) and
[EXTRAPOLATION_ORDER](#CP2K_INPUT.FORCE_EVAL.DFT.QS.EXTRAPOLATION_ORDER) keywords under the
`FORCE_EVAL/DFT/QS` section.

### Notes about k-point sampling

The usage of [k-point sampling in Brillouin zone](./dft/k-points), which is enabled with the section
[FORCE_EVAL/DFT/KPOINTS](#CP2K_INPUT.FORCE_EVAL.DFT.KPOINTS), makes a significant difference in the
file type of unformatted wavefunction restart files:

- If the section is not explicitly set, referred to as an "implicit gamma" formalism below, the
  `.wfn` suffix will be used (which is not to be confused with the plaintext AIM wavefunction file
  also with the `.wfn` suffix that several other quantum chemistry programs can output) and the
  contents include atomic basis sets and per-spin eigenvalues and occupation numbers;
- If the section is explicitly set, referred to as an "explicit k-point" formalism below, regardless
  of the actual number of k-points sampled in the scheme, the `.kp` suffix will be used and the
  content include density matrix for each spin and cell image.

The entire workflow of writing and reading wavefunction restart files are handled with separate
machinery, and thus it is not possible to swap them for restarting purposes. Also compare the
treatment of mismatching system configurations between the restart file and the input file:

- Number of atoms (`natom`) - both implicit gamma and explicit k-point will reject a mismatch and
  fallback to another method of initial guess;
- Number of spins (`nspin`) - explicit k-point will reject a mismatch, whereas implicit gamma may
  allow one but only under the restriction that (1) the number of spins in the restart is fewer than
  requested in the input, and (2) the input requests a multiplicity of 1;
- Number of basis functions in atomic orbitals (`nao`) - explicit k-point will reject a mismatch,
  implicit gamma can allow one;
- Number of molecular orbitals (`nmo`) - specific to implicit gamma, some truncation or padding with
  zeroes will be carried out in case of a MO mismatch, but if the number of occupied orbitals in the
  restart as indicated by the index of HOMO is higher than total number of available MOs from input
  due to more electrons of a spin in the restart than input, then restarting is rejected.
- Number of cell images (`nimages`) - specific to explicit k-point, a mismatch is allowed.

For explicit k-point sampling, the [FULL_GRID](#CP2K_INPUT.FORCE_EVAL.DFT.KPOINTS.FULL_GRID) keyword
in the `KPOINTS` section controls whether all k-points (including symmetry-equivalent ones) or only
irreducible k-points should be used and written to the restart. `FULL_GRID ON` uses all, which has
maximum restart flexibility; the default `FULL_GRID OFF` is faster but restricted to the same
symmetry settings. K-point restart files generated with `FULL_GRID ON` can be used with
`FULL_GRID OFF` in the restart calculation (or vice versa). Both converge to the same energy,
typically in 1 SCF step.

A complete example of k-point samping and restart files is demonstrated below. Despite distinctions
described above, the relevant input keywords are otherwise identical with the implicit gamma case.

```text
&GLOBAL
  PROJECT H_SYM
  RUN_TYPE ENERGY
&END GLOBAL

&FORCE_EVAL
  &DFT
    BASIS_SET_FILE_NAME GTH_BASIS_SETS
    POTENTIAL_FILE_NAME POTENTIAL
    &KPOINTS
      FULL_GRID ON
      SCHEME MONKHORST-PACK 2 2 2
      SYMMETRY ON
    &END KPOINTS
    &SCF
      SCF_GUESS ATOMIC
      &PRINT
        &RESTART ON
        &END RESTART
      &END PRINT
    &END SCF
    &XC
      &XC_FUNCTIONAL PADE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC 3.56683 3.56683 3.56683
    &END CELL
    &COORD
      SCALED
      H     0.100000    0.000000    0.000000
      H     0.500000    0.500000    0.000000
      ...
    &END COORD
    &KIND H
      BASIS_SET SZV-GTH
      POTENTIAL GTH-PADE-q1
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
```

This generates `H_SYM-RESTART.kp`. Restarting with a different `FULL_GRID` setting:

```text
&DFT
  WFN_RESTART_FILE_NAME H_SYM-RESTART.kp
  &KPOINTS
    FULL_GRID OFF  ! Different from initial calculation
    SCHEME MONKHORST-PACK 2 2 2
    SYMMETRY ON
  &END KPOINTS
  &SCF
    SCF_GUESS RESTART
    MAX_SCF 3
    &PRINT
      &RESTART OFF
      &END RESTART
    &END PRINT
  &END SCF
  &XC
    &XC_FUNCTIONAL PADE
    &END XC_FUNCTIONAL
  &END XC
&END DFT
```

The CP2K test suite includes k-point restart tests in `tests/QS/regtest-kp-1/`.

### Harris Chain

The Harris functional converts a k-point wavefunction into a gamma-point wavefunction, enabling
k-point-converged densities to seed gamma-point calculations. The three-step chain is: `k-point SCF`
→ `k-point Harris functional` → `gamma-point DFT`.

**Step 1: K-point SCF** — generates `{PROJECT}-1_0.kp`:

```none
&GLOBAL
  PROJECT Carbon
  RUN_TYPE ENERGY
&END GLOBAL

&FORCE_EVAL
  &DFT
    BASIS_SET_FILE_NAME BASIS_SET
    POTENTIAL_FILE_NAME GTH_POTENTIALS
    &KPOINTS
      FULL_GRID ON
      SCHEME MONKHORST-PACK 2 2 2
      SYMMETRY ON
    &END KPOINTS
    &SCF
      SCF_GUESS ATOMIC
      &PRINT
        &RESTART ON
        &END RESTART
      &END PRINT
    &END SCF
    &XC
      &XC_FUNCTIONAL PADE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC 3.56683 3.56683 3.56683
    &END CELL
    &COORD ...  
    &KIND C
      BASIS_SET ORB DZVP-GTH-PADE
      POTENTIAL GTH-PADE-q4
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
```

This generates `Carbon-1_0.kp`.

**Step 2: Harris functional with k-points** — reads the `.kp` file and writes a gamma-point
wavefunction to `{PROJECT}-Harris-1_0.kp` via `HARRIS_OUTPUT_WFN`:

```none
&GLOBAL
  PROJECT Carbon
  RUN_TYPE ENERGY
&END GLOBAL

&FORCE_EVAL
  &DFT
    BASIS_SET_FILE_NAME BASIS_SET
    POTENTIAL_FILE_NAME GTH_POTENTIALS
    WFN_RESTART_FILE_NAME ./Carbon-1_0.kp
    &ENERGY_CORRECTION
      ENERGY_FUNCTIONAL HARRIS
      HARRIS_BASIS HARRIS
      &KPOINTS
        FULL_GRID ON
        SCHEME MONKHORST-PACK 2 2 2
        SYMMETRY ON
      &END KPOINTS
      &PRINT
        &HARRIS_OUTPUT_WFN
        &END HARRIS_OUTPUT_WFN
      &END PRINT
      &XC
        &XC_FUNCTIONAL PBE
        &END XC_FUNCTIONAL
      &END XC
    &END ENERGY_CORRECTION
    &SCF
      SCF_GUESS RESTART
      &PRINT
        &RESTART OFF
      &END RESTART
    &END PRINT
    &XC
      &XC_FUNCTIONAL PADE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
  &SUBSYS
    &KIND C
      BASIS_SET ORB DZVP-GTH-PADE
      BASIS_SET HARRIS SZV-GTH-PADE
      POTENTIAL GTH-PADE-q4
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
```

Key points:

- `ENERGY_CORRECTION` + `ENERGY_FUNCTIONAL HARRIS` enables the Harris functional
- `HARRIS_BASIS HARRIS` uses the smaller Harris basis set
- `KPOINTS` in `ENERGY_CORRECTION` must match the step 1 grid
- `HARRIS_OUTPUT_WFN` writes `{PROJECT}-Harris-1_0.kp`
- The `KIND` section needs both the orbital and the `HARRIS` basis set

This generates `Carbon-Harris-1_0.kp`.

**Step 3: Gamma-point DFT** — reads the Harris wavefunction and optionally writes a `.wfn` file:

```none
&GLOBAL
  PROJECT Carbon_gamma
  RUN_TYPE ENERGY
&END GLOBAL

&FORCE_EVAL
  &DFT
    BASIS_SET_FILE_NAME BASIS_SET
    POTENTIAL_FILE_NAME GTH_POTENTIALS
    WFN_RESTART_FILE_NAME ./Carbon-Harris-1_0.kp
    &SCF
      SCF_GUESS RESTART
      &PRINT
        &RESTART ON
      &END RESTART
    &END PRINT
    &XC
      &XC_FUNCTIONAL PADE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
  &SUBSYS
    &KIND C
      BASIS_SET ORB DZVP-GTH-PADE
      POTENTIAL GTH-PADE-q4
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
```

- No `KPOINTS` section — this is a gamma-point calculation
- `&RESTART ON` writes `Carbon_gamma-RESTART.wfn`

The CP2K test suite includes these Harris chain tests in `tests/QS/regtest-harris-kp/`
(`cc_kp_01.inp`, `cc_kp_02.inp`, `cc_kp_03.inp`) which demonstrate the three-step k-point restart
process.

### Band and NEB Calculations

A band calculation writes a separate wavefunction restart for every replica. CP2K derives the
replica project names from the project in `GLOBAL`, for example:

```none
neb-BAND1-RESTART.wfn
neb-BAND2-RESTART.wfn
neb-BAND3-RESTART.wfn
```

With k-point sampling, the corresponding files have the `.kp` extension. To reuse these files, set
`SCF_GUESS RESTART` but normally omit `WFN_RESTART_FILE_NAME`:

```none
&GLOBAL
  PROJECT neb
  RUN_TYPE BAND
&END GLOBAL

&EXT_RESTART
  RESTART_FILE_NAME neb-1.restart
&END EXT_RESTART

&FORCE_EVAL
  &DFT
    # Do not set WFN_RESTART_FILE_NAME here.
    &SCF
      SCF_GUESS RESTART
    &END SCF
  &END DFT
  ...
&END FORCE_EVAL
```

CP2K then selects `neb-BAND<N>-RESTART.wfn` or `neb-BAND<N>-RESTART.kp` after assigning the current
replica project name. In contrast, an explicitly specified `WFN_RESTART_FILE_NAME` is interpreted
literally and is shared by all replicas; CP2K does not insert the `-BAND<N>` part. An explicit name
is therefore useful when all replicas should deliberately start from the same wavefunction, but a
name such as `neb-RESTART.kp` will not select the per-replica files shown above.

The `EXT_RESTART` file restores the band coordinates and optimizer state. This is independent of the
per-replica wavefunction restarts, so use both `EXT_RESTART` and `SCF_GUESS RESTART` when both parts
of an interrupted band calculation should be continued. Keep the original `PROJECT` name and run
from the directory containing the restart files, or preserve their relative paths.

### CDFT

For mixed CDFT calculations, each constrained state requires its own wavefunction restart file.
Specify `SCF_GUESS RESTART` and `WFN_RESTART_FILE_NAME` per `FORCE_EVAL` subblock:

```none
&FORCE_EVAL
  METHOD MIXED
  &MIXED
    MIXING_TYPE MIXED_CDFT
    &MIXED_CDFT
      WFN_OVERLAP TRUE
      WFN_RESTART_FILE_NAME reference_wavefunction.wfn
    &END MIXED_CDFT
  &END MIXED
&END FORCE_EVAL
```

Mixed CDFT calculations require separate wavefunction restart files for each constrained state. The
CP2K test suite includes examples in `tests/QS/regtest-cdft-3/`: single-state CDFT calculations
generate wavefunction files that are then used as restarts for mixed CDFT calculations.
