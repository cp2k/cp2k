# Restarting CP2K Calculations

This guide explains how to restart CP2K calculations from previously saved wavefunction files.
Restarting is useful for continuing calculations that were interrupted, extending molecular dynamics
simulations, or using converged wavefunctions as starting points for new calculations.

## Basic

### Restarting a Calculation

To restart a CP2K calculation, you need a previously saved wavefunction file (typically with `.wfn`
extension), and two keywords in your input:

```none
&DFT
  WFN_RESTART_FILE_NAME your_restart_file.wfn
  &SCF
    SCF_GUESS RESTART
  &END SCF
&END DFT
```

`WFN_RESTART_FILE_NAME` points to the wavefunction file, and `SCF_GUESS RESTART` instructs CP2K to
use it as the initial guess.

### Generating Restart Files

To produce a wavefunction restart file during a calculation, configure the `RESTART` print section:

```none
  &PRINT
    &RESTART
      FILENAME ./your_project_name
      BACKUP_COPIES 3
      COMMON_ITERATION_LEVELS 1
      &EACH
        JUST_ENERGY 1  ! Write at end of each SCF (produces your_project_name-RESTART.wfn)
        QS_SCF 0
      &END EACH
    &END RESTART
  &END PRINT
```

The project name is set in the GLOBAL section.

```none
&GLOBAL
  PROJECT your_project_name
&END GLOBAL
```

### Key Keywords

- **`FILENAME`**: Base name for restart files (suffix is `-RESTART.wfn` or `-{step}_0.wfn`)
- **`BACKUP_COPIES`**: Number of backup copies to retain (default: 1)
- **`COMMON_ITERATION_LEVELS`**: How many iteration levels share a common filename (default: 0)
- **`JUST_ENERGY 1`**: Write at end of SCF → produces `{FILENAME}-RESTART.wfn`
- **`QS_SCF N`**: Write every N SCF steps → produces `{FILENAME}-{i}_0.wfn`

For geometry optimizations or MD, write at each step:

```none
&EACH
  GEO_OPT 1  ! or MD 1
  QS_SCF 0
&END EACH
```

Best practices:

1. Use descriptive filenames that include the project name and step number.
1. Verify the original calculation converged properly before restarting.
1. Ensure the restart uses the same basis sets, functionals, and parameters.

### Wavefunction File Format

CP2K wavefunction restart files (`*.wfn`) are binary and architecture-specific. The exact format may
vary between CP2K versions — use restart files with the same CP2K version that generated them.

### Troubleshooting

- **Incompatible Parameters**: Ensure the restart uses identical settings to the original (basis
  sets, k-points, XC functional, etc.).
- **File Not Found**: Verify the restart file path is correct and the file exists.
- **Corrupted Files**: Try an earlier backup copy.

## K-Points

For periodic systems using k-point sampling, CP2K stores wavefunction restart files with the `.kp`
extension. The procedure is the same as the basic case but uses a different file extension.

The key parameter is `FULL_GRID` in the `KPOINTS` section:

- **`FULL_GRID ON`**: Writes all k-points (including symmetry-equivalent ones). Use this for maximum
  restart flexibility.
- **`FULL_GRID OFF`** (default): Writes only irreducible k-points. Faster but restricted to the same
  symmetry settings.

### Generating K-Point Restart Files

```none
&KPOINTS
  FULL_GRID ON
  SCHEME MONKHORST-PACK 2 2 2
  SYMMETRY ON
&END KPOINTS
```

Using a k-point restart file:

```none
&DFT
  WFN_RESTART_FILE_NAME your_project-RESTART.kp
  &KPOINTS
    FULL_GRID ON
    SCHEME MONKHORST-PACK 2 2 2
    SYMMETRY ON
  &END KPOINTS
  &SCF
    SCF_GUESS RESTART
  &END PRINT
  ...
&END DFT
```

K-point restart files generated with `FULL_GRID ON` can be used with `FULL_GRID OFF` in the restart
calculation (or vice versa). Both converge to the same energy, typically in 1 SCF step.

### Complete Example

```none
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

```none
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
  &XC
    &XC_FUNCTIONAL PADE
    &END XC_FUNCTIONAL
  &END XC
&END DFT
```

The CP2K test suite includes k-point restart tests in `tests/QS/regtest-kp-1/`.

## Harris Chain

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

## Molecular Dynamics

For continuing molecular dynamics simulations, CP2K stores positions, velocities, cell parameters,
and thermostat/barostat state in `.restart` files. Use the `EXT_RESTART` section for fine-grained
control over which components to restore.

Keywords:

- **`RESTART_FILE_NAME`** / **`EXTERNAL_FILE`**: Restart file to read
- **`RESTART_DEFAULT`**: Set all `RESTART_*` options at once (default: TRUE)
- **`RESTART_POS`**: Restart positions (bare keyword = TRUE)
- **`RESTART_VEL`**: Restart velocities (bare keyword = TRUE)
- **`RESTART_CELL`**: Restart cell parameters (bare keyword = TRUE)
- **`RESTART_THERMOSTAT`**: Restart thermostat state (bare keyword = TRUE)
- **`RESTART_BAROSTAT`**: Restart barostat state (bare keyword = TRUE)
- **`RESTART_RANDOMG`**: Restart RNG state (bare keyword = TRUE)
- **`RESTART_COUNTERS`**: Restart step counter and walltime (bare keyword = TRUE)
- **`RESTART_BAROSTAT_THERMOSTAT`**: Restart barostat thermostat (bare keyword = TRUE)
- **`RESTART_SHELL_POS` / `RESTART_CORE_POS`**: Restart shell-model positions
- **`RESTART_SHELL_VELOCITY` / `RESTART_CORE_VELOCITY`**: Restart shell-model velocities
- **`RESTART_SHELL_THERMOSTAT`**: Restart shell thermostat

Common scenarios:

1. **Restart positions only (fresh velocity distribution):**

   ```none
   &EXT_RESTART
     EXTERNAL_FILE equilibration.restart
     RESTART_POS
     RESTART_VEL FALSE
   &END EXT_RESTART
   ```

1. **Full restart:**

   ```none
   &EXT_RESTART
     EXTERNAL_FILE previous_md.restart
     RESTART_POS
     RESTART_VEL
     RESTART_CELL
   &END EXT_RESTART
   ```

1. **NPT restart:**

   ```none
   &EXT_RESTART
     EXTERNAL_FILE npt_equil.restart
     RESTART_POS
     RESTART_VEL
     RESTART_CELL
     RESTART_BAROSTAT .FALSE.
   &END EXT_RESTART


   ```

Best practices: check parameter consistency, keep backup copies for long simulations, and validate
energies and temperatures after restarting.

## CDFT

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

## References

For more information about restarting calculations, see:

- The [CP2K Input Reference](https://manual.cp2k.org/trunk/CP2K_INPUT.html) for detailed keyword
  documentation
- Example input files in the `tests` directory of the CP2K distribution
- The [CP2K Forum](https://groups.google.com/group/cp2k) for community support
