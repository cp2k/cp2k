# Restarting CP2K Calculations

This guide explains how to restart CP2K calculations from previously saved wavefunction files.
Restarting is useful for continuing calculations that were interrupted, extending molecular dynamics
simulations, or using converged wavefunctions as starting points for new calculations.

## Basic Restart Procedure

To restart a CP2K calculation, you need to:

1. Have a previously saved wavefunction file (typically with `.wfn` extension)
1. Configure your input file to read this wavefunction
1. Set the appropriate SCF guess method

## Input File Configuration

### Wavefunction Restart File

The main keyword for restarting is `WFN_RESTART_FILE_NAME` in the `DFT` section:

```none
&DFT
  ...
  WFN_RESTART_FILE_NAME your_restart_file.wfn
  ...
&END DFT
```

### SCF Guess Method

You should also set the SCF guess method to `RESTART`:

```none
&SCF
  SCF_GUESS RESTART
  ...
&END SCF
```

## Complete Example

Here's a complete example showing how to restart a calculation:

```none
@SET RESTART_WFN          TRUE
@SET WFN_FILE             previous_calculation-1_0.wfn

&GLOBAL
  PROJECT my_restarted_calculation
  RUN_TYPE ENERGY
&END GLOBAL

&FORCE_EVAL
  METHOD QS
  &DFT
    BASIS_SET_FILE_NAME  BASIS_MOLOPT
    POTENTIAL_FILE_NAME  POTENTIAL
    @IF ( ${RESTART_WFN} == TRUE )
      WFN_RESTART_FILE_NAME ${WFN_FILE}
    @ENDIF
    CHARGE 0
    &MGRID
      CUTOFF 100
      NGRIDS 5
    &END MGRID
    &SCF
      @IF ( ${RESTART_WFN} == TRUE )
        SCF_GUESS RESTART
      @ENDIF
      @IF ( ${RESTART_WFN} == FALSE )
        SCF_GUESS ATOMIC
      @ENDIF
      EPS_SCF 1.0E-5
      MAX_SCF 20
      &OT ON
        MINIMIZER CG
        PRECONDITIONER FULL_ALL
      &END OT
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
&END FORCE_EVAL
```

## Restart File Generation

To generate restart files during a calculation, you need to configure the `RESTART` section within
the `PRINT` subsection of `SCF`. This will create wavefunction files that can be used later for
restarting calculations.

### Basic Restart File Generation

The most common approach is to write a restart file at the end of a successful SCF calculation:

```none
&SCF
  ...
  &PRINT
    &RESTART
      FILENAME ./your_project_name
      BACKUP_COPIES 0
      COMMON_ITERATION_LEVELS 1
      &EACH
        JUST_ENERGY 1  ! Write restart file when energy is computed (at end of SCF)
        QS_SCF 0
      &END EACH
    &END RESTART
  &END PRINT
  ...
&END SCF
```

### Key Keywords for Restart File Generation

1. **`FILENAME`**: Sets the base filename for restart files (wavefunction files will have this
   prefix, plus `-RESTART.wfn` or `-{step}_0.wfn` suffix)
1. **`BACKUP_COPIES`**: Controls how many backup copies to keep (default: 1; 0 = no backups)
1. **`COMMON_ITERATION_LEVELS`**: Determines how many iteration levels are merged into a common
   filename (default: 0 = each level gets its own filename)
1. **`JUST_ENERGY`**: When set to 1, writes restart file at the end of SCF when energy is computed
   (produces `{FILENAME}-RESTART.wfn`)
1. **`QS_SCF N`**: Write restart file every N SCF iterations (produces `{FILENAME}-{i}_0.wfn`)

### Complete Example with File Generation

Here's a complete example that generates a restart file named `previous_calculation-RESTART.wfn`
(matching the filename expected by the restart examples above):

```none
@SET PROJECT_NAME previous_calculation

&GLOBAL
  PROJECT ${PROJECT_NAME}
  RUN_TYPE ENERGY
&END GLOBAL

&FORCE_EVAL
  METHOD QS
  &DFT
    BASIS_SET_FILE_NAME  BASIS_MOLOPT
    POTENTIAL_FILE_NAME  POTENTIAL
    CHARGE 0
    &MGRID
      CUTOFF 100
      NGRIDS 5
    &END MGRID
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF 1.0E-5
      MAX_SCF 20
      &OT ON
        MINIMIZER CG
        PRECONDITIONER FULL_ALL
      &END OT
      &PRINT
        &RESTART
          FILENAME ./${PROJECT_NAME}
          BACKUP_COPIES 0
          COMMON_ITERATION_LEVELS 1
          &EACH
            JUST_ENERGY 1  ! This will generate ${PROJECT_NAME}-1_0.wfn
            QS_SCF 0
          &END EACH
        &END RESTART
      &END PRINT
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
&END FORCE_EVAL
```

This configuration will generate a restart file named `previous_calculation-RESTART.wfn` (the
`-RESTART.wfn` suffix is produced when `JUST_ENERGY 1` is set). This file can be used in the restart
examples shown earlier in this documentation.

### Controlling Restart File Frequency

You can control how often restart files are written:

```none
&EACH
  QS_SCF 10  ! Write restart file every 10 SCF steps
  JUST_ENERGY 0
&END EACH
```

Or write at specific intervals:

```none
&EACH
  QS_SCF 0
  JUST_ENERGY 1  ! Only write at the end when energy is computed
&END EACH
```

### Multiple Restart Files

For geometry optimizations or molecular dynamics, you might want to write restart files at each
optimization step:

```none
&SCF
  &PRINT
    &RESTART
      FILENAME ./${PROJECT_NAME}
      BACKUP_COPIES 3  ! Keep 3 backup copies
      &EACH
        GEO_OPT 1  ! Write restart file at each geometry optimization step
      &END EACH
    &END RESTART
  &END PRINT
&END SCF
```

## Best Practices

1. **File Naming**: Use descriptive names for your restart files that include the project name and
   possibly the step number.

1. **Backup Copies**: Consider keeping backup copies of restart files by setting `BACKUP_COPIES` to
   a positive number.

1. **Convergence Check**: Before restarting, verify that the original calculation converged
   properly.

1. **Consistency**: Ensure that the restart calculation uses the same basis sets, functionals, and
   other parameters as the original calculation.

1. **Performance**: Restarting from a converged wavefunction can significantly reduce the number of
   SCF iterations needed for convergence.

## Common Use Cases

### Continuing a Previous Calculation

When you want to continue a calculation that was stopped (either intentionally or due to an
interruption):

```bash
cp2k.sopt previous_input.inp --restart previous_output.wfn
```

### Using as Initial Guess

When you want to use a converged wavefunction as an initial guess for a similar system or slightly
modified calculation:

```none
&DFT
  WFN_RESTART_FILE_NAME good_initial_guess.wfn
  SCF_GUESS RESTART
  ...
&END DFT
```

### Molecular Dynamics Restarts

For continuing molecular dynamics simulations, the restart procedure is similar but typically
involves additional files for velocities and positions.

## Troubleshooting

- **Incompatible Parameters**: If you get errors about incompatible parameters, ensure your restart
  calculation uses identical settings to the original.

- **File Not Found**: Verify the restart file path is correct and the file exists.

- **Corrupted Files**: If a restart file appears corrupted, try using an earlier backup copy if
  available.

- **Performance Issues**: If restarting doesn't improve performance, the wavefunction may not be a
  good guess for the new system.

## Advanced Topics

### Selective Restarting

In some cases, you may want to restart only certain parts of a calculation. CP2K provides options
for this through various restart-related keywords.

### Mixed CDFT Restarts

For mixed CDFT calculations, each state can have its own restart file:

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
  ...
&END FORCE_EVAL
```

Mixed CDFT calculations require separate wavefunction restart files for each constrained state. The
CP2K test suite includes such tests in `tests/QS/regtest-cdft-3/`, where single-state CDFT
calculations (e.g., `HeH-cdft-state-1.inp`, `HeH-cdft-state-2.inp`) generate wavefunction files that
are then used as restarts for mixed CDFT calculations (e.g., `HeH-mixed-cdft-1.inp`):

```none
@SET WFN_FILE_1  HeH-cdft-state-1-1_0.wfn  ! Restart file for state 1
@SET WFN_FILE_2  HeH-cdft-state-2-1_0.wfn  ! Restart file for state 2
```

Each force eval subblock uses `WFN_RESTART_FILE_NAME ${WFN_FILE_1}` or `${WFN_FILE_2}` accordingly,
with `SCF_GUESS RESTART` to restart from the respective state's wavefunction.

### Restart File Format

CP2K restart files are binary files containing the wavefunction coefficients and other essential
information. The exact format may vary between CP2K versions, so it's generally recommended to use
restart files with the same version of CP2K that generated them.

### External Restart Files for Molecular Dynamics

For molecular dynamics simulations, CP2K provides the `EXT_RESTART` section which allows you to
selectively restart specific components from external restart files. This is particularly useful
when you want to continue a trajectory with different parameters (e.g., switching from NVT to NVE
ensemble).

#### Basic EXT_RESTART Usage

The `EXT_RESTART` section allows fine-grained control over what to restart:

```none
&EXT_RESTART
  EXTERNAL_FILE your_restart_file.restart
  RESTART_POS                    ! Restart atomic positions
  RESTART_VEL FALSE             ! Don't restart velocities (start from 0)
  RESTART_CELL                  ! Restart cell parameters
  RESTART_THERMOSTAT .FALSE.    ! Don't restart thermostat state
  RESTART_BAROSTAT .FALSE.      ! Don't restart barostat state
  RESTART_RANDOMG .FALSE.       ! Don't restart random number generator
&END EXT_RESTART
```

#### Common EXT_RESTART Scenarios

1. **Restarting positions only (for new velocity distribution):**

   ```none
   &EXT_RESTART
     EXTERNAL_FILE equilibration.restart
     RESTART_POS
     RESTART_VEL FALSE
   &END EXT_RESTART
   ```

1. **Full restart for continuing MD:**

   ```none
   &EXT_RESTART
     EXTERNAL_FILE previous_md.restart
     RESTART_POS
     RESTART_VEL
     RESTART_CELL
   &END EXT_RESTART
   ```

1. **Restarting cell for NPT simulations:**

   ```none
   &EXT_RESTART
     EXTERNAL_FILE npt_equil.restart
     RESTART_POS
     RESTART_VEL
     RESTART_CELL
     RESTART_BAROSTAT .FALSE.
   &END EXT_RESTART
   ```

#### EXT_RESTART Keywords

- **`RESTART_FILE_NAME`** (alias: `EXTERNAL_FILE`): Specifies the restart file to read from
- **`RESTART_DEFAULT`**: Sets all `RESTART_*` options to the same value at once (default: TRUE)
- **`RESTART_COUNTERS`**: Restart MD counters (step number, walltime, etc.) (bare keyword = TRUE)
- **`RESTART_POS`**: Restart atomic positions (bare keyword = TRUE)
- **`RESTART_VEL`**: Restart atomic velocities (bare keyword = TRUE)
- **`RESTART_RANDOMG`**: Restart random number generator state (bare keyword = TRUE)
- **`RESTART_BAROSTAT`**: Restart barostat state (bare keyword = TRUE)
- **`RESTART_BAROSTAT_THERMOSTAT`**: Restart barostat thermostat state (bare keyword = TRUE)
- **`RESTART_THERMOSTAT`**: Restart thermostat state (bare keyword = TRUE)
- **`RESTART_CELL`**: Restart cell parameters (bare keyword = TRUE)
- **`RESTART_SHELL_POS`**: Restart shell positions (for shell-model) (bare keyword = TRUE)
- **`RESTART_CORE_POS`**: Restart core positions (for shell-model) (bare keyword = TRUE)
- **`RESTART_SHELL_VELOCITY`**: Restart shell velocities (bare keyword = TRUE)
- **`RESTART_CORE_VELOCITY`**: Restart core velocities (bare keyword = TRUE)
- **`RESTART_SHELL_THERMOSTAT`**: Restart shell thermostat (bare keyword = TRUE)
- **`RESTART_OPTIMIZE_INPUT_VARIABLES`**: Restart optimize input variables (bare keyword = TRUE)

#### Example: Equilibration to Production Workflow

This example shows how to use EXT_RESTART to continue from an equilibration run to a production run:

```none
&GLOBAL
  PROJECT production_run
  RUN_TYPE MD
&END GLOBAL

&EXT_RESTART
  EXTERNAL_FILE equilibration-1.restart
  RESTART_POS
  RESTART_VEL
  RESTART_CELL
  RESTART_THERMOSTAT .FALSE.    ! Start with fresh thermostat
  RESTART_BAROSTAT .FALSE.     ! Start with fresh barostat
&END EXT_RESTART

&MOTION
  &MD
    ENSEMBLE NVE
    STEPS 10000
    TEMPERATURE 300
    TIMESTEP 2.0
  &END MD
  &PRINT
    &RESTART
      &EACH
        MD 100  ! Write restart file every 100 steps
      &END EACH
    &END RESTART
  &END PRINT
&END MOTION

&FORCE_EVAL
  METHOD QS
  &DFT
    ! Your DFT settings
    BASIS_SET_FILE_NAME BASIS_MOLOPT
    POTENTIAL_FILE_NAME POTENTIAL
    &SCF
      SCF_GUESS RESTART
      &PRINT
        &RESTART
          FILENAME ./production_run
          BACKUP_COPIES 3
        &END RESTART
      &END PRINT
    &END SCF
  &END DFT
&END FORCE_EVAL
```

### Architecture-Specific Wavefunction Files

As mentioned in the CP2K.org documentation, wavefunction restart files (`*.wfn` files) are binary
and architecture-specific. If you need to transfer restart files between different architectures or
CP2K versions, you may need to use the tools provided in `cp2k/tools/RestartTools`.

### Best Practices for MD Restarts

1. **Consistency Check**: When restarting MD simulations, ensure the new run uses compatible
   parameters with the original simulation.

1. **Selective Restarting**: Use `EXT_RESTART` to carefully control which components to restart,
   especially when changing simulation parameters.

1. **File Organization**: Keep restart files organized with clear naming conventions that include
   the simulation type and step number.

1. **Backup Strategy**: For long MD simulations, maintain multiple backup copies of restart files
   using `BACKUP_COPIES`.

1. **Validation**: After restarting, always validate that the continuation produces reasonable
   results by checking energies, temperatures, and other key properties.

## References

For more information about restarting calculations, see:

- The [CP2K Input Reference](https://manual.cp2k.org/trunk/CP2K_INPUT.html) for detailed keyword
  documentation
- Example input files in the `tests` directory of the CP2K distribution
- The [CP2K Forum](https://groups.google.com/group/cp2k) for community support
