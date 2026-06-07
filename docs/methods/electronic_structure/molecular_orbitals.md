# Molecular orbitals output

CP2K can write Kohn--Sham molecular orbital (MO) information in several formats. These outputs serve
different purposes: some are convenient for quick inspection in the standard output, some are
intended for visualization, and some provide enough basis-set information for external
post-processing.

The most useful choices are:

| Output              | Print section                                              | Main use                                                                 | Notes                                                                                                                           |
| ------------------- | ---------------------------------------------------------- | ------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------- |
| Text MO information | [&MO](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO)                 | Inspect MO energies, occupations, and AO coefficients in the output file | Good for small systems or debugging; files can become large if eigenvectors are printed.                                        |
| Molden file         | [&MO_MOLDEN](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO_MOLDEN)   | Visualize Gamma-only orbitals with external tools                        | Not available for k-point calculations. CP2K can optionally write cell and pseudopotential metadata.                            |
| K-point MO dump     | [&MO_KP](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO_KP)           | Export MO information for k-point calculations                           | CP2K-specific `.mokp` text format. It includes k-points, eigenvalues, occupations, complex MO coefficients, and AO information. |
| MO cube files       | [&MO_CUBES](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO_CUBES)     | Visualize selected orbitals on a real-space grid                         | Useful for frontier-orbital plots; may produce many large files.                                                                |
| MO openPMD files    | [&MO_OPENPMD](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO_OPENPMD) | Grid-based orbital output in openPMD format                              | More suitable for structured, possibly large, grid data than plain cube files.                                                  |

## Text MO information

The [&MO](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO) print section writes MO information to a text file.
It is most useful when you want to inspect a limited number of orbitals, their energies,
occupations, and optionally their AO coefficients.

```text
&FORCE_EVAL
  &DFT
    &PRINT
      &MO
        EIGENVALUES T
        OCCUPATION_NUMBERS T
        COEFFICIENTS F
        MO_INDEX_RANGE 1 -1
        NDIGITS 6
      &END MO
    &END PRINT
  &END DFT
&END FORCE_EVAL
```

Useful controls include:

- [EIGENVALUES](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO.ENERGIES),
  [OCCUPATION_NUMBERS](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO.OCCUPATION_NUMBERS), and
  [COEFFICIENTS](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO.COEFFICIENTS): choose which parts of the MO
  information are printed.
- [MO_INDEX_RANGE](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO.MO_INDEX_RANGE): restrict the printed orbital
  range. Use `-1` as the last index to print all available orbitals.
- [CARTESIAN](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO.CARTESIAN): print coefficients in the Cartesian
  basis instead of the default spherical basis.
- [CARTESIAN_OVERLAP](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO.CARTESIAN_OVERLAP): print the Cartesian
  overlap matrix together with Cartesian MO coefficients.

For diagonalization calculations, unoccupied orbitals must be available through
[ADDED_MOS](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.ADDED_MOS). For OT calculations, CP2K can generate the
requested virtual orbitals for printing through the OT eigensolver.

## Molden file for Gamma-only calculations

[&MO_MOLDEN](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO_MOLDEN) writes the orbitals in Molden format. This
is the most convenient format for visualizing molecular orbitals from Gamma-only calculations.

```text
&FORCE_EVAL
  &DFT
    &PRINT
      &MO_MOLDEN
        NDIGITS 9
        UNIT ANGSTROM
        WRITE_CELL T
        WRITE_PSEUDO T
        MARK_GHOST T
        NLUMO 5
      &END MO_MOLDEN
    &END PRINT
  &END DFT
&END FORCE_EVAL
```

The Molden format was originally designed for molecular systems. CP2K therefore adds optional
metadata that is useful for periodic or pseudopotential calculations:

- [WRITE_CELL](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO_MOLDEN.WRITE_CELL): write a `[Cell]` block with
  the simulation-cell vectors.
- [WRITE_PSEUDO](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO_MOLDEN.WRITE_PSEUDO): write a `[Pseudo]` block
  with effective nuclear charges.
- [MARK_GHOST](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO_MOLDEN.MARK_GHOST): mark ghost atoms in the
  `[Atoms]` block by setting their atomic number to zero.
- [NLUMO](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO_MOLDEN.NLUMO): include a selected number of unoccupied
  orbitals. `0` means no virtual orbitals, and `-1` means all available virtual orbitals.
- [GTO_KIND](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO_MOLDEN.GTO_KIND): choose spherical or Cartesian
  Gaussian functions. The default is spherical.

```{note}
`WRITE_PSEUDO` already writes the effective nuclear charge of ghost atoms as zero. Therefore,
`MARK_GHOST` is mainly needed when a visualization or post-processing tool relies on the atomic
number in the `[Atoms]` block rather than the `[Pseudo]` block.
```

```{note}
Molden output is not available for k-point calculations. Use
[&MO_KP](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO_KP) instead.
```

## K-point MO output: `.mokp`

From version 2026.2, CP2K starts to write `.mokp` text format for k-point calculations through
[&MO_KP](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO_KP). This format is intended for post-processing tools
that need the k-point-resolved MO coefficients together with enough AO information to interpret
them.

```text
&FORCE_EVAL
  &DFT
    &PRINT
      &MO_KP
        AO_EXPORT_TYPE GTO_BASIS
        NDIGITS 9
        UNIT ANGSTROM
      &END MO_KP
    &END PRINT
  &END DFT
&END FORCE_EVAL
```

A `.mokp` file contains:

- cell vectors, atomic positions, nuclear charges, effective nuclear charges, and AO ranges;
- the k-point list and k-point weights;
- eigenvalues and occupations for each k-point and spin channel;
- real and, when needed, imaginary parts of the MO coefficient matrix;
- either a GTO basis definition or explicit overlap matrices, controlled by
  [AO_EXPORT_TYPE](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO_KP.AO_EXPORT_TYPE).

The two AO export modes target different post-processing workflows:

- `AO_EXPORT_TYPE GTO_BASIS` writes a compact GTO basis definition. The post-processing code
  reconstructs the overlap matrix, including its k-point dependence.
- `AO_EXPORT_TYPE OVERLAP_MATRIX` writes the explicit overlap matrix $S(k)$ for each k-point. This
  produces larger files but avoids reconstructing the AO overlap externally.

```{note}
`&MO_KP` is only active for k-point calculations. For Gamma-only calculations, use `&MO_MOLDEN`,
`&MO`, or `&MO_CUBES` depending on the intended use.
```

## MO cube files

[&MO_CUBES](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO_CUBES) writes selected orbitals on a real-space grid
in cube format. This is mainly useful for visualizing frontier orbitals, defect states, adsorbate
states, or localized molecular orbitals.

```text
&FORCE_EVAL
  &DFT
    &PRINT
      &MO_CUBES
        NHOMO 2
        NLUMO 5
        STRIDE 2 2 2
      &END MO_CUBES
    &END PRINT
  &END DFT
&END FORCE_EVAL
```

Useful controls include:

- [NHOMO](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO_CUBES.NHOMO) and
  [NLUMO](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO_CUBES.NLUMO): choose how many occupied and unoccupied
  orbitals are written.
- [HOMO_LIST](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO_CUBES.HOMO_LIST): request specific occupied
  orbitals and override `NHOMO`.
- [STRIDE](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO_CUBES.STRIDE): reduce the grid resolution to make
  files smaller.
- [MAX_FILE_SIZE_MB](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO_CUBES.MAX_FILE_SIZE_MB): let CP2K choose a
  suitable stride to limit cube-file size.
- [WRITE_CUBE](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO_CUBES.WRITE_CUBE): compute eigenvalues while
  suppressing cube-file writing when set to false.

Cube files can become very large, especially for large cells or when many orbitals are requested.
For large grid-based orbital data, [&MO_OPENPMD](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO_OPENPMD) may be
more appropriate.

## MO openPMD file

[&MO_OPENPMD](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO_OPENPMD) writes selected molecular orbitals on the
real-space grid in the openPMD format. Similar to
[&MO_CUBES](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO_CUBES), it is mainly intended for grid-based
visualization or post-processing of selected orbitals rather than for exporting AO-basis MO
coefficients. Compared with plain cube files, openPMD is better suited for structured and
potentially large data, supports parallel I/O through the openPMD backend, and can store the data in
formats such as ADIOS2 BP files depending on the openPMD configuration.

```text
&FORCE_EVAL
  &DFT
    &PRINT
      &MO_OPENPMD
        NHOMO 1
        NLUMO 5
        STRIDE 2 2 2
      &END MO_OPENPMD
    &END PRINT
  &END DFT
&END FORCE_EVAL
```

The most relevant options are [NHOMO](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO_OPENPMD.NHOMO),
[NLUMO](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO_OPENPMD.NLUMO), and
[HOMO_LIST](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO_OPENPMD.HOMO_LIST), which select the orbitals to be
written. [STRIDE](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO_OPENPMD.STRIDE) reduces the grid resolution
and therefore the output size. More advanced openPMD runtime options can be passed through
[OPENPMD_CFG](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO_OPENPMD.OPENPMD_CFG) or
[OPENPMD_CFG_FILE](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO_OPENPMD.OPENPMD_CFG_FILE).

```{important}
Printing openPMD file requires CP2K to be built with
[openPMD](../../technologies/libraries.md#openpmd-structured-output) support.
```

## Choosing an output format

Use [&MO_MOLDEN](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO_MOLDEN) when the calculation is Gamma-only and
the goal is visualization in common molecular-orbital viewers. Use
[&MO_KP](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO_KP) when the calculation uses k-points and a
post-processing tool needs k-point-resolved coefficients. Use
[&MO_CUBES](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO_CUBES) for real-space visualization of a small
number of selected orbitals. Use [&MO](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO) for compact textual
inspection and debugging.
