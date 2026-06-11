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
| TREXIO file         | [&TREXIO](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.TREXIO)         | Export orbital information in TREXIO HDF5 format                         | Useful for interoperability with quantum-chemistry post-processing tools. Requires TREXIO support in the CP2K build.            |
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
        ENERGIES T
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

- [ENERGIES](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO.ENERGIES),
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

### Interpreting MO energies and energy gaps

The MO energies printed by CP2K are one-particle eigenvalues of the chosen electronic-structure
method. They are useful for comparing orbital orderings, occupations, frontier orbitals, and
k-point-resolved band structures. However, they should not in general be interpreted as directly
measurable quantities. Experimental observables such as ionization potentials, electron affinities,
redox potentials, optical excitation energies, and quasiparticle gaps involve total-energy
differences, many-body effects, structural or environmental relaxation, or other physical
contributions that are not captured by a simple one-particle eigenvalue alone.

For finite systems, the HOMO--LUMO gap is the difference between the lowest unoccupied and highest
occupied MO eigenvalues. For periodic k-point calculations, the relevant band gap is obtained by
comparing the valence-band maximum and conduction-band minimum over all k-points and spin channels;
the minimum gap may be indirect. In metallic or smeared calculations, occupations around the Fermi
level can make a HOMO--LUMO or band-gap interpretation ambiguous.

The numerical value of an eigenvalue gap depends on the exchange-correlation functional, basis set,
pseudopotential, spin treatment, k-point sampling, and the number and quality of available
unoccupied states.

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

A shortened Molden file written with `WRITE_CELL` and `WRITE_PSEUDO` looks like:

```text
[Molden Format]
[Cell] Angs
    3.574965       0.000000       0.000000
    0.000000       3.574965       0.000000
    0.000000       0.000000       3.574965
[Atoms] Angs
C      1     6       0.000000      -0.000000       0.000000
C      2     6       2.681237       0.893728       2.681237
...
[Pseudo]
C      1     4
C      2     4
...
[GTO]
       1       0
                        s       5    1.00
                                                        5.605330752        0.111532882
                                                        2.113016391        0.153142194
...
[MO]
Ene=   -3.5996287670E-01
Spin= Alpha
Occup=   2.0000000
     1  2.155619005E-01
     2  1.904995482E-03
...
```

Here `[Cell]` and `[Pseudo]` are CP2K extensions to the original Molden format. The `[GTO]` block
contains the Gaussian basis information, and the `[MO]` block contains the orbital energy, spin,
occupation, and AO coefficients of each printed orbital.

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

The two AO export modes target different post-processing workflows:

- `AO_EXPORT_TYPE GTO_BASIS` writes a compact GTO basis definition. The post-processing code
  reconstructs the overlap matrix, including its k-point dependence.
- `AO_EXPORT_TYPE OVERLAP_MATRIX` writes the explicit overlap matrix $S(k)$ for each k-point. This
  produces larger files but avoids reconstructing the AO overlap externally.

A `.mokp` file contains:

- cell vectors, atomic positions, nuclear charges, effective nuclear charges, and AO ranges;
- the k-point list and k-point weights;
- eigenvalues and occupations for each k-point and spin channel;
- real and, when needed, imaginary parts of the MO coefficient matrix;
- either a GTO basis definition or explicit overlap matrices, controlled by
  [AO_EXPORT_TYPE](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO_KP.AO_EXPORT_TYPE).

A shortened `.mokp` file looks like:

```text
# CP2K_KPOINT_MO_DUMP, Version 2.0
# All energies in Hartree, lengths in Angstrom, k-points in fractional reciprocal coords
# DIMENSIONS: natom nspins nao nkp =       8       1     104      14
# NMO =      20
# USE_REAL_WFN = F
# AO_EXPORT_TYPE = GTO_BASIS
# SPARSE_THRESHOLD = 1.00000E-09
# CELL_VECTORS [Angstrom]
    3.574965       0.000000       0.000000
    0.000000       3.574965       0.000000
    0.000000       0.000000       3.574965
# ATOM_LIST: Atom_ID  Element  Z_nuc  Z_eff  X [Ang]  Y [Ang]  Z [Ang]  First_AO  Last_AO
     1  C      6     4       0.000000      -0.000000       0.000000         1        13
...
# KPOINT_LIST: ikp  kx  ky  kz  weight
     1   -3.33333333E-01   -3.33333333E-01   -3.33333333E-01    7.40740741E-02
...
# GTO_BASIS (spherical, MOLDEN convention for contraction coefficients)
     1   0
s      5    1.00
    5.605330752E+00    1.115328819E-01
...
# BEGIN_KPOINT_SPIN ikp ispin =     1     1
# EIGENVALUES
   -2.634414483E-01   -1.634039635E-01   -1.633988780E-01   -1.633988754E-01
...
# OCCUPATIONS
    2.000000000E+00    2.000000000E+00    2.000000000E+00    2.000000000E+00
...
# MO_COEFF_RE (Sparse: mo_index  ao_index  value)
     1     1     1.152519336E-01
     1     2     6.519198415E-03
...
# MO_COEFF_IM (Sparse: mo_index  ao_index  value)
     1     1    -1.607824559E-01
     1     2    -8.895941146E-03
...
```

```{note}
`&MO_KP` is only active for k-point calculations. For Gamma-only calculations, use `&MO_MOLDEN`,
`&MO`, or `&MO_CUBES` depending on the intended use.
```

## TREXIO HDF5 file

```{important}
Printing TREXIO file requires CP2K to be built with
[TREXIO](../../technologies/libraries.md#trexio-unified-computational-chemistry-format) and
[HDF5](../../technologies/libraries.md#hdf5) support.
```

[&TREXIO](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.TREXIO) writes orbital information in the
[TREXIO](https://trex-coe.github.io/trexio/) format. CP2K writes the HDF5 backend and appends the
`.h5` suffix to the requested filename. If
[FILENAME](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.TREXIO.FILENAME) is not set explicitly, the output file
is named `<PROJECT>-TREXIO.h5`.

```text
&FORCE_EVAL
  &DFT
    &PRINT
      &TREXIO
        FILENAME orbitals
        CARTESIAN F
      &END TREXIO
    &END PRINT
  &END DFT
&END FORCE_EVAL
```

Useful controls include:

- [FILENAME](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.TREXIO.FILENAME): choose the body of the `.h5` output
  filename.
- [CARTESIAN](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.TREXIO.CARTESIAN): store the MOs in the Cartesian
  basis instead of the default spherical basis.
- [FULL_KPOINT_GRID](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.TREXIO.FULL_KPOINT_GRID): for symmetry-reduced
  k-point calculations, export MOs on the full unreduced k-point grid instead of only the
  irreducible grid.
- [REUSE_SCF_MOS](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.TREXIO.REUSE_SCF_MOS): when `FULL_KPOINT_GRID` is
  active, try to reconstruct the full-grid MO coefficients from the symmetry-reduced SCF orbitals
  before falling back to a full-grid diagonalization.

The resulting file, here `orbitals.h5`, stores the information in a structured binary format rather
than as a plain text file. It contains system metadata, nuclear coordinates and effective charges,
cell and periodic-boundary information, electron numbers, basis and AO metadata, and canonical MO
information such as MO coefficients, energies, occupations, and spin labels. For k-point
calculations, CP2K also writes the k-point information and, when needed, the imaginary part of the
MO coefficients.

A TREXIO file is binary, so the useful "file output example" is its HDF5 layout rather than raw text
content. For example, inspecting the file with an HDF5 tool shows groups and datasets such as:

```text
$ h5ls orbitals.h5
metadata                 Group
nucleus                  Group
cell                     Group
pbc                      Group
electron                 Group
state                    Group
basis                    Group
ao                       Group
ecp                      Group
mo                       Group

$ h5ls orbitals.h5/mo
coefficient              Dataset {ao_num, mo_num}
energy                   Dataset {mo_num}
occupation               Dataset {mo_num}
spin                     Dataset {mo_num}
type                     Dataset {...}
...
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

Each selected orbital is written to a separate cube file. The file name contains the MO index and
spin channel, for example `PROJECT-WFN_00016_1-1_0.cube` for MO 16, spin 1. A shortened cube file
looks like:

```text
-Quickstep-
 WAVEFUNCTION    16 spin 1 i.e. HOMO - 0
    8    0.000000    0.000000    0.000000
   36    0.187500    0.000000    0.000000
   36    0.000000    0.187500    0.000000
   36    0.000000    0.000000    0.187500
    6    0.000000    0.000000   -0.000000    0.000000
    6    0.000000    2.681237    0.893728    2.681237
...
  1.23456E-07 -2.34567E-07  3.45678E-07 -4.56789E-07  5.67890E-07 -6.78901E-07
...
```

The first two lines are comments. The following lines define the grid origin, grid dimensions, and
grid vectors, followed by the atom list and then the real-space values of the selected orbital.

```{note}
CP2K writes the real-space orbital values in fixed-width `E13.5` fields following the
[Gaussian cube file convention](https://gaussian.com/cubegen/). Some downstream cube-file parsers
assume this fixed-width layout rather than only whitespace-separated floating-point values, so the
formatting should be preserved when cube files are post-processed or rewritten.
```

Cube files can become very large, especially for large cells or when many orbitals are requested.
For large grid-based orbital data, [&MO_OPENPMD](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO_OPENPMD) may be
more appropriate.

## MO openPMD file

```{important}
Printing openPMD file requires CP2K to be built with
[openPMD](../../technologies/libraries.md#openpmd-structured-output) support.
```

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

## Choosing an output format

- Use [&MO_MOLDEN](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO_MOLDEN) when the calculation is Gamma-only
  and the goal is visualization in common molecular-orbital viewers.
- Use [&MO_KP](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO_KP) when the calculation uses k-points and a
  post-processing tool needs k-point-resolved coefficients.
- Use [&TREXIO](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.TREXIO) when a structured HDF5 orbital file is
  needed for interoperability with external quantum-chemistry tools.
- Use [&MO_CUBES](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO_CUBES) or
  [&MO_OPENPMD](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO_OPENPMD) for real-space visualization of a small
  number of selected orbitals.
- Use [&MO](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.MO) for compact textual inspection and debugging.

```{note}
The overall phase of an individual molecular orbital is arbitrary. Therefore, a sign flip of a
single MO, or a phase rotation in a complex k-point calculation, is not by itself a physical
difference. Post-processing tools and regression tests should compare phase-invariant quantities
when appropriate.
```
