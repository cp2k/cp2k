# Wannier90 interface

CP2K can generate the input and matrix files required by
[Wannier90](https://wannier.org/) through
[&WANNIER90](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.WANNIER90). The interface is experimental.
It prepares the k-point mesh, eigenvalues, and overlap matrices from a periodic Quickstep
calculation; the subsequent construction and use of Wannier functions are performed by Wannier90.

Wannier90 requires a complete, uniformly weighted k-point mesh with its nearest-neighbour
connectivity. A high-symmetry band path is not a suitable input mesh. See [](../dft/k-points) for
k-point sampling and convergence.

## Basic workflow

1. Run a periodic, diagonalization-based SCF calculation with a converged k-point mesh and enough
   bands for the intended Wannierization.
2. Enable `&DFT%PRINT%WANNIER90`. CP2K writes the Wannier90 input and data files.
3. Add the Wannier90 settings specific to the calculation, such as projections, disentanglement
   windows, or post-processing options, to the generated `.win` file, then run Wannier90 with the
   same seed name.

For example, the following exports a four-function Wannierization using the k-point mesh already
used by the SCF calculation:

```text
&FORCE_EVAL
  &DFT
    &KPOINTS
      SCHEME MONKHORST-PACK 6 6 6
    &END KPOINTS
    &PRINT
      &WANNIER90
        SEED_NAME silicon
        KPOINTS_SOURCE SCF
        WANNIER_FUNCTIONS 4
        ADDED_MOS 4
      &END WANNIER90
    &END PRINT
  &END DFT
&END FORCE_EVAL
```

[WANNIER_FUNCTIONS](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.WANNIER90.WANNIER_FUNCTIONS) sets the number
of Wannier functions. [ADDED_MOS](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.WANNIER90.ADDED_MOS) provides
additional bands for the export, and
[EXCLUDE_BANDS](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.WANNIER90.EXCLUDE_BANDS) can remove selected bands
from it. Choose the exported band window and the subsequent Wannier90 settings for the particular
material and target property.

## Generated files

With `SEED_NAME silicon`, CP2K writes the following Wannier90 files:

- `silicon.win`: a starting Wannier90 input file containing the cell, atomic positions, exported
  band count, and k-point mesh;
- `silicon.mmn`: overlap matrices between neighbouring k-points;
- `silicon.eig`: eigenvalues for the exported bands; and
- `silicon.amn`: an identity projection matrix, only when `USE_BLOCH_PHASES T` is used.

CP2K regenerates these files when the calculation is run. Preserve a separate copy of a completed
Wannier90 input file, or add project-specific settings after the CP2K export has finished.

## Selecting the k-point source

[KPOINTS_SOURCE](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.WANNIER90.KPOINTS_SOURCE) selects how CP2K builds
the export mesh.

### Use the SCF mesh

`KPOINTS_SOURCE SCF` uses the k-point mesh from `&DFT%KPOINTS`. It supports explicit Gamma,
Monkhorst--Pack, MacDonald, and equally weighted `GENERAL` meshes, and is the preferred choice when
the Wannier90 export should match the SCF calculation:

```text
&DFT
  &KPOINTS
    SCHEME MONKHORST-PACK 6 6 6
  &END KPOINTS
  &PRINT
    &WANNIER90
      KPOINTS_SOURCE SCF
      ...
    &END WANNIER90
  &END PRINT
&END DFT
```

`KPOINTS_SOURCE SCF` requires an active `&DFT%KPOINTS` section. For an explicit `GENERAL` mesh, the
points must have equal weights and must form a complete mesh from which Wannier90 connectivity can
be constructed. If CP2K cannot infer the mesh dimensions, set
[MP_GRID](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.WANNIER90.MP_GRID) explicitly.

When the SCF calculation used k-point symmetry reduction, CP2K regenerates the corresponding full
mesh for the export. Wannier90 needs that full mesh even though the SCF calculation solved only its
irreducible subset.

### Use a separate Monkhorst--Pack mesh

`KPOINTS_SOURCE MP_GRID` is the historical default. It builds a full Monkhorst--Pack mesh from
[MP_GRID](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.WANNIER90.MP_GRID), independently of the SCF k-point
setup:

```text
&PRINT
  &WANNIER90
    KPOINTS_SOURCE MP_GRID
    MP_GRID 6 6 6
    ...
  &END WANNIER90
&END PRINT
```

This path performs the necessary full-mesh diagonalizations for the Wannier90 files. It can be used
when the SCF calculation is Gamma-only or when the export mesh intentionally differs from the SCF
mesh, but using a separately chosen mesh requires its own convergence assessment.

## Reusing SCF orbitals

With `KPOINTS_SOURCE SCF`,
[REUSE_SCF_MOS](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.WANNIER90.REUSE_SCF_MOS) is enabled by default.
CP2K reuses the SCF orbital coefficients directly when the SCF mesh is already complete. It can also
reconstruct some time-reversal and atomic-symmetry-related points from a symmetry-reduced SCF mesh.

The reconstruction is only used when the exported band window is suitable. In particular, a
symmetry-reconstructed window must contain complete degenerate subspaces. If it cuts through a
degenerate subspace, CP2K falls back to a full-mesh diagonalization before writing the Wannier90
files. This fallback is intentional and preserves a well-defined exported subspace.

[VALIDATE_REUSE_SCF_MOS](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.WANNIER90.VALIDATE_REUSE_SCF_MOS) builds a
full-mesh reference and compares it with the reconstructed orbitals. It is expensive and intended
for development and diagnostic use, not routine production calculations.

## Bloch phases and projections

[USE_BLOCH_PHASES](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.WANNIER90.USE_BLOCH_PHASES) applies the CP2K
Bloch-phase gauge and writes an identity `.amn` projection file. It is valid only when
`WANNIER_FUNCTIONS` equals the number of exported bands. Disentanglement calculations, or any case
with fewer Wannier functions than exported bands, still require explicit Wannier90 projections.

## Limitations

The CP2K Wannier90 interface is experimental. In particular:

- use a complete k-point mesh rather than a band path;
- verify the convergence of the SCF and export meshes for the target quantity;
- inspect CP2K output when exporting from a symmetry-reduced SCF mesh, since CP2K may reconstruct
  the missing orbitals or fall back to full-mesh diagonalization; and
- consult the Wannier90 documentation for localization, projection, disentanglement, interpolation,
  and post-processing settings that are not controlled by CP2K.
