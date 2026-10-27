# CC4S interface

CP2K can export a converged restricted Hartree–Fock reference for MP2, CCSD and CCSD(T) calculations
with [CC4S](https://github.com/cc4s/cc4s). CP2K writes canonical orbital energies, a Gaussian RI
Coulomb vertex and a CC4S input file. Run CC4S separately to calculate correlation energies. CP2K
does not need to be linked against CC4S.

The interface has been tested with **CC4S v1.1, revision `8204fb21`**. The underlying methods are
described in the [RI-HFX paper](#Bussy2023) and the
[auxiliary-field representation paper](#Hummel2017).

## Reference calculation

Use [GLOBAL/RUN_TYPE ENERGY](#CP2K_INPUT.GLOBAL.RUN_TYPE) and a converged, restricted closed-shell
HF calculation with integer occupations, at least one occupied and one virtual orbital, and a
positive HOMO–LUMO gap. Set `XC_FUNCTIONAL NONE`, use exactly one `HF` section with `FRACTION 1`,
and select [HF/RI/RI_FLAVOR MO](#CP2K_INPUT.FORCE_EVAL.DFT.XC.HF.RI.RI_FLAVOR). Supply a suitable
Gaussian fitting basis through `BASIS_SET RI_HFX` for each atomic kind. These settings belong to
RI-HFX; the `WF_CORRELATION` integral settings do not configure this exporter. ADMM and smearing are
unsupported. `ENERGY_FORCE` is also accepted, but any forces are HF forces.

## Enable export

Add the following inside [FORCE_EVAL/DFT/PRINT](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT) to enable export.
All six optional keywords are shown with their default values:

```text
&CC4S ON
  METHOD CCSD
  N_FROZEN_CORE 0
  MAX_ITER 100
  ENERGY_CONVERGENCE 1e-7
  AMPLITUDES_CONVERGENCE 1e-6
  MAX_MEMORY 1024
&END CC4S
```

| Keyword                                                                                | Meaning                                                                                                                                                                                 |
| -------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [METHOD](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.CC4S.METHOD)                                 | `MP2` generates MP2; `CCSD` adds CCSD; `CCSD_T` also adds perturbative triples.                                                                                                         |
| [N_FROZEN_CORE](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.CC4S.N_FROZEN_CORE)                   | Number of lowest occupied spatial orbitals excluded from correlation; must leave at least one occupied orbital. With pseudopotentials, excludes additional explicitly treated orbitals. |
| [MAX_ITER](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.CC4S.MAX_ITER)                             | Maximum CC4S coupled-cluster iterations.                                                                                                                                                |
| [ENERGY_CONVERGENCE](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.CC4S.ENERGY_CONVERGENCE)         | CC4S energy convergence threshold, in Hartree.                                                                                                                                          |
| [AMPLITUDES_CONVERGENCE](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.CC4S.AMPLITUDES_CONVERGENCE) | CC4S amplitude convergence threshold.                                                                                                                                                   |
| [MAX_MEMORY](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.CC4S.MAX_MEMORY)                         | Estimated export workspace limit, in MB per MPI rank; excludes existing SCF storage and CC4S memory.                                                                                    |

The section is OFF by default and independent of `GLOBAL/PRINT_LEVEL`. The iteration and convergence
settings control CC4S, independently of CP2K's SCF and integral thresholds.

Complete small input examples are provided in
[`tests/QS/regtest-cc4s`](https://github.com/cp2k/cp2k/tree/master/tests/QS/regtest-cc4s): `h2o.inp`
with its `BASIS_CC4S` file, and `he-gamma.inp`. These are regression examples, not numerically
converged production settings.

## Boundary conditions and interactions

| System                                 | CELL and POISSON periodicity | HF potential | RI metric |
| -------------------------------------- | ---------------------------- | ------------ | --------- |
| Isolated molecule                      | `NONE` in both               | `COULOMB`    | `COULOMB` |
| Neutral three-dimensional Γ-point cell | `XYZ` in both                | `TRUNCATED`  | `HFX`     |

Periodic export supports Γ-point calculations only. `CUTOFF_RADIUS` must be positive and smaller
than half the minimum perpendicular spacing between opposite cell faces.

## Exported files

Export occurs once after converged SCF in a single-point calculation. Geometry optimization and MD
are unsupported.

Each export writes the following files in the CP2K working directory:

| File                     | Contents                                                                           |
| ------------------------ | ---------------------------------------------------------------------------------- |
| `EigenEnergies.yaml`     | Tensor header, Fermi level and complete orbital-energy list in `metaData.energies` |
| `EigenEnergies.elements` | Plain-text orbital energies in Hartree                                             |
| `CoulombVertex.yaml`     | Tensor header, dimensions and interaction metadata                                 |
| `CoulombVertex.elements` | Raw binary real 64-bit vertex elements, without Fortran record markers             |
| `cc4s.in`                | CC4S algorithm sequence and solver settings                                        |

The vertex has index order `(auxiliary, orbital, orbital)`, with the auxiliary index varying
fastest. Its units are the square root of Hartree. The headers use tensor format version `100`,
`Real64` and `unit: 1.0`. No Fock-matrix file is written. See the CC4S
[CoulombVertex specification](https://manuals.cc4s.org/objects/CoulombVertex.html) for the tensor
convention; CP2K supplies Gaussian fitting factors rather than the plane-wave construction discussed
there.
