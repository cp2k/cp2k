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
The filename prefix and CC4S-specific keywords are shown with their default values:

```text
&CC4S ON
  FILENAME cc4s
  METHOD CCSD
  N_FROZEN_CORE 0
  MAX_ITER 100
  ENERGY_CONVERGENCE 1e-7
  AMPLITUDES_CONVERGENCE 1e-6
  OUTPUT_BUFFER_SIZE 8
  MAX_MEMORY 1024
&END CC4S
```

| Keyword                                                                                | Meaning                                                                                                                                                                                 |
| -------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [FILENAME](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.CC4S.FILENAME)                             | Common prefix for the generated files, using standard CP2K print-key naming.                                                                                                            |
| [METHOD](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.CC4S.METHOD)                                 | `MP2` generates MP2; `CCSD` adds CCSD; `CCSD_T` also adds perturbative triples.                                                                                                         |
| [N_FROZEN_CORE](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.CC4S.N_FROZEN_CORE)                   | Number of lowest occupied spatial orbitals excluded from correlation; must leave at least one occupied orbital. With pseudopotentials, excludes additional explicitly treated orbitals. |
| [MAX_ITER](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.CC4S.MAX_ITER)                             | Maximum CC4S coupled-cluster iterations.                                                                                                                                                |
| [ENERGY_CONVERGENCE](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.CC4S.ENERGY_CONVERGENCE)         | CC4S energy convergence threshold, in Hartree.                                                                                                                                          |
| [AMPLITUDES_CONVERGENCE](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.CC4S.AMPLITUDES_CONVERGENCE) | CC4S amplitude convergence threshold.                                                                                                                                                   |
| [OUTPUT_BUFFER_SIZE](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.CC4S.OUTPUT_BUFFER_SIZE)         | Maximum output buffer per MPI rank, in MiB (1 MiB = 1048576 bytes). Positive real value; defaults to 8.                                                                                 |
| [MAX_MEMORY](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.CC4S.MAX_MEMORY)                         | Estimated export workspace limit, in MB per MPI rank; includes output and MPI reduction buffers, but excludes existing SCF storage and CC4S memory.                                     |

`OUTPUT_BUFFER_SIZE` is rounded down to whole 8-byte values and must hold at least one value. Its
element count must fit a default `INTEGER` (at most 2,147,483,647 with 32-bit integers). The actual
buffer is capped by the exported vertex size and printed in MiB. This keyword controls output
gathering and writing; it does not change transformation batches or numerical thresholds.
`MAX_MEMORY` remains in decimal MB (1 MB = 1,000,000 bytes).

Export is OFF by default. `ON` enables it independently of `GLOBAL/PRINT_LEVEL`; standard print
levels are also accepted. The iteration and convergence settings above control CC4S, independently
of CP2K's SCF and integral thresholds. Generic `EACH` and `ADD_LAST` controls need not be specified;
they select output only when the post-SCF exporter is called.

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

With `PROJECT water` and the defaults above, the files in the CP2K working directory are:

| File                                | Contents                                                                           |
| ----------------------------------- | ---------------------------------------------------------------------------------- |
| `water-cc4s-EigenEnergies.yaml`     | Tensor header, Fermi level and complete orbital-energy list in `metaData.energies` |
| `water-cc4s-EigenEnergies.elements` | Plain-text orbital energies in Hartree                                             |
| `water-cc4s-CoulombVertex.yaml`     | Tensor header, dimensions and interaction metadata                                 |
| `water-cc4s-CoulombVertex.elements` | Raw binary real 64-bit vertex elements, without Fortran record markers             |
| `water-cc4s.in`                     | CC4S algorithm sequence and solver settings                                        |

With `PROJECT water`:

- `FILENAME sample` produces `water-sample-EigenEnergies.yaml`, `water-sample-CoulombVertex.yaml`
  and `water-sample.in`, with matching `.elements` files.
- `FILENAME ./sample` produces `sample-EigenEnergies.yaml`, `sample-CoulombVertex.yaml` and
  `sample.in`, with matching `.elements` files and no project prefix.

Paths may refer to existing directories; CP2K does not create them. Exact-name (`=name`) and stdout
forms are unsupported because this export writes multiple distinct files. The generated input
references the matching tensor filenames; run CC4S from the directory containing these files.
Standard iteration suffixes, if selected, apply to the whole file set.

The vertex has index order `(auxiliary, orbital, orbital)`, with the auxiliary index varying
fastest. Its units are the square root of Hartree. The headers use tensor format version `100`,
`Real64` and `unit: 1.0`. No Fock-matrix file is written. See the CC4S
[CoulombVertex specification](https://manuals.cc4s.org/objects/CoulombVertex.html) for the tensor
convention; CP2K supplies Gaussian fitting factors rather than the plane-wave construction discussed
there.
