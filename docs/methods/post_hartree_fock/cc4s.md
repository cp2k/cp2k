# CC4S interface

CP2K can export data from a converged, restricted closed-shell RI-HF calculation for subsequent MP2,
CCSD and CCSD(T) calculations with [CC4S](https://github.com/cc4s/cc4s). The export contains
canonical orbital energies, a Gaussian RI Coulomb vertex and a CC4S input file. Run CC4S separately
to calculate correlation energies. CP2K does not need to be linked against CC4S.

The interface has been tested with **CC4S v1.1, revision `8204fb21`**. The underlying methods are
described in the [RI-HFX paper](#Bussy2023) and the
[auxiliary-field representation paper](#Hummel2017).

## Reference calculation

Use [GLOBAL/RUN_TYPE ENERGY](#CP2K_INPUT.GLOBAL.RUN_TYPE) and converge a restricted closed-shell HF
reference with integer occupations, at least one occupied and one virtual orbital, and a positive
HOMO–LUMO gap. For an isolated molecule, add the following inside `FORCE_EVAL/DFT`:

```text
&XC
  &XC_FUNCTIONAL NONE
  &END XC_FUNCTIONAL
  &HF
    FRACTION 1.0
    &INTERACTION_POTENTIAL
      POTENTIAL_TYPE COULOMB
    &END INTERACTION_POTENTIAL
    &RI
      RI_FLAVOR MO
      RI_METRIC COULOMB
    &END RI
  &END HF
&END XC
```

Use exactly one `HF` section. Supply an orbital basis and a suitable Gaussian fitting basis through
`BASIS_SET RI_HFX` for each atomic kind. The `RI` settings above configure RI-HFX; the
`WF_CORRELATION` integral settings do not configure this exporter. For a periodic Γ-point cell, use
the interaction potential and RI metric shown [below](#boundary-conditions-and-interactions). ADMM
is not supported; do not include the `AUXILIARY_DENSITY_MATRIX_METHOD` section. Smearing is also
unsupported. `ENERGY_FORCE` is accepted, but any forces are HF forces.

A wavefunction restart can make the final HF calculation cheaper or help it converge. For example,
first converge a GGA calculation or an HF calculation with different SCF settings. Then set
`WFN_RESTART_FILE_NAME` in `DFT` to its saved wavefunction and `SCF_GUESS RESTART` in `DFT/SCF`. Use
the final RI-HF settings above and enable CC4S export only in this final calculation. The restarted
HF calculation must converge before export; it may take just one SCF iteration, but this is not
guaranteed. See [restarting CP2K calculations](../restarting.md) for wavefunction-restart setup.

## Enable export

Add the following [CC4S](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.CC4S) section inside
`FORCE_EVAL/DFT/PRINT`. It shows the standard `FILENAME` control and all CC4S-specific keywords at
their default values:

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

| Keyword                                                                                | Meaning                                                                                                                                                                              |
| -------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| [FILENAME](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.CC4S.FILENAME)                             | Common prefix for the exported files, following CP2K print-key naming.                                                                                                               |
| [METHOD](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.CC4S.METHOD)                                 | `MP2` generates MP2; `CCSD` adds CCSD; `CCSD_T` also adds perturbative triples.                                                                                                      |
| [N_FROZEN_CORE](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.CC4S.N_FROZEN_CORE)                   | Number of lowest occupied spatial orbitals omitted from export; at least one occupied orbital must remain. With a pseudopotential, these are additional explicitly treated orbitals. |
| [MAX_ITER](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.CC4S.MAX_ITER)                             | Maximum number of CC4S coupled-cluster iterations.                                                                                                                                   |
| [ENERGY_CONVERGENCE](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.CC4S.ENERGY_CONVERGENCE)         | CC4S correlation-energy convergence threshold in Hartree.                                                                                                                            |
| [AMPLITUDES_CONVERGENCE](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.CC4S.AMPLITUDES_CONVERGENCE) | CC4S amplitude convergence threshold.                                                                                                                                                |
| [OUTPUT_BUFFER_SIZE](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.CC4S.OUTPUT_BUFFER_SIZE)         | Maximum output buffer per MPI rank, in MiB.                                                                                                                                          |
| [MAX_MEMORY](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.CC4S.MAX_MEMORY)                         | Estimated export workspace budget per MPI rank, in decimal MB.                                                                                                                       |

`OUTPUT_BUFFER_SIZE` controls gathering and writing the vertex. It must be large enough to hold at
least one double-precision floating-point element of the Coulomb vertex; the actual buffer may be
smaller than the requested maximum. `MAX_MEMORY` selects orbital and auxiliary transformation batch
sizes. Its estimate includes output and MPI reduction buffers, but excludes existing SCF storage and
CC4S memory. It is not a limit on total process memory. CP2K stops if the budget cannot accommodate
canonicalization or the smallest transformation batch.

Export is off by default. `ON` enables it independently of `GLOBAL/PRINT_LEVEL`. The CC4S iteration
and convergence settings above do not change CP2K's SCF or integral thresholds. The exporter checks
MPI and tensor-library count limits independently of `MAX_MEMORY`.

Complete small input examples are provided in
[`tests/QS/regtest-cc4s`](https://github.com/cp2k/cp2k/tree/master/tests/QS/regtest-cc4s): `h2o.inp`
with its `BASIS_CC4S` file, and `he-gamma.inp`. These are regression examples, not numerically
converged production settings.

## Boundary conditions and interactions

| System                                 | `CELL` and `POISSON` periodicity | HF potential | RI metric |
| -------------------------------------- | -------------------------------- | ------------ | --------- |
| Isolated molecule                      | `NONE` in both                   | `COULOMB`    | `COULOMB` |
| Neutral three-dimensional Γ-point cell | `XYZ` in both                    | `TRUNCATED`  | `HFX`     |

Periodic export is limited to Γ-only calculations without k-point sampling. `CUTOFF_RADIUS` must be
positive and smaller than half the minimum perpendicular spacing between opposite cell faces.

## Exported files

Export occurs once after converged SCF in a single-point calculation. Geometry optimization and MD
are unsupported.

With `PROJECT water` and the default `FILENAME cc4s`, CP2K writes these files in its working
directory:

| File                                | Contents                                                                            |
| ----------------------------------- | ----------------------------------------------------------------------------------- |
| `water-cc4s-EigenEnergies.yaml`     | Tensor header, Fermi energy, HOMO index and orbital energies in `metaData.energies` |
| `water-cc4s-EigenEnergies.elements` | Plain-text orbital energies in Hartree                                              |
| `water-cc4s-CoulombVertex.yaml`     | Tensor header, dimensions and interaction metadata                                  |
| `water-cc4s-CoulombVertex.elements` | Raw binary Real64 vertex elements, without Fortran record markers                   |
| `water-cc4s.in`                     | CC4S algorithm sequence and solver settings                                         |

With `PROJECT water`:

- `FILENAME sample` produces `water-sample-EigenEnergies.yaml`, `water-sample-CoulombVertex.yaml`
  and `water-sample.in`, with matching `.elements` files.
- `FILENAME ./sample` produces `sample-EigenEnergies.yaml`, `sample-CoulombVertex.yaml` and
  `sample.in`, with matching `.elements` files and no project prefix.

`FILENAME` can include an existing directory; CP2K does not create directories for the export.
Exact-name (`=name`) and stdout forms are unsupported because the export writes multiple files. The
generated CC4S input references the matching tensor filenames. Run CC4S from the directory
containing the files. CP2K reports the expected and written binary vertex sizes in its output.

The vertex has index order `(auxiliary, orbital, orbital)`, with the auxiliary index varying
fastest. Its units are the square root of Hartree. The headers use tensor format version `100`,
`Real64` and `unit: 1.0`. No Fock-matrix file is written. See the CC4S
[CoulombVertex specification](https://manuals.cc4s.org/objects/CoulombVertex.html) for the tensor
convention; CP2K supplies Gaussian fitting factors rather than the plane-wave construction discussed
there.
