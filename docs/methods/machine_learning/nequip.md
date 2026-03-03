# NequIP and Allegro

NequIP and Allegro are frameworks designed for developing interatomic potentials for molecular
dynamics simulations using deep equivariant neural networks. The methodology and recent
high-performance upgrades are described in detail in the literature ([](#Batzner2022),
[](#Musaelian2023), [](#Tan2025)).

The CP2K interface is compatible with models trained and compiled with NequIP version >= 0.7.0, and
it is consistent with the LAMMPS `pair_nequip_allegro` (v0.7.0) integration.

## Input Section

Inference in CP2K has been unified and is configured entirely through the
[NEQUIP](#CP2K_INPUT.FORCE_EVAL.MM.FORCEFIELD.NONBONDED.NEQUIP) section within the `&NONBONDED`
forcefield parameters.

An example of the input configuration:

```text
&FORCEFIELD
  &NONBONDED
    &NEQUIP
      MODEL_TYPE  NEQUIP # possible choices are NEQUIP or ALLEGRO
      ATOMS H O
      POT_FILE_NAME NequIP/waterscan-neq0.16.nequip.pth
      UNIT_ENERGY eV
      UNIT_FORCES eV*angstrom^-1
      UNIT_LENGTH angstrom
    &END NEQUIP
  &END NONBONDED
&END FORCEFIELD

```

- [MODEL_TYPE](#CP2K_INPUT.FORCE_EVAL.MM.FORCEFIELD.NONBONDED.NEQUIP.MODEL_TYPE): Specifies the
  architecture of the loaded model (`NEQUIP` or `ALLEGRO`).
- [ATOMS](#CP2K_INPUT.FORCE_EVAL.MM.FORCEFIELD.NONBONDED.NEQUIP.ATOMS): Expects a list of
  elements/kinds.
- [POT_FILE_NAME](#CP2K_INPUT.FORCE_EVAL.MM.FORCEFIELD.NONBONDED.NEQUIP.POT_FILE_NAME): The path to
  the NequIP/Allegro model.
- **`UNIT_*`**: These tags explicitly define the units for the model's internal lengths, energies,
  and forces.

Full example input files demonstrating production-ready molecular dynamics setups can be found in
the CP2K regression tests directory:

- `tests/Fist/regtest-nequip/water-bulk.inp`
- `tests/Fist/regtest-allegro/water-bulk.inp`

## Compiling CP2K with LibTorch

Running NequIP or Allegro requires compiling CP2K with the LibTorch library. Versions of LibTorch
(2.4 through 2.7) are supported.

For CP2K binaries running on CPUs, installing the toolchain using the flag `--with-libtorch` is
sufficient.

To benefit from GPU acceleration, either compile LibTorch from scratch or download the precompiled
LibTorch library for CUDA from [PyTorch](https://pytorch.org) and provide the appropriate path to
the toolchain script:

```shell
./install_cp2k_toolchain.sh --with-libtorch=<path-to-libtorch-cuda>

```

## Validation & Reproducibility

- **Comparison with LAMMPS:** We have verified that this implementation numerically reproduces the
  results of the LAMMPS `pair_nequip_allegro` plugin.
- **Data:** The training datasets, model files inside `data/NequIP` and `data/Allegro`, input
  scripts, and parity plots used for validation are available on Zenodo:
  [doi:10.5281/zenodo.18848354](https://doi.org/10.5281/zenodo.18848354).

## Further Resources

For additional references on NequIP, Allegro, and equivariant neural networks (e3nn), see:

- **High-Performance Upgrades:** Paper [](#Tan2025) and source code at
  [github.com/mir-group/pair_nequip_allegro](https://github.com/mir-group/pair_nequip_allegro).
- **Allegro:** Paper [](#Musaelian2023) and source code at
  [github.com/mir-group/allegro](https://github.com/mir-group/allegro).
- **NequIP:** Paper [](#Batzner2022) and source code at
  [github.com/mir-group/nequip](https://github.com/mir-group/nequip).
- **e3nn:** For an introduction to Euclidean neural networks, visit [e3nn.org](https://e3nn.org) and
  [doi:10.5281/zenodo.7430260](https://dx.doi.org/10.5281/zenodo.7430260).
