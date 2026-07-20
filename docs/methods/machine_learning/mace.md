# MACE

MACE is a framework for building interatomic potentials using higher-order equivariant
message-passing neural networks. The methodology is described in detail in the literature by Batatia
et al. (2022).

**Note:** Running MACE requires a CP2K build with LibTorch support

Like the [NequIP/Allegro](https://manual.cp2k.org/trunk/methods/machine_learning/nequip.html)
interface, CP2K runs MACE models through the generic LibTorch interface: a trained MACE model is
exported once to a self-contained TorchScript file (`.pth`), which CP2K then loads and evaluates at
runtime. No Python interpreter is involved during the simulation.

## Exporting a MACE model for CP2K

Wraps a trained MACE model (`.model`) and compiles it to CP2K-loadable TorchScrip file with the
helper script `cp2k/tools/mace/create_cp2k_model.py`:

```shell
python create_cp2k_model.py my_mace.model --dtype float64
# -> writes my_mace.model-cp2k.pth
```

This conversion only needs to be performed once on a machine with both `torch` and `mace` installed.
The resulting `*.pth` file uses the same tensor and metadata format as the NequIP interface (inputs
`pos`, `edge_index`, `edge_cell_shift`, `cell`, `atom_types`; outputs `atomic_energy`, `forces`,
`virial`) and embeds the metadata (`num_types`, `r_max`, `type_names`, `model_dtype`) that CP2K
reads to build the neighbour graph. The `torch` version used for export must be compatible with the
LibTorch version linked into CP2K.

## Input Section

Inference is configured through the [MACE](#CP2K_INPUT.FORCE_EVAL.MM.FORCEFIELD.NONBONDED.MACE)
section within the `&NONBONDED` forcefield parameters:

```text
&FORCEFIELD
  &NONBONDED
    &MACE
      ATOMS Cu
      POT_FILE_NAME MACE/my_mace.model-cp2k.pth
    &END MACE
  &END NONBONDED
&END FORCEFIELD
```

- [ATOMS](#CP2K_INPUT.FORCE_EVAL.MM.FORCEFIELD.NONBONDED.MACE.ATOMS): a list of elements/kinds; the
  mapping to the model type list must be consistent with the coordinates in `&COORDS`/`&TOPOLOGY`.
- [POT_FILE_NAME](#CP2K_INPUT.FORCE_EVAL.MM.FORCEFIELD.NONBONDED.MACE.POT_FILE_NAME): path to the
  exported MACE model.

MACE is a message-passing model with a non-local receptive field. As with NequIP, the interface
evaluates the full system on every MPI rank and divides the energy, forces, and virial by the number
of ranks.

## Further Resources

- **MACE:** Paper [](#Batatia2022) and source code at
  [github.com/ACEsuit/mace](https://github.com/ACEsuit/mace).
- **e3nn:** For an introduction to Euclidean neural networks, visit [e3nn.org](https://e3nn.org).
