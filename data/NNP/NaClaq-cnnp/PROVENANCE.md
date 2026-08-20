# NaCl(aq) RPA C-NNP, one committee member

Source: https://github.com/niamhon/nacl-water Commit: 459f90fdce1e68647f1acfa9ef43f4e8b00e69e3 Path
in source: `models/RPA/final_model/train_001/`

License: CC-BY-SA 4.0. The license notice and its canonical URIs travel with the model in `LICENSE`
beside this file, as that license requires. This is data, distributed alongside CP2K rather than
combined with it, so it does not affect the GPL-2.0-or-later terms of the source tree.

Paper: N. O'Neill, B. X. Shi, K. Fong, A. Michaelides, C. Schran, "To Pair or not to Pair?
Machine-Learned Explicitly-Correlated Electronic Structure for NaCl in Water", *J. Phys. Chem.
Lett.* **15**, 6081 (2024). DOI: 10.1021/acs.jpclett.4c01030 arXiv: 2311.01527

Architecture (from input.nn):

- 4 elements: O, H, Na, Cl
- 2 hidden layers, 25 nodes each, tanh-tanh-linear
- cutoff_type 2 (tanh cubed), 12 Bohr (~6.35 Å) symmetry-function cutoff
- scale_symmetry_functions + center_symmetry_functions

## Why this model is in the repository

The published model is an 8-member committee. Only `train_001` is committed here, because the
regression test that needs it, `tests/NNP/regtest-1/NaCl-dilute-hetero.inp`, depends on the
symmetry-function layout rather than on the committee: it is the smallest available model with four
elements, two of which sit beyond the symmetry-function cutoff from each other in a dilute box. That
gives heterogeneous angular groups whose second neighbour list is empty, which no water-only model
can reach, and which is the case the test asserts is skipped consistently by both the descriptor and
the force-scatter stage.

Files (under `train_001/`):

```
weights.001.data   H weights (Z=1)
weights.008.data   O weights (Z=8)
weights.011.data   Na weights (Z=11)
weights.017.data   Cl weights (Z=17)
scaling.data
input.nn
```
