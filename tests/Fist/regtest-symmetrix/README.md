# symmetrix regression test

This test evaluates the [symmetrix](https://github.com/wcwitt/symmetrix) MACE backend on two rigid
water molecules and checks the total energy against a stored reference.

It only runs when CP2K is built with `-DCP2K_USE_SYMMETRIX=ON` (it is gated on the `symmetrix`
cp2kflag in `tests/TEST_DIRS`), so a standard build skips it.

## The model

`data/Symmetrix/toy_mace-1-8.json` is a tiny untrained MACE model (H and O, 2 channels, seeded
weights) exported with symmetrix's `extract_mace` tool. It carries no physics; it exists so the
backend can be exercised without a large download. It was produced with `mace-torch` 0.3.16 as
follows:

```python
import numpy as np
import torch
from e3nn import o3
from mace import modules

torch.manual_seed(42)
torch.set_default_dtype(torch.float64)

model = modules.ScaleShiftMACE(
    r_max=4.0,
    num_bessel=4,
    num_polynomial_cutoff=5,
    max_ell=3,
    interaction_cls=modules.interaction_classes["RealAgnosticResidualInteractionBlock"],
    interaction_cls_first=modules.interaction_classes["RealAgnosticInteractionBlock"],
    num_interactions=2,
    num_elements=2,
    hidden_irreps=o3.Irreps("2x0e"),
    MLP_irreps=o3.Irreps("16x0e"),
    atomic_energies=np.array([-0.5, -75.0]),
    avg_num_neighbors=8.0,
    atomic_numbers=[1, 8],
    correlation=3,
    gate=torch.nn.functional.silu,
    atomic_inter_scale=1.0,
    atomic_inter_shift=0.0,
)
torch.save(model.double(), "toy_mace.model")
```

followed by an export at 128 spline points, written as compact JSON to keep the file below 300 kB.
`make_toy_model.py` in this directory does both steps and is the definitive recipe:

```
python make_toy_model.py -o ../../../data/Symmetrix/toy_mace-1-8.json
```

Note this deliberately does not use `python -m symmetrix.cli.extract_mace`: that CLI has no option
for the spline resolution (it always uses its own default) and writes indented JSON, so it cannot
reproduce the committed file. Regenerating with different `mace-torch`/`torch` versions may shift
the weights, which would change the reference values in `TEST_FILES.toml`; the two have to be
updated together.

Trained models are exported the same way; see `docs/methods/machine_learning/symmetrix.md`.

`sym_h2o.inp` resolves `Symmetrix/toy_mace-1-8.json` relative to `CP2K_DATA_DIR` (the `data/`
directory by default).
