# Methods

```{toctree}
---
titlesonly:
maxdepth: 2
---
dft/index
post_hartree_fock/index
semiempiricals/index
machine_learning/index
embedding/index
qm_mm/index
sampling/index
optimization/index
properties/index
```

CP2K combines several levels of molecular and materials simulation in one input framework. The
method sections below cover electronic-structure methods, machine-learned and semi-empirical models,
embedding and QM/MM workflows, sampling, optimization, and property calculations.

For most new atomistic simulations, start with Quickstep density functional theory. Choose a
consistent basis-set and potential pair, converge the real-space grid and SCF settings, and then add
the property, sampling, or higher-level method needed for the scientific question. More specialized
methods such as MP2, RPA, GW/BSE, QM/MM, and machine-learning potentials often reuse the same
structural ideas but add method-specific basis sets, auxiliary spaces, coupling terms, or
post-processing steps.

## Selected Videos

The videos below provide CP2K-focused introductions to the main method families and workflows.

```{youtube} 3Cw4h3MrZ8k
---
align: center
privacy_mode:
---
```

```{youtube} wyux20qVlck
---
align: center
privacy_mode:
---
```

```{youtube} zSt8KQ2Hf3c
---
align: center
privacy_mode:
---
```
