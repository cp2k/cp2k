# Density Functional Theory

```{toctree}
---
titlesonly:
maxdepth: 2
---
gpw
gapw
hartree-fock/index
local_ri
constrained
cneo
linear_scaling
k-points
basis_sets
pseudopotentials
cutoff
```

Density functional theory in CP2K is primarily provided by the Quickstep module. Most production
calculations use the Gaussian and plane waves (GPW) method with Gaussian basis sets,
pseudopotentials, and real-space grids for densities and potentials. The Gaussian augmented plane
waves (GAPW) method extends the same framework to all-electron and more core-sensitive calculations.

For new inputs, first choose a consistent basis-set and potential pair, then converge the MGRID
cutoffs and the SCF settings for the target property. The pages in this section collect the main
Quickstep ingredients: GPW/GAPW, hybrid functionals and ADMM, local RI, constraints, k-points, basis
sets, pseudopotentials, and grid convergence.

## References

- [](#K%C3%BChne2020)
- [](#Iannuzzi2026)

```{youtube} 3Cw4h3MrZ8k
---
align: center
privacy_mode:
---
```
