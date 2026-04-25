# Gaussian Augmented Plane Waves

The Gaussian augmented plane wave ({term}`GAPW`) method extends GPW so that all-electron
calculations and calculations with very small-core pseudopotentials become practical in CP2K. The
central idea is to keep the smooth part of the density on the regular GPW grids while treating the
rapidly varying density close to the nuclei with atom-centered contributions.

GAPW is useful when the core electron density matters, for example in all-electron calculations,
core-level spectroscopy, magnetic properties, and some small-core pseudopotential setups. For
standard valence-only pseudopotential DFT calculations, GPW is usually simpler and faster.

## Activating GAPW

GAPW is activated in the [QS](#CP2K_INPUT.FORCE_EVAL.DFT.QS) section:

```text
&FORCE_EVAL
  METHOD Quickstep
  &DFT
    &QS
      METHOD GAPW
    &END QS
  &END DFT
&END FORCE_EVAL
```

All-electron GAPW calculations also require all-electron basis sets and `POTENTIAL ALL` for the
corresponding atomic kinds:

```text
&KIND O
  BASIS_SET SVP-MOLOPT-GGA-ae
  POTENTIAL ALL
  LEBEDEV_GRID 110
  RADIAL_GRID 80
&END KIND
```

A complete tested water example is available as [](gapw_h2o.inp). It is intentionally small and is
meant as a starting point rather than as a production-quality benchmark.

## Accuracy Parameters

Several GAPW-specific tolerances control the split between soft grid-based and hard atom-centered
contributions:

- [EPSFIT](#CP2K_INPUT.FORCE_EVAL.DFT.QS.EPSFIT) controls how Gaussian exponents are split into the
  hard and soft parts. Lowering it includes harder functions in the soft density and usually
  requires a larger [CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.CUTOFF).
- [EPSRHO0](#CP2K_INPUT.FORCE_EVAL.DFT.QS.EPSRHO0) controls the range used for the hard compensation
  density contribution.
- [EPSSVD](#CP2K_INPUT.FORCE_EVAL.DFT.QS.EPSSVD) controls the singular value decomposition tolerance
  used for projector matrices.

The atom-centered integration grid is controlled per kind with
[LEBEDEV_GRID](#CP2K_INPUT.FORCE_EVAL.SUBSYS.KIND.LEBEDEV_GRID) and
[RADIAL_GRID](#CP2K_INPUT.FORCE_EVAL.SUBSYS.KIND.RADIAL_GRID). Increasing these values can improve
the electron count and the accuracy of properties that depend on the near-core density, but it also
increases cost.

## Practical Guidance

When setting up a GAPW calculation:

- Use an all-electron basis set for `POTENTIAL ALL`, or a basis set designed for the chosen
  small-core pseudopotential.
- Inspect the electron count printed by CP2K after SCF convergence. It is a useful diagnostic for
  the quality of the hard/soft density split.
- Tighten `EPSFIT`, `EPSRHO0`, `EPSSVD`, and the atomic grids only as much as needed for the target
  property.
- Increase the density [CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.CUTOFF) when harder Gaussian
  exponents are included in the soft density.
- Prefer GPW when the calculation does not need all-electron or near-core accuracy.

## See Also

- [](#Lippert1999)
- [](#VandeVondele2006)
- [The CP2K Program Package Made Simple](https://doi.org/10.1021/acs.jpcb.5c05851)

```{youtube} L0hKLjvjIFU
---
url_parameters: ?start=4
align: center
privacy_mode:
---
```
