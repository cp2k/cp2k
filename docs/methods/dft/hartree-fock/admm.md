# HFX with ADMM

The auxiliary density matrix method ({term}`ADMM`) reduces the cost of Hartree-Fock exchange in
hybrid DFT calculations by projecting the density matrix from the primary orbital basis onto a
smaller auxiliary basis. CP2K evaluates exact exchange in the auxiliary basis and adds a correction
term for the difference between the primary and auxiliary exchange descriptions.

ADMM is most useful when exact exchange is the bottleneck, especially with larger or more diffuse
Gaussian basis sets. It is commonly used for hybrid DFT, and the ADMM2 variant is also supported by
several post-SCF methods that reuse exact-exchange machinery.

## Basic Setup

An ADMM calculation needs three pieces of input:

- a hybrid functional or another setup that evaluates Hartree-Fock exchange,
- an auxiliary basis set for each atomic kind, specified with `BASIS_SET AUX_FIT`,
- an [AUXILIARY_DENSITY_MATRIX_METHOD](#CP2K_INPUT.FORCE_EVAL.DFT.AUXILIARY_DENSITY_MATRIX_METHOD)
  section that selects the ADMM variant and correction functional.

For example:

```text
&DFT
  BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH
  BASIS_SET_FILE_NAME BASIS_ADMM_UZH
  POTENTIAL_FILE_NAME POTENTIAL_UZH
  &AUXILIARY_DENSITY_MATRIX_METHOD
    ADMM_TYPE ADMMS
    EXCH_CORRECTION_FUNC PBEX
  &END AUXILIARY_DENSITY_MATRIX_METHOD
  &XC
    &XC_FUNCTIONAL PBE
    &END XC_FUNCTIONAL
    &HF
      FRACTION 0.25
    &END HF
  &END XC
&END DFT

&SUBSYS
  &KIND O
    BASIS_SET ccGRB-D-q6
    BASIS_SET AUX_FIT admm-dz-q6
    POTENTIAL GTH-HYB-q6
  &END KIND
&END SUBSYS
```

## Choosing the Auxiliary Basis

The auxiliary basis should be chosen for the primary basis family and for the intended accuracy.
For MOLOPT-style calculations, the `BASIS_ADMM_MOLOPT` family provides compact auxiliary bases. The
newer UZH basis-set collection includes `BASIS_ADMM_UZH` and related basis files for correlation
consistent setups. All-electron calculations can use all-electron auxiliary basis sets when
available.

The auxiliary basis is part of the approximation. A too small auxiliary basis can make the exchange
correction large and reduce accuracy; a too large one gives back less speedup. For production work,
test at least one larger auxiliary basis or compare against a smaller reference system without ADMM.

## Choosing the ADMM Variant

[ADMM_TYPE](#CP2K_INPUT.FORCE_EVAL.DFT.AUXILIARY_DENSITY_MATRIX_METHOD.ADMM_TYPE) is a shortcut that
sets the projection, purification, and scaling options consistently. `ADMM1` and `ADMM2` are the
original variants, while `ADMMS`, `ADMMP`, and `ADMMQ` use additional models introduced later.
`ADMM2` is often the most broadly supported variant for workflows beyond ground-state hybrid DFT.

The
[EXCH_CORRECTION_FUNC](#CP2K_INPUT.FORCE_EVAL.DFT.AUXILIARY_DENSITY_MATRIX_METHOD.EXCH_CORRECTION_FUNC)
keyword selects the exchange functional used for the ADMM correction. It should be chosen
consistently with the exchange part of the main exchange-correlation setup; `PBEX` is a common
choice for PBE-based hybrid calculations.

## Practical Checks

When using ADMM:

- keep the same primary basis and potential convergence checks that would be used without ADMM,
- check the sensitivity to the auxiliary basis size,
- compare total energies, forces, or target properties against a non-ADMM reference for a small
  representative system,
- remember that ADMM accelerates the exchange calculation but does not replace convergence of the
  primary Gaussian basis, real-space grid, or SCF thresholds.

## See Also

- [](#Guidon2009)
- [](#Guidon2010)
- [](#Merlot2014)
- [The CP2K Program Package Made Simple](https://doi.org/10.1021/acs.jpcb.5c05851)

```{youtube} snG4fbpI0_g
---
url_parameters: ?start=1673
align: center
privacy_mode:
---
```
