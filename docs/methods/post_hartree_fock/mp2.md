# Møller–Plesset Perturbation Theory

The MP2 method computes

$$
E^{(2)} = - \sum_{ijab} \frac{(ia|jb)[2(ia|jb)-(ib|ja)]}{\epsilon_{a}+\epsilon_{b}-\epsilon_{i}-\epsilon_{j}}
$$

with the GPW method as described in [](#DelBen2012) and [](#DelBen2013). CP2K also implements double-hybrid functionals.

CP2K provides three different implementations: A canonical implementation, a GPW-based implementation and a RI-based implementation. Different functionals require rescaling of different contributions (singulet, triplet). This is achieved using the keywords [SCALE_S](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.SCALE_S) and [SCALE_T](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.SCALE_T) for all implementations.

## Canonical Implementation

The canonical implementation is available for GPW and GAPW reference calculations and thus allows for core-corrections. It is the most costly implementation. Gradients (Nuclear forces, stress tensors) are not analytically available. A simple WF_CORRELATION section looks like

```none
&WF_CORRELATION
  # See prelimaries
  MEMORY    1200
  NUMBER_PROC  1
  # Only if not the Coulomb potential is required
  &INTEGRALS
  &END
  &MP2
    METHOD DIRECT_CANONICAL
  &END
&END
```

## GPW-based Implementation

The MP2-based implementation is available only for GPW calculations.[](#DelBen2012) It is the cheaper than the canonical implementation. Gradients (Nuclear forces, stress tensors) are not analytically available. A simple WF_CORRELATION section looks like

```none
&WF_CORRELATION
  # See prelimaries
  MEMORY    1200
  NUMBER_PROC  1
  # Only if not the Coulomb potential is required
  &INTEGRALS
    &WFC_GPW
      # To be optimized
      CUTOFF 200
      REL_CUTOFF 50
    &END
  &END
  &MP2
    METHOD MP2_GPW
  &END
&END
```

## RI-based Implementation

The RI-based implementation is the most affordable MP2-implementation.[](#DelBen2013) Analytical gradients are available for GPW-based integrals (not for GAPW reference calculations).[](#DelBen2015b) MME-integration does not support stress tensor calculations, Obara-Saika integration does not implement any gradients. A simple WF_CORRELATION section looks like

```none
&WF_CORRELATION
  # See prelimaries
  MEMORY    1200
  NUMBER_PROC  1
  # Only if not the Coulomb potential is required
  &INTEGRALS
    &WFC_GPW
      # To be optimized
      CUTOFF 200
      REL_CUTOFF 50
    &END
  &END
  &RI_MP2
    # Larger block sizes require more memory but reduce communication, -1 (default) let CP2K choose it
    BLOCK_SIZE 2
    # This keyword determines how many copies of the rank-three tensor are kept in the memory of all CPUs. Large number of groups increasethe memory demands but reduce communication (default: -1, automatic determination)
    NUMBER_INTEGRATION_GROUPS 2
  &END
&END
```

## Performance Considerations

MP2-calculations are generally very expensive for large systems due to their quintical scaling with respect to the number of atoms. Thus, for larger systems, one should switch to the cheaper RPA or SOS-MP2 methods.
CP2K employs local matrix-matrix multiplications for contractions. These can be accelerated with the SpLA library if CP2K was linked against it. In that way, GPU acceleration is possible if no GPU-accelerated DGEMM implementation is available and SpLA was configured accordingly.
