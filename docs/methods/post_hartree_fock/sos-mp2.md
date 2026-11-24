# Laplace-Transformed Scaled-Opposite-Spin-MP2

The Laplace-transformed Scaled-Opposite-Spin-MP2 (LT-RI-SOS-MP2) method computes

$$
E^{LT-RI-SOS-MP2} &= - \int_{0}^{\infty}d\tau\text{Tr}\left(\overline{Q}(\tau)^2\right)\\
\overline{Q}_{RS}(\tau) &= \sum_{ia}B_{iaR}e^{-\left(\epsilon_a-\epsilon_i\right)\tau}B_{iaS}
$$

as implemented in [](#DelBen2013). The integration is performed numerically. Double-hybrid
functionals are available in CP2K.

CP2K provides two different implementations: a quartically-scaling and a low-scaling
cubically-scaling implementation (see [](low-scaling.md)). Different functionals require rescaling
which is achieved using the keyword [SCALE_S](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.SCALE_S).
Here, we will focus on the quartically-scaling implementation only. The closely related random-phase
approximation (RPA) is described in [](rpa.md).

## LT-RI-SOS-MP2

CP2K implements two quadrature schemes only the Minimax scheme requiring usually 6-8 quadrature
points.[](#DelBen2013) Minimax quadrature rules have to be preoptimized such that not all possible
numbers of quadrature points are available. Scaling of the energy contribution is possible with the
[SCALE_S](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.SCALE_S) keyword which is also used by MP2
calculations. The input is simpler than in the RPA-case

```none
  &RI_SOS_MP2
    # Works similar than in RPA
    NUM_INTEG_GROUPS -1
    # Larger values improve the accuracy but increase the costs
    QUADRATURE_POINTS 6
  &END
```

## LT-RI-SOS-MP2 Gradient Calculations

Analytical gradients are available and work similarly to RI-MP2 calculations and RI-RPA
calculations. As in case of RI-MP2, it makes use of the
[EPS_CANONICAL](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.CANONICAL_GRADIENTS.EPS_CANONICAL)
keyword. Larger values improve the numerical accuracy but increase the computational costs
significantly because CP2K assumes that the number of relevant pairs is negligible.

## Performance Considerations

Without employing the low-scaling implementation, RI-SOS-MP2 calculations scale quartically with
respect to the number of atoms. The implementation relies on parallel matrix-matrix multiplications.
Their costs can be reduced with the COSMA library and accelerated on GPUs if COSMA was configured
accordingly.
