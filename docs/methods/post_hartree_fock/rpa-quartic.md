# RPA for molecules and large supercells

This page describes the RI-RPA implementation of Del Ben *et al.* \[[](#DelBen2013),
[](#DelBen2015)\] with nuclear gradients by Stein and Hutter \[[](#Stein2024)\]. It works at the
$\Gamma$ point only (no *k*-point sampling) and is intended for molecules and for large supercells,
for example liquid water ([Del Ben *et al.* 2015](https://doi.org/10.1063/1.4927325)). The
computation time increases with the fourth power of the system size. For an introduction to RPA and
for the other RPA implementations in CP2K, see [](rpa.md).

The direct RI-RPA method computes the RPA correlation energy

$$
E^{RI-dRPA} &= - \frac{1}{4\pi}\int_{-\infty}^{\infty}d\omega\text{Tr}\left(\text{log}(1+Q(\omega))-Q(\omega)\right)\\
Q_{RS}(\omega) &= 2\sum_{ia}B_{iaR}\frac{\epsilon_a-\epsilon_i}{\left(\epsilon_a-\epsilon_i\right)^2+\omega^2}B_{iaS}
$$

as implemented in [](#DelBen2013). The integration is performed numerically. Double-hybrid
functionals are available in CP2K. Different functionals require rescaling which is achieved using
the keyword [SCALE_RPA](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.SCALE_RPA). The closely
related Laplace-transformed scaled-opposite-spin MP2 method is described in [](sos-mp2.md).

## RI-dRPA

CP2K implements two quadrature schemes: Clenshaw-Curtis and Minimax. The first requires 30-40
quadrature points, the latter 6-8 quadrature points. Minimax quadrature rules have to be
preoptimized such that not all possible numbers of quadrature points are available.

RPA correlation energies are usually combined with exact exchange energies. In CP2K, this is
available by activating the [HF](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.HF) section
which is setup like an ordinary HF section.

Several Beyond-RPA schemes are available: The Renormalized Screened Exchange (RSE) correction, the
Approximate Exchange Kernel (AXK) and the Second-Order Screened Exchange corrections (SOSEX). These
are enabled with the [RSE](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.RSE) keyword and the
[EXCHANGE_CORRECTION](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.EXCHANGE_CORRECTION)
section. This results into the following possible RPA section

```none
  &RI_RPA
    # Choose it as large as necessary and as small as possible
    QUADRATURE_POINTS 6
    # Choose it as large as possible (must be a divisor of the number of quadrature points and the number of processes)
    # -1 is default and let CP2K decide on that value
    # Larger values increase the memory demands but reduce communication
    NUM_INTEG_GROUPS -1
    # The RSE correction is only relevant for a non-HF reference
    RSE .TRUE.
    # Exchange corrections may be quite costly
    &EXCHANGE_CORRECTION [NONE|AXK|SOSEX]
      # The Hartree-Fock-based implementation scales better for larger systems but introduces more noise
      USE_HFX_IMPLEMENTATION F
      # This parameter is ignored if USE_HFX_IMPLEMENTATION is set to T
      # Larger values improve performance but increase the memory demands
      BLOCK_SIZE 16
    &END
  &END
```

## RI-dRPA Gradient Calculations

Analytical gradients are only available for RPA calculations using a minimax grid, but not for the
beyond-RPA methods (RSE, AXK, SOSEX).[](#Stein2024) The general setup is similar to RI-MP2
gradients. In addition, there are two more relevant keywords:
[DOT_PRODUCT_BLKSIZE](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.CANONICAL_GRADIENTS.DOT_PRODUCT_BLKSIZE)
and
[MAX_PARALLEL_COMM](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.CANONICAL_GRADIENTS.MAX_PARALLEL_COMM).
The first splits the contraction along the auxiliary index to improve numerical stability. By
default, this feature is turned off. The second keyword determines the number of parallel
communication channels for non-blocking communication. Larger numbers allow more overlap but
increase the memory requirements. Larger values than 3 are commonly not necessary.

## Performance Considerations

Without employing the low-scaling implementation, RI-RPA calculations scale quartically with respect
to the number of atoms. The implementation relies on parallel matrix-matrix multiplications. Their
costs can be reduced with the COSMA library and accelerated on GPUs if COSMA was configured
accordingly.

HF calculations may be accelerated using ADMM which is recommended if very diffuse basis functions
are employed. The RSE correction employs the same HF section as for the calculation of the HF
exchange energy. The costs of AXK and SOSEX depend on the implementation in use. The HF-based
implementation is cheaper for larger systems. Its performance depends on the available memory per
process group (see [GROUP_SIZE](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.GROUP_SIZE) and
[MAX_MEMORY](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.HF.MEMORY.MAX_MEMORY)). The costs
becomes neglegible for large systems. Without the HF implementation, CP2K relies on SpLA-accelerated
local matrix-matrix multiplications. These become the bottleneck for larger systems but introduce
less numerical noise than the HF-based implementation.
