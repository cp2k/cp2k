# Preliminaries

Here, we explain the relevant preparation steps and things to consider with all post-Hartree-Fock
methods regarding choice of basis sets and pseudopotentials. Because of the different building
blocks, we have to distinguish at least two kinds of basis sets: the primary basis set (PBS) to
represent the orbital functions and the RI basis set to expand products of PBS functions into a set
of auxiliary functions such that the four-center electron repulsion integrals (in Mulliken-notation)
are approximated by

$$
(ia|jb)=\sum_{PQ}(ia|P) (P|Q)^{-1} (Q|jb)
$$

On top of the PBS and the RI basis set, the ADMM approximation of HF may introduce a third kind of
basis set (see there for further information).

## Step 1: Choice of Primary Basis Sets

Because of the slow convergence of the calculated properties from post-Hartree-Fock calculations
with respect to the size of the PBS, post-Hartree-Fock methods demand for larger PBSs than ordinary
hybrid DFT or HF calculations and for very accurate results a basis set convergence scheme employing
consistent PBS for accurate results. For this purpose, post-Hartree-Fock calculations are commonly
performed with correlation-consistent (cc)-basis sets. Without an extrapolation scheme, we recommend
basis sets of at least triple-zeta (cc-TZ) or augmented double-zeta (aug-cc-DZ) quality. With an
extrapolation scheme, one needs a second basis set for each atom kind of different size. The usual
extrapolation formula for correlation energies is given by

$$
E(X) = A + \frac{B}{X^3}
$$

with $X$ being 2, 3 or 4 in case of DZ, TZ or QZ zeta basis sets, respectively. If possible, one
should refrain from non-augmented DZ basis sets as their results are usually not converged enough
for extrapolation. Suitable basis sets can be found in the `BASIS_RI_cc-TZ` and `BASIS_ccGRB` basis
set files in the data directory.

## Step 2: Choice of RI Basis Sets (Only for RI-based Methods)

The RI approximation asks for an additional auxiliary basis set. They are usually specified as
[BASIS_SET RI_AUX](#CP2K_INPUT.FORCE_EVAL.SUBSYS.KIND.BASIS_SET) in the input file. CP2K provides
several options to find suitable RI basis sets:

- For a few main group elements, the `BASIS_RI_cc-TZ` file contains RI basis sets for the given
  PBSs. The literature may help with missing elements.
- Optimize RI basis sets yourself. This is recommended if a larger number of calculations is
  required to reduce the computational costs significantly. Check the
  [OPT_RI_BASIS](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI.OPT_RI_BASIS) and the regtest suite
  (QS/regtest-ri-opt) fur further information.
- In all other cases, CP2K can generate RI basis sets automatically (see the
  [AUTO_BASIS RI_AUX](#CP2K_INPUT.FORCE_EVAL.DFT.AUTO_BASIS). These basis sets provide a decent
  accuracy but are roughly twice as large as optimized basis sets.

## Step 3: Pseudopotentials

Pseudopotentials should reflect the underlying hybrid or HF calculation used to determine orbitals
and orbital energies (see the respective sections).

## Step 4: Preoptimize Reference Orbitals

It is recommended to restart the orbital calculation step from preoptimized HF or hybrid orbitals.
This reduces computational costs because DFT and HF calculations do not perform well in case of the
large number of CPU cores required for post-Hartree-Fock methods.

## Step 5: Setup of Integral Calculation Step

The integration step is setup with the
[INTEGRALS](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.INTEGRALS) section. Because most densities
occurring in the formalisms of the post-Hartree-Fock methods have a zero net-charge, such that the
Coulomb operator does not need to be truncated it is required by HF. CP2K has its own integration
routines for post-Hartree-Fock methods exploiting this circumstance. The available integration
methods (compare the [ERI_METHOD](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.INTEGRALS.ERI_METHOD)
keyword) are Gaussian-Plane-Wave integration (GPW), Minimax-Ewald integration (MME) and Obara Saika
integration (OS).

GPW is fully supported, MME does not implement stress tensors and does not support short-ranged
operators (truncated Coulomb, erfc-Coulomb, ...), OS implements gradients only in case of low
scaling methods and does not support automatically generated RI basis sets. The GPW and MME
integration methods are further configured in the respective
[WFC_GPW](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.INTEGRALS.WFC_GPW) and
[ERI_MME](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.INTEGRALS.ERI_MME) sections. In case of the
GPW integration method, the respective CUTOFF parameters should be tuned as known from ordinary
G(A)PW integrations in CP2K although the primary cutoff parameter can be chosen much smaller
(150-300 Ry) then usual.

The RI metric can be set by the user with the [RI](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI)
section. In case of low scaling methods (with activated
[LOW_SCALING](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.LOW_SCALING) section), the metric is by
default the overlap metric, else it will be the potential operator will be used. If the default is
not requested, the RI section should be set.

## General Structure of Post-Hartree-Fock Methods

The general structure is very similar. Each post-Hartree-Fock method has its own section which has
to be activated to start the respective calculation. Each section has subsections to further
configure the respective kind of calculation. Most post-Hartree-Fock calculations require the
following structure

```none
&WF_CORRELATION
  # Determine the available memory
  MEMORY 1000
  # Determines the size of sub groups used in most methods
  # The lower the number, the less communication but the more memory (only relevant for very large systems)
  GROUP_SIZE 1
  &INTEGRALS
    # Here, we show the setup with the GPW integration
    ERI_METHOD GPW
    # Tune the accuracy and performance of the integral calculation step (default is usually sufficient)
    SIZE_LATTICE_SUM 5
    &WFC_GPW
      # Set it as large as necessary and as small as possible (check energies)
      CUTOFF 200
      REL_CUTOFF 50
    &END
    # Only necessary for special functionals
    &INTERACTION_POTENTIAL
    &END
    # Relevant only if low-scaling methods or a potential operator without short-range contributions is requested
    &RI
      # Set up the RI metric
      &RI_METRIC
      &END
    &END
  &END
  # Set the section corresponding to the requested method (mixtures of different methods are not possible)
  &<NAME_OF_METHOD>
    <OPTIONS_TO_CONFIGURE_METHOD>
  &END
&END
```

## Gradient Calculations

Gradient calculations are available for RI-MP2, RI-RPA and RI-SOS-MP2 calculations. In case of
low-scaling implementations, consult the
[CPHF](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.LOW_SCALING.CPHF) section in the LOW_SCALING
section, otherwise the
[CANONICAL_GRADIENTS](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.CANONICAL_GRADIENTS) section for
further information. The setup of the CANONICAL_GRADIENTS section islooks like

```none
&CANONICAL_GRADIENTS
  # This threshold switches to explicit contractions for almost degenerate orbital pairs (excluding diagonal elements)
  EPS_CANONICAL 1.0E-6
  # Try to turn it on, but it may crash in the execution step
  FREE_HFX_BUFFER .FALSE.
  &CPHF
    # This threshold is to be optimized (lower values increase accuracy but require more time)
    EPS_CONV 1.0E-6
    MAX_ITER 20
  &END
&END
```

## Basis Set Superposition Error (BSSE)

Energy calculations using post-Hartree-Fock methods suffer from severe basis set superposition error
(BSSE). In CP2K, this is enhanced by the many different basis sets in use. Check the
[BSSE](#CP2K_INPUT.FORCE_EVAL.BSSE) section for further information on the automatic setup of BSSE
calculations and the [GHOST](#CP2K_INPUT.FORCE_EVAL.SUBSYS.KIND.GHOST) for the manual setup of BSSE
calculations.
