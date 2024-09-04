# Møller–Plesset Perturbation Theory

The MP2 method computes

$$
E^{(2)} = - \sum_{ijab} \frac{(ia|jb)[2(ia|jb)-(ib|ja)]}{\epsilon_{a}+\epsilon_{b}-\epsilon_{i}-\epsilon_{j}}
$$

with the GPW method as described in [](#DelBen2012) and [](#DelBen2013). CP2K also implements
double-hybrid functionals.

CP2K provides three different implementations: A canonical implementation, a GPW-based
implementation and a RI-based implementation. Different functionals require rescaling of different
contributions (singulet, triplet). This is achieved using the keywords
[SCALE_S](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.SCALE_S) and
[SCALE_T](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.SCALE_T) for all implementations.

## Canonical Implementation

The canonical implementation is available for GPW and GAPW reference calculations and thus allows
for core-corrections. It is the most costly implementation. Gradients (Nuclear forces, stress
tensors) are not analytically available. A simple WF_CORRELATION section looks like

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

The MP2-based implementation is available only for GPW calculations.[](#DelBen2012) It is the
cheaper than the canonical implementation. Gradients (Nuclear forces, stress tensors) are not
analytically available. A simple WF_CORRELATION section looks like

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

The RI-based implementation is the most affordable MP2-implementation.[](#DelBen2013) Analytical
gradients are available for GPW-based integrals (not for GAPW reference
calculations).[](#DelBen2015b) MME-integration does not support stress tensor calculations,
Obara-Saika integration does not implement any gradients. A simple WF_CORRELATION section looks like

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
    # This keyword determines how many copies of the rank-three tensor are kept in the memory of all CPUs. Large number of groups increase the memory demands but reduce communication (default: -1, automatic determination)
    NUMBER_INTEGRATION_GROUPS 2
  &END
&END
```

## RI-MP2 Gradient Calculations

Analytical gradients (nuclear gradients, stress tensors) are available for RI-MP2 and related
double-hybrid functionals.[](#DelBen2015b),[](#Stein2022) HF calculations may also be accelerated by
the ADMM approximation which is recommended if very diffuse basis functions are included. The
computational costs of the MP2 part is 2-3 times as high as for energy-only calculations. A typical
input file of a RI-MP2 calculation using ADMM-accelerated HF is

```none
&GLOBAL
  PRINT_LEVEL LOW
  PROJECT example
  RUN_TYPE FORCE_EVAL
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
  &DFT
    # This file contains basis sets for common elements (H, C, N, O; for F, Al, Cl, Si, P, S only TZ basis sets)
    BASIS_SET_FILE_NAME BASIS_RI_cc-TZ
    BASIS_SET_FILE_NAME BASIS_ADMM
    POTENTIAL_FILE_NAME POTENTIAL
    # Plug in your preconverged wfn file here
    WFN_RESTART_FILE_NAME ./example_HF.wfn
    # Use if no RI basis set is available (possible values: SMALL, MEDIUM, LARGE, HUGE)
    AUTO_BASIS RI_AUX LARGE

    &AUXILIARY_DENSITY_MATRIX_METHOD
      # Purification is not implemented
      ADMM_PURIFICATION_METHOD NONE
      # Try different options (check with reference values or use the default)
      EXCH_CORRECTION_FUNC PBEX
      # other methods are not implemented
      METHOD BASIS_PROJECTION
    &END AUXILIARY_DENSITY_MATRIX_METHOD
    &MGRID
      # Adjust as usual
      CUTOFF 600
      REL_CUTOFF 50
    &END MGRID
    &POISSON
      PERIODIC XYZ
      POISSON_SOLVER WAVELET
    &END POISSON
    &QS
      # Choose as tight as necessary
      EPS_DEFAULT 1.0E-10
      METHOD GPW
    &END QS
    &SCF
      # Choose as tight as necessary
      EPS_SCF 1.0E-8
      MAX_SCF 100
      SCF_GUESS RESTART
    &END SCF
    &XC
      &HF
        FRACTION 1.0000000
        &INTERACTION_POTENTIAL
          # Adjust the cutoff radius according to your cell
          CUTOFF_RADIUS 1.5
          POTENTIAL_TYPE TRUNCATED
          T_C_G_DATA t_c_g.dat
        &END INTERACTION_POTENTIAL
        &SCREENING
          # Tune these parameters
          EPS_SCHWARZ 1.0E-10
          EPS_SCHWARZ_FORCES 1.0E-5
          SCREEN_ON_INITIAL_P .FALSE.
        &END SCREENING
      &END HF
      &WF_CORRELATION
        MEMORY 500
        NUMBER_PROC 1
        &CANONICAL_GRADIENTS
          # The default is usually good enough
          EPS_CANONICAL 1E-6
          # This is the option which should always work
          # Try to set it to .TRUE. (may segfault)
          FREE_HFX_BUFFER .FALSE.
          &CPHF
            # Choose as tight as necessary
            EPS_CONV 1.0E-6
            # Smaller values of EPS_CONV may require more iterations
            MAX_ITER 10
          &END CPHF
        &END CANONICAL_GRADIENTS
        &INTEGRALS
          &WFC_GPW
            # Adjust these parameters to your needs of accuracy
            CUTOFF 200
            EPS_FILTER 1.0E-12
            EPS_GRID 1.0E-8
            REL_CUTOFF 50
          &END WFC_GPW
        &END INTEGRALS
        &RI_MP2
          # Determine an automatic block size, if memory is low, set it to 1
          BLOCK_SIZE -1
        &END RI_MP2
      &END WF_CORRELATION
      &XC_FUNCTIONAL NONE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
  &PRINT
    &FORCES
    &END FORCES
  &END PRINT
  &SUBSYS
    # Adjust cell information and coordinations as necessary
    &CELL
      ABC [angstrom] 5.0 5.0 5.0
      PERIODIC XYZ
    &END CELL
    &COORD
      O       0.000000    0.000000    -0.211000
      H       0.000000   -0.844000     0.495000
      H       0.000000    0.744000     0.495000
    &END COORD
    &KIND H
      # orbital and RI basis set should match
      BASIS_SET cc-TZ
      BASIS_SET RI_AUX RI_TZ
      # Use the largest affordable one
      BASIS_SET AUX_FIT cpFIT3
      POTENTIAL GTH-HF-q1
    &END KIND
    &KIND O
      BASIS_SET cc-TZ
      BASIS_SET RI_AUX RI_TZ
      BASIS_SET AUX_FIT cpFIT3
      POTENTIAL GTH-HF-q6
    &END KIND
    &TOPOLOGY
      &CENTER_COORDINATES
      &END CENTER_COORDINATES
    &END TOPOLOGY
  &END SUBSYS
&END FORCE_EVAL
```

A few more words are required regarding systems with many degenerate occupied orbital pairs as they
are found in structures with many symmetry-equivalent atoms. In that case, non-diagonal elements of
the density matrix of almost degenerate occupied-orbital pairs have to be calculated explicitly for
numerical reason. This is significantly more expensive than for non-degenerate orbital pairs. This
behavior can be tuned using the
[EPS_CANONICAL](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.CANONICAL_GRADIENTS.EPS_CANONICAL)
keyword.

## Performance Considerations

MP2-calculations are generally very expensive for large systems due to their quintical scaling with
respect to the number of atoms. Thus, for larger systems, one should switch to the cheaper RPA or
SOS-MP2 methods. CP2K employs local matrix-matrix multiplications for contractions. These can be
accelerated with the SpLA library if CP2K was linked against it. In that way, GPU acceleration is
possible if no GPU-accelerated DGEMM implementation is available and SpLA was configured
accordingly.
