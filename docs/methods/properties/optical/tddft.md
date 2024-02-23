# Time-Dependent DFT

This is a short tutorial on how to run linear-response time-dependent density functional theory
(LR-TDDFT) computations for absorption and emission spectroscopy. The TDDFT module enables a
description of excitation energies and excited-state computations within the Tamm-Dancoff
approximation (TDA) featuring GGA and hybrid functionals as well as semi-empirical simplified TDA
kernels. The details of the implementation can be found in [Strand2019] and in [Hehn2022] for
corresponding excited-state gradients. Note that the current module is based on an earlier TDDFT
implementation [](#Iannuzzi2005). Please cite these papers if you were to use the TDDFT module for
the computation of excitation energies ([Strand2019], [](#Iannuzzi2005)) or excited-state gradients
([Hehn2022]).

## Brief theory recap

The implementation in CP2K is based on the Tamm-Dancoff approximation (TDA), which describes each
excited state $p$ with the excitation energy $\Omega^p$ and the corresponding excited-state
eigenvectors $\mathbf{X}^p$ as an Hermitian eigenvalue problem

$$
\mathbf{A} \mathbf{X}^p &= \Omega^p \mathbf{S} \mathbf{X}^p \, , \\
\sum_{\kappa k} [ F_{\mu \kappa \sigma} \delta_{ik} - F_{ik \sigma} S_{\mu \kappa} ] X^p_{\kappa k \sigma} + \sum_{\lambda} K_{\mu \lambda \sigma} [\mathbf{D}^{{\rm{\tiny{X}}}p}] C_{\lambda i \sigma} &=  \sum_{\kappa} \Omega^p S_{\mu \kappa} X^p_{\kappa i \sigma} \, . 
$$

The Hermitian matrix $\mathbf{A}$ contains as zeroth-order contributions the difference in the
Kohn-Sham (KS) orbital energies $\mathbf{F}$, and to first order kernel contributions $\mathbf{K}$
which comprise --depending on the chosen density functional approximation -- Coulomb $\mathbf{J}$
and exact exchange $\mathbf{K}^{\rm{\tiny{EX}}}$ contributions as well as contributions due to the
exchange-correlation (XC) potential $\mathbf{V}^{\rm{\tiny{XC}}}$ and kernel
$\mathbf{f}^{\rm{\tiny{XC}}}$,

$$
F_{\mu \nu \sigma} [\mathbf{D}] &= h_{\mu \nu} + J_{\mu \nu \sigma} [\mathbf{D}] - a_{\rm{\tiny{EX}}}K^{\rm{\tiny{EX}}}_{\mu \nu \sigma} [\mathbf{D}] + V_{\mu \nu \sigma}^{\rm{\tiny{XC}}} \, , \\
K_{\mu \nu \sigma} [\mathbf{D}^{{\rm{\tiny{X}}}p}] &=  J_{\mu \nu \sigma} [\mathbf{D}^{{\rm{\tiny{X}}}p}] - a_{\rm{\tiny{EX}}} K^{\rm{\tiny{EX}}}_{\mu \nu \sigma}[\mathbf{D}^{{\rm{\tiny{X}}}p}] + \sum_{\kappa \lambda \sigma'} f^{\rm{\tiny{XC}}}_{\mu \nu \sigma,\kappa \lambda \sigma'} D_{\kappa \lambda \sigma'}^{{\rm{\tiny{X}}}p} \, .
$$

$\mathbf{S}$ denotes the atomic-orbital overlap matrix, $\mathbf{C}$ the occupied ground-state KS
orbitals and $\mathbf{D}$ and $\mathbf{D}^{\rm{\tiny{X}}}$ ground-state and response density
matrices,

$$
D_{\mu \nu \sigma} &= \sum_k C_{\mu k \sigma} C_{\nu k \sigma}^{\rm{T}} \, , \\
D_{\mu \nu \sigma}^{{\rm{\tiny{X}}}p} &= \frac{1}{2} \sum_{k} ( X^p_{\mu k \sigma} C_{\nu k \sigma}^{\rm{T}} + C_{\mu k \sigma} (X^p_{\nu k \sigma})^{\rm{T}} ) \, .
$$

Within the current implementation, symmetrization and orthogonalization of the response density
matrix is ensured at each step of the Davidson algorithm. The current implementation features to
approximate the exact exchange contribution of hybrid functionals using the auxiliary density matrix
method (ADMM). Furthermore, the standard kernel can be approximated using the semi-empirical
simplified Tamm-Dancoff approximation (sTDA), neglecting in this case XC contributions and
approximating both Coulomb and exchange contributions $\mathbf{J}$ and $\mathbf{K}$ using
semi-empirical operators $\boldsymbol{\gamma}^{\rm{\tiny{J}}}$ and
$\boldsymbol{\gamma}^{\rm{\tiny{K}}}$ depending on the interatomic distance $R_{AB}$ of atoms $A$
and $B$,

$$
\gamma^{\rm{\tiny{J}}}(A,B) &= \left ( \frac{1}{(R_{AB})^{\alpha} + \eta^{-\alpha}}  \right)^{1/\alpha} \, , \\
\gamma^{\rm{\tiny{K}}}(A,B) &= \left ( \frac{1}{(R_{AB})^{\beta} +( a_{\rm{\tiny{EX}}}\eta)^{- \beta} } \right )^{1/\beta} \, ,
$$

that depend on the chemical hardness $\eta$, the Fock-exchange mixing parameter $a_{\rm{\tiny{EX}}}$
and powers of $\alpha$ and $\beta$ for either Coulomb and exchange interactions.

Within the current implementation, oscillator strengths can be calculated for molecular systems in
the length form and for periodic systems using the velocity form (see [Strand2019]).

Based on the TDA eigenvalue problem, excited-state gradients can be formulated based on a
variational Lagrangian for each excited state $p$,

$$
L [\mathbf{X}, \mathbf{C}, \Omega, \bar{\mathbf{W}}^{\rm{\tiny{X}}}, \bar{\mathbf{Z}}, \bar{\mathbf{W}}^{\rm{\tiny{C}}} ] &= \Omega -  \sum_{\kappa \lambda k l \sigma} \Omega ( X_{\kappa k \sigma }^{\rm{T}} S_{\kappa \lambda } X_{\lambda l \sigma} - \delta_{kl} ) \\
&-  \sum_{kl \sigma} ( \bar{W}_{kl \sigma}^{\rm{\tiny{X}}} )^{\rm{T}}  \sum_{\kappa \lambda} \frac{1}{2} ( C_{\kappa k \sigma}^{\rm{T}} S_{\kappa \lambda} X_{\lambda l \sigma}  + X_{\kappa k \sigma}^{\rm{T}}S_{\kappa \lambda} C_{\lambda l \sigma}) \\
&+ \sum_{\kappa k \sigma}( \bar{Z}_{\kappa k \sigma})^{\rm{T}}  \sum_{\lambda} ( F_{\kappa \lambda \sigma}C_{\lambda k \sigma} - S_{\kappa \lambda } C_{\lambda k \sigma} \varepsilon_{k \sigma}) \\
&- \sum_{kl\sigma} (\bar{W}^{\rm{\tiny{C}}}_{kl \sigma})^{\rm{T}}  ( S_{kl \sigma} - \delta_{kl})\, .
$$

introducing Lagrange multipliers $\bar{\mathbf{W}}^{\rm{\tiny{X}}}$,
$\bar{\mathbf{W}}^{\rm{\tiny{C}}}$, and $\bar{\mathbf{Z}}$ to ensure stationarity of the
corresponding ground-state (GS) equations and to account for the geometric dependence of the
Gaussian orbitals and thus requiring to solve the Z vector equation iteratively.

## The LR-TDDFT input section

To compute absorption spectra, parameters defining the LR-TDDFT computation have to be specified in
the [TDDFPT](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.TDDFPT) subsection. Furthermore,
[RUN_TYPE](#CP2K_INPUT.GLOBAL.RUN_TYPE) has to be set to `ENERGY` and the underlying KS ground-state
reference has to be specified in the [DFT](#CP2K_INPUT.FORCE_EVAL.DFT) section.

The most important keywords and subsections of [TDDFPT](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.TDDFPT)
are:

- [KERNEL](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.TDDFPT.KERNEL): option for the kernel matrix
  $\mathbf{K}$ to choose between the full kernel for GGA or hybrid functionals and the simplified
  TDA kernel
- [NSTATES](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.TDDFPT.NSTATES): number of excitation energies to be
  computed
- [CONVERGENCE](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.TDDFPT.CONVERGENCE): threshold for the convergence
  of the Davidson algorithm
- [RKS_TRIPLETS](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.TDDFPT.RKS_TRIPLETS): option to switch from the
  default computation of singlet excitation energies to triplet excitation energies
- [RESTART](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.TDDFPT.RESTART): the keyword enables the restart of
  the TDDFPT computation if a corresponding restart file (.tdwfn) exists
- [WFN_RESTART_FILE_NAME](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.TDDFPT.WFN_RESTART_FILE_NAME): for a
  restart of the TDDFPT computation, the name of the restart file has to be specified using this
  keyword

To compute excited-state gradients and thus corresponding fluorescence spectra, the excited state to
be optimized furthermore has to be specified by adding the subsection `EXCITED_STATES` of the
section `DFT`.

## Simple examples

### Excitation energies for acetone

The following input is a standard input for calculating excitation energies with the hybrid
functional PBE0. It should be noted that it is possible to compute the exact exchange integrals for
hybrid functionals analytically at however high computational costs. It is therefore recommended to
use the Auxiliary Density Matrix Method (ADMM) [](#Guidon2010) to approximate the exact exchange
contribution. For GGAs, the corresponding ADMM sections are not required and it is sufficient to
specify the chosen functional in the subsection `&XC_FUNCTIONAL` of the `&XC` section.

```none
&GLOBAL
  PROJECT S20Acetone
  RUN_TYPE ENERGY
  PREFERRED_DIAG_LIBRARY SL
  PRINT_LEVEL medium
&END GLOBAL
&FORCE_EVAL
  METHOD Quickstep
    &PROPERTIES
      &TDDFPT                          ! input section for TDDFPT
       KERNEL FULL                       ! specification of the underlying kernel matrix K
                                         ! FULL kernel is for GGA and hybrid functional computations
                                         ! sTDA kernel is referring to a semi-empirical sTDA computation
       NSTATES 10                      ! specifies the number of excited states to be computed
       MAX_ITER   100                  ! number of iterations for the Davidson algorithm
       CONVERGENCE [eV] 1.0e-7         ! convergence threshold in eV
       RKS_TRIPLETS F                  ! Keyword to choose between singlet and triplet excitations
     !  &XC                            ! If choosing kernel FULL, the underlying functional can be
     !   &XC_FUNCTIONAL PBE0             ! specified by adding an XC section
     !   &END XC_FUNCTIONAL              ! The functional can be chosen independently from the chosen
     !  &END XC                        ! GS functional except when choosing ADMM 
     !  &MGRID                         ! It is also possible to choose a separate grid for the real-space 
     !    CUTOFF 800                   ! integration of the response density in the TDDFT part,
     !    REL_CUTOFF 80                ! however, in general a consistent setup for GS and ES is recommended
     !  &END MGRID 
      &END TDDFPT
    &END PROPERTIES
  &DFT
    &QS
      METHOD GPW
      EPS_DEFAULT 1.0E-17
      EPS_PGF_ORB 1.0E-20
    &END QS
    &SCF
      SCF_GUESS restart
      &OT
         PRECONDITIONER FULL_ALL
         MINIMIZER DIIS
      &END OT
      &OUTER_SCF
         MAX_SCF 900
         EPS_SCF 1.0E-7
      &END OUTER_SCF
      MAX_SCF 10
      EPS_SCF 1.0E-7
    &END SCF
    POTENTIAL_FILE_NAME POTENTIAL_UZH
    BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH
    BASIS_SET_FILE_NAME BASIS_ADMM_UZH
    &MGRID
      CUTOFF 800
      REL_CUTOFF 80
    &END MGRID
    &AUXILIARY_DENSITY_MATRIX_METHOD       ! For hybrid functionals, it is recommended to choose ADMM 
      METHOD BASIS_PROJECTION              ! the ADMM environment for ground and excited state has to be 
      EXCH_SCALING_MODEL NONE              ! identical
      EXCH_CORRECTION_FUNC NONE            ! Triple-zeta auxiliary basis sets are recommended (see below)
      ADMM_PURIFICATION_METHOD NONE        ! For periodic systems (see below), only specific ADMM options
    &END AUXILIARY_DENSITY_MATRIX_METHOD   ! are available
    &POISSON
       PERIODIC NONE
       POISSON_SOLVER WAVELET
    &END
    &XC
     &XC_FUNCTIONAL PBE0
     &END XC_FUNCTIONAL
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC [angstrom] 14.0 14.0 14.0
      PERIODIC NONE
    &END CELL
    &COORD
       C 0.000000 1.282877 -0.611721
       C 0.000000 -1.282877 -0.611721
       C 0.000000 0.000000 0.185210
       O 0.000000 0.000000 1.392088
       H 0.000000 2.133711 0.059851
       H -0.876575 1.319344 -1.256757
       H 0.876575 1.319344 -1.256757
       H 0.000000 -2.133711 0.059851
       H 0.876575 -1.319344 -1.256757
       H -0.876575 -1.319344 -1.256757
    &END COORD
    &TOPOLOGY
     &CENTER_COORDINATES T
     &END
    &END
    &KIND H
      BASIS_SET ORB DZVP-MOLOPT-PBE0-GTH-q1  ! in general it is recommended to use larger basis  sets
      BASIS_SET AUX_FIT admm-dzp-q1          ! for the primary and auxiliary basis (TZVP/tzp)
      POTENTIAL GTH-PBE0-q1
    &END KIND
    &KIND O
      BASIS_SET ORB DZVP-MOLOPT-PBE0-GTH-q6
      BASIS_SET AUX_FIT admm-dzp-q6
      POTENTIAL GTH-PBE0-q6
    &END KIND
    &KIND C
      BASIS_SET ORB DZVP-MOLOPT-PBE0-GTH-q4
      BASIS_SET AUX_FIT admm-dzp-q4
      POTENTIAL GTH-PBE0-q4
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
```

In the resulting output file, there is a `TDDFPT` section reporting the steps of the calculations.
The initial guess is referring to the zeroth-order KS energy differences with the printout also
listing the transition from the corresponding occupied to virtual orbital.

```none
 -------------------------------------------------------------------------------
 -                            TDDFPT Initial Guess                             -
 -------------------------------------------------------------------------------
          State         Occupied      ->      Virtual          Excitation
          number         orbital              orbital          energy (eV)
 -------------------------------------------------------------------------------
             1               12                   13              6.63336
             2               11                   13              9.62185
             3               10                   13             10.34045
             4                9                   13             10.78846
             5                8                   13             11.05691
             6               12                   14             12.07194
             7                7                   13             12.42370
             8                6                   13             12.63557
             9               12                   15             12.65312
            10                5                   13             12.68633
```

After convergence of the iterative Davidson algorithm, CP2K is printing for each of the calculated
excited states the excitation energy in eV, the corresponding transition dipole as well as the
oscillator strength. The form of the dipole transition integrals can be chosen by modifying the
keyword [DIPOLE_FORM](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.TDDFPT.DIPOLE_MOMENTS.DIPOLE_FORM). Possible
options are `BERRY` (valid for fully periodic systems only), `LENGTH` (valid for molecular systems
only) and `VELOCITY` (both). When referring to the length form, the reference point to calculate
electric dipole moments can be chosen in the subsection by specifying the coordinates of the
reference point
[REFERENCE_POINT](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.TDDFPT.DIPOLE_MOMENTS.REFERENCE_POINT) or by
choosing one of the options `COM` (center of mass), `COAC` (center of atomic charges),
`USER_DEFINED` (user-defined point) and `ZERO` (origin of the coordinate system) for the
[REFERENCE](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.TDDFPT.DIPOLE_MOMENTS.REFERENCE) keyword.

```none
 -------------------------------------------------------------------------------
 -  TDDFPT run converged in 10 iteration(s) 
 -------------------------------------------------------------------------------

 R-TDDFPT states of multiplicity 1
 Transition dipoles calculated using velocity formulation

         State    Excitation        Transition dipole (a.u.)        Oscillator
         number   energy (eV)       x           y           z     strength (a.u.)
         ------------------------------------------------------------------------
 TDDFPT|      1       4.67815   2.4840E-08 -1.9187E-08  5.9673E-08   5.21038E-16
 TDDFPT|      2       8.73074  -8.1610E-09  1.0502E-08  9.1823E-08   1.84132E-15
 TDDFPT|      3       9.17033   1.6661E-02  2.4911E-08 -9.2612E-08   6.23670E-05
 TDDFPT|      4       9.70094  -2.3716E-09  7.2323E-07  7.0660E-01   1.18663E-01
 TDDFPT|      5       9.88658  -3.1922E-08 -4.1476E-01  1.0768E-06   4.16669E-02
 TDDFPT|      6      10.78269   3.5483E-08  4.2695E-01 -1.3123E-07   4.81550E-02
 TDDFPT|      7      10.82893   2.3885E-01 -3.0329E-08  4.5882E-09   1.51356E-02
 TDDFPT|      8      10.92800  -2.4034E-07 -2.4299E-08  5.8160E-08   1.65293E-14
 TDDFPT|      9      11.32208   8.3808E-08 -6.8844E-01  1.9166E-07   1.31469E-01
 TDDFPT|     10      11.97282   1.2605E-08  1.7029E-08  3.3026E-01   3.19930E-02
```

A more detailed analysis of the excitations is given subsequently, listing orbital contributions for
each transition with the corresponding excitation amplitudes. The keyword
[MIN_AMPLITUDE](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.TDDFPT.MIN_AMPLITUDE) regulates the threshold for
the smallest amplitude to print.

```none
-------------------------------------------------------------------------------
 -                            Excitation analysis                              -
 -------------------------------------------------------------------------------
        State             Occupied              Virtual             Excitation
        number             orbital              orbital             amplitude
 -------------------------------------------------------------------------------
             1   4.67815 eV
                                12                   13               0.998438
             2   8.73074 eV
                                10                   13               0.995560
                                 6                   13               0.087643
             3   9.17033 eV
                                 9                   13               0.993356
                                 7                   13               0.075589
             4   9.70094 eV
                                11                   13              -0.893130
                                12                   16               0.297568
                                 5                   13               0.269023
                                12                   28              -0.108047
                                 9                   15               0.066145
                                 9                   30              -0.050159
```

### Excited-state gradient for acetone

To perform an excited-state optimization for emission spectroscopy, the run type has to be set to
`GEO_OPT` and the state to be optimized has to be specified in the
[EXCITED_STATES](#CP2K_INPUT.FORCE_EVAL.DFT.EXCITED_STATES) section. Note that the number of excited
states chosen in the TDDFPT section should be larger or at least equal to the number of the chosen
excited state.

```none
 &GLOBAL
  PROJECT S20Acetone
  RUN_TYPE ENERGY_FORCE                ! The run type has to be changed to ENERGY_FORCE of GEO_OPT
  PREFERRED_DIAG_LIBRARY SL
  PRINT_LEVEL medium
 &END GLOBAL
 &PROPERTIES
    &TDDFPT
       NSTATES 10
       MAX_ITER   100
       CONVERGENCE [eV] 1.0e-7
       RKS_TRIPLETS F
       ADMM_KERNEL_CORRECTION_SYMMETRIC T   ! required keyword when using hybrid functionals and ADMM for
    &END TDDFPT                             ! the exact exchange contribution for ES gradients
 &END PROPERTIES
 ...
 &DFT
    &QS
      METHOD GPW
      EPS_DEFAULT 1.0E-17
      EPS_PGF_ORB 1.0E-20
    &END QS
    &EXCITED_STATES T
     STATE 1                     ! the excited state to be optimized has to be specified in this section
    &END EXCITED_STATES
    &SCF
      SCF_GUESS restart
      &OT
         PRECONDITIONER FULL_ALL
         MINIMIZER DIIS
      &END OT
      &OUTER_SCF
         MAX_SCF 900
         EPS_SCF 1.0E-7
      &END OUTER_SCF
      MAX_SCF 10
      EPS_SCF 1.0E-7
    &END SCF
 ...
 &END DFT
```

The resulting output file contains for each step of the geometry optimization the already explained
output of the TDDFPT computation and additionally information on the excited state that is optimized
as well as the iterative solution of the Z vector equations.

```none
 !--------------------------- Excited State Energy ----------------------------!
 Excitation Energy [Hartree]                                        0.1719188002
 Total Energy [Hartree]                                           -36.3332484083
 !-----------------------------------------------------------------------------!
 !--------------------------- Excited State Forces ----------------------------!

  Iteration    Method   Restart      Stepsize      Convergence         Time
  ------------------------------------------------------------------------------
        1        PCG       F         0.00E+00      0.0231470476        0.02
        2        PCG       F         0.30E+00      0.0003972862        1.56
        3        PCG       F         0.60E+00      0.0000361012        3.11
        4        PCG       F         0.62E+00      0.0000031865        4.65
        5        PCG       F         0.63E+00      0.0000002183        6.19
        6        PCG       F         0.73E+00      0.0000000147        7.73
        7        PCG       F         0.70E+00      0.0000000014        9.27
        8        PCG       F         0.72E+00      0.0000000001       10.82
        9        PCG       F         0.71E+00      0.0000000000       12.35
 !-----------------------------------------------------------------------------!

 ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]:              -36.333248408276823

 -------------------------------------------------------------------------------
```

### Choosing a semi-empirical kernel

To speed up computation times for broad-band absorption spectra, a semi-empirical simplified
Tamm-Dancoff (sTDA) kernel can be chosen by setting
[KERNEL](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.TDDFPT.KERNEL) to `sTDA`. The semi-empirical electron
repulsion operators depend on several empirical parameters, which need to be adjusted depending on
the system under investigation. Most importantly, the amount of exact exchange, scaled by adjusting
the parameter $a_{\rm{\tiny{EX}}}$ needs to be chosen carefully and is really crucial to ensure a
well-balanced treatment of exact exchange in the GS and ES potential energy surfaces. Too large
fractions of exchange in the excited state can lead to negative excitation energies. In general, a
relatively small amount of exchange $a_{\rm{\tiny{EX}}}=0.2 / 0.1$ is therefore recommended.

For details on the STDA method see also <https://github.com/grimme-lab/stda>.

The most important keywords of the subsection [sTDA](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.TDDFPT.STDA)
are:

- [FRACTION](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.TDDFPT.STDA.FRACTION): fraction of exact exchange
  $a_{\rm{\tiny{EX}}}$
- [MATAGA_NISHIMOTO_CEXP](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.TDDFPT.STDA.MATAGA_NISHIMOTO_CEXP):
  keyword to modify the parameter $\alpha$ of $\gamma^{\rm{\tiny{J}}}$
- [MATAGA_NISHIMOTO_XEXP](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.TDDFPT.STDA.MATAGA_NISHIMOTO_XEXP):
  keyword to modify the parameter $\beta$ of $\gamma^{\rm{\tiny{K}}}$
- [DO_EWALD](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.TDDFPT.STDA.DO_EWALD): keyword to switch on Ewald
  summation for the Coulomb contributions, required when treating periodic systems. (Exact exchange
  is treated within the minimum image convention and does not require any adjustment.)

```none
&PROPERTIES
  &TDDFPT
   KERNEL sTDA     ! switches on the semi-empirical kernel sTDA
   &sTDA
     FRACTION 0.2  ! it is crucial to adjust the fraction of exact exchange
   &END sTDA
   NSTATES 10
   MAX_ITER   100
   CONVERGENCE [eV] 1.0e-7
   RKS_TRIPLETS F
  &END TDDFPT
&END PROPERTIES
```

In the output, it can be checked that the sTDA kernel was switched on and the correct parameters for
$\alpha$, $\beta$ and $a_{\rm{\tiny{EX}}}$ were chosen:

```none
 sTDA| HFX Fraction                                                       0.2000
 sTDA| Mataga-Nishimoto exponent (C)                                      1.5160
 sTDA| Mataga-Nishimoto exponent (X)                                      0.5660
 sTDA| TD matrix filter                                           0.10000000E-09

 -------------------------------------------------------------------------------
 -                      sTDA Kernel: Create Matrix SQRT(S)                     -
 -------------------------------------------------------------------------------
```

The subsequent printout is then identical to the printout for GGA and hybrid functional kernels:

```none
 R-TDDFPT states of multiplicity 1
 Transition dipoles calculated using velocity formulation

         State    Excitation        Transition dipole (a.u.)        Oscillator
         number   energy (eV)       x           y           z     strength (a.u.)
         ------------------------------------------------------------------------
 TDDFPT|      1       4.88982   8.2927E-09  8.6085E-09 -2.2495E-08   7.77366E-17
 TDDFPT|      2       9.05888   4.5649E-09 -4.9483E-09 -8.2739E-08   1.52940E-15
 TDDFPT|      3       9.23700   1.4100E-02 -1.2802E-08  3.2962E-08   4.49914E-05
 TDDFPT|      4       9.52129   1.1297E-08 -4.7865E-08  8.3174E-01   1.61372E-01
 TDDFPT|      5      10.01790  -9.4987E-09 -4.8324E-01 -5.0387E-08   5.73140E-02
 TDDFPT|      6      10.86728  -1.6175E-08 -3.6799E-01 -7.5353E-09   3.60529E-02
 TDDFPT|      7      10.97656  -3.7898E-01 -1.0287E-08 -3.3941E-09   3.86231E-02
 TDDFPT|      8      11.21831   4.5569E-08 -4.2529E-09  1.9115E-08   6.76117E-16
 TDDFPT|      9      11.40817   2.7777E-08 -6.8033E-01 -3.7411E-08   1.29362E-01
 TDDFPT|     10      11.82525  -1.2576E-09  1.3562E-08  3.8592E-01   4.31486E-02
```

It should be noted that it is possible to combine sTDA with the semi-empirical ground-state
reference GFN1-xTB. However, it is then recommended to adjust all parameters [](#Grimme2016) and to
apply corrections to shift the virtual KS orbital eigenvalues. A shift can be applied by adding the
keyword [EV_SHIFT](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.TDDFPT.EV_SHIFT) (for open-shell systems
[EOS_SHIFT](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.TDDFPT.EOS_SHIFT)).

### Periodic systems

Adjustments are required when switching to periodic boundary conditions. As already mentioned above,
oscillator strengths should be calculated using the `BERRY` or `VELOCITY` formulation. When choosing
the `sTDA` kernel, [DO_EWALD](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.TDDFPT.STDA.DO_EWALD) has to be
switched to true to activate Ewald summation for Coulomb contributions. When computing ES gradients
using [KERNEL](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.TDDFPT.KERNEL) `FULL` in combination with hybrid
functionals and ADMM, only the ADMM2 method relying on basis projection is implemented in
combination with the default for the exchange functional for the first-order GGA correction term.

```none
&AUXILIARY_DENSITY_MATRIX_METHOD
  METHOD BASIS_PROJECTION
  ADMM_PURIFICATION_METHOD NONE
  EXCH_SCALING_MODEL NONE
  EXCH_CORRECTION_FUNC PBEX
&END
```

### Natural transition orbitals

Natural transition orbitals (NTOs) are printed when choosing `PRINT_LEVEL medium` in the `GLOBAL`
section or when enabling the
[NTO_ANALYSIS](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.TDDFPT.PRINT.NTO_ANALYSIS) section. For this
purpose it is required to generate unoccupied orbitals and the keyword
[LUMO](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.TDDFPT.NLUMO) enables to adjust the number of unoccupied
orbitals. It is possible to print the NTOs as `CUBE_FILES` or in Molden format, with the latter
being activated with the keyword
[MOS_MOLDEN](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.TDDFPT.PRINT.MOS_MOLDEN).

[hehn2022]: https://doi.org/10.1021/acs.jctc.2c00144
[strand2019]: https://doi.org/10.1063/1.5078682
