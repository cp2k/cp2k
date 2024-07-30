# X-Ray Ab-Initio Correction Scheme

As mentioned in [](./tddft), XAS LR-TDDFT results need to be rigidly shifted to match experiments.
This is due to self-interaction error and the lack of orbital relaxation upon the creation of the
core hole. An *ab-initio* correction scheme was developed to address these issues. Theory and
benchmarks were published in [](#Bussy2021b). Please cite this paper if you were to use this method.

## Brief theory recap

XAS LR-TDDFT yield excitation energies as correction to ground state Kohn-Sham orbital energy
differences, namely:

$$
\omega = \varepsilon_a - \varepsilon_I + \Delta_{xc},
$$

where $\varepsilon_a$ is the orbital energy of a virtual MO and $\varepsilon_I$ the energy of the
donor core MO. Under Koopman's condition, these energies are interpreted as the electron affinity
and and the ionization potential (IP). However, DFT is notoriously bad at predicting accurate
absolute orbital eigenvalues. Therefore, and because $|\varepsilon_I| >> |\varepsilon_a|$,
excitation energies are expected to be widely improved if the DFT energy $\varepsilon_I$ were to be
replace by an accurate value of the IP.

The IP can be accurately calculated using the second-order electron propagator equation:

$$
\text{IP}_I = -\varepsilon_I - \frac{1}{2} \sum_{ajk}\frac{|\langle Ia||jk\rangle|^2}{-\text{IP}_I + \varepsilon_a -\varepsilon_j -\varepsilon_k} - \frac{1}{2}\sum_{abj}\frac{|\langle Ij||ab\rangle|^2}{-\text{IP}_I + \varepsilon_j - \varepsilon_a - \varepsilon_b}
$$

where $a, b$ refer to virtual Hartree-Fock spin-orbitals and $j,k$ to occupied HF spin-orbitals. The
DFT generalization of this theory is known as GW2X ([](#Shigeta2001)). It involves calculating the
Generalized Fock matrix and the rotation of the occupied and virtual DFT orbitals separately, such
that they become pseudocanonical. Alternatively, the diagonal elements of the generalized Fock
matrix can be used as approximations for the orbital energies (thus saving on the orbital rotation).
This is known as the GW2X\* method.

## The GW2X input subsection

The parameters defining the GW2X correction to XAS LR-TDDFT are found in the
[GW2X](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP.GW2X) subsection of
[XAS_TDP](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP). GW2X will only work with hybrid functionals (or full
Hartree-Fock), as the machinery necessary for the calculation of the generalized Fock matrix is not
available otherwise.

There are not many parameters to set for the GW2X correction. Simply adding an empty
[GW2X](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP.GW2X) subsection is usually enough. The electron
propagator equation is solved iteratively with a Newton-Raphson scheme.
[EPS_GW2X](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP.GW2X.EPS_GW2X) controls the convergence threshold and
[MAX_GW2X_ITER](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP.GW2X.MAX_GW2X_ITER) the maximum number of
iterations allowed. The [PSEUDO_CANONICAL](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP.GW2X.PSEUDO_CANONICAL)
keyword controls whether the original GW2X scheme or its simplified GW2X\* version is run (by
default, the original GW2X is on). [C_SS](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP.GW2X.C_SS) and
[C_OS](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP.GW2X.C_OS) allow to scale the same- and opposite-spin
components (as in SOS- and SCS-MP2). Finally, if
[XPS_ONLY](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP.GW2X.XPS_ONLY) is set, only the core IP is calculated
and the XAS LR-TDDFT calculation is skipped.

## Simple examples

### OCS molecule (L-edge + SOC)

This example covers GW2X corrected L-edge spectroscopy with spin-orbit coupling.

```none
&GLOBAL
  PROJECT OCS
  PRINT_LEVEL MEDIUM
  RUN_TYPE ENERGY
&END GLOBAL
&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME BASIS_GW2X
    POTENTIAL_FILE_NAME POTENTIAL
    AUTO_BASIS RI_XAS MEDIUM

    &MGRID
      CUTOFF 800
      REL_CUTOFF 50
      NGRIDS 5
    &END MGRID
    &QS
      METHOD GAPW
    &END QS

    &POISSON
      PERIODIC NONE
      PSOLVER MT
    &END

    &SCF
      EPS_SCF 1.0E-8
      MAX_SCF 50
    &END SCF

    &XC
      &XC_FUNCTIONAL                ! The PBEh(45%) functional
         &LIBXC
            FUNCTIONAL GGA_C_PBE
         &END LIBXC
         &LIBXC
            FUNCTIONAL GGA_X_PBE
            SCALE 0.55
         &END LIBXC
      &END XC_FUNCTIONAL

      &HF
         FRACTION 0.45
      &END HF
    &END XC

    &XAS_TDP
      &DONOR_STATES
         DEFINE_EXCITED BY_KIND
         KIND_LIST S
         STATE_TYPES 2p          ! Need to look for the S 2p states within the 7 MOs with lowest energy;
         N_SEARCH 7              ! one S 1s, one S 2s, three S 2s , one C 1s and one O 1s
         LOCALIZE                ! Localization is required
      &END DONOR_STATES

      EXCITATIONS RCS_SINGLET
      EXCITATIONS RCS_TRIPLET
      SOC

      GRID S  300 500

      N_EXCITED 150
      TAMM_DANCOFF

      &GW2X                      ! This is the only difference in the input file with respect to a
      &END GW2X                  ! standard XAS_TDP calculation (defaults parameters are used)

      &KERNEL
         RI_REGION 3.0
         &XC_FUNCTIONAL
            &LIBXC
               FUNCTIONAL GGA_C_PBE
            &END LIBXC
            &LIBXC
               FUNCTIONAL GGA_X_PBE
               SCALE 0.55
            &END LIBXC
         &END XC_FUNCTIONAL
         &EXACT_EXCHANGE
            FRACTION 0.45
         &END EXACT_EXCHANGE
      &END KERNEL

    &END XAS_TDP
  &END DFT
  &SUBSYS
    &CELL
      ABC 10.0 10.0 10.0
      PERIODIC NONE
    &END CELL
    &COORD
      C         5.0000000209        4.9999999724        5.2021372095
      O         5.0000000094        5.0000000290        6.3579624316
      S         5.0000000207        5.0000000007        3.6399034216
    &END COORD
    &KIND C
      BASIS_SET aug-pcX-2
      POTENTIAL ALL
    &END KIND
    &KIND O
      BASIS_SET aug-pcX-2
      POTENTIAL ALL
    &END KIND
    &KIND S
      BASIS_SET aug-pcX-2
      POTENTIAL ALL
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
```

The only difference between the above input file and that of a standard XAS LR-TDDFT calculation is
the addition of the [GW2X](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP.GW2X) subsection. In this case, only
default parameters are used, which corresponds to the original GW2X scheme with a convergence
threshold of 0.01 eV. Note that the core specific all-electron aug-pcX-2 basis set is used (triple
zeta quality). This inputs corresponds to an entry of table II in [](#Bussy2021b), although slacker
parameters are used here (in order to make this tutorial cheap and easy to run, this particular
calculations takes ~2 minutes on 4 cores).

In the output file, the correction for each S $2p$ is displayed. Note that the correction amounts to
a shift of 1.9 eV compared to standard XAS LR-TDDFT, leading to a first singlet excitation energy of
164.4 eV (at the L$_3$ edge). This fits
[experimental results](<https://doi.org/10.1016/s0301-0104(97)00111-0>) within 0.1 eV. thus clearly
improving the XAS LR-TDDFT result. Note that the core IPs, including spin-orbit coupling effects,
are also provided. These can be directly used to produce a XPS spectrum. The content of the
`OCS.spectrum` file yields the corrected spectrum directly.

```none
    - GW2X correction for donor MO with spin  1 and MO index    5:
                             iteration                convergence (eV)
                                     1                       10.047536
                                     2                        1.237503
                                     3                        0.014416
                                     4                       -0.000000

      Final GW2X shift for this donor MO (eV):   1.927146


    - GW2X correction for donor MO with spin  1 and MO index    6:
                             iteration                convergence (eV)
                                     1                        6.197650
                                     2                        4.963008
                                     3                        0.241838
                                     4                        0.000439

      Final GW2X shift for this donor MO (eV):   1.907648


    - GW2X correction for donor MO with spin  1 and MO index    7:
                             iteration                convergence (eV)
                                     1                        6.197650
                                     2                        4.963008
                                     3                        0.241838
                                     4                        0.000439

      Final GW2X shift for this donor MO (eV):   1.907648


    Calculations done:

    First singlet XAS excitation energy (eV):                165.014087
    First triplet XAS excitation energy (eV):                164.681850
    First SOC XAS excitation energy (eV):                    164.396537

    Ionization potentials for XPS (GW2X + SOC):              170.602279
                                                             169.457339
                                                             169.367465

```

### Solid NH$_3$ (K-edge, periodic)

This is a much larger example of a periodic system, namely solid ammonia. This example is much
heavier to run (~45 minutes on 24 cores).

```none

&GLOBAL
  PROJECT NH3
  RUN_TYPE ENERGY
  PRINT_LEVEL MEDIUM
&END GLOBAL
&FORCE_EVAL
  METHOD QS
  &DFT
    BASIS_SET_FILE_NAME BASIS_GW2X
    BASIS_SET_FILE_NAME BASIS_ADMM
    BASIS_SET_FILE_NAME BASIS_MOLOPT
    POTENTIAL_FILE_NAME POTENTIAL
    AUTO_BASIS RI_XAS MEDIUM

    &QS
      METHOD GAPW
    &END QS

    &MGRID
      CUTOFF 600
      REL_CUTOFF 50
      NGRIDS 5
    &END MGRID

    &SCF
      SCF_GUESS RESTART
      EPS_SCF 1.0E-8
      MAX_SCF 30

      &OT
         MINIMIZER CG
         PRECONDITIONER FULL_ALL
      &END OT

      &OUTER_SCF
         MAX_SCF 6
         EPS_SCF 1.0E-8
      &END OUTER_SCF

    &END SCF

    &AUXILIARY_DENSITY_MATRIX_METHOD
      ADMM_PURIFICATION_METHOD NONE
    &END AUXILIARY_DENSITY_MATRIX_METHOD

    &XC
      &XC_FUNCTIONAL
         &LIBXC
            FUNCTIONAL GGA_X_PBE
            SCALE 0.55
         &END
         &LIBXC
            FUNCTIONAL GGA_C_PBE
         &END
      &END XC_FUNCTIONAL
      &HF
         FRACTION 0.45
         &INTERACTION_POTENTIAL
            POTENTIAL_TYPE TRUNCATED
            CUTOFF_RADIUS 5.0
         &END INTERACTION_POTENTIAL
      &END HF
    &END XC

    &XAS_TDP
      &DONOR_STATES
         DEFINE_EXCITED BY_KIND
         KIND_LIST Nx
         STATE_TYPES 1s
         N_SEARCH 1
         LOCALIZE
      &END DONOR_STATES

      TAMM_DANCOFF
      GRID Nx 300 500
      E_RANGE 30.0

      &GW2X
      &END

      &KERNEL
         &XC_FUNCTIONAL
            &LIBXC
               FUNCTIONAL GGA_X_PBE
               SCALE 0.55
            &END
            &LIBXC
               FUNCTIONAL GGA_C_PBE
            &END
         &END XC_FUNCTIONAL
         &EXACT_EXCHANGE
            OPERATOR TRUNCATED
            CUTOFF_RADIUS 5.0
            FRACTION 0.45
         &END EXACT_EXCHANGE
      &END KERNEL

    &END XAS_TDP
  &END DFT
  &SUBSYS
    &CELL
      ABC   10.016118  10.016118  10.016118
    &END CELL
    &TOPOLOGY
      COORD_FILE_FORMAT XYZ
      COORD_FILE_NAME NH3.xyz
    &END TOPOLOGY
    &KIND H
      BASIS_SET DZVP-MOLOPT-SR-GTH
      BASIS_SET AUX_FIT FIT3
      POTENTIAL GTH-PBE
    &END KIND
    &KIND N
      BASIS_SET DZVP-MOLOPT-SR-GTH
      BASIS_SET AUX_FIT FIT3
      POTENTIAL GTH-PBE
    &END KIND
    &KIND Nx
      ELEMENT N
      BASIS_SET aug-pcseg-2
      BASIS_SET AUX_FIT aug-admm-2
      POTENTIAL ALL
    &END KIND
  &END SUBSYS
&END FORCE_EVAL

```

Again, the only difference with respect to a standard XAS-LRTDDFT input file is the
[GW2X](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP.GW2X) subsection. This input file corresponds exactly to
figure 3 a) of [](#Bussy2021b). In this case, the GW2X correction amounts to a blue shift of 3.7 eV,
aligning the calculated spectrum to the experimental one remarkably well. The NH3 structure file as
well as the necessary basis set file (also for the OCS example) are available
[here](https://github.com/cp2k/cp2k-examples/tree/master/x-ray/gw2x).

## FAQ

### How can I make the GW2X correction run faster ?

The GW2X correction scheme scales cubically with the number of MOs in the system. Therefore, the
best way to improve performance is to reduce that number. Because an accurate description of the
core region is only necessary for the exited atoms, all other atoms can be described with
pseudopotentials. This drastically reduces the number of MOs since only valence states are kept. In
the solid NH3 example above, all nitrogen atoms are equivalent under symmetry. Therefore, their
individual contribution to the XAS spectrum is bound to be the same. This allows for the description
of a single nitrogen atom at the all-electron level, while all others (and the hydrogens) use
pseudopotentials. Note that the [ADMM](#CP2K_INPUT.FORCE_EVAL.DFT.AUXILIARY_DENSITY_MATRIX_METHOD)
approximation is also utilized. This greatly reduces the cost of the underlying hybrid DFT
calculation, as well as the evaluation of the generalized Fock matrix as required by GW2X.

### Why don't I get the absolute core IP in periodic systems ?

For molecules in non-periodic boundary conditions, the potential is such that it is zero far away.
In the periodic case, the zero is ill defined. As a consequence, all Kohn-Sham eigenvalues end up
shifted by some unknown, constant amount. Therefore, their absolute values and that of the
calculated IP cannot be interpreted in a physical manner. However, the correction scheme depends on
the difference $|\varepsilon_a-\varepsilon_I|$, where the shift cancels out.

### Why is the LOCALIZE keyword required ?

In order to efficiently evaluate the antisymmetric integrals of the type $\langle Ia || jk \rangle$,
the same local RI scheme as XAS_TDP is used. Therefore, the core state $I$ needs to be local in
space. However, the rotation required to get the pseudocanonical orbitals needed for the original
GW2X scheme may break this localization, provided that there are other equivalent atoms in the
system. To prevent that from happening, all core states localized on other atoms are ignored for the
rotation and the subsequent IP calculation. This has negligible impact since core states belonging
to different atoms only weakly interact. It is however important to keep the value of the
[LOCALIZE](#CP2K_INPUT.FORCE_EVAL.DFT.XAS.LOCALIZE) keyword to a minimum to insure that only core
states are ignored.
