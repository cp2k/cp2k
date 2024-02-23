# X-Ray Absorption from TDDFT

This a a short tutorial on how to run near-edge X-ray absorption spectroscopy calculations using
linear-response TDDFT. The method is implemented in CP2K under the XAS_TDP name. It relies on
core-level specific approximations that enables efficient calculations of large and periodic
systems. Both K- and L-edge are available. The details of the method can be found in [](#Bussy2021).
Please cite this paper if you were to use the XAS_TDP method for work you publish.

```{note}
The XAS LR-TDDFT method comes with a correction scheme that is described in [](./correction_scheme).
```

## Brief theory recap

The method is based on 3 main core-specific approximations that boost the calculation efficiency.
The first one is the core-valence separation. Due to large differences in energy and localization,
core and valence states only weakly couple. Thus, when dealing with XAS, it is customary to simply
ignore excitations from valence states.

The second approximation is the sudden approximation, in which the relaxation of electrons beyond
the core region is neglected upon excitation of a core electron. Combined with the localized nature
of core states, this allows to treat excitations one at a time rather than all at once. This is more
efficient in the sense that diagonaling a series of small matrices scales better than diagonalizing
a single much larger one.

Finally, a lot of 4-center 2-electron integrals have to be computed. Thanks to the core-valence
separation and the sudden approximation, all required integrals involve the core state from which
the exciation takes place. This allows for a core-specific resolution of the identity scheme (RI).
For the Coulomb integrals:

$$
  (pI|Jq) \approx \sum_{\mu, \nu} \ (pI|\mu) \ (\mu|\nu)^{-1} \ (\nu|Jq)
$$

where $p,q$ represent atomic orbitals (Gaussian type orbitals/GTOs) and $I, J$ core orbtials. For
non-zero integrals, $p,I$ and $q,J$ have to overlap. Since $I,J$ are localized core orbitals
centered on the same atom, it is sufficient to take a RI basis only made GTOs centered on the
excited atom. Note that in case of K-edge spectroscopy, there is only one core state to consider and
$I=J$. For L-edge, $I,J$ span all three degenerate $2p$ states. This leads to particularly efficient
integral evaluations.

For the exchange-correlation kernel, the RI scheme reads:

$$
(pI|f_{xc}|Jq) \approx \sum_{\kappa, \lambda, \mu, \nu, } \ (pI|\kappa) \ (\kappa|\lambda)^{-1} \ (\lambda|f_{xc}|\mu) \ (\mu|\nu)^{-1} (\nu|Jq)
$$

where all integrals but $(\lambda|f_{xc}|\mu)$ are the same as for the Coulomb kernel above. Since
the RI basis elements $\lambda, \mu$ are centered on the excited atoms, we only need the density in
its vicinity. For this, we use a simple projection:

$$
n(\mathbf{r}) &=\sum_\sigma\sum_{pq} P^\sigma_{pq} \ \varphi_p(\mathbf{r}) \varphi_q(\mathbf{r})\\
%\pause
&\approx \sum_\sigma \sum_{pq}\sum_{\mu\nu} P^\sigma_{pq} \ (pq\mu) \ S_{\mu\nu}^{-1} \ \chi_\nu(\mathbf{r})\\
&= \sum_\nu d_\nu \ \chi_\nu(\mathbf{r}),
$$

which turns the density into a linear combination of RI basis elements. This allows for easy and
simple numerical integration of $(\lambda|f_{xc}|\mu)$. Note that the quality of the projection may
suffer if there are (heavy) atoms close by since their core states may not be well described (GTOs
are only sharp at their center). This can be addressed by either using pseudopotentials for the
neighbors or adding their RI basis function for the projection.

For these approximations to work, the core states to be excited need to be identified among the
Kohn-Sham oritals. They need to have a strong $1s$, $2s$ or $2p$ nature and be well localized. Not
fullfilling these conditions will lead to wrong results.

## The XAS_TDP input section

The parameters defining XAS LR-TDDFT calculations are found in the
[XAS_TDP](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP) subsection of [DFT](#CP2K_INPUT.FORCE_EVAL.DFT). Some
external parameters also need to be set to specific values. In particular,
[RUN_TYPE](#CP2K_INPUT.GLOBAL.RUN_TYPE) should be set to `ENERGY` and
[METHOD](#CP2K_INPUT.FORCE_EVAL.DFT.QS.METHOD) to `GAPW`. The combination of GAPW and all-electron
basis sets allow for an accurate description of core states.

The most important keywords and subsections of [XAS_TDP](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP) are:

- [DONOR_STATES](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP.DONOR_STATES): which define which core states
  need to be excited (and where to look for them)
- [KERNEL](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP.KERNEL): where the XC functional and exact exchange
  interaction (for hybrid TDDFT) are defined
- [GRID](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP.GRID): which defines the integration grids for the xc
  kernel $(\lambda|f_{xc}|\mu)$

The defaults value of all other keywords are in principle good enough.

Note that the first requirement for XAS LR-TDDFT is that the ground state calculation on which it is
based is of good quality.

## Simple examples

Illustrative examples usually tell more than long texts. Below, some typical input examples are
displayed with explanations. They should cover most common use cases.

### CO$_2$ molecule (K-edge)

This is a simple C and O K-edge calculation of the CO$_2$ molecule in the gas phase. The annotated
input file is displayed below:

```none
&GLOBAL
  PROJECT CO2
  RUN_TYPE ENERGY
&END GLOBAL

&FORCE_EVAL
  &DFT
    BASIS_SET_FILE_NAME EMSL_BASIS_SETS
    POTENTIAL_FILE_NAME POTENTIAL
    AUTO_BASIS RI_XAS MEDIUM              ! size of automatically generated RI basis

    &MGRID
      CUTOFF 500
      REL_CUTOFF 40
      NGRIDS 5
    &END MGRID

    &QS
      METHOD GAPW                         ! It is necesary to use the GAPW method for
    &END QS                               ! accurate description of core states

    &POISSON
      PERIODIC NONE
      PSOLVER MT
    &END

    &SCF
      EPS_SCF 1.0E-8
      MAX_SCF 30
    &END SCF

    &XC
      &XC_FUNCTIONAL
         &LIBXC
            FUNCTIONAL HYB_GGA_XC_BHandHLYP
         &END LIBXC
      &END XC_FUNCTIONAL
      &HF
         FRACTION 0.5                    ! BHandHLYP functional requires 50% exact exchange
      &END HF
    &END XC

    &XAS_TDP
      &DONOR_STATES
         DEFINE_EXCITED BY_INDEX         ! We look for states by atom index:
         ATOM_LIST 1 2                   ! we want to excite atoms 1 and 2
         STATE_TYPES 1s 1s               ! from their 1s core state.
         N_SEARCH 3                      ! The 3 lowest energy MOs need to be searched (C1s, O1s, O1s)
         LOCALIZE                        ! States need to be actively localized because O atoms are
      &END DONOR_STATES                  ! equivalent under symmetry

      GRID C 250 500                     ! Integration grid dimensions for C and O excited atoms
      GRID O 250 500                     ! there are 250 angular points (Lebedev grid) and 500
                                         ! radial points

      &KERNEL
         RI_REGION 2.0                   ! Include RI basis elements from atoms within a 2.0 Ang
                                         ! sphere radius around the excited atom for the density projection
         &XC_FUNCTIONAL
            &LIBXC
               FUNCTIONAL HYB_GGA_XC_BHandHLYP
            &END LIBXC
         &END XC_FUNCTIONAL
         &EXACT_EXCHANGE
            FRACTION 0.5                 ! Definition of the functional for the TDDFT kernel
         &END EXACT_EXCHANGE             ! Here (and usually) taken to be the same as the ground state
      &END KERNEL
    &END XAS_TDP

  &END DFT
  &SUBSYS
    &CELL
      ABC 10.0 10.0 10.0
      PERIODIC NONE
    &END CELL
    &COORD
      C 5.00 5.00 5.00
      O 5.00 5.00 6.16
      O 5.00 5.00 3.84
    &END COORD
    &KIND C
      BASIS_SET 6-311G**                ! Using all-electron basis sets and potential is necessary
      POTENTIAL ALL                     ! for the correct description of core states
    &END KIND
    &KIND O
      BASIS_SET 6-311G**
      POTENTIAL ALL
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
```

There are a few points of particular interest in the [XAS_TDP](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP)
section of this input file. When defining the
[DONOR_STATES](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP.DONOR_STATES) subsection, we specify that excited
atoms are defined by their index and proceed to list atoms 1 and 2. This is such that excitations
take place from the Carbon 1s level and from **one** of the Oxygen 1s levels. Since both O atoms are
equivalent under symmetry, there is no need to run the calculation for both of them (this
calculation is small enough and can run on a single processor, this is mostly for illustration
purposes).

In the [KERNEL](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP.KERNEL) subsection, the
[RI_REGION](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP.KERNEL.RI_REGION) keyword is set to 2.0. This is such
that all atoms within a 2.0 Angstrom radius of the current atom provide RI basis functions for the
projection of the density (see last equation of theory recap). This leads to a better quality
description, especially for the Carbon atom which has 2 heavier Oxygens around. Note that if the O
atoms were described with the use of a pseudopotential, there would be no need for
[RI_REGION](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP.KERNEL.RI_REGION).

In the resulting output file, there is a [XAS_TDP](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP) section
reporting the steps of the calculations. Especially important are the parts related to donor core
state identification. For each excited atom/core level combination, a small report about Mulliken
population analysis and overlap with pure donor level is printed. For the Carbon 1s:

```none
# Start of calculations for donor state of type 1s for atom   1 of kind C

    The following localized MO(s) have been associated with the donor state(s)
    based on the overlap with the components of a minimal STO basis:
                                             Spin   MO index     overlap(sum)
                                                1          3          1.00113

    The next best overlap for spin 1 is 0.00000 for MO with index    1

    Mulliken population analysis retricted to the associated MO(s) yields:
                                                  Spin  MO index     charge
                                                     1         3      1.002
```

Both the overlap and the Mulliken charge should be as close to 1.0 as possible. This ensure that the
molecular orbital selected is of the correct type (here projection on a C 1s Slater type orbital)
and properly localized (there is a full electron associated to this MO, on this atom). If those
numbers are lower, something went wrong with the core level identification. This is usually solved
by increasing [N_SEARCH](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP.DONOR_STATES.N_SEARCH) or/and by usinge
the [LOCALIZE](#CP2K_INPUT.FORCE_EVAL.DFT.XAS.LOCALIZE) keyword in case some atoms are equivalent
under symmetry.

Spectral information are given in a separate file named `CO2.spectrum`, where excitation energies
are listed with corresponding oscillator strengths, for each excited core level.

### Tetrahedral NaAlO$_2$ (K-edge, periodic)

This example is about crystalline sodium aluminate and illustrates how large periodic structures can
be efficiently simulated.

```none
&GLOBAL
   PROJECT  sodal
   RUN_TYPE ENERGY
   PRINT_LEVEL LOW
&END GLOBAL

&FORCE_EVAL
   METHOD QS
   &DFT
      BASIS_SET_FILE_NAME  BASIS_ADMM
      ! the pcseg-n and admm-n basis set families can be downloaded at https://www.basissetexchange.org
      BASIS_SET_FILE_NAME  BASIS_PCSEG
      BASIS_SET_FILE_NAME  BASIS_MOLOPT
      POTENTIAL_FILE_NAME  POTENTIAL
      AUTO_BASIS RI_XAS MEDIUM

      &QS
         METHOD GAPW                         ! GAPW is necessary for core states
      &END QS

      &AUXILIARY_DENSITY_MATRIX_METHOD       ! The ADMM methog greatly accelerated the ground state calculation
         ADMM_PURIFICATION_METHOD NONE       ! This is the simplest ADMM scheme and has proven to work well
      &END AUXILIARY_DENSITY_MATRIX_METHOD

      &SCF
         MAX_SCF    30
         EPS_SCF    1.0E-06

         &OT
           MINIMIZER DIIS
           PRECONDITIONER FULL_ALL
         &END OT
         &OUTER_SCF
            MAX_SCF    6
            EPS_SCF    1.0E-06
         &END
      &END SCF

      &MGRID
         CUTOFF 400
         REL_CUTOFF 40
         NGRIDS 5
      &END

      &XC
         &XC_FUNCTIONAL PBE               ! This is the PBEh functional with 45% HFX
            &PBE                          ! Large fraction of HFX are ususally needed for XAS LR-TDDFT
               SCALE_X 0.55
            &END
         &END XC_FUNCTIONAL

         &HF
            FRACTION 0.45
            &INTERACTION_POTENTIAL
               POTENTIAL_TYPE TRUNCATED   ! The tuncated Coulomb potential has to be used in PBCs
               CUTOFF_RADIUS 5.0          ! with a cutoff radius lower than half the cell size
            &END INTERACTION_POTENTIAL
            &SCREENING
               EPS_SCHWARZ 1.0E-6         ! Screening HFX integrals boosts performance
            &END SCREENING
         &END HF
      &END XC

      &XAS_TDP
         &DONOR_STATES
            DEFINE_EXCITED BY_KIND        ! We define the excited atoms by kind, which is named Alx here
            KIND_LIST Alx                 ! There is only one Alx atom in the coordinates since all Al
            STATE_TYPES 1s                ! atoms are equivalent under symmetry. The Alx atom is the only
            N_SEARCH 1                    ! one decribed at all-electron level, which is why we use
         &END DONOR_STATES                ! N_SEARCH = 1. There is also no need to LOCALIZE

         TAMM_DANCOFF                     ! TDA is turned on by default, but we make it explicit here
         GRID Alx 150 300
         ENERGY_RANGE 20.0                ! This means that we onluy solve for excitation energies that are
                                          ! up to 20.0 eV above the first energy
         &OT_SOLVER
            MINIMIZER DIIS                ! The iterative OT solver is typically more efficient than
            EPS_ITER 1.0E-4               ! full diagonalization for large systems
         &END OT_SOLVER

         &KERNEL
            &XC_FUNCTIONAL PBE
               &PBE
                  SCALE_X 0.55
               &END
            &END XC_FUNCTIONAL

            &EXACT_EXCHANGE
               OPERATOR TRUNCATED
               RANGE  5.0
               SCALE 0.45
            &END EXACT_EXCHANGE
         &END KERNEL
      &END XAS_TDP
   &END DFT

   &SUBSYS
      &CELL
         ABC 10.467947   10.651128   14.393541
      &END CELL
      &TOPOLOGY
         COORD_FILE_NAME sodal.xyz
         COORD_FILE_FORMAT xyz
      &END TOPOLOGY
    &KIND O
      BASIS_SET DZVP-MOLOPT-SR-GTH
      BASIS_SET AUX_FIT FIT3
      POTENTIAL GTH-PBE
    &END KIND
    &KIND Na
      ELEMENT Na
      BASIS_SET DZVP-MOLOPT-SR-GTH
      BASIS_SET AUX_FIT FIT3
      POTENTIAL GTH-PBE
    &END
    &KIND Al
      BASIS_SET DZVP-MOLOPT-SR-GTH
      BASIS_SET AUX_FIT FIT3
      POTENTIAL GTH-PBE
    &END
    &KIND Alx                          ! All atoms but the single Alx are described using pseudopotentials
      ELEMENT Al                       ! This greatly reduces the number of basis function and the cost of
      BASIS_SET pcseg-2                ! the calculation in general. AUX_FIT basis sets are for ADMM
      BASIS_SET AUX_FIT admm-2
      POTENTIAL ALL
    &END
   &END SUBSYS
&END FORCE_EVAL
```

There are many performance oriented keywords and subsection in the above input. Most importantly,
only one atom is treated at the all-electron level (the one atom from which the excitation takes
place), all other are described using pseudopotentials. Also quite important is the usage of the
[ADMM](#CP2K_INPUT.FORCE_EVAL.DFT.AUXILIARY_DENSITY_MATRIX_METHOD) method. This allows for very
efficient evaluation of the HFX energy in the ground state calculation. Finally, the
[OT](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.DIAGONALIZATION.OT) iterative solver is used. Since only a
handful of eigenvalues are calculated (those within 20.0 eV of the first excitation energy), this
scales much better than a full digonalization. Note that the
[RI_REGION](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP.KERNEL.RI_REGION) keyword is absent (it is set to 0
by default). Since the neighbors of the excited Al atom are described with pseudopotentials, there
is no need for extra RI basis function for the projection of the density.

This input file would generate a spectrum such as the one visible on figure 4 of [](#Bussy2021).
This is a much larger calculation than the first example though and would require a few hours on
20-30 processors (mostly to converge the SCF). In you are interested in reproducing this result,
input, geometry and pcseg-2/admm-2 basis sets are available
[here](https://www.cp2k.org/_media/howto:sodal.zip).

### TiCl$_4$ molecule (L-edge + SOC)

This example covers L-edge spectroscopy with the addition of spin-orbit coupling.

```none
&GLOBAL
  PROJECT TiCl4
  PRINT_LEVEL LOW
  RUN_TYPE ENERGY
&END GLOBAL
&FORCE_EVAL
  &DFT
    BASIS_SET_FILE_NAME  BASIS_DEF2-TZVPD
    POTENTIAL_FILE_NAME  POTENTIAL
    AUTO_BASIS RI_XAS    LARGE

    &POISSON
      PERIODIC NONE
      PSOLVER MT
    &END POISSON
    &QS
      METHOD GAPW
    &END QS

    &MGRID
      CUTOFF 800
      REL_CUTOFF 50
      NGRIDS 5
    &END

    &SCF
      EPS_SCF 1.0E-8
      MAX_SCF 200
      &MIXING
         METHOD BROYDEN_MIXING
         ALPHA 0.2
         BETA 1.5
         NBROYDEN 8
      &END MIXING
    &END SCF

    &XC
      &XC_FUNCTIONAL
         &LIBXC
            FUNCTIONAL HYB_GGA_XC_B3LYP
         &END LIBXC
      &END XC_FUNCTIONAL
      &HF
         FRACTION 0.2
      &END HF
    &END XC

    &XAS_TDP
      &DONOR_STATES
         DEFINE_EXCITED BY_KIND
         KIND_LIST Ti
         STATE_TYPES 2p             ! 2p core state for L-edge
      &END DONOR_STATES             ! No need to LOCALIZE since only one Ti atom

      TAMM_DANCOFF FALSE            ! TDA is on by default, get full TDDFT like this
      DIPOLE_FORM LENGTH

      GRID Ti 500 1000              ! This is a fairly dense grid

      EXCITATIONS RCS_SINGLET       ! For SOC calculations in closed-shell system, these 3 keywords
      EXCITATIONS RCS_TRIPLET       ! are required. Singlet and triplet excitation are coupled together
      SOC                           ! with the SOC hamiltonian

      &KERNEL
         RI_REGION 5.0              ! To get the best possible density projection
      &XC_FUNCTIONAL
         &LIBXC
            FUNCTIONAL HYB_GGA_XC_B3LYP
         &END LIBXC
      &END XC_FUNCTIONAL
         &EXACT_EXCHANGE
            FRACTION 0.2
         &END EXACT_EXCHANGE
      &END KERNEL
    &END XAS_TDP

  &END DFT
  &SUBSYS
    &KIND Cl
      BASIS_SET def2-TZVPD
      POTENTIAL ALL
      RADIAL_GRID 80                ! The GAPW grids are also used to evaluate the SOC operator
      LEBEDEV_GRID 120              ! it is good practice to use sligthly larger ones than the default
    &END KIND
    &KIND Ti
      BASIS_SET def2-TZVPD
      POTENTIAL ALL
      RADIAL_GRID 80
      LEBEDEV_GRID 120
    &END KIND
    &CELL
      ABC 10.0 10.0 10.0
      PERIODIC NONE
    &END CELL
    &TOPOLOGY
      COORD_FILE_FORMAT XYZ
      COORD_FILE_NAME TiCl4.xyz
      &CENTER_COORDINATES
      &END CENTER_COORDINATES
    &END TOPOLOGY
  &END SUBSYS
&END FORCE_EVAL
```

The structure of the input file is not very different from the CO$_2$ example. Notable differences
are the [DONOR_STATES](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP.DONOR_STATES) subsection where 2p states
are specified and the combinations of the
[EXCITATIONS](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP.EXCITATIONS) and
[SOC](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP.SPIN_ORBIT_COUPLING) keywords. Indeed, the way spin-orbit
coupling is treated in [XAS_TDP](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP) is by coupling together singlet
and triplet excitation via the ZORA SOC Hamiltionian.

Note that this calculation is meant to be a benchmark calculation, hence the overall larger
[GRID](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP.GRID) and
[CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.CUTOFF) values. The result can be seen in [](#Bussy2021),
figure 1 a). The calculation takes about 10 minutes on 4 cores. All necessary files are available
[here](https://www.cp2k.org/_media/howto:ticl4.zip).

In the output file, the donor state identification yields overlaps that are greater than one. This
is due to the degenerate nature of 2p states. The candidate Kohn-Sham orbital is projected on 3 STOs
for 2px, 2py and 2pz. To avoid cancelling contributions, the sum of the absolute overlap is taken.

```none
  # Start of calculations for donor state of type 2p for atom   1 of kind Ti

    The following canonical MO(s) have been associated with the donor state(s)
    based on the overlap with the components of a minimal STO basis:
                                             Spin   MO index     overlap(sum)
                                                1          7          1.36751
                                                1          8          1.36751
                                                1          9          0.99786

    The next best overlap for spin 1 is 0.06653 for MO with index   27

    Mulliken population analysis retricted to the associated MO(s) yields:
                                                  Spin  MO index     charge
                                                     1         7      1.000
                                                     1         8      1.000
                                                     1         9      1.000
```

## FAQ

### Which functional and basis sets to use ?

Hybrid functionals with high fraction of Hartree-Fock exchange are know to perform well for core
spectroscopy. PBEh($\alpha=0.45$) and BHandHLYP have had success with this particular
implementation. In periodic boundary conditions, the truncated Coulomb operator should be used (with
truncation radius \< half cell parameter).

For appropriate description of core states, all-electron basis sets should be used for the excited
atom(s). MOLOPT basis sets and pseudopotentials can be used on all other atoms. There exist core
specific basis such as pcX-n and cc-pCVXZ, but their usage is not necessary (based on basis set
convergence studies on small molecules). Note that the pcseg-n basis sets are nice to use as they
come with their own ADMM basis.

### How do I make my calculation more accurate ?

The first necessity is to have a good ground state calculation. Thus, any change of paramter
improving the `SCF` will reflect on the quality of the LR-TDDFT calculation.

Within [XAS_TDP](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP), a few parameters may play a role:

- [EPS_FILTER](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP.EPS_FILTER) and
  [EPS_PGF_XAS](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP.EPS_PGF_XAS) are used for screening. Lowering
  those will result in slower but more accurate calculations
- Increasing the [GRID](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP.GRID) dimensions will improve the quality
  of the numerical integration of the XC kernel. The upper limit for the number of angular points is
  974\. There is not upper limit for the radial point, but performance may suffer if too large.
  Usually, something like `GRID C 250 500` is sufficient.
- Increasing the [RI_REGION](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP.KERNEL.RI_REGION) in the
  [KERNEL](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP.KERNEL) subsection will lead to a more accurate
  projection of the density on the RI basis. Basis functions centered on atoms within the region
  (defined by a sphere around the excited atom) are added for the projection. Increasing this
  parameter should come together with denser [GRID](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP.GRID).

By default, the RI basis used for the integral and the projection is autamatically generated. The
quality of the RI basis can be changed via the [AUTO_BASIS](#CP2K_INPUT.FORCE_EVAL.DFT.AUTO_BASIS)
keyword in the [DFT](#CP2K_INPUT.FORCE_EVAL.DFT) section. To improve from the default `MEDIUM` size,
one can used: `AUTO_BASIS RI_XAS LARGE/HUGE`. Note that an external RI basis set can also be
provided.

### How do I make my calculation faster ?

All points mentioned above in the accuracy section can also be tweaked for performence. In general,
lowering accuracy will lead to faster calculations.

For large systems, it is recommanded to use the
[OT](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.DIAGONALIZATION.OT) iterative solver rather than the default
full digonalization. This should improve the scaling of the method. See the NaAlO$_2$ example.

The Tamm-Dancoff approximation is well established and generally yields results as good as full
TDDFT. It is moreover much cheaper than the latter. It is turned on by default, but you may want to
make sure it is enabled.

The use of [ADMM](#CP2K_INPUT.FORCE_EVAL.DFT.AUXILIARY_DENSITY_MATRIX_METHOD) is highly recommanded
for large systems, where the ground state HFX evaluation is the main bottleneck. It is also
recommanded to use pseudopotentials on all atoms that are not excited as all-electron basis set tend
to be large. If there exist no proper ADMM basis for the all-electron basis used for the excited
atom, you may use the full basis as `AUX_FIT`. If the latter is very diffuse, it may be beneficial
to remove the most diffuse elements.

The code is also both MPI and OMP parallelized. Using more core will, to a certain degree, speedup
your calculations as well.

### How do I plot a spectrum from the \*.spectrum output file

For each donor state in the system, the \*.spectrum file contains a list of excitation energies and
corresponding oscillator strengths. This yields a stick spectrum which needs to be artificially
broadened to match experiments. This is typically done using Gaussian of Lorentzian functions. Note
that in case of spin-orbit calculation at the L-edge, results for the singlet, triplet and SOC
excitations are given.

Remember that XAS LR-TDDFT produces an accurate spectrum, but it is usually wrongly positioned on
the energy axis. Again, to match experiment, a rigid shift must be applied to the result.

### My calculation yields a wrong/unphysical result, what do I do ?

Assuming that the ground state calculation is sound and well converged, there are two main causes
for failure.

The proper core state was not found/identified. In the [XAS_TDP](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP)
part of the output file, look for the Mulliken population analysis and the overlap with a STO basis.
Both quantities should close to 1.0. If they are not, there is a problem with the selected core
state. You may want to increase the value of `N_SEARCH` in
[DONOR_STATES](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP.DONOR_STATES) to scan additional Kohn-Sham
orbitals. If there are multiple atoms that are equivalent under symmetry, make sure to use the
`LOCALIZE` keyword of [DONOR_STATES](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP.DONOR_STATES) as well.

The numerical integration of the XC kernel $(\lambda|f_{xc}|\mu)$ is not accurate enough. In this
case, you may want to increase the density of the integration grid with the
[GRID](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP.GRID) keyword, increase the
[RI_REGION](#CP2K_INPUT.FORCE_EVAL.DFT.XAS_TDP.KERNEL.RI_REGION) and/or increase the quality of the
genenerated RI basis set. Note that using an external RI basis set may also help as the basis
generation scheme may fail (For example: def2-QZVP for Zn, use def2-QZVP-RIFIT in this case).

Finally, keep in mind that calculated spectra need to be rigidly shifted by some energy to match
experiment.
