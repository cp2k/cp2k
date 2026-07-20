# Extended Tight Binding

This is a short tutorial on how to run GFN-xTB simulations with CP2K. The details on the theory and
the original implementation can be found in [](#Grimme2017) for GFN1-xTB and in [](#Bannwarth2019)
for GFN2-xTB.

The GFN-xTB methods are implemented either natively in CP2K or through the tblite library as
described in [](#Katbashev2025) and
[Alizadeh2026](https://chemrxiv.org/doi/full/10.26434/chemrxiv.15003939). The natively available
methods are GFN0-xTB and GFN1-xTB, while the tblite library provides GFN1-xTB, IPEA1-xTB and
GFN2-xTB.

In this tutorial after a brief theory recap, there is an example on how to run GFN2-xTB simulation
with tblite and later there is an example on how to run GFN1-xTB with the native implementation.

## Brief theory recap

The semi-empirical GFN1-xTB energy expression comprises contributions due to electronic (EL),
isotropic electrostatic (IES), atom-pairwise repulsion (REP), dispersion (DISP), and halogen-bonding
(XB) terms,

$$
E_{\rm{\tiny{GFN1-xTB}}} = E_{\rm{\tiny{EL}}} + E_{\rm{\tiny{IES}}} + E_{\rm{\tiny{REP}}} + E_{\rm{\tiny{DISP}}} + E_{\rm{\tiny{XB}}} \, .
$$

The GFN2-xTB, is defined similarly using anisotropic electrostatic (AES) in addition to isotropic
contributions,

$$
E_{\rm{\tiny{GFN2-xTB}}} = E_{\rm{\tiny{EL}}} + E_{\rm{\tiny{AES}}} + E_{\rm{\tiny{REP}}} + E_{\rm{\tiny{DISP}}} \, .
$$

### 1. Electronic energy

The electronic energy contribution,

$$
E_{\rm{\tiny{EL}}} =  \sum_i^{\rm{\tiny{occ}}} n_i \langle \Psi_i | h_0 | \Psi_i \rangle - T_{\rm{\tiny{el}}} S_{\rm{\tiny{el}}} \, ,
$$

contains zeroth-order contributions based on a zeroth-order Hamiltonian $h_0$, the valence molecular
orbitals $\Psi_i$, occupation numbers $n_i$ as well as the electronic temperature times entropy term
$T_{\rm{\tiny{el}}}S_{\rm{\tiny{el}}}$ from fractional orbital occupations.

### 2. Electrostatic energy

The isotropic electrostatic energy contains a second-order contributions which are optimized
self-consistently as well as third-order diagonal contributions,

$$
E_{\rm{\tiny{IES}}} = \frac{1}{2} \sum_{A,B} \sum_{{l}^A}\sum_{{l'}^B} p_l^A p_{{l'}}^B \gamma_{AB,ll'} + \frac{1}{3}\sum_{A} \Gamma_A q_A^3 \, .
$$

The second order contributions are described using the semi-empirical electron repulsion operator
$\gamma_{AB,ll'}$ which depends on the interatomic distance of atoms $A$ and $B$ as well as further
empirical parameters that are specific for different angular momenta $l$ and $l'$. The monopole
charges of the second-order expression are optimized self-consistently,

$$
p_l^A = p_l^{A_0} - \sum_{\nu}^{N_{\rm{\tiny{AO}}}} \sum_{\mu \in A, \mu \in l} S_{\mu \nu } P_{\mu \nu} \, ,
$$

referring to the atomic orbital overlap matrix $\mathbf{S}$ and the density matrix $\mathbf{P}$.

The remaining diagonal terms represent a cubic charge correction based on the Mulliken charge $q_A$
of atom $A$ and the charge derivative $\Gamma_A$ of the atomic Hubbard parameter $\eta_A$.

In GFN2-xTB the electrostatic energy is extended to anisotropic electrostatic contributions,

$$
E_{\rm{\tiny{AES}}} = E_{\rm{\tiny{IES}}}
+ \sum_{A,B}\sum_{i} \mu_i q_B T_{AB,i}
+ \frac12\sum_{A,B}\sum_{ij} \mu_{A,i} \mu_{B,j} T_{AB,ij}
+ \sum_{A,B}\sum_{ij} \theta_{A,ij} q_B T_{AB,ij}
\, ,
$$

including the atomic dipole moments $\mu_A$ and atomic quadrupole moments $\theta_A$. The
interaction tensors $T_{AB}$ are defined by

$$
T_{AB,i} = -f(R_{AB})\frac{R_{AB,i}}{R^3_{AB}}
\, ,
$$

for the charge-dipole interaction and

$$
T_{AB,ij} = f(R_{AB})\frac{3R_{AB,i}R_{AB,j} - R^2_{AB}}{R^5_{AB}}
\, ,
$$

for the dipole-dipole and charge-dipole interaction, both interaction tensors include a short-ranged
damping function $f$.

### 3. Repulsion

Repulsion is described via an atom-pairwise potential,

$$
E_{\rm{\tiny{REP}}} = \sum_{AB} \frac{Z_A^{\rm{\tiny{eff}}} Z_B^{\rm{\tiny{eff}}} }{R_{AB}} \exp^{- (\alpha_A \alpha_B)^{1/2} (R_{AB})^{k_f}} \, ,
$$

with the effective nuclear charge $\mathbf{Z}^{\rm{\tiny{eff}}}$ as well as the global or
element-specific parameters $k_f$ and $\alpha$.

### 4. Dispersion

Dispersion is included by the well-established D3 method in the BJ-damping scheme [](#Grimme2010)
for GFN1-xTB and using the D4 method with self-consistent charges for GFN2-xTB.

### 5. Corrections

Corrections for element-specific interactions are possible using either a halogen-bonding correction
term (XB) or a generic nonbonding potential correction (NONBOND). Note that the generic nonbonding
potential correction is CP2K specific and thus the so-obtained energy differs from the original
GFN1-xTB method,

$$
E_{\rm{\tiny{GFN1-xTB+NONBOND}}} = E_{\rm{\tiny{GFN1-xTB}}} + E_{\rm{\tiny{NONBOND}}} \, .
$$

## The GFN-xTB input section

In this section the keywords for the native GFN-xTB implementation are described. The most important
keywords and subsections of section [XTB](#CP2K_INPUT.FORCE_EVAL.DFT.QS.XTB) are:

- [CHECK_ATOMIC_CHARGES](#CP2K_INPUT.FORCE_EVAL.DFT.QS.XTB.CHECK_ATOMIC_CHARGES): the cubic charge
  diagonal contribution is checked to be numerically stable by switching the keyword to true.

- [USE_HALOGEN_CORRECTION](#CP2K_INPUT.FORCE_EVAL.DFT.QS.XTB.USE_HALOGEN_CORRECTION): keyword to
  switch off contribution $E_{\rm{\tiny{XB}}}$ to correct halogen interactions, default is to
  include this correction

- [DO_NONBONDED](#CP2K_INPUT.FORCE_EVAL.DFT.QS.XTB.DO_NONBONDED): add a generic correction potential
  to correct bond- or atomic-specific interactions

- [PARAMETER](#CP2K_INPUT.FORCE_EVAL.DFT.QS.XTB.PARAMETER): it is possible to add this section with
  corresponding keywords to modify xTB parameters

- [DO_EWALD](#CP2K_INPUT.FORCE_EVAL.DFT.QS.XTB.DO_EWALD): keyword to activate Ewald summation for
  periodic boundary conditions (PBC); has to be switched to true in case of PBC. Starting from CP2K
  2026.2 the periodicity is directly taken from the cell definition, which makes this keyword
  obsolete.

The additional keywords
[COULOMB_INTERACTION](#CP2K_INPUT.FORCE_EVAL.DFT.QS.XTB.COULOMB_INTERACTION),
[COULOMB_LR](#CP2K_INPUT.FORCE_EVAL.DFT.QS.XTB.COULOMB_LR) and
[TB3_INTERACTION](#CP2K_INPUT.FORCE_EVAL.DFT.QS.XTB.TB3_INTERACTION) are for debugging purposes only
and it is recommended to use the default options here.

## Input section for tblite

Available methods in the [METHOD](#CP2K_INPUT.FORCE_EVAL.DFT.QS.XTB.TBLITE.METHOD) keyword are

- GFN1: for GFN1-xTB method
- GFN2: for GFN2-xTB method
- IPEA1: for IPEA1-xTB version of GFN1-xTB
- PARAM: read parameter file from filename provided via
  [PARAM](#CP2K_INPUT.FORCE_EVAL.DFT.QS.XTB.TBLITE.PARAM) keyword, the format for the parameter file
  is following tblite structure as described in the
  [documentation](https://tblite.readthedocs.io/en/latest/spec/parameter.html).

For enabling GFN-xTB methods via the tblite library, set the
[METHOD](#CP2K_INPUT.FORCE_EVAL.DFT.QS.METHOD) keyword to xTB and the
[GFN_TYPE](#CP2K_INPUT.FORCE_EVAL.DFT.QS.XTB.GFN_TYPE) in the
[XTB](#CP2K_INPUT.FORCE_EVAL.DFT.QS.XTB) block to tblite. The
[TBLITE](#CP2K_INPUT.FORCE_EVAL.DFT.QS.XTB.TBLITE) section allows to configure the input for the
library, like selecting the actual GFN-xTB method via the
[METHOD](#CP2K_INPUT.FORCE_EVAL.DFT.QS.XTB.TBLITE.METHOD) keyword.

```
&QS
  METHOD xTB
  &XTB
    GFN_TYPE TBLITE
    SCC_MIXER AUTO
    &TBLITE
      METHOD GFN2
      ACCURACY 1.0
    &END TBLITE
  &END XTB
&END QS
```

The [ACCURACY](#CP2K_INPUT.FORCE_EVAL.DFT.QS.XTB.TBLITE.ACCURACY) keyword controls the convergence
thresholds of the electronic mixer in the tblite library. Lower values correspond to tighter
convergence.

## Simple examples

### GFN-xTB ground-state energy based on Tblite

The following input is an example for calculating a single point calculation of ice-Ih crystal with
GFN2-xTB method. Please note that k-points are fully supported for tblite in CP2K.

```
&GLOBAL
  PRINT_LEVEL LOW
  PROJECT ice_Ih_GFN2_k333
  RUN_TYPE ENERGY
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
  &DFT
    &QS
      EPS_DEFAULT 1.0E-12
      METHOD xTB
      &XTB
        GFN_TYPE TBLITE
        &TBLITE
          METHOD GFN2
          ACCURACY 0.1
        &END TBLITE
      &END XTB
    &END QS
    &KPOINTS
      SCHEME MACDONALD 3 3 3 0.0 0.0 0.0
      FULL_GRID T
    &END KPOINTS
    &SCF
      EPS_SCF 1.0E-9
      MAX_SCF 300
      SCF_GUESS MOPAC
      &MIXING
        METHOD DIRECT_P_MIXING
        ALPHA 0.2
      &END MIXING
      &PRINT
        &RESTART OFF
        &END RESTART
      &END PRINT
    &END SCF
  &END DFT
  &SUBSYS
    &CELL
      PERIODIC XYZ
      A 7.678093000000 0.000000000000 0.000000000000
      B 3.839046000000 6.649423000000 0.000000000000
      C 0.000000000000 0.000000000000 7.234567000000
    &END CELL
    &COORD
      SCALED
      H   0.000007000000  0.334718000000  0.199799000000
      H   0.665252000000  0.000010000000  0.199780000000
      H   0.334714000000  0.665261000000  0.199791000000
      H   0.334760000000  0.999976000000  0.699786000000
      H   0.999980000000  0.665239000000  0.699800000000
      H   0.665240000000  0.334755000000  0.699799000000
      H   0.544461000000  0.000006000000  0.019584000000
      H  -0.000011000000  0.455520000000  0.019605000000
      H   0.455505000000  0.544480000000  0.019594000000
      H   0.455543000000  0.999980000000  0.519584000000
      H   0.999998000000  0.544446000000  0.519600000000
      H   0.544461000000  0.455524000000  0.519592000000
      H   0.332240000000  0.879481000000  0.984486000000
      H   0.211731000000  1.120507000000  0.984481000000
      H   0.879462000000  0.788265000000  0.984489000000
      H   0.788248000000  0.332255000000  0.984498000000
      H   0.667731000000  0.211748000000  0.984486000000
      H   0.120526000000  0.211703000000  0.484497000000
      H   0.667762000000  1.120506000000  0.484488000000
      H   0.788272000000  0.879477000000  0.484480000000
      H   0.211722000000  0.667743000000  0.484495000000
      H   0.332240000000  0.788252000000  0.484486000000
      H  -0.120503000000  0.332222000000  0.484499000000
      H   1.120486000000  0.667750000000  0.984491000000
      O   0.000009000000  0.331330000000  0.061616000000
      O   0.668637000000  0.000012000000  0.061595000000
      O   0.331326000000  0.668660000000  0.061608000000
      O   0.331370000000  0.999972000000  0.561601000000
      O   0.668637000000  0.331348000000  0.561616000000
      O   0.335676000000  0.999992000000  0.936746000000
      O   0.664294000000  0.335703000000  0.936768000000
      O   0.000012000000  0.335655000000  0.436771000000
      O   0.664333000000  0.999995000000  0.436740000000
      O   0.335670000000  0.664304000000  0.436760000000
      O   0.999979000000  0.668634000000  0.561617000000
      O   0.999974000000  0.664306000000  0.942299010000
    &END COORD
  &END SUBSYS
&END FORCE_EVAL
```

### Unrestricted calculation with spGFN2-xTB

In case of open-shell calculations, a spin-polarization term can be enabled with the
[LSD](#CP2K_INPUT.FORCE_EVAL.DFT.UKS) keyword in CP2K. In this case, tblite automatically allows the
usage of spGFN2-xTB for calculations as described in [Neugebauer2023](#Neugebauer2023). An example
for triplet oxygen is shown here.

```
&FORCE_EVAL
  &DFT
    LSD
    MULTIPLICITY 3
    &QS
      EPS_DEFAULT 1.00E-12
      METHOD xTB
      &XTB
        GFN_TYPE TBLITE
        SCC_MIXER TBLITE
        &TBLITE
          METHOD GFN2
        &END TBLITE
      &END XTB
    &END QS
    &SCF
      ADDED_MOS 1 3
      EPS_SCF 1.e-10
      MAX_SCF 200
      SCF_GUESS MOPAC
      &PRINT
        &RESTART OFF
        &END RESTART
      &END PRINT
    &END SCF
  &END DFT
  &SUBSYS
    &CELL
      ABC 20.0 20.0 20.0
      PERIODIC NONE
    &END CELL
    &COORD
      O     0.000000    0.000000    0.000000
      O     1.208000    0.000000    0.000000
    &END COORD
  &END SUBSYS
&END FORCE_EVAL
```

### GFN1-xTB ground-state energy based on native Quickstep

The following input is a standard example for calculating GFN1-xTB ground-state energies using the
native quickstep implementation [XTB](#CP2K_INPUT.FORCE_EVAL.DFT.QS.XTB).

```none
&GLOBAL
  RUN_TYPE  ENERGY
  PROJECT_NAME xtb
  PRINT_LEVEL  MEDIUM
  PREFERRED_DIAG_LIBRARY SL
&END GLOBAL
&FORCE_EVAL
 METHOD QS
&DFT
  &QS
   METHOD XTB
   &XTB
    CHECK_ATOMIC_CHARGES F    ! Keyword to check if Mulliken charges are physically reasonable
    DO_EWALD  T               ! Ewald summation is required for periodic structures
    USE_HALOGEN_CORRECTION T  ! Element-specific correction for halogen interactions (Cl, Br) with (O, N)
   &END XTB
  &END QS
  &SCF
   SCF_GUESS RESTART
   MAX_SCF 50
   EPS_SCF 1.E-6
   &OT ON
     PRECONDITIONER FULL_SINGLE_INVERSE
     MINIMIZER DIIS
   &END
   &OUTER_SCF
     MAX_SCF 200
     EPS_SCF 1.E-6
   &END OUTER_SCF
  &END SCF
 &END DFT
  &SUBSYS
    &TOPOLOGY
      COORD_FILE_FORMAT  xyz
      COORD_FILE_NAME  input.xyz
      CONNECTIVITY OFF
      &CENTER_COORDINATES
      &END CENTER_COORDINATES
     &END TOPOLOGY
    &CELL
      ABC  21.64 21.64 21.64
      ALPHA_BETA_GAMMA 90.0 90.0 90.0
      PERIODIC XYZ
    &END CELL
  &END SUBSYS
&END FORCE_EVAL
```

The so-obtained output is listing information on the chosen system-specific parameters. Note that
parameters can be changed manually by adding a
[PARAMETER](#CP2K_INPUT.FORCE_EVAL.DFT.QS.XTB.PARAMETER) section to the
[XTB](#CP2K_INPUT.FORCE_EVAL.DFT.QS.XTB) section and specifying corresponding keywords for the
specific parameters with the adjusted values or by giving the path to the modified parameter file,
adding the keywords [PARAM_FILE_PATH](#CP2K_INPUT.FORCE_EVAL.DFT.QS.XTB.PARAMETER.PARAM_FILE_PATH)
and [PARAM_FILE_NAME](#CP2K_INPUT.FORCE_EVAL.DFT.QS.XTB.PARAMETER.PARAM_FILE_NAME).

```none
                  #####   #####        #          ####### ######   
                 #     # #     #      #              #    #     #  
                 #     # #           #    ##   ##    #    #     #  
                 #     #  #####     #      ## ##     #    ######   
                 #   # #       #   #        ###      #    #     #  
                 #    #  #     #  #        ## ##     #    #     #  
                  #### #  #####  #        ##   ##    #    ######   
                                                                   

 xTB| Parameter file                                              xTB_parameters
 xTB| Basis expansion STO-NG                                                   6
 xTB| Basis expansion STO-NG for Hydrogen                                      4
 xTB| Halogen interaction potential                                            F
 xTB| Halogen interaction potential cutoff radius                         20.000
 xTB| Nonbonded interactions                                                   F
 xTB| D3 Dispersion: Parameter                                         dftd3.dat
 xTB| Huckel constants ks kp kd                        1.850     2.250     2.000
 xTB| Huckel constants ksp k2sh                                  2.080     2.850
 xTB| Mataga-Nishimoto exponent                                            2.000
 xTB| Repulsion potential exponent                                         1.500
 xTB| Coordination number scaling kcn(s) kcn(p) kc     0.006    -0.003    -0.005
 xTB| Electronegativity scaling                                           -0.007
 xTB| Halogen potential scaling kxr kx2                          1.300     0.440
```

Analogously to any other self-consistent field optimization (SCF) method, the output also includes
the energy and convergence during the SCF steps with the finally converged GFN1-xTB energy.

```none
 SCF WAVEFUNCTION OPTIMIZATION

  ----------------------------------- OT ---------------------------------------
  Minimizer      : DIIS                : direct inversion
                                         in the iterative subspace
                                         using   7 DIIS vectors
                                         safer DIIS on
  Preconditioner : FULL_SINGLE_INVERSE : inversion of 
                                         H + eS - 2*(Sc)(c^T*H*c+const)(Sc)^T
  Precond_solver : DEFAULT
  stepsize       :    0.08000000                  energy_gap     :    0.08000000
  eps_taylor     :   0.10000E-15                  max_taylor     :             4
  ----------------------------------- OT ---------------------------------------

  Step     Update method      Time    Convergence         Total energy    Change
  ------------------------------------------------------------------------------
     1 OT DIIS     0.80E-01    0.5     0.01213502      -947.7483409153 -9.48E+02
     2 OT DIIS     0.80E-01    0.3     0.00675007      -951.5762826800 -3.83E+00
     3 OT DIIS     0.80E-01    0.3     0.00092877      -953.2164544959 -1.64E+00
     4 OT DIIS     0.80E-01    0.3     0.00034159      -953.2591478247 -4.27E-02
     5 OT DIIS     0.80E-01    0.3     0.00018348      -953.2687102329 -9.56E-03
     6 OT DIIS     0.80E-01    0.3     0.00009265      -953.2707750500 -2.06E-03
     7 OT DIIS     0.80E-01    0.3     0.00005495      -953.2714236504 -6.49E-04
     8 OT DIIS     0.80E-01    0.3     0.00002612      -953.2716704946 -2.47E-04
     9 OT DIIS     0.80E-01    0.3     0.00001585      -953.2717390500 -6.86E-05
    10 OT DIIS     0.80E-01    0.3     0.00001020      -953.2717664315 -2.74E-05
    11 OT DIIS     0.80E-01    0.3     0.00000564      -953.2717774258 -1.10E-05
    12 OT DIIS     0.80E-01    0.3     0.00000354      -953.2717818198 -4.39E-06
    13 OT DIIS     0.80E-01    0.3     0.00000206      -953.2717839406 -2.12E-06
    14 OT DIIS     0.80E-01    0.3     0.00000127      -953.2717844831 -5.42E-07
    15 OT DIIS     0.80E-01    0.3     0.00000077      -953.2717846336 -1.51E-07

  *** SCF run converged in    15 steps ***


  Core Hamiltonian energy:                                   -962.45147378153547
  Repulsive potential energy:                                   8.84897617161771
  Electronic energy:                                            0.76461561909348
  DFTB3 3rd order energy:                                       0.33228335538302
  Dispersion energy:                                           -0.76618599817727

  Total energy:                                              -953.27178463361872

  outer SCF iter =    1 RMS gradient =   0.77E-06 energy =       -953.2717846336
  outer SCF loop converged in   1 iterations or   15 steps
```

### Adding a generic correction potential

It is possible to add a generic non bonded correction potential. The potential form can be chosen
freely and needs to be specified by adding the keyword
[FUNCTION](#CP2K_INPUT.FORCE_EVAL.DFT.QS.XTB.NONBONDED.GENPOT.FUNCTION). Included parameters and
variables have to be specified using the keywords
[VARIABLES](#CP2K_INPUT.FORCE_EVAL.DFT.QS.XTB.NONBONDED.GENPOT.VARIABLES) and
[PARAMETERS](#CP2K_INPUT.FORCE_EVAL.DFT.QS.XTB.NONBONDED.GENPOT.PARAMETERS). The section can be
repeated as often as required and enables to include pairwise, element-specific correction
potentials. The implementation also features analytic gradients for this option.

```none

&GLOBAL
  RUN_TYPE  ENERGY
  PROJECT_NAME xtb
  PRINT_LEVEL  MEDIUM
  PREFERRED_DIAG_LIBRARY SL
&END GLOBAL
&FORCE_EVAL
 METHOD QS
&DFT
  &QS
   METHOD XTB
   &XTB
    CHECK_ATOMIC_CHARGES F
    DO_EWALD  T
    USE_HALOGEN_CORRECTION T
    DO_NONBONDED T               ! Possible option to include a generic non-bonded potential
     &NONBONDED                  ! Specification of the potential, keyword can be repeated
      &GENPOT
       ATOMS Kr Br
       FUNCTION Aparam*exp(-Bparam*r)-Cparam/r**8  ! Potential formula has to be specified
       PARAMETERS Aparam Bparam Cparam             ! Parameters included in the formula above
       VALUES 70.0 1.0 0.0                         ! Explicit values for the parameters
       VARIABLES r
       RCUT 40.5
      &END GENPOT
     &END NONBONDED
   &END XTB
  &END QS
  &SCF
   SCF_GUESS RESTART
   MAX_SCF 50
   EPS_SCF 1.E-6
   &OT ON
     PRECONDITIONER FULL_SINGLE_INVERSE
     MINIMIZER DIIS
   &END
   &OUTER_SCF
     MAX_SCF 200
     EPS_SCF 1.E-6
   &END OUTER_SCF
  &END SCF
 &END DFT
  &SUBSYS
    &TOPOLOGY
      COORD_FILE_FORMAT  xyz
      COORD_FILE_NAME  input.xyz
      CONNECTIVITY OFF
      &CENTER_COORDINATES
      &END CENTER_COORDINATES
     &END TOPOLOGY
    &CELL
      ABC  21.64 21.64 21.64
      ALPHA_BETA_GAMMA 90.0 90.0 90.0
      PERIODIC XYZ
    &END CELL
  &END SUBSYS
&END FORCE_EVAL
```
