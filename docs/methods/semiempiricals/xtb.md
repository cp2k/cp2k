# eXtended Tight Binding

This is a short tutorial on how to run GFN1-xTB computations. The details on the theory and the
original implementation by Grimme can be found in [](#Grimme2017). Please cite this paper if you
were to use the GFN1-xTB module.

## Brief theory recap

The semi-empirical GFN1-xTB energy expression comprises contributions due to electronic (EL),
atom-pairwise repulsion (REP), dispersion (DISP), and halogen-bonding (XB) terms,

$$
E_{\rm{\tiny{GFN1-xTB}}} = E_{\rm{\tiny{EL}}} + E_{\rm{\tiny{REP}}} + E_{\rm{\tiny{DISP}}} + E_{\rm{\tiny{XB}}} \, .
$$

### 1. Electronic energy

The electronic energy contribution,

$$
E_{\rm{\tiny{EL}}} =  \sum_i^{\rm{\tiny{occ}}} n_i \langle \Psi_i | h_0 | \Psi_i \rangle + \frac{1}{2} \sum_{A,B} \sum_{{l}^A}\sum_{{l'}^B} p_l^A p_{{l'}}^B \gamma_{AB,ll'} + \frac{1}{3}\sum_{A} \Gamma_A q_A^3 - T_{\rm{\tiny{el}}} S_{\rm{\tiny{el}}} \, ,
$$

contains zeroth-order contributions based on a zeroth-order Hamiltonian $h_0$, the valence molecular
orbitals $\Psi_i$, occupation numbers $n_i$ as well as second-order contributions which are
optimized self-consistently as well as third-order diagonal contributions. The second order
contributions are described using the semi-empirical electron repulsion operator $\gamma_{AB,ll'}$
which depends on the interatomic distance of atoms $A$ and $B$ as well as further empirical
parameters that are specific for different angular momenta $l$ and $l'$. The monopole charges of the
second-order expression are optimized self-consistently,

$$
p_l^A = p_l^{A_0} - \sum_{\nu}^{N_{\rm{\tiny{AO}}}} \sum_{\mu \in A, \mu \in l} S_{\mu \nu } P_{\mu \nu} \, ,
$$

referring to the atomic orbital overlap matrix $\mathbf{S}$ and the density matrix $\mathbf{P}$.

The remaining diagonal terms represent a cubic charge correction based on the Mulliken charge $q_A$
of atom $A$ and the charge derivative $\Gamma_A$ of the atomic Hubbard parameter $\eta_A$.
Furthermore, the electronic temperature times entropy term $T_{\rm{\tiny{el}}}S_{\rm{\tiny{el}}}$
enables fractional orbital occupations.

### 2. Repulsion

Repulsion is described via an atom-pairwise potential,

$$
E_{\rm{\tiny{REP}}} = \sum_{AB} \frac{Z_A^{\rm{\tiny{eff}}} Z_B^{\rm{\tiny{eff}}} }{R_{AB}} \exp^{- (\alpha_A \alpha_B)^{1/2} (R_{AB})^{k_f}} \, ,
$$

with the effective nuclear charge $\mathbf{Z}^{\rm{\tiny{eff}}}$ as well as the global or
element-specific parameters $k_f$ and $\alpha$.

### 3. Dispersion

Dispersion is included by the well-established D3 method in the BJ-damping scheme [](#Grimme2010).

### 4. Corrections

Corrections for element-specific interactions are possible using either a halogen-bonding correction
term (XB) or a generic nonbonding potential correction (NONBOND). Note that the generic nonbonding
potential correction is CP2K specific and thus the so-obtained energy differs from the original
GFN1-xTB method,

$$
E_{\rm{\tiny{GFN1-xTB+NONBOND}}} = E_{\rm{\tiny{GFN1-xTB}}} + E_{\rm{\tiny{NONBOND}}} \, .
$$

## The GFN1-xTB input section

The most important keywords and subsections of section [XTB](#CP2K_INPUT.FORCE_EVAL.DFT.QS.XTB) are:

- [DO_EWALD](#CP2K_INPUT.FORCE_EVAL.DFT.QS.XTB.DO_EWALD): keyword to activate Ewald summation for
  periodic boundary conditions (PBC); has to be switched to true in case of PBC
- [USE_HALOGEN_CORRECTION](#CP2K_INPUT.FORCE_EVAL.DFT.QS.XTB.USE_HALOGEN_CORRECTION): keyword to
  switch off contribution $E_{\rm{\tiny{XB}}}$ to correct halogen interactions, default is to
  include this correction
- [CHECK_ATOMIC_CHARGES](#CP2K_INPUT.FORCE_EVAL.DFT.QS.XTB.CHECK_ATOMIC_CHARGES): the cubic charge
  diagonal contribution is checked to be numerically stable by switching the keyword to true.
- [DO_NONBONDED](#CP2K_INPUT.FORCE_EVAL.DFT.QS.XTB.DO_NONBONDED): add a generic correction potential
  to correct bond- or atomic-specific interactions
- [PARAMETER](#CP2K_INPUT.FORCE_EVAL.DFT.QS.XTB.PARAMETER): it is possible to add this section with
  corresponding keywords to modify xTB parameters

The additional keywords
[COULOMB_INTERACTION](#CP2K_INPUT.FORCE_EVAL.DFT.QS.XTB.COULOMB_INTERACTION),
[COULOMB_LR](#CP2K_INPUT.FORCE_EVAL.DFT.QS.XTB.COULOMB_LR) and
[TB3_INTERACTION](#CP2K_INPUT.FORCE_EVAL.DFT.QS.XTB.TB3_INTERACTION) are for debugging purposes only
and it is recommended to use the default options here.

## Simple examples

### GFN1-xTB ground-state energy for

The following input is an examplary standard input for calculating GFN1-xTB ground-state energies.

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
