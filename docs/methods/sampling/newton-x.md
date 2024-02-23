# Surface Hopping with NEWTON-X

This is a short tutorial on how to use the CP2K-NEWTONX interface to a) generate initial conditions
to compute photoabsorption spectra and b) to run non-adiabatic dynamics simulations using orbital
derivative couplings. A more comprehensive tutorial on all NEWTONX features, including a
documentation of the required specifications for the CP2K interface, can be found on the NEWTONX
homepage, <https://newtonx.org/documentation-tutorials/>.

## Brief theory recap

The interface enables to use electronic-structure data from CP2K and combine it with the surface
hopping module of NEWTONX. Excitation energies $\Omega^M$ and excited-state eigenvectors
$\mathbf{X}^M$ to describe the excited state $M$ are provided by CP2K, relying on the Tamm-Dancoff
eigenvalue problem,

$$
\mathbf{A} \mathbf{X}^M &= \Omega^M \mathbf{S} \mathbf{X}^M \, , \\
\sum_{\kappa k} [ F_{\mu \kappa \sigma} \delta_{ik} - F_{ik \sigma} S_{\mu \kappa} ] X^M_{\kappa k \sigma} + \sum_{\lambda} K_{\mu \lambda \sigma} [\mathbf{D}^{{\rm{\tiny{X}}}M}] C_{\lambda i \sigma} &=  \sum_{\kappa} \Omega^M S_{\mu \kappa} X^M_{\kappa i \sigma} \, ,
$$

with $\mathbf{S}$ representing the conventional atomic-orbital overlap matrix, $\mathbf{F}$ the
Kohn-Sham matrix, $\mathbf{K}$ the kernel comprising -- depending on the chosen functional --
Coulomb, exchange and exchange-correlation contributions, and $\mathbf{C}$ the molecular orbital
coefficients. $\mu, \nu, \dots$ denote atomic orbitals, $i, j, \dots$ occupied molecular orbitals.
The corresponding excited-state gradient is obtained setting up a variational Lagrangian and taking
the derivative with respect to the nuclear coordinates $\mathbf{R}$ (see also
[](../properties/optical/tddft)).

By performing a TDDFPT computation, excitation energies $\Omega^M (\mathbf{R}(t))$, excited-state
eigenvectors $\mathbf{X}^M (\mathbf{R}(t))$ and corresponding excited-state gradients
$\nabla \Omega^M (\mathbf{R}(t))$ are provided by CP2K. On the so-defined potential energy surfaces,
the nuclei are propagated classically relying on the surface hopping code of NEWTONX,

$$
\mathbf{R}(t + \Delta t) &= \mathbf{R} (t) + \mathbf{v} (t) \Delta t + \frac{1}{2} \mathbf{a}(t) \Delta t^2  \, ,\\
\mathbf{v} (t + \Delta t) &= \mathbf{v} (t) + \frac{1}{2} (\mathbf{a} (t) + \mathbf{a} (t+ \Delta t) ) \Delta t  \, , \\
\mathbf{a} (t) &= - \frac{1}{m} \nabla \Omega^M (\mathbf{R}(t)) \, .
$$

The coefficients $c^M (t)$ of the total wave function $\Psi (\mathbf{R}(t))$ over all excited states
$M$ are obtained implying hopping probabilities $P_{M\rightarrow N}$ of Tully's surface hopping,

$$
\Psi (\mathbf{R}(t)) &= \sum_{M} c^{M} (t) \Psi^M (\mathbf{R}(t)) \\
i \frac{{\rm{d}} c^M (t)}{{\rm{d}}t} &= \sum_N c^N (t) \left ( \delta_{MN} E_N (\mathbf{R}(t)) - i \sigma_{MN} (t) \right ) \, , \\
P_{M \rightarrow N} &= {\rm{max}} \left [ 0, \frac{-2 \Delta t}{| c^M|^2} {\rm{Re}} (c^M c^{N \ast}) \sigma_{MN} \right ] \, .
$$

The therefore required non-adiabatic time derivative couplings $\sigma_{MN}$ can be obtained relying
on semi-empirical models (Baeck-An; please cite
[Barbatti et al., Open Research Europe 1, 49 (2021)](https://doi.org/10.12688/openreseurope.13624.1).)
or as numerical time derivative couplings (orbital time derivative (OD); please cite
[Ryabinkin et al., J. Phys. Chem. Lett. 6, 4200 (2015)](https://doi.org/10.1021/acs.jpclett.5b02062);
[Barbatti et al., Molecules 21, 1603 (2016)](https://doi.org/10.3390/molecules21111603).), with the
corresponding molecular orbital overlap matrix $\mathbf{S}^{{\rm{\tiny{t-\Delta t,t}}}}$ being
provided by CP2K,

$$
\sigma_{MN}^{{\rm{\tiny{OD}}}} &= \sum_{ia} X_{ia}^{M} \frac{\partial }{\partial t} X_{ia}^N + \sum_{iab} X_{ia}^M X_{ib}^N  S_{ab}^{{\rm{\tiny{t-\Delta t,t}}}} - \sum_{ija} P_{ij} X_{ia}^M X_{ja}^N
 S_{ji}^{{\rm{\tiny{t-\Delta t,t}}}} \\
S_{pq}^{{\rm{\tiny{t - \Delta t , t}}}} &= \frac{\langle \phi_i (\mathbf{R}(t- \Delta t )) | \phi_j (\mathbf{R} (t)) \rangle}{\Delta t} \, .
$$

$a,b, \dots$ denote virtual molecular orbitals.

## General input setup

The input sections for TDDFPT energy and gradient computations are described in
[](../properties/optical/tddft). To furthermore provide the required CP2K output, subsequently read
in by NEWTONX, the following print statements have to be added to the CP2K input files:

- [FORCE_EVAL.PRINT.FORCES](#CP2K_INPUT.FORCE_EVAL.PRINT.FORCES): prints the excited-state forces
- [TDDFPT.PRINT.NAMD_PRINT](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.TDDFPT.PRINT.NAMD_PRINT) with keyword
  option [PRINT_PHASES](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.TDDFPT.PRINT.NAMD_PRINT.PRINT_PHASES):
  prints the excited-state eigenvectors in MO format as well as the corresponding phases.
- [VIBRATIONAL_ANALYSIS.PRINT.NAMD_PRINT](#CP2K_INPUT.VIBRATIONAL_ANALYSIS.PRINT.NAMD_PRINT): prints
  normal modes to generate initial conditions It should furthermore be noted that cartesian
  coordinates have to be provided in terms of the external file `coord.cp2k` and that the number of
  atoms has to be specified in the CP2K input file in the [SUBSYS](#CP2K_INPUT.FORCE_EVAL.SUBSYS)
  section.

## A) Initial conditions and photoabsorption spectra

The following tutorial to obtain photoabsorption spectra is based on section 2 of
<https://vdv.dcf.mybluehost.me/nx/wp-content/uploads/2020/02/tutorial-2_2.pdf>. For the
electronic-structure calculation with CP2K, a `cp2k.inp` and `cp2k.par` file as well as a coordinate
file named `coord.cp2k` has to be provided in a subdirectory called `JOB_AD`. Furthermore, a
vibrational analysis computation has to be performed to provide cartesian normal modes, with the
input file including the corresponding `NAMD print` section.

Examplary input files for computing the absorption spectrum as well as for performing a vibrational
analysis for a single water molecule with CP2K are given below:

```none
&GLOBAL
  PROJECT excited_states_for_h2o 
  RUN_TYPE ENERGY
  PREFERRED_DIAG_LIBRARY SL
  PRINT_LEVEL medium
&END GLOBAL
&FORCE_EVAL
  &PRINT                      # print statement for ground-state or excited-state forces
    &FORCES
    &END FORCES
  &END PRINT
  METHOD Quickstep
  &PROPERTIES
    &TDDFPT                    # TDDFPT input section to compute 10 excited states
      &DIPOLE_MOMENTS
        DIPOLE_FORM LENGTH
      &END DIPOLE_MOMENTS
      KERNEL FULL
      NSTATES 10
      MAX_ITER   100
      MAX_KV 20
      CONVERGENCE [eV] 1.0e-5
      RKS_TRIPLETS F
      &PRINT                     # NAMD print section to print excited-state eigenvectors
        &NAMD_PRINT
          PRINT_VIRTUALS T
          PRINT_PHASES T
        &END NAMD_PRINT
      &END PRINT
    &END TDDFPT
  &END PROPERTIES
  &DFT
    &QS
      METHOD GAPW
      EPS_DEFAULT 1.0E-17
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
    POTENTIAL_FILE_NAME POTENTIAL
    BASIS_SET_FILE_NAME EMSL_BASIS_SETS
    &MGRID
      CUTOFF 1000
      REL_CUTOFF 100
      NGRIDS 5
    &END MGRID
    &POISSON
      PERIODIC NONE
      PSOLVER MT
    &END
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC 8.0 8.0 8.0
      PERIODIC NONE
    &END CELL
                                    # Coordinates are provided externally for the interface
    &COORD
      @include coord.cp2k
    &END COORD
    &TOPOLOGY
      &CENTER_COORDINATES T
      &END
      NATOMS 3                       # specifying number of atoms for NEWTONX
      CONNECTIVITY OFF
    &END TOPOLOGY
    &KIND H
      BASIS_SET 6-311Gxx
      POTENTIAL ALL
    &END KIND
    &KIND O
      BASIS_SET 6-311Gxx
      POTENTIAL ALL
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
```

```none
&GLOBAL
  PROJECT normal_modes_for_h2o
  RUN_TYPE VIBRATIONAL_ANALYSIS      #computing normal modes to generate initial conditions
  PREFERRED_DIAG_LIBRARY SL
  PRINT_LEVEL medium
&END GLOBAL
&FORCE_EVAL
  &PRINT
    &FORCES
    &END FORCES
  &END PRINT
  METHOD Quickstep
  &DFT
    &QS
      METHOD GAPW                   # GAPW enables comparison with all-electron molecular program codes like Turbomole
      EPS_DEFAULT 1.0E-17
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
    POTENTIAL_FILE_NAME POTENTIAL
    BASIS_SET_FILE_NAME EMSL_BASIS_SETS
    &MGRID
      CUTOFF 1000
      REL_CUTOFF 100
      NGRIDS 5
    &END MGRID
    &POISSON
      PERIODIC NONE
      PSOLVER MT
    &END
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC 8.0 8.0 8.0
      PERIODIC NONE
    &END CELL
                                    # coordinates must be provided as external file for NEWTONX
    &COORD
      @include coord.cp2k
    &END COORD
    &TOPOLOGY
      &CENTER_COORDINATES T
      &END
      NATOMS 3
      CONNECTIVITY OFF
    &END TOPOLOGY
    &KIND H
      BASIS_SET 6-311Gxx
      POTENTIAL ALL
    &END KIND
    &KIND O
      BASIS_SET 6-311Gxx
      POTENTIAL ALL
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
&VIBRATIONAL_ANALYSIS
  &PRINT
    &NAMD_PRINT                      # keyword to enable printing of cartesian normal modes
    &END NAMD_PRINT
  &END PRINT
  DX 0.001
&END VIBRATIONAL_ANALYSIS
```

The input file `cp2k.par` includes all specifications regarding the executable and parallelization
setup.

```none
 parallel = 16
 exec = cp2k.psmp
```

Furthermore, a `initqp_input` file has to be generated for NEWTONX following the instructions given
in the NEWTONX tutorial. Specifications for CP2K in the `initqp_input` file are the following:

- The file comprising the normal modes of the CP2K frequency computation -- for the above input
  provided as `normal_modes_for_h2o-VIBRATIONS-1.eig`-- has to be specified as
  `file_nmodes = normal_modes_for_h2o-VIBRATIONS-1.eig`.
- The electronic structure program has to be specified as CP2K by defining `iprog = 10`.

```none
&dat
 nact = 2
 iprog = 10
 numat = 3
 npoints = 500
 file_geom = geom
 file_nmodes = normal_modes_for_h2o-VIBRATIONS-1.eig
 anh_f = 1
 rescale = n
 temp = 0
 ics_flg = n
 chk_e = 1
 nis = 1
 nfs = 11
 kvert = 1
 de = 100
 prog = 14
 iseed = 0
 lvprt = 1
/
```

After providing the excited-state CP2K computation based on input file `h2o_cp2k.inp` in the
subdirectory `JOB_AD`, the normal modes `normal_modes_for_h2o-VIBRATIONS-1.eig` of the frequency
computation and the `initqp_input` file for NEWTONX, the script initcond.pl of NEWTONX can be
executed to generate initial conditions. The resulting initcond-output file of NEWTONX, it is first
stated that the read-in cartesian normal modes are transferred to mass-weighted normal modes.

```none
Cartesian normal modes (1/sqrt(amu))

        0.00        0.00        0.00        0.00        0.00        0.00     1523.92     3851.12

      0.0000     -0.0492      0.0001     -0.1268      0.5632     -0.0083      0.0000     -0.0000
     -0.0886      0.0000     -0.0000     -0.0169      0.0047      0.5777      0.0000     -0.0000
     -0.0000     -0.0000     -0.0000      0.5630      0.1269      0.0155     -0.0715      0.0487
      0.0001      0.3905     -0.0004     -0.1267      0.5632     -0.0082     -0.4184     -0.5910
      0.7043      0.0008      0.7071     -0.0162      0.0040      0.5768      0.0000      0.0000
     -0.0001     -0.5885      0.0007      0.5630      0.1270      0.0155      0.5678     -0.3867
      0.0000      0.3905     -0.0004     -0.1267      0.5632     -0.0083      0.4184      0.5910
      0.7043     -0.0009     -0.7071     -0.0170      0.0051      0.5768      0.0000      0.0000
     -0.0000      0.5885     -0.0007      0.5630      0.1269      0.0154      0.5678     -0.3867

     3986.44

      0.0712
     -0.0000
      0.0000
     -0.5650
      0.0000
     -0.4222
     -0.5650
      0.0000
      0.4222

Mass weighted normal modes
Frequencies will be multiplied by ANH_F =    1.00000

        0.00        0.00        0.00        0.00        0.00        0.00     1523.92     3851.12

      0.0001     -0.1967      0.0006     -0.5069      2.2526     -0.0330      0.0000     -0.0000
     -0.3543      0.0000     -0.0000     -0.0677      0.0186      2.3104      0.0000     -0.0000
     -0.0001     -0.0000     -0.0002      2.2517      0.5077      0.0619     -0.2861      0.1949
      0.0001      0.3920     -0.0004     -0.1272      0.5654     -0.0083     -0.4200     -0.5933
      0.7071      0.0008      0.7099     -0.0162      0.0040      0.5791      0.0000      0.0000
     -0.0001     -0.5908      0.0007      0.5652      0.1275      0.0155      0.5700     -0.3882
      0.0000      0.3921     -0.0004     -0.1272      0.5654     -0.0083      0.4200      0.5933
      0.7071     -0.0009     -0.7099     -0.0171      0.0051      0.5790      0.0000      0.0000
     -0.0000      0.5908     -0.0007      0.5652      0.1274      0.0155      0.5700     -0.3882

     3986.44

      0.2847
     -0.0000
      0.0000
     -0.5672
      0.0000
     -0.4238
     -0.5672
      0.0000
      0.4238
```

The thereon based initial conditions are summarized in external output files for each state, dubbed
"final_output_XXX", comprising information on the various geometries and velocities as examplarily
given below:

```none
 Initial condition =     1
 Geometry in COLUMBUS and NX input format:
 o     8.0    5.00630777    5.00000001    4.46399957   15.99491464
 h     1.0    6.37684065    5.00000128    5.50815661    1.00782504
 h     1.0    3.52303474    5.00000149    5.58297278    1.00782504
 Velocity in NX input format:
   -0.000089112    0.000000000   -0.000020915
    0.000417197    0.000000002    0.000694479
    0.000997296    0.000000013   -0.000362483
 Epot of initial state (eV):    0.0865  Epot of final state (eV):     19.0799
 Vertical excitation (eV):     18.9935  Is Ev in the required range? YES
 Ekin of initial state (eV):    0.0479  Etot of initial state (eV):    0.1343
 Oscillator strength:           0.1221
 State:                         10
```

Based on the initial conditions, the broadened photoabsorption spectrum can be computed with the
nxinp script. As outlined in section 2.7 of the cited NEWTONX tutorial, the so-obtained output file
`cross-section.dat` comprises the data points of the computed photoabsorption spectrum as visualized
below:

## B) Non-adiabatic dynamics using orbital determinant derivatives
