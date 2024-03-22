# Bethe-Salpeter Equation

In this section, we discuss the basics for computing linear optical properties of molecules within
the Bethe-Salpeter formalism. The Bethe-Salpeter equation (BSE) enables the computation of
electronic excitations with and without the Tamm-Dancoff approximation (TDA) based on $G_0W_0$
eigenvalues (cf. [GW]-subsection for details on the GW method), i.e. we perform BSE@$G_0W_0$@DFT.
The corresponding input options can be found in the [BSE]-section.

## Theory

The current implementation in CP2K includes the full diagonalization of the full BSE in the static
approximation, which reads in matrix notation

$$\left( \begin{array}{cc}A &  B\\B &  A\end{array} \right)\left( \begin{array}{cc}\mathbf{X}^{(n)}\\\mathbf{Y}^{(n)}\end{array} \right) = \Omega^{(n)}\left(\begin{array}{cc}1&0\\0&-1\end{array}\right)\left(\begin{array}{cc}\mathbf{X}^{(n)}\\\mathbf{Y}^{(n)}\end{array}\right) \quad ,$$

as well as the solution of the reduced problem in the TDA, i.e. $B = 0$

$$ A \mathbf{X}^{(n)}_\mathbf{TDA} = \Omega^{(n)}_\mathbf{TDA} \mathbf{X}^{(n)}_\mathbf{TDA} \quad .$$

In this notation, the $n$-th excited state is characterized by the excitation energy $\Omega^{(n)}$
and the eigenvector $\mathbf{X}^{(n)}$. The excitation character, i.e. transitions from occupied to
unoccupied molecular orbitals (MO's), introduces a combined index structure, e.g. the component
$X_{ia}^{(n)}$ of the eigenvector describes the transition amplitude from occupied orbital $i$ to
unoccupied orbital $a$ contributing to the $n$-th excitation.

The matrices $A$ and $B$ in this product space of occupied ($i,j$) and unoccupied ($a,b$) orbitals
read

$$ \begin{align}
    A_{ia,jb} &= (\varepsilon_a-\varepsilon_i)\delta_{ij}\delta_{ab} + \alpha^\mathrm{(S/T)}
    v_{ia,jb} - W_{ij,ab}(\omega=0) \quad ,\\
    B_{ia,jb} &= \alpha^\mathrm{(S/T)} v_{ia,bj} - W_{ib,aj}(\omega=0) \quad .
\end{align}$$

The scalar factor $\alpha^\mathrm{(S/T)}$ describes the spin configuration of the computed
excitation (cf.
[SPIN_CONFIG](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.SPIN_CONFIG)). Further, the
Kronecker delta $\delta_{p,q}$ as well as several quantities from GW enter (cf. [](#Wilhelm2016)):
the GW quasiparticle energies $\varepsilon_{p}$, the bare Coulomb interaction $v_{pq,rs}$ and the
statically screened Coulomb interaction $W_{pq,rs}(\omega=0)$, where $(p,q,r,s)$ denote the indices
from either subspace.

Until now, there is no iterative diagonalization scheme implemented, i.e. the runtime scales with
$\mathcal{O}(N^6)$ and memory consumption scales with $\mathcal{O}(N^4)$, where $N$ is the number of
electrons.

A manuscript with detailed definitions and formulas is currently in preparation.

## BSE input section

```{note}
The accuracy of the BSE relies heavily on well-converged settings in the prior DFT and GW steps (BSE@$G_0W_0$@DFT). 
For example, the chosen [XC_FUNCTIONAL](#CP2K_INPUT.FORCE_EVAL.DFT.XC.XC_FUNCTIONAL.SECTION_PARAMETERS) and the parameters for the analytic continutation in GW ([QUADRATURE_POINTS](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.QUADRATURE_POINTS), [NPARAM_PADE](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.NPARAM_PADE) and [OMEGA_MAX_FIT](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.OMEGA_MAX_FIT)) can have a profound influence on the excitation energies.
In particular, all MO's included in the BSE have to be corrected in GW by setting [CORR_MOS_OCC](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.CORR_MOS_OCC) and [CORR_MOS_VIRT](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.CORR_MOS_VIRT) to a sufficiently large number. By default, one should set it to a negative number, which defaults to the correction of all orbitals, as specified in the manual sections given above.
```

```{note}
The [BSE] is implemented on top of the old [GW] code. Please make sure that you are editing the correct subsections (cf. [GW]-subsection) of the input.
```

The parameters defining the BSE calculation have to be specified in the [BSE]-subsection when
running a [RUN_TYPE](#CP2K_INPUT.GLOBAL.RUN_TYPE) `ENERGY` calculation. As highlighted above, ensure
a converged DFT and [GW] calculation before. The most important keywords are:

- [BSE_APPROX](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.BSE_APPROX): Option to
  switch on/off the TDA, i.e. $B=0$. Can either be the full BSE, TDA or both.
- [ENERGY_CUTOFF_OCC](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.ENERGY_CUTOFF_OCC):
  Cutoff for the occupied MO's. If the energy of MO with index $i$ is more than `ENERGY_CUTOFF_OCC`
  away from the HOMO (highest occupied molecular orbital) energy, i.e.
  $\epsilon_{\mathrm{HOMO}}-\epsilon_i>\mathtt{ENERGY\_CUTOFF\_OCC}$, then all occupied MO's with
  index $\le i$ are not included in the matrices $A_{ia,jb}$ and $B_{ia,jb}$.
- [ENERGY_CUTOFF_VIRT](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.ENERGY_CUTOFF_VIRT):
  Cutoff for the unoccupied MO's. If the energy of MO with index $a$ is more than
  `ENERGY_CUTOFF_VIRT` away from the LUMO (lowest unoccupied molecular orbital) energy, i.e.
  $\epsilon_a-\epsilon_{\mathrm{LUMO}}>\mathtt{ENERGY\_CUTOFF\_VIRT}$, then all unoccupied MO's with
  index $\ge a$ are not included in the matrices $A_{ia,jb}$ and $B_{ia,jb}$.

```{note}
Usage of [ENERGY_CUTOFF_OCC](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.ENERGY_CUTOFF_OCC) and [ENERGY_CUTOFF_VIRT](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.ENERGY_CUTOFF_VIRT) requires careful checks of convergence!
```

- [EPS_X](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.EPS_X): Threshold for printing
  entries $X_{ia}^{(n)}$ of the eigenvector. Only entries with absolute value larger than `EPS_X`
  are printed.
- [NUM_PRINT_EXC](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.NUM_PRINT_EXC): Number
  of excitation energies $\Omega^{(n)}$ and eigenvectors $\mathbf{X}^{(n)}$ to be printed. Does not
  affect computation time.
- [SPIN_CONFIG](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.SPIN_CONFIG): Specifies
  the desired spin configuration of the excitation. This determines the scalar factor
  $\alpha^\mathrm{(S/T)}$, which is 2 for Singlet excitations and 0 for Triplet excitations.

## Minimal example for a BSE calculation

## Input

The following input for cp2k is a minimal example, which can be run quickly by

```none
mpirun -n 1 cp2k.psmp BSE_H2.inp
```

within 5 minutes on 1 core (yet with a considerable memory consumption). Input and output files are
also available [here](https://www.cp2k.org/_media/howto:bse_example_h2.zip).

```none
&GLOBAL
  PROJECT  H2
  RUN_TYPE ENERGY
&END GLOBAL
&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME BASIS-aug       ! Custom Basis set file (aug-cc-pVDZ and aug-cc-pVDZ-RIFIT from Basis Set exchange)
    POTENTIAL_FILE_NAME POTENTIAL     
    &QS
      METHOD GAPW                       ! All electron calculation
      EPS_DEFAULT 1.0E-40
      EPS_PGF_ORB 1.0E-40
    &END QS
    &POISSON
      PERIODIC NONE
      PSOLVER MT
    &END
    &SCF
      SCF_GUESS RESTART
      EPS_SCF 1e-7
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE                ! Choice of functional has a profound influence on the results
      &END XC_FUNCTIONAL
      &WF_CORRELATION
        &RI_RPA
          QUADRATURE_POINTS 500         ! Convergence parameter: Check carefully
          &GW
            CORR_OCC   -1               ! Correct all occupied orbitals
            CORR_VIRT  -1               ! Correct all unoccupied orbitals
            RI_SIGMA_X
            NPARAM_PADE 128             ! Convergence parameter: Check carefully
            OMEGA_MAX_FIT 3.675         ! Convergence parameter: Check carefully
            &BSE
              BSE_APPROX BOTH           ! In this case, full BSE and TDA are calculated
            &END BSE
          &END GW
        &END RI_RPA
      &END
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC 25 25 25
      PERIODIC NONE
    &END CELL
    &COORD
      H 0.0000 0.0000 0.0000            ! H2 molecule geometry from GW100 Paper
      H 0.0000 0.0000 0.74144
    &END COORD
    &TOPOLOGY
      &CENTER_COORDINATES
      &END
    &END TOPOLOGY
    &KIND H
      BASIS_SET ORB    aug-cc-pVDZ      ! Basis sets also require convergence checks
      BASIS_SET RI_AUX aug-cc-pVDZ-RIFIT
      POTENTIAL ALL
    &END KIND
    &KIND O
      BASIS_SET ORB    aug-cc-pVDZ
      BASIS_SET RI_AUX aug-cc-pVDZ-RIFIT
      POTENTIAL ALL
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
```

The basis sets `aug-cc-pVDZ` and `aug-cc-pVDZ-RIFIT` in `BASIS-aug` can be obtained from the Basis
Set Exchange Library:
<a href="https://www.basissetexchange.org/basis/aug-cc-pvdz/format/cp2k/?version=1&elements=1" target="_blank">aug-cc-pVDZ</a>,
<a href="https://www.basissetexchange.org/basis/aug-cc-pvdz-rifit/format/cp2k/?version=1&elements=1" target="_blank">aug-cc-pVDZ-RIFIT</a>.
The geometry for $H_2$ was taken from the [GW100](https://doi.org/10.1021/acs.jctc.5b00453)-Paper.

## Output

In the resulting output file (cf. [Download](https://www.cp2k.org/_media/howto:bse_example_h2.zip)),
the [BSE]-section itself starts with a banner after the [GW]-section. Therein, all lines a formatted
with a trailing `BSE|`.

### Characteristics of the BSE run

We examine the first lines of the BSE output, which contain a number of characteristics for the run,
line by line:

```none
 BSE| No cutoff given for occupied orbitals
 BSE| No cutoff given for virtual orbitals
 BSE| First occupied index                                                     1
 BSE| Last virtual index (not MO index!)                                      17
 BSE| Energy of first occupied index [eV]                                -15.501
 BSE| Energy of last virtual index [eV]                                  101.702
 BSE| Energy difference of first occupied index to HOMO [eV]              -0.000
 BSE| Energy difference of last virtual index to LUMO [eV]                99.580
 BSE| Number of GW-corrected occupied MOs                                      1
 BSE| Number of GW-corrected virtual MOs                                      17
 BSE|
 BSE| No truncation of BSE matrices applied
 BSE|
 BSE| Total peak memory estimate from BSE [GB]                             0.000
 BSE| Peak memory estimate per MPI rank from BSE [GB]                      0.000
 BSE|
 BSE| Diagonalizing aux. matrix with size of A. This will take around 0.E+00 s.
 BSE| Diagonalizing A. This will take around 0.E+00 s.
```

From the first two lines, we see that there was no cutoff applied, i.e.
[ENERGY_CUTOFF_OCC](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.ENERGY_CUTOFF_OCC)
and
[ENERGY_CUTOFF_VIRT](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.ENERGY_CUTOFF_VIRT)
have not been invoked. Therefore, all indices from GW enter, which is one occupied orbital and 17
virtual/unoccupied orbitals. As already denoted in the output, the last virtual index is not the MO
index, but the number within the unoccupied orbitals. The following four lines contain information
about the energy structure of the included orbitals. </br> In line 9 and 10, the actual number of
GW-corrected orbitals are listed (cf.
[CORR_MOS_OCC](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.CORR_MOS_OCC) and
[CORR_MOS_VIRT](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.CORR_MOS_VIRT)). In case, the
buildup of the BSE-matrices would include MOs with uncorrected energies, the run will abort.

The characteristics section is then finalized by a message about the truncation (without cutoffs,
there are no cutoffs applied) and estimates for the peak memory
($\mathcal{O}((N_\mathrm{occupied} N_\mathrm{unoccupied})^2)$ $=\mathcal{O}(N^4)$) and runtime per
diagonalization ($\mathcal{O}((N_\mathrm{occupied} N_\mathrm{unoccupied})^3)$ $=\mathcal{O}(N^6)$).

### Results

Depending on the chosen
[BSE_APPROX](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.BSE_APPROX), a banner
signalizes the start of the respective results section. Afterwards, the most important formulas and
quantities are summarized before the excitation energies and the single-particle transitions, i.e.
the eigenvector elements $X_{ia}^{(n)}$, are printed up to the given
[EPS_X](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.EPS_X). For the full solution (no
TDA applied), the first lines of the output for the energies of the requested singlet excitation
look like

```none
 BSE| Excitation energies from solving the BSE without the TDA:
 BSE|
 BSE|       Excitation        Multiplet  TDA/full BSE     Excitation energy (eV)
 BSE|                1    Singlet State        -full-                    11.4625
 BSE|                2    Singlet State        -full-                    12.4855
 BSE|                3    Singlet State        -full-                    15.2853
```

with the columns excitation index $n$, the requested multiplet, the
[BSE_APPROX](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.BSE_APPROX) and the energy
in eV. The single-particle transitions are displayed like that (for the first three excitations):

```none
 BSE| Excitations are built up by the following single-particle transitions,
 BSE| neglecting contributions where |X_ia^n| <  0.10 :
 BSE|         -- Quick reminder: HOMO i =    1 and LUMO a =    2 --
 BSE|
 BSE| n-th exc.    i =>     a            TDA/full BSE                   |X_ia^n|
 BSE|
 BSE|         1    1 =>     2                  -full-                     0.6681
 BSE|         1    1 =>     4                  -full-                     0.2465
 BSE|
 BSE|         2    1 =>     3                  -full-                     0.7060
 BSE|
 BSE|         3    1 =>     5                  -full-                     0.7078
```

Here, some reminders for [EPS_X](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.EPS_X)
and the indices of HOMO and LUMO are printed. The columns display the excitation index $n$, the
single-particle indices $i$ and $a$ contributing to the $n$-th excitation, the
[BSE_APPROX](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.BSE_APPROX) and the absolute
value of the eigenvector entry $|X_{ia}^{(n)}|$.

In the case of the $H_2$, the first excitation is mainly built up by a transition from the first MO
to the second MO, i.e. the LUMO, but also contains a considerable contribution from the 1=>4
(HOMO=>LUMO+2) transition. The remaining contributions of the normalized $\mathbf{X}^{(n)}$ are
smaller than `0.10` and are therefore not printed.

[bse]: #CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE
[gw]: #CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW
