# *GW* + Bethe-Salpeter equation

In this section, we discuss the basics for computing optical properties of molecules using the
Bethe-Salpeter equation (BSE) in CP2K \[[Graml2024b](#Graml2024b)\]. The BSE enables the computation
of electronic excitation energies and optical absorption spectra, for a review, see
\[[Blase2018](#Blase2018), [Blase2020](#Blase2020), [Bruneval2015](#Bruneval2015),
[Sander2015](#Sander2015)\]. In this howto, we describe in Sec.
[1](#theory-and-implementation-of-bse) the theory and implementation of BSE, in Sec. [2](#bse-input)
the BSE input keywords and in Sec. [3](#minimal-example-for-a-bse-calculation) a full CP2K input
file of a BSE calculation and the corresponding output.

## 1. Theory and implementation of BSE

A central goal of a BSE calculation is to compute electronic excitation energies
$\Omega^{(n)}, n=1,2,\ldots$ (cf. Refs. \[[Blase2018](#Blase2018), [Blase2020](#Blase2020),
[Graml2024b](#Graml2024b)\] for more usecases and details).

The following ingredients are necessary for computing $\Omega^{(n)}$:

- Occupied Kohn-Sham (KS) orbitals $\varphi_i(\mathbf{r})$ and empty KS orbitals
  $\varphi_a(\mathbf{r})$ from a DFT calculation, where $i=1,\ldots,N_\mathrm{occ}$ and
  $a=N_\mathrm{occ}+1,\ldots,N_\mathrm{occ}+N_\mathrm{empty}$,
- $GW$ eigenvalues $\varepsilon_i^{GW}$ and $\varepsilon_a^{GW}$ of corresponding KS orbitals.

In CP2K, it is possible to use $G_0W_0$, ev$GW_0$ or ev$GW$ eigenvalues, see details in [GW] and in
Ref. \[[Golze2019](#Golze2019)\], i.e. we perform BSE@$G_0W_0$/ev$GW_0$/ev$GW$@DFT. Thus, also input
parameters for a DFT and $GW$ calculation can be given (see full input in Sec. [3.1](#input-file)).
We obtain optical properties from BSE solving the following generalized eigenvalue problem that
involves the block matrix $ABBA$:

$$\left( \begin{array}{cc}A &  B\\B &  A\end{array} \right)\left( \begin{array}{cc}\mathbf{X}^{(n)}\\\mathbf{Y}^{(n)}\end{array} \right) = \Omega^{(n)}\left(\begin{array}{cc}1&0\\0&-1\end{array}\right)\left(\begin{array}{cc}\mathbf{X}^{(n)}\\\mathbf{Y}^{(n)}\end{array}\right) \quad .$$

We abbreviate $A$ and $B$ as matrices with index $A_{ia,jb}$, i.e. they have
$N_\mathrm{occ}N_\mathrm{empty}$ rows and $N_\mathrm{occ}N_\mathrm{empty}$ columns. The entries of
$A$ and $B$ are given by \[[Blase2018](#Blase2018)\]

$$ \begin{align}
    A_{ia,jb} &= (\varepsilon_a^{GW}-\varepsilon_i^{GW})\delta_{ij}\delta_{ab} + \alpha^\mathrm{S/T}
    v_{ia,jb} - W_{ij,ab}(\omega=0) \quad ,\\
    B_{ia,jb} &= \alpha^\mathrm{(S/T)} v_{ia,bj} - W_{ib,aj}(\omega=0) \quad .
\end{align}$$

where $\delta_{ij}$ is the Kronecker delta. The user sets $\alpha^S=2$ for computing singlet
excitations and $\alpha^T=0$ for computing triplet excitations. $v_{pq,rs}$ is the bare Coulomb
interaction and $W_{pq,rs}(\omega=0)$ the statically ($\omega=0$) screened Coulomb interaction ,
where $p,q,r,s \in [ 1, N_\mathrm{occ}+N_\mathrm{empty}]$ are KS orbital indices.
$(\mathbf{X}^{(n)},\mathbf{Y}^{(n)})$ with elements $X_{ia}^{(n)}$ and $Y_{ia}^{(n)}$ are the
eigenvectors of the excitation $n$ which relate to the wavefunction of the electronic excitation
\[[Blase2020](#Blase2020)\],

$$ \begin{align}
\Psi_\text{excitation}^{(n)}(\mathbf{r}_e,\mathbf{r}_h) = \sum_{ia} X_{ia}^{(n)} \varphi_i(\mathbf{r}_h) \varphi_a(\mathbf{r}_e) + Y_{ia}^{(n)} \varphi_i(\mathbf{r}_e) \varphi_a(\mathbf{r}_h) \quad ,
\end{align}$$

i.e. $X_{ia}^{(n)}$ and $Y_{ia}^{(n)}$ describe the transition amplitude between occupied orbital
$\varphi_i$ and empty orbital $\varphi_a$ of the $n$-th excitation.

the Tamm-Dancoff approximation (TDA) is also implemented, which constrains $B=0$. In case $A$ is
positive definite, and excitation energies $\Omega^{(n)}>0$, we have $\mathbf{Y}=0$ and $\mathbf{X}$
can be computed from the Hermitian eigenvalue problem

$$ A \mathbf{X}^{(n)}_\mathbf{TDA} = \Omega^{(n)}_\mathbf{TDA} \mathbf{X}^{(n)}_\mathbf{TDA} \quad .$$

Diagonalizing $A$ in TDA, or the full block-matrix $ABBA$, takes in the order of
$(N_\mathrm{occ} N_\mathrm{empty})^3$ floating point operations. This translates to a computational
scaling of $O(N^6)$ in the system size $N$, e.g. the number of electrons.

## 2. BSE input

For starting a BSE calculation one needs to set the [RUN_TYPE](#CP2K_INPUT.GLOBAL.RUN_TYPE) to
`ENERGY` and the following sections for [GW] and [BSE]:

```
&GW
  SELF_CONSISTENCY      G0W0      ! can be changed to EV_GW0 or EV_GW
  &BSE  
    TDA                 TDA+ABBA  ! Diagonalizing ABBA and A
    SPIN_CONFIG         SINGLET   ! or TRIPLET
    NUM_PRINT_EXC       20        ! Number of printed excitations
    ENERGY_CUTOFF_OCC   -1        ! Set to positive numbers (eV) to
    ENERGY_CUTOFF_EMPTY -1        ! truncate matrices A_ia,jb and B_ia,jb
  &END BSE
&END GW
```

In the upper GW/BSE section, the following keywords have been used:

- [SELF_CONSISTENCY](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.SELF_CONSISTENCY):
  Determines which GW self-consistency ($G_0W_0$, ev$GW_0$ or ev$G_0W_0$) is used to calculate the
  single-particle GW energies $\varepsilon_p^{GW}$ needed in the BSE calculation.

- [TDA](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.TDA): Three options available:

  - `ON` diagonalize of $A$,
  - `OFF` Generalized diagonalization of $ABBA$,
  - `TDA+ABBA` CP2K diagonalizes $ABBA$ as well as $A$.

- [SPIN_CONFIG](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.SPIN_CONFIG): Two options
  available:

  - `SINGLET` for computing singlet excitation energies ($\alpha^S=2$),
  - `TRIPLET` for computing triplet excitation energies ($\alpha^T=0$).

- [NUM_PRINT_EXC](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.NUM_PRINT_EXC): Number
  of excitations $N_\text{exc}^\text{print}$ to be printed.

- [ENERGY_CUTOFF_OCC](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.ENERGY_CUTOFF_OCC)
  $E_\text{cut}^\text{occ}$: Restrict occupied molecular orbital (MO) indices $i$ and only use
  occupied MOs with
  $\varepsilon_i\in[\varepsilon_{i=\text{HOMO}}^{GW}-E_\text{cut}^\text{occ},\varepsilon_{i=\text{HOMO}}^{GW}]$.
  Setting a small `ENERGY_CUTOFF_OCC` drastically reduces the computation time and the memory
  consumption, but also might affect the computed excitation energies $\Omega^{(n)}$. Recommended to
  use for large systems with more than 30 atoms, but we recommend a careful convergence test by
  increasing `ENERGY_CUTOFF_OCC` and observing the effect on $\Omega^{(n)}$ \[[Liu2020](#Liu2020)\].

- [ENERGY_CUTOFF_EMPTY](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.ENERGY_CUTOFF_EMPTY)
  $E_\text{cut}^\text{empty}$: Analogous to `ENERGY_CUTOFF_OCC`, but for the empty states, i.e. only
  empty states in the interval
  $\varepsilon_a\in[\varepsilon_{a=\text{LUMO}}^{GW},\varepsilon_{a=\text{LUMO}}^{GW}+E_\text{cut}^\text{empty}]$.

In addition to these keywords, the following settings from DFT will have an influence on the BSE
excitation energies:

- [XC_FUNCTIONAL](#CP2K_INPUT.FORCE_EVAL.DFT.XC.XC_FUNCTIONAL): Choose between one of the available
  xc-functionals. The starting point can have a profound influence on the excitation energies
  \[[Bruneval2015](#Bruneval2015), [Gui2018](#Gui2018), [Jacquemin2017](#Jacquemin2017)\]. We either
  recommend the PBE functional as DFT starting point when using BSE@ev$GW_0$@PBE or the PBE0
  functional, when using BSE@$G_0W_0$@PBE0 (see also
  [SELF_CONSISTENCY](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.SELF_CONSISTENCY)).

- [BASIS_SET](#CP2K_INPUT.FORCE_EVAL.SUBSYS.KIND.BASIS_SET): Specify the basis set, which affects
  $N_\mathrm{empty}$ and thus the size of the matrices $A_{ia,jb}$ and $B_{ia,jb}$. The
  <a href="https://www.basissetexchange.org/basis/aug-cc-pvdz/format/cp2k/?version=1&elements=1" target="_blank">`aug-cc-pVDZ`</a>
  basis set should be sufficient for most calculations, but needs to be checked regarding
  convergence, e.g. using
  <a href="https://www.basissetexchange.org/basis/aug-cc-pvtz/format/cp2k/?version=1&elements=1" target="_blank">`aug-cc-pVTZ`</a>.

The memory consumption of the BSE algorithm is large, it is approximately
$100 \cdot N_\mathrm{occ}^2 N_\mathrm{empty}^2$ Bytes. You can see $N_\mathrm{occ}$,
$N_\mathrm{empty}$ and the estimated memory consumption from the BSE output. The BSE implementation
is well parallelized, i.e. you can use several nodes that can provide the memory.

We have benchmarked the numerical precision of our BSE implementation and we found excellent
agreement within only 10 meV compared to the BSE implementation in FHI aims
\[[Liu2020](#Liu2020),[Graml2024b](#Graml2024b)\].

The current BSE implementation in CP2K works for molecules. The inclusion of periodic boundary
conditions in a $\Gamma$-only approach and with full $k$-point sampling is work in progress.

## 3. Minimal example for a BSE calculation

### 3.1 Input file

In this section, we provide a minimal example of a BSE calculation on $\mathrm{H}_2$. For the
calculation you need the input file BSE_H2.inp and the aug-cc-DZVP basis
([Download](https://www.cp2k.org/_media/howto:bse_example_h2.zip)).

Please copy both files into your working directory and run CP2K by

```none
mpirun -n 1 cp2k.psmp BSE_H2.inp
```

which requires 12 GB RAM and takes roughly 2 minutes on 1 core. You can find the output file also in
the [Download](https://www.cp2k.org/_media/howto:bse_example_h2.zip).

```none
&GLOBAL
  PROJECT  H2
  RUN_TYPE ENERGY
&END GLOBAL
&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME BASIS-aug               ! Custom Basis set file (aug-cc-pVDZ and aug-cc-pVDZ-RIFIT 
    POTENTIAL_FILE_NAME POTENTIAL               ! from the Basis Set Exchange library - www.basissetexchange.org/)
    &QS
      METHOD GAPW                               ! All electron calculation
      EPS_DEFAULT 1.0E-16
      EPS_PGF_ORB 1.0E-16
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
      &XC_FUNCTIONAL PBE0                       ! Choice of functional has a profound influence on the results
      &END XC_FUNCTIONAL
      &WF_CORRELATION
        &RI_RPA                                 ! In the RI_RPA and the GW section, additional numerical parameters, e.g.
          &GW
            SELF_CONSISTENCY      G0W0          ! can be changed to EV_GW0 or EV_GW
            &BSE  
              TDA                 TDA+ABBA      ! Diagonalizing ABBA and A
              SPIN_CONFIG         SINGLET       ! or TRIPLET
              NUM_PRINT_EXC       20            ! Number of printed excitations
              ENERGY_CUTOFF_OCC   -1            ! Set to positive numbers (eV) to
              ENERGY_CUTOFF_EMPTY -1            ! truncate matrices A_ia,jb and B_ia,jb
            &END BSE
          &END GW
        &END RI_RPA
      &END WF_CORRELATION
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC 20 20 20
      PERIODIC NONE
    &END CELL
    &COORD
      H 0.0000 0.0000 0.0000                    ! H2 molecule geometry from GW100 Paper
      H 0.0000 0.0000 0.74144
    &END COORD
    &TOPOLOGY
      &CENTER_COORDINATES
      &END
    &END TOPOLOGY
    &KIND H
      BASIS_SET ORB    aug-cc-pVDZ              ! For production runs, the basis set should be checked for convergence.
      BASIS_SET RI_AUX aug-cc-pVDZ-RIFIT        ! In general, pVDZ should be a solid choice.
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
The geometry for $\mathrm{H}_2$ was taken from the
[GW100](https://doi.org/10.1021/acs.jctc.5b00453)-Paper.

### 3.2 Output

The BSE calculation outputs the excitation energies $\Omega^{(n)}$ for *n* = 1, ...,
$N_\text{exc}^\text{print}$:

```none
 BSE| Excitation energies from solving the BSE without the TDA:
 BSE|
 BSE|     Excitation n        Multiplet  TDA/full BSE   Excitation energy Ω (eV)
 BSE|                1    Singlet State        -full-                    11.8621
 BSE|                2    Singlet State        -full-                    12.8264
 BSE|                3    Singlet State        -full-                    15.6415
 BSE|                4    Singlet State        -full-                    15.6415
```

Afterwards, the single-particle transitions, i.e. the eigenvector elements $X_{ia}^{(n)}$, with
$|X_{ia}^{(n)}|$ > [EPS_X](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.EPS_X), are
printed:

```none
 BSE| Excitations are built up by the following single-particle transitions,
 BSE| neglecting contributions where |X_ia^n| <  0.10 :
 BSE|         -- Quick reminder: HOMO i =    1 and LUMO a =    2 --
 BSE|
 BSE| Excitation n      i =>     a               TDA/full BSE           |X_ia^n|
 BSE|
 BSE|            1      1 =>     2                     -full-             0.6630
 BSE|            1      1 =>     4                     -full-             0.2549
 BSE|
 BSE|            2      1 =>     3                     -full-             0.7059
 BSE|
 BSE|            3      1 =>     6                     -full-             0.7076
 BSE|
 BSE|            4      1 =>     5                     -full-             0.7076
```

In the case of the $\mathrm{H}_2$, the lowest excitation *n* = 1 is mainly built up by a transition
from the HOMO (i=1) to the LUMO (a=2), what is apparent from
$X_{i=\text{HOMO},a=\text{LUMO}}^{(n=1)}=0.663$, and also contains a considerable contribution from
the 1=>4 (HOMO=>LUMO+2) transition ($X_{i=\text{HOMO},a=\text{LUMO+2}}^{(n=1)}=0.2549$ ).

### 3.3 Large scale calculations

Going to larger systems is a challenge for a $GW$+BSE-calculation, since the memory consumption
increases with $N_\mathrm{occ}^2 N_\mathrm{empty}^2$. The used $N_\mathrm{occ}$, $N_\mathrm{empty}$
and the required memory of a calculation are printed in the output file to help you estimating the
memory consumption and distributing the calculation on several nodes to provide the memory. In the
following, we provide a sample output of a BSE calculation on a nanographene with 206 atoms:

```none
 BSE| Cutoff occupied orbitals [eV]                                       80.000
 BSE| Cutoff empty orbitals [eV]                                          10.000
 BSE| First occupied index                                                   155
 BSE| Last empty index (not MO index!)                                       517
 BSE| Energy of first occupied index [eV]                                -24.774
 BSE| Energy of last empty index [eV]                                      7.863
 BSE| Energy difference of first occupied index to HOMO [eV]              19.487
 BSE| Energy difference of last empty index to LUMO [eV]                   9.639
 BSE| Number of GW-corrected occupied MOs                                    400
 BSE| Number of GW-corrected empty MOs                                       600
 BSE|
 BSE| Total peak memory estimate from BSE [GB]                          2221.591
 BSE| Peak memory estimate per MPI rank from BSE [GB]                      4.628
```

To enable BSE calculations on large molecules, we recommend to use a large supercomputer and setting
the keywords
[ENERGY_CUTOFF_OCC](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.ENERGY_CUTOFF_OCC)
and
[ENERGY_CUTOFF_EMPTY](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.ENERGY_CUTOFF_EMPTY),
see details given above.

## 4. References

<a id="Blase2018">\[Blase2018\]</a> X. Blase, I. Duchemin, D. Jacquemin, *The Bethe–Salpeter
equation in chemistry: relations with TD-DFT, applications and challenges*,
[Chem. Soc. Rev., **47**, 1022 (2018)](https://doi.org/10.1039/c7cs00049a).

<a id="Blase2020">\[Blase2020\]</a> X. Blase, I. Duchemin, D. Jacquemin, P.-F. Loos, *The
Bethe−Salpeter Equation Formalism: From Physics to Chemistry*,
[J. Phys. Chem. Lett., **11**, 7371−7382 (2020)](https://doi.org/10.1021/acs.jpclett.0c01875).

<a id="Bruneval2015">\[Bruneval2015\]</a> F. Bruneval, S. M. Hamed, J. B. Neaton, *A systematic
benchmark of the ab initio Bethe-Salpeter equation approach for low-lying optical excitations of
small organic molecules*,
[J. Chem. Phys. **142**, 244101 (2015)](https://doi.org/10.1063/1.4922489).

<a id="Jacquemin2017">\[Jacquemin2017\]</a> D. Jacquemin, I. Duchemin, X. Blase, *Is the
Bethe–Salpeter Formalism Accurate for Excitation Energies? Comparisons with TD-DFT, CASPT2, and
EOM-CCSD*,
[J. Phys. Chem. Lett. **8**, 1524–1529 (2017)](https://doi.org/10.1021/acs.jpclett.7b00381).

<a id="Golze2019">\[Golze2019\]</a> D. Golze, M. Dvorak, P. Rinke, *The GW Compendium: A Practical
Guide to Theoretical Photoemission Spectroscopy*,
[Front. Chem. **7**, 377 (2019)](https://doi.org/10.3389/fchem.2019.00377).

<a id="Graml2024b">\[Graml2024b\]</a> M. Graml, J. Wilhelm, Manuscript in preparation

<a id="Gui2018">\[Gui2018\]</a> X. Gui, C. Holzer, W. Klopper, *Accuracy Assessment of GW Starting
Points for Calculating Molecular Excitation Energies Using the Bethe–Salpeter Formalism*,
[J. Chem. Theory Comput., **14**, 2127-2136 (2018)](https://doi.org/10.1021/acs.jctc.8b00014).

<a id="Liu2020">\[Liu2020\]</a> C. Liu, J. Kloppenburg, Y. Yao, X. Ren, H. Appel, Y. Kanai, V. Blum,
*All-electron ab initio Bethe-Salpeter equation approach to neutral excitations in molecules with
numeric atom-centered orbitals*,
[J. Chem. Phys. **152**, 044105 (2020)](https://doi.org/10.1063/1.5123290).

<a id="Sander2015">\[Sander2015\]</a> T. Sander, E. Maggio, G. Kresse, *Beyond the Tamm-Dancoff
approximation for extended systems using exact diagonalization*,
[Phys. Rev. B **92**, 045209 (2015)](https://doi.org/10.1103/PhysRevB.92.045209).

<a id="Schreiber2008">\[Schreiber2008\]</a> M. Schreiber, M. R. Silva-Junior, S. P. A. Sauer, W.
Thiel, *Benchmarks for electronically excited states: CASPT2, CC2, CCSD, and CC3*,
[J. Chem. Phys. **7**, 134110 (2008)](https://doi.org/10.1063/1.2889385).

[bse]: #CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE
[gw]: #CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW
