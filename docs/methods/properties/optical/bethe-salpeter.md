# *GW* + Bethe-Salpeter equation

The Bethe-Salpeter equation (BSE) is a method for computing electronic excitation energies and
optical absorption spectra. We describe in Sec. [1](#header-theory) the theory and implementation of
BSE, in Sec. [2](#header-input) the BSE input keywords and in Sec. [3](#header-example) a full CP2K
input file of a BSE calculation and the corresponding output. For reviews on BSE, see
\[[](#Blase2018), [](#Blase2020), [](#Bruneval2015), [](#Sander2015)\].

(header-theory)=

## 1. Theory and implementation of BSE

A central goal of a BSE calculation is to compute electronic excitation energies
$\Omega^{(n)}, n=1,2,\ldots$.

The following ingredients are necessary for computing $\Omega^{(n)}$:

- Occupied Kohn-Sham (KS) orbitals $\varphi_i(\mathbf{r})$ and empty KS orbitals
  $\varphi_a(\mathbf{r})$ from a DFT calculation, where $i=1,\ldots,N_\mathrm{occ}$ and
  $a=N_\mathrm{occ}+1,\ldots,N_\mathrm{occ}+N_\mathrm{empty}$,
- *GW* eigenvalues $\varepsilon_i^{GW}$ and $\varepsilon_a^{GW}$ of corresponding KS orbitals.

It is possible to use *G*<sub>0</sub>*W*<sub>0</sub>, ev*GW*<sub>0</sub> or ev*GW* eigenvalues, see
details in [GW] and in Ref. \[[](#Golze2019)\], i.e. we perform
BSE@*G*<sub>0</sub>*W*<sub>0</sub>/ev*GW*<sub>0</sub>/ev*GW*@DFT. Thus, also input parameters for a
DFT and *GW* calculation influence the BSE calcualation (see full input in Sec.
[3.1](#header-input-file)). We obtain optical properties from BSE solving the following generalized
eigenvalue problem that involves the block matrix $ABBA$:

$$\left( \begin{array}{cc}A &  B\\B &  A\end{array} \right)\left( \begin{array}{cc}\mathbf{X}^{(n)}\\\mathbf{Y}^{(n)}\end{array} \right) = \Omega^{(n)}\left(\begin{array}{cc}1&0\\0&-1\end{array}\right)\left(\begin{array}{cc}\mathbf{X}^{(n)}\\\mathbf{Y}^{(n)}\end{array}\right) \quad .$$

We abbreviate $A$ and $B$ as matrices with index $A_{ia,jb}$, i.e. they have
$N_\mathrm{occ}N_\mathrm{empty}$ rows and $N_\mathrm{occ}N_\mathrm{empty}$ columns. The entries of
$A$ and

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
\[[](#Blase2020)\],

$$ \begin{align}
\Psi_\text{excitation}^{(n)}(\mathbf{r}_e,\mathbf{r}_h) = \sum_{ia} X_{ia}^{(n)} \varphi_i(\mathbf{r}_h) \varphi_a(\mathbf{r}_e) + Y_{ia}^{(n)} \varphi_i(\mathbf{r}_e) \varphi_a(\mathbf{r}_h) \quad ,
\end{align}$$

i.e. $X_{ia}^{(n)}$ and $Y_{ia}^{(n)}$ describe the transition amplitude between occupied orbital
$\varphi_i$ and empty orbital $\varphi_a$ of the $n$-th excitation.

The Tamm-Dancoff approximation (TDA) is also implemented, which constrains $B=0$. In case $A$ is
positive definite, and excitation energies $\Omega^{(n)}>0$, we have $\mathbf{Y}=0$ and $\mathbf{X}$
can be computed from the Hermitian eigenvalue problem

$$ A \mathbf{X}^{(n)}_\mathbf{TDA} = \Omega^{(n)}_\mathbf{TDA} \mathbf{X}^{(n)}_\mathbf{TDA} \quad .$$

Diagonalizing $A$ in TDA, or the full block-matrix $ABBA$, takes in the order of
$(N_\mathrm{occ} N_\mathrm{empty})^3$ floating point operations. This translates to a computational
scaling of $O(N^6)$ in the system size $N$.

(header-input)=

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
  Determines which *GW* self-consistency (*G*<sub>0</sub>*W*<sub>0</sub>, ev*GW*<sub>0</sub> or
  ev*GW*) is used to calculate the single-particle *GW* energies $\varepsilon_p^{GW}$ needed in the
  BSE calculation.

- [TDA](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.TDA): Three options available:

  - `ON` diagonalize $A$,
  - `OFF` generalized diagonalization of $ABBA$,
  - `TDA+ABBA` CP2K diagonalizes $ABBA$ as well as $A$.

- [SPIN_CONFIG](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.SPIN_CONFIG): Two options
  available: Choose between `SINGLET` for computing singlet excitation energies $(\alpha^S = 2)$ and
  `TRIPLET` for computing triplet excitation energies $(\alpha^T=0)$. Standard is `SINGLET` as an
  electronic excitation directly after photoexcitation is a singlet due to angular momentum
  conservation; triplet excited states can form by intersystem crossing.

- [NUM_PRINT_EXC](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.NUM_PRINT_EXC): Number
  of excitations $N_\text{exc}^\text{print}$ to be printed.

- [ENERGY_CUTOFF_OCC](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.ENERGY_CUTOFF_OCC)
  $E_\text{cut}^\text{occ}$: Restrict occupied molecular orbital (MO) indices $i$ and only use
  occupied MOs with
  $\varepsilon_i\in[\varepsilon_{i=\text{HOMO}}^{GW}-E_\text{cut}^\text{occ},\varepsilon_{i=\text{HOMO}}^{GW}]$.
  Setting a small `ENERGY_CUTOFF_OCC` drastically reduces the computation time and the memory
  consumption, but also might affect the computed excitation energies $\Omega^{(n)}$. Recommended to
  use for large systems with more than 30 atoms, but we recommend a careful convergence test by
  increasing `ENERGY_CUTOFF_OCC` and observing the effect on $\Omega^{(n)}$ \[[](#Liu2020)\].

- [ENERGY_CUTOFF_EMPTY](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.ENERGY_CUTOFF_EMPTY)
  $E_\text{cut}^\text{empty}$: Analogous to `ENERGY_CUTOFF_OCC`, but for the empty states, i.e. only
  empty states in the interval
  $\varepsilon_a\in[\varepsilon_{a=\text{LUMO}}^{GW},\varepsilon_{a=\text{LUMO}}^{GW}+E_\text{cut}^\text{empty}]$.

The following settings from DFT will also have an influence on the BSE excitation energies:

- [XC_FUNCTIONAL](#CP2K_INPUT.FORCE_EVAL.DFT.XC.XC_FUNCTIONAL): Choose between one of the available
  xc-functionals. The starting point can have a profound influence on the excitation energies
  \[[Bruneval2015](#Bruneval2015), [](#Gui2018), [](#Jacquemin2017)\]. We either recommend the PBE
  functional as DFT starting point when using BSE@ev*GW*<sub>0</sub>@PBE or the PBE0 functional,
  when using BSE@*G*<sub>0</sub>*W*<sub>0</sub>@PBE0 (see also
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
agreement within only 10 meV compared to the BSE implementation in FHI aims \[[](#Liu2020)\].

The current BSE implementation in CP2K works for molecules. The inclusion of periodic boundary
conditions in a Γ-only approach and with full *k*-point sampling is work in progress.

(header-example)=

## 3. Minimal example for a BSE calculation

(header-input-file)=

### 3.1 Input file

In this section, we provide a minimal example of a BSE calculation on H<sub>2</sub>. For the
calculation you need the input file BSE_H2.inp and the aug-cc-DZVP basis
([Download](https://github.com/cp2k/cp2k-examples/tree/master/bethe-salpeter/H2)).

Please copy both files into your working directory and run CP2K by

```none
mpirun -n 1 cp2k.psmp BSE_H2.inp
```

which requires 12 GB RAM and takes roughly 2 minutes on 1 core. You can find the input and output
file [here](https://github.com/cp2k/cp2k-examples/tree/master/bethe-salpeter/H2). We use the basis
sets `aug-cc-pVDZ` and `aug-cc-pVDZ-RIFIT` from the file `BASIS-aug`. These basis sets can be
obtained from the Basis Set Exchange Library:
<a href="https://www.basissetexchange.org/basis/aug-cc-pvdz/format/cp2k/?version=1&elements=1" target="_blank">aug-cc-pVDZ</a>,
<a href="https://www.basissetexchange.org/basis/aug-cc-pvdz-rifit/format/cp2k/?version=1&elements=1" target="_blank">aug-cc-pVDZ-RIFIT</a>.
The geometry for H<sub>2</sub> was taken from \[[](#vanSetten2015)\].

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

In the case of H<sub>2</sub>, the lowest excitation *n* = 1 is mainly built up by a transition from
the HOMO (i=1) to the LUMO (a=2), what is apparent from
$X_{i=\text{HOMO},a=\text{LUMO}}^{(n=1)}=0.663$, and also contains a considerable contribution from
the 1=>4 (HOMO=>LUMO+2) transition ($X_{i=\text{HOMO},a=\text{LUMO+2}}^{(n=1)}=0.2549$ ).

### 3.3 Large scale calculations

Going to larger systems is a challenge for a *GW*+BSE-calculation, since the memory consumption
increases with $N_\mathrm{occ}^2 N_\mathrm{empty}^2$. The used $N_\mathrm{occ}$, $N_\mathrm{empty}$
and the required memory of a calculation are printed in the output file to estimate the memory
consumption. We recommend to use several nodes to provide the required memory. In the following, we
provide a sample output of a BSE calculation on a nanographene with 206 atoms which has a memory
requirement of 2.2 TB RAM:

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

[bse]: #CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE
[gw]: #CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW
