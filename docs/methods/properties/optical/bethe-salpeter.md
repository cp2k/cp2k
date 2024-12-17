# *GW* + Bethe-Salpeter equation

The Bethe-Salpeter equation (BSE) is a method for computing electronic excitation energies and
optical absorption spectra. We describe in Sec. [1](#header-theory) the theory and implementation of
BSE, in Sec. [2](#header-input) the BSE input keywords and in Sec. [3](#header-example) a full CP2K
input file of a BSE calculation and the corresponding output. For reviews on BSE, see
\[[](#Blase2018), [](#Blase2020), [](#Bruneval2015), [](#Sander2015)\].

(header-theory)=

## 1. Theory and implementation of BSE

### 1.1 Electronic excitation energies and the excitation wave function

A central goal of a BSE calculation is to compute excited state properties such as the electronic
excitation energies $\Omega^{(n)} (n=1,2,\ldots)$.

The following ingredients are necessary for such a computation:

- Occupied Kohn-Sham (KS) orbitals $\varphi_i(\mathbf{r})$ and empty KS orbitals
  $\varphi_a(\mathbf{r})$ from a DFT calculation, where $i=1,\ldots,N_\mathrm{occ}$ and
  $a=N_\mathrm{occ}+1,\ldots,N_\mathrm{occ}+N_\mathrm{empty}$,
- *GW* eigenvalues $\varepsilon_i^{GW}$ and $\varepsilon_a^{GW}$ of corresponding KS orbitals.

It is possible to use *G*<sub>0</sub>*W*<sub>0 </sub>, ev*GW*<sub>0</sub> or ev*GW* eigenvalues, see
details in [GW] and in Ref. \[[](#Golze2019)\], i.e. we perform BSE@*G* <sub>0
</sub>*W*<sub>0</sub>/ev*GW*<sub>0</sub>/ev*GW*@DFT. Thus, also input parameters for a DFT and *GW*
calculation influence the BSE calculation (see full input in Sec. [3.1](#header-input-file)). We
obtain the properties of the excited states from the BSE by solving the following generalized
eigenvalue problem that involves the block matrix $ABBA$:

$$
\left( \begin{array}{cc}A &  B\\B &  A\end{array} \right)\left( \begin{array}{cc}\mathbf{X}^{(n)}\\\mathbf{Y}^{(n)}\end{array} \right) = \Omega^{(n)}\left(\begin{array}{cc}1&0\\0&-1\end{array}\right)\left(\begin{array}{cc}\mathbf{X}^{(n)}\\\mathbf{Y}^{(n)}\end{array}\right) \quad .
$$

We abbreviate $A$ and $B$ as matrices with index $A_{ia,jb}$, i.e. they have
$N_\mathrm{occ}N_\mathrm{empty}$ rows and $N_\mathrm{occ}N_\mathrm{empty}$ columns. The entries of
$A$ and

$$
\begin{align}
    A_{ia,jb} &= (\varepsilon_a^{GW}-\varepsilon_i^{GW})\delta_{ij}\delta_{ab} + \alpha^\mathrm{S/T}
    v_{ia,jb} - W_{ij,ab}(\omega=0) \quad ,\\
    B_{ia,jb} &= \alpha^\mathrm{(S/T)} v_{ia,bj} - W_{ib,aj}(\omega=0) \quad .
\end{align}
$$

where $\delta_{ij}$ is the Kronecker delta. The user sets $\alpha^S=2$ for computing singlet
excitations and $\alpha^T=0$ for computing triplet excitations. $v_{pq,rs}$ is the bare Coulomb
interaction and $W_{pq,rs}(\omega=0)$ the statically ($\omega=0$) screened Coulomb interaction,
where $p,q,r,s \in [ 1, N_\mathrm{occ}+N_\mathrm{empty}]$ are KS orbital indices. Note here, that
the screened Coulomb interaction is always computed from DFT quantitites, i.e. $W_{0}(\omega=0)$
enters the BSE. $(\mathbf{X}^{(n)},\mathbf{Y}^{(n)})$ with elements $X_{ia}^{(n)}$ and
$Y_{ia}^{(n)}$ are the eigenvectors of the excitation $n$, which relate to the wave function of the
electronic excitation \[[](#Blase2020)\],

$$
\begin{align}
\Psi_\text{excitation}^{(n)}(\mathbf{r}_e,\mathbf{r}_h) = \sum_{ia} X_{ia}^{(n)} \varphi_i(\mathbf{r}_h) \varphi_a(\mathbf{r}_e) + Y_{ia}^{(n)} \varphi_i(\mathbf{r}_e) \varphi_a(\mathbf{r}_h) \quad ,
\end{align}
$$

i.e. $X_{ia}^{(n)}$ and $Y_{ia}^{(n)}$ describe the transition amplitude between occupied orbital
$\varphi_i$ and empty orbital $\varphi_a$ of the $n$-th excitation.

The Tamm-Dancoff approximation (TDA) is also implemented, which constrains $B=0$ and
$\mathbf{Y}^{(n)}=0$. In TDA, $\Omega^{(n)}_\mathrm{TDA}$ and $\mathbf{X}^{(n)}_\mathrm{TDA}$ can be
computed from the Hermitian eigenvalue problem

$$
A \mathbf{X}^{(n)}_\mathrm{TDA} = \Omega^{(n)}_\mathrm{TDA} \mathbf{X}^{(n)}_\mathrm{TDA} \quad .
$$

Diagonalizing $A$ in TDA, or the full block-matrix $ABBA$, takes in the order of
$(N_\mathrm{occ} N_\mathrm{empty})^3$ floating point operations. This translates to a computational
scaling of $O(N^6)$ in the system size $N$.

### 1.2 Optical absorption spectrum

The BSE further allows the investigation of optical properties. For example, the optical absorption
spectrum can be computed as the imaginary part of the dynamical dipole polarizability tensor
$\alpha_{\mu,\mu'}(\omega) $ with $(\mu,\mu'\in\{x,y,z\})$:

$$
\begin{align}
\alpha_{\mu,\mu'}(\omega) 
= - \sum_n \frac{2 \Omega^{(n)} d^{(n)}_{\mu} d^{(n)}_{\mu'}}{(\omega+i\eta)^2-\left(\Omega^{(n)}\right)^2}
\quad ,
\end{align}
$$

where we have introduced an artificial broadening $\eta$. The transition moments $d^{(n)}_{\mu}$ are
computed in the length gauge $(\mu\in\{x,y,z\})$ as

$$
\begin{align}
d^{(n)}_{\mu} = \sqrt{2} \sum_{i,a} \langle \varphi_i|\hat{\mu}| \varphi_a \rangle (X_{ia}^{(n)} + Y_{ia}^{(n)}) 
\quad .
\end{align}
$$

When the molecules are not aligned, e.g. for gas phase and liquids, the spatial average is
sufficient, i.e. the optical absorption spectrum can be computed as

$$
\begin{align}
\mathrm{Im}\left[\bar{\alpha}(\omega)\right] = \frac{1}{3} \sum_{\mu\in\{x,y,z\}} \mathrm{Im}\left[\alpha_{\mu,\mu}(\omega)\right]
\quad .
\end{align}
$$

We can rewrite the last equation as

$$
\begin{align}
\mathrm{Im}\left[\bar{\alpha}(\omega)\right] 
= - \mathrm{Im}\left[
  \sum_n \frac{f^{(n)}}{(\omega+i\eta)^2-\left(\Omega^{(n)}\right)^2}
  \right]
\quad .
\end{align}
$$

where we introduced the oscillator strengths $f^{(n)}$, which are defined by

$$
\begin{align}
f^{(n)} = \frac{2}{3} \Omega^{(n)} \sum_{\mu\in\{x,y,z\}} | d^{(n)}_{\mu} |^2
\quad .
\end{align}
$$

Additionally, the photoabsorption cross section tensor

$$
\begin{align}
\sigma_{\mu,\mu'}(\omega)  = \frac{4 \pi \omega}{c} \mathrm{Im}\left[\alpha_{\mu,\mu'}(\omega) \right]
\quad .
\end{align}
$$

is printed, where $c$ denotes the speed of light.

### 1.3 Visualizing excited states using Natural Transition Orbitals (NTOs)

In order to analyse the excitation wave function independent of a specific choice of the molecular
orbitals $\varphi_p(\mathbf{r})$, we can rewrite it as

$$
\begin{align}
\Psi_\text{excitation}^{(n)}(\mathbf{r}_e,\mathbf{r}_h) = 
\sum_I {\lambda_I^{(n)}} \phi_I^{(n)}(\mathbf{r}_e) \chi_I^{(n)}(\mathbf{r}_h)
\quad .
\end{align}
$$

in terms of the natural transitions orbitals (NTOs) $\phi_I^{(n)}(\mathbf{r}_e) $ and
$\chi_I^{(n)}(\mathbf{r}_h)$. Here, we introduce the idea of electron and holes: In the excitation
process, an electron leaves a hole behind in the (formerly) occupied orbitals
$\chi_I^{(n)}(\mathbf{r}_h)$ and is afterwards located in the empty orbitals
$\phi_I^{(n)}(\mathbf{r}_e) $. We call $\phi_I^{(n)}(\mathbf{r}_e) $ and
$\chi_I^{(n)}(\mathbf{r}_h)$ a NTO pair of electron and hole with a weight $\lambda_I^{(n)}$, where
$I=1,\cdots,N_b$ with $N_b=N_\mathrm{occ}+N_\mathrm{empty}$.

The weights $\lambda_I^{(n)}$ are sorted in decreasing order, i.e.:

$$\lambda_1^{(n)} \ge \lambda_2^{(n)} \ge \cdots \ge 0$$

For many cases, we observe $\lambda_1^{(n)} \gg \lambda_2^{(n)}$, which indicates a dominant
single-particle excitation character for this $n$ and allows a visual analysis of a small number of
NTO pairs.

Assuming $\lambda_1^{(n)} = 1$ and $\lambda_{I\neq 1}=0$, the excitation wave function simply is
given as a product

$$
\begin{align}
\Psi_\text{excitation}^{(n)}(\mathbf{r}_e,\mathbf{r}_h) = 
\phi_1^{(n)}(\mathbf{r}_e) \chi_1^{(n)}(\mathbf{r}_h)
\quad .
\end{align}
$$

In this case, the electron is excited from the occupied NTO $\chi_1^{(n)}(\mathbf{r}_h)$ to the
empty state $\phi_1^{(n)}(\mathbf{r}_e)$. This process leaves a hole at $\mathbf{r}_h$ and creates
an excited electron at $\mathbf{r}_e$.

Computationally, the NTO pairs $\phi_I^{(n)}(\mathbf{r}_e) $ and $\chi_I^{(n)}(\mathbf{r}_h)$ as
well as $\lambda_I^{(n)}$ are obtained by a singular value decomposition of the transition density
matrix

$$
T = \begin{pmatrix}
    0 &  {X}^{(n)}\\
    \left({Y}^{(n)} \right)^T &  0
    \end{pmatrix} \quad ,
$$

i.e.:

$$
\begin{align}
    {T}^{(n)} &=  
    {U}^{(n)} 
    {\Lambda^{(n)}}
    \left({V}^{(n)}\right)^T
    \\
    \phi_I^{(n)}(\mathbf{r}_e) &= \sum_{p=1}^{N_b} \varphi_p(\mathbf{r}_e) V_{p,I}^{(n)} \quad ,
    \\
    \chi_I^{(n)}(\mathbf{r}_h) &= \sum_{q=1}^{N_b} \varphi_q(\mathbf{r}_h) U_{q,I}^{(n)} \quad .
\end{align}
$$

### 1.4 Measures for the size of an excited state

In this subsection, we introduce several measures, which allow us to analyse the wave function of an
excited state $\Psi_\text{excitation}^{(n)}(\mathbf{r}_e,\mathbf{r}_h)$ and its spatial properties,
following Ref. \[[](#Mewes2018)\].

To that end, we define the exciton expectation value with respect to a generic operator $\hat{O}$ as

$$
\begin{align}
\langle \hat{O} \rangle _\text{exc}^{(n)}
=
\frac{ 
 \langle \Psi_\text{excitation}^{(n)} | \hat{O} | \Psi_\text{excitation}^{(n)}\rangle 
}{
 \langle \Psi_\text{excitation}^{(n)} | \Psi_\text{excitation}^{(n)}\rangle 
}
\quad ,
\end{align}
$$

where we drop the excitation index $n$ from now on for better readability.

For each excitation level $n$, we have then several quantities, which allow us to quantify the
spatial extent of the electron, the hole and their combined two-particle character as combined
electron-hole pair, i.e. as an exciton. By that, we can often determine the type of the excited
state, i.e. distinguish between, e.g., valence, Rydberg or charge-transfer states
\[[](#Mewes2018)\].

First, we define the distance between electron and hole as

$$
\begin{align}
d_{h \rightarrow e} = | \langle \mathbf{r}_h - \mathbf{r}_e \rangle_\mathrm{exc} | \quad ,
\end{align}
$$

which can be used to distinguish different classes of excitations: For example in a charge-transfer
state, electron and hole sit on different parts of the molecule and therefore have a non-vanishing
electron-hole distance $d_{h \rightarrow e}$.

Further, we can measure the size of electron and hole, respectively, as

$$
\begin{align}
\sigma_{e/h} = \sqrt{ 
  \langle \mathbf{r}_{e/h}^2 \rangle_\mathrm{exc} 
  - \langle \mathbf{r}_{e/h} \rangle_\mathrm{exc} ^2
  } \quad ,
\end{align}
$$

which allow us to distinguish between Rydberg states, where $\sigma_h \ll \sigma_e$, and valence
states, where $\sigma_h \approx \sigma_e$.

Closely related to these quantities, we can also define the exciton size

$$
\begin{align}
d_\mathrm{exc} = \sqrt{ \langle |\mathbf{r}_h - \mathbf{r}_e|^2 \rangle_\mathrm{exc} } \quad .
\end{align}
$$

which quantifies the spatial extent of the combined electron-hole pair. As one would expect, the
exciton size $d_\mathrm{exc}$ increases when $d_{h \rightarrow e}$, $\sigma_{e}$ or $\sigma_{h}$
increase.

Finally, we quantify the correlation of electron and hole by the electron-hole correlation
coefficient

$$
\begin{align}
R_{eh} = \frac{1}{\sigma_e \sigma_h} \left( \langle \mathbf{r}_h \cdot \mathbf{r}_e \rangle_\mathrm{exc}
- \langle \mathbf{r}_h \rangle_\mathrm{exc} \cdot \langle \mathbf{r}_e \rangle_\mathrm{exc} \right)
\quad ,
\end{align}
$$

which allows us to distinguish between correlated ($R_{eh}>0$) motion, i.e. bound excitons, and
anticorrelated ($R_{eh}<0$) motion, where electron and hole try to avoid each other.

(header-input)=

## 2. BSE input

For starting a BSE calculation one needs to set the [RUN_TYPE](#CP2K_INPUT.GLOBAL.RUN_TYPE) to
`ENERGY` and the following sections for [GW] and [BSE]:

```
&GW
  SELF_CONSISTENCY      evGW0         ! We strongly recommend to use evGW0 (and PBE) for BSE runs
  &BSE  
    TDA                 TDA+ABBA      ! Diagonalizing ABBA and A
    SPIN_CONFIG         SINGLET       ! or TRIPLET
    NUM_PRINT_EXC       15            ! Number of printed excitations
    ENERGY_CUTOFF_OCC   -1            ! Set to positive numbers (eV) to
    ENERGY_CUTOFF_EMPTY -1            ! truncate matrices A_ia,jb and B_ia,jb
    NUM_PRINT_EXC_DESCR -1            ! Number of printed exciton descriptors, -1 to default to NUM_PRINT_EXC
    &BSE_SPECTRUM                     ! Activates computation and output of optical absorption spectrum
      ETA_LIST 0.01 0.02              ! Multiple broadenings can be specified within one run
    &END BSE_SPECTRUM
    &NTO_ANALYSIS
      STATE_LIST 1 2                  ! List of states for which NTOs are computed
      CUBE_FILES T                    ! Write cube files for NTOs
      STRIDE 6 6 6                    ! Coarse grid for smaller example cube files
    &END NTO_ANALYSIS
  &END BSE
&END GW
```

In the upper GW/BSE section, the following keywords have been used:

- [SELF_CONSISTENCY](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.SELF_CONSISTENCY):
  Determines which *GW* self-consistency (*G*<sub>0</sub>*W*<sub>0</sub>, ev*GW*<sub>0</sub> or
  ev*GW*) is used to calculate the single-particle *GW* energies $\varepsilon_p^{GW}$ needed in the
  BSE calculation. We recommend to always use ev*GW*<sub>0</sub> for BSE runs. The screened coulomb
  interaction is always computed from the DFT level, i.e. $W_0(\omega=0)$ is used for the
  $ABBA$-matrix, including BSE@ev*GW*@DFT.

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

- [NUM_PRINT_EXC_DESCR](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.NUM_PRINT_EXC_DESCR):
  Number of excitations, for which the exciton descriptors are printed.

- [BSE_SPECTRUM](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.BSE_SPECTRUM): Activates
  computation and printing of the optical absorption spectrum. For each chosen option in
  [TDA](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.TDA), one file is printed with
  columns frequency $\omega$, $\mathrm{Im}\left[\bar{\alpha}(\omega)\right]$ and the imaginary part
  of the elements of the dynamical dipole polarizability tensor
  $\mathrm{Im}\left[{\alpha_{\mu,\mu'}}(\omega)\right]$ as well as another file with the respective
  entries for the photoabsorption cross section tensor ${\sigma_{\mu,\mu'}}(\omega)$. The frequency
  range, step size and the broadening $\eta$ can be specified by the user (cf. keywords in
  [BSE_SPECTRUM](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.BSE_SPECTRUM)). Further,
  multiple broadenings $\eta$ can be given for one cp2k run
  ([ETA_LIST](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.BSE_SPECTRUM.ETA_LIST)),
  resulting in one file per specified $\eta$ and
  [TDA](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.TDA) (see Sec.
  [3.2](#header-output)).

- [NTO_ANALYSIS](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.NTO_ANALYSIS): Activates
  computation and printing of the natural transition orbitals. Within this section, the user can
  specify the NTO pairs, which should be printed, e.g. by the excitation indices of the desired
  levels $n$, their minimum oscillator strength $f^{(n)}$ or by the minimum weight of their NTO
  coefficient $\lambda_I$. For each NTO pair $\phi_I^{(n)}(\mathbf{r}_e) $ and
  $\chi_I^{(n)}(\mathbf{r}_h)$ of a certain excitation $n$, which fulfills these requirements, and
  the chosen options in [TDA](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.TDA), the
  NTO coefficients $\lambda_I$ are computed and the corresponding NTOs are printed to separate files
  if
  [CUBE_FILES](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.NTO_ANALYSIS.CUBE_FILES)
  is activated.

The following settings from DFT will also have an influence on the BSE excitation energies and
optical properties:

- [XC_FUNCTIONAL](#CP2K_INPUT.FORCE_EVAL.DFT.XC.XC_FUNCTIONAL): Choose between one of the available
  xc-functionals. The starting point can have a profound influence on the excitation energies
  \[[](#Knysh2024)\]. Motivated by the discussion in \[[](#Schambeck2024)\], we strongly recommend
  to use BSE@ev*GW*<sub>0</sub>@PBE, i.e. the PBE functional as DFT starting point (see also
  [SELF_CONSISTENCY](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.SELF_CONSISTENCY)).
- [BASIS_SET](#CP2K_INPUT.FORCE_EVAL.SUBSYS.KIND.BASIS_SET): Specify the basis set, which affects
  $N_\mathrm{empty}$ and thus the size of the matrices $A_{ia,jb}$ and $B_{ia,jb}$. The
  <a href="https://www.basissetexchange.org/basis/aug-cc-pvdz/format/cp2k/?version=1&elements=1" target="_blank">`aug-cc-pVDZ`</a>
  basis set should be sufficient for most calculations, but needs to be checked regarding
  convergence, e.g. using
  <a href="https://www.basissetexchange.org/basis/aug-cc-pvtz/format/cp2k/?version=1&elements=1" target="_blank">`aug-cc-pVTZ`</a>
  .

The memory consumption of the BSE algorithm is large, it is approximately
$100 \cdot N_\mathrm{occ}^2 N_\mathrm{empty}^2$ Bytes. You can see $N_\mathrm{occ}$,
$N_\mathrm{empty}$ and the estimated memory consumption from the BSE output. The BSE implementation
is well parallelized, i.e. you can use several nodes that can provide the memory.

We have benchmarked the numerical precision of our BSE implementation and compared its results to
the BSE implementation in FHI aims \[[](#Liu2020)\]. For our recommended settings, i.e.
BSE@ev*GW*<sub>0</sub>@PBE with the aug-cc-pVDZ basis set, we have found excellent agreement with
less than 5 meV mean absolute deviation averaged over the first 10 excitation levels and the 28
molecules in *Thiel's set* for Singlet excitations.

The current BSE implementation in CP2K works for molecules. The inclusion of periodic boundary
conditions in a Γ-only approach and with full *k*-point sampling is work in progress.

(header-example)=

## 3. Minimal example for a BSE calculation

(header-input-file)=

### 3.1 Input file

In this section, we provide a minimal example of a BSE calculation on H<sub>2</sub>. For the
calculation you need the input file BSE_H2.inp and the aug-cc-pVDZ basis
([Download](https://github.com/cp2k/cp2k-examples/tree/master/bethe-salpeter/H2)).

Please copy both files into your working directory and run CP2K by

```none
mpirun -n 1 cp2k.psmp BSE_H2.inp
```

which requires 5 GB RAM and takes roughly 45 seconds on 1 core. You can find the input and output
file [here](https://github.com/cp2k/cp2k-examples/tree/master/bethe-salpeter/H2). We use the basis
sets `aug-cc-pVDZ` and `aug-cc-pVDZ-RIFIT` from the file `BASIS-aug`. These basis sets can be
obtained from the Basis Set Exchange Library:
<a href="https://www.basissetexchange.org/basis/aug-cc-pvdz/format/cp2k/?version=1&elements=1" target="_blank">`aug-cc-pVDZ`</a>,
<a href="https://www.basissetexchange.org/basis/aug-cc-pvdz-rifit/format/cp2k/?version=1&elements=1" target="_blank">`aug-cc-pVDZ-RIFIT`</a>.
The geometry for H<sub>2</sub> was taken from \[[](#vanSetten2015)\].

(header-output)=

### 3.2 Output

The BSE calculation outputs the excitation energies $\Omega^{(n)}$ for *n* = 1, ...,
$N_\text{exc}^\text{print}$:

```none
BSE| Excitation energies from solving the BSE without the TDA:
BSE|
BSE|     Excitation n   Spin Config       TDA/ABBA   Excitation energy Ω^n (eV)
BSE|                1       Singlet         -ABBA-                      12.0361
BSE|                2       Singlet         -ABBA-                      13.0327
BSE|                3       Singlet         -ABBA-                      15.8729
BSE|                4       Singlet         -ABBA-                      15.8729
```

Afterwards, the single-particle transitions, i.e. the eigenvector elements $X_{ia}^{(n)}$ and
$Y_{ia}^{(n)}$, with $|X_{ia}^{(n)}|$ >
[EPS_X](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.EPS_X) or $|Y_{ia}^{(n)}|$ >
[EPS_X](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.EPS_X), are printed. The third
column is either

- $\Rightarrow$: for excitations, i.e. entries of X_ia^n, or
- $\Leftarrow$: for deexcitations, i.e. entries of Y_ia^n.

```none
BSE| The following single-particle transitions have significant contributions,
BSE| i.e. |X_ia^n| > 0.100 or |Y_ia^n| > 0.100, respectively :
BSE|         -- Quick reminder: HOMO i =    1 and LUMO a =    2 --
BSE|
BSE| Excitation n           i =>/<=     a      TDA/ABBA       |X_ia^n|/|Y_ia^n|
BSE|
BSE|            1           1    =>     2        -ABBA-                  0.6687
BSE|            1           1    =>     4        -ABBA-                  0.2431
BSE|
BSE|            2           1    =>     3        -ABBA-                  0.7061
BSE|
BSE|            3           1    =>     5        -ABBA-                  0.7077
BSE|
BSE|            4           1    =>     6        -ABBA-                  0.7077
```

In the case of H<sub>2</sub>, the lowest excitation *n* = 1 is mainly built up by a transition from
the HOMO (i=1) to the LUMO (a=2), what is apparent from
$X_{i=\text{HOMO},a=\text{LUMO}}^{(n=1)}= 0.6687$, and also contains a considerable contribution
from the 1=>4 (HOMO=>LUMO+2) transition ($X_{i=\text{HOMO},a=\text{LUMO+2}}^{(n=1)}=0.2431$ ).

In the next section, the optical properties, i.e. the transition moments $d^{(n)}_{r}$ and the
oscillator strengths $f^{(n)}$, are printed:

```none
BSE| Optical properties from solving the BSE without the TDA:
BSE|
BSE|  Excitation n  TDA/ABBA        d_x^n     d_y^n     d_z^n Osc. strength f^n
BSE|             1    -ABBA-        0.000    -0.000    -1.191             0.418
BSE|             2    -ABBA-        0.000     0.000     0.000             0.000
BSE|             3    -ABBA-       -0.836     0.713    -0.000             0.469
BSE|             4    -ABBA-       -0.713    -0.836    -0.000             0.469
```

We can see that $n=1,3,4$ are optically active since they have a nonzero oscillator strength
$f^{(n)}$.

In case [BSE_SPECTRUM](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.BSE_SPECTRUM) is
invoked, the frequency-resolved optical absorption spectrum is printed in addition. For each
specified $\eta$ in
[ETA_LIST](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.BSE_SPECTRUM.ETA_LIST), a
separate file, e.g. `BSE-ABBA-eta=0.010.spectrum`, is created, which starts with

```none
# Imaginary part of polarizability α_{µ µ'}(ω) =  \sum_n [2 Ω^n d_µ^n d_µ'^n] / [(ω+iη)²- (Ω^n)²] from Bethe Salpeter equation for method -ABBA-
#     Frequency (eV)       Im{α_{avg}(ω)}          Im{α_xx(ω)}          Im{α_xy(ω)}          Im{α_xz(ω)}          Im{α_yx(ω)}          Im{α_yy(ω)}          Im{α_yz(ω)}          Im{α_zx(ω)}          Im{α_zy(ω)}          Im{α_zz(ω)}
          0.00000000           0.00000000           0.00000000           0.00000000           0.00000000           0.00000000           0.00000000           0.00000000           0.00000000           0.00000000           0.00000000
          0.10000000           0.00005283           0.00003298          -0.00000000           0.00000000          -0.00000000           0.00003298           0.00000000           0.00000000           0.00000000           0.00009251
```

The next table summarizes some of the exciton descriptors of each excitation level $n$:

```none
BSE| Exciton descriptors from solving the BSE without the TDA:
BSE|
BSE|    n      c_n     d_eh [Å]     σ_e [Å]     σ_h [Å]   d_exc [Å]        R_eh
BSE|    1    1.026       0.0000      2.0542      0.8576      2.2521     -0.0176
BSE|    2    1.004       0.0000      2.9922      0.8578      3.1127     -0.0000
BSE|    3    1.005       0.0000      1.5859      0.8568      1.8134     -0.0077
BSE|    4    1.005       0.0000      1.5859      0.8568      1.8134     -0.0077
```

where we obtain the normalization factor
$c_n = \langle \Psi_\text{excitation}^{(n)} | \Psi_\text{excitation}^{(n)}\rangle $, the vectorial
distance between electron and hole $d_{e\rightarrow h}$, sizes of electron $\sigma_e$ and hole
$\sigma_h$, the exciton size $d_\mathrm{exc}$ and the electron-hole correlation coefficient
$R_{eh}$.

Further, the expectation values $\langle \mathbf{r}_e \rangle_\mathrm{exc}$ and
$\langle \mathbf{r}_h \rangle_\mathrm{exc}$ as well as the reference point of origin and the atomic
geometry are printed.

The last section summarizes the characteristics of the NTO analysis:

```none
BSE| Number of excitations, for which NTOs are computed                       2
BSE| Threshold for oscillator strength f^n                                  ---
BSE| Threshold for NTO weights (λ_I)^2                                    0.010
BSE|
BSE| NTOs from solving the BSE without the TDA:
BSE|
BSE|  Excitation n  TDA/ABBA   Index of NTO I               NTO weights (λ_I)^2
BSE|
BSE|             1    -ABBA-                1                           1.01320
BSE|             1    -ABBA-                2                           0.01320
BSE|
BSE|             2    -ABBA-                1                           1.00210
```

In the input file, we have specified a list of three excitation levels, no constraint on the
oscillator strengths $f^{(n)}$ and used the default threshold on the weights, given by
$\lambda_I^2 \leq 0.010$, which is repeated at the beginning of the section. Therefore, we obtain
two NTO pairs for $n=1$ with significant weight and one NTO pair for $n=2$. For $n=1$, we can verify
that the NTO weights satisfy
$1.01320 + 0.01320 \approx \sum_I {\left( \lambda_I^{(n)}\right)}^2 = c_n \approx 1.026$ when
comparing to the output of the exciton descriptor section.

### 3.3 Large scale calculations

Going to larger systems is a challenge for a *GW*+BSE-calculation, since the memory consumption
increases with $N_\mathrm{occ}^2 N_\mathrm{empty}^2$. The used $N_\mathrm{occ}$, $N_\mathrm{empty}$
and the required memory of a calculation are printed in the output file to estimate the memory
consumption. [Here](https://github.com/cp2k/cp2k-examples/tree/master/bethe-salpeter/Nanographene),
you can find a sample output of a BSE@evGW0@PBE calculation on a nanographene with 206 atoms, which
has a peak memory requirement of 2.5 TB RAM:

```none
BSE| Cutoff occupied orbitals [eV]                                       80.000
BSE| Cutoff empty orbitals [eV]                                          10.000
BSE| First occupied index                                                   155
BSE| Last empty index (not MO index!)                                       324
BSE| Energy of first occupied index [eV]                                -25.080
BSE| Energy of last empty index [eV]                                      9.632
BSE| Energy difference of first occupied index to HOMO [eV]              20.289
BSE| Energy difference of last empty index to LUMO [eV]                  10.549
BSE| Number of GW-corrected occupied MOs                                    334
BSE| Number of GW-corrected empty MOs                                       324
BSE|
BSE| Total peak memory estimate from BSE [GB]                         8.725E+02
BSE| Peak memory estimate per MPI rank from BSE [GB]                      5.453
```

To enable BSE calculations on large molecules, we recommend to use large clusters with increased RAM
and explicitly setting the keywords
[ENERGY_CUTOFF_OCC](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.ENERGY_CUTOFF_OCC)
and
[ENERGY_CUTOFF_EMPTY](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE.ENERGY_CUTOFF_EMPTY),
see details given above.

[bse]: #CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.BSE
[gw]: #CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW
