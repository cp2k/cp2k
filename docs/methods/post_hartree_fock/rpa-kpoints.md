# RPA with *k*-points (RI-RS)

For periodic systems with a small unit cell and full *k*-point sampling, CP2K computes the RPA total
energy per unit cell,

$$
E^\text{RPA} = E^\text{DFT} - E^\text{xc} + E^\text{EXX} + E_\text{c}^\text{RPA}\,,
$$

where $E^\text{DFT}$ is the DFT total energy, $E^\text{xc}$ the exchange-correlation energy of the
DFT functional, $E^\text{EXX}$ the exact exchange energy and $E_\text{c}^\text{RPA}$ the RPA
correlation energy. $E^\text{EXX}$ and $E_\text{c}^\text{RPA}$ are evaluated with the SCF orbitals
on the *k*-point mesh of [DFT/KPOINTS](#CP2K_INPUT.FORCE_EVAL.DFT.KPOINTS). The implementation
requires 3D periodicity. Nuclear gradients and the stress tensor are in preparation, as is a
publication on the method. For an introduction to RPA and for the other RPA implementations in CP2K,
see [](rpa.md).

## Real-space resolution of the identity (RI-RS)

An RPA calculation needs products of two Gaussian basis functions,
$\phi_\mu(\mathbf{r})\phi_\nu(\mathbf{r})$. In the resolution of the identity (RI), such a product
is expanded in an auxiliary basis $\{\varphi_P\}$, the RI basis. The expansion needs the
three-center integrals $(\mu\nu|P)$, and with *k*-point sampling, there are very many of them. The
real-space resolution of the identity (RI-RS) \[[](#Duchemin2019), [](#Duchemin2021)\] replaces
these three-center integrals by the values of the basis functions on a small set of real-space grid
points $\mathbf{r}_\ell$ around every atom:

$$
(\mu\nu|P) \approx \sum_\ell \phi_\mu(\mathbf{r}_\ell)\,\phi_\nu(\mathbf{r}_\ell)\, Z_{\ell P}\,.
$$

The fit coefficients $Z_{\ell P}$ are computed once. Afterwards, the density response $\chi$ and the
RPA correlation energy follow from fast matrix multiplications. CP2K builds the RI basis and the
grid points automatically, see [Input](#input). RI-RS is also used for
[*GW* band structures with *k*-points](../electronic_structure/band/gw.md#4-gw-band-structures-of-solids-with-k-point-sampling).

## Theory

The RPA correlation energy per unit cell is

$$
E_\text{c}^\text{RPA} = \frac{1}{2\pi N_k}\sum_{\mathbf{q}}\int_0^\infty d\omega\,
\text{Tr}\left[\ln\varepsilon(\mathbf{q},i\omega) - \varepsilon(\mathbf{q},i\omega) + 1\right]
$$

with the number of *k*-points $N_k$, the wavevectors $\mathbf{q}$ of the *k*-point mesh and the
dielectric matrix $\varepsilon(\mathbf{q},i\omega)$ in the RI basis at imaginary frequency
$i\omega$. As for *GW*, CP2K computes the Green's function $G^\mathbf{R}_{\ell\ell'}(\pm i\tau)$ on
the RI-RS grid points $\mathbf{r}_\ell$ for lattice vectors $\mathbf{R}$ in imaginary time $i\tau$
from the SCF orbitals. The density response is then an element-wise product, which is transformed to
the wavevectors $\mathbf{q}$ and projected onto the RI basis functions $P, Q$ with the RI-RS
coefficients $Z_{\ell P}$:

$$
\chi^\mathbf{R}_{\ell\ell'}(i\tau) = -2\, G^\mathbf{R}_{\ell\ell'}(i\tau)\,
G^{-\mathbf{R}}_{\ell'\ell}(-i\tau)\,,
$$

$$
\chi_{PQ}(\mathbf{q},i\tau) = \sum_{\ell\ell'} Z^*_{\ell P}(\mathbf{q})
\left[\sum_\mathbf{R} e^{i\mathbf{q}\cdot\mathbf{R}}\, \chi^\mathbf{R}_{\ell\ell'}(i\tau)\right]
Z_{\ell' Q}(\mathbf{q})\,.
$$

With this sign convention, $\chi_{PQ}(\mathbf{q},i\tau)$ is a positive semidefinite matrix. A cosine
transform on the minimax grid gives $\chi_{PQ}(\mathbf{q},i\omega)$ and the dielectric matrix

$$
\varepsilon(\mathbf{q},i\omega) = 1 + V^{1/2}(\mathbf{q})\, M^{-1}(\mathbf{q})\,
\chi(\mathbf{q},i\omega)\, M^{-1}(\mathbf{q})\, V^{1/2}(\mathbf{q})
$$

with the Coulomb matrix $V_{PQ}(\mathbf{q})$ and the RI metric matrix $M_{PQ}(\mathbf{q})$ in the RI
basis. The trace of the logarithm follows from the Cholesky decomposition $\varepsilon = LL^\dagger$
without diagonalization, and the frequency integral is a sum over the minimax frequencies $\omega_j$
with weights $w_j$:

$$
E_\text{c}^\text{RPA} = \frac{1}{2\pi N_k}\sum_{j}w_j\sum_{\mathbf{q}}\sum_P
\left[2\ln L_{PP}(\mathbf{q},i\omega_j) - \varepsilon_{PP}(\mathbf{q},i\omega_j) + 1\right]\,.
$$

$E^\text{EXX}$ is summed over the same $\mathbf{q}$ mesh, and the $1/q^2$ singularity of the Coulomb
interaction is treated with the correction of Gygi and Baldereschi \[Phys. Rev. B **34**, 4405
(1986)\]. The RI-RS fit needs the quantity
$d_{\ell P} = \sum_{\mu\nu} \phi_\mu(\mathbf{r}_\ell)\phi_\nu(\mathbf{r}_\ell)(\mu\nu|P)$ with
Gaussian basis functions $\phi_\mu$. CP2K computes $d_{\ell P}$ either from analytic three-center
integrals with lattice sums and space-group symmetry or by numerical integration on a uniform
real-space grid.

## Input

The calculation is activated via the
[RI_RS_RPA](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RS_RPA) section:

```none
&DFT
  ...
  &KPOINTS
    SCHEME MONKHORST-PACK 8 8 8
    PARALLEL_GROUP_SIZE   -1
  &END KPOINTS
  &XC
    &XC_FUNCTIONAL PBE
    &END XC_FUNCTIONAL
    &WF_CORRELATION
      &RI_RS_RPA
        NUM_TIME_FREQ_POINTS 10
        &AUTO_RI
          RI_AO_RATIO 4
        &END AUTO_RI
        &RI_RS
          &GRID_OPTIMIZATION
            RS_AO_RATIO 18
          &END GRID_OPTIMIZATION
        &END RI_RS
      &END RI_RS_RPA
    &END WF_CORRELATION
  &END XC
&END DFT
```

The keywords and convergence parameters are:

- [DFT/KPOINTS](#CP2K_INPUT.FORCE_EVAL.DFT.KPOINTS): The *k*-point mesh of DFT is also the
  $\mathbf{q}$ mesh of $E^\text{EXX}$ and $E_\text{c}^\text{RPA}$. It is a convergence parameter.
- [NUM_TIME_FREQ_POINTS](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RS_RPA.NUM_TIME_FREQ_POINTS):
  Number of points of the minimax imaginary-time and imaginary-frequency grids. 8 to 12 points are
  usually sufficient for semiconductors and insulators (the default of 30 is needed for *GW*; RPA
  usually needs much fewer points). For Si, $E_\text{c}^\text{RPA}$ is identical within 0.03 meV for
  12, 20 and 30 points.
- [AUTO_RI](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RS_RPA.AUTO_RI) with
  [RI_AO_RATIO](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RS_RPA.AUTO_RI.RI_AO_RATIO):
  Automatic RI basis. `RI_AO_RATIO` is the number of RI basis functions per Gaussian basis function
  (default: 3). Use 3 to 4; the benchmark below was computed with 4. Values above 5 do not improve
  the result and can lead to numerical problems, because the RI basis functions become linearly
  dependent.
- [RI_RS/GRID_OPTIMIZATION](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RS_RPA.RI_RS.GRID_OPTIMIZATION)
  with
  [RS_AO_RATIO](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RS_RPA.RI_RS.GRID_OPTIMIZATION.RS_AO_RATIO):
  Automatic RI-RS grid. `RS_AO_RATIO` is the number of RI-RS grid points per Gaussian basis
  function. As a rule of thumb, use 18 to 20; the lattice constants and bulk moduli of the benchmark
  below were computed with 18 (Si: 20). Do not increase it without need: values above about 25 can
  lead to numerical problems.
- [D_LP_COMPUTATION](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RS_RPA.RI_RS.D_LP_COMPUTATION):
  `ANALYTIC` or `NUMERIC` computation of $d_{\ell P}$. With the default `AUTO`, CP2K estimates the
  computation time of both and uses the faster one. `NUMERIC` is not feasible for all-electron basis
  sets.

The output lists the three energies in Hartree:

```none
RI-RS-RPA| Exact exchange energy E^EXX:
RI-RS-RPA| RPA correlation energy E_c^RPA:
RI-RS-RPA| Total energy E^DFT-E^xc+E^EXX+E_c^RPA:
```

## Pseudopotentials and basis sets

We recommend the GTH-HF pseudopotentials for RPA calculations. The exchange energy of RPA is the
exact exchange energy, and a pseudopotential that was generated for Hartree-Fock describes the
exchange interaction between core and valence electrons consistently with it. For Si, the RPA
lattice constant is 5.412 Å with the GTH-HF pseudopotential and 5.398 Å with the GTH-PBE
pseudopotential, while plane-wave RPA calculations give 5.414 Å to 5.438 Å (see the benchmark
below). The SCF is a PBE calculation also with the GTH-HF pseudopotential.

As basis sets, we recommend the correlation-consistent basis sets cc-pVXZ for solids of Ye and
Berkelbach ([Ye and Berkelbach 2022](https://doi.org/10.1021/acs.jctc.1c01245)), which are optimized
for the GTH-HF pseudopotentials (basis set names `cc-pVTZ-GTH-HF-q4`, `cc-pVQZ-GTH-HF-q4`, ... with
the number of valence electrons after `q`). The RPA correlation energy converges slowly with the
basis-set size: for Si, the RPA lattice constant decreases by 0.017 Å from cc-pVTZ to cc-pVQZ. All
calculations of the benchmark below use cc-pVQZ. The basis set and pseudopotential files are part of
the benchmark repository.

## Preliminary benchmark: lattice constants and bulk moduli

The tables compare lattice constants $a_0$ and bulk moduli $B_0$ of six crystals from CP2K with
all-electron and plane-wave results from other codes and with experiment (corrected for zero-point
vibrations, except the bulk modulus of MgO). The last column is the mean absolute error (MAE) with
respect to experiment. The CP2K calculations use a *k*-point mesh of 8x8x8 and the cc-pVQZ basis set
without basis-set extrapolation; PBE with GTH-PBE pseudopotentials, RPA@PBE with GTH-HF
pseudopotentials, `RI_AO_RATIO 4`, `RS_AO_RATIO 18` (Si: 20) and 20 time and frequency points. $a_0$
and $B_0$ follow from a third-order Birch-Murnaghan fit to the total energy at nine lattice
constants. One RPA calculation of Si with cc-pVQZ takes at most 24 minutes on a single compute node
with 128 cores.

PBE lattice constant (Å), MAE in mÅ:

|                       | Si        | C         | SiC       | BN        | MgO       | LiF       | MAE    |
| --------------------- | --------- | --------- | --------- | --------- | --------- | --------- | ------ |
| **CP2K, GTH-PBE**     | **5.462** | **3.573** | **4.371** | **3.629** | **4.260** | **4.072** | **49** |
| all-electron [1, 2]   | 5.476     | 3.575     | 4.384     | 3.628     | 4.254     | 4.069     | 52     |
| plane waves, PAW [1]  | 5.469     | 3.574     | 4.380     | 3.626     | –         | 4.068     | 47     |
| plane waves, ONCV [3] | 5.479     | 3.567     | –         | 3.621     | 4.255     | 4.062     | 51     |
| experiment [1, 4]     | 5.421     | 3.553     | 4.346     | 3.592     | 4.189     | 3.972     |        |

PBE bulk modulus (GPa), MAE in GPa:

|                      | Si     | C       | SiC     | BN      | MgO     | LiF    | MAE    |
| -------------------- | ------ | ------- | ------- | ------- | ------- | ------ | ------ |
| **CP2K, GTH-PBE**    | **88** | **430** | **211** | **369** | **148** | **68** | **20** |
| all-electron [1, 2]  | 89     | 435     | 213     | 374     | 149     | 67     | 18     |
| plane waves, PAW [1] | 88     | 431     | 210     | 370     | –       | 67     | 21     |
| experiment [1, 4]    | 101    | 455     | 229     | 410     | 165     | 76     |        |

RPA lattice constant (Å), MAE in mÅ:

|                         | Si        | C         | SiC       | BN        | MgO       | LiF       | MAE    |
| ----------------------- | --------- | --------- | --------- | --------- | --------- | --------- | ------ |
| **CP2K, GTH-HF**        | **5.412** | **3.580** | **4.358** | **3.631** | **4.191** | **4.003** | **20** |
| plane waves, PAW [4, 5] | 5.432     | 3.572     | 4.365     | 3.621     | 4.225     | 3.998     | 23     |
| norm-cons. PAW [6]      | 5.428     | 3.562     | 4.353     | –         | 4.200     | 3.994     | 11     |
| standard PAW [6]        | 5.438     | 3.571     | 4.367     | –         | 4.217     | 3.993     | 21     |
| plane waves, ONCV [3]   | 5.414     | 3.567     | –         | 3.617     | 4.200     | 3.996     | 16     |
| experiment [1, 4]       | 5.421     | 3.553     | 4.346     | 3.592     | 4.189     | 3.972     |        |

RPA bulk modulus (GPa), MAE in GPa:

|                         | Si     | C       | SiC     | BN      | MgO     | LiF    | MAE    |
| ----------------------- | ------ | ------- | ------- | ------- | ------- | ------ | ------ |
| **CP2K, GTH-HF**        | **97** | **430** | **218** | **376** | **173** | **89** | **16** |
| plane waves, PAW [4, 5] | 99     | 441     | 223     | 378     | 168     | 76     | 10     |
| plane waves, GPAW [7]   | 95     | 435     | 220     | –       | 163     | 71     | 8      |
| experiment [1, 4]       | 101    | 455     | 229     | 410     | 165     | 76     |        |

- [1] Tran and Blaha, Phys. Rev. B **83**, 235118 (2011).
- [2] Bosoni *et al.*, Nat. Rev. Phys. **6**, 45 (2024).
- [3] Pitts, Contant and Hellgren, Phys. Rev. B **112**, 085137 (2025).
- [4] Harl, Schimka and Kresse, Phys. Rev. B **81**, 115126 (2010).
- [5] Humer *et al.*, arXiv:2602.14923 (2026).
- [6] Klimeš, Kaltak, Maggio and Kresse, J. Chem. Phys. **143**, 102816 (2015).
- [7] Patrick and Thygesen, J. Chem. Phys. **143**, 102802 (2015).

PBE from CP2K agrees with the all-electron and plane-wave PBE results within 0.017 Å and 5 GPa. RPA
reduces the MAE of the lattice constant from 0.049 Å (PBE) to 0.020 Å, which is in the range of the
plane-wave RPA calculations (0.011 Å to 0.023 Å). The RPA lattice constants from CP2K deviate by at
most 0.034 Å (MgO) from the plane-wave RPA results. The RPA bulk modulus of LiF from CP2K is 13 GPa
above the plane-wave RPA value of [4, 5]; this deviation is under investigation.

These numbers are preliminary. Inputs, outputs and a summary of all calculations are available in
the GitHub repository
[CP2K_PR6627_RI_RS_RPA_benchmarks](https://github.com/JWilhelm/CP2K_PR6627_RI_RS_RPA_benchmarks).

## Practical recommendations

- ***k*-point mesh:** On an even *k*-point mesh, $E^\text{EXX}$ needs at least four *k*-points per
  direction. A 2x2x2 mesh gives a wrong $E^\text{EXX}$, and CP2K stops with an error. On a mesh that
  contains the $\Gamma$ point (odd mesh), the $\mathbf{q}=0$ term of $E^\text{EXX}$ is incomplete,
  which causes an error that decreases as $1/N_k$, and CP2K prints a note in the output.
- **Equation of state:** The RI-RS grid is geometry-specific, i.e. it is optimized for every
  geometry. The energy as a function of the lattice constant therefore scatters around a smooth
  curve, in the benchmark by less than 1 meV for Si, C, SiC, BN and LiF and by 3 meV for MgO. The
  RPA lattice constants changed by at most 0.003 Å between different RI and RI-RS settings. The bulk
  modulus changed by a few percent.
- **All-electron calculations:** The automatic RI basis and the automatic RI-RS grid do not work for
  all-electron basis sets because of the large range of the Gaussian exponents. All-electron
  calculations need an `RI_AUX` basis set in `&KIND` and a tabulated RI-RS grid
  ([GRID_SELECT](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RS_RPA.RI_RS.GRID_SELECT)).
  All-electron RPA energies can be computed, but $E^\text{EXX}$ scatters by 26–160 meV between
  geometries, so all-electron RPA equations of state do not work at present.
