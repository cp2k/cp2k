# Band structure from *GW*

*GW* is a method for computing the electronic band structure of solids and molecular orbital
energies of molecules. In this tutorial, we describe the band gap problem of DFT to motivate the
usage of *GW*, we briefly discuss the theoretical framework of *GW* and the *GW* implementation in
CP2K. We also provide input and output files of *GW* calculations performed with CP2K. The *GW*
theory of this tutorial is taken from \[[](#Golze2019)\] and many useful details on *GW* can be
found in this review.

## 1. The band gap problem of DFT

When considering a solid under periodic boundary conditions, the Kohn-Sham (KS) equations read

$$\left(- \frac{\nabla^2}{2m} + v_\text{ext}(\mathbf{r}) + v_\text{Hartree}(\mathbf{r}) + 
v_\text{xc}(\mathbf{r}) \right) \psi_{n\mathbf{k}}(\mathbf{r}) = 
\varepsilon_{n\mathbf{k}}^\text{DFT} \psi_{n\mathbf{k}}(\mathbf{r}). $$

The KS equations are solved to obtain the KS orbitals $\psi_{n\mathbf{k}}(\mathbf{r})$ and the KS
eigenvalues $\varepsilon_{n\mathbf{k}}^\text{DFT}$ with band index $n$ and crystal momentum
$\mathbf{k}$. For a molecule, the KS orbitals $\psi_n(\mathbf{r})$ and the KS eigenvalues
$\varepsilon_n^\text{DFT}$ do not carry a crystal momentum $\mathbf{k}$, only the molecular-orbital
index $n$.

The KS eigenvalues $\varepsilon_{n\mathbf{k}}^\text{DFT}$ are often used to approximate the
electronic band structure of a solid. This approximation comes with limitations:

- When using one of the common GGA exchange-correlation (xc) functionals, the band gap in the KS-DFT
  band structure $\varepsilon_{n\mathbf{k}}^\text{DFT}$ is much too small compared to experimental
  band gaps (Fig. 26 in \[[](#Golze2019)\]). Even with the exact xc functional, the band gap in the
  KS-DFT band structure $\varepsilon_{n\mathbf{k}}^\text{DFT}$ will be too small due to the
  derivative discontinuity.

- The GGA band structure $\varepsilon_{n\mathbf{k}}^\text{DFT}$ is insensitive to screening by the
  environment. As an example, the GGA eigenvalue gap
  $\varepsilon_{n=\text{LUMO}}^\text{DFT}-\varepsilon_{n=\text{HOMO}}^\text{DFT}$ of a molecule is
  almost identical for the molecule in gas phase and the molecule being on a surface. In experiment
  however, the surface induces an image-charge effect which can reduce the HOMO-LUMO gap of the
  molecule by several eV compared to gas phase. In more general terms, this is a non-local screening
  effect by the surface which is absent in common approximate GGA xc functionals. A similar band gap
  reduction due to non-local screening is present when bringing two materials close to each other,
  for example when stacking two sheets of atomically thin materials on top of each other.

- One might use hybrid xc functionals like HSE06 to obtain band gaps from KS eigenvalues that align
  more closely with experimental values than GGA band gaps. However, the issue remains that
  non-local screening effects by the environment are not included in hybrid functionals.

The above issues are known as the band gap problem of DFT.

## 2. Theory of *GW* band structure calculations

Green's function theory offers a framework for calculating electron removal and addition energies,
known as quasiparticle energies. Hedin’s equations provide an exact method for computing these
quasiparticle energies within Green's function theory. The *GW* approximation simplifies Hedin's
equations by approximating the self-energy Σ as the product of the Green's function *G* and the
screened Coulomb interaction *W*,

$$\Sigma^{GW}(\mathbf{r}_1,\mathbf{r}_2,t)= iG(\mathbf{r}_1,\mathbf{r}_2,t)W(\mathbf{r}_1,\mathbf{r}_2,t). $$

A big success of *GW* is that it captures non-local screening effects on the electronic band
structure as previously discussed and that band gaps computed from *GW* can be in excellent
agreement with experiment.

*GW* calculations in CP2K start from a KS-DFT calculation, i.e., we assume that the above KS
equations have been solved. In the *G*<sub>0</sub>*W*<sub>0</sub> approach, we use KS orbitals and
KS eigenvalues to compute the Green's function $G_0$ of non-interacting electrons and the screened
Coulomb interaction $W_0$ in the random-phase approximation. One then computes the
*G*<sub>0</sub>*W*<sub>0</sub> self-energy $\Sigma^{G_0W_0}(t)$ by replacing in the above equation
$G\rightarrow G_0$ and $W\rightarrow W_0$ and Fourier transforms $\Sigma^{G_0W_0}(t)$ to
frequency/energy $\varepsilon$, $\Sigma^{G_0W_0}(\varepsilon)$. We further approximate in
*G*<sub>0</sub>*W*<sub>0</sub> that KS orbitals are the quasiparticle wavefunctions. Then, the
*G*<sub>0</sub>*W*<sub>0</sub> quasiparticle energies follow,

$$ \varepsilon_{n\mathbf{k}}^{G_0W_0} = \varepsilon_{n\mathbf{k}}^\text{DFT} + \braket{\psi_{n\mathbf{k}}|
\Sigma^{G_0W_0}(\varepsilon_{n\mathbf{k}}^{G_0W_0}) - v_\text{xc}|\psi_{n\mathbf{k}}} .$$

We might interpret that we remove the spurious xc contribution from DFT,
$\braket{\psi_{n\mathbf{k}}|
v_\text{xc}|\psi_{n\mathbf{k}}}$, from the DFT eigenvalue $\varepsilon_{n\mathbf{k}}^\text{DFT}$ and
we add the xc contribution from *G*<sub>0</sub>*W*<sub>0</sub>, $\braket{\psi_{n\mathbf{k}}|
\Sigma^{G_0W_0}(\varepsilon_{n\mathbf{k}}^{G_0W_0}) |\psi_{n\mathbf{k}}}$.

CP2K also allows to perform eigenvalue-selfconsistency in $G$ (ev*GW*<sub>0</sub>). In general, we
recommend to use ev*GW*<sub>0</sub>@PBE for both molecules and solids, motivated by the discussion
in \[[](#Schambeck2024)\].

Please note that *G*<sub>0</sub>*W*<sub>0</sub> and ev*GW*<sub>0</sub> use the KS orbitals as
quasiparticle wavefunctions. This requires that the KS orbitals are reasonable. For systems with a
complicated electronic structure, for example transition-metal complexes or a defect with unpaired
spins in a solid, a hybrid functional may be needed to obtain a qualitatively correct electron
density and thus qualitatively correct quasiparticle wavefunctions. In such cases, one might need to
start the *GW* calculation from a hybrid functional.

CP2K contains two different *GW* implementations:

- *GW* for molecules \[[](#Wilhelm2016)\] (Sec. 3)
- *GW* for computing the band structure of a solid with a small unit cell with *k*-point sampling
  (\[[](#Pasquier2025)\] and publication in preparation, Sec. 4)

In the following, we will discuss details and usage of these *GW* implementations.

## 3. *GW* for molecules

For starting a *G*<sub>0</sub>*W*<sub>0</sub>, ev*GW*<sub>0</sub> or ev*GW* calculation for a
molecule, one needs to set the [RUN_TYPE](#CP2K_INPUT.GLOBAL.RUN_TYPE) to `ENERGY`. The section
[DFT/XC/WF_CORRELATION/RI_RPA/GW](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW) should be
defined:

```
&XC
 &XC_FUNCTIONAL PBE
 &END XC_FUNCTIONAL
 ! GW for molecules is part of the WF_CORRELATION section
 &WF_CORRELATION
  &RI_RPA
   ! use 100 points to perform the frequency integration in GW
   QUADRATURE_POINTS 100
   &GW
    SELF_CONSISTENCY  G0W0   ! can be changed to EV_GW0 or EV_GW
   &END GW
  &END RI_RPA
 &END WF_CORRELATION
&END XC
```

In the upper *GW* section, the following keywords have been used:

- [QUADRATURE_POINTS](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.QUADRATURE_POINTS): Number
  of imaginary-frequency points used for computing the self-energy (Eq. (21) in
  \[[](#Wilhelm2016)\]). 100 points are usually enough for converging quasiparticle energies within
  10 meV.

- [SELF_CONSISTENCY](#CP2K_INPUT.FORCE_EVAL.DFT.XC.WF_CORRELATION.RI_RPA.GW.SELF_CONSISTENCY):
  Determines which *GW* self-consistency (*G*<sub>0</sub>*W*<sub>0</sub>, ev*GW*<sub>0</sub> or
  ev*GW*) is used to calculate the *GW* quasiparticle energies. For molecules, CP2K also offers
  eigenvalue self-consistency in both $G$ and $W$ (ev*GW*): the quasiparticle energies replace the
  KS eigenvalues in $G$ and in $W$ until the quasiparticle energies no longer change.

The numerical precision of the *GW* implementation is 10 meV compared to reference calculations, for
example the *GW*100 test set \[[](#vanSetten2015)\]. To help new users get familiar with the *GW*
implementation in CP2K, we recommend starting by reproducing the HOMO and LUMO
*G*<sub>0</sub>*W*<sub>0</sub>@PBE quasiparticle energies of the H<sub>2</sub>O molecule from the
GW100 test set. The reference values are $\varepsilon_\text{HOMO}^{G_0W_0\text{@PBE}}$ = -11.97 eV
and $\varepsilon_\text{LUMO}^{G_0W_0\text{@PBE}}$ = 2.37 eV; CP2K input and output files for the
*G*<sub>0</sub>*W*<sub>0</sub>@PBE calculation of the H<sub>2</sub>O molecule are available
[here](https://github.com/cp2k/cp2k-examples/tree/master/gw/1_H2O_GW100).

The following settings from DFT will also have an influence on *GW* quasiparticle energies:

- [XC_FUNCTIONAL](#CP2K_INPUT.FORCE_EVAL.DFT.XC.XC_FUNCTIONAL): Starting xc functional for the
  *G*<sub>0</sub>*W*<sub>0</sub>, ev*GW*<sub>0</sub> or ev*GW* calculation. We recommend to use
  ev*GW*<sub>0</sub>@PBE, as discussed in \[[](#Schambeck2024)\]. For further guidance on selecting
  an appropriate DFT starting functional and self-consistency scheme for your system, you may
  consult \[[](#Golze2019)\].

- [BASIS_SET](#CP2K_INPUT.FORCE_EVAL.SUBSYS.KIND.BASIS_SET): The basis set is of Gaussian type and
  can have strong influence on the quasiparticle energies. For computing quasiparticle energies, a
  basis set extrapolation is necessary, see Fig. 2a in \[[](#Wilhelm2016)\] and we recommend
  all-electron GAPW calculations with correlation-consistent basis sets from the Basis Set Exchange
  Library. For computing the HOMO-LUMO gap from *GW*, we recommend augmented basis sets, for example
  [`aug-cc-pVDZ`](https://www.basissetexchange.org/basis/aug-cc-pvdz/format/cp2k/?version=1&elements=1,6,8)
  and
  [`aug-cc-pVTZ`](https://www.basissetexchange.org/basis/aug-cc-pvtz/format/cp2k/?version=1&elements=1,6,8).
  As for `RI_AUX` basis set, we recommend the corresponding RIFIT basis sets, for example
  [`aug-cc-pVDZ-RIFIT`](https://www.basissetexchange.org/basis/aug-cc-pvdz-rifit/format/cp2k/?version=1&elements=1,6,8).

The computational effort of the *GW* calculation increases with *N*<sup>4</sup> in system size *N*.
The memory requirement increases with *N*<sup>3</sup>. For running large-scale calculations, we
recommend starting with a small molecule. After successfully completing the *GW* calculation for the
small molecule, you can gradually increase the molecule size. The computational resources needed for
larger molecules can then be estimated using the *N*<sup>4</sup> scaling for computation time and
*N*<sup>3</sup> scaling for memory. The output provides a useful lower limit of the required memory:

```
  RI_INFO| Total memory for (ia|K) integrals:                       90555.51 MiB
  RI_INFO| Total memory for G0W0-(nm|K) integrals:                   4625.99 MiB
```

When facing out-of-memory, please increase the number of nodes of your calculation.

Input and output of a large-scale *GW* calculation on a nanographene with 200 atoms is available
[here](https://github.com/cp2k/cp2k-examples/tree/master/gw/2_Nanographene_206_atoms_N%5E4_scaling).

## 4. *GW* band structures of solids with *k*-point sampling

CP2K computes *GW* band structures of 2D materials (`PERIODIC XY`) and of 3D bulk crystals
(`PERIODIC XYZ`) with a small unit cell. The calculation starts from a DFT calculation with
*k*-point sampling and then computes the *G*<sub>0</sub>*W*<sub>0</sub> or ev*GW*<sub>0</sub>
quasiparticle energies $\varepsilon_{n\mathbf{k}}^{GW}$ along a path in the Brillouin zone.

### 4.1 Real-space resolution of the identity (RI-RS)

A *GW* calculation needs products of two Gaussian basis functions, $\phi_\mu(\mathbf{r})
\phi_\nu(\mathbf{r})$. In the resolution of the identity (RI), such a product is expanded in an
auxiliary basis $\{\varphi_P\}$, the RI basis. The expansion needs the three-center integrals
$(\mu\nu|P)$, and with *k*-point sampling, there are very many of them. The real-space resolution of
the identity (RI-RS) \[[](#Duchemin2019), [](#Duchemin2021)\] replaces these three-center integrals
by the values of the basis functions on a small set of real-space grid points $\mathbf{r}_\ell$
around every atom:

$$ (\mu\nu|P) \approx \sum_\ell \phi_\mu(\mathbf{r}_\ell)\,\phi_\nu(\mathbf{r}_\ell)\,Z_{\ell P}\,. $$

The fit coefficients $Z_{\ell P}$ are computed once. Afterwards, the density response $\chi$, the
screened Coulomb interaction $W$ and the self-energy $\Sigma$ follow from fast matrix
multiplications.

Users do not need to construct the RI basis or the grid points themselves. CP2K builds both
automatically for the material:

- **Automatic RI basis** ([AUTO_RI](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.BANDSTRUCTURE.GW.AUTO_RI)):
  CP2K starts from a large RI basis and compresses it for every atom to the
  [RI_AO_RATIO](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.BANDSTRUCTURE.GW.AUTO_RI.RI_AO_RATIO)-fold number
  of the atom's basis functions. We recommend `RI_AO_RATIO 3` (this is also the default).
- **Automatic RI-RS grid**
  ([RI_RS/GRID_OPTIMIZATION](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.BANDSTRUCTURE.GW.RI_RS.GRID_OPTIMIZATION)):
  CP2K optimizes the positions of the grid points around every atom such that the RI-RS formula
  above reproduces the exact three-center integrals. Every atom gets
  [RS_AO_RATIO](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.BANDSTRUCTURE.GW.RI_RS.GRID_OPTIMIZATION.RS_AO_RATIO)
  grid points per basis function. We recommend `RS_AO_RATIO 7` (this is also the default).

### 4.2 *k*-point sampling

*GW* uses the *k*-point mesh of [DFT/KPOINTS](#CP2K_INPUT.FORCE_EVAL.DFT.KPOINTS) with $N_k$
*k*-points; there is no separate *k*-point setting for *GW*. From the KS orbitals
$\psi_{n\mathbf{k}}$ and KS eigenvalues $\varepsilon_{n\mathbf{k}}^\text{DFT}$ on this mesh, CP2K
computes the Green's function on the RI-RS grid points $\mathbf{r}_\ell$ in imaginary time $i\tau$
(here for $\tau>0$, sum over empty bands $n$, Fermi energy $\varepsilon_\text{F}$) and transforms it
to the lattice vectors $\mathbf{R}$:

$$
G_{\ell\ell'}(\mathbf{k},i\tau) = -\sum_n^\text{empty} \psi_{n\mathbf{k}}(\mathbf{r}_\ell)\,
\psi^*_{n\mathbf{k}}(\mathbf{r}_{\ell'})\,
e^{-\tau|\varepsilon_{n\mathbf{k}}^\text{DFT}-\varepsilon_\text{F}|}\,,
$$

$$
G^\mathbf{R}_{\ell\ell'}(i\tau) = \frac{1}{N_k}\sum_\mathbf{k} e^{-i\mathbf{k}\cdot\mathbf{R}}\,
G_{\ell\ell'}(\mathbf{k},i\tau)\,.
$$

For $\tau<0$, the sum runs over the occupied bands and the sign is positive. In this representation,
the density response $\chi$ and the correlation self-energy $\Sigma^\text{c}$ are element-wise
products, without any sum over grid points or *k*-points:

$$
\chi^\mathbf{R}_{\ell\ell'}(i\tau) = -2\, G^\mathbf{R}_{\ell\ell'}(i\tau)\,
G^{-\mathbf{R}}_{\ell'\ell}(-i\tau)\,,
$$

$$
\Sigma^{\text{c},\mathbf{R}}_{\ell\ell'}(i\tau) = - G^\mathbf{R}_{\ell\ell'}(i\tau)\,
W^{\text{c},\mathbf{R}}_{\ell\ell'}(i\tau)\,.
$$

The factor 2 in $\chi$ accounts for the two spin states of a closed-shell system. Between these two
steps, $\chi^\mathbf{R}_{\ell\ell'}(i\tau)$ is transformed to the wavevectors $\mathbf{q}$ of the
*k*-point mesh and to imaginary frequency $i\omega$, and is projected onto the RI basis with the
RI-RS coefficients $Z_{\ell P}$ from Sec. 4.1. The dielectric matrix
$\varepsilon_{PQ}(\mathbf{q},i\omega)$ is inverted for every $\mathbf{q}$ and $i\omega$, which gives
the correlation part $W^\text{c}=W-V$ of the screened Coulomb interaction $W$ ($V$: bare Coulomb
interaction). $W^\text{c}$ is transformed back to the RI-RS grid, to lattice vectors and to
imaginary time:

$$
\chi_{PQ}(\mathbf{q},i\tau) = \sum_{\ell\ell'} Z^*_{\ell P}(\mathbf{q})\,
\chi_{\ell\ell'}(\mathbf{q},i\tau)\, Z_{\ell' Q}(\mathbf{q})
\,,
$$

$$
\chi_{PQ}(\mathbf{q},i\omega) \;\rightarrow\; \varepsilon_{PQ}(\mathbf{q},i\omega)
\;\rightarrow\; W^\text{c}_{PQ}(\mathbf{q},i\omega)
\;\rightarrow\; W^{\text{c},\mathbf{R}}_{\ell\ell'}(i\tau)\,.
$$

The self-energy in the lattice representation can be transformed to any $\mathbf{k}$, so the
quasiparticle energies are computed directly at the *k*-points of the band structure path, which do
not need to be part of the *k*-point mesh ($\Sigma^\text{x}$: exchange self-energy):

$$
\Sigma_{n\mathbf{k}}(i\omega) = \sum_{\ell\ell'} \psi^*_{n\mathbf{k}}(\mathbf{r}_\ell)
\,\Sigma_{\ell\ell'}(\mathbf{k},i\omega)\,\psi_{n\mathbf{k}}(\mathbf{r}_{\ell'})\,,
$$

$$
\Sigma_{\ell\ell'}(\mathbf{k},i\omega) = \sum_\mathbf{R} e^{i\mathbf{k}\cdot\mathbf{R}}
\left[\Sigma^{\text{x},\mathbf{R}}_{\ell\ell'} + \Sigma^{\text{c},\mathbf{R}}_{\ell\ell'}(i\omega)\right].
$$

Because of the element-wise products, the computation time of $\chi$ and $\Sigma$ increases only
linearly with the number of *k*-points. The $\mathbf{q}\rightarrow 0$ contribution of $W$, which a
finite *k*-point mesh does not sample, is added by a correction for 3D crystals and for 2D
materials. A publication with all details of the method is in preparation.

### 4.3 Input

*k*-point sampling is defined in the [DFT/KPOINTS](#CP2K_INPUT.FORCE_EVAL.DFT.KPOINTS) section. The
*GW* calculation uses the same *k*-point mesh as DFT, i.e. the density response $\chi$, the screened
Coulomb interaction $W$ and the self-energy $\Sigma$ are computed on this mesh:

```
&DFT
  ...
  &KPOINTS
    SCHEME MONKHORST-PACK 12 12 12   ! 3D crystal; use e.g. 16 16 1 for a 2D material
    PARALLEL_GROUP_SIZE   -1
  &END KPOINTS
&END DFT
```

The *GW* band structure calculation is activated via the
[PROPERTIES/BANDSTRUCTURE/GW](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.BANDSTRUCTURE.GW) section. A minimal
input for ev*GW*<sub>0</sub>@PBE with spin-orbit coupling reads:

```
&PROPERTIES
  &BANDSTRUCTURE
    &GW
      &AUTO_RI
        RI_AO_RATIO 3
      &END AUTO_RI
      &RI_RS
        &GRID_OPTIMIZATION
          RS_AO_RATIO 7
        &END GRID_OPTIMIZATION
      &END RI_RS
      &EVGW0
      &END EVGW0
    &END GW
    &SOC
    &END SOC
    &BANDSTRUCTURE_PATH
      NPOINTS 20
      SPECIAL_POINT GAMMA 0.000 0.000 0.000
      SPECIAL_POINT X     0.000 0.500 0.500
      SPECIAL_POINT L     0.500 0.500 0.500
    &END BANDSTRUCTURE_PATH
  &END BANDSTRUCTURE
&END PROPERTIES
```

All other parameters have well-tested default values. The sections and keywords from above are:

- [&AUTO_RI](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.BANDSTRUCTURE.GW.AUTO_RI): Automatic RI basis, see
  Sec. 4.1. No `RI_AUX` basis set needs to be specified in `&KIND`.
- [&RI_RS](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.BANDSTRUCTURE.GW.RI_RS) with
  [&GRID_OPTIMIZATION](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.BANDSTRUCTURE.GW.RI_RS.GRID_OPTIMIZATION):
  Activates RI-RS with the automatic grid, see Sec. 4.1.
- [&EVGW0](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.BANDSTRUCTURE.GW.EVGW0): Activates ev*GW*<sub>0</sub>.
  Without this section, CP2K performs *G*<sub>0</sub>*W*<sub>0</sub>.
- [&SOC](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.BANDSTRUCTURE.SOC): Activates spin-orbit coupling (SOC
  from Hartwigsen-Goedecker-Hutter pseudopotentials \[[](#Hartwigsen1998)\]). SOC also needs
  `POTENTIAL_FILE_NAME  GTH_SOC_POTENTIALS`.
- [&BANDSTRUCTURE_PATH](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.BANDSTRUCTURE.BANDSTRUCTURE_PATH): Specify
  the *k*-path in the Brillouin zone for computing the band structure. Relative *k*-coordinates are
  needed which you can retrieve for your crystal structure from \[[](#Setyawan2010)\].

The *k*-point mesh is a convergence parameter for both DFT and *GW*. For the benchmark below, we
used 12×12×12 for 3D crystals and 16×16×1 for 2D materials.

### 4.4 Basis sets

We recommend the augmented basis sets `aug-SZV-MOLOPT-GTH-tier-2`, `aug-DZVP-MOLOPT-GTH-tier-2` and
`aug-TZVP-MOLOPT-GTH-tier-2` from the file `BASIS_AUG_MOLOPT_GTH` together with GTH
pseudopotentials. `aug-SZV-MOLOPT-GTH-tier-2` is a cheap choice for a first band structure,
`aug-TZVP-MOLOPT-GTH-tier-2` is the accurate choice.

**Why augmented basis sets?** The MOLOPT basis sets are optimized for ground-state total energies,
which depend only on the occupied orbitals. The occupied orbitals and their energies are therefore
well converged, but the empty orbitals are not: no basis function was optimized for them. By the
variational principle, an orbital that is poorly described by the basis set has a too high energy.
With a MOLOPT basis set, the conduction bands are thus too high in energy and the band gap is too
large. In *GW*, the empty orbitals enter twice, as conduction bands of the band structure and in the
sums over empty bands in $G$ and $\chi$ (Sec. 4.2), so that *GW* band gaps converge slowly with the
size of a ground-state basis set. Additional diffuse basis functions describe the empty orbitals,
lower the conduction bands towards their basis-set limit and make the *GW* band gap converge much
faster with the basis-set size.

**Construction.** An augmented basis set contains the unchanged MOLOPT basis set (`SZV-MOLOPT-GTH`,
`DZVP-MOLOPT-GTH` or `TZVP-MOLOPT-GTH`) plus one additional contracted Gaussian function per angular
momentum $l$:

- `tier-1`: one additional function for every $l$ from 0 to the highest angular momentum
  $l_\text{max}$ of the MOLOPT basis set. Example carbon: `SZV-MOLOPT-GTH` has one *s* and one *p*
  function, `aug-SZV-MOLOPT-GTH-tier-1` has two *s* and two *p* functions.
- `tier-2`: as `tier-1`, and one additional function with $l_\text{max}+1$. Example carbon:
  `aug-SZV-MOLOPT-GTH-tier-2` has two *s*, two *p* and one *d* function.

The additional functions were optimized on a set of molecules such that the five lowest
*GW*+Bethe-Salpeter excitation energies agree with the complete-basis-set limit. At the same time,
the condition number of the overlap matrix was kept small, so that the SCF of solids and large
molecules converges without problems. The total energy is of similar quality as with the underlying
MOLOPT basis set. The same basis sets are suited for other excited-state methods (TDDFT,
Bethe-Salpeter). With *k*-point sampling, the RI basis is generated automatically (`&AUTO_RI`, Sec.
4.1), no RI basis set file is needed.

### 4.5 Analytic continuation

CP2K computes the *GW* self-energy on the imaginary frequency axis, $\Sigma_{n\mathbf{k}}(i\omega)$.
The quasiparticle equation from Sec. 2, however, needs the self-energy at real frequencies,
$\Sigma_{n\mathbf{k}}(\omega)$. CP2K obtains it by analytic continuation: a Padé approximant is
fitted to the values $\Sigma_{n\mathbf{k}}(i\omega_j)$ on the imaginary frequency grid, and this
fitted function is then evaluated at real frequencies,

$$ \Sigma_{n\mathbf{k}}(i\omega_j) \;\xrightarrow{\;\text{fit}\;}\; \Sigma_{n\mathbf{k}}^\text{Padé}(z)
\;\xrightarrow{\;z\,=\,\omega\;}\; \Sigma_{n\mathbf{k}}(\omega)\,. $$

By default, CP2K uses a Nevanlinna function for this fit
([ANALYTIC_CONTINUATION](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.BANDSTRUCTURE.GW.ANALYTIC_CONTINUATION)
`NEVANLINNA`, \[[](#Fei2021)\]), named after the Finnish mathematician Rolf Nevanlinna, i.e. only
functions with a positive spectral function on the real-frequency axis, as for the exact
self-energy. The quasiparticle energy is the median over several such fits to subsets of the
imaginary frequencies, and their spread serves as an error estimate: in
`bandstructure_SCF_and_G0W0`, a level is marked `*` if its analytic continuation is unstable (the
fits disagree, or a level between it and the HOMO or LUMO is marked), and `m` if the quasiparticle
equation has several solutions of comparable weight.

```{warning}
With analytic continuation, *GW* quasiparticle energies are only reliable in the energy window from
about 3 eV below the valence band maximum to about 3 eV above the conduction band minimum,
$\varepsilon_\text{VBM}-3\,\text{eV} \lesssim \varepsilon_{n\mathbf{k}}^{GW} \lesssim
\varepsilon_\text{CBM}+3\,\text{eV}$. Deeper occupied levels and higher unoccupied levels can be very
unreliable. Increasing the number of imaginary-frequency points does not fix this, because the
continuation from the imaginary to the real frequency axis is a numerically ill-posed problem. We are
working on more robust methods; until then, please keep this limitation in mind when you interpret
*GW* band structures. With `NEVANLINNA`, levels outside the reliable window are marked `*`; the marks
need enough imaginary-frequency points to be meaningful, we recommend NUM_TIME_FREQ_POINTS ≥ 26.
```

### 4.6 Output, accuracy and computation time

The *GW* band structure is written to the files `bandstructure_SCF_and_G0W0` and
`bandstructure_SCF_and_G0W0_plus_SOC`. The direct and indirect band gaps are also listed in the CP2K
output file.

Inputs, outputs and *GW* band structures of 2D materials and 3D crystals with five different basis
sets are available in the GitHub repository
[CP2K_GW_RIRS_Bandstructures_3D_2D](https://github.com/JWilhelm/CP2K_GW_RIRS_Bandstructures_3D_2D).
We recommend to start from one of these inputs. The table summarizes the mean computational cost and
the mean absolute error (MAE) of the *G*<sub>0</sub>*W*<sub>0</sub>@PBE band gap:

| Materials                                | Basis set                    | Cost (core-h) | Gap MAE to ref. (eV) |
| ---------------------------------------- | ---------------------------- | ------------: | -------------------: |
| 2D materials<br>(MoS₂, MoSe₂, WS₂, WSe₂) | `aug-SZV-MOLOPT-GTH-tier-2`  |           2.2 |                0.074 |
| 2D materials<br>(MoS₂, MoSe₂, WS₂, WSe₂) | `aug-TZVP-MOLOPT-GTH-tier-2` |          22.8 |                0.040 |
| 3D crystals                              | `aug-SZV-MOLOPT-GTH-tier-2`  |          42.7 |                0.205 |
| 3D crystals                              | `aug-TZVP-MOLOPT-GTH-tier-2` |         973.4 |                0.137 |

The cost is the mean over all materials. As an example, 22.8 core-hours correspond to about 11
minutes on a single compute node with 128 cores. For 2D materials, the *GW* calculation is very
cheap. For 3D crystals, the *GW* calculation is only roughly 2–10 times more expensive than the
preceding PBE calculation.

### 4.7 Preliminary benchmark: band gaps compared to VASP

The table lists preliminary *G*<sub>0</sub>*W*<sub>0</sub>@PBE band gaps (in eV) of 29 3D crystals
computed with CP2K (indirect gap, without spin-orbit coupling, *k*-point mesh 12×12×12, all *GW*
parameters at their default values) together with plane-wave *G*<sub>0</sub>*W*<sub>0</sub>@PBE band
gaps from VASP ([Klimeš 2014](https://doi.org/10.1103/PhysRevB.90.075125) and
[Varrassi 2025](https://doi.org/10.1038/s41524-025-01833-w)). The CP2K columns are labeled by the
basis set, where aug-SZV tier-2 stands for `aug-SZV-MOLOPT-GTH-tier-2` etc. The last row is the mean
absolute deviation (MAD) of CP2K from VASP over the 28 crystals with a VASP value; where two VASP
values are available, their mean is the reference.

| Material         | CP2K<br>aug-SZV<br>tier-2 | CP2K<br>aug-DZVP<br>tier-2 | CP2K<br>aug-TZVP<br>tier-2 | VASP<br>[Klimeš 2014](https://doi.org/10.1103/PhysRevB.90.075125) | VASP<br>[Varrassi 2025](https://doi.org/10.1038/s41524-025-01833-w) |
| ---------------- | ------------------------: | -------------------------: | -------------------------: | ----------------------------------------------------------------: | ------------------------------------------------------------------: |
| AlAs (As: q23)   |                      2.01 |                       2.07 |                       2.16 |                                                              2.31 |                                                                2.36 |
| AlAs (As: q5)    |                      2.00 |                       1.97 |                       2.10 |                                                              2.31 |                                                                2.36 |
| AlN              |                      6.50 |                       5.80 |                       5.86 |                                                                 – |                                                                6.07 |
| AlP              |                      2.28 |                       2.36 |                       2.52 |                                                              2.48 |                                                                2.46 |
| AlSb (Sb: q23)   |                      1.75 |                       1.65 |                       1.77 |                                                              1.87 |                                                                1.88 |
| AlSb (Sb: q5)    |                      1.66 |                       1.56 |                       1.68 |                                                              1.87 |                                                                1.88 |
| BN (cubic)       |                      6.28 |                       6.29 |                       6.32 |                                                              6.40 |                                                                6.39 |
| BN (hexagonal)   |                      5.95 |                       5.94 |                       6.00 |                                                                 – |                                                                6.17 |
| BP               |                      2.17 |                       2.02 |                       2.09 |                                                                 – |                                                                2.13 |
| C                |                      5.58 |                       5.61 |                       5.59 |                                                                 – |                                                                5.72 |
| CdS (wurtzite)   |                      2.01 |                       2.08 |                       2.17 |                                                                 – |                                                                2.19 |
| CdS (zincblende) |                      1.93 |                       1.99 |                       2.08 |                                                              2.15 |                                                                   – |
| CdSe             |                      1.14 |                       1.23 |                       1.33 |                                                              1.52 |                                                                1.71 |
| CdTe             |                      1.49 |                       1.35 |                       1.41 |                                                              1.57 |                                                                1.68 |
| GaAs             |                      1.15 |                       1.06 |                       1.13 |                                                              1.23 |                                                                1.34 |
| GaP              |                      2.18 |                       2.18 |                       2.34 |                                                              2.30 |                                                                2.37 |
| Ge               |                      0.49 |                       0.56 |                       0.59 |                                                              0.67 |                                                                   – |
| InP              |                      1.17 |                       1.07 |                       1.20 |                                                              1.23 |                                                                1.29 |
| LiCl             |                      8.49 |                       8.93 |                       9.17 |                                                                 – |                                                                9.21 |
| LiF              |                     14.30 |                      13.72 |                      13.72 |                                                                 – |                                                               13.74 |
| MgO              |                      7.66 |                       7.44 |                       7.39 |                                                              7.55 |                                                                7.51 |
| NaCl             |                      7.66 |                       7.96 |                       8.28 |                                                                 – |                                                                9.27 |
| Si               |                      1.11 |                       1.08 |                       1.17 |                                                                 – |                                                                1.30 |
| SiC              |                      2.35 |                       2.28 |                       2.35 |                                                              2.42 |                                                                2.42 |
| TiO₂ (anatase)   |                      3.10 |                       3.65 |                       3.74 |                                                                 – |                                                                   – |
| TiO₂ (rutile)    |                      2.80 |                       3.32 |                       3.41 |                                                                 – |                                                                3.39 |
| ZnS              |                      3.32 |                       3.30 |                       3.37 |                                                              3.46 |                                                                3.48 |
| ZnSe             |                      2.22 |                       2.20 |                       2.29 |                                                              2.55 |                                                                2.65 |
| ZnTe             |                      2.46 |                       2.05 |                       2.10 |                                                              2.27 |                                                                2.33 |
| **MAD to VASP**  |                  **0.30** |                   **0.24** |                   **0.15** |                                                                   |                                                                     |

For AlAs and AlSb, two calculations with different GTH pseudopotentials are listed: As and Sb with 5
valence electrons (q5) and with 23 valence electrons (q23). Zn and Cd have 20 valence electrons, Ga
and In have 21. For NaCl, CP2K is 1.0 eV below the VASP value; all-electron
*G*<sub>0</sub>*W*<sub>0</sub>@PBE calculations with FHI-aims give 7.80 eV and 7.92 eV
([Ren *et al.*, Phys. Rev. Materials **5**, 013807 (2021)](https://doi.org/10.1103/PhysRevMaterials.5.013807)).
For TiO₂ anatase, no *G*<sub>0</sub>*W*<sub>0</sub>@PBE band gap is available in the two VASP
references.

ZnO is not part of the table. The *G*<sub>0</sub>*W*<sub>0</sub>@PBE band gap of ZnO has been
debated for a long time, because it converges very slowly with the number of empty bands
([Shih *et al.* 2010](https://doi.org/10.1103/PhysRevLett.105.146401),
[Friedrich *et al.* 2011](https://doi.org/10.1103/PhysRevB.83.081101)). With the present basis sets,
CP2K cannot yet compute a converged *GW* band gap of ZnO.

For the four 2D materials MoS₂, MoSe₂, WS₂ and WSe₂, the table lists the
*G*<sub>0</sub>*W*<sub>0</sub>@PBE band gap at K with spin-orbit coupling (in eV, *k*-point mesh
16×16×1) together with VASP ([arXiv:2507.18411](https://arxiv.org/abs/2507.18411)):

| Material        | CP2K<br>aug-SZV<br>tier-2 | CP2K<br>aug-DZVP<br>tier-2 | CP2K<br>aug-TZVP<br>tier-2 | VASP |
| --------------- | ------------------------: | -------------------------: | -------------------------: | ---: |
| MoS₂            |                      2.36 |                       2.34 |                       2.33 | 2.29 |
| MoSe₂           |                      1.92 |                       1.98 |                       1.97 | 2.01 |
| WS₂             |                      2.41 |                       2.37 |                       2.39 | 2.37 |
| WSe₂            |                      1.93 |                       1.98 |                       1.97 | 2.02 |
| **MAD to VASP** |                  **0.07** |                   **0.03** |                   **0.03** |      |

These numbers are preliminary; a publication is in preparation. All raw data (inputs, outputs, band
structures and compute times, also for the basis sets `aug-DZVP-MOLOPT-GTH-tier-1` and
`aug-TZVP-MOLOPT-GTH-tier-1`) are available in the GitHub repository
[CP2K_GW_RIRS_Bandstructures_3D_2D](https://github.com/JWilhelm/CP2K_GW_RIRS_Bandstructures_3D_2D).

In case anything does not work, please feel free to contact jan.wilhelm (at) ur.de.
