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

- *GW* for molecules \[[](#Wilhelm2016a)\] (Sec. 3)
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
  \[[](#Wilhelm2016a)\]). 100 points are usually enough for converging quasiparticle energies within
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
  basis set extrapolation is necessary, see Fig. 2a in \[[](#Wilhelm2016a)\] and we recommend
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

### 4.2 Input

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
used 12×12×12 for 3D crystals and 16×16×1 for 2D materials. We recommend the augmented basis sets
`aug-SZV-MOLOPT-GTH-tier-2`, `aug-DZVP-MOLOPT-GTH-tier-2` and `aug-TZVP-MOLOPT-GTH-tier-2` from the
file `BASIS_AUG_MOLOPT_GTH` together with GTH pseudopotentials. `aug-SZV-MOLOPT-GTH-tier-2` is a
cheap choice for a first band structure, `aug-TZVP-MOLOPT-GTH-tier-2` is the accurate choice.

### 4.3 Analytic continuation

CP2K computes the *GW* self-energy on the imaginary frequency axis, $\Sigma_{n\mathbf{k}}(i\omega)$.
The quasiparticle equation from Sec. 2, however, needs the self-energy at real frequencies,
$\Sigma_{n\mathbf{k}}(\omega)$. CP2K obtains it by analytic continuation: a Padé approximant is
fitted to the values $\Sigma_{n\mathbf{k}}(i\omega_j)$ on the imaginary frequency grid, and this
fitted function is then evaluated at real frequencies,

$$ \Sigma_{n\mathbf{k}}(i\omega_j) \;\xrightarrow{\;\text{fit}\;}\; \Sigma_{n\mathbf{k}}^\text{Padé}(z)
\;\xrightarrow{\;z\,=\,\omega\;}\; \Sigma_{n\mathbf{k}}(\omega)\,. $$

```{warning}
With analytic continuation, *GW* quasiparticle energies are only reliable in the energy window from
about 3 eV below the valence band maximum to about 3 eV above the conduction band minimum,
$\varepsilon_\text{VBM}-3\,\text{eV} \lesssim \varepsilon_{n\mathbf{k}}^{GW} \lesssim
\varepsilon_\text{CBM}+3\,\text{eV}$. Deeper occupied levels and higher unoccupied levels can be very
unreliable. Increasing the number of imaginary-frequency points does not fix this, because the
continuation from the imaginary to the real frequency axis is a numerically ill-posed problem. We are
working on more robust methods; until then, please keep this limitation in mind when you interpret
*GW* band structures.
```

### 4.4 Output, accuracy and computation time

The *GW* band structure is written to the files `bandstructure_SCF_and_G0W0` and
`bandstructure_SCF_and_G0W0_plus_SOC`. The direct and indirect band gaps are also listed in the CP2K
output file.

Inputs, outputs and *GW* band structures of 2D materials and 3D crystals with five different basis
sets are available in the GitHub repository
[CP2K_GW_RIRS_Bandstructures_3D_2D](https://github.com/JWilhelm/CP2K_GW_RIRS_Bandstructures_3D_2D).
We recommend to start from one of these inputs. The table summarizes the mean computational cost and
the mean absolute error (MAE) of the *G*<sub>0</sub>*W*<sub>0</sub>@PBE band gap:

| Materials                             | Basis set                    | Cost (core-h) | Gap MAE to ref. (eV) | Gap MAE to prev. CP2K (eV) |
| ------------------------------------- | ---------------------------- | ------------: | -------------------: | -------------------------: |
| 2D materials (MoS₂, MoSe₂, WS₂, WSe₂) | `aug-SZV-MOLOPT-GTH-tier-2`  |           2.2 |                0.074 |                      0.040 |
| 2D materials (MoS₂, MoSe₂, WS₂, WSe₂) | `aug-TZVP-MOLOPT-GTH-tier-2` |          22.8 |                0.040 |                      0.034 |
| 3D crystals                           | `aug-SZV-MOLOPT-GTH-tier-2`  |          42.7 |                0.205 |                          – |
| 3D crystals                           | `aug-TZVP-MOLOPT-GTH-tier-2` |         973.4 |                0.137 |                          – |

The cost is the mean over all materials. As an example, 22.8 core-hours correspond to about 11
minutes on a single compute node with 128 cores. For 2D materials, the *GW* calculation is very
cheap. For 3D crystals, the *GW* calculation is only roughly 2–10 times more expensive than the
preceding PBE calculation.

In case anything does not work, please feel free to contact jan.wilhelm (at) ur.de.
