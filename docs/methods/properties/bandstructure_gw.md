# Electronic band structure from *GW*

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
$\mathbf{k}$. For a molecule, the KS orbitals $\psi_n(\mathbf{r})$ and KS eigenvalues
$\varepsilon_n^\text{DFT}$ only carry a single quantum number, the molecular-orbital index $n$.

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
v_\text{xc}|\psi_{n\mathbf{k}}}$, from the DFT eigenvalue
$\varepsilon_{n\mathbf{k}}^\text{DFT}$ and we add the xc contribution from
*G*<sub>0</sub>*W*<sub>0</sub>,
$\braket{\psi_{n\mathbf{k}}|
\Sigma^{G_0W_0}(\varepsilon_{n\mathbf{k}}^{G_0W_0}) |\psi_{n\mathbf{k}}}$.

The DFT xc start functional of *G*<sub>0</sub>*W*<sub>0</sub> can influence the
*G*<sub>0</sub>*W*<sub>0</sub> quasiparticle energies. For molecules, it is recommended to start the
*G*<sub>0</sub>*W*<sub>0</sub> calculation from the PBE0 functional and for solids, from the PBE
functional.

CP2K also allows to perform eigenvalue-selfconsistency in $G$ (ev*GW*<sub>0</sub>) and
eigenvalue-selfconsistency in $G$ and in $W$ (ev*GW*). For solids, it has been shown that band gaps
from ev*GW*$_0$@PBE can be in better agreement with experimental band gaps than band gaps from
*G*<sub>0</sub>*W*<sub>0</sub>@PBE.

CP2K contains three different *GW* implementations:

- *GW* for molecules \[[](#Wilhelm2016)\] (Sec. 3)
- *GW* for computing the band structure of a solid with small unit cell with *k*-point sampling in
  DFT (publication in preparation TODO: insert reference, Sec. 4)
- *GW* for computing the band structure of a large cell in a Γ-only approach \[[](#Graml2024)\]
  (Sec. 5)

In the following, we will discuss details and usage of these *GW* implementations.

## 3. *GW* for molecules

For starting a *G*<sub>0</sub>*W*<sub>0</sub>, ev*GW*<sub>0</sub> or ev*GW* calculation for a
molecule, one needs to set the [RUN_TYPE](#CP2K_INPUT.GLOBAL.RUN_TYPE) to `ENERGY` and one needs to
put the following section:

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
  ev*GW*) is used to calculate the *GW* quasiparticle energies.

The numerical precision of the *GW* implementation is 10 meV compared to reference calculations, for
example the *GW*100 test set \[[](#vanSetten2015)\]. To help new users get familiar with the *GW*
implementation in CP2K, we recommend starting by reproducing the HOMO and LUMO
*G*<sub>0</sub>*W*<sub>0</sub>@PBE quasiparticle energies of the H<sub>2</sub>O molecule from the
GW100 test set. The reference values are $\varepsilon_\text{HOMO}^{G_0W_0\text{@PBE}}$ = -11.97 eV
and $\varepsilon_\text{LUMO}^{G_0W_0\text{@PBE}}$ = 2.37 eV; CP2K input and output files for the
*G*<sub>0</sub>*W*<sub>0</sub>@PBE calculation of the H<sub>2</sub>O molecule are available
[here](https://github.com/cp2k/cp2k-examples/tree/master/bandstructure_gw/1_H2O_GW100).

The following settings from DFT will also have an influence on *GW* quasiparticle energies:

- [XC_FUNCTIONAL](#CP2K_INPUT.FORCE_EVAL.DFT.XC.XC_FUNCTIONAL): Starting xc functional for the
  *G*<sub>0</sub>*W*<sub>0</sub>, ev*GW*<sub>0</sub> or ev*GW* calculation. For molecules, we
  recommend either ev*GW*<sub>0</sub>@PBE or *G*<sub>0</sub>*W*<sub>0</sub>@PBE0. For further
  guidance on selecting an appropriate DFT starting functional and self-consistency scheme for your
  system, you may consult \[[](#Golze2019)\].

- [BASIS_SET](#CP2K_INPUT.FORCE_EVAL.SUBSYS.KIND.BASIS_SET): The basis set is of Gaussian type and
  can have strong influence on the quasiparticle energies. For computing quasiparticle energies, a
  basis set extrapolation is necessary, see Fig. 2a in \[[](#Wilhelm2016)\] and we recommend
  all-electron GAPW calculations with correlation-consistent basis sets from the
  <a href="https://www.basissetexchange.org/" target="_blank">EMSL database</a>. For computing the
  HOMO-LUMO gap from *GW*, we recommend augmented basis sets, for example
  <a href="https://www.basissetexchange.org/basis/aug-cc-pvdz/format/cp2k/?version=1&elements=1,6,8" target="_blank">`aug-cc-pVDZ`</a>
  and
  <a href="https://www.basissetexchange.org/basis/aug-cc-pvtz/format/cp2k/?version=1&elements=1,6,8" target="_blank">`aug-cc-pVTZ`</a>.
  As `RI_AUX` basis set, we recommend the RIFIT basis sets from the EMSL database, for example
  <a href="https://www.basissetexchange.org/basis/aug-cc-pvdz-rifit/format/cp2k/?version=1&elements=1,6,8" target="_blank">`aug-cc-pVDZ-RIFIT`</a>.

The computational effort of the *GW* calculation increases with *N*<sup>4</sup> in system size *N*.
The memory requirement increases with *N*<sup>3</sup>. For running large-scale calculations, we
recommend starting with a small molecule. After successfully completing the *GW* calculation for the
small molecule, you can gradually increase the molecule size. The computational resources needed for
larger molecules can then be estimated using the *N*<sup>4</sup> scaling for computation time and
*N*<sup>3</sup> scaling for memory. The output provides a useful lower limit of the required memory
is given: (TODO: will be replaced by large-scale calculation)

```
  RI_INFO| Total memory for (ia|K) integrals:                           1.55 MiB
  RI_INFO| Total memory for G0W0-(nm|K) integrals:                      8.08 MiB
```

When facing out-of-memory, please increase the number of nodes of your calculation.

Input and output of a large-scale *GW* calculation on a nanographene with 200 atoms is available
<a href="https://github.com/cp2k/cp2k-examples/tree/master/bandstructure_gw/2_200_atom_molecule_nanographene" target="_blank">here</a>.

## 4. *GW* for small unit cells with *k*-point sampling

For a periodic *GW* calculation, *k*-point sampling is required. *k*-point sampling is included in
the DFT section, see below a *k*-point section for a 2D-periodic cell:

```
&DFT
  ...
  &KPOINTS
    SCHEME MONKHORST-PACK 32 32 1
    PARALLEL_GROUP_SIZE   -1
  &END KPOINTS
&END DFT
```

The *k*-point mesh size is a convergence parameter, 32x32 for a 2D material is expected to reach
convergence of the *GW* bandgap within 10 meV.

The periodic *GW* calculation is activated via the bandstructure section:

```
 &PROPERTIES
    &BANDSTRUCTURE
      &GW
        NUM_TIME_FREQ_POINTS         30
        MEMORY_PER_PROC              10
        EPS_FILTER              1.0E-11
        REGULARIZATION_RI        1.0E-2
        CUTOFF_RADIUS_RI            7.0
      &END GW
      &SOC
      &END SOC
      &BANDSTRUCTURE_PATH
        NPOINTS 19
        SPECIAL_POINT K     0.33 0.33 0.00
        SPECIAL_POINT GAMMA 0.00 0.00 0.00
        SPECIAL_POINT M     0.00 0.50 0.00
      &END
    &END BANDSTRUCTURE
  &END PROPERTIES
```

All parameters from above have been chosen to converge the *GW* bandgap within 10 meV, see also
details in forthcoming publication (TODO Link):

- [NUM_TIME_FREQ_POINTS](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.BANDSTRUCTURE.GW.NUM_TIME_FREQ_POINTS):
  Number of imaginary-time and imaginary-frequency points used for computing the self-energy.
  Between 20 and 30 points are usually enough for converging quasiparticle energies within 10 meV.
  Grids up to 34 points are available.
- [MEMORY_PER_PROC](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.BANDSTRUCTURE.GW.MEMORY_PER_PROC): Specifies
  the available memory per MPI process in GB which is specific to your hardware. A larger
  `MEMORY_PER_PROC` can increase performance but also the memory requirement increases.
- [EPS_FILTER](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.BANDSTRUCTURE.GW.EPS_FILTER): Filter for
  three-center integrals, 10<sup>-11</sup> should be well-converged.
- [REGULARIZATION_RI](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.BANDSTRUCTURE.GW.REGULARIZATION_RI):
  Regularization parameter for resolution-of-the-identity (RI) basis set. For big RI basis set (> 50
  RI function per atom) we recommend 10<sup>-2</sup> to prevent for linear dependencies. For small
  RI basis set, one can turn RI regularization off by setting 0.0.
- [CUTOFF_RADIUS_RI](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.BANDSTRUCTURE.GW.CUTOFF_RADIUS_RI): Cutoff
  radius of truncated Coulomb metric in Å. A larger cutoff leads to faster the RI basis set
  convergence, but also computational cost increases. A cutoff of 7 Å is an accurate choice.
- [&SOC](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.BANDSTRUCTURE.SOC): Activates spin-orbit coupling (SOC
  from Hartwigsen-Goedecker-Hutter pseudopotentials \[[Hartwigsen1998](#Hartwigsen1998)\]). SOC also
  needs `POTENTIAL_FILE_NAME  GTH_SOC_POTENTIALS`.
- [&BANDSTRUCTURE_PATH](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.BANDSTRUCTURE.BANDSTRUCTURE_PATH): Specify
  the *k*-path in the Brillouin zone for computing the band structure. Relative *k*-coordinates are
  needed which you can retrieve for your crystal structure from \[[](#Setyawan2010)\].

We recommend TZVP-MOLOPT basis sets together with GTH/HGH pseudopotentials, see basis set
convergence study in (TODOref). At present, 2d-periodic boundary conditions are supported, 1d- and
3d-periodic boundary conditions are work in progress.

The *GW* band structure is written to the files `bandstructure_SCF_and_G0W0` and
`bandstructure_SCF_and_G0W0_plus_SOC`. The direct and indirect band gaps are also listed in the CP2K
output file. When facing an out-of-memory crash, please increase `MEMORY_PER_PROC`. An input and
output for a *G*<sub>0</sub>*W*<sub>0</sub>@PBE band structure calculation of the 2d material
WSe<sub>2</sub> can be found
\[<a href="https://github.com/cp2k/cp2k-examples/tree/master/bandstructure_gw/3_3-atom_unit_cell_2d_WSe2_loose_parameters" target="_blank">here</a>\]
using loose parameters (*G*<sub>0</sub>*W*<sub>0</sub>@PBE band gap: 2.30 eV, computation time: 3
hours on 3 large-memory nodes) and
\[<a href="https://github.com/cp2k/cp2k-examples/tree/master/bandstructure_gw/4_3-atom_unit_cell_2d_WSe2_tight_parameters" target="_blank">here</a>\]
using tight parameters (*G*<sub>0</sub>*W*<sub>0</sub>@PBE band gap: 2.30 eV, computation time: 12
hours on 20 large-memory nodes).

## 5. *GW* for large cells in Γ-only approach

For a large unit cell, a Γ-only *GW* algorithm is available in CP2K \[[](#Graml2024)\]. The
requirement on the cell is that elements of the density matrix decay by several orders of magnitude
when the two basis functions of the matrix element have a distance similar to the cell size. As rule
of thumb, for a 2d material, a 9x9 unit cell is large enough for the Γ-only algorithm, see Fig. 1 in
\[[](#Graml2024)\].

The input file for a Γ-only *GW* calculation is identical as for *GW* for small cells with *k*-point
sampling except that the `&KPOINTS` section in DFT needs to be removed. An exemplary input and
output is available
\[<a href="https://github.com/cp2k/cp2k-examples/tree/master/bandstructure_gw/5_9x9_supercell_2d_MoS2" target="_blank">here</a>\].
Running the input file requires access to a large computer (the calculation took 2.5 hours on 32
nodes on Noctua2 cluster in Paderborn). The computational parameters from this input file reach
numerical convergence of the band gap within ~ 50 meV (TZVP basis set, 10 time and frequency
points). Detailed convergence tests are available in the SI, Table S1 of \[[](#Graml2024)\] We
recommend the numerical parameters from the input file for large-scale *GW* calculations. The code
prints restart files with ending .matrix that can be used to restart a crashed calculation.

In case anything does not work, please feel free to contact jan.wilhelm (at) ur.de.
