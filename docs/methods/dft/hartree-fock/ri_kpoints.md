# HFX-RI with k-Points

The resolution-of-the-identity (RI) technique is implemented for various methods in CP2K, each time
in a slightly different flavor. Hartree-Fock exchange (HFX) has two distinct RI implementations: One
for [$\Gamma$-point calculations](./ri_gamma), and one for k-point sampling, which is covered on
this page.

The implementation for RI-HFXk (with k-point sampling) is described in
[doi:10.1063/5.0189659](https://doi.org/10.1063/5.0189659). All files necessary to run the examples
discussed below can be found [here](https://www.cp2k.org/_media/howto:ri_hfx_examples.zip).

## Application domain

The RI-HFX implementation with k-point sampling (RI-HFXk) is optimized for the simulation of small
unit cells and dense k-point meshes. It is much more efficient than equivalent $\Gamma$-point
supercell calculations.

## Brief theory recap

In HFX calculations with k-point sampling, the exact-exchange contribution to the KS matrix at a
given k-point is expressed as

$$
K_{\mu\nu}^\mathbf{k} = \sum_\mathbf{k'} P_{\sigma\lambda}^{\mathbf{k'}} (\mu^\mathbf{k}\sigma^{\mathbf{k'}} | \nu^\mathbf{k}\lambda^{\mathbf{k'}})
$$

In CP2K, most k-space matrices are obtained via a Fourier transform of their real-space counterpart.
This is also the case for the exact-exchange matrix:

$$
K_{\mu\nu}^{\mathbf{k}} = \sum_\mathbf{R} e^{i\mathbf{k}\cdot\mathbf{R}}\ K^\mathbf{R}_{\mu\nu}
$$

where $\mathbf{R}$ is the translation vector for a periodic image of the simulation cell. The
real-space matrix is expressed as:

$$
K^\mathbf{b}_{\mu\nu} = \sum_{\mathbf{a},\mathbf{c}}\ P_{\sigma,\lambda}^{\mathbf{c}}\ (\mu^\mathbf{0}\sigma^\mathbf{a}| \nu^\mathbf{b}\lambda^\mathbf{a+c})
$$

where $\mathbf{a, b, c}$ are periodic image indices. Note that in $\Gamma$-point calculations, there
is a single k-point with $\mathbf{k} = 0$, leading to the same density matrix for all $\mathbf{c}$.
It can then be taken out of the sum, and the $\Gamma$-point formula from far above is recovered.

In the RI-HFXk method, the real-space exact-exchange matrices are calculated with a local
atom-specific RI basis:

$$
\begin{aligned}
    K_{\mu i,\nu j}^\textbf{b} = \sum_{\mathbf{a}, \mathbf{c}} \ &P^\mathbf{c}_{\sigma,\lambda} \ (\mu^\mathbf{0}_i\sigma^\mathbf{a}\lfloor P^\mathbf{0}_i)\ (P^\mathbf{0}_i\lfloor R^\mathbf{0}_i)^{-1}\ (R^\mathbf{0}_i | S^\mathbf{b}_j)\\
    &(S^\mathbf{b}_j\lfloor Q^\mathbf{b}_j)^{-1}\ (Q^\mathbf{b}_j\lfloor \nu^\mathbf{b}_j \lambda^{\mathbf{a}+\mathbf{c}})
\end{aligned}
$$

where indices $i,j$ refer to atoms, and $K^\mathbf{b}_{\mu i,\nu j}$ correspond to the $\mu,\nu$ AO
pair in the $i,j$ atomic block of the matrix for periodic image $\mathbf{b}$. The local RI basis
{$P^\mathbf{0}_i$} for atom $i$ in the reference cell is composed of RI basis elements coming from
all atoms within a sphere a radius $R_\text{max}$ centered on atom $i$. $R_\text{max}$ is the extent
of the most diffuse AO in the system. You can find more details in the paper accepted, not yet
published paper.

Because building the real-space exchange matrices is the most expensive part of the calculation, the
method has a constant cost regardless of the k-point mesh. When the number of k-point becomes high
(~1000), that scaling becomes linear because a matrix has to be diagonalized at each k-point. The
ADMM method can be seamlessly used with RI-HFXk.

## Simple example (graphene band structure)

In this example, the band structure of graphene at the ADMM-PBE0 level of theory is calculated. We
use the `pob-TZVP-rev2` basis set because it is not diffuse, and therefore efficient. Note that more
diffuse basis sets, such as the `ccGRB` family distributed with CP2K, seem to produce more robust
results. As the ADMM auxiliary basis set, we use the lower quality `pob-DZVP-rev2`.

To calculate a band structure, we first converge a SCF cycle with a dense, general k-point mesh. In
this case, we use a 19x19x1 Monkhorst-Pack mesh. Note that the special point K where the Dirac cone
is located is not explicitly in the mesh. This allows for a simpler calculation since the system is
treated as a semi-conductor. The mesh is however dense enough that the physics is well captured.
After the SCF is converged, the collection of real-space KS matrices $F^\textbf{b}_{\mu\nu}$ is back
Fourier transformed to $F^\textbf{k}_{\mu\nu}$ along the k-point path of interest specified in the
input. The matrices are diagonalized, and their eigenvalues used as band energies. You can used the
[cp2k_bs2csv](https://github.com/cp2k/cp2k-output-tools) to transform the resulting CP2K output to a
more workable CSV file.

Note that in this input file, we select a value of $1.0\times 10^{-6}$ for
[EPS_PGF_ORB](#CP2K_INPUT.ATOM.METHOD.XC.HF.RI.EPS_PGF_ORB). This parameter controls the range of
AOs, and threfore the extent of the local atom-specific RI basis sets. This value leads to
particularly high accuracy. The default of $1.0\times 10^{-5}$ is typically enough.

The HFX potential was selected as the Truncated Coulomb operator with a cutoff radius of $R_C = 5.0$
Angstroms. As for all periodic HFX calculations, a limited range potential is required. In case of
k-point calculations, a sphere with the radius of $R_C$ must fit the in the BvK supercell
(equivalent of L/2 requirement in $\Gamma$-point calculations). If the requirement is not met, a
warning is issued. The ideal truncation radius to be used for converged results is system dependend,
and careful testing should be done. In practice, most calculations are converged from $R_C = 6.0$
Angstroms on.

Note that there is no specification for the RI metric in the input file below. The default is to use
the same operator for the HFX potential and the RI metric (TC with $R_C = 5.0$ Angstroms, in this
case). Contrary to $\Gamma$-point calculations, using longer ranged RI metric does not dramatically
affect speed in RI-HFXk, while insuring best possible accuracy.

As for most HFX calculations, the SCF convergence can be spedup by restarting from a converged PBE
wavefunction. Note that in k-point restart files, the real-space density matrices are dumped.
However, many more images are required for HFX calculation than for PBE, due to the non-locality of
exact-exchange. Using a very tight value of
[EPS_PGF_ORB](#CP2K_INPUT.ATOM.METHOD.XC.HF.RI.EPS_PGF_ORB) (e.g. $1.0\times 10^{-12}$) in the
initial PBE calculation leads to a lot of images there as well. An example is provided in the
example file [bundle](https://www.cp2k.org/_media/howto:ri_hfx_examples.zip). The example bellow
takes about 5 minutes to run on 32 CPUs if restarted from a PBE wavefunction, and 10 minutes
otherwise.

```none
&GLOBAL
  PROJECT graphene_kp
  RUN_TYPE ENERGY
&END GLOBAL
&FORCE_EVAL
  &DFT
    BASIS_SET_FILE_NAME BASIS_pob
    POTENTIAL_FILE_NAME POTENTIAL
    SORT_BASIS EXP
    AUTO_BASIS RI_HFX MEDIUM
    !restarting from a converged PBE calculation lead to less SCF steps
    WFN_RESTART_FILE_NAME graphene_pbe-RESTART.kp

    !Turning on the ADMM approximation
    &AUXILIARY_DENSITY_MATRIX_METHOD
      ADMM_TYPE ADMMS
    &END AUXILIARY_DENSITY_MATRIX_METHOD

    &QS
      !sometimes necessary when running small systems with a lot of CPUs
      PW_GRID_BLOCKED FALSE
      METHOD GAPW
      !needs to be the same value as that in RI%EPS_PGF_ORB
      EPS_PGF_ORB 1.0E-6
    &END  QS

    &MGRID
      CUTOFF 600
      REL_CUTOFF 60
      NGRIDS 5
    &END MGRID

    &SCF
      EPS_SCF 1.0E-06
      MAX_SCF 50
      !typically need lower threshold to start DIIS with k-points
      EPS_DIIS 0.05
      SCF_GUESS RESTART
    &END SCF

    &XC
      &XC_FUNCTIONAL
        &PBE
          SCALE_X 0.75
        &END
      &END XC_FUNCTIONAL
      &HF
        FRACTION 0.25
        &RI
          KP_NGROUPS 16
          !using a smaller than default EPS_PGF_ORB allows for a
          !more accurate calculation with a larger local RI basis
          EPS_PGF_ORB 1.0E-6
        &END RI
        &INTERACTION_POTENTIAL
          !Always use a limited ranged potential in PBCs
          POTENTIAL_TYPE TRUNCATED
          CUTOFF_RADIUS 5.0
        &END INTERACTION_POTENTIAL
      &END HF
    &END XC
    &KPOINTS
      SCHEME MONKHORST-PACK 19 19 1
    &END KPOINTS
    &PRINT
      &BAND_STRUCTURE
        ADDED_MOS 5
        &KPOINT_SET
          NPOINTS 50
          SPECIAL_POINT GAMMA 0.0000000000 0.0000000000 0.0000000000
          SPECIAL_POINT M 0.5000000000 0.0000000000 0.0000000000
          SPECIAL_POINT K 0.3333333333 0.3333333333 0.0000000000
          SPECIAL_POINT GAMMA 0.0000000000 0.0000000000 0.0000000000
        &END  KPOINT_SET
        FILE_NAME graphene_kp.bs
      &END BAND_STRUCTURE
    &END PRINT
  &END DFT
  &SUBSYS
    &CELL
      !enough space between 2 sheets of graphene not to interact
      ABC 2.46 2.46 20.000
      ALPHA_BETA_GAMMA 90.0 90.0 120.0
    &END CELL
    &COORD
      SCALED
      C 0.3333333 0.6666667 0.000
      C 0.6666667 0.3333333 0.000
    &END COORD
    &KIND C
      BASIS_SET pob-TZVP-rev2
      BASIS_SET AUX_FIT pob-DZVP-rev2
      POTENTIAL ALL
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
```

## Important input parameters

There are a few important input parameters for RI-HFXk calculations:

- [EPS_FILTER](#CP2K_INPUT.ATOM.METHOD.XC.HF.RI.EPS_FILTER): the filtering threshold for sparse
  tensors. Works the same way as for $\Gamma$-point calculations (see above).
- [RI_METRIC](#CP2K_INPUT.ATOM.METHOD.XC.HF.RI.RI_METRIC): using the default value for the RI
  metric, which correspond to the choice of HFX potential, is the way to go. It insures best
  possible accuracy, while only marginally increasing the costs.
- [KP_NGROUPS](#CP2K_INPUT.ATOM.METHOD.XC.HF.RI.KP_NGROUPS): this is a performance keyword. During
  the calculation of the real-space exact-exchange matrices, the work is split among MPI
  subcommunicators. Using more groups drastically speeds up the calculation (efficienctly up to 16
  groups, reasonably up to 32). This comes with a memory overhead though, as some data must be
  replicated on each subgroup. The total number of MPI ranks must be divisible by the number of
  groups.
- [EPS_PGF_ORB](#CP2K_INPUT.ATOM.METHOD.XC.HF.RI.EPS_PGF_ORB): generally determines the range of
  GTOs in the AO basis. As such, it also determines the extent of the local atom-specific RI basis
  used for RI-HFXk. The default value of $1.0\times 10^{-5}$ has proven to be accurate and fast.
- [KP_USE_DELTA_P](#CP2K_INPUT.ATOM.METHOD.XC.HF.RI.KP_USE_DELTA_P): when set to .TRUE. (default
  value), the next SCF step is calculated using the density matrix difference, rather than the full
  new density matrix. This helps with computational efficiency by increasing sparsity. If your
  calculation struggles to converge, you can try to turn this off.

The rest of the [HF/RI](#CP2K_INPUT.ATOM.METHOD.XC.HF.RI) input parameters related to k-point
sampling (all with a KP\_ prefix) have little to no impact, and their default values are good
enough.
