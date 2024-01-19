# RI Hartree-Fock Exchange (incl. k-points)

The resolution-of-the-identity (RI) technique is implemented for various methods in CP2K, each time
in a slightly different flavor. Hartree-Fock exchange (HFX) has two distinct RI implementations, one
for $\Gamma$-point (and non-periodic) calculations, and one for k-point sampling. This page contains
a brief description of each method, with illustrative examples. The implementations are described in
two publications; [](#Bussy2023) for standard RI-HFX, and an accepted, but not yet published paper
for RI-HFXk (with k-point sampling). All files necessary to run the examples discussed below can be
found [here](https://www.cp2k.org/_media/howto:ri_hfx_examples.zip).

## $\Gamma$-point and non-periodic RI-HFX

### Application domain

The CP2K RI-HFX implementation has proven to be particularly efficient for the calculation of
molecules with large basis sets, and small/medium size solids. For large and sparse periodic systems
(such as water), the original 4-center HFX implementation is more efficient. The ADMM method, known
for greatly speeding up HFX calculations, can be seamlessly used with RI-HFX. Nuclear forces and
stress tensors are available.

### Brief theory recap

For $\Gamma$-point calculations, the exact-exchange contribution to the KS matrix is obtained by
contraction of the density matrix with 2-electron 4-center electron repulsion integrals (ERIs):

$$
K_{\mu\nu} = P_{\sigma\lambda} \sum_{\mathbf{a, b, c}}(\mu^\mathbf{0}\sigma^\mathbf{a}|\nu^\mathbf{b}\lambda^\mathbf{a+c})
$$

where the AO indices are implicitly summed (Einstein convention). $\mathbf{a, b, c}$ correspond to
periodic images of the simulation cell. The sum over periodic images can be done once and for all,
resulting in a single 4-center quantity $(\mu\sigma|\nu\lambda)$. For molecules in non-periodic
boundary conditions, $\mathbf{a}=\mathbf{b}=\mathbf{c}=\mathbf{0}$. In the RI approximation, 2- and
3-center ERIs are used instead:

$$
K_{\mu\nu} = P_{\sigma\lambda} \Big[ (\mu\sigma\lfloor P)\ (P\lfloor Q)^{-1}\ (Q|R)\ (R\lfloor S)^{-1}\ (S\lfloor \nu\lambda) \Big]
$$

with $P,Q,R,S$ GTOs forming a global RI basis, ideally covering the same space as all possible AO
products. $\lfloor$ represents the RI metric, typically taken to have a smaller range than the HFX
potential. The shortest range RI metric is the overlap, resulting in the well known RI-SVS
approximation.

When using RI-HFX, the computational bottleneck shifts from calculating 4-center 2-electron ERIs to
sparse-tensor contraction. Typically, RI-HFX is less memory intensive and GPU accelerated. However,
it's asymptotical scaling is $\sim\mathcal{O}(N^3)$, while 4-center HFX is linear. Therefore, RI-HFX
is more suitable for smaller dense systems and molecules.

### Simple examples

The following two examples illustrate the use of RI-HFX where it is most efficient: for a molecule
with a large basis, and a small dense system. Different input options are shown for maximum coverage
(ADMM, hybdrid-DFT, short range RI metric, automatically generated RI basis, etc.)

#### Molecule

In this example, a geometry optimization of the glycine molecule is calculated at the Hartree-Fock
level of theory, with the cc-pVQZ basis set (taken from https://www.basissetexchange.org). The RI
basis set, referred to as `RI_HFX` in the `&KIND` input section, is explicitly provided by the user
in this example. It takes about 70 seconds to run on 16 CPUs, per geometry optimization step.

There are a couple take home messages with this example:

- If optimized RI basis sets (such as cc-pVQZ-JKFIT in this case) exist, one should use them for
  more efficiency. They are smaller than their automatically generated counter parts, albeit
  sometimes not as accurate.
- In non-periodic calculations, use the 1/r default Coulomb interaction for the HFX potential and
  the RI metric. Short range RI metric introduce sparsity and reduces the cost of calculating
  $(\mu\sigma\lfloor P)$, which is only really useful in PBCs.

In this specific example, RI-HFX is much more efficient than the original 4-center implementation.
You can try running the same input with the `&HF%RI` section commented out for a demonstration.

```none
&GLOBAL  
  PROJECT glycine  
  PRINT_LEVEL MEDIUM  
  RUN_TYPE GEO_OPT  
&END GLOBAL  
&FORCE_EVAL  
  METHOD Quickstep  
  &DFT  
    BASIS_SET_FILE_NAME BASIS_cc-pVQZ  
    POTENTIAL_FILE_NAME POTENTIAL  
    !Sort basis function accoring to their exponent for more sparsity 
    SORT_BASIS EXP  
  
    &MGRID  
      CUTOFF 500  
      REL_CUTOFF 50  
      NGRIDS 5  
    &END MGRID  
  
    &QS  
      !all-electron calculations require GAPW  
      METHOD GAPW  
    &END QS  
  
    &POISSON  
      !non-periodic calculation for this molecule  
      PERIODIC NONE  
      PSOLVER WAVELET  
    &END  
  
    &SCF  
      EPS_SCF 1.0E-6  
      MAX_SCF 50  
      &END SCF  
    &XC  
      &XC_FUNCTIONAL NONE  
      &END XC_FUNCTIONAL  
      &HF  
        !Pure RI Hartree-Fock calculation using only defaults:  
        ! -HFX potential is the 1/r Coulomb interaction  
        ! -RI metric is also the 1/r Coulomb interaction  
        ! -Default accuracy parameters (good in most cases)  
        &RI  
        &END RI  
      &END HF  
    &END XC  
  &END DFT  
  &SUBSYS  
    &CELL  
      ABC 10.0 10.0 10.0  
      PERIODIC NONE  
    &END CELL  
    &COORD  
      C -0.04879702 -0.00000000 1.40419128  
      N -1.35021542 0.00000000 2.04225544  
      C -0.04354337 0.00000000 -0.12235209  
      O -1.02422569 -0.00000000 -0.83489570  
      O 1.22983691 0.00000000 -0.61028238  
      H 1.14837668 -0.00000000 -1.58391528  
      H 0.53209836 -0.87421885 1.73662058  
      H 0.53209836 0.87421885 1.73662058  
      H -1.88873390 0.81280508 1.74087629  
      H -1.88873390 -0.81280508 1.7408762
    &END COORD
    &TOPOLOGY  
      !Always a good idea to put the molecule in the middle of the simulation cell in non-PBCs  
      &CENTER_COORDINATES  
      &END CENTER_COORDINATES  
    &END TOPOLOGY  
    &KIND C  
      BASIS_SET cc-pVQZ  
      BASIS_SET RI_HFX cc-pVQZ-JKFIT  
      POTENTIAL ALL  
    &END KIND  
   &KIND O  
      BASIS_SET cc-pVQZ  
      BASIS_SET RI_HFX cc-pVQZ-JKFIT  
      POTENTIAL ALL  
    &END KIND  
    &KIND N  
      BASIS_SET cc-pVQZ  
      BASIS_SET RI_HFX cc-pVQZ-JKFIT  
      POTENTIAL ALL  
    &END KIND  
    &KIND H  
      BASIS_SET cc-pVQZ  
      BASIS_SET RI_HFX cc-pVQZ-JKFIT  
      POTENTIAL ALL  
    &END KIND  
  &END SUBSYS  
&END FORCE_EVAL
```

#### Solid

In this example, the energy of a 64 atoms cell of bulk silicon is calculated at the ADMM-PBE0 level.
Because there exists no pre-optimized RI basis set corresponding to `admm-dzp`, the RI basis is
generated on the fly. Here, the 25% exact-exchange fraction of the PBE0 functional is calculated
with RI-HFX and ADMM. The HFX potential is taken to be the truncated Coulomb potential with a cutoff
radius of 5.4 Angstroms, which is smaller than half the cell size (as necessary for all periodic HFX
calculations). In this case, the RI metric is selected to be the overlap, for optimal computational
efficiency. Note that sparser systems might require a longer ranged RI metric (see discussion in
[](#Bussy2023)). This example runs in about 10 minutes on 32 CPUs.

There are a handful of take home messages with this example:

- Always use HFX potentials that decay within half the simulation cell in periodic HFX calculations
- Using short range RI metrics is crucial for performance
- ADMM is available in all its flavors, and it helps reducing computational costs
- Using the keyword `SORT_BASIS EXP` helps create sparsity and improves performance
- Automatically generated RI basis sets are very accurate, although larger and less efficient than
  pre-optimized ones

```none
&GLOBAL  
  PROJECT Si64  
  RUN_TYPE ENERGY  
&END GLOBAL  
&FORCE_EVAL  
  &DFT  
    BASIS_SET_FILE_NAME BASIS_ccGRB_UZH  
    BASIS_SET_FILE_NAME BASIS_ADMM_UZH  
    POTENTIAL_FILE_NAME POTENTIAL_UZH
    !sort the basis function according to their exponent for more sparsity  
    SORT_BASIS EXP
    !generate the RI_HFX basis set on the fly
    AUTO_BASIS RI_HFX SMALL  
    
    !turn on ADMM
    &AUXILIARY_DENSITY_MATRIX_METHOD  
      ADMM_TYPE ADMMS  
    &END AUXILIARY_DENSITY_MATRIX_METHOD  
  
    &MGRID  
      CUTOFF 600  
      REL_CUTOFF 50  
      NGRIDS 5  
    &END MGRID  
  
    &SCF  
      EPS_SCF 1.0E-6  
      MAX_SCF 40  
    &END SCF  
  
    &XC  
      &XC_FUNCTIONAL  
        &PBE  
          SCALE_X 0.75  
        &END PBE 
            &END XC_FUNCTIONAL  
            &HF  
                FRACTION 0.25  
                &INTERACTION_POTENTIAL
                  !Important to use a limited range potential in periodic HFX  
                    POTENTIAL_TYPE TRUNCATED
                    !5.4 < half cell dimension
                    CUTOFF_RADIUS 5.4  
                &END  
                &RI
                  !overlap metric for maximal efficiency  
                    RI_METRIC IDENTITY  
                &END RI  
            &END HF  
        &END XC  
    &END DFT  
    &SUBSYS  
        &CELL  
            ABC 10.861395 10.861395 10.861395  
        &END CELL  
        &TOPOLOGY  
            COORD_FILE_FORMAT XYZ  
            COORD_FILE_NAME ./Si64.xyz  
        &END TOPOLOGY  
        &KIND Si  
            BASIS_SET ccGRB-D-q4  
            BASIS_SET AUX_FIT admm-dzp-q4  
            POTENTIAL GTH-PBE0-q4  
        &END KIND  
    &END SUBSYS  
&END FORCE_EVAL
```

### Important input parameters

There are some few input parameters of particular importance when running RI-HFX, mostly concerning
the efficiency to accuracy balance. Here is a list:

- `RI_METRIC`: the choice of RI metric is crucial for performance in periodic calculations. A
  shorter ranged RI metric will generally be more efficient, while a longer ranged metric (e.g. TC
  with a cutoff radius of 1.5 Angstrom) can be more accurate. See [](#Bussy2023) for a full
  discussion.
- `EPS_FILTER`: the bottleneck of the RI-HFX calculation consists of sparse tensor contractions.
  This parameter is the threshold for block sparsity during calculations, with a safe default of
  $1.0\times10^{-9}$. A looser threshold might significantly speedup calculations (not recommended
  to go higher than $1.0\times10^{-8}$).
- `MEMORY_CUT`: this keyword influences the batching strategy for large tensor contractions. In
  RI-HFX, some tensor contractions are done in multiple steps, with large intermediate results.
  Storing these intermediates can lead to memory shortage. The value of `MEMORY_CUT` (default of 3)
  indicates how large tensors are split into batches, such that smaller intermediates can be stored.
  Note that MEMORY_CUT 3 does not mean that the total memory consumption of the program is divided
  by three (initial and final tensors are not affected, as well as the data concerning the rest of
  the calculation). A higher value reduces the memory footprint, but leads to performance overheads.

Most other keywords have little to no importance, and their default values are fine.

## RI-HFXk with k-point sampling

### Application domain

The RI-HFX implementation with k-point sampling (RI-HFXk) is optimized for the simulation of small
unit cells and dense k-point meshes. It is much more efficient than equivalent $\Gamma$-point
supercell calculations.

### Brief theory recap

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

### Simple example (graphene band structure)

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

Note that in this input file, we select a value of $1.0\times 10^{-6}$ for `EPS_PGF_ORB`. This
parameter controls the range of AOs, and threfore the extent of the local atom-specific RI basis
sets. This value leads to particularly high accuracy. The default of $1.0\times 10^{-5}$ is
typically enough.

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
exact-exchange. Using a very tight value of `EPS_PGF_ORB` (e.g. $1.0\times 10^{-12}$) in the initial
PBE calculation leads to a lot of images there as well. An example is provided in the example file
[bundle](https://www.cp2k.org/_media/howto:ri_hfx_examples.zip). The example bellow takes about 5
minutes to run on 32 CPUs if restarted from a PBE wavefunction, and 10 minutes otherwise.

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

### Important input parameters

There are a few important input parameters for RI-HFXk calculations:

- `EPS_FILTER`: the filtering threshold for sparse tensors. Works the same way as for $\Gamma$-point
  calculations (see above).
- `RI_METRIC`: using the default value for the RI metric, which correspond to the choice of HFX
  potential, is the way to go. It insures best possible accuracy, while only marginally increasing
  the costs.
- `NGROUPS`: this is a performance keyword. During the calculation of the real-space exact-exchange
  matrices, the work is split among MPI subcommunicators. Using more groups drastically speeds up
  the calculation (efficienctly up to 16 groups, reasonably up to 32). This comes with a memory
  overhead though, as some data must be replicated on each subgroup. The total number of MPI ranks
  must be divisible by the number of groups.
- `EPS_PGF_ORB`: generally determines the range of GTOs in the AO basis. As such, it also determines
  the extent of the local atom-specific RI basis used for RI-HFXk. The default value of
  $1.0\times 10^{-5}$ has proven to be accurate and fast.
- `KP_USE_DELTA_P`: when set to .TRUE. (default value), the next SCF step is calculated using the
  density matrix difference, rather than the full new density matrix. This helps with computational
  efficiency by increasing sparsity. If your calculation struggles to converge, you can try to turn
  this off.

The rest of the `&HF%RI` input parameters related to k-point sampling (all with a KP\_ prefix) have
little to no impact, and their default values are good enough.
