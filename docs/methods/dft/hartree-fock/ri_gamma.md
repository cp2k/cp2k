# HFX-RI for Γ-Point (non-periodic)

The resolution-of-the-identity (RI) technique is implemented for various methods in CP2K, each time
in a slightly different flavor. Hartree-Fock exchange (HFX) has two distinct RI implementations: One
for [k-point sampling](./ri_kpoints) and one for $\Gamma$-point calculations, which is covered on
this page.

The implementations for standard RI-HFX is described in [](#Bussy2023). All files necessary to run
the examples discussed below can be found
[here](https://www.cp2k.org/_media/howto:ri_hfx_examples.zip).

## Application domain

The CP2K RI-HFX implementation has proven to be particularly efficient for the calculation of
molecules with large basis sets, and small/medium size solids. For large and sparse periodic systems
(such as water), the original 4-center HFX implementation is more efficient. The ADMM method, known
for greatly speeding up HFX calculations, can be seamlessly used with RI-HFX. Nuclear forces and
stress tensors are available.

## Brief theory recap

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

## Simple examples

The following two examples illustrate the use of RI-HFX where it is most efficient: for a molecule
with a large basis, and a small dense system. Different input options are shown for maximum coverage
(ADMM, hybdrid-DFT, short range RI metric, automatically generated RI basis, etc.)

### Molecule

In this example, a geometry optimization of the glycine molecule is calculated at the Hartree-Fock
level of theory, with the cc-pVQZ basis set (taken from https://www.basissetexchange.org). The RI
basis set, referred to as `RI_HFX` in the [KIND](#CP2K_INPUT.FORCE_EVAL.SUBSYS.KIND) input section,
is explicitly provided by the user in this example. It takes about 70 seconds to run on 16 CPUs, per
geometry optimization step.

There are a couple take home messages with this example:

- If optimized RI basis sets (such as cc-pVQZ-JKFIT in this case) exist, one should use them for
  more efficiency. They are smaller than their automatically generated counter parts, albeit
  sometimes not as accurate.
- In non-periodic calculations, use the 1/r default Coulomb interaction for the HFX potential and
  the RI metric. Short range RI metric introduce sparsity and reduces the cost of calculating
  $(\mu\sigma\lfloor P)$, which is only really useful in PBCs.

In this specific example, RI-HFX is much more efficient than the original 4-center implementation.
You can try running the same input with the [HF/RI](#CP2K_INPUT.ATOM.METHOD.XC.HF.RI) section
commented out for a demonstration.

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

### Solid

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
- Setting the keyword [SORT_BASIS](#CP2K_INPUT.FORCE_EVAL.DFT.SORT_BASIS) to `EXP` helps create
  sparsity and improves performance
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

## Important input parameters

There are some few input parameters of particular importance when running RI-HFX, mostly concerning
the efficiency to accuracy balance. Here is a list:

- [RI_METRIC](#CP2K_INPUT.ATOM.METHOD.XC.HF.RI.RI_METRIC): the choice of RI metric is crucial for
  performance in periodic calculations. A shorter ranged RI metric will generally be more efficient,
  while a longer ranged metric (e.g. TC with a cutoff radius of 1.5 Angstrom) can be more accurate.
  See [](#Bussy2023) for a full discussion.
- [EPS_FILTER](#CP2K_INPUT.ATOM.METHOD.XC.HF.RI.EPS_FILTER): the bottleneck of the RI-HFX
  calculation consists of sparse tensor contractions. This parameter is the threshold for block
  sparsity during calculations, with a safe default of $1.0\times10^{-9}$. A looser threshold might
  significantly speedup calculations (not recommended to go higher than $1.0\times10^{-8}$).
- [MEMORY_CUT](#CP2K_INPUT.ATOM.METHOD.XC.HF.RI.MEMORY_CUT): this keyword influences the batching
  strategy for large tensor contractions. In RI-HFX, some tensor contractions are done in multiple
  steps, with large intermediate results. Storing these intermediates can lead to memory shortage.
  The value of [MEMORY_CUT](#CP2K_INPUT.ATOM.METHOD.XC.HF.RI.MEMORY_CUT) (default of 3) indicates
  how large tensors are split into batches, such that smaller intermediates can be stored. Note that
  `MEMORY_CUT 3` does not mean that the total memory consumption of the program is divided by three
  (initial and final tensors are not affected, as well as the data concerning the rest of the
  calculation). A higher value reduces the memory footprint, but leads to performance overheads.

Most other keywords have little to no importance, and their default values are fine.
