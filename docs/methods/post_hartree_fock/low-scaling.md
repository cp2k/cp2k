# Low-scaling post Hartree-Fock

## Brief theory recap

The low-scaling post-HF methods implemented in CP2K are SOS-MP2 and (direct) RPA. To achieve the
desired low-scaling (sub-cubic), all equations are formulated in the atomic orbital (AO) basis,
using the resolution-of-the-identity technique with a short range metric, and a numerical Laplace
transform. A central quantity of post-HF is the collection of 4-center 2-electron electron repulsion
integrals (ERIs). In the AO basis and within the RI approximation:

$$
(\mu\sigma|\nu\lambda) = (\mu\sigma\lfloor P)\ (P\lfloor Q)^{-1}\ (Q|R)\ (R\lfloor S)^{-1}\ (S\lfloor \nu\lambda)
$$

where $\lfloor$ denotes the short range RI metric, and $P,Q,R,S$ are RI basis elements. A short
range metric increases the sparsity of the 3-center RI integrals $(\mu\sigma\lfloor P)$, resulting
in more efficient sparse tensor contractions. The shortest possible ranged RI metric is the overlap.
Using the truncated Coulomb with a short truncation radius (1.5-2.0 Angstroms) is recommended in
sparser systems (e.g. water), while the overlap is typically sufficient for dense solids. The 3- and
2-center RI ERIs are calculated analytically with the libint library. The 2-center potential ERIs,
$(Q|R)$, are calculated numerically with the GPW scheme, because it uses the (very) long range 1/r
Coulomb interaction.

During post-HF calculations, the ERIs discussed above are contracted with various other quantities,
with a specific order of operations for maximum performance. The details can be found in
[](#Wilhelm2016b) for the initial RPA energy implementation. [](#Bussy2023) covers the detailed
implementation of SOS-MP2 and RPA nuclear gradients. Calculations for the energy, forces, MD, cell
and gemoetry optimization are possible. While these low-scaling methods have more favorable scaling
than their standard implementations (sub-cubic vs quartic), they suffer from a heavy prefactor.
Therefore, they are only interesting for large to very large systems (hundreds of atoms), and will
still consume a lot of resources. For smaller system, the standard MP2 and RPA implementations are
recomended (they also have gradients implemented).

## Simple examples

We cover two simple examples of low-scaling post-HF calculations: SOS-MP2 water MD, and titania
RPA@PBE0 cell optimization. These examples are meant to be representative of the available options,
and the parameters to pay attention to. They are also relatively afforable (although the titania
input might require a few nodes on an HPC). To maintain a resonable cost, the systems are small and
with low quality basis sets. In practice, such small systems should be simulated with the standard
quartic scaling implementation. Moreover, it is never recommended to run post-HF calculations with
double-zeta quality basis sets. All files required to run these examples can be found
[here](https://www.cp2k.org/_media/howto:posthf_examples.zip).

### SOS-MP2 liquid water

This example is a short 3 steps MD liquid water calculation (32 molecules). The RI-MP2 optimized
cc-TZ and RI_TZ basis sets are used. Whenever available, one should use pre-optimized RI basis sets
for better performance. By default, the RI basis is generated on the fly. While it yields accurate
results, the automatically generated RI basis sets are larger, and therefore less efficient. Using a
short range RI metric improves efficiency by introducing sparsity. However, in somewhat sparse
systems like water, the overlap metric lacks accuracy ([](#Bussy2023)). Using the truncated Coulomb
operator with a cutoff radius of 1.5-2.0 Angstroms as a RI metric is safer.

All post-HF calculations need an accurate and well converged SCF as a starting point. In this case,
a pure Hartree-Fock calculation is performed. Note the use of the truncated Coulomb operator as HFX
potential, with cutoff radius \< L/2, as always required in periodic HF. Here, the value of 4.5
Angstroms is somewhat small, and simulating a larger water box would be advised.

A numerical Laplace transform is used to reduce the scaling of the method. The integral is replaced
by a weighted sum over a few grid points, placed according to the MINIMAX algortihm. Typically,
using 6-8 points is enough. CP2K can go up to 20 points, although more points do not necessarily
increase accuracy (and sometime bring numerical instability). The MINIMAX quadrature depends on the
band gap of the system and the total spread of the SCF eigenvalues. If the seleceted number of grid
point is innapropriate, a warning is issued.

```none
&GLOBAL  
  PROJECT water32  
  RUN_TYPE MD  
  PRINT_LEVEL MEDIUM  
&END GLOBAL  
&MOTION  
  &MD  
    STEPS 3  
  &END MD  
&END MOTION  
&FORCE_EVAL  
  &DFT  
    !cc-TZ: RI-MP2 optimized basis sets  
    BASIS_SET_FILE_NAME BASIS_RI_cc-TZ  
    POTENTIAL_FILE_NAME POTENTIAL  
    SORT_BASIS EXP  
  
    &MGRID  
      CUTOFF 600  
      REL_CUTOFF 50  
      NGRIDS 5  
    &END MGRID  
  
    &SCF  
      SCF_GUESS RESTART  
      EPS_SCF 1.0E-6  
      MAX_SCF 40  
    &END SCF  
  
    &XC  
      &XC_FUNCTIONAL NONE  
      &END XC_FUNCTIONAL  
      &HF  
        FRACTION 1.0  
        &INTERACTION_POTENTIAL  
          !TC potential with cutoff < L/2 for periodic HFX  
          POTENTIAL_TYPE TRUNCATED  
          CUTOFF_RADIUS 4.5  
        &END INTERACTION_POTENTIAL  
        &MEMORY  
          !maximum memory allocated to HFX ERI storage, per MPI rank  
          !the optimal number depends on the specifics of the computer  
          MAX_MEMORY 4000  
        &END MEMORY  
      &END HF
      &WF_CORRELATION  
        !explicit opposite-spin scaling  
        SCALE_S 1.3  
        &RI_SOS_MP2  
          !MINIMAX quadrature by default  
          QUADRATURE_POINTS 6  
        &END RI_SOS_MP2
        !Enabling low-scaling  
        &LOW_SCALING  
          MEMORY_CUT 3  
        &END LOW_SCALING  
        &RI  
          &RI_METRIC  
            !Short range RI metric for SOS-MP2  
            POTENTIAL_TYPE TRUNCATED  
            CUTOFF_RADIUS 1.5  
          &END RI_METRIC  
        &END RI  
        &INTEGRALS  
          ERI_METHOD GPW  
          &WFC_GPW  
            !Safe yet faster than default values  
            CUTOFF 200  
            REL_CUTOFF 40  
          &END  WFC_GPW
        &END INTEGRALS
      &END WF_CORRELATION 
    &END XC  
  &END DFT  
  &SUBSYS  
    &CELL  
      ABC 9.8528 9.8528 9.8528  
    &END CELL  
    &TOPOLOGY  
      COORD_FILE_FORMAT XYZ  
      COORD_FILE_NAME ./H2O-32.xyz  
    &END TOPOLOGY  
    &KIND H  
      BASIS_SET cc-DZ  
      BASIS_SET RI_AUX RI_DZ  
      POTENTIAL GTH-HF  
    &END KIND  
    &KIND O  
      BASIS_SET cc-DZ  
      BASIS_SET RI_AUX RI_DZ  
      POTENTIAL GTH-HF  
    &END KIND  
  &END SUBSYS  
&END FORCE_EVAL
```

### RPA bulk titania

This second example is 2 steps of a cell optimization of titania at the RPA@PBE0 level (72 atoms).
All low-scaling post-HF methods implemented in CP2K can be used with the ADMM (ADMM2 flavor)
approximation to speed up the Hartree-Fock calculations for the SCF, the response forces, and EXX
(in case of RPA). Here, the admm-dzp basis is used as an auxiliary ADMM basis set, and ccGRB-D as
the primary basis. The RI basis is generated automatically.

The overlap RI metric is used in this example. For dense systems, this choice of metric offers a
good tradeoff between accuracy and efficiency. Note that the `&HF` input section is repeated; it
appears once for the SCF specification, and once for the RPA. In RPA@PBE0, the SCF is first
converged at the PBE0 level of theory. Then, the RPA correlation energy is computed using the PBE0
orbitals and eigenvalues, and finally, the PBE0 exchange-correlation energy is replaced by the exact
exchange energy (EXX) calculated with the same orbitals. The `&HF` section in RPA refers to that
last step. If taken to be the same as the SCF (except for `FRACTION` and `MEMORY`), the stored ERIs
can be reused and rescaled instead of recomputed, thus saving precious computational time. The
`ADMM` keyword in RPA specifies that the EXX should be calculated within the ADMM approximation.
Note the truncated Coulomb potential used for HFX, and its cutoff radius of 4.5 Angstroms (\< L/2),
necessary for periodic HFX calculations. This cutoff radius is on the shorter size, and using a
larger simulation cell would be advised.

The `MINIMAX_QUADRATURE` is explicitely requested in the input. This is a necessary step for
efficiency, as the default grid follows a Clenshaw-Curtis scheme, which requires many more
integration points for the same accuracy.

```none
&GLOBAL  
  PROJECT TiO2  
  RUN_TYPE CELL_OPT  
  PRINT_LEVEL MEDIUM  
&END GLOBAL  
&MOTION  
  &CELL_OPT  
    MAX_ITER 2  
  &END CELL_OPT  
&END MOTION  
&FORCE_EVAL  
  STRESS_TENSOR ANALYTICAL  
  &DFT  
    BASIS_SET_FILE_NAME BASIS_ccGRB_UZH  
    BASIS_SET_FILE_NAME BASIS_ADMM_UZH  
    POTENTIAL_FILE_NAME POTENTIAL_UZH  
    SORT_BASIS EXP  
    !automatically generated RI basis set  
    AUTO_BASIS RI_AUX SMALL  
  
    !enabling the ADMM2 approximation  
    &AUXILIARY_DENSITY_MATRIX_METHOD  
      METHOD BASIS_PROJECTION  
      ADMM_PURIFICATION_METHOD NONE  
      EXCH_CORRECTION_FUNC PBEX  
    &END AUXILIARY_DENSITY_MATRIX_METHOD   
  
    &MGRID  
      CUTOFF 600  
      REL_CUTOFF 50  
      NGRIDS 5  
    &END MGRID  
  
    &SCF  
      EPS_SCF 1.0E-6  
      MAX_SCF 20  
      &OT  
        PRECONDITIONER FULL_ALL  
        MINIMIZER DIIS  
      &END OT  
      &OUTER_SCF  
        MAX_SCF 5  
        EPS_SCF 1.0E-6  
      &END OUTER_SCF  
    &END SCF  
  
    &XC  
      &XC_FUNCTIONAL  
        &PBE  
          SCALE_X 0.75  
        &END PBE  
      &END XC_FUNCTIONAL  
      &HF  
        FRACTION 0.25  
        !always use a short range potential < L/2 in periodic HFX  
        &INTERACTION_POTENTIAL  
          POTENTIAL_TYPE TRUNCATED  
          CUTOFF_RADIUS 4.5  
        &END INTERACTION_POTENTIAL  
        &MEMORY  
          !maximum memory allocated to HFX ERI storage, per MPI rank  
          !the optimal number depends on the specifics of the computer  
          MAX_MEMORY 4000  
        &END MEMORY  
      &END HF
      &WF_CORRELATION  
        &RI_RPA  
          MINIMAX_QUADRATURE  
          QUADRATURE_POINTS 6  
          !calculate EXX with ADMM, using the same HF section as in SCF  
          !so that the integrals can be resued without recomputing  
          ADMM  
          &HF  
            FRACTION 1.0  
            &INTERACTION_POTENTIAL  
              POTENTIAL_TYPE TRUNCATED  
              CUTOFF_RADIUS 4.5  
            &END INTERACTION_POTENTIAL  
          &END HF  
        &END RI_RPA  
        !enabling low-scaling  
        &LOW_SCALING  
          MEMORY_CUT 3  
        &END LOW_SCALING
        &RI  
          !overlap RI metric is appropriate for dense systems  
          &RI_METRIC  
            POTENTIAL_TYPE IDENTITY  
          &END RI_METRIC 
        &END RI  
        &INTEGRALS  
          ERI_METHOD GPW  
          !Safe yet faster than default values  
          &WFC_GPW  
            CUTOFF 200  
            REL_CUTOFF 40  
          &END WFC_GPW  
        &END INTEGRALS  
      &END WF_CORRELATION
    &END XC  
  &END DFT  
  &SUBSYS  
    &CELL  
      ABC 9.330 9.330 9.107  
    &END CELL  
    &TOPOLOGY  
      COORD_FILE_FORMAT XYZ  
      COORD_FILE_NAME ./TiO2.xyz  
    &END TOPOLOGY  
    &KIND Ti  
      BASIS_SET ccGRB-D-q12  
      BASIS_SET AUX_FIT admm-dz-q12  
      POTENTIAL GTH-PBE0-q12  
    &END KIND  
    &KIND O  
      BASIS_SET ccGRB-D-q6  
      BASIS_SET AUX_FIT admm-dz-q6  
      POTENTIAL GTH-PBE0-q6  
    &END KIND  
  &END SUBSYS  
&END FORCE_EVAL
```

## Important input parameters

There are a few importnat input parameters for low-scaling post-HF calculations:

- `SORT_BASIS`, in the DFT input section, should be set to EXP. This way, basis elements are sorted
  according to their exponents, leading to increased sparsity and performance.
- `EPS_FILTER`, in the LOW_SCALING input section. The bottleneck of the low-scaling post-HF
  calculations consists of sparse tensor contractions. This parameter is the threshold for block
  sparsity during calculations, with a safe default of $1.0\times10^{-9}$. A looser threshold might
  significantly speedup calculations (not recommended to go higher than $1.0\times10^{-8}$).
- `MEMORY_CUT`, in the LOW_SCALING input section. This keyword influences the batching strategy for
  large tensor contractions. In post-HF, some tensor contractions are done in multiple steps, with
  large intermediate results. Storing these intermediates can lead to memory shortage. The value of
  `MEMORY_CUT` (default of 5) indicates how large tensors are split into batches, such that smaller
  intermediates can be stored. Note that MEMORY_CUT 3 does not mean that the total memory
  consumption of the program is divided by three (initial and final tensors are not affected, as
  well as the data concerning the rest of the calculation). A higher value reduces the memory
  footprint, but leads to performance overheads.
- `RI METRIC`: the choice of RI metric is crucial for performance in periodic calculations.
  A\
  shorter ranged RI metric will generally be more efficient, while a longer ranged metric (e.g.
  TC\
  with a cutoff radius of 1.5 Angstrom) can be more accurate. See [](#Bussy2023) for a full
  discussion.
