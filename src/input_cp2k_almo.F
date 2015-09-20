!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief input for the ALMO SCF section
!> \author Rustam Khaliullin
! *****************************************************************************
MODULE input_cp2k_almo
  USE bibliography,                    ONLY: Khaliullin2007,&
                                             Khaliullin2008,&
                                             Khaliullin2013
  USE input_constants,                 ONLY: &
       almo_deloc_none, almo_deloc_scf, almo_deloc_x, almo_deloc_x_then_scf, &
       almo_deloc_xalmo_1diag, almo_deloc_xalmo_scf, almo_deloc_xalmo_x, &
       almo_frz_crystal, almo_frz_none, almo_scf_diag, almo_scf_pcg, &
       atomic_guess, cg_dai_yuan, cg_fletcher, cg_fletcher_reeves, &
       cg_hager_zhang, cg_hestenes_stiefel, cg_liu_storey, cg_polak_ribiere, &
       cg_zero, molecular_guess, optimizer_diis, optimizer_pcg, prec_default, &
       prec_zero
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_add_subsection,&
                                             section_create,&
                                             section_release,&
                                             section_type
  USE kinds,                           ONLY: dp
  USE string_utilities,                ONLY: s2a
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_almo'

  INTEGER,  PARAMETER, PRIVATE :: optimizer_block_diagonal_diis = 1 
  INTEGER,  PARAMETER, PRIVATE :: optimizer_block_diagonal_pcg = 2 
  INTEGER,  PARAMETER, PRIVATE :: optimizer_xalmo_pcg = 3
  
  PUBLIC :: create_almo_scf_section

CONTAINS

! *****************************************************************************
!> \brief create the almo scf section
!> \param section ...
!> \par History
!>       2011.05 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
SUBROUTINE create_almo_scf_section(section)
    TYPE(section_type), POINTER              :: section

    CHARACTER(len=*), PARAMETER :: routineN = 'create_almo_scf_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,failure)
    CALL section_create(section,"ALMO_SCF",&
         description="Settings for a class of efficient linear scaling methods based "//&
         "on absolutely localized orbitals"//&
         " (ALMOs). ALMO methods are currently restricted to closed-shell molecular systems.",&
         n_keywords=4, n_subsections=3, repeats=.FALSE.,&
         citations=(/Khaliullin2013/))

    NULLIFY (keyword)

    CALL keyword_create(keyword, name="EPS_FILTER",&
         description="Threshold for the matrix sparsity filter",&
         usage="EPS_FILTER 1.e-6", default_r_val=1.e-7_dp)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)
    
    CALL keyword_create(keyword, name="ALMO_SCF_GUESS",&
         description="The method to generate initial ALMOs.",&
         usage="ALMO_SCF_GUESS MOLECULAR",&
         default_i_val=atomic_guess,&
         enum_c_vals=s2a("MOLECULAR", "ATOMIC"),&
         enum_desc=s2a("SCF calculations on single molecules controled by the regular SCF "//&
                       "keywords outside ALMO options. This kind of calculation is expensive "//&
                       "and only recommended if ALMO SCF does not converge from the ATOMIC guess.",&
                       "Superpoisiton of atomic densities."),&
         enum_i_vals=(/molecular_guess,atomic_guess/))
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="ALMO_ALGORITHM",&
         description="Specifies the algorithm to update block-diagonal ALMOs during the SCF procedure.",&
         usage="ALMO_ALGORITHM DIAG",&
         default_i_val=almo_scf_diag,&
         !enum_c_vals=s2a("DIAG", "DM_SIGN","PCG"),&
         enum_c_vals=s2a("DIAG", "PCG"),&
         enum_desc=s2a("DIIS-accelerated diagonalization controlled by ALMO_OPTIMIZER_DIIS. "//&
                       "Recommended for large systems containing small fragments.",&
                       !"Update the density matrix using linear scaling routines. "//&
                       !"Recommended if large fragments are present.",&
                       "Energy minimization with a PCG algorithm controlled by ALMO_OPTIMIZER_PCG."),&
         !enum_i_vals=(/almo_scf_diag,almo_scf_dm_sign,almo_scf_pcg/),&
         enum_i_vals=(/almo_scf_diag,almo_scf_pcg/))
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="DELOCALIZE_METHOD",&
         description="Methods to reintroduce electron delocalization, which is excluded "//&
                     "with the block-diagonal ALMO reference. Electron delocalization can "//&
                     "be computed using either fully delocalized MOs or spatially restricted "//&
                     "ALMOs (called extended ALMOs or XALMOs). All methods below use a PCG "//&
                     "optimizer controlled by XALMO_OPTIMIZER_PCG. The only exception is "//&
                     "the non-iterative XALMO_1DIAG.",&
         usage="DELOCALIZE_METHOD XALMO_X",&
         default_i_val=almo_deloc_xalmo_x,&
         enum_c_vals=s2a("NONE","XALMO_1DIAG","XALMO_X","XALMO_SCF","FULL_X","FULL_SCF","FULL_X_THEN_SCF"),&
         enum_desc=s2a("Neglect electron delocalization",&
         "Correction based on one diagonalization of the spatially projected Hamiltonian (XALMO)",&
         "Single excitation correction (no Hamiltonian re-build) with spatial restrictions (XALMO)",&
         "Self-consistent treatment of delocalization with spatial restrictions (XALMO)",&
         "Single excitation correction (no Hamiltonian re-build) without spatial restrictions",&
         "Self-consistent treatment of delocalization without spatial restrictions",&
         "Single excitation correction followed by full SCF procedure, both without spatial restrictions"),&
         enum_i_vals=(/almo_deloc_none,almo_deloc_xalmo_1diag,almo_deloc_xalmo_x,almo_deloc_xalmo_scf,&
         almo_deloc_x,almo_deloc_scf,almo_deloc_x_then_scf/))
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="XALMO_R_CUTOFF_FACTOR",&
         description="Controls the localization radius of XALMOs: "//&
         !"r0 = r0_factor*(radius(at1)+radius(at2)) + r0_shift",&
         "R_cutoff = XALMO_R_CUTOFF_FACTOR*(radius(at1)+radius(at2))",&
         usage="XALMO_R_CUTOFF_FACTOR 1.6", default_r_val=1.50_dp)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    !CALL keyword_create(keyword, name="DELOCALIZE_EPS_ITER",&
    !     description="Obsolete and to be deleted: use EPS_ERROR in XALMO_OPTIMIZER_PCG",&
    !     usage="DELOCALIZE_EPS_ITER 1.e-5", default_r_val=1.e-5_dp)
    !CALL section_add_keyword(section,keyword)
    !CALL keyword_release(keyword)

    !CALL keyword_create(keyword, name="DELOCALIZE_EPS_LIN_SEARCH",&
    !     description="Obsolete and to be deleted: use EPS_ERROR_LIN_SEARCH in XALMO_OPTIMIZER_PCG",&
    !     usage="DELOCALIZE_EPS_LIN_SEARCH 1.e-6", default_r_val=1.e-7_dp)
    !CALL section_add_keyword(section,keyword)
    !CALL keyword_release(keyword)

    !CALL keyword_create(keyword, name="BLOCKED_MAX_ITER",&
    !     description="Obsolete and to be deleted: use MAX_ITER in ALMO_OPTIMIZER_DIIS or ALMO_OPTIMIZER_PCG",&
    !     usage="BLOCKED_MAX_ITER 200", default_i_val=100)
    !CALL section_add_keyword(section,keyword)
    !CALL keyword_release(keyword)

    !CALL keyword_create(keyword, name="BLOCKED_EPS_ITER",&
    !     description="Obsolete and to be deleted: use EPS_ERROR in ALMO_OPTIMIZER_DIIS or ALMO_OPTIMIZER_PCG",&
    !     usage="BLOCKED_EPS_ITER 1.e-5", default_r_val=1.e-5_dp)
    !CALL section_add_keyword(section,keyword)
    !CALL keyword_release(keyword)

    !CALL keyword_create(keyword, name="DELOCALIZE_MAX_ITER",&
    !     description="Obsolete and to be deleted: use MAX_ITER in XALMO_OPTIMIZER_PCG",&
    !     usage="DELOCALIZE_MAX_ITER 200", default_i_val=100)
    !CALL section_add_keyword(section,keyword)
    !CALL keyword_release(keyword)

    !CALL keyword_create(keyword, name="DOMAIN_LAYOUT_MOS",&
    !     description="Each electron in the system is constrained to its own delocalization domain."//&
    !     " This keyword creates groups of electrons that share the same domain.",&
    !     usage="DOMAIN_LAYOUT_MOS MOLECULAR",&
    !     default_i_val=almo_domain_layout_molecular,&
    !     enum_c_vals=s2a( "ORBITAL", "ATOMIC", "MOLECULAR"),&
    !     enum_desc=s2a("Each electron can have its own delocalization domain",&
    !     "All electrons of an atom are delocalized over the same domain",&
    !     "All electrons of a molecule are delocalized over the same domain."//&
    !     " This is the only implemented option"),&
    !     enum_i_vals=(/almo_domain_layout_orbital,almo_domain_layout_atomic,almo_domain_layout_molecular/))
    !CALL section_add_keyword(section,keyword)
    !CALL keyword_release(keyword)

    !CALL keyword_create(keyword, name="DOMAIN_LAYOUT_AOS",&
    !     description="Atomic orbitals or groups of atomic orbitals represent domains over which electrons "//&
    !     "can be delocalized. This keyword imposes constraints on the structure of domains",&
    !     usage="DOMAIN_LAYOUT_MOS MOLECULAR",&
    !     default_i_val=almo_domain_layout_molecular,&
    !     enum_c_vals=s2a("ATOMIC", "MOLECULAR"),&
    !     enum_desc=s2a("Atomic blocks represent domains. That is, if a basis function on an atom is in"//&
    !     " domain A then all basis functions on this atom are in domain A.",&
    !     "Molecular subsets represent domains. That is, if a basis function on a molecule is"//&
    !     " in domain A then all basis functions on this molecule are in domain A. "//&
    !     "This is the only implemented option"),&
    !     enum_i_vals=(/almo_domain_layout_atomic,almo_domain_layout_molecular/))
    !CALL section_add_keyword(section,keyword)
    !CALL keyword_release(keyword)

    !CALL keyword_create(keyword, name="MATRIX_CLUSTERING_MOS",&
    !     description="Blocks of matrices in the MO basis set are distributed for parallel computations. "//&
    !     "This keywords specifies the type of matrix blocks.",&
    !     usage="MATRIX_CLUSTERING_MOS MOLECULAR",&
    !     default_i_val=almo_mat_distr_molecular,&
    !     enum_c_vals=s2a("ATOMIC", "MOLECULAR"),&
    !     enum_desc=s2a("Not recommended. ZZZ Maybe used for the PAO-based methods in the future",&
    !                   "All molecular orbitals of a molecule belong to the same block."),&
    !     enum_i_vals=(/almo_mat_distr_atomic,almo_mat_distr_molecular/))
    !CALL section_add_keyword(section,keyword)
    !CALL keyword_release(keyword)

    !CALL keyword_create(keyword, name="MATRIX_CLUSTERING_AOS",&
    !     description="Blocks of matrices in the AO basis set are distributed for parallel computations."//&
    !     " This keywords specifies the type of matrix blocks.",&
    !     usage="MATRIX_CLUSTERING_AOS MOLECULAR",&
    !     default_i_val=almo_mat_distr_molecular,&
    !     enum_c_vals=s2a("ATOMIC", "MOLECULAR"),&
    !     enum_desc=s2a("All atomic orbitals of an atom belong to the "//&
    !     "same block. Use only if there are very large molecules in the system. "//&
    !     "ZZZ Maybe used for the PAO-based methods in the future",&
    !     "All atomic orbitals of a molecule belong to the same block."),&
    !     enum_i_vals=(/almo_mat_distr_atomic,almo_mat_distr_molecular/))
    !CALL section_add_keyword(section,keyword)
    !CALL keyword_release(keyword)

    !CALL keyword_create(keyword, name="CONSTRAINT_TYPE",&
    !     description="Determines the type of ALMO constraints",&
    !     usage="CONSTRAINT_TYPE BLOCK_DIAGONAL",&
    !     default_i_val=almo_constraint_block_diagonal,&
    !     enum_c_vals=s2a("BLOCK_DIAGONAL", "DISTANCE", "AO_OVERLAP"),&
    !     enum_desc=s2a("MO coefficient matrix is block-diagonal",&
    !     "MO coefficients are quenched according to the distance criterion",&
    !     "MO coefficients are quenched according to the AO overlap criterion"),&
    !     enum_i_vals=(/almo_constraint_block_diagonal,almo_constraint_distance,&
    !     almo_constraint_ao_overlap/))
    !CALL section_add_keyword(section,keyword)
    !CALL keyword_release(keyword)
 
    !!CALL keyword_create(keyword, name="OUTER_MAX_SCF_Q",&
    !!     description="Maximum number of the outer loop SCF iterations for optimization of quenched ALMOs",&
    !!     usage="OUTER_MAX_SCF_Q 10", default_i_val=0)
    !!CALL section_add_keyword(section,keyword)
    !!CALL keyword_release(keyword)

    !!CALL keyword_create(keyword, name="EPS_USE_PREV_AS_GUESS",&
    !!     description="SCF convergence below which quantities from previous iterations"//&
    !!                 " can be used as initial guess for the current iteration.",&
    !!     usage="EPS_USE_PREV_AS_GUESS 0.01", default_r_val=0.001_dp)
    !!CALL section_add_keyword(section,keyword)
    !!CALL keyword_release(keyword)

    !!CALL keyword_create(keyword, name="MIXING_FRACTION",&
    !!     description="Weight of the new KS matrix in the mixing during the SCF procedure.",&
    !!     usage="MIXING_FRACTION 0.45", default_r_val=0.45_dp)
    !!CALL section_add_keyword(section,keyword)
    !!CALL keyword_release(keyword)

    !!CALL keyword_create(keyword, name="FIXED_MU",&
    !!     description="Fix chemical potential or optimize it to get "//&
    !!                 "the correct number of electrons",&
    !!     usage="FIXED_MU .TRUE.", default_l_val=.FALSE.,&
    !!     lone_keyword_l_val=.TRUE.)
    !!CALL section_add_keyword(section,keyword)
    !!CALL keyword_release(keyword)

    !!CALL keyword_create(keyword, name="MU",&
    !!     description="Value (or initial guess) for the chemical potential."//&
    !!                 " Provide energy between HOMO and LUMO energy.",&
    !!     usage="MU 0.0", default_r_val=-0.1_dp)
    !!CALL section_add_keyword(section,keyword)
    !!CALL keyword_release(keyword)

    !CALL keyword_create(keyword, name="XALMO_ALGORITHM",&
    !     description="Specifies the algorithm to update ALMOs on eXtended domains (XALMOs).",&
    !     usage="XALMO_ALGORITHM PCG",&
    !     default_i_val=almo_scf_diag,&
    !     enum_c_vals=s2a("DDIAG", "PCG"),&
    !     enum_desc=s2a("Domain diagonalization",&
    !                   "Preconditioned conjugate gradient"),&
    !     enum_i_vals=(/almo_scf_diag,almo_scf_pcg/))
    !CALL section_add_keyword(section,keyword)
    !CALL keyword_release(keyword)

    !CALL keyword_create(keyword, name="QUENCHER_RADIUS_TYPE",&
    !     description="Determines the type of atomic radii used for imposing the ALMO constraints",&
    !     usage="QUENCHER_RADIUS_TYPE VDW",&
    !     default_i_val=do_bondparm_vdw,&
    !     enum_c_vals=s2a("COVALENT", "VDW"),&
    !     enum_desc=s2a("Covalent atomic radii",&
    !     "Van der Waals atomic radii"),&
    !     enum_i_vals=(/do_bondparm_covalent,do_bondparm_vdw/))
    !CALL section_add_keyword(section,keyword)
    !CALL keyword_release(keyword)
 
    !CALL keyword_create(keyword, name="QUENCHER_R0_FACTOR",&
    !     description="Parameter to calculate the inner soft cutoff radius: "//&
    !     !"r0 = r0_factor*(radius(at1)+radius(at2)) + r0_shift",&
    !     "r0 = r0_factor*(radius(at1)+radius(at2))",&
    !     usage="QUENCHER_R0_FACTOR 1.05", default_r_val=1.05_dp)
    !CALL section_add_keyword(section,keyword)
    !CALL keyword_release(keyword)

    !!CALL keyword_create(keyword, name="QUENCHER_R0_SHIFT",&
    !!     description="Parameter to calculate the inner soft cutoff radius (in Angstrom): "//&
    !!     "r0 = r0_factor*(radius(at1)+radius(at2)) + r0_shift",&
    !!     usage="QUENCHER_R0_SHIFT 0.0", default_r_val=0.0_dp)
    !!     
    !!CALL section_add_keyword(section,keyword)
    !!CALL keyword_release(keyword)

    !CALL keyword_create(keyword, name="QUENCHER_R1_FACTOR",&
    !     description="Parameter to calculate the outer soft cutoff radius: "//&
    !     !"r1 = r1_factor*(radius(at1)+radius(at2)) + r1_shift",&
    !     "r1 = r1_factor*(radius(at1)+radius(at2))",&
    !     usage="QUENCHER_R1_FACTOR 1.55", default_r_val=1.55_dp)
    !CALL section_add_keyword(section,keyword)
    !CALL keyword_release(keyword)

    !!CALL keyword_create(keyword, name="QUENCHER_R1_SHIFT",&
    !!     description="Parameter to calculate the outer soft cutoff radius (in Angstrom): "//&
    !!     "r1 = r1_factor*(radius(at1)+radius(at2)) + r1_shift",&
    !!     usage="QUENCHER_R1_SHIFT 0.0", default_r_val=0.0_dp)
    !!     
    !!CALL section_add_keyword(section,keyword)
    !!CALL keyword_release(keyword)

    !CALL keyword_create(keyword, name="QUENCHER_AO_OVERLAP_0",&
    !     description="Overlap value of the inner soft cutoff",&
    !     usage="QUENCHER_AO_OVERLAP_0 1.0E-4", default_r_val=1.0E-4_dp)
    !CALL section_add_keyword(section,keyword)
    !CALL keyword_release(keyword)

    !CALL keyword_create(keyword, name="QUENCHER_AO_OVERLAP_1",&
    !     description="Overlap value of the outer soft cutoff",&
    !     usage="QUENCHER_AO_OVERLAP_1 1.0E-6", default_r_val=1.0E-6_dp)
    !CALL section_add_keyword(section,keyword)
    !CALL keyword_release(keyword)

    !CALL keyword_create(keyword, name="ENVELOPE_AMPLITUDE",&
    !     description="Defines an upper bound on the maximum norm of the MO coefficients",&
    !     usage="ENVELOPE_AMPLITUDE 1.0", default_r_val=1.0_dp)
    !CALL section_add_keyword(section,keyword)
    !CALL keyword_release(keyword)

    !CALL keyword_create(keyword, name="DELOC_CAYLEY_TENSOR_TYPE",&
    !     description="Tensor properties of occupied and virtual indices",&
    !     usage="DELOC_CAYLEY_TENSOR_TYPE ORTHOGONAL",&
    !     default_i_val=tensor_orthogonal,&
    !     enum_c_vals=s2a("ORTHOGONAL", "BIORTHOGONAL"),&
    !     enum_desc=s2a("Orthogonalize both occupied and virtual orbitals",&
    !     "Contravariant virtual (MOs or AOs) and covariant occupied orbitals"),&
    !     enum_i_vals=(/tensor_orthogonal,tensor_up_down/))
    !CALL section_add_keyword(section,keyword)
    !CALL keyword_release(keyword)
 
    !CALL keyword_create(keyword, name="DELOC_CAYLEY_CONJUGATOR",&
    !     description="Various methods to compute step directions in the CG algorithm",&
    !     usage="DELOC_CAYLEY_CONJUGATOR POLAK_RIBIERE",&
    !     default_i_val=cg_hager_zhang,&
    !     enum_c_vals=s2a("ZERO", "POLAK_RIBIERE", "FLETCHER_REEVES",&
    !     "HESTENES_STIEFEL", "FLETCHER", "LIU_STOREY", "DAI_YUAN","HAGER_ZHANG"),&
    !     enum_desc=s2a("Steepest descent","Polak and Ribiere",&
    !     "Fletcher and Reeves","Hestenes and Stiefel",&
    !     "Fletcher (Conjugate descent)","Liu and Storey",&
    !     "Dai and Yuan","Hager and Zhang"),&
    !     enum_i_vals=(/cg_zero,cg_polak_ribiere,cg_fletcher_reeves,&
    !                   cg_hestenes_stiefel,cg_fletcher,cg_liu_storey,&
    !                   cg_dai_yuan,cg_hager_zhang/))
    !CALL section_add_keyword(section,keyword)
    !CALL keyword_release(keyword)

    !CALL keyword_create(keyword, name="DELOC_CAYLEY_MAX_ITER",&
    !     description="Maximum number of CG iterations to solve Ricatti equations",&
    !     usage="DELOC_CAYLEY_MAX_ITER 100",default_i_val=50)
    !CALL section_add_keyword(section,keyword)
    !CALL keyword_release(keyword)
 
    !CALL keyword_create(keyword, name="DELOC_CAYLEY_EPS_CONVERGENCE",&
    !     description="Convergence criterion of the CG algorithm",&
    !     usage="DELOC_CAYLEY_EPS_CONVERGENCE 1.e-6", default_r_val=1.e-7_dp)
    !CALL section_add_keyword(section,keyword)
    !CALL keyword_release(keyword)

    !CALL keyword_create(keyword, name="DELOC_CAYLEY_VIR_PRECOND",&
    !     description="Use preconditioner for the virtual subspace",&
    !     usage="DELOC_CAYLEY_VIR_PRECOND .TRUE.", default_l_val=.TRUE.,&
    !     lone_keyword_l_val=.TRUE.)
    !CALL section_add_keyword(section,keyword)
    !CALL keyword_release(keyword)

    !CALL keyword_create(keyword, name="DELOC_CAYLEY_OCC_PRECOND",&
    !     description="Use preconditioner for the occupied subspace",&
    !     usage="DELOC_CAYLEY_OCC_PRECOND .TRUE.", default_l_val=.TRUE.,&
    !     lone_keyword_l_val=.TRUE.)
    !CALL section_add_keyword(section,keyword)
    !CALL keyword_release(keyword)

    !CALL keyword_create(keyword, name="DELOC_TRUNCATE_VIRTUALS",&
    !     description="Truncation of the virtual subspace",&
    !     usage="DELOC_TRUNCATE_VIRTUALS MINIMAL",&
    !     default_i_val=virt_full,&
    !     enum_c_vals=s2a("FULL", "MINIMAL","OCC_SIZE", "EXACT_NUMBER_PER_DOMAIN"),&
    !     enum_desc=s2a("Keep all virtual orbitals","Retained virtuals "//&
    !     "complement occupied orbitals to form the minimal basis set",&
    !     "Number of virtuals is equal to the number of occupied orbitals",&
    !     "Specify exact number of virtuals per domain with DELOC_VIRT_PER_DOMAIN"),&
    !     enum_i_vals=(/virt_full,virt_minimal,virt_occ_size,&
    !                   virt_number/))
    !CALL section_add_keyword(section,keyword)
    !CALL keyword_release(keyword)

    !CALL keyword_create(keyword, name="DELOC_VIRT_PER_DOMAIN",&
    !     description="Number of virtual orbitals (per domain, atom or molecule) "//&
    !     "retained to obtain the delocalization correction",&
    !     usage="DELOC_VIRT_PER_DOMAIN",default_i_val=-1)
    !CALL section_add_keyword(section,keyword)
    !CALL keyword_release(keyword)
 
    !CALL keyword_create(keyword, name="DELOC_USE_OCC_ORBS",&
    !     description="Use occupied orbitals (as opposed to density matrix) "//&
    !     "to calculate correction for electron delocalization",&
    !     usage="DELOC_USE_OCC_ORBS .TRUE.", default_l_val=.TRUE.,&
    !     lone_keyword_l_val=.TRUE.)
    !CALL section_add_keyword(section,keyword)
    !CALL keyword_release(keyword)

    !CALL keyword_create(keyword, name="DELOC_CAYLEY_USE_VIRT_ORBS",&
    !     description="Use virtual orbitals (as opposed to the 1-P projector) "//&
    !     "to calculate correction for electron delocalization. Works only if "//&
    !     "DELOC_USE_OCC_ORBS is set to TRUE",&
    !     usage="DELOC_CAYLEY_USE_VIRT_ORBS .TRUE.", default_l_val=.FALSE.,&
    !     lone_keyword_l_val=.TRUE.)
    !CALL section_add_keyword(section,keyword)
    !CALL keyword_release(keyword)

    !CALL keyword_create(keyword, name="DELOC_CAYLEY_LINEAR",&
    !     description="Neglect the quadratic term in the Riccati equations. "//&
    !     "Equivalent to the first order correction to the occupied orbitals "//&
    !     "(second order correction to the energy)",&
    !     usage="DELOC_CAYLEY_LINEAR .FALSE.", default_l_val=.FALSE.,&
    !     lone_keyword_l_val=.TRUE.)
    !CALL section_add_keyword(section,keyword)
    !CALL keyword_release(keyword)

    !CALL keyword_create(keyword, name="OPT_K_OUTER_MAX_ITER",&
    !     description="Maximum number of outer loop iterations to optimize retained virtual orbitals",&
    !     usage="OPT_K_OUTER_MAX_ITER 10",default_i_val=1)
    !CALL section_add_keyword(section,keyword)
    !CALL keyword_release(keyword)
 
    !CALL keyword_create(keyword, name="OPT_K_MAX_ITER",&
    !     description="Maximum number of iterations to optimize retained virtual orbitals",&
    !     usage="OPT_K_MAX_ITER 100",default_i_val=100)
    !CALL section_add_keyword(section,keyword)
    !CALL keyword_release(keyword)
 
    !CALL keyword_create(keyword, name="OPT_K_EPS_CONVERGENCE",&
    !     description="Convergence criterion of the optimization algorithm",&
    !     usage="OPT_K_EPS_CONVERGENCE 1.e-5", default_r_val=1.e-5_dp)
    !CALL section_add_keyword(section,keyword)
    !CALL keyword_release(keyword)

    !CALL keyword_create(keyword, name="OPT_K_TRIAL_STEP_SIZE",&
    !     description="Size of the trial step along the gradient",&
    !     usage="OPT_K_TRIAL_STEP_SIZE 0.05", default_r_val=0.05_dp)
    !CALL section_add_keyword(section,keyword)
    !CALL keyword_release(keyword)

    !CALL keyword_create(keyword, name="OPT_K_TRIAL_STEP_SIZE_MULTIPLIER",&
    !     description="The trial step size is obtained by multiplying the optimal step size "//&
    !     "from the previous iteration",&
    !     usage="OPT_K_TRIAL_STEP_SIZE_multiplier 1.0", default_r_val=1.4_dp)
    !CALL section_add_keyword(section,keyword)
    !CALL keyword_release(keyword)

    !CALL keyword_create(keyword, name="OPT_K_CONJ_ITER_START",&
    !     description="Iteration for switching from the steepest descent algorithm "//&
    !     "to conjugate gradient",&
    !     usage="OPT_K_CONJ_ITER_START 5",default_i_val=0)
    !CALL section_add_keyword(section,keyword)
    !CALL keyword_release(keyword)
 
    !CALL keyword_create(keyword, name="OPT_K_CONJ_ITER_FREQ_RESET",&
    !     description="Reset frequency of the conjugate gradient direction",&
    !     usage="OPT_K_CONJ_ITER_FREQ_RESET 20",default_i_val=1000000)
    !CALL section_add_keyword(section,keyword)
    !CALL keyword_release(keyword)
 
    !CALL keyword_create(keyword, name="OPT_K_CONJUGATOR",&
    !     description="Various methods to compute step directions in the CG algorithm",&
    !     usage="OPT_K_CONJUGATOR POLAK_RIBIERE",&
    !     default_i_val=cg_hager_zhang,&
    !     enum_c_vals=s2a("ZERO", "POLAK_RIBIERE", "FLETCHER_REEVES",&
    !     "HESTENES_STIEFEL", "FLETCHER", "LIU_STOREY", "DAI_YUAN","HAGER_ZHANG"),&
    !     enum_desc=s2a("Steepest descent","Polak and Ribiere",&
    !     "Fletcher and Reeves","Hestenes and Stiefel",&
    !     "Fletcher (Conjugate descent)","Liu and Storey",&
    !     "Dai and Yuan","Hager and Zhang"),&
    !     enum_i_vals=(/cg_zero,cg_polak_ribiere,cg_fletcher_reeves,&
    !                   cg_hestenes_stiefel,cg_fletcher,cg_liu_storey,&
    !                   cg_dai_yuan,cg_hager_zhang/))
    !CALL section_add_keyword(section,keyword)
    !CALL keyword_release(keyword)

    !CALL keyword_create(keyword, name="OPT_K_PREC_ITER_START",&
    !     description="Start using the preconditioner (approximate preconditioners "//&
    !     "might not be valid on early iterations)",&
    !     usage="OPT_K_PREC_ITER_START 2",default_i_val=0)
    !CALL section_add_keyword(section,keyword)
    !CALL keyword_release(keyword)
 
    !CALL keyword_create(keyword, name="OPT_K_PREC_ITER_FREQ_UPDATE",&
    !     description="Frequency for updating the preconditioner",&
    !     usage="OPT_K_PREC_ITER_FREQ_UPDATE 10",default_i_val=1)
    !CALL section_add_keyword(section,keyword)
    !CALL keyword_release(keyword)
 
    NULLIFY(subsection)
    CALL create_optimizer_section(subsection,optimizer_block_diagonal_diis)
    CALL section_add_subsection(section, subsection)
    CALL section_release(subsection)

    NULLIFY(subsection)
    CALL create_optimizer_section(subsection,optimizer_block_diagonal_pcg)
    CALL section_add_subsection(section, subsection)
    CALL section_release(subsection)

    NULLIFY(subsection)
    CALL create_optimizer_section(subsection,optimizer_xalmo_pcg)
    CALL section_add_subsection(section, subsection)
    CALL section_release(subsection)

    !NULLIFY(subsection)
    !CALL create_almo_eda_section(subsection)
    !CALL section_add_subsection(section, subsection)
    !CALL section_release(subsection)


  END SUBROUTINE create_almo_scf_section

! *****************************************************************************
!> \brief The ALMO EDA section controls decomposition analysis based on ALMOs
!> \param section ...
!> \par History
!>       2014.10 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
SUBROUTINE create_almo_eda_section(section)

    TYPE(section_type), POINTER              :: section

    CHARACTER(len=*), PARAMETER :: routineN = 'create_almo_eda_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,failure)
   
    CALL section_create(section,"ALMO_DA",&
         description="Controls decomposition analysis based on ALMOs and XALMOs. "//&
         "Not yet fully implemented.",&
         n_keywords=1, n_subsections=0, repeats=.FALSE.,&
         citations=(/Khaliullin2007,Khaliullin2008/))

    NULLIFY (keyword)

    CALL keyword_create(keyword, name="FRZ_TERM",&
         description="Perform calculations on single molecules to compute the frozen density term",&
         usage="FRZ_TERM ISOLATED", default_i_val=0,&
         !enum_c_vals=s2a("SKIP", "ISOLATED", "CRYSTAL"),&
         enum_c_vals=s2a("SKIP", "CRYSTAL"),&
         enum_desc=s2a("Do not compute the frozen energy term.",&
         !"Use isolated molecules as the reference.",&
         "Use molecules in their positions in the crystal cell as the reference. "//&
         !"Substantially more expensive than ISOLATED. "//&
         "Be careful interpreting this term for systems with charged fragmetns."),&
         enum_i_vals=(/almo_frz_none,&
         !              almo_frz_isolated,&
                       almo_frz_crystal/))
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)


  END SUBROUTINE create_almo_eda_section

! *****************************************************************************
!> \brief The optimizer section is a collection of keywords that are similar
!>        to all optimization methods (e.g. target error, number of iterations)
!> \param section ...
!> \param optimizer_id - allows to adapt the standard section for specific needs 
!> \par History
!>       2012.03 created [Rustam Z Khaliullin]
!>       2014.10 fully integrated [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
SUBROUTINE create_optimizer_section(section,optimizer_id)

    TYPE(section_type), POINTER              :: section
    INTEGER, INTENT(IN)                      :: optimizer_id

    CHARACTER(len=*), PARAMETER :: routineN = 'create_optimizer_section', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: optimizer_type
    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,failure)
   
    ! choose the name of the section
    SELECT CASE(optimizer_id)
    CASE(optimizer_block_diagonal_diis)
       CALL section_create(section,"ALMO_OPTIMIZER_DIIS",&
            description="Controls the iterative DIIS-accelerated optimization of block-diagonal ALMOs.",&
            n_keywords=3, n_subsections=0, repeats=.FALSE.)
       optimizer_type=optimizer_diis
    CASE(optimizer_block_diagonal_pcg)
       CALL section_create(section,"ALMO_OPTIMIZER_PCG",&
            description="Controls the PCG optimization of block-diagonal ALMOs.",&
            n_keywords=6, n_subsections=0, repeats=.FALSE.)
       optimizer_type=optimizer_pcg
    CASE(optimizer_xalmo_pcg)
       CALL section_create(section,"XALMO_OPTIMIZER_PCG",&
            description="Controls the PCG optimization of extended ALMOs.",&
            n_keywords=6, n_subsections=0, repeats=.FALSE.)
       optimizer_type=optimizer_pcg
    CASE DEFAULT
       CPErrorMessage(cp_failure_level,routineP,"No default values allowed")
       CPPrecondition(.FALSE.,cp_failure_level,routineP,failure)
    END SELECT

    NULLIFY (keyword)

    ! add common keywords
    CALL keyword_create(keyword, name="MAX_ITER",&
         description="Maximum number of iterations",&
         usage="MAX_ITER 100", default_i_val=20)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="EPS_ERROR",&
         description="Target value of the MAX norm of the error",&
         usage="EPS_ERROR 1.E-6", default_r_val=1.0E-5_dp)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    ! add keywords specific to each type
    SELECT CASE(optimizer_type)
    CASE(optimizer_diis)

       CALL keyword_create(keyword, name="N_DIIS",&
            description="Number of error vectors to be used in the DIIS "//&
            "optimization procedure",&
            usage="N_DIIS 5", default_i_val=6)
       CALL section_add_keyword(section,keyword)
       CALL keyword_release(keyword)

    CASE(optimizer_pcg)
       
       CALL keyword_create(keyword, name="LIN_SEARCH_EPS_ERROR",&
            description="Target value of the gradient norm during the linear search",&
            usage="LIN_SEARCH_EPS_ERROR 1.E-2", default_r_val=1.0E-3_dp)
       CALL section_add_keyword(section,keyword)
       CALL keyword_release(keyword)

       CALL keyword_create(keyword, name="LIN_SEARCH_STEP_SIZE_GUESS",&
            description="The size of the first step in the linear search",&
            usage="LIN_SEARCH_STEP_SIZE_GUESS 0.1", default_r_val=1.0_dp)
       CALL section_add_keyword(section,keyword)
       CALL keyword_release(keyword)

       CALL keyword_create(keyword, name="MAX_ITER_OUTER_LOOP",&
            description="Maximum number of iterations in the outer loop. "//&
            "Use the outer loop to update the preconditioner and reset the conjugator. "//&
            "This can speed up convergence significantly.",&
            usage="MAX_ITER 10", default_i_val=0)
       CALL section_add_keyword(section,keyword)
       CALL keyword_release(keyword)

       CALL keyword_create(keyword, name="PRECONDITIONER",&
            description="Select a preconditioner for the conjugate gradient optimization",&
            usage="PRECONDITIONER NONE",&
            default_i_val=-1,&
            enum_c_vals=s2a("DEFAULT", "NONE"),&
            enum_desc=s2a("Default preconditioner","Do not use preconditioner"),&
            enum_i_vals=(/prec_default,prec_zero/))
       CALL section_add_keyword(section,keyword)
       CALL keyword_release(keyword)

       CALL keyword_create(keyword, name="CONJUGATOR",&
            description="Various methods to compute step directions in the PCG optimization",&
            usage="CONJUGATOR POLAK_RIBIERE",&
            default_i_val=cg_hager_zhang,&
            enum_c_vals=s2a("ZERO", "POLAK_RIBIERE", "FLETCHER_REEVES",&
            "HESTENES_STIEFEL", "FLETCHER", "LIU_STOREY", "DAI_YUAN","HAGER_ZHANG"),&
            enum_desc=s2a("Steepest descent","Polak and Ribiere",&
            "Fletcher and Reeves","Hestenes and Stiefel",&
            "Fletcher (Conjugate descent)","Liu and Storey",&
            "Dai and Yuan","Hager and Zhang"),&
            enum_i_vals=(/cg_zero,cg_polak_ribiere,cg_fletcher_reeves,&
                          cg_hestenes_stiefel,cg_fletcher,cg_liu_storey,&
                          cg_dai_yuan,cg_hager_zhang/))
       CALL section_add_keyword(section,keyword)
       CALL keyword_release(keyword)

    END SELECT 


  END SUBROUTINE create_optimizer_section

!!! *****************************************************************************
!!!> \brief The developer section allows uninterrupted development of the code
!!!>        by keeping untested keywords neatly separated from those publicly
!!!>        available
!!!> \param section ...
!!!> \param error ...
!!!> \par History
!!!>       2014.10 created [Rustam Z Khaliullin]
!!!> \author Rustam Z Khaliullin
!!! *****************************************************************************
!!  SUBROUTINE create_developer_section(section)
!!
!!    TYPE(section_type), POINTER              :: section
!!
!!    CHARACTER(len=*), PARAMETER :: routineN = 'create_developer_section', &
!!      routineP = moduleN//':'//routineN
!!
!!    LOGICAL                                  :: failure
!!    TYPE(keyword_type), POINTER              :: keyword
!!
!!    failure=.FALSE.
!!
!!    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,failure)
!!    IF (.NOT. failure) THEN
!!      
!!       CALL section_create(section,"DEVELOPER",&
!!            description="Developer section for all ALMO-based methods. "//&
!!                        "Allows uninterrupted development of the code "//&
!!                        "by keeping untested keywords neatly separated from those publicly available",&
!!            n_keywords=3, n_subsections=0, repeats=.FALSE.)
!!             
!!
!!       NULLIFY (keyword)
!!   
!!
!!    END IF
!!
!!  END SUBROUTINE create_developer_section

END MODULE input_cp2k_almo
