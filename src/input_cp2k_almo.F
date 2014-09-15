!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief input for the ALMO SCF section
!> \author Rustam Khaliullin
! *****************************************************************************
MODULE input_cp2k_almo
  USE bibliography,                    ONLY: Khaliullin2013
  USE f77_blas
  USE input_constants,                 ONLY: &
       almo_deloc_full_scf, almo_deloc_none, almo_deloc_qscf, almo_deloc_qx, &
       almo_deloc_x, almo_deloc_x_then_full_scf, almo_scf_diag, almo_scf_pcg, &
       cg_dai_yuan, cg_fletcher, cg_fletcher_reeves, cg_hager_zhang, &
       cg_hestenes_stiefel, cg_liu_storey, cg_polak_ribiere, cg_zero
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_create,&
                                             section_type
  USE kinds,                           ONLY: dp
  USE string_utilities,                ONLY: s2a
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_almo'

  PUBLIC :: create_almo_scf_section

! RZK-MERGE
!  USE cp_output_handling,              ONLY: cp_print_key_section_create
!  USE input_section_types,             ONLY: section_add_keyword,&
!                                             section_add_subsection,&
!                                             section_create,&
!                                             section_release,&
!                                             section_type
!  USE input_val_types,                 ONLY: integer_t,&
!                                             real_t

CONTAINS

! *****************************************************************************
!> \brief create the almo scf section
!> \param section ...
!> \param error ...
!> \par History
!>       2011.05 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
SUBROUTINE create_almo_scf_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_almo_scf_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

!TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"ALMO_SCF",&
            description="Settings for methods based on absolutely localized orbitals"//&
            " (ALMOs). Restricted to closed-shell molecular systems.",&
            n_keywords=11, n_subsections=0, repeats=.FALSE.,&
            citations=(/Khaliullin2013/),&
            required=.FALSE., error=error)

       NULLIFY (keyword)

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
       !     enum_i_vals=(/almo_domain_layout_orbital,almo_domain_layout_atomic,almo_domain_layout_molecular/),&
       !     error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

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
       !     enum_i_vals=(/almo_domain_layout_atomic,almo_domain_layout_molecular/),&
       !     error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="MATRIX_CLUSTERING_MOS",&
       !     description="Blocks of matrices in the MO basis set are distributed for parallel computations. "//&
       !     "This keywords specifies the type of matrix blocks.",&
       !     usage="MATRIX_CLUSTERING_MOS MOLECULAR",&
       !     default_i_val=almo_mat_distr_molecular,&
       !     enum_c_vals=s2a("ATOMIC", "MOLECULAR"),&
       !     enum_desc=s2a("Not recommended. ZZZ Maybe used for the PAO-based methods in the future",&
       !                   "All molecular orbitals of a molecule belong to the same block."),&
       !     enum_i_vals=(/almo_mat_distr_atomic,almo_mat_distr_molecular/),&
       !     error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

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
       !     enum_i_vals=(/almo_mat_distr_atomic,almo_mat_distr_molecular/),&
       !     error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="CONSTRAINT_TYPE",&
       !     description="Determines the type of ALMO constraints",&
       !     usage="CONSTRAINT_TYPE BLOCK_DIAGONAL",&
       !     default_i_val=almo_constraint_block_diagonal,&
       !     enum_c_vals=s2a("BLOCK_DIAGONAL", "DISTANCE", "AO_OVERLAP"),&
       !     enum_desc=s2a("MO coefficient matrix is block-diagonal",&
       !     "MO coefficients are quenched according to the distance criterion",&
       !     "MO coefficients are quenched according to the AO overlap criterion"),&
       !     enum_i_vals=(/almo_constraint_block_diagonal,almo_constraint_distance,&
       !     almo_constraint_ao_overlap/),&
       !     error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)
    
       CALL keyword_create(keyword, name="BLOCKED_MAX_ITER",&
            description="Maximum number of SCF iterations for optimization of block-diagonal ALMOs",&
            usage="BLOCKED_MAX_ITER 200", default_i_val=100,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="BLOCKED_N_DIIS",&
            description="Number of error vectors to be used in the DIIS "//&
            "optimization procedure of block-diagonal ALMOs",&
            usage="BLOCKED_N_DIIS 5", default_i_val=7,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="BLOCKED_EPS_ITER",&
            description="Target accuracy for the convergence of optimization of block-diagonal ALMOs",&
            usage="BLOCKED_EPS_ITER 1.e-5", default_r_val=1.e-5_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DELOCALIZE_MAX_ITER",&
            description="Maximum number of SCF iterations for optimization of ALMOs",&
            usage="DELOCALIZE_MAX_ITER 200", default_i_val=100,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="OUTER_MAX_SCF_Q",&
       !     description="Maximum number of the outer loop SCF iterations for optimization of quenched ALMOs",&
       !     usage="OUTER_MAX_SCF_Q 10", default_i_val=0,error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DELOCALIZE_EPS_ITER",&
            description="Target accuracy for the convergence of optimization of ALMOs",&
            usage="DELOCALIZE_EPS_ITER 1.e-5", default_r_val=1.e-5_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DELOCALIZE_EPS_LIN_SEARCH",&
            description="Target accuracy of the line search for the PCG optimization of ALMOs",&
            usage="DELOCALIZE_EPS_LIN_SEARCH 1.e-6", default_r_val=1.e-7_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_FILTER",&
            description="Threshold for the sparsity filter",&
            usage="EPS_FILTER 1.e-6", default_r_val=1.e-7_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       
       CALL keyword_create(keyword, name="DELOCALIZE_CONJUGATOR",&
            description="Various methods to compute step directions in the PCG optimization of ALMOs",&
            usage="DELOCALIZE_CONJUGATOR POLAK_RIBIERE",&
            default_i_val=cg_hager_zhang,&
            enum_c_vals=s2a("ZERO", "POLAK_RIBIERE", "FLETCHER_REEVES",&
            "HESTENES_STIEFEL", "FLETCHER", "LIU_STOREY", "DAI_YUAN","HAGER_ZHANG"),&
            enum_desc=s2a("Steepest descent","Polak and Ribiere",&
            "Fletcher and Reeves","Hestenes and Stiefel",&
            "Fletcher (Conjugate descent)","Liu and Storey",&
            "Dai and Yuan","Hager and Zhang"),&
            enum_i_vals=(/cg_zero,cg_polak_ribiere,cg_fletcher_reeves,&
                          cg_hestenes_stiefel,cg_fletcher,cg_liu_storey,&
                          cg_dai_yuan,cg_hager_zhang/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="EPS_USE_PREV_AS_GUESS",&
       !     description="SCF convergence below which quantities from previous iterations"//&
       !                 " can be used as initial guess for the current iteration.",&
       !     usage="EPS_USE_PREV_AS_GUESS 0.01", default_r_val=0.001_dp,error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="MIXING_FRACTION",&
       !     description="Weight of the new KS matrix in the mixing during the SCF procedure.",&
       !     usage="MIXING_FRACTION 0.45", default_r_val=0.45_dp,error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="FIXED_MU",&
       !     description="Fix chemical potential or optimize it to get "//&
       !                 "the correct number of electrons",&
       !     usage="FIXED_MU .TRUE.", default_l_val=.FALSE.,&
       !     lone_keyword_l_val=.TRUE.,  error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="MU",&
       !     description="Value (or initial guess) for the chemical potential."//&
       !                 " Provide energy between HOMO and LUMO energy.",&
       !     usage="MU 0.0", default_r_val=-0.1_dp,error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="ALMO_UPDATE_ALGORITHM_BD",&
       !     description="Specifies the algorithm to update block-digonal ALMOs during the SCF procedure.",&
       !     usage="ALMO_UPDATE_ALGORITHM_BD DIAG",&
       !     default_i_val=almo_scf_diag,&
       !     enum_c_vals=s2a("DIAG", "DM_SIGN","PCG"),&
       !     enum_desc=s2a("Straightforward diagonalization. Recommended for systems "//&
       !                   "comprising large number of small fragments.",&
       !                   "Update the density matrix using linear scaling routines. "//&
       !                   "Recommended if large fragments are present.",&
       !                   "Direct minimization using the preconditioned cojugate gradient algorithm. "//&
       !                   "It has not been extensively tested."),&
       !     enum_i_vals=(/almo_scf_diag,almo_scf_dm_sign,almo_scf_pcg/),&
       !     error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DELOCALIZE_ALGORITHM",&
            description="Specifies the algorithm to update ALMOs.",&
            usage="DELOCALIZE_ALGORITHM BDIAG",&
            default_i_val=almo_scf_diag,&
            enum_c_vals=s2a("DDIAG", "PCG"),&
            enum_desc=s2a("Domain diagonalization",&
                          "Preconditioned conjugate gradient"),&
            enum_i_vals=(/almo_scf_diag,almo_scf_pcg/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="QUENCHER_RADIUS_TYPE",&
       !     description="Determines the type of atomic radii used for imposing the ALMO constraints",&
       !     usage="QUENCHER_RADIUS_TYPE VDW",&
       !     default_i_val=do_bondparm_vdw,&
       !     enum_c_vals=s2a("COVALENT", "VDW"),&
       !     enum_desc=s2a("Covalent atomic radii",&
       !     "Van der Waals atomic radii"),&
       !     enum_i_vals=(/do_bondparm_covalent,do_bondparm_vdw/),&
       !     error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)
    
       CALL keyword_create(keyword, name="DELOCALIZE_R_CUTOFF_FACTOR",&
            description="Parameter to calculate the localization radius: "//&
            !"r0 = r0_factor*(radius(at1)+radius(at2)) + r0_shift",&
            "R_cutoff = DELOCALIZE_R_CUTOFF_FACTOR*(radius(at1)+radius(at2))",&
            usage="DELOCALIZE_R_CUTOFF_FACTOR 1.6", default_r_val=1.60_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="QUENCHER_R0_FACTOR",&
       !     description="Parameter to calculate the inner soft cutoff radius: "//&
       !     !"r0 = r0_factor*(radius(at1)+radius(at2)) + r0_shift",&
       !     "r0 = r0_factor*(radius(at1)+radius(at2))",&
       !     usage="QUENCHER_R0_FACTOR 1.05", default_r_val=1.05_dp,&
       !     error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !!CALL keyword_create(keyword, name="QUENCHER_R0_SHIFT",&
       !!     description="Parameter to calculate the inner soft cutoff radius (in Angstrom): "//&
       !!     "r0 = r0_factor*(radius(at1)+radius(at2)) + r0_shift",&
       !!     usage="QUENCHER_R0_SHIFT 0.0", default_r_val=0.0_dp,&
       !!     error=error)
       !!CALL section_add_keyword(section,keyword,error=error)
       !!CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="QUENCHER_R1_FACTOR",&
       !     description="Parameter to calculate the outer soft cutoff radius: "//&
       !     !"r1 = r1_factor*(radius(at1)+radius(at2)) + r1_shift",&
       !     "r1 = r1_factor*(radius(at1)+radius(at2))",&
       !     usage="QUENCHER_R1_FACTOR 1.55", default_r_val=1.55_dp,&
       !     error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !!CALL keyword_create(keyword, name="QUENCHER_R1_SHIFT",&
       !!     description="Parameter to calculate the outer soft cutoff radius (in Angstrom): "//&
       !!     "r1 = r1_factor*(radius(at1)+radius(at2)) + r1_shift",&
       !!     usage="QUENCHER_R1_SHIFT 0.0", default_r_val=0.0_dp,&
       !!     error=error)
       !!CALL section_add_keyword(section,keyword,error=error)
       !!CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="QUENCHER_AO_OVERLAP_0",&
       !     description="Overlap value of the inner soft cutoff",&
       !     usage="QUENCHER_AO_OVERLAP_0 1.0E-4", default_r_val=1.0E-4_dp,&
       !     error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="QUENCHER_AO_OVERLAP_1",&
       !     description="Overlap value of the outer soft cutoff",&
       !     usage="QUENCHER_AO_OVERLAP_1 1.0E-6", default_r_val=1.0E-6_dp,&
       !     error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="ENVELOPE_AMPLITUDE",&
       !     description="Defines an upper bound on the maximum norm of the MO coefficients",&
       !     usage="ENVELOPE_AMPLITUDE 1.0", default_r_val=1.0_dp,&
       !     error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DELOCALIZE_METHOD",&
            description="Methods to treat electron delocalization after block-diagonal ALMOs are computed",&
            usage="DELOCALIZE_METHOD X_R",&
            default_i_val=almo_deloc_qx,&
            enum_c_vals=s2a("NONE","X","SCF","X_R","SCF_R","X_THEN_SCF"),&
            enum_desc=s2a("Neglect electron delocalization",&
            "Single excitation delocalization energy correction",&
            "Fully self-consistent treatment of delocalization",&
            "Spacially restricted single excitation delocalization energy correction",&
            "Spacially restricted self-consistent treatment of delocalization",&
            "Single excitation delocalization energy correction followed by full SCF procedure"),&
            enum_i_vals=(/almo_deloc_none,almo_deloc_x,almo_deloc_full_scf,&
            almo_deloc_qx,almo_deloc_qscf,&
            !almo_deloc_xk,&
            almo_deloc_x_then_full_scf/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="DELOC_CAYLEY_TENSOR_TYPE",&
       !     description="Tensor properties of occupied and virtual indices",&
       !     usage="DELOC_CAYLEY_TENSOR_TYPE ORTHOGONAL",&
       !     default_i_val=tensor_orthogonal,&
       !     enum_c_vals=s2a("ORTHOGONAL", "BIORTHOGONAL"),&
       !     enum_desc=s2a("Orthogonalize both occupied and virtual orbitals",&
       !     "Contravariant virtual (MOs or AOs) and covariant occupied orbitals"),&
       !     enum_i_vals=(/tensor_orthogonal,tensor_up_down/),&
       !     error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)
    
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
       !                   cg_dai_yuan,cg_hager_zhang/),&
       !     error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="DELOC_CAYLEY_MAX_ITER",&
       !     description="Maximum number of CG iterations to solve Ricatti equations",&
       !     usage="DELOC_CAYLEY_MAX_ITER 100",default_i_val=50,error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)
    
       !CALL keyword_create(keyword, name="DELOC_CAYLEY_EPS_CONVERGENCE",&
       !     description="Convergence criterion of the CG algorithm",&
       !     usage="DELOC_CAYLEY_EPS_CONVERGENCE 1.e-6", default_r_val=1.e-7_dp,&
       !     error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="DELOC_CAYLEY_VIR_PRECOND",&
       !     description="Use preconditioner for the virtual subspace",&
       !     usage="DELOC_CAYLEY_VIR_PRECOND .TRUE.", default_l_val=.TRUE.,&
       !     lone_keyword_l_val=.TRUE., error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="DELOC_CAYLEY_OCC_PRECOND",&
       !     description="Use preconditioner for the occupied subspace",&
       !     usage="DELOC_CAYLEY_OCC_PRECOND .TRUE.", default_l_val=.TRUE.,&
       !     lone_keyword_l_val=.TRUE., error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

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
       !                   virt_number/),&
       !     error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="DELOC_VIRT_PER_DOMAIN",&
       !     description="Number of virtual orbitals (per domain, atom or molecule) "//&
       !     "retained to obtain the delocalization correction",&
       !     usage="DELOC_VIRT_PER_DOMAIN",default_i_val=-1,error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)
    
       !CALL keyword_create(keyword, name="DELOC_USE_OCC_ORBS",&
       !     description="Use occupied orbitals (as opposed to density matrix) "//&
       !     "to calculate correction for electron delocalization",&
       !     usage="DELOC_USE_OCC_ORBS .TRUE.", default_l_val=.TRUE.,&
       !     lone_keyword_l_val=.TRUE., error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="DELOC_CAYLEY_USE_VIRT_ORBS",&
       !     description="Use virtual orbitals (as opposed to the 1-P projector) "//&
       !     "to calculate correction for electron delocalization. Works only if "//&
       !     "DELOC_USE_OCC_ORBS is set to TRUE",&
       !     usage="DELOC_CAYLEY_USE_VIRT_ORBS .TRUE.", default_l_val=.FALSE.,&
       !     lone_keyword_l_val=.TRUE.,  error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="DELOC_CAYLEY_LINEAR",&
       !     description="Neglect the quadratic term in the Riccati equations. "//&
       !     "Equivalent to the first order correction to the occupied orbitals "//&
       !     "(second order correction to the energy)",&
       !     usage="DELOC_CAYLEY_LINEAR .FALSE.", default_l_val=.FALSE.,&
       !     lone_keyword_l_val=.TRUE.,  error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="OPT_K_OUTER_MAX_ITER",&
       !     description="Maximum number of outer loop iterations to optimize retained virtual orbitals",&
       !     usage="OPT_K_OUTER_MAX_ITER 10",default_i_val=1,error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)
    
       !CALL keyword_create(keyword, name="OPT_K_MAX_ITER",&
       !     description="Maximum number of iterations to optimize retained virtual orbitals",&
       !     usage="OPT_K_MAX_ITER 100",default_i_val=100,error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)
    
       !CALL keyword_create(keyword, name="OPT_K_EPS_CONVERGENCE",&
       !     description="Convergence criterion of the optimization algorithm",&
       !     usage="OPT_K_EPS_CONVERGENCE 1.e-5", default_r_val=1.e-5_dp,&
       !     error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="OPT_K_TRIAL_STEP_SIZE",&
       !     description="Size of the trial step along the gradient",&
       !     usage="OPT_K_TRIAL_STEP_SIZE 0.05", default_r_val=0.05_dp,&
       !     error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="OPT_K_TRIAL_STEP_SIZE_MULTIPLIER",&
       !     description="The trial step size is obtained by multiplying the optimal step size "//&
       !     "from the previous iteration",&
       !     usage="OPT_K_TRIAL_STEP_SIZE_multiplier 1.0", default_r_val=1.4_dp,&
       !     error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="OPT_K_CONJ_ITER_START",&
       !     description="Iteration for switching from the steepest descent algorithm "//&
       !     "to conjugate gradient",&
       !     usage="OPT_K_CONJ_ITER_START 5",default_i_val=0,error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)
    
       !CALL keyword_create(keyword, name="OPT_K_CONJ_ITER_FREQ_RESET",&
       !     description="Reset frequency of the conjugate gradient direction",&
       !     usage="OPT_K_CONJ_ITER_FREQ_RESET 20",default_i_val=1000000,error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)
    
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
       !                   cg_dai_yuan,cg_hager_zhang/),&
       !     error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="OPT_K_PREC_ITER_START",&
       !     description="Start using the preconditioner (approximate preconditioners "//&
       !     "might not be valid on early iterations)",&
       !     usage="OPT_K_PREC_ITER_START 2",default_i_val=0,error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)
    
       !CALL keyword_create(keyword, name="OPT_K_PREC_ITER_FREQ_UPDATE",&
       !     description="Frequency for updating the preconditioner",&
       !     usage="OPT_K_PREC_ITER_FREQ_UPDATE 10",default_i_val=1,error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)
    
       !CALL keyword_create(keyword, name="INTEGER01",&
       !     description="Development keyword",&
       !     usage="INTEGER01 1",default_i_val=0,error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)
    
       !CALL keyword_create(keyword, name="INTEGER02",&
       !     description="Development keyword",&
       !     usage="INTEGER02 1",default_i_val=0,error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)
    
       !CALL keyword_create(keyword, name="INTEGER03",&
       !     description="Development keyword",&
       !     usage="INTEGER03 1",default_i_val=0,error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)
    
       !CALL keyword_create(keyword, name="INTEGER04",&
       !     description="Development keyword",&
       !     usage="INTEGER04 1",default_i_val=0,error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)
    
       !CALL keyword_create(keyword, name="INTEGER05",&
       !     description="Development keyword",&
       !     usage="INTEGER05 1",default_i_val=0,error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)
    
       !CALL keyword_create(keyword, name="REAL01",&
       !     description="Development keyword",&
       !     usage="REAL01", default_r_val=0.0_dp,&
       !     error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="REAL02",&
       !     description="Development keyword",&
       !     usage="REAL02", default_r_val=0.0_dp,&
       !     error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="REAL03",&
       !     description="Development keyword",&
       !     usage="REAL03", default_r_val=0.0_dp,&
       !     error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="REAL04",&
       !     description="Development keyword",&
       !     usage="REAL04", default_r_val=0.0_dp,&
       !     error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="REAL05",&
       !     description="Development keyword",&
       !     usage="REAL05", default_r_val=0.0_dp,&
       !     error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="LOGICAL01",&
       !     description="Development keyword",&
       !     usage="LOGICAL01 .FALSE.", default_l_val=.FALSE.,&
       !     lone_keyword_l_val=.TRUE.,  error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="LOGICAL02",&
       !     description="Development keyword",&
       !     usage="LOGICAL02 .FALSE.", default_l_val=.FALSE.,&
       !     lone_keyword_l_val=.TRUE.,  error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="LOGICAL03",&
       !     description="Development keyword",&
       !     usage="LOGICAL03 .FALSE.", default_l_val=.FALSE.,&
       !     lone_keyword_l_val=.TRUE.,  error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="LOGICAL04",&
       !     description="Development keyword",&
       !     usage="LOGICAL04 .FALSE.", default_l_val=.FALSE.,&
       !     lone_keyword_l_val=.TRUE.,  error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="LOGICAL05",&
       !     description="Development keyword",&
       !     usage="LOGICAL05 .FALSE.", default_l_val=.FALSE.,&
       !     lone_keyword_l_val=.TRUE.,  error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)
       !
       !CALL keyword_create(keyword, name="INT_LIST01",&
       !     description="Development keyword",&
       !     usage="INT_LIST01 {integer} {integer} .. {integer}", required=.FALSE.,&
       !     n_var=-1, type_of_var=integer_t, error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="INT_LIST02",&
       !     description="Development keyword",&
       !     usage="INT_LIST02 {integer} {integer} .. {integer}", required=.FALSE.,&
       !     n_var=-1, type_of_var=integer_t, error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="INT_LIST03",&
       !     description="Development keyword",&
       !     usage="INT_LIST03 {integer} {integer} .. {integer}", required=.FALSE.,&
       !     n_var=-1, type_of_var=integer_t, error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="INT_LIST04",&
       !     description="Development keyword",&
       !     usage="INT_LIST04 {integer} {integer} .. {integer}", required=.FALSE.,&
       !     n_var=-1, type_of_var=integer_t, error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="INT_LIST05",&
       !     description="Development keyword",&
       !     usage="INT_LIST05 {integer} {integer} .. {integer}", required=.FALSE.,&
       !     n_var=-1, type_of_var=integer_t, error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword,name="REAL_LIST01",&
       !     description="Development keyword",&
       !     usage="REAL_LIST01 {real} {real} .. {real}", required=.FALSE.,&
       !     type_of_var=real_t, n_var=-1, error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword,name="REAL_LIST02",&
       !     description="Development keyword",&
       !     usage="REAL_LIST02 {real} {real} .. {real}", required=.FALSE.,&
       !     type_of_var=real_t, n_var=-1, error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword,name="REAL_LIST03",&
       !     description="Development keyword",&
       !     usage="REAL_LIST03 {real} {real} .. {real}", required=.FALSE.,&
       !     type_of_var=real_t, n_var=-1, error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword,name="REAL_LIST04",&
       !     description="Development keyword",&
       !     usage="REAL_LIST04 {real} {real} .. {real}", required=.FALSE.,&
       !     type_of_var=real_t, n_var=-1, error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword,name="REAL_LIST05",&
       !     description="Development keyword",&
       !     usage="REAL_LIST05 {real} {real} .. {real}", required=.FALSE.,&
       !     type_of_var=real_t, n_var=-1, error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !NULLIFY(subsection)
       !CALL create_optimizer_section(subsection,error)
       !CALL section_add_subsection(section, subsection, error=error)
       !CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_almo_scf_section

! *****************************************************************************
!> \brief create the optimizer section
!> \param section ...
!> \param error ...
!> \par History
!>       2012.03 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
SUBROUTINE create_optimizer_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_optimizer_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       !CALL section_create(section,"ORTHOGONALIZATION_OPT",&
       !     description="Controls iterative orthogonalization of ALMOs.",&
       !     n_keywords=6, n_subsections=0, repeats=.FALSE.,&
       !     required=.FALSE., error=error)

       !NULLIFY (keyword)
   
       !CALL keyword_create(keyword, name="MAX_ITER",&
       !     description="Maximum number of iterations",&
       !     usage="MAX_ITER 100", default_i_val=100,error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="PRECONDITIONER",&
       !     description="Various methods to compute preconditioner during the CG optimization",&
       !     usage="PRECONDITIONER 0",&
       !     default_i_val=0,&
       !     enum_c_vals=s2a("NONE", "SIGMA_ZERO"),&
       !     enum_i_vals=(/0,1/),&
       !     error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="CONJUGATOR",&
       !     description="Various methods to compute conjugated step directions during the CG optimization",&
       !     usage="CONJUGATOR POLAK_RIBIERE",&
       !     default_i_val=cg_hager_zhang,&
       !     enum_c_vals=s2a("ZERO", "POLAK_RIBIERE", "FLETCHER_REEVES",&
       !     "HESTENES_STIEFEL", "FLETCHER", "LIU_STOREY", "DAI_YUAN","HAGER_ZHANG"),&
       !     enum_desc=s2a("Steepest descent","Polak and Ribiere",&
       !     "Fletcher and Reeves","Hestenes and Stiefel",&
       !     "Fletcher (Conjugate descent)","Liu and Storey",&
       !     "Dai and Yuan","Hager and Zhang"),&
       !     enum_i_vals=(/cg_zero,cg_polak_ribiere,cg_fletcher_reeves,&
       !                   cg_hestenes_stiefel,cg_fletcher,cg_liu_storey,&
       !                   cg_dai_yuan,cg_hager_zhang/),&
       !     error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="GRAD_TARGET",&
       !     description="Target value of the MAX gradient norm",&
       !     usage="GRAD_TARGET 1.E-6", default_r_val=1.0E-4_dp,error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="LIN_SEARCH_GRAD_TARGET",&
       !     description="Target value of the gradient for the linear search",&
       !     usage="LIN_SEARCH_GRAD_TARGET 1.E-6", default_r_val=1.0E-3_dp,error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

       !CALL keyword_create(keyword, name="STEP_SIZE_GUESS",&
       !     description="Target value of the MAX gradient norm",&
       !     usage="STEP_SIZE_GUESS 0.5", default_r_val=1.0_dp,error=error)
       !CALL section_add_keyword(section,keyword,error=error)
       !CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_optimizer_section

END MODULE input_cp2k_almo
