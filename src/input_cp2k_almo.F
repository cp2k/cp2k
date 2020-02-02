!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief input for the ALMO SCF section
!> \author Rustam Khaliullin
! **************************************************************************************************
MODULE input_cp2k_almo
   USE bibliography,                    ONLY: Khaliullin2007,&
                                              Khaliullin2008,&
                                              Khaliullin2013,&
                                              Scheiber2018
   USE cp_output_handling,              ONLY: cp_print_key_section_create,&
                                              low_print_level
   USE input_constants,                 ONLY: &
        almo_deloc_none, almo_deloc_scf, almo_deloc_x, almo_deloc_x_then_scf, &
        almo_deloc_xalmo_1diag, almo_deloc_xalmo_scf, almo_deloc_xalmo_x, almo_frz_crystal, &
        almo_frz_none, almo_scf_diag, almo_scf_pcg, almo_scf_skip, almo_scf_trustr, atomic_guess, &
        cg_dai_yuan, cg_fletcher, cg_fletcher_reeves, cg_hager_zhang, cg_hestenes_stiefel, &
        cg_liu_storey, cg_polak_ribiere, cg_zero, molecular_guess, op_loc_berry, op_loc_pipek, &
        optimizer_diis, optimizer_lin_eq_pcg, optimizer_pcg, optimizer_trustr, &
        spd_inversion_dense_cholesky, spd_inversion_ls_hotelling, spd_inversion_ls_taylor, &
        trustr_cauchy, trustr_dogleg, trustr_steihaug, xalmo_prec_domain, xalmo_prec_full, &
        xalmo_prec_zero, xalmo_trial_r0_out, xalmo_trial_simplex
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
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_almo'

   INTEGER, PARAMETER, PRIVATE :: optimizer_block_diagonal_diis = 1
   INTEGER, PARAMETER, PRIVATE :: optimizer_block_diagonal_pcg = 2
   INTEGER, PARAMETER, PRIVATE :: optimizer_xalmo_pcg = 3
   INTEGER, PARAMETER, PRIVATE :: optimizer_xalmo_trustr = 4
   INTEGER, PARAMETER, PRIVATE :: optimizer_newton_pcg_solver = 5
   INTEGER, PARAMETER, PRIVATE :: optimizer_nlmo_pcg = 6

   PUBLIC :: create_almo_scf_section

CONTAINS

! **************************************************************************************************
!> \brief create the almo scf section
!> \param section ...
!> \par History
!>       2011.05 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE create_almo_scf_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_almo_scf_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="ALMO_SCF", &
                          description="Settings for a class of efficient linear scaling methods based "// &
                          "on absolutely localized orbitals"// &
                          " (ALMOs). ALMO methods are currently restricted to closed-shell molecular systems.", &
                          n_keywords=10, n_subsections=6, repeats=.FALSE., &
                          citations=(/Khaliullin2013, Scheiber2018/))

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_FILTER", &
                          description="Threshold for the matrix sparsity filter", &
                          usage="EPS_FILTER 1.e-6", default_r_val=1.e-7_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ALMO_SCF_GUESS", &
                          description="The method to generate initial ALMOs.", &
                          usage="ALMO_SCF_GUESS MOLECULAR", &
                          default_i_val=atomic_guess, &
                          enum_c_vals=s2a("MOLECULAR", "ATOMIC"), &
                          enum_desc=s2a("SCF calculations on single molecules controlled by the regular SCF "// &
                                        "keywords outside ALMO options. This kind of calculation is expensive "// &
                                        "and only recommended if ALMO SCF does not converge from the ATOMIC guess.", &
                                        "Superposition of atomic densities."), &
                          enum_i_vals=(/molecular_guess, atomic_guess/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MO_OVERLAP_INV_ALG", &
                          description="Algorithm to invert MO overlap matrix.", &
                          usage="MO_OVERLAP_INV_ALG LS_HOTELLING", &
                          default_i_val=spd_inversion_ls_hotelling, &
                          enum_c_vals=s2a("LS_HOTELLING", "LS_TAYLOR", "DENSE_CHOLESKY"), &
                          enum_desc=s2a("Linear scaling iterative Hotelling algorithm. Fast for large sparse matrices.", &
                                        "Linear scaling algorithm based on Taylor expansion of (A+B)^(-1).", &
                                        "Stable but dense Cholesky algorithm. Cubically scaling."), &
                          enum_i_vals=(/spd_inversion_ls_hotelling, spd_inversion_ls_taylor, spd_inversion_dense_cholesky/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      !CALL keyword_create(keyword, __LOCATION__, name="STOP_SCF_EARLY",&
      !     description="Stops SCF using EPS_ERROR_EARLY or MAX_ITER_EARLY", &
      !     usage="STOP_SCF_EARLY .TRUE.", default_l_val=.FALSE.,&
      !     lone_keyword_l_val=.TRUE.)
      !CALL section_add_keyword(section,keyword)
      !CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ALMO_EXTRAPOLATION_ORDER", &
                          description="Number of previous states used for the ASPC extrapolation of the ALMO "// &
                          "initial guess. 0 implies that the guess is given by ALMO_SCF_GUESS at each step.", &
                          usage="ALMO_EXTRAPOLATION_ORDER 3", default_i_val=3)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="XALMO_EXTRAPOLATION_ORDER", &
                          description="Number of previous states used for the ASPC extrapolation of the initial guess "// &
                          "for the delocalization correction.", &
                          usage="XALMO_EXTRAPOLATION_ORDER 1", default_i_val=0)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ALMO_ALGORITHM", &
                          description="Specifies the algorithm to update block-diagonal ALMOs during the SCF procedure.", &
                          usage="ALMO_ALGORITHM DIAG", &
                          default_i_val=almo_scf_diag, &
                          !enum_c_vals=s2a("DIAG", "DM_SIGN","PCG"),&
                          !enum_c_vals=s2a("DIAG", "PCG"), &
                          enum_c_vals=s2a("DIAG", "PCG", "SKIP"), &
                          enum_desc=s2a("DIIS-accelerated diagonalization controlled by ALMO_OPTIMIZER_DIIS. "// &
                                        "Recommended for large systems containing small fragments.", &
                                        !"Update the density matrix using linear scaling routines. "//&
                                        !"Recommended if large fragments are present.",&
                                        "Energy minimization with a PCG algorithm controlled by ALMO_OPTIMIZER_PCG.", &
                                        "Skip optimization of block-diagonal ALMOs."), &
                          !enum_i_vals=(/almo_scf_diag,almo_scf_dm_sign,almo_scf_pcg/),&
                          enum_i_vals=(/almo_scf_diag, almo_scf_pcg, almo_scf_skip/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="XALMO_ALGORITHM", &
                          description="Specifies the algorithm to update ALMOs on eXtended domains (XALMOs).", &
                          usage="XALMO_ALGORITHM TRUST_REGION", &
                          default_i_val=almo_scf_pcg, &
                          enum_c_vals=s2a("DIAG", "PCG", "TRUST_REGION"), &
                          enum_desc=s2a("DIIS-accelerated diagonalization.", &
                                        "PCG algorithm controlled by XALMO_OPTIMIZER_PCG.", &
                                        "Trust-region optimizer controlled by XALMO_OPTIMIZER_TRUSTR"), &
                          enum_i_vals=(/almo_scf_diag, almo_scf_pcg, almo_scf_trustr/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="XALMO_TRIAL_WF", &
                          description="Determines the form of the trial XALMOs.", &
                          usage="XALMO_TRIAL_WF SIMPLE", &
                          default_i_val=xalmo_trial_r0_out, &
                          enum_c_vals=s2a("SIMPLE", "PROJECT_R0_OUT"), &
                          enum_desc=s2a("Straightforward AO-basis expansion.", &
                                        "Block-diagonal ALMOs plus the XALMO term projected onto the unoccupied "// &
                                        "ALMO-subspace."), &
                          enum_i_vals=(/xalmo_trial_simplex, xalmo_trial_r0_out/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="DELOCALIZE_METHOD", &
         description="Methods to reintroduce electron delocalization, which is excluded "// &
         "with the block-diagonal ALMO reference. Electron delocalization can "// &
         "be computed using either fully delocalized MOs or spatially restricted "// &
         "ALMOs (called eXtended ALMOs or XALMOs). All methods below use either a PCG or trust-region"// &
         "optimizers controlled by XALMO_OPTIMIZER_* subsections. The only exception is "// &
         "the non-iterative XALMO_1DIAG.", &
         usage="DELOCALIZE_METHOD XALMO_X", &
         default_i_val=almo_deloc_xalmo_x, &
         enum_c_vals=s2a("NONE", "XALMO_1DIAG", "XALMO_X", "XALMO_SCF", "FULL_X", "FULL_SCF", "FULL_X_THEN_SCF"), &
         enum_desc=s2a("Neglect electron delocalization", &
                       "Correction based on one diagonalization of the spatially projected Hamiltonian (XALMO)", &
                       "Single excitation correction (no Hamiltonian re-build) with spatial restrictions (XALMO)", &
                       "Self-consistent treatment of delocalization with spatial restrictions (XALMO)", &
                       "Single excitation correction (no Hamiltonian re-build) without spatial restrictions", &
                       "Self-consistent treatment of delocalization without spatial restrictions", &
                       "Single excitation correction followed by full SCF procedure, both without spatial restrictions"), &
         enum_i_vals=(/almo_deloc_none, almo_deloc_xalmo_1diag, almo_deloc_xalmo_x, almo_deloc_xalmo_scf, &
                       almo_deloc_x, almo_deloc_scf, almo_deloc_x_then_scf/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="XALMO_R_CUTOFF_FACTOR", &
                          description="Controls the localization radius of XALMOs: "// &
                          !"r0 = r0_factor*(radius(at1)+radius(at2)) + r0_shift",&
                          "R_cutoff = XALMO_R_CUTOFF_FACTOR*(radius(at1)+radius(at2))", &
                          usage="XALMO_R_CUTOFF_FACTOR 1.6", default_r_val=1.50_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RETURN_ORTHOGONALIZED_MOS", &
                          description="Orthogonalize final ALMOs before they are returned"// &
                          "to Quickstep (i.e. for calculation of properties)", &
                          usage="RETURN_ORTHOGONALIZED_MOS .TRUE.", default_l_val=.TRUE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CONSTRUCT_NLMOS", &
                          description="Turns on post-SCF construction of NLMOs", &
                          usage="CONSTRUCT_NLMOS .TRUE.", default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      !CALL keyword_create(keyword, __LOCATION__, name="DOMAIN_LAYOUT_MOS",&
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

      !CALL keyword_create(keyword, __LOCATION__, name="DOMAIN_LAYOUT_AOS",&
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

      !CALL keyword_create(keyword, __LOCATION__, name="MATRIX_CLUSTERING_MOS",&
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

      !CALL keyword_create(keyword, __LOCATION__, name="MATRIX_CLUSTERING_AOS",&
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

      !CALL keyword_create(keyword, __LOCATION__, name="CONSTRAINT_TYPE",&
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

      !CALL keyword_create(keyword, __LOCATION__, name="QUENCHER_RADIUS_TYPE",&
      !     description="Determines the type of atomic radii used for imposing the ALMO constraints",&
      !     usage="QUENCHER_RADIUS_TYPE VDW",&
      !     default_i_val=do_bondparm_vdw,&
      !     enum_c_vals=s2a("COVALENT", "VDW"),&
      !     enum_desc=s2a("Covalent atomic radii",&
      !     "Van der Waals atomic radii"),&
      !     enum_i_vals=(/do_bondparm_covalent,do_bondparm_vdw/))
      !CALL section_add_keyword(section,keyword)
      !CALL keyword_release(keyword)

      !CALL keyword_create(keyword, __LOCATION__, name="QUENCHER_R0_FACTOR",&
      !     description="Parameter to calculate the inner soft cutoff radius: "//&
      !     !"r0 = r0_factor*(radius(at1)+radius(at2)) + r0_shift",&
      !     "r0 = r0_factor*(radius(at1)+radius(at2))",&
      !     usage="QUENCHER_R0_FACTOR 1.05", default_r_val=1.05_dp)
      !CALL section_add_keyword(section,keyword)
      !CALL keyword_release(keyword)

      !!CALL keyword_create(keyword, __LOCATION__, name="QUENCHER_R0_SHIFT",&
      !!     description="Parameter to calculate the inner soft cutoff radius (in Angstrom): "//&
      !!     "r0 = r0_factor*(radius(at1)+radius(at2)) + r0_shift",&
      !!     usage="QUENCHER_R0_SHIFT 0.0", default_r_val=0.0_dp)
      !!
      !!CALL section_add_keyword(section,keyword)
      !!CALL keyword_release(keyword)

      !CALL keyword_create(keyword, __LOCATION__, name="QUENCHER_R1_FACTOR",&
      !     description="Parameter to calculate the outer soft cutoff radius: "//&
      !     !"r1 = r1_factor*(radius(at1)+radius(at2)) + r1_shift",&
      !     "r1 = r1_factor*(radius(at1)+radius(at2))",&
      !     usage="QUENCHER_R1_FACTOR 1.55", default_r_val=1.55_dp)
      !CALL section_add_keyword(section,keyword)
      !CALL keyword_release(keyword)

      !!CALL keyword_create(keyword, __LOCATION__, name="QUENCHER_R1_SHIFT",&
      !!     description="Parameter to calculate the outer soft cutoff radius (in Angstrom): "//&
      !!     "r1 = r1_factor*(radius(at1)+radius(at2)) + r1_shift",&
      !!     usage="QUENCHER_R1_SHIFT 0.0", default_r_val=0.0_dp)
      !!
      !!CALL section_add_keyword(section,keyword)
      !!CALL keyword_release(keyword)

      !CALL keyword_create(keyword, __LOCATION__, name="QUENCHER_AO_OVERLAP_0",&
      !     description="Overlap value of the inner soft cutoff",&
      !     usage="QUENCHER_AO_OVERLAP_0 1.0E-4", default_r_val=1.0E-4_dp)
      !CALL section_add_keyword(section,keyword)
      !CALL keyword_release(keyword)

      !CALL keyword_create(keyword, __LOCATION__, name="QUENCHER_AO_OVERLAP_1",&
      !     description="Overlap value of the outer soft cutoff",&
      !     usage="QUENCHER_AO_OVERLAP_1 1.0E-6", default_r_val=1.0E-6_dp)
      !CALL section_add_keyword(section,keyword)
      !CALL keyword_release(keyword)

      !CALL keyword_create(keyword, __LOCATION__, name="ENVELOPE_AMPLITUDE",&
      !     description="Defines an upper bound on the maximum norm of the MO coefficients",&
      !     usage="ENVELOPE_AMPLITUDE 1.0", default_r_val=1.0_dp)
      !CALL section_add_keyword(section,keyword)
      !CALL keyword_release(keyword)

      !CALL keyword_create(keyword, __LOCATION__, name="DELOC_CAYLEY_TENSOR_TYPE",&
      !     description="Tensor properties of occupied and virtual indices",&
      !     usage="DELOC_CAYLEY_TENSOR_TYPE ORTHOGONAL",&
      !     default_i_val=tensor_orthogonal,&
      !     enum_c_vals=s2a("ORTHOGONAL", "BIORTHOGONAL"),&
      !     enum_desc=s2a("Orthogonalize both occupied and virtual orbitals",&
      !     "Contravariant virtual (MOs or AOs) and covariant occupied orbitals"),&
      !     enum_i_vals=(/tensor_orthogonal,tensor_up_down/))
      !CALL section_add_keyword(section,keyword)
      !CALL keyword_release(keyword)

      !CALL keyword_create(keyword, __LOCATION__, name="DELOC_CAYLEY_CONJUGATOR",&
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

      !CALL keyword_create(keyword, __LOCATION__, name="DELOC_CAYLEY_MAX_ITER",&
      !     description="Maximum number of CG iterations to solve Ricatti equations",&
      !     usage="DELOC_CAYLEY_MAX_ITER 100",default_i_val=50)
      !CALL section_add_keyword(section,keyword)
      !CALL keyword_release(keyword)

      !CALL keyword_create(keyword, __LOCATION__, name="DELOC_CAYLEY_EPS_CONVERGENCE",&
      !     description="Convergence criterion of the CG algorithm",&
      !     usage="DELOC_CAYLEY_EPS_CONVERGENCE 1.e-6", default_r_val=1.e-7_dp)
      !CALL section_add_keyword(section,keyword)
      !CALL keyword_release(keyword)

      !CALL keyword_create(keyword, __LOCATION__, name="DELOC_CAYLEY_VIR_PRECOND",&
      !     description="Use preconditioner for the virtual subspace",&
      !     usage="DELOC_CAYLEY_VIR_PRECOND .TRUE.", default_l_val=.TRUE.,&
      !     lone_keyword_l_val=.TRUE.)
      !CALL section_add_keyword(section,keyword)
      !CALL keyword_release(keyword)

      !CALL keyword_create(keyword, __LOCATION__, name="DELOC_CAYLEY_OCC_PRECOND",&
      !     description="Use preconditioner for the occupied subspace",&
      !     usage="DELOC_CAYLEY_OCC_PRECOND .TRUE.", default_l_val=.TRUE.,&
      !     lone_keyword_l_val=.TRUE.)
      !CALL section_add_keyword(section,keyword)
      !CALL keyword_release(keyword)

      !CALL keyword_create(keyword, __LOCATION__, name="DELOC_TRUNCATE_VIRTUALS",&
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

      !CALL keyword_create(keyword, __LOCATION__, name="DELOC_VIRT_PER_DOMAIN",&
      !     description="Number of virtual orbitals (per domain, atom or molecule) "//&
      !     "retained to obtain the delocalization correction",&
      !     usage="DELOC_VIRT_PER_DOMAIN",default_i_val=-1)
      !CALL section_add_keyword(section,keyword)
      !CALL keyword_release(keyword)

      !CALL keyword_create(keyword, __LOCATION__, name="DELOC_USE_OCC_ORBS",&
      !     description="Use occupied orbitals (as opposed to density matrix) "//&
      !     "to calculate correction for electron delocalization",&
      !     usage="DELOC_USE_OCC_ORBS .TRUE.", default_l_val=.TRUE.,&
      !     lone_keyword_l_val=.TRUE.)
      !CALL section_add_keyword(section,keyword)
      !CALL keyword_release(keyword)

      !CALL keyword_create(keyword, __LOCATION__, name="DELOC_CAYLEY_USE_VIRT_ORBS",&
      !     description="Use virtual orbitals (as opposed to the 1-P projector) "//&
      !     "to calculate correction for electron delocalization. Works only if "//&
      !     "DELOC_USE_OCC_ORBS is set to TRUE",&
      !     usage="DELOC_CAYLEY_USE_VIRT_ORBS .TRUE.", default_l_val=.FALSE.,&
      !     lone_keyword_l_val=.TRUE.)
      !CALL section_add_keyword(section,keyword)
      !CALL keyword_release(keyword)

      !CALL keyword_create(keyword, __LOCATION__, name="DELOC_CAYLEY_LINEAR",&
      !     description="Neglect the quadratic term in the Riccati equations. "//&
      !     "Equivalent to the first order correction to the occupied orbitals "//&
      !     "(second order correction to the energy)",&
      !     usage="DELOC_CAYLEY_LINEAR .FALSE.", default_l_val=.FALSE.,&
      !     lone_keyword_l_val=.TRUE.)
      !CALL section_add_keyword(section,keyword)
      !CALL keyword_release(keyword)

      !CALL keyword_create(keyword, __LOCATION__, name="OPT_K_OUTER_MAX_ITER",&
      !     description="Maximum number of outer loop iterations to optimize retained virtual orbitals",&
      !     usage="OPT_K_OUTER_MAX_ITER 10",default_i_val=1)
      !CALL section_add_keyword(section,keyword)
      !CALL keyword_release(keyword)

      !CALL keyword_create(keyword, __LOCATION__, name="OPT_K_MAX_ITER",&
      !     description="Maximum number of iterations to optimize retained virtual orbitals",&
      !     usage="OPT_K_MAX_ITER 100",default_i_val=100)
      !CALL section_add_keyword(section,keyword)
      !CALL keyword_release(keyword)

      !CALL keyword_create(keyword, __LOCATION__, name="OPT_K_EPS_CONVERGENCE",&
      !     description="Convergence criterion of the optimization algorithm",&
      !     usage="OPT_K_EPS_CONVERGENCE 1.e-5", default_r_val=1.e-5_dp)
      !CALL section_add_keyword(section,keyword)
      !CALL keyword_release(keyword)

      !CALL keyword_create(keyword, __LOCATION__, name="OPT_K_TRIAL_STEP_SIZE",&
      !     description="Size of the trial step along the gradient",&
      !     usage="OPT_K_TRIAL_STEP_SIZE 0.05", default_r_val=0.05_dp)
      !CALL section_add_keyword(section,keyword)
      !CALL keyword_release(keyword)

      !CALL keyword_create(keyword, __LOCATION__, name="OPT_K_TRIAL_STEP_SIZE_MULTIPLIER",&
      !     description="The trial step size is obtained by multiplying the optimal step size "//&
      !     "from the previous iteration",&
      !     usage="OPT_K_TRIAL_STEP_SIZE_multiplier 1.0", default_r_val=1.4_dp)
      !CALL section_add_keyword(section,keyword)
      !CALL keyword_release(keyword)

      !CALL keyword_create(keyword, __LOCATION__, name="OPT_K_CONJ_ITER_START",&
      !     description="Iteration for switching from the steepest descent algorithm "//&
      !     "to conjugate gradient",&
      !     usage="OPT_K_CONJ_ITER_START 5",default_i_val=0)
      !CALL section_add_keyword(section,keyword)
      !CALL keyword_release(keyword)

      !CALL keyword_create(keyword, __LOCATION__, name="OPT_K_CONJ_ITER_FREQ_RESET",&
      !     description="Reset frequency of the conjugate gradient direction",&
      !     usage="OPT_K_CONJ_ITER_FREQ_RESET 20",default_i_val=1000000)
      !CALL section_add_keyword(section,keyword)
      !CALL keyword_release(keyword)

      !CALL keyword_create(keyword, __LOCATION__, name="OPT_K_CONJUGATOR",&
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

      !CALL keyword_create(keyword, __LOCATION__, name="OPT_K_PREC_ITER_START",&
      !     description="Start using the preconditioner (approximate preconditioners "//&
      !     "might not be valid on early iterations)",&
      !     usage="OPT_K_PREC_ITER_START 2",default_i_val=0)
      !CALL section_add_keyword(section,keyword)
      !CALL keyword_release(keyword)

      !CALL keyword_create(keyword, __LOCATION__, name="OPT_K_PREC_ITER_FREQ_UPDATE",&
      !     description="Frequency for updating the preconditioner",&
      !     usage="OPT_K_PREC_ITER_FREQ_UPDATE 10",default_i_val=1)
      !CALL section_add_keyword(section,keyword)
      !CALL keyword_release(keyword)

      NULLIFY (subsection)
      CALL create_optimizer_section(subsection, optimizer_block_diagonal_diis)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (subsection)
      CALL create_optimizer_section(subsection, optimizer_block_diagonal_pcg)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (subsection)
      CALL create_optimizer_section(subsection, optimizer_xalmo_pcg)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (subsection)
      CALL create_optimizer_section(subsection, optimizer_xalmo_trustr)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (subsection)
      CALL create_optimizer_section(subsection, optimizer_nlmo_pcg)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (subsection)
      CALL create_matrix_iterate_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (subsection)
      CALL create_almo_analysis_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_almo_scf_section

! **************************************************************************************************
!> \brief The optimizer section is a collection of keywords that are similar
!>        to all optimization methods (e.g. target error, number of iterations)
!> \param section ...
!> \param optimizer_id   allows to adapt the standard section for specific needs
!> \par History
!>       2012.03 created [Rustam Z Khaliullin]
!>       2014.10 fully integrated [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   RECURSIVE SUBROUTINE create_optimizer_section(section, optimizer_id)

      TYPE(section_type), POINTER                        :: section
      INTEGER, INTENT(IN)                                :: optimizer_id

      CHARACTER(len=*), PARAMETER :: routineN = 'create_optimizer_section', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: optimizer_type
      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      NULLIFY (section)

      ! choose the name of the section
      SELECT CASE (optimizer_id)
      CASE (optimizer_block_diagonal_diis)
         CALL section_create(section, __LOCATION__, name="ALMO_OPTIMIZER_DIIS", &
                             description="Controls the iterative DIIS-accelerated optimization of block-diagonal ALMOs.", &
                             n_keywords=5, n_subsections=0, repeats=.FALSE.)
         optimizer_type = optimizer_diis
      CASE (optimizer_block_diagonal_pcg)
         CALL section_create(section, __LOCATION__, name="ALMO_OPTIMIZER_PCG", &
                             description="Controls the PCG optimization of block-diagonal ALMOs.", &
                             n_keywords=9, n_subsections=1, repeats=.FALSE.)
         optimizer_type = optimizer_pcg
      CASE (optimizer_nlmo_pcg)
         CALL section_create(section, __LOCATION__, name="NLMO_OPTIMIZER_PCG", &
                             description="Controls the PCG optimization of nonorthogonal localized MOs.", &
                             n_keywords=9, n_subsections=1, repeats=.FALSE.)
         optimizer_type = optimizer_pcg
         NULLIFY (subsection)
         CALL create_penalty_section(subsection)
         CALL section_add_subsection(section, subsection)
         CALL section_release(subsection)
      CASE (optimizer_xalmo_pcg)
         CALL section_create(section, __LOCATION__, name="XALMO_OPTIMIZER_PCG", &
                             description="Controls the PCG optimization of extended ALMOs.", &
                             n_keywords=10, n_subsections=2, repeats=.FALSE.)
         NULLIFY (subsection)
         CALL create_optimizer_section(subsection, optimizer_newton_pcg_solver)
         CALL section_add_subsection(section, subsection)
         CALL section_release(subsection)
         optimizer_type = optimizer_pcg
      CASE (optimizer_xalmo_trustr)
         CALL section_create(section, __LOCATION__, name="XALMO_OPTIMIZER_TRUSTR", &
                             description="Controls the trust-region optimization of extended ALMOs. "// &
                             "Trust radius is varied in the outer loop. Once the trust radius is "// &
                             "chosen (and fixed) the model function can be minized using varius "// &
                             "approaches. Currently, an iterative conjugate-gradient approach is "// &
                             "used and controled by the inner loop", &
                             n_keywords=10, n_subsections=0, repeats=.FALSE.)
         optimizer_type = optimizer_trustr
      CASE (optimizer_newton_pcg_solver)
         CALL section_create(section, __LOCATION__, name="XALMO_NEWTON_PCG_SOLVER", &
                             description="Controls an iterative solver of the Newton-Raphson linear equation.", &
                             n_keywords=4, n_subsections=0, repeats=.FALSE.)
         optimizer_type = optimizer_lin_eq_pcg
      CASE DEFAULT
         CPABORT("No default values allowed")
      END SELECT

      NULLIFY (keyword)

      ! add common keywords
      CALL keyword_create(keyword, __LOCATION__, name="MAX_ITER", &
                          description="Maximum number of iterations", &
                          usage="MAX_ITER 100", default_i_val=20)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_ERROR", &
                          description="Target value of the MAX norm of the error", &
                          usage="EPS_ERROR 1.E-6", default_r_val=1.0E-5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! add common keywords
      CALL keyword_create(keyword, __LOCATION__, name="MAX_ITER_EARLY", &
                          description="Maximum number of iterations for truncated SCF "// &
                          "(e.g. Langevin-corrected MD). Negative values mean that this keyword is not used.", &
                          usage="MAX_ITER_EARLY 5", default_i_val=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_ERROR_EARLY", &
                          description="Target value of the MAX norm of the error for truncated SCF "// &
                          "(e.g. Langevin-corrected MD). Negative values mean that this keyword is not used.", &
                          usage="EPS_ERROR_EARLY 1.E-2", default_r_val=-1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! add keywords specific to each type
      IF (optimizer_type .EQ. optimizer_diis) THEN

         CALL keyword_create(keyword, __LOCATION__, name="N_DIIS", &
                             description="Number of error vectors to be used in the DIIS "// &
                             "optimization procedure", &
                             usage="N_DIIS 5", default_i_val=6)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

      ENDIF

      IF (optimizer_type .EQ. optimizer_pcg .OR. &
          optimizer_type .EQ. optimizer_lin_eq_pcg .OR. &
          optimizer_type .EQ. optimizer_trustr) THEN

         CALL keyword_create(keyword, __LOCATION__, name="MAX_ITER_OUTER_LOOP", &
                             description="Maximum number of iterations in the outer loop. "// &
                             "Use the outer loop to update the preconditioner and reset the conjugator. "// &
                             "This can speed up convergence significantly.", &
                             usage="MAX_ITER 10", default_i_val=0)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

         CALL keyword_create(keyword, __LOCATION__, name="PRECONDITIONER", &
                             description="Select a preconditioner for the conjugate gradient optimization", &
                             usage="PRECONDITIONER DOMAIN", &
                             default_i_val=xalmo_prec_domain, &
                             enum_c_vals=s2a("NONE", "DEFAULT", "DOMAIN", "FULL"), &
                             enum_desc=s2a("Do not use preconditioner", &
                                           "Same as DOMAIN preconditioner", &
                                           "Invert preconditioner domain-by-domain."// &
                                           " The main component of the linear scaling algorithm", &
                                           "Solve linear equations step=-H.grad on the entire space"), &
                             enum_i_vals=(/xalmo_prec_zero, xalmo_prec_domain, &
                                           xalmo_prec_domain, xalmo_prec_full/))
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

      ENDIF

      IF (optimizer_type .EQ. optimizer_pcg) THEN

         CALL keyword_create(keyword, __LOCATION__, name="LIN_SEARCH_EPS_ERROR", &
                             description="Target value of the gradient norm during the linear search", &
                             usage="LIN_SEARCH_EPS_ERROR 1.E-2", default_r_val=1.0E-3_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

         CALL keyword_create(keyword, __LOCATION__, name="LIN_SEARCH_STEP_SIZE_GUESS", &
                             description="The size of the first step in the linear search", &
                             usage="LIN_SEARCH_STEP_SIZE_GUESS 0.1", default_r_val=1.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

         CALL keyword_create(keyword, __LOCATION__, name="PRECOND_FILTER_THRESHOLD", &
                             description="Select eigenvalues of the preconditioner "// &
                             "that are smaller than the threshold and project out the "// &
                             "corresponding eigenvectors from the gradient. No matter "// &
                             "how large the threshold is the maximum number of projected "// &
                             "eienvectors for a fragment equals to the number of occupied "// &
                             "orbitals of fragment's neighbors.", &
                             usage="PRECOND_FILTER_THRESHOLD 0.1", default_r_val=0.5_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

      ENDIF

      IF (optimizer_type .EQ. optimizer_pcg .OR. &
          optimizer_type .EQ. optimizer_trustr) THEN

         CALL keyword_create(keyword, __LOCATION__, name="CONJUGATOR", &
                             description="Various methods to compute step directions in the PCG optimization", &
                             usage="CONJUGATOR POLAK_RIBIERE", &
                             default_i_val=cg_hager_zhang, &
                             enum_c_vals=s2a("ZERO", "POLAK_RIBIERE", "FLETCHER_REEVES", &
                                             "HESTENES_STIEFEL", "FLETCHER", "LIU_STOREY", "DAI_YUAN", "HAGER_ZHANG"), &
                             enum_desc=s2a("Steepest descent", "Polak and Ribiere", &
                                           "Fletcher and Reeves", "Hestenes and Stiefel", &
                                           "Fletcher (Conjugate descent)", "Liu and Storey", &
                                           "Dai and Yuan", "Hager and Zhang"), &
                             enum_i_vals=(/cg_zero, cg_polak_ribiere, cg_fletcher_reeves, &
                                           cg_hestenes_stiefel, cg_fletcher, cg_liu_storey, &
                                           cg_dai_yuan, cg_hager_zhang/))
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

      ENDIF

      IF (optimizer_type .EQ. optimizer_trustr) THEN

         CALL keyword_create(keyword, __LOCATION__, name="ALGORITHM", &
                             description="Selects an algorithm to solve the fixed-radius subproblem", &
                             usage="ALGORITHM CG", &
                             default_i_val=trustr_cauchy, &
                             enum_c_vals=s2a("CG", "CAUCHY", "DOGLEG"), &
                             enum_desc=s2a("Steihaug's iterative CG algorithm that does not invert model Hessian", &
                                           "Compute simple Cauchy point", &
                                           "Dogleg optimizer"), &
                             enum_i_vals=(/trustr_steihaug, trustr_cauchy, trustr_dogleg/))
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

         CALL keyword_create(keyword, __LOCATION__, name="ETA", &
                             description="Must be between 0.0 and 0.25. Rho value below which the "// &
                             "optimization of the model function is not accepted and the "// &
                             "optimization is restarted from the same point but decreased "// &
                             "trust radius. Rho is the ratio of the actual over predicted "// &
                             "change in the objective function", &
                             usage="ETA 0.1", default_r_val=0.25_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

         CALL keyword_create(keyword, __LOCATION__, name="MODEL_GRAD_NORM_RATIO", &
                             description="Stop the fixed-trust-radius (inner) loop optimization "// &
                             "once the ratio of the current norm of the model gradient over the "// &
                             "initial norm drops below this threshold", &
                             usage="MODEL_GRAD_NORM_RATIO 1.E-2", default_r_val=0.01_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

         CALL keyword_create(keyword, __LOCATION__, name="INITIAL_TRUST_RADIUS", &
                             description="Initial trust radius", &
                             usage="INITIAL_TRUST_RADIUS 0.1", default_r_val=0.1_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

         CALL keyword_create(keyword, __LOCATION__, name="MAX_TRUST_RADIUS", &
                             description="Maximum allowed trust radius", &
                             usage="MAX_TRUST_RADIUS 1.0", default_r_val=2.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

      ENDIF

   END SUBROUTINE create_optimizer_section

! **************************************************************************************************
!> \brief The section controls iterative matrix operations like SQRT or inverse
!> \param section ...
!> \par History
!>       2017.05 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE create_matrix_iterate_section(section)

      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      NULLIFY (section)

      CALL section_create(section, __LOCATION__, name="MATRIX_ITERATE", &
                          description="Controls linear scaling iterative procedure on matrices: inversion, sqrti, etc. "// &
                          "High-order Lanczos accelerates convergence provided it can estimate the eigenspectrum correctly.", &
                          n_keywords=4, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_TARGET_FACTOR", &
                          description="Multiplication factor that determines acceptable error in the iterative procedure. "// &
                          "Acceptable error = EPS_TARGET_FACTOR * EPS_FILTER", &
                          usage="EPS_TARGET_FACTOR 100.0", default_r_val=10.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_LANCZOS", &
                          description="Threshold for Lanczos eigenvalue estimation.", &
                          usage="EPS_LANCZOS 1.0E-4", default_r_val=1.0E-3_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ORDER_LANCZOS", &
                          description="Order of the Lanczos estimator. Use 0 to turn off. Do not use 1.", &
                          usage="ORDER_LANCZOS 5", default_i_val=3)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_ITER_LANCZOS", &
                          description="Maximum number of Lanczos iterations.", &
                          usage="MAX_ITER_LANCZOS 64", default_i_val=128)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_matrix_iterate_section

! **************************************************************************************************
!> \brief The section controls penalty methods
!> \param section ...
!> \par History
!>       2018.01 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE create_penalty_section(section)

      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      NULLIFY (section)

      CALL section_create(section, __LOCATION__, name="PENALTY", &
                          description="Add penalty terms to the energy functional.", &
                          n_keywords=3, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="OPERATOR", &
         description="Type of opertator which defines the spread functional", &
         usage="OPERATOR PIPEK", &
         enum_c_vals=s2a("BERRY", "PIPEK"), &
         enum_i_vals=(/op_loc_berry, op_loc_pipek/), &
         default_i_val=op_loc_berry)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PENALTY_STRENGTH", &
                          description="Strength of the orthogonalization penalty", &
                          usage="PENALTY_STRENGTH 1.1", default_r_val=1.1_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PENALTY_STRENGTH_DECREASE_FACTOR", &
                          description="Factor that decreases the strength of the orthogonalization penalty.", &
                          usage="PENALTY_STRENGTH_DECREASE_FACTOR 1.1", default_r_val=1.1_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DETERMINANT_TOLERANCE", &
         description="Stop the optimization of the penalty strength if the determinant of the overlap changes less than this tolerance threshold.", &
                          usage="DETERMINANT_TOLERANCE 1.0E-4", default_r_val=1.0E-3_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FINAL_DETERMINANT", &
                          description="The final determinant that obtained after optimization.", &
                          usage="FINAL_DETERMINANT 0.1", default_r_val=0.1_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="COMPACTIFICATION_FILTER_START", &
                          description="Set orbital coefficients with absolute value smaller than this value to zero.", &
                          usage="COMPACTIFICATION_FILTER_START 1.e-6", default_r_val=-1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="VIRTUAL_NLMOS", &
                          description="Localize virtual oribtals", &
                          usage="VIRTUAL_NLMOS .TRUE.", default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_penalty_section

! **************************************************************************************************
!> \brief The section controls electronic structure analysis based on ALMOs
!> \param section ...
!> \par History
!>       2014.10 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE create_almo_analysis_section(section)

      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_almo_analysis_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      NULLIFY (section)

      CALL section_create(section, __LOCATION__, name="ANALYSIS", &
                          description="Controls electronic structure analysis based on ALMOs and XALMOs.", &
                          n_keywords=2, n_subsections=1, repeats=.FALSE., &
                          citations=(/Khaliullin2007, Khaliullin2008/))

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="Activation of ALMO-based electronic structure analysis.", &
                          usage="&ANALYSIS T", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FROZEN_MO_ENERGY_TERM", &
                          description="Perform calculations on single molecules to compute the frozen density term", &
                          usage="FROZEN_MO_ENERGY_TERM SUBLATTICE", default_i_val=almo_frz_none, &
                          !enum_c_vals=s2a("SKIP", "ISOLATED", "SUBLATTICE"),&
                          enum_c_vals=s2a("SKIP", "SUBLATTICE"), &
                          enum_desc=s2a("Do not compute the frozen MO energy term.", &
                                        !"Use isolated gas-phase molecules as the reference.",&
                                        "Use energies of single molecules in their positions in the crystal "// &
                                        "cell as the reference. "// &
                                        "This term does not have an interpretation if fragmetns are charged."), &
                          enum_i_vals=(/almo_frz_none, &
                                        !almo_frz_isolated,&
                                        almo_frz_crystal/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      NULLIFY (subsection)
      CALL create_almo_print_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_almo_analysis_section

! *****************************************************************************
!> \brief Create a section with all print keys for ALMO methods
!> \param section ...
!> \par History
!>       2016.05 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE create_almo_print_section(section)

      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_almo_print_section', &
         routineP = moduleN//':'//routineN

      TYPE(section_type), POINTER                        :: print_key

      CPASSERT(.NOT. ASSOCIATED(section))
      NULLIFY (section, print_key)

      CALL section_create(section, __LOCATION__, name="PRINT", &
                          description="Controls printing of detailed ALMO decomosition analysis results", &
                          n_keywords=0, n_subsections=2, repeats=.TRUE.)

      ! create print keys within this subsection
      CALL cp_print_key_section_create(print_key, __LOCATION__, "ALMO_EDA_CT", &
                                       description="Controls printing of the electron transfer terms in "// &
                                       "ALMO energy decomposition analysis. "// &
                                       "File format: first column is the index of the electron acceptor "// &
                                       "fragment, second -- is the index of the electron donor fragment, "// &
                                       "third - the energy change (a.u.) associated with the electron transfer.", &
                                       print_level=low_print_level, common_iter_levels=1, &
                                       filename="EDA")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "ALMO_CTA", &
                                       description="Controls printing of the electron transfer terms in "// &
                                       "ALMO charge transfer analysis. "// &
                                       "File format: first column is the index of the electron acceptor "// &
                                       "fragment, second -- is the index of the electron donor fragment, "// &
                                       "third - the change (a.u.) transferred between the two fragments.", &
                                       print_level=low_print_level, common_iter_levels=1, &
                                       filename="CTA")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

   END SUBROUTINE create_almo_print_section

END MODULE input_cp2k_almo
