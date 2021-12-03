!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief function that build the dft section of the input
!> \par History
!>      10.2005 moved out of input_cp2k [fawzi]
!> \author fawzi
! **************************************************************************************************
MODULE input_cp2k_ec
   USE bibliography,                    ONLY: Niklasson2003,&
                                              VandeVondele2012,&
                                              brehm2018,&
                                              brehm2020,&
                                              brehm2021,&
                                              rycroft2009,&
                                              thomas2015
   USE cp_output_handling,              ONLY: cp_print_key_section_create,&
                                              debug_print_level,&
                                              high_print_level
   USE input_constants,                 ONLY: &
        bqb_opt_exhaustive, bqb_opt_normal, bqb_opt_off, bqb_opt_patient, bqb_opt_quick, &
        ec_diagonalization, ec_functional_harris, ec_ls_solver, ec_matrix_sign, ec_matrix_tc2, &
        ec_matrix_trs4, ec_mo_solver, ec_ot_atomic, ec_ot_diag, ec_ot_gs, kg_cholesky, &
        ls_cluster_atomic, ls_cluster_molecular, ls_s_inversion_hotelling, &
        ls_s_inversion_sign_sqrt, ls_s_preconditioner_atomic, ls_s_preconditioner_molecular, &
        ls_s_preconditioner_none, ls_s_sqrt_ns, ls_s_sqrt_proot, ls_scf_sign_ns, &
        ls_scf_sign_proot, ot_precond_full_all, ot_precond_full_kinetic, ot_precond_full_single, &
        ot_precond_full_single_inverse, ot_precond_none, ot_precond_s_inverse, precond_mlp, &
        voro_radii_cov, voro_radii_unity, voro_radii_user, voro_radii_vdw
   USE input_cp2k_mm,                   ONLY: create_dipoles_section
   USE input_cp2k_xc,                   ONLY: create_xc_section
   USE input_keyword_types,             ONLY: keyword_create,&
                                              keyword_release,&
                                              keyword_type
   USE input_section_types,             ONLY: section_add_keyword,&
                                              section_add_subsection,&
                                              section_create,&
                                              section_release,&
                                              section_type
   USE input_val_types,                 ONLY: char_t,&
                                              integer_t,&
                                              real_t
   USE kinds,                           ONLY: dp
   USE string_utilities,                ONLY: s2a
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_ec'

   PUBLIC :: create_ec_section

CONTAINS

! **************************************************************************************************
!> \brief creates the ENERGY CORRECTION section
!> \param section ...
!> \author JGH
! **************************************************************************************************
   SUBROUTINE create_ec_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))

      NULLIFY (keyword)
      CALL section_create(section, __LOCATION__, name="ENERGY_CORRECTION", &
                          description="Sets the various options for the Energy Correction", &
                          n_keywords=0, n_subsections=2, repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="Controls the activation of the energy_correction", &
                          usage="&ENERGY_CORRECTION T", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! add a special XC section
      NULLIFY (subsection)
      CALL create_xc_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      ! add a section for solver keywords
      NULLIFY (subsection)
      CALL create_ec_solver_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      ! add a print section for properties
      NULLIFY (subsection)
      CALL create_ec_print_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL keyword_create(keyword, __LOCATION__, name="ENERGY_FUNCTIONAL", &
                          description="Functional used in energy correction", &
                          usage="ENERGY_FUNCTIONAL HARRIS", &
                          default_i_val=ec_functional_harris, &
                          enum_c_vals=s2a("HARRIS"), &
                          enum_desc=s2a("Harris functional"), &
                          enum_i_vals=(/ec_functional_harris/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="HARRIS_BASIS", &
                          description="Specifies the type of basis to be used for the KG energy correction."// &
                          "Options are: (1) the default orbital basis (ORBITAL);"// &
                          "(2) the primitive functions of the default orbital basis (PRIMITIVE);"// &
                          "(3) the basis set labeled in Kind section (HARRIS)", &
                          usage="HARRIS_BASIS ORBITAL", &
                          type_of_var=char_t, default_c_val="ORBITAL", n_var=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SKIP_EC", &
                          description="Skip EC calculation if ground-state calculation has not converged.", &
                          usage="SKIP_EC T", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAO", &
                          description="Use modified atomic orbitals (MAO) to solve Harris equation", &
                          usage="MAO T", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAO_MAX_ITER", &
                          description="Maximum iterations in MAO optimization. ", &
                          usage="MAO_MAX_ITER 100 ", default_i_val=0)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAO_EPS_GRAD", &
                          description="Threshold used for MAO iterations. ", &
                          usage="MAO_EPS_GRAD 1.0E-4 ", default_r_val=1.0E-5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAO_EPS1", &
                          description="Occupation threshold used to determine number of MAOs."// &
                          " KIND section MAO keyword sets the minimum.", &
                          usage="MAO_EPS1 0.001 ", default_r_val=1000.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAO_IOLEVEL", &
                          description="Verbosity of MAO output: (0) no output ... (3) fully verbose", &
                          usage="MAO_IOLEVEL 0 ", default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ALGORITHM", &
                          description="Algorithm used to solve KS equation", &
                          usage="ALGORITHM DIAGONALIZATION", &
                          default_i_val=ec_diagonalization, &
                          enum_c_vals=s2a("DIAGONALIZATION", "MATRIX_SIGN", &
                                          "TRS4", "TC2", "OT_DIAG"), &
                          enum_desc=s2a("Diagonalization of KS matrix.", &
                                        "Matrix Sign algorithm", &
                                        "Trace resetting trs4 algorithm", &
                                        "Trace resetting tc2 algorithm", &
                                        "OT diagonalization"), &
                          enum_i_vals=(/ec_diagonalization, ec_matrix_sign, &
                                        ec_matrix_trs4, ec_matrix_tc2, ec_ot_diag/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FACTORIZATION", &
                          description="Algorithm used to calculate factorization of overlap matrix", &
                          usage="FACTORIZATION CHOLESKY", &
                          default_i_val=kg_cholesky, &
                          enum_c_vals=s2a("CHOLESKY"), &
                          enum_desc=s2a("Cholesky factorization of overlap matrix"), &
                          enum_i_vals=(/kg_cholesky/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_DEFAULT", &
                          description="Threshold used for accuracy estimates within energy correction. ", &
                          usage="EPS_DEFAULT 1.0E-7 ", default_r_val=1.0E-7_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! Keywords for LS solver of Harris functional
      CALL keyword_create(keyword, __LOCATION__, name="EPS_FILTER", &
                          description="Threshold used for filtering matrix operations.", &
                          usage="EPS_FILTER 1.0E-12", default_r_val=1.0E-12_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_LANCZOS", &
                          description="Threshold used for lanczos estimates.", &
                          usage="EPS_LANCZOS 1.0E-4", default_r_val=1.0E-3_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_ITER_LANCZOS", &
                          description="Maximum number of lanczos iterations.", &
                          usage="MAX_ITER_LANCZOS ", default_i_val=128)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MU", &
                          description="Value (or initial guess) for the chemical potential,"// &
                          " i.e. some suitable energy between HOMO and LUMO energy.", &
                          usage="MU 0.0", default_r_val=-0.1_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FIXED_MU", &
                          description="Should the calculation be performed at fixed chemical potential,"// &
                          " or should it be found fixing the number of electrons", &
                          usage="FIXED_MU .TRUE.", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="S_PRECONDITIONER", &
                          description="Preconditions S with some appropriate form.", &
                          usage="S_PRECONDITIONER MOLECULAR", &
                          default_i_val=ls_s_preconditioner_atomic, &
                          enum_c_vals=s2a("NONE", "ATOMIC", "MOLECULAR"), &
                          enum_desc=s2a("No preconditioner", &
                                        "Using atomic blocks", &
                                        "Using molecular sub-blocks. Recommended if molecules are defined and not too large."), &
                          enum_i_vals=(/ls_s_preconditioner_none, ls_s_preconditioner_atomic, ls_s_preconditioner_molecular/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="S_SQRT_METHOD", &
                          description="Method for the caclulation of the sqrt of S.", &
                          usage="S_SQRT_METHOD NEWTONSCHULZ", &
                          default_i_val=ls_s_sqrt_ns, &
                          enum_c_vals=s2a("NEWTONSCHULZ", "PROOT"), &
                          enum_desc=s2a("Using a Newton-Schulz-like iteration", &
                                        "Using the p-th root method."), &
                          enum_i_vals=(/ls_s_sqrt_ns, ls_s_sqrt_proot/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="S_SQRT_ORDER", &
                          variants=s2a("SIGN_SQRT_ORDER"), &
                          description="Order of the iteration method for the calculation of the sqrt of S.", &
                          usage="S_SQRT_ORDER 3", default_i_val=3)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SIGN_METHOD", &
                          description="Method used for the computation of the sign matrix.", &
                          usage="SIGN_METHOD NEWTONSCHULZ", &
                          default_i_val=ls_scf_sign_ns, &
                          citations=(/VandeVondele2012, Niklasson2003/), &
                          enum_c_vals=s2a("NEWTONSCHULZ", "PROOT"), &
                          enum_desc=s2a("Newton-Schulz iteration.", &
                                        "p-th order root iteration"), &
                          enum_i_vals=(/ls_scf_sign_ns, ls_scf_sign_proot/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SIGN_ORDER", &
                          description="Order of the method used for the computation of the sign matrix.", &
                          usage="SIGN_ORDER 2", &
                          default_i_val=2)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DYNAMIC_THRESHOLD", &
                          description="Should the threshold for the purification be chosen dynamically", &
                          usage="DYNAMIC_THRESHOLD .TRUE.", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NON_MONOTONIC", &
                          description="Should the purification be performed non-monotonically. Relevant for TC2 only.", &
                          usage="NON_MONOTONIC .TRUE.", default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="MATRIX_CLUSTER_TYPE", &
         description="Specify how atomic blocks should be clustered in the used matrices, in order to improve flop rate, "// &
         "and possibly speedup the matrix multiply. Note that the atomic s_preconditioner can not be used."// &
         "Furthermore, since screening is on matrix blocks, "// &
         "slightly more accurate results can be expected with molecular.", &
         usage="MATRIX_CLUSTER_TYPE MOLECULAR", &
         default_i_val=ls_cluster_atomic, &
         enum_c_vals=s2a("ATOMIC", "MOLECULAR"), &
         enum_desc=s2a("Using atomic blocks", &
                       "Using molecular blocks."), &
         enum_i_vals=(/ls_cluster_atomic, ls_cluster_molecular/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SINGLE_PRECISION_MATRICES", &
                          description="Matrices used within the LS code can be either double or single precision.", &
                          usage="SINGLE_PRECISION_MATRICES", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="S_INVERSION", &
                          description="Method used to compute the inverse of S.", &
                          usage="S_PRECONDITIONER MOLECULAR", &
                          default_i_val=ls_s_inversion_sign_sqrt, &
                          enum_c_vals=s2a("SIGN_SQRT", "HOTELLING"), &
                          enum_desc=s2a("Using the inverse sqrt as obtained from sign function iterations.", &
                                        "Using the Hotellign iteration."), &
                          enum_i_vals=(/ls_s_inversion_sign_sqrt, ls_s_inversion_hotelling/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="REPORT_ALL_SPARSITIES", &
                          description="Run the sparsity report at the end of the SCF", &
                          usage="REPORT_ALL_SPARSITIES", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CHECK_S_INV", &
                          description="Perform an accuracy check on the inverse/sqrt of the s matrix.", &
                          usage="CHECK_S_INV", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OT_INITIAL_GUESS", &
                          description="Initial guess of density matrix used for OT Diagonalization", &
                          usage="OT_INITIAL_GUESS ATOMIC", &
                          default_i_val=ec_ot_atomic, &
                          enum_c_vals=s2a("ATOMIC", "GROUND_STATE"), &
                          enum_desc=s2a("Generate an atomic density using the atomic code", &
                                        "Using the ground-state density."), &
                          enum_i_vals=(/ec_ot_atomic, ec_ot_gs/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_ec_section

! **************************************************************************************************
!> \brief creates the linear scaling solver section
!> \param section ...
!> \author Joost VandeVondele [2010-10], JGH [2019-12]
! **************************************************************************************************
   SUBROUTINE create_ec_solver_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="RESPONSE_SOLVER", &
                          description="Specifies the parameters of the linear scaling solver routines", &
                          n_keywords=24, n_subsections=3, repeats=.FALSE., &
                          citations=(/VandeVondele2012/))

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS", &
                          description="Target accuracy for the convergence of the conjugate gradient.", &
                          usage="EPS 1.e-10", default_r_val=1.e-12_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_FILTER", &
                          description="Threshold used for filtering matrix operations.", &
                          usage="EPS_FILTER 1.0E-10", default_r_val=1.0E-10_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_LANCZOS", &
                          description="Threshold used for lanczos estimates.", &
                          usage="EPS_LANCZOS 1.0E-4", default_r_val=1.0E-3_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_ITER", &
                          description="Maximum number of conjugate gradient iteration"// &
                          "to be performed for one optimization.", &
                          usage="MAX_ITER 200", default_i_val=50)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_ITER_LANCZOS", &
                          description="Maximum number of lanczos iterations.", &
                          usage="MAX_ITER_LANCZOS 128", default_i_val=128)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="METHOD", &
                          description="Algorithm used to solve response equation."// &
                          "Both solver are conjugate gradient based, but use either a vector (MO-coefficient)"// &
                          "or density matrix formalism in the orthonormal AO-basis to obtain response density", &
                          usage="METHOD SOLVER", &
                          default_i_val=ec_ls_solver, &
                          enum_c_vals=s2a("MO_SOLVER", "AO_ORTHO"), &
                          enum_desc=s2a("Solver based on MO (vector) formalism", &
                                        "Solver based on density matrix formalism"), &
                          enum_i_vals=(/ec_mo_solver, ec_ls_solver/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="PRECONDITIONER", &
         description="Type of preconditioner to be used with MO conjugate gradient solver. "// &
         "They differ in effectiveness, cost of construction, cost of application. "// &
         "Properly preconditioned minimization can be orders of magnitude faster than doing nothing."// &
         "Only multi-level conjugate gradient preconditioner (MULTI_LEVEL) available for AO response solver (AO_ORTHO).", &
         usage="PRECONDITIONER FULL_ALL", &
         default_i_val=precond_mlp, &
         enum_c_vals=s2a("FULL_ALL", "FULL_SINGLE_INVERSE", "FULL_SINGLE", "FULL_KINETIC", "FULL_S_INVERSE", &
                         "MULTI_LEVEL", "NONE"), &
         enum_desc=s2a("Most effective state selective preconditioner based on diagonalization, "// &
                       "requires the ENERGY_GAP parameter to be an underestimate of the HOMO-LUMO gap. "// &
                       "This preconditioner is recommended for almost all systems, except very large systems where "// &
                       "make_preconditioner would dominate the total computational cost.", &
                       "Based on H-eS cholesky inversion, similar to FULL_SINGLE in preconditioning efficiency "// &
                       "but cheaper to construct, "// &
                       "might be somewhat less robust. Recommended for large systems.", &
                       "Based on H-eS diagonalisation, not as good as FULL_ALL, but somewhat cheaper to apply. ", &
                       "Cholesky inversion of S and T, fast construction, robust, and relatively good, "// &
                       "use for very large systems.", &
                       "Cholesky inversion of S, not as good as FULL_KINETIC, yet equally expensive.", &
                       "Based on same CG as AO-solver itself, but uses cheaper linear transformation", &
                       "skip preconditioning"), &
         enum_i_vals=(/ot_precond_full_all, ot_precond_full_single_inverse, ot_precond_full_single, &
                       ot_precond_full_kinetic, ot_precond_s_inverse, precond_mlp, ot_precond_none/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="S_PRECONDITIONER", &
                          description="Preconditions S with some appropriate form.", &
                          usage="S_PRECONDITIONER MOLECULAR", &
                          default_i_val=ls_s_preconditioner_atomic, &
                          enum_c_vals=s2a("NONE", "ATOMIC", "MOLECULAR"), &
                          enum_desc=s2a("No preconditioner", &
                                        "Using atomic blocks", &
                                        "Using molecular sub-blocks. Recommended if molecules are defined and not too large."), &
                          enum_i_vals=(/ls_s_preconditioner_none, ls_s_preconditioner_atomic, ls_s_preconditioner_molecular/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="S_SQRT_METHOD", &
                          description="Method for the caclulation of the sqrt of S.", &
                          usage="S_SQRT_METHOD NEWTONSCHULZ", &
                          default_i_val=ls_s_sqrt_ns, &
                          enum_c_vals=s2a("NEWTONSCHULZ", "PROOT"), &
                          enum_desc=s2a("Using a Newton-Schulz-like iteration", &
                                        "Using the p-th root method."), &
                          enum_i_vals=(/ls_s_sqrt_ns, ls_s_sqrt_proot/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="S_SQRT_ORDER", &
                          variants=s2a("SIGN_SQRT_ORDER"), &
                          description="Order of the iteration method for the calculation of the sqrt of S.", &
                          usage="S_SQRT_ORDER 3", default_i_val=3)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="MATRIX_CLUSTER_TYPE", &
         description="Specify how atomic blocks should be clustered in the used matrices, in order to improve flop rate, "// &
         "and possibly speedup the matrix multiply. Note that the atomic s_preconditioner can not be used."// &
         "Furthermore, since screening is on matrix blocks, "// &
         "slightly more accurate results can be expected with molecular.", &
         usage="MATRIX_CLUSTER_TYPE MOLECULAR", &
         default_i_val=ls_cluster_atomic, &
         enum_c_vals=s2a("ATOMIC", "MOLECULAR"), &
         enum_desc=s2a("Using atomic blocks", &
                       "Using molecular blocks."), &
         enum_i_vals=(/ls_cluster_atomic, ls_cluster_molecular/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SINGLE_PRECISION_MATRICES", &
                          description="Matrices used within the LS code can be either double or single precision.", &
                          usage="SINGLE_PRECISION_MATRICES", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="S_INVERSION", &
                          description="Method used to compute the inverse of S.", &
                          usage="S_PRECONDITIONER MOLECULAR", &
                          default_i_val=ls_s_inversion_sign_sqrt, &
                          enum_c_vals=s2a("SIGN_SQRT", "HOTELLING"), &
                          enum_desc=s2a("Using the inverse sqrt as obtained from sign function iterations.", &
                                        "Using the Hotellign iteration."), &
                          enum_i_vals=(/ls_s_inversion_sign_sqrt, ls_s_inversion_hotelling/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RESTART", &
                          description="Restart the response calculation if the restart file exists", &
                          usage="RESTART", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RESTART_EVERY", &
                          description="Restart the conjugate gradient after the specified number of iterations.", &
                          usage="RESTART_EVERY 50", default_i_val=50)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_ec_solver_section

! **************************************************************************************************
!> \brief Create the print dft section
!> \param section the section to create
!> \author fbelle - from create_print_dft_section
! **************************************************************************************************
   SUBROUTINE create_ec_print_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="PRINT", &
                          description="Section of possible print options in EC code.", &
                          n_keywords=0, n_subsections=1, repeats=.FALSE.)

      NULLIFY (print_key, keyword)

      ! Output of BQB volumetric files
      CALL cp_print_key_section_create(print_key, __LOCATION__, name="E_DENSITY_BQB", &
                                       description="Controls the output of the electron density to the losslessly"// &
                                       " compressed BQB file format, see [Brehm2018] (via LibBQB"// &
                                       " see <a href=""https://brehm-research.de/bqb"" target=""_blank"">"// &
                                       " https://brehm-research.de/bqb</a> ). Currently does not work"// &
                                       "  with changing cell vector (NpT ensemble).", &
                                       print_level=debug_print_level + 1, filename="", &
                                       citations=(/Brehm2018/))

      CALL keyword_create(keyword, __LOCATION__, name="SKIP_FIRST", &
                          description="Skips the first step of a MD run (avoids duplicate step if restarted).", &
                          usage="SKIP_FIRST T", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="STORE_STEP_NUMBER", &
                          description="Stores the step number and simulation time in the comment line of each BQB"// &
                          " frame. Switch it off for binary compatibility with original CP2k CUBE files.", &
                          usage="STORE_STEP_NUMBER F", default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CHECK", &
                          description="Performs an on-the-fly decompression of each compressed BQB frame to check"// &
                          " whether the volumetric data exactly matches, and aborts the run if not so.", &
                          usage="CHECK T", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OVERWRITE", &
                          description="Specify this keyword to overwrite the output BQB file if"// &
                          " it already exists. By default, the data is appended to an existing file.", &
                          usage="OVERWRITE T", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="HISTORY", &
                          description="Controls how many previous steps are taken into account for extrapolation in"// &
                          " compression. Use a value of 1 to compress the frames independently.", &
                          usage="HISTORY 10", n_var=1, default_i_val=10, type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PARAMETER_KEY", &
                          description="Allows to supply previously optimized compression parameters via a"// &
                          " parameter key (alphanumeric character sequence starting with 'at')."// &
                          " Just leave away the 'at' sign here, because CP2k will otherwise"// &
                          " assume it is a variable name in the input", &
                          usage="PARAMETER_KEY <KEY>", n_var=1, default_c_val="", type_of_var=char_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OPTIMIZE", &
                          description="Controls the time spent to optimize the parameters for compression efficiency.", &
                          usage="OPTIMIZE {OFF,QUICK,NORMAL,PATIENT,EXHAUSTIVE}", repeats=.FALSE., n_var=1, &
                          default_i_val=bqb_opt_quick, &
                          enum_c_vals=s2a("OFF", "QUICK", "NORMAL", "PATIENT", "EXHAUSTIVE"), &
                          enum_desc=s2a("No optimization (use defaults)", "Quick optimization", &
                                        "Standard optimization", "Precise optimization", "Exhaustive optimization"), &
                          enum_i_vals=(/bqb_opt_off, bqb_opt_quick, bqb_opt_normal, bqb_opt_patient, bqb_opt_exhaustive/))
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      ! Voronoi Integration via LibVori
      CALL cp_print_key_section_create(print_key, __LOCATION__, name="VORONOI", &
                                       description="Controls the Voronoi integration of the total electron density"// &
                                       " for the computation of electromagnetic moments, see [Thomas2015],"// &
                                       " [Brehm2020], and [Brehm2021] (via LibVori,"// &
                                       " see <a href=""https://brehm-research.de/voronoi"" target=""_blank"">"// &
                                       " https://brehm-research.de/voronoi</a> ).", &
                                       print_level=debug_print_level + 1, filename="", &
                                       citations=(/Rycroft2009, Thomas2015, Brehm2018, Brehm2020, Brehm2021/))

      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="Appends frames to already existing .voronoi file.", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SANITY_CHECK", &
                          description="Performs a sanity check before each Voronoi integration, i.e.,"// &
                          " checks if every grid point is located in exactly one Voronoi cell.", &
                          usage="SANITY_CHECK T", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OVERWRITE", &
                          description="Specify this keyword to overwrite any existing ""properties.emp"" file if"// &
                          " it already exists. By default, the data is appended to an existing .emp file.", &
                          usage="OVERWRITE T", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SKIP_FIRST", &
                          description="Skips the first step of a MD run (avoids duplicate step if restarted).", &
                          usage="SKIP_FIRST T", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="VERBOSE", &
                          description="Switches on verbose screen output of the Voronoi integration.", &
                          usage="VERBOSE T", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OUTPUT_EMP", &
                          description="Writes the resulting electromagnetic moments to a binary file ""properties.emp""."// &
                          " The file name cannot be changed.", &
                          usage="OUTPUT_EMP T", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OUTPUT_TEXT", &
                          description="Writes the resulting electromagnetic moments to text files (*.voronoi)."// &
                          " The file name is specified via FILENAME.", &
                          usage="OUTPUT_TEXT T", default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="REFINEMENT_FACTOR", &
                          description="Sets the refinement factor for the Voronoi integration.", &
                          usage="REFINEMENT 2", n_var=1, default_i_val=1, type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="VORONOI_RADII", &
                          description="Which atomic radii to use for the radical Voronoi tessellation.", &
                          usage="VORONOI_RADII {UNITY,VDW,COVALENT,USER}", repeats=.FALSE., n_var=1, &
                          default_i_val=voro_radii_vdw, &
                          enum_c_vals=s2a("UNITY", "VDW", "COVALENT", "USER"), &
                          enum_desc=s2a("Use unity radii (non-radical Voronoi tessellation)", "Use VdW atom radii", &
                                        "Use covalent atom radii", "Use user-specified atom radii"), &
                          enum_i_vals=(/voro_radii_unity, voro_radii_vdw, voro_radii_cov, voro_radii_user/))
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="USER_RADII", &
                          description="Defines user atom radii for the radical Voronoi tessellation (one per atom).", &
                          usage="USER_RADII {real} {real} {real}", repeats=.FALSE., &
                          unit_str="angstrom", &
                          type_of_var=real_t, n_var=-1)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      !Printing of Moments
      CALL create_dipoles_section(print_key, "MOMENTS", high_print_level)
      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="MAX_MOMENT", &
         description="Maximum moment to be calculated. Values higher than 1 not implemented under periodic boundaries.", &
         usage="MAX_MOMENT {integer}", &
         repeats=.FALSE., &
         n_var=1, &
         type_of_var=integer_t, &
         default_i_val=1)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, &
                          name="MAGNETIC", &
                          description="Calculate also magnetic moments, only implemented without periodic boundaries", &
                          usage="MAGNETIC yes", &
                          repeats=.FALSE., &
                          n_var=1, &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

   END SUBROUTINE create_ec_print_section

END MODULE input_cp2k_ec
