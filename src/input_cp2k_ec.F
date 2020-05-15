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
                                              VandeVondele2012
   USE input_constants,                 ONLY: &
        ec_diagonalization, ec_functional_harris, ec_ls_solver, ec_matrix_sign, ec_matrix_tc2, &
        ec_matrix_trs4, ec_mo_solver, ec_ot_atomic, ec_ot_diag, ec_ot_gs, kg_cholesky, &
        ls_cluster_atomic, ls_cluster_molecular, ls_s_inversion_hotelling, &
        ls_s_inversion_sign_sqrt, ls_s_preconditioner_atomic, ls_s_preconditioner_molecular, &
        ls_s_preconditioner_none, ls_s_sqrt_ns, ls_s_sqrt_proot, ls_scf_sign_ns, &
        ls_scf_sign_proot, ot_precond_full_all, ot_precond_full_kinetic, ot_precond_full_single, &
        ot_precond_full_single_inverse, ot_precond_none, ot_precond_s_inverse, precond_mlp
   USE input_cp2k_xc,                   ONLY: create_xc_section
   USE input_keyword_types,             ONLY: keyword_create,&
                                              keyword_release,&
                                              keyword_type
   USE input_section_types,             ONLY: section_add_keyword,&
                                              section_add_subsection,&
                                              section_create,&
                                              section_release,&
                                              section_type
   USE input_val_types,                 ONLY: char_t
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

END MODULE input_cp2k_ec
