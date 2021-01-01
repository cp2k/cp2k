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
        ec_curvy_steps, ec_diagonalization, ec_functional_harris, ec_matrix_sign, ec_matrix_tc2, &
        ec_matrix_trs4, kg_cholesky, ls_cluster_atomic, ls_cluster_molecular, &
        ls_s_inversion_hotelling, ls_s_inversion_sign_sqrt, ls_s_preconditioner_atomic, &
        ls_s_preconditioner_molecular, ls_s_preconditioner_none, ls_s_sqrt_ns, ls_s_sqrt_proot, &
        ls_scf_sign_ns, ls_scf_sign_proot
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
                          n_keywords=0, n_subsections=1, repeats=.FALSE.)

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
                          enum_c_vals=s2a("DIAGONALIZATION", "CURVY_STEPS", &
                                          "MATRIX_SIGN", "TRS4", "TC2"), &
                          enum_desc=s2a("Diagonalization of KS matrix.", &
                                        "Head-Grodon curvy step algorithm", &
                                        "Matrix Sign algorithm", &
                                        "Trace resetting trs4 algorithm", &
                                        "Trace resetting tc2 algorithm"), &
                          enum_i_vals=(/ec_diagonalization, ec_curvy_steps, ec_matrix_sign, &
                                        ec_matrix_trs4, ec_matrix_tc2/))
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
                          usage="EPS_DEFAULT 1.0E-6 ", default_r_val=1.0E-12_dp)
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
      CALL section_create(section, __LOCATION__, name="LS_SOLVER", &
                          description="Specifies the parameters of the linear scaling solver routines", &
                          n_keywords=24, n_subsections=3, repeats=.FALSE., &
                          citations=(/VandeVondele2012/))

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="EPS_FILTER", &
                          description="Threshold used for filtering matrix operations.", &
                          usage="EPS_FILTER 1.0E-7", default_r_val=1.0E-6_dp)
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

      CALL keyword_create(keyword, __LOCATION__, name="EXTRAPOLATION_ORDER", &
                          description="Number of previous matrices used for the ASPC extrapolation of the initial guess. "// &
                          "0 implies that an atomic guess is used at each step. "// &
                          "low (1-2) will result in a drift of the constant of motion during MD. "// &
                          "high (>5) might be somewhat unstable, leading to more SCF iterations.", &
                          usage="EXTRAPOLATION_ORDER 3", default_i_val=4)
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
                          usage="REPORT_ALL_SPARSITIES", default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CHECK_S_INV", &
                          description="Perform an accuracy check on the inverse/sqrt of the s matrix.", &
                          usage="CHECK_S_INV", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

!     CALL create_ls_curvy_section(subsection)
!     CALL section_add_subsection(section, subsection)
!     CALL section_release(subsection)

!     CALL create_chebyshev_section(subsection)
!     CALL section_add_subsection(section, subsection)
!     CALL section_release(subsection)

   END SUBROUTINE create_ec_solver_section

END MODULE input_cp2k_ec
