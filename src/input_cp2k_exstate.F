!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Excited state input section
!> \par History
!>      01.2020 created
!> \author jgh
! **************************************************************************************************
MODULE input_cp2k_exstate
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

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_exstate'

   PUBLIC :: create_exstate_section

CONTAINS

! **************************************************************************************************
!> \brief creates the EXCITED ENERGY section
!> \param section ...
!> \author JGH
! **************************************************************************************************
   SUBROUTINE create_exstate_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_exstate_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))

      NULLIFY (keyword)
      CALL section_create(section, __LOCATION__, name="EXCITED_STATES", &
                          description="Sets the various options for Excited State Potential Energy Calculations", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="Controls the activation of the excited states", &
                          usage="&EXCITED_STATES T", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="STATE", &
                          description="Excited state to be used in calculation.", &
                          usage="STATE 2", &
                          default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_exstate_section

END MODULE input_cp2k_exstate
