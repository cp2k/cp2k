!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Input definition and setup for EEQ model
!> \author JGH [2024]
! **************************************************************************************************
MODULE eeq_input
   USE input_keyword_types,             ONLY: keyword_create,&
                                              keyword_release,&
                                              keyword_type
   USE input_section_types,             ONLY: section_add_keyword,&
                                              section_create,&
                                              section_type,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   TYPE eeq_solver_type
      LOGICAL                                :: direct = .FALSE.
      LOGICAL                                :: sparse = .FALSE.
      REAL(KIND=dp)                          :: eps_diis = 1.0E-09_dp
      REAL(KIND=dp)                          :: alpha = 0.75_dp
      INTEGER                                :: mdiis = 12
      INTEGER                                :: sdiis = 3
      INTEGER                                :: max_diis = 500
   END TYPE eeq_solver_type

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'eeq_input'

   PUBLIC :: eeq_solver_type
   PUBLIC :: read_eeq_param
   PUBLIC :: create_eeq_control_section

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_eeq_control_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="EEQ", &
                          description="Parameters needed for EEQ method and solver", &
                          n_keywords=1, n_subsections=1, repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="DIRECT", &
                          description="Use a direct method to solve the EEQ equations in PBC (matrix solver)", &
                          usage="DIRECT", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SPARSE", &
                          description="Use a sparse method to solve the EEQ equations. (NYA)", &
                          usage="SPARSE", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_DIIS", &
                          description="Accuracy for the iterative solver.", &
                          usage="EPS_DIIS 1.0E-10", default_r_val=1.0e-10_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ALPHA", &
                          description="Step length of initial steepest descent steps.", &
                          usage="ALPHA 1.0", default_r_val=0.75_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_DIIS", &
                          description="Max. number of iterations for EEQ solver.", &
                          usage="MAX_DIIS  100", default_i_val=500)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MDIIS", &
                          description="Max. number of DIIS vectors used.", &
                          usage="MDIIS 10", default_i_val=12)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SDIIS", &
                          description="Number of vectors accumulated before starting DIIS.", &
                          usage="SDIIS 4", default_i_val=3)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_eeq_control_section

! **************************************************************************************************
!> \brief ...
!> \param eeq_section ...
!> \param eeq_sparam ...
! **************************************************************************************************
   SUBROUTINE read_eeq_param(eeq_section, eeq_sparam)

      TYPE(section_vals_type), POINTER                   :: eeq_section
      TYPE(eeq_solver_type), INTENT(INOUT)               :: eeq_sparam

      CALL section_vals_val_get(eeq_section, "DIRECT", l_val=eeq_sparam%direct)
      CALL section_vals_val_get(eeq_section, "SPARSE", l_val=eeq_sparam%sparse)
      CALL section_vals_val_get(eeq_section, "EPS_DIIS", r_val=eeq_sparam%eps_diis)
      CALL section_vals_val_get(eeq_section, "ALPHA", r_val=eeq_sparam%alpha)
      CALL section_vals_val_get(eeq_section, "MAX_DIIS", i_val=eeq_sparam%max_diis)
      CALL section_vals_val_get(eeq_section, "MDIIS", i_val=eeq_sparam%mdiis)
      CALL section_vals_val_get(eeq_section, "SDIIS", i_val=eeq_sparam%sdiis)

   END SUBROUTINE read_eeq_param

END MODULE eeq_input
