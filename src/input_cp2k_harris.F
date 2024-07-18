!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Harris input section
! **************************************************************************************************
MODULE input_cp2k_harris
   USE input_constants,                 ONLY: hden_atomic,&
                                              hfun_harris,&
                                              horb_default
   USE input_keyword_types,             ONLY: keyword_create,&
                                              keyword_release,&
                                              keyword_type
   USE input_section_types,             ONLY: section_add_keyword,&
                                              section_create,&
                                              section_type
   USE string_utilities,                ONLY: s2a
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_harris'

   PUBLIC :: create_harris_section

CONTAINS

! **************************************************************************************************
!> \brief creates the HARRIS_METHOD section
!> \param section ...
!> \author JGH
! **************************************************************************************************
   SUBROUTINE create_harris_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))

      NULLIFY (keyword)
      CALL section_create(section, __LOCATION__, name="HARRIS_METHOD", &
                          description="Sets the various options for the Harris method", &
                          n_keywords=5, n_subsections=0, repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="Controls the activation of the Harris method", &
                          usage="&HARRIS_METHOD T", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ENERGY_FUNCTIONAL", &
                          description="Functional used in energy correction", &
                          usage="ENERGY_FUNCTIONAL HARRIS", &
                          default_i_val=hfun_harris, &
                          enum_c_vals=s2a("HARRIS"), &
                          enum_desc=s2a("Harris functional"), &
                          enum_i_vals=(/hfun_harris/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DENSITY_SOURCE", &
                          description="Method to create the input density", &
                          usage="DENSITY_SOURCE ATOMIC", &
                          default_i_val=hden_atomic, &
                          enum_c_vals=s2a("ATOMIC"), &
                          enum_desc=s2a("Atomic densities"), &
                          enum_i_vals=(/hden_atomic/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ORBITAL_BASIS", &
                          description="Specifies the type of basis to be used for the energy functional. ", &
                          default_i_val=horb_default, &
                          enum_c_vals=s2a("ATOMIC_KIND_BASIS"), &
                          enum_desc=s2a("Atomic kind orbital basis"), &
                          enum_i_vals=(/horb_default/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DEBUG_FORCES", &
                          description="Additional output to debug Harris method forces.", &
                          usage="DEBUG_FORCES T", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="DEBUG_STRESS", &
                          description="Additional output to debug Harris method stress.", &
                          usage="DEBUG_STRESS T", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_harris_section

END MODULE input_cp2k_harris
