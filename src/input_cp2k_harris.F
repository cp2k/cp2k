!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Harris input section
! **************************************************************************************************
MODULE input_cp2k_harris
   USE input_constants,                 ONLY: hden_atomic,&
                                              hden_cube,&
                                              hden_cube_fit,&
                                              hfun_harris,&
                                              horb_default
   USE input_keyword_types,             ONLY: keyword_create,&
                                              keyword_release,&
                                              keyword_type
   USE input_section_types,             ONLY: section_add_keyword,&
                                              section_create,&
                                              section_type
   USE input_val_types,                 ONLY: char_t
   USE kinds,                           ONLY: dp
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
                          n_keywords=11, n_subsections=0, repeats=.FALSE.)

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
                          enum_i_vals=[hfun_harris])
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DENSITY_SOURCE", &
                          description="Method to create the input density", &
                          usage="DENSITY_SOURCE ATOMIC", &
                          default_i_val=hden_atomic, &
                          enum_c_vals=s2a("ATOMIC", "CUBE", "CUBE_FIT"), &
                          enum_desc=s2a("Atomic densities", &
                                        "Electron density read directly from a Gaussian cube file", &
                                        "Cube density fitted by a constrained AO density matrix. "// &
                                        "This dense reference method currently requires GPW, a "// &
                                        "spin-restricted calculation, and the Gamma point"), &
                          enum_i_vals=[hden_atomic, hden_cube, hden_cube_fit])
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FILE_DENSITY", &
                          description="Filename of the electron density in Gaussian cube format. "// &
                          "Used when DENSITY_SOURCE is CUBE or CUBE_FIT. The cube grid has to coincide with "// &
                          "the finest CP2K real-space grid.", &
                          usage="FILE_DENSITY <FILENAME>", &
                          type_of_var=char_t, default_c_val="", n_var=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FIT_MAX_ITER", &
                          description="Maximum number of projected-gradient iterations for CUBE_FIT.", &
                          usage="FIT_MAX_ITER 50", default_i_val=50)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FIT_EPS", &
                          description="Target RMS error of the fitted density for CUBE_FIT.", &
                          usage="FIT_EPS 1.0E-3", default_r_val=1.0E-3_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FIT_STEP_SIZE", &
                          description="Initial spectral projected-gradient step size for CUBE_FIT.", &
                          usage="FIT_STEP_SIZE 1.0", default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FIT_MAX_BACKTRACK", &
                          description="Maximum number of line-search reductions in each CUBE_FIT iteration.", &
                          usage="FIT_MAX_BACKTRACK 20", default_i_val=20)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ORBITAL_BASIS", &
                          description="Specifies the type of basis to be used for the energy functional. ", &
                          default_i_val=horb_default, &
                          enum_c_vals=s2a("ATOMIC_KIND_BASIS"), &
                          enum_desc=s2a("Atomic kind orbital basis"), &
                          enum_i_vals=[horb_default])
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
