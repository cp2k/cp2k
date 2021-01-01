!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief manage control variables for the maximum overlap method
! **************************************************************************************************
MODULE qs_mom_types
   USE bibliography,                    ONLY: Barca2018,&
                                              Gilbert2008
   USE input_constants,                 ONLY: momproj_norm,&
                                              momproj_sum,&
                                              momtype_imom,&
                                              momtype_mom
   USE input_keyword_types,             ONLY: keyword_create,&
                                              keyword_release,&
                                              keyword_type
   USE input_section_types,             ONLY: section_add_keyword,&
                                              section_create,&
                                              section_type
   USE input_val_types,                 ONLY: integer_t
   USE string_utilities,                ONLY: s2a
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_mom_types'

   PUBLIC :: create_mom_section

CONTAINS

! **************************************************************************************************
!> \brief Create CP2K input section for variable occupancy using the Maximum
!>        Overlap Method. Only with diagonalization methods, i.e. not with OT
!> \param section section to create
!> \date       20.06.2013
!> \par History
!>      06.2013 created [MattW]
!>      01.2016 (DE)OCC_ALPHA and (DE)OCC_BETA keywords accept a list of
!>              molecular orbitals. Added two extra keywords: START_ITER
!>              and PROJ_FORMULA [Sergey Chulkov]
!> \author     MattW
!> \version    1.0
! **************************************************************************************************
   SUBROUTINE create_mom_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))

      CALL section_create(section, __LOCATION__, &
                          name="MOM", &
                          description="Define type and parameters for the maximum overlap method (MOM) "// &
                          "to determine orbital occupancies. "// &
                          "The MOM procedures activated by this section are only active for diagonalization "// &
                          "methods, i.e. not with minimization methods based on OT. "// &
                          "Incompatible with core-level excitation spectra (XAS).", &
                          n_keywords=8, n_subsections=0, repeats=.FALSE., &
                          citations=(/Gilbert2008, Barca2018/))

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="_SECTION_PARAMETERS_", &
                          description="Controls the activation of the MOM procedure", &
                          usage="MOM ON", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MOM_TYPE", &
                          description="Revision of the maximum overlap method to be used", &
                          usage="MOM_TYPE mom", default_i_val=momtype_mom, &
                          enum_c_vals=s2a("MOM", "IMOM"), &
                          enum_desc=s2a( &
                          "Originally proposed MOM protocol which uses molecular orbitals"// &
                          " from the previous SCF cycle as reference", &
                          "Initial-MOM which uses molecular orbitals of the initial guess as reference"), &
                          enum_i_vals=(/momtype_mom, momtype_imom/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="START_ITER", &
                          description="SCF iteration cycle to start the MOM procedure. "// &
                          "Could be used for ground state calculations only "// &
                          "in order to stabilise oscillating SCF iterations", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=integer_t, &
                          default_i_val=0, &
                          usage="START_ITER 2")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="DEOCC_ALPHA", &
                          description="Alpha orbitals to be deoccupied", &
                          repeats=.FALSE., &
                          n_var=-1, &
                          type_of_var=integer_t, &
                          default_i_val=0, &
                          usage="DEOCC_ALPHA 10 11 ...")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="DEOCC_BETA", &
                          description="Beta orbitals to be deoccupied", &
                          repeats=.FALSE., &
                          n_var=-1, &
                          type_of_var=integer_t, &
                          default_i_val=0, &
                          usage="DEOCC_BETA 10 11 ...")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="OCC_ALPHA", &
                          description="Alpha orbitals to be occupied", &
                          repeats=.FALSE., &
                          n_var=-1, &
                          type_of_var=integer_t, &
                          default_i_val=0, &
                          usage="OCC_ALPHA 12 15 ...")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="OCC_BETA", &
                          description="Beta orbitals to be occupied", &
                          repeats=.FALSE., &
                          n_var=-1, &
                          type_of_var=integer_t, &
                          default_i_val=0, &
                          usage="OCC_BETA 12 15 ...")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PROJ_FORMULA", &
                          description="Projection formula to be used", &
                          usage="PROJ_FORMULA norm", default_i_val=momproj_norm, &
                          enum_c_vals=s2a("NORM", "SUM"), &
                          enum_desc=s2a( &
                          "The one which ignores the phase of molecular orbitals: proj_j = \sqrt(\sum_i overlap(i, j)^2)", &
                          "The one proposed in the original paper: proj_j = |\sum_i overlap(i, j)|"), &
                          enum_i_vals=(/momproj_norm, momproj_sum/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_mom_section

END MODULE qs_mom_types
