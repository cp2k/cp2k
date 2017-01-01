!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

MODULE qs_fb_input
   USE bibliography,                    ONLY: Rayson2009
   USE cp_units,                        ONLY: cp_unit_to_cp2k
   USE input_keyword_types,             ONLY: keyword_create,&
                                              keyword_release,&
                                              keyword_type
   USE input_section_types,             ONLY: section_add_keyword,&
                                              section_create,&
                                              section_type
   USE input_val_types,                 ONLY: logical_t,&
                                              real_t
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_fb_input'

   PUBLIC :: create_filtermatrix_section

CONTAINS

! **************************************************************************************************
!> \brief Input section for filter matrix diagonalisation method
!> \param section : section to be created
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! **************************************************************************************************
   SUBROUTINE create_filtermatrix_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_filtermatrix_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))

      CALL section_create(section, "FILTER_MATRIX", &
                          description=" ", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)

      CALL keyword_create(keyword, &
                          name="FILTER_TEMPERATURE", &
                          description="Temperature used for the filter function used "// &
                          "to construct the filter matrix.", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=cp_unit_to_cp2k(value=10000.0_dp, &
                                                        unit_str="K"), &
                          unit_str="K", &
                          usage="FILTER_TEMPERATURE [K] 10000", &
                          citations=(/Rayson2009/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, &
                          name="AUTO_CUTOFF_SCALE", &
                          description="Scalar constant multiplied to maximum orbital "// &
                          "size of each atom, used for automatically "// &
                          "creating cutoff radii for atomic matrices", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=0.5_dp, &
                          usage="AUTO_CUTOFF_SCALE 0.5_dp", &
                          citations=(/Rayson2009/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, &
                          name="EPS_FB", &
                          description="Default tolerance used in generating the filter "// &
                          "matrix. Anything less than EPS_FB will be "// &
                          "regarded as zero", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=1.e-12_dp, &
                          usage="EPS_FB 1.e-12", &
                          citations=(/Rayson2009/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, &
                          name="COLLECTIVE_COMMUNICATION", &
                          description="If set to TRUE, then all MPI communications "// &
                          "required for the construction of the "// &
                          "filter matrix is done at the start and end "// &
                          "of each filter matrix calculation. This "// &
                          "makes communications more efficient, at "// &
                          "the expense of using more memory. If you "// &
                          "find the fb_fltrmat_add_blkcol_mpi times "// &
                          "at the end of CP2K output is high, then "// &
                          "run again with this option set to .TRUE.", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=logical_t, &
                          default_l_val=.FALSE., &
                          usage="COLLECTIVE_COMMUNICATION T")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_filtermatrix_section

END MODULE qs_fb_input
