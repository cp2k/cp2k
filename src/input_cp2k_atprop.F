!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief input section for atomic properties
!> \par History
!>      07.2011 created
!> \author JHU
! *****************************************************************************
MODULE input_cp2k_atprop
  USE bibliography,                    ONLY: Kikuchi2009
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_create,&
                                             section_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_atprop'

  PUBLIC :: create_atprop_section

CONTAINS

! *****************************************************************************
!> \brief Creates the ATOMIC section
!> \param section the section to create
!> \author JHU
! *****************************************************************************
  SUBROUTINE create_atprop_section(section)
    TYPE(section_type), POINTER              :: section

    CHARACTER(len=*), PARAMETER :: routineN = 'create_atprop_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,failure)
    CALL section_create(section,name="ATOMIC",&
         description="Controls the calculation of atomic properties. "//&
                     "Printing is controled by FORCE_EVAL / PRINT / PROGRAM_RUN_INFO",&
         repeats=.FALSE., &
         citations=(/Kikuchi2009/))

    NULLIFY(keyword)

    CALL keyword_create(keyword, name="ENERGY",&
         description="Calculate atomic energies ",&
         usage="ENERGY {logical}",&
         repeats=.FALSE.,&
         n_var=1,&
         default_l_val=.FALSE.,&
         lone_keyword_l_val=.TRUE.)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="PRESSURE",&
         description="Calculate atomic pressure tensors ",&
         usage="PRESSURE {logical}",&
         repeats=.FALSE.,&
         n_var=1,&
         default_l_val=.FALSE.,&
         lone_keyword_l_val=.TRUE.)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

  END SUBROUTINE create_atprop_section

END MODULE input_cp2k_atprop
