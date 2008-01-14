!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****h* cp2k/input_cp2k_rsgrid *
!!
!!   NAME
!!     input_cp2k_rsgrid
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Joost
!!
!!   MODIFICATION HISTORY
!!     01.2008 Created
!!
!!   SOURCE
!****************************************************************************
MODULE input_cp2k_rsgrid
  USE f77_blas
  USE input_constants
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_add_subsection,&
                                             section_create,&
                                             section_release,&
                                             section_type
  USE kinds,                           ONLY: dp
  USE string_utilities,                ONLY: s2a
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_rsgrid'

PUBLIC :: create_rsgrid_section

!***
!****************************************************************************
CONTAINS

!***************************************************************************
!!****f* input_cp2k_rsgrid/create_rsgrid_section *
!!
!!   NAME
!!     create_rsgrid_section
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Joost
!!
!!*** **********************************************************************
  SUBROUTINE create_rsgrid_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_rsgrid_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN  
       CALL section_create(section,name="RS_GRID",&
            description="Set options that influence how the realspace grids are being distributed in parallel runs.",&
            n_keywords=5, n_subsections=0, repeats=.TRUE., required=.FALSE.,error=error)

       CALL keyword_create(keyword, name="DISTRIBUTION_TYPE",&
            description="Parallelization strategy.",&
            usage="DISTRIBUTION_TYPE DISTRIBUTED",&
            enum_c_vals=s2a("AUTOMATIC","DISTRIBUTED","REPLICATED"),&
            enum_i_vals=(/rsgrid_automatic,rsgrid_distributed,rsgrid_replicated /),&
            enum_desc=s2a("Use heuristic rules to decide between distributed and replicated", &
                          "Force a distributed setup if possible",&
                          "Force a replicated setup"),&
            default_i_val=rsgrid_automatic, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)        

       CALL keyword_create(keyword, name="DISTRIBUTION_LAYOUT",&
            description="Specifies the number of slices in the x, y and z directions."//&
                        "-1 specifies that any number of slices is OK."//&
                        "If a given distribution can not be satisfied, a replicated grid will result.",&
            usage="DISTRIBUTION_LAYOUT",&
            repeats=.FALSE.,required=.FALSE.,n_var=3,&
            default_i_vals=(/-1,-1,-1/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


    END IF
  END SUBROUTINE create_rsgrid_section

END MODULE input_cp2k_rsgrid
