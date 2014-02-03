!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      01.2008 Created
!> \author Joost
! *****************************************************************************
MODULE input_cp2k_rsgrid
  USE f77_blas
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_create,&
                                             section_type
  USE kinds,                           ONLY: dp
  USE realspace_grid_types,            ONLY: rsgrid_automatic,&
                                             rsgrid_distributed,&
                                             rsgrid_replicated
  USE string_utilities,                ONLY: s2a
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_rsgrid'

PUBLIC :: create_rsgrid_section

!***
CONTAINS

! *****************************************************************************
!> \author Joost
! *****************************************************************************
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

       NULLIFY(keyword)
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
                        "If a given distribution can not be satisfied, a replicated grid will result."//&
                        "Also see LOCK_DISTRIBUTION.",&
            usage="DISTRIBUTION_LAYOUT",&
            repeats=.FALSE.,required=.FALSE.,n_var=3,&
            default_i_vals=(/-1,-1,-1/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAX_DISTRIBUTED_LEVEL",&
            description="If the multigrid-level of a grid is larger than the parameter,"//&
                        " it will not be distributed in the automatic scheme.",&
            usage="MAX_DISTRIBUTED_LEVEL 1",&
            default_i_val=2, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="LOCK_DISTRIBUTION",&
            description="Expert use only, only basic QS deals correctly with a non-default value."//&
                        "If the distribution is locked, a grid will have the same distribution as"//&
                        "the next finer multigrid (provided it is distributed)."//&
                        "If unlocked, all grids can be distributed freely.",&
            usage="LOCK_DISTRIBUTION TRUE",&
            default_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MEMORY_FACTOR",&
            description="A grid will only be distributed if the memory usage for that grid (including halo) "//&
                        "is smaller than a replicated grid by this parameter.",&
            usage="MEMORY_FACTOR 4.0",&
            default_r_val=2.0_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="HALO_REDUCTION_FACTOR",&
            description="Can be used to reduce the halo of the distributed grid (experimental features).",&
            usage="HALO_REDUCTION_FACTOR 0.5",&
            default_r_val=1.0_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
    END IF
  END SUBROUTINE create_rsgrid_section

END MODULE input_cp2k_rsgrid
