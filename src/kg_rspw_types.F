!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C ) 2002 - 2004  CP2K developers group                         !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****** cp2k/kg_rspw_types [1.0] *
!!
!!   NAME
!!     kg_rspw_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     gloria,30.10.2004, previously in kg_rho_types 
!!
!!   MODIFICATION HISTORY
!!     JGH (22-Feb-03) PW grid options added
!!     gt 16-nov-03 moved initialization in this new module
!!
!!   SOURCE
!******************************************************************************

MODULE kg_rspw_types
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_rs_pool_types,                ONLY: cp_rs_pool_p_type,&
                                             cp_rs_pool_type,&
                                             rs_pools_dealloc
  USE cube_utils,                      ONLY: cube_info_type,&
                                             destroy_cube_info
  USE gaussian_gridlevels,             ONLY: destroy_gaussian_gridlevel,&
                                             gridlevel_info_type
  USE kinds,                           ONLY: dbl
  USE l_utils,                         ONLY: destroy_l_info,&
                                             l_info_type
  USE pw_pool_types,                   ONLY: pw_pool_release,&
                                             pw_pool_p_type,&
                                             pw_pool_type,&
                                             pw_pools_dealloc
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='kg_rspw_types'
  PUBLIC :: kg_rspw_type, kg_rspw_release, kg_rspw_get, &
            kg_rspw_retain

  TYPE kg_rspw_type
     INTEGER :: auxbas_grid
     INTEGER :: ref_count
     TYPE(gridlevel_info_type),POINTER              :: gridlevel_info
     TYPE(l_info_type),        POINTER              :: l_info
     TYPE(cube_info_type), DIMENSION(:), POINTER    :: cube_info
     TYPE(pw_pool_p_type), DIMENSION(:), POINTER    :: pw_pools
     TYPE(cp_rs_pool_p_type), DIMENSION(:), POINTER :: rs_pools
  END TYPE kg_rspw_type
!-----------------------------------------------------------------------------!

CONTAINS

!-----------------------------------------------------------------------------!
!!***
!******************************************************************************
SUBROUTINE kg_rspw_get(kg_rspw, l_info, cube_info, gridlevel_info,&
     auxbas_pw_pool,auxbas_grid,rs_pools,pw_pools,auxbas_rs_pool,&
     error)
    TYPE(kg_rspw_type), POINTER               :: kg_rspw
    TYPE(pw_pool_p_type), DIMENSION(:), &
     OPTIONAL, POINTER                      :: pw_pools
    TYPE(l_info_type), OPTIONAL, POINTER     :: l_info
    TYPE(cube_info_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: cube_info
    TYPE(gridlevel_info_type), OPTIONAL, &
      POINTER                                :: gridlevel_info
    TYPE(pw_pool_type), OPTIONAL, POINTER    :: auxbas_pw_pool
    INTEGER, INTENT(out), OPTIONAL           :: auxbas_grid
    TYPE(cp_rs_pool_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: rs_pools
    TYPE(cp_rs_pool_type), OPTIONAL, POINTER :: auxbas_rs_pool
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_rspw_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(kg_rspw),cp_failure_level,routineP,error,failure)
  CPPrecondition(kg_rspw%ref_count>0,cp_failure_level,routineP,error,failure)

  IF (.NOT. failure) THEN
     IF (PRESENT(pw_pools)) pw_pools => kg_rspw%pw_pools
     IF (PRESENT(rs_pools)) rs_pools => kg_rspw%rs_pools
     IF (PRESENT(l_info)) l_info =>  kg_rspw%l_info
     IF (PRESENT(cube_info)) cube_info => kg_rspw%cube_info
     IF (PRESENT(gridlevel_info)) gridlevel_info => kg_rspw%gridlevel_info
     IF (PRESENT(auxbas_pw_pool)) THEN
        auxbas_pw_pool => kg_rspw%pw_pools(kg_rspw%auxbas_grid)%pool
     END IF
     IF (PRESENT(auxbas_rs_pool)) THEN
        auxbas_rs_pool => kg_rspw%rs_pools(kg_rspw%auxbas_grid)%pool
     END IF
     IF (PRESENT(auxbas_grid)) auxbas_grid = kg_rspw%auxbas_grid
  END IF
END SUBROUTINE kg_rspw_get
!***************************************************************************
SUBROUTINE kg_rspw_retain(kg_rspw,error)
    TYPE(kg_rspw_type), POINTER               :: kg_rspw
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_rspw_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
 
  CPPrecondition(ASSOCIATED(kg_rspw),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(kg_rspw%ref_count>0,cp_failure_level,routineP,error,failure)
     kg_rspw%ref_count=kg_rspw%ref_count+1
  END IF
END SUBROUTINE kg_rspw_retain
!***************************************************************************

!!****f* kg_force/release_kg_rspw [1.0] *
!!
!!   NAME
!!     release_kg_rspw
!!
!!   FUNCTION
!!     releases the memory used by the kg_rspw and kg_rspw 
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************
  SUBROUTINE kg_rspw_release(kg_rspw,error)               


    TYPE(kg_rspw_type), POINTER               :: kg_rspw
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_rspw_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, igrid_level, stat
    LOGICAL                                  :: failure

!---------------------------------------------------------------------------

   CALL timeset(routineN,"I","",handle)
   failure=.FALSE.  
   IF (ASSOCIATED(kg_rspw)) THEN
     CPPrecondition(kg_rspw%ref_count>0,cp_failure_level,routineP,error,failure)
     kg_rspw%ref_count=kg_rspw%ref_count-1
     IF (kg_rspw%ref_count<1) THEN
       IF (ASSOCIATED(kg_rspw%gridlevel_info)) THEN
         CALL destroy_gaussian_gridlevel(kg_rspw%gridlevel_info,error=error)
         DEALLOCATE(kg_rspw%gridlevel_info,stat=stat)
         CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
       IF (ASSOCIATED(kg_rspw%l_info)) THEN
         CALL destroy_l_info(kg_rspw%l_info)
         DEALLOCATE(kg_rspw%l_info, stat=stat)
         CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
       IF (ASSOCIATED(kg_rspw%cube_info)) THEN
         DO igrid_level=1,SIZE(kg_rspw%cube_info)
           CALL destroy_cube_info(kg_rspw%cube_info(igrid_level))
         END DO
        DEALLOCATE (kg_rspw%cube_info,STAT=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF   
      CALL pw_pools_dealloc(kg_rspw%pw_pools,error=error)
      CALL rs_pools_dealloc(kg_rspw%rs_pools,error=error)
      DEALLOCATE(kg_rspw, stat=stat)
      CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
  END IF
  NULLIFY(kg_rspw)
  CALL timestop(0.0_dbl,handle)
  END SUBROUTINE kg_rspw_release

!***************************************************************************

END MODULE kg_rspw_types

!******************************************************************************
