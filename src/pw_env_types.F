!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/pw_env_types [1.0] *
!!
!!   NAME
!!     pw_env_types
!!
!!   FUNCTION
!!     container for various plainwaves related things
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     10.2002 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE pw_env_types
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_rs_pool_types,                ONLY: cp_rs_pool_p_type,&
                                             cp_rs_pool_type,&
                                             rs_pools_dealloc,&
                                             rs_pools_flush_cache
  USE cube_utils,                      ONLY: cube_info_type,&
                                             destroy_cube_info
  USE gaussian_gridlevels,             ONLY: destroy_gaussian_gridlevel,&
                                             gridlevel_info_type
  USE kinds,                           ONLY: dp
  USE l_utils,                         ONLY: destroy_l_info,&
                                             l_info_release,&
                                             l_info_type
  USE pw_pool_types,                   ONLY: pw_pool_p_type,&
                                             pw_pool_type,&
                                             pw_pools_dealloc,&
                                             pw_pools_flush_cache
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='pw_env_types'

  PUBLIC :: pw_env_type
  PUBLIC :: pw_env_retain, pw_env_release, pw_env_get, pw_env_flush_cache
!***
!****************************************************************************

!!****s* m/pw_env_types [1.0] *
!!
!!   NAME
!!     pw_env_types
!!
!!   FUNCTION
!!     contained for different pw related things
!!
!!   NOTES
!!     be careful in parallel nsmax for the rs_grids is choosen with
!!     multigrid in mind! (well not as of this writing but...)
!!
!!   ATTRIBUTES
!!     - pw_pools: pools for multigrids
!!     - rs_pools: pools for the realspace grids
!!     - gridlevel_info: gaussian gridlevel info
!!     - l_info: informations for l calc
!!     - cube_info: info on cube and max sphere size for multigrids
!!     - aux_bas: which multigrid is the auxiliar basis
!!     
!!     readonly attributes
!!     - auxbas_pw_pool: a pool that allocates grids in the auxiliary basis
!!     - auxbas_rs_pool: real space grid in the auxiliary basis, be careful
!!       in parallel nsmax is choosen with multigrid in mind!
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     10.2002 created [fawzi]
!!     04.2003 added rs grid pools [fawzi]
!!
!!   SOURCE
  !***************************************************************************
  TYPE pw_env_type
     INTEGER :: ref_count
     TYPE(pw_pool_p_type), DIMENSION(:), POINTER :: pw_pools
     TYPE(cp_rs_pool_p_type), DIMENSION(:), POINTER :: rs_pools
     TYPE(gridlevel_info_type), POINTER :: gridlevel_info
     TYPE(l_info_type), POINTER :: l_info ! right place for it?
     TYPE(l_info_type), POINTER :: l_info_rho0 ! right place for it?
     TYPE(cube_info_type), DIMENSION(:), POINTER :: cube_info
     INTEGER :: auxbas_grid
  END TYPE pw_env_type
!!***
  !****************************************************************************

CONTAINS

!!****f* pw_env_types/pw_env_get [1.0] *
!!
!!   NAME
!!     pw_env_get
!!
!!   SYNOPSIS
!!     Subroutine pw_env_get(pw_env, pw_pools, l_info, cube_info,&
!!         gridlevel_info, auxbas_pw_pool, auxbas_cube_info, auxbas_grid,&
!!         rs_pools, auxbas_rs_pool, error)
!!       Type(pw_env_type), Pointer:: pw_env
!!       Type(pw_pool_p_type), Dimension(:), Pointer, Optional:: pw_pools
!!       Type(l_info_type), Pointer, Optional:: l_info
!!       Type(cube_info_type), Dimension(:), Pointer, Optional:: cube_info
!!       Type(gridlevel_info_type), Pointer, Optional:: gridlevel_info
!!       Type(pw_pool_type), Pointer, Optional:: auxbas_pw_pool
!!       Type(cube_info_type), Pointer, Optional:: auxbas_cube_info
!!       Integer, Intent (OUT), Optional:: auxbas_grid
!!       Type(cp_rs_pool_p_type), Dimension(:), Pointer, Optional::&
!!         rs_pools
!!       Type(cp_rs_pool_type), Pointer, Optional:: auxbas_rs_pool
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pw_env_get
!!
!!   FUNCTION
!!     returns the various attributes of the pw env
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pw_env: the einvironment you want info about
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!     For the other see the attributes of pw_env_type
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     10.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE pw_env_get(pw_env, pw_pools,l_info,l_info_rho0, cube_info, gridlevel_info,&
     auxbas_pw_pool,auxbas_cube_info,auxbas_grid,rs_pools,auxbas_rs_pool,&
     error)
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: pw_pools
    TYPE(l_info_type), OPTIONAL, POINTER     :: l_info
    TYPE(l_info_type), OPTIONAL, POINTER     :: l_info_rho0
    TYPE(cube_info_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: cube_info
    TYPE(gridlevel_info_type), OPTIONAL, &
      POINTER                                :: gridlevel_info
    TYPE(pw_pool_type), OPTIONAL, POINTER    :: auxbas_pw_pool
    TYPE(cube_info_type), OPTIONAL, POINTER  :: auxbas_cube_info
    INTEGER, INTENT(out), OPTIONAL           :: auxbas_grid
    TYPE(cp_rs_pool_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: rs_pools
    TYPE(cp_rs_pool_type), OPTIONAL, POINTER :: auxbas_rs_pool
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_env_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(pw_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     IF (PRESENT(pw_pools)) pw_pools => pw_env%pw_pools
     IF (PRESENT(rs_pools)) rs_pools => pw_env%rs_pools
     IF (PRESENT(l_info)) l_info => pw_env%l_info
     IF (PRESENT(l_info_rho0)) l_info_rho0 => pw_env%l_info_rho0
     IF (PRESENT(cube_info)) cube_info => pw_env%cube_info
     IF (PRESENT(gridlevel_info)) gridlevel_info => pw_env%gridlevel_info
     IF (PRESENT(auxbas_pw_pool)) THEN
        auxbas_pw_pool => pw_env%pw_pools(pw_env%auxbas_grid)%pool
     END IF
     IF (PRESENT(auxbas_pw_pool)) THEN
        auxbas_pw_pool => pw_env%pw_pools(pw_env%auxbas_grid)%pool
     END IF
     IF (PRESENT(auxbas_rs_pool)) THEN
        auxbas_rs_pool => pw_env%rs_pools(pw_env%auxbas_grid)%pool
     END IF
     IF (PRESENT(auxbas_grid)) auxbas_grid = pw_env%auxbas_grid
  END IF
END SUBROUTINE pw_env_get
!***************************************************************************

!!****f* pw_env_types/pw_env_retain [1.0] *
!!
!!   NAME
!!     pw_env_retain
!!
!!   SYNOPSIS
!!     Subroutine pw_env_retain(pw_env, error)
!!       Type(pw_env_type), Pointer:: pw_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pw_env_retain
!!
!!   FUNCTION
!!     retains the pw_env (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pw_env: the pw_env to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     10.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE pw_env_retain(pw_env,error)
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(pw_env%ref_count>0,cp_failure_level,routineP,error,failure)
     pw_env%ref_count=pw_env%ref_count+1
  END IF
END SUBROUTINE pw_env_retain
!***************************************************************************

!!****f* pw_env_types/pw_env_release [1.0] *
!!
!!   NAME
!!     pw_env_release
!!
!!   SYNOPSIS
!!     Subroutine pw_env_release(pw_env, error)
!!       Type(pw_env_type), Pointer:: pw_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pw_env_release
!!
!!   FUNCTION
!!     releases the given pw_env (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pw_env: the pw_env to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     10.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE pw_env_release(pw_env, error)
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: igrid_level, stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (ASSOCIATED(pw_env)) THEN
     CPPrecondition(pw_env%ref_count>0,cp_failure_level,routineP,error,failure)
     pw_env%ref_count=pw_env%ref_count-1
     IF (pw_env%ref_count<1) THEN
        CALL pw_pools_dealloc(pw_env%pw_pools,error=error)
        IF (ASSOCIATED(pw_env%gridlevel_info)) THEN
           CALL destroy_gaussian_gridlevel(pw_env%gridlevel_info,error=error)
           DEALLOCATE(pw_env%gridlevel_info,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(pw_env%l_info)) THEN
           CALL l_info_release(pw_env%l_info)
           IF(pw_env%l_info%ref_count == 0) THEN
              DEALLOCATE(pw_env%l_info, stat=stat)
              CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
           END IF
        END IF
        IF (ASSOCIATED(pw_env%l_info_rho0)) THEN
           CALL l_info_release(pw_env%l_info_rho0)
           IF(pw_env%l_info_rho0%ref_count == 0) THEN
              DEALLOCATE(pw_env%l_info_rho0, stat=stat)
              CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
           END IF
        END IF
        IF (ASSOCIATED(pw_env%cube_info)) THEN
           DO igrid_level=1,SIZE(pw_env%cube_info)
              CALL destroy_cube_info(pw_env%cube_info(igrid_level))
           END DO
           DEALLOCATE(pw_env%cube_info,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        CALL rs_pools_dealloc(pw_env%rs_pools, error=error)
        DEALLOCATE(pw_env, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(pw_env)
END SUBROUTINE pw_env_release
!***************************************************************************

!!****f* pw_env_types/pw_env_flush_cache [1.0] *
!!
!!   NAME
!!     pw_env_flush_cache
!!
!!   SYNOPSIS
!!     Subroutine pw_env_flush_cache(pw_env, error)
!!       Type(pw_env_type), Pointer:: pw_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pw_env_flush_cache
!!
!!   FUNCTION
!!     flushes the cached pws
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pw_env: the pw_env to be flushed
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     11.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE pw_env_flush_cache(pw_env,error)
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_env_flush_cache', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     IF (ASSOCIATED(pw_env%pw_pools)) THEN
        CALL pw_pools_flush_cache(pw_env%pw_pools,error=error)
     END IF
     IF (ASSOCIATED(pw_env%rs_pools)) THEN
        CALL rs_pools_flush_cache(pw_env%rs_pools, error=error)
     END IF
  END IF
END SUBROUTINE pw_env_flush_cache
!***************************************************************************

END MODULE 
