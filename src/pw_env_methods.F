!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002 -2003 CP2K developers group                            !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/pw_env_methods [1.0] *
!!
!!   NAME
!!     pw_env_methods
!!
!!   FUNCTION
!!     methods of pw_env that have dependence on qs_env
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     10.2002 created [fawzi]
!!     JGH (22-Feb-03) PW grid options added
!!     04.2003 added rs grid pools [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE pw_env_methods
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_rs_pool_types,                ONLY: cp_rs_pool_p_type,&
                                             cp_rs_pool_type,&
                                             rs_pool_create,&
                                             rs_pool_release,&
                                             rs_pool_retain,&
                                             rs_pools_copy,&
                                             rs_pools_dealloc
  USE cube_utils,                      ONLY: cube_info_type,&
                                             destroy_cube_info,&
                                             init_cube_info
  USE dft_types,                       ONLY: dft_control_type
  USE gaussian_gridlevels,             ONLY: destroy_gaussian_gridlevel,&
                                             gridlevel_info_type,&
                                             init_gaussian_gridlevel
  USE global_types,                    ONLY: global_environment_type
  USE greens_fn,                       ONLY: greens_function_type,&
                                             pw_green_fn_init
  USE kinds,                           ONLY: wp=>dp
  USE l_utils,                         ONLY: destroy_l_info,&
                                             init_l_info,&
                                             l_info_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_release,&
                                             pw_env_retain,&
                                             pw_env_type
  USE pw_grid_types,                   ONLY: FULLSPACE,&
                                             HALFSPACE,&
                                             pw_grid_type
  USE pw_grids,                        ONLY: pw_find_cutoff,&
                                             pw_grid_construct,&
                                             pw_grid_release,&
                                             pw_grid_setup
  USE pw_pool_types,                   ONLY: pw_pool_create,&
                                             pw_pool_init_coeff,&
                                             pw_pool_p_type,&
                                             pw_pool_release,&
                                             pw_pool_retain,&
                                             pw_pool_type,&
                                             pw_pools_copy,&
                                             pw_pools_dealloc
  USE pws,                             ONLY: init_pw_poisson_solver
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_mo_types,                     ONLY: get_mo_set
  USE simulation_cell,                 ONLY: cell_type
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='pw_env_methods'

  PUBLIC :: pw_env_create, pw_env_rebuild
!***
!****************************************************************************
CONTAINS

!!****f* pw_env_methods/pw_env_create [1.0] *
!!
!!   NAME
!!     pw_env_create
!!
!!   SYNOPSIS
!!     Subroutine pw_env_create(pw_env, qs_env, error)
!!       Type(pw_env_type), Pointer:: pw_env
!!       Type(qs_environment_type), Intent (IN), Target, Optional:: qs_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pw_env_create
!!
!!   FUNCTION
!!     creates a pw_env, if qs_env is given calls pw_env_rebuild
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pw_env: the pw_env that gets created
!!     - qs_env: if given it is used to initialize the pw_env
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
SUBROUTINE pw_env_create(pw_env,qs_env,error)
  TYPE(pw_env_type), POINTER :: pw_env
  TYPE(qs_environment_type), INTENT(in), TARGET, OPTIONAL :: qs_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='pw_env_create',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat
  failure=.FALSE.
  
  ALLOCATE(pw_env, stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     NULLIFY(pw_env%pw_pools, pw_env%gridlevel_info,&
          pw_env%l_info, pw_env%cube_info, pw_env%rs_pools)
     pw_env%auxbas_grid=-1
     pw_env%ref_count=1
     IF (PRESENT(qs_env)) CALL pw_env_rebuild(pw_env,qs_env=qs_env,error=error)
  END IF
END SUBROUTINE pw_env_create
!***************************************************************************

!!****f* pw_env_methods/pw_env_rebuild [1.0] *
!!
!!   NAME
!!     pw_env_rebuild
!!
!!   SYNOPSIS
!!     Subroutine pw_env_rebuild(pw_env, qs_env, error)
!!       Type(pw_env_type), Pointer:: pw_env
!!       Type(qs_environment_type), Intent (IN), Target:: qs_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pw_env_rebuild
!!
!!   FUNCTION
!!     rebuilds the pw_env data (necessary at the beginning or if the
!!     cell or cutoffs change)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pw_env: the environment to rebuild
!!     - qs_env: the qs_env where to get the cell, cutoffs,...
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
SUBROUTINE pw_env_rebuild(pw_env, qs_env, error)
  TYPE(pw_env_type), POINTER :: pw_env
  TYPE(qs_environment_type), INTENT(in), TARGET :: qs_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  INTEGER, PARAMETER :: maxgridpoints = 100 ! should be be avoided (see later)
  LOGICAL :: failure, pw_grid_information
  CHARACTER(len=*), PARAMETER :: routineN='pw_env_rebuild',&
       routineP=moduleN//':'//routineN
  REAL(wp)          :: rel_cutoff, maxradius, cutilev
  INTEGER :: handle,i,igrid_level,stat, ngrid_level, cmax,cmaxl, istat,maxl
  INTEGER, DIMENSION(3) :: nsmax
  REAL(wp), DIMENSION(:), POINTER :: cutoff
  TYPE(cell_type), POINTER          :: cell
  TYPE(pw_grid_type), POINTER :: pw_grid, old_pw_grid
  TYPE(dft_control_type), POINTER :: dft_control
  TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
  TYPE(pw_pool_p_type),DIMENSION(:), POINTER :: pw_pools
  TYPE(cp_rs_pool_p_type),DIMENSION(:), POINTER :: rs_pools
  TYPE(cp_para_env_type), POINTER :: para_env
  TYPE(cp_logger_type), POINTER :: logger

  CALL timeset(routineN,"I","",handle)

  failure=.FALSE.
  NULLIFY(cutoff,cell,pw_grid,old_pw_grid,dft_control,atomic_kind_set,&
       pw_pools,rs_pools, para_env)
  logger => cp_error_get_logger(error)
  pw_grid_information=logger%print_keys%pw_grid_information
 
  CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(pw_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL pw_pools_dealloc(pw_env%pw_pools,error=error)
     CALL rs_pools_dealloc(pw_env%rs_pools,error=error)
     IF (ASSOCIATED(pw_env%gridlevel_info)) THEN
        CALL destroy_gaussian_gridlevel(pw_env%gridlevel_info,error=error)
     ELSE
        ALLOCATE(pw_env%gridlevel_info,stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     END IF
     IF (ASSOCIATED(pw_env%l_info)) THEN
        CALL destroy_l_info(pw_env%l_info)
     ELSE
        ALLOCATE(pw_env%l_info,stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     END IF
     IF (ASSOCIATED(pw_env%cube_info)) THEN
        DO igrid_level=1,SIZE(pw_env%cube_info)
           CALL destroy_cube_info(pw_env%cube_info(igrid_level))
        END DO
        DEALLOCATE(pw_env%cube_info,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
     NULLIFY(pw_env%pw_pools, pw_env%cube_info)
  END IF
  IF (.NOT.failure) THEN
     CALL get_qs_env(qs_env=qs_env,&
          dft_control=dft_control,atomic_kind_set=atomic_kind_set,&
          cell=cell, para_env=para_env)
     
     rel_cutoff = dft_control%qs_control%relative_cutoff
     cutoff => dft_control%qs_control%e_cutoff
     ngrid_level = SIZE(cutoff)
     
     ! init gridlevel_info
     CALL init_gaussian_gridlevel(pw_env%gridlevel_info,&
          ngrid_levels=ngrid_level,cutoff=cutoff,rel_cutoff=rel_cutoff)

     ! init pw_grids and pools
     ALLOCATE(pw_pools(ngrid_level),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     ALLOCATE(rs_pools(ngrid_level),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

     DO igrid_level=1,ngrid_level ! keep the pw_grids in a separate array?
        ALLOCATE(pw_grid,stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        CALL pw_grid_construct(pw_grid)
        pw_grid%para%rs_dims ( 1 ) = para_env%num_pe
        pw_grid%para%rs_dims ( 2 ) = 1
        
        IF ( dft_control % qs_control % pw_grid_opt % spherical ) THEN
           cutilev = cutoff(igrid_level)
           pw_grid%grid_span = HALFSPACE
        ELSE IF ( dft_control % qs_control % pw_grid_opt % fullspace ) THEN
           cutilev = -cutoff(igrid_level)
           pw_grid%grid_span = FULLSPACE
        ELSE
           cutilev = -cutoff(igrid_level)
           pw_grid%grid_span = HALFSPACE
        END IF

        IF (igrid_level == 1) THEN
           IF (pw_grid_information) THEN
              CALL pw_grid_setup(cell,pw_grid,cutilev,&
                   pe_group=para_env%group,&
                   info= cp_logger_get_default_unit_nr(logger,&
                   local=.FALSE.), blocked=.FALSE.)
           ELSE
              CALL pw_grid_setup(cell,pw_grid,cutilev,&
                   pe_group=para_env%group, blocked=.FALSE.)
           END IF
           old_pw_grid => pw_grid
        ELSE
           IF (pw_grid_information) THEN
              CALL pw_grid_setup(cell,pw_grid,cutilev,&
                   pe_group=para_env%group,&
                   info=cp_logger_get_default_unit_nr(logger,&
                   local=.FALSE.), blocked=.FALSE., &
                   ref_grid=old_pw_grid)
           ELSE
              CALL pw_grid_setup(cell,pw_grid,cutilev,&
                   pe_group=para_env%group, blocked=.FALSE., &
                   ref_grid=old_pw_grid)
           END IF
        END IF
        
        ! init pw_pools
        NULLIFY(pw_pools(igrid_level)%pool)
        CALL pw_pool_create(pw_pools(igrid_level)%pool,&
             pw_grid=pw_grid,error=error)

        ! init rs_pools
        ! nsmax hard coded for the time being
        ! gives the minimum thickness of a slice (or related quantity)
        ! is not yet used in rs_grid_setup
        ! this is dependent from the progression factor of the multigrids
        ! and the softest gaussian
        nsmax(:)=1000
        NULLIFY(rs_pools(igrid_level)%pool)
        CALL rs_pool_create(rs_pools(igrid_level)%pool,&
             el_struct=pw_grid,nsmax=nsmax,error=error)

        CALL pw_grid_release(pw_grid,error=error)
     END DO
     pw_env%pw_pools => pw_pools
     pw_env%rs_pools => rs_pools

     ! init cube info
     ALLOCATE (pw_env%cube_info(ngrid_level),STAT=istat)
     IF (istat /= 0) CALL stop_memory(routineP,"cube_info",0)

     cmax = 0
     DO igrid_level=1,ngrid_level
        maxradius = maxgridpoints*MAXVAL(pw_pools(igrid_level)%pool%pw_grid%dr)
        CALL init_cube_info(pw_env%cube_info(igrid_level),&
             pw_pools(igrid_level)%pool%pw_grid%dr(:),&
             maxradius,cmaxl)
        cmax = MAX(cmax,cmaxl)
     END DO

     ! init l_info
     CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
          maxl=maxl)
     ! maxl+1 to allow for the forces
     CALL init_l_info(pw_env%l_info,maxl+1,cmax)
     
     ! init auxbas_grid
     DO i=1,ngrid_level
        IF (cutoff(i) == dft_control%qs_control%cutoff) pw_env%auxbas_grid=i
     END DO

     ! init green's function
     CALL pw_green_fn_init(cell%green,&
          pw_env%pw_pools(pw_env%auxbas_grid)%pool%pw_grid)
  END IF
  CALL timestop(0.0_wp,handle)
END SUBROUTINE pw_env_rebuild
!***************************************************************************

END MODULE pw_env_methods
