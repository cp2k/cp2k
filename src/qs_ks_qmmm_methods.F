
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2003  CP2K developers group                            !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****** cp2k/qs_ks_qmmm_methods [1.0] *
!!
!!   NAME
!!     qs_ks_qmmm_methods
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     T.Laino
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE qs_ks_qmmm_methods
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_env_retain,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qmmm_types,                      ONLY: qmmm_env_qm_type
  USE cube_utils,                      ONLY: init_cube_info,&
                                             cube_info_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type, pw_env_retain
  USE pw_pool_types,                   ONLY: pw_pool_give_back_coeff,&
                                             pw_pool_init_coeff,&
                                             pw_pool_p_type,&
                                             pw_pool_retain,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type
  USE qs_ks_qmmm_types,                ONLY: qs_ks_qmmm_env_type,&
                                             qs_ks_qmmm_release
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE input_constants, ONLY: do_qmmm_grid, do_qmmm_gauss, do_qmmm_none
  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qs_ks_qmmm_methods'
  INTEGER, SAVE, PRIVATE :: last_ks_qmmm_nr=0

  PUBLIC :: ks_qmmm_env_rebuild, qs_ks_qmmm_create
CONTAINS
  !!****f* force_env_methods/ks_qmmm_env_rebuild [1.0] *
  !!
  !!   NAME
  !!     ks_qmmm_env_rebuild
  !!
  !!   FUNCTION
  !!     Initialize the ks_qmmm_env
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     Teodoro Laino
  !!
  !!   MODIFICATION HISTORY
  !!     05.2004 created [tlaino]
  !!
  !!*** **********************************************************************
  SUBROUTINE ks_qmmm_env_rebuild(qs_env,qmmm_env,error)
    ! Arguments
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error
    TYPE(qs_environment_type), OPTIONAL, &
         POINTER                                :: qs_env 
    TYPE(qs_ks_qmmm_env_type), POINTER       :: ks_qmmm_env
    TYPE(qmmm_env_qm_type), POINTER             :: qmmm_env
    ! Local Variables
    
    NULLIFY(ks_qmmm_env)
    CALL get_qs_env(qs_env=qs_env,&
         ks_qmmm_env=ks_qmmm_env)
    
    !   *** allocate the ks_qmmm env if not allocated yet!**
    IF (.NOT.ASSOCIATED(ks_qmmm_env)) THEN
       CALL qs_ks_qmmm_create(ks_qmmm_env=ks_qmmm_env,qs_env=qs_env,&
            qmmm_env=qmmm_env,error=error)
       CALL set_qs_env(qs_env=qs_env, ks_qmmm_env=ks_qmmm_env)
       CALL qs_ks_qmmm_release(ks_qmmm_env=ks_qmmm_env,error=error)
    END IF
  END SUBROUTINE ks_qmmm_env_rebuild
  !***************************************************************************

!!****f* qs_ks_methods/qs_ks_qmmm_create [1.0] *
!!
!!   NAME
!!     qs_ks_qmmm_create
!!
!!   FUNCTION
!!     allocates and initializes the given ks_qmmm_env.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - ks_qmmm_env: the ks_qmmm env to be initialized
!!     - qs_env: the qs environment
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     05.2004 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE qs_ks_qmmm_create(ks_qmmm_env, qs_env, qmmm_env, error)
    TYPE(qs_ks_qmmm_env_type), POINTER            :: ks_qmmm_env
    TYPE(qmmm_env_qm_type), pointer :: qmmm_env
    TYPE(qs_environment_type), POINTER            :: qs_env
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ks_qmmm_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat, cmax, igrid
    LOGICAL                                  :: failure
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(pw_pool_p_type), DIMENSION (:), pointer :: pools
    TYPE(cube_info_type), DIMENSION(:), pointer :: cube_info

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
    CPPreconditionNoFail(.NOT.ASSOCIATED(ks_qmmm_env),cp_failure_level,routineP,error)
    ALLOCATE(ks_qmmm_env, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    NULLIFY( ks_qmmm_env%pw_env,&
             ks_qmmm_env%cube_info)
    IF (.not.failure) THEN
       NULLIFY( auxbas_pw_pool)
       CALL get_qs_env(qs_env=qs_env,&
            pw_env=ks_qmmm_env%pw_env)
       CALL pw_env_get(ks_qmmm_env%pw_env,auxbas_pw_pool=auxbas_pw_pool)
       CALL pw_env_retain(ks_qmmm_env%pw_env,error=error)

       ks_qmmm_env%n_evals=0
       ks_qmmm_env%ref_count=1
       last_ks_qmmm_nr=last_ks_qmmm_nr+1
       ks_qmmm_env%id_nr=last_ks_qmmm_nr

       CALL pw_pool_init_coeff(auxbas_pw_pool,ks_qmmm_env%v_qmmm_rspace,&
            use_data=REALDATA3D, in_space=REALSPACE)
  
       IF (qmmm_env%qmmm_coupl_type.EQ.do_qmmm_gauss) THEN
          CALL pw_env_get(ks_qmmm_env%pw_env,pw_pools=pools,error=error)
          ALLOCATE(cube_info(SIZE(pools)), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO igrid=1,SIZE(pools)
             CALL init_cube_info(cube_info(igrid),&
               pools(igrid)%pool%pw_grid%dr(:),&
               qmmm_env%maxRadius(igrid),cmax)
          END DO
          ks_qmmm_env%cube_info => cube_info
       END IF
       !
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE qs_ks_qmmm_create
!***************************************************************************

END MODULE qs_ks_qmmm_methods
