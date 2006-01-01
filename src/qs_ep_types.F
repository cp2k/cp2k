!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


!!****h* cp2k/qs_ep_types [1.0] *
!!
!!   NAME
!!     qs_ep_types
!!
!!   FUNCTION
!!     the types needed for the calculation of the enrgy with perturbation
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     11.2002 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE qs_ep_types
  USE cp_fm_types,                     ONLY: cp_fm_p_type
  USE cp_fm_vect,                      ONLY: cp_fm_vect_dealloc
  USE kinds,                           ONLY: dp
  USE preconditioner,                  ONLY: destroy_preconditioner,&
                                             preconditioner_p_type
  USE qs_environment_types,            ONLY: qs_env_release,&
                                             qs_environment_p_type,&
                                             qs_environment_type
  USE qs_p_env_types,                  ONLY: p_env_release,&
                                             qs_p_env_type
  USE qs_p_sparse_psi,                 ONLY: p_proj_release,&
                                             qs_p_projection_p_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qs_ep_types'

  PUBLIC :: qs_ep_env_type, qs_ep_env_p_type
  PUBLIC :: ep_env_retain, ep_env_release, ep_env_get
!!***
!****************************************************************************

!!****s* qs_ep_types/qs_ep_env_type [1.0] *
!!
!!   NAME
!!     qs_ep_env_type
!!
!!   FUNCTION
!!     environement for an ep calculation
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - id_nr: identification number (unique)
!!     - ref_count: reference count (see doc/ReferenceCounting.html)
!!     - sub_qs_env: qs_env that represent the pieces of the full system
!!     - sub_proj: projection from the full system to the pieces
!!     - main_qs_env: qs_env of the full system
!!     - main_p_env: perturbation environment for the full system
!!     - m_pi_Hrho_psi0d: the rhs of the linear system to solve for
!!       the perturbation ( - pi Hrho psi0d )
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     11.2002 created [fawzi]
!!
!!   SOURCE
  !***************************************************************************
  TYPE qs_ep_env_type
     INTEGER :: id_nr, ref_count
     TYPE(qs_environment_p_type), DIMENSION(:), POINTER :: sub_qs_env
     TYPE(qs_p_projection_p_type), DIMENSION(:), POINTER :: sub_proj
     TYPE(qs_environment_type), POINTER :: main_qs_env
     TYPE(qs_p_env_type), POINTER :: main_p_env
     TYPE(cp_fm_p_type), DIMENSION(:), POINTER :: m_pi_Hrho_psi0d
     TYPE(cp_fm_p_type), DIMENSION(:), POINTER :: psi1
     TYPE(preconditioner_p_type), DIMENSION(:), POINTER :: precond
  END TYPE qs_ep_env_type
!!***
  !****************************************************************************

!!****s* qs_ep_types/qs_ep_env_p_type [1.0] *
!!
!!   NAME
!!     qs_ep_env_p_type
!!
!!   FUNCTION
!!     type to have array of pointers
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - ep_env: the pointer to the ep env
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     12.2002 created [fawzi]
!!
!!   SOURCE
  !***************************************************************************
  TYPE qs_ep_env_p_type
     TYPE(qs_ep_env_type), POINTER :: ep_env
  END TYPE qs_ep_env_p_type
!!***
  !****************************************************************************

CONTAINS

!!****f* qs_ep_types/ep_env_retain [1.0] *
!!
!!   NAME
!!     ep_env_retain
!!
!!   SYNOPSIS
!!     Subroutine ep_env_retain(ep_env, error)
!!       Type(qs_ep_env_type), Pointer:: ep_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine ep_env_retain
!!
!!   FUNCTION
!!     retains an ep environment (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - ep_env: the environment to retain
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
SUBROUTINE ep_env_retain(ep_env, error)
    TYPE(qs_ep_env_type), POINTER            :: ep_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(ep_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(ep_env%ref_count>0,cp_failure_level,routineP,error,failure)
     ep_env%ref_count=ep_env%ref_count+1
  END IF
END SUBROUTINE ep_env_retain
!***************************************************************************

!!****f* qs_ep_types/ep_env_release [1.0] *
!!
!!   NAME
!!     ep_env_release
!!
!!   SYNOPSIS
!!     Subroutine ep_env_release(ep_env, error)
!!       Type(qs_ep_env_type), Pointer:: ep_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine ep_env_release
!!
!!   FUNCTION
!!     releases the given ep_env (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - ep_env: the ep_env to release
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
SUBROUTINE ep_env_release(ep_env, error)
    TYPE(qs_ep_env_type), POINTER            :: ep_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (ASSOCIATED(ep_env)) THEN
     CPPrecondition(ep_env%ref_count>0,cp_failure_level,routineP,error,failure)
     ep_env%ref_count=ep_env%ref_count-1
     IF (ep_env%ref_count<1) THEN
        IF (ASSOCIATED(ep_env%sub_qs_env)) THEN
           DO i=1,SIZE(ep_env%sub_qs_env)
              CALL qs_env_release(ep_env%sub_qs_env(i)%qs_env,error=error)
           END DO
           DEALLOCATE(ep_env%sub_qs_env, stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(ep_env%sub_proj)) THEN
           DO i=1,SIZE(ep_env%sub_proj)
              CALL p_proj_release(ep_env%sub_proj(i)%projection,error=error)
           END DO
           DEALLOCATE(ep_env%sub_proj, stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        CALL qs_env_release(ep_env%main_qs_env, error=error)
        CALL p_env_release(ep_env%main_p_env,error=error)
        CALL cp_fm_vect_dealloc(ep_env%m_pi_Hrho_psi0d,error=error)
        CALL cp_fm_vect_dealloc(ep_env%psi1,error=error)
        IF (ASSOCIATED(ep_env%precond)) THEN
           DO i=1,1 !SIZE(ep_env%precond) ! same precond for all spins
              CALL destroy_preconditioner(ep_env%precond(i)%preconditioner)
           END DO
           DEALLOCATE(ep_env%precond,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        DEALLOCATE(ep_env, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
END SUBROUTINE ep_env_release
!***************************************************************************

!!****f* qs_ep_methods/ep_env_get [1.0] *
!!
!!   NAME
!!     ep_env_get
!!
!!   SYNOPSIS
!!     Subroutine ep_env_get(ep_env, id_nr, sub_qs_env, sub_proj,&
!!         main_qs_env, main_p_env, preconditioner, error)
!!       Type(qs_ep_env_type), Pointer:: ep_env
!!       Integer, Intent (OUT), Optional:: id_nr
!!       Type(qs_environment_p_type), Dimension(:), Pointer, Optional::&
!!         sub_qs_env
!!       Type(qs_p_projection_p_type), Dimension(:), Pointer, Optional::&
!!         sub_proj
!!       Type(qs_environment_type), Pointer, Optional:: main_qs_env
!!       Type(qs_p_env_type), Pointer, Optional:: main_p_env
!!       Type(preconditioner_p_type), Dimension(:), Pointer, Optional::&
!!         preconditioner
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine ep_env_get
!!
!!   FUNCTION
!!     returns the various attributes of the ep_env
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - ep_env: the ep_env to get the info from
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!     For the other arguments see the attributes of qs_ep_env_type
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     12.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE ep_env_get(ep_env,id_nr, sub_qs_env, sub_proj, main_qs_env,&
     main_p_env, preconditioner, error)
    TYPE(qs_ep_env_type), POINTER            :: ep_env
    INTEGER, INTENT(out), OPTIONAL           :: id_nr
    TYPE(qs_environment_p_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: sub_qs_env
    TYPE(qs_p_projection_p_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: sub_proj
    TYPE(qs_environment_type), OPTIONAL, &
      POINTER                                :: main_qs_env
    TYPE(qs_p_env_type), OPTIONAL, POINTER   :: main_p_env
    TYPE(preconditioner_p_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: preconditioner
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_env_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(ep_env),cp_failure_level,routineP,error,failure)
  IF (.not.failure) THEN
     CPPrecondition(ep_env%ref_count>0,cp_failure_level,routineP,error,failure)
  END IF

  IF (.NOT. failure) THEN
     IF (PRESENT(id_nr)) id_nr=ep_env%id_nr
     IF (PRESENT(sub_qs_env)) sub_qs_env => ep_env%sub_qs_env
     IF (PRESENT(sub_proj)) sub_proj => ep_env%sub_proj
     IF (PRESENT(main_qs_env)) main_qs_env => ep_env%main_qs_env
     IF (PRESENT(main_p_env)) main_p_env => ep_env%main_p_env
!FM     IF (PRESENT(nspins)) THEN
!FM        CPPrecondition(ASSOCIATED(ep_env%main_qs_env),cp_failure_level,routineP,error,failure)
!FM        IF (.not.failure) THEN
!FM           nspins=ep_env%main_qs_env%dft_control%nspins
!FM        END IF
!FM     END IF
     IF (PRESENT(preconditioner)) THEN
        preconditioner => ep_env%precond
     END IF
  END IF
END SUBROUTINE ep_env_get
!***************************************************************************

END MODULE qs_ep_types
