!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief A common interface (wrapper) for a callback into the md_run loop.
!>        Currently this is only used by the glbopt machinery, but its meant
!>        to be extended if others need to controll the md_run loop, too.
!>
!> \par History
!>      11.2012 created [Ole]
!> \author Ole
! *****************************************************************************
MODULE mdctrl_methods
  USE glbopt_callback,                 ONLY: glbopt_md_callback
  USE md_environment_types,            ONLY: md_environment_type
  USE mdctrl_types,                    ONLY: mdctrl_type
  USE message_passing,                 ONLY: mp_abort
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

 PUBLIC :: mdctrl_callback

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mdctrl_methods'

CONTAINS

! *****************************************************************************
!> \brief This is called by md_run for each step during during its main-loop.
!> \param mdctrl data which is passed on to the wrapped client-routine
!> \param md_env contains the current state of the md_run
!> \param should_stop can be used to abort the md_run
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE mdctrl_callback(mdctrl, md_env, should_stop, error)
    TYPE(mdctrl_type), POINTER               :: mdctrl
    TYPE(md_environment_type), POINTER       :: md_env
    LOGICAL, INTENT(inout)                   :: should_stop
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mdctrl_callback', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(md_env), cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(mdctrl), cp_failure_level,routineP,error,failure)

    IF(ASSOCIATED(mdctrl%glbopt)) THEN
      CALL glbopt_md_callback(mdctrl%glbopt, md_env, should_stop, error)

    !ELSE IF(ASSOCIATED(mdctrl%your_own_hook)) THEN ...

    ELSE
      CALL mp_abort("mdctrl_callback: No hook found.")
    ENDIF

  END SUBROUTINE mdctrl_callback

END MODULE mdctrl_methods

