!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****m* cp2k/qs_arpack_tddfpt *
!!
!!   NAME
!!     qs_arpack_tddfpt
!!
!!   FUNCTION
!!     Provides an environment to calculate the action of the
!!     tddfpt operator on a vector. This is the routine that
!!     arpack needs.
!!
!!   SYNOPSIS
!!     
!!
!!   MODIFICATION HISTORY
!!     09-2002, TCH, created
!!
!!   SEE ALSO
!!     qs_arpack_diag_n, qs_arpack_full
!!
!!****
MODULE qs_arpack_tddfpt

  USE global_types, ONLY : global_environment_type
  USE kinds,        ONLY : dbl
  USE qs_p_types,   ONLY : qs_p_env_type
  USE qs_p_utils,   ONLY : p_op_l1, p_op_l2

  IMPLICIT NONE

  TYPE op_tddfpt_type
     TYPE(qs_p_env_type), POINTER :: p_env
     CHARACTER                    :: mode
  END TYPE op_tddfpt_type

  CHARACTER(len=*), PARAMETER :: moduleN="qs_arpack_tddfpt"

  PRIVATE
  PUBLIC :: arpack_op_tddfpt, op_tddfpt_type, setup_op_tddfpt
  
CONTAINS

! *****************************************************************************

!!****f* cp2k/qs_arpack_tddfpt/arpack_op_tddfpt *
!!
!!   NAME
!!     arpack_op_tddfpt
!!
!!   FUNCTION
!!     Calculate the action of the tddfpt operator on a
!!     input vector vin. Store output in vout
!!
!!   ARGUMENTS
!!     - op       : operator containing needed info
!!     - vin      : vector to operate on
!!     - vout     : result
!!     - glob_env : global environment
!!
!!   NOTES
!!     Doesn't work yet.
!!
!!   MODIFICATION HISTORY
!!     09-2002, TCH, created
!!
!!****
  SUBROUTINE arpack_op_tddfpt(op, vin, vout, glob_env)
  
    IMPLICIT NONE

  	! arguments
    TYPE(op_tddfpt_type), INTENT(inout)       :: op
    REAL(dbl), DIMENSION(:), POINTER          :: vin, vout
    TYPE(global_environment_type), INTENT(in) :: glob_env
  
  	! locals
    CHARACTER(len=*), PARAMETER               :: routineN='arpack_op_tddfpt',&
                                                 routineP=moduleN//':'//routineN
  
    ! code
    IF (glob_env%ionode) THEN
       WRITE (glob_env%scr, '(A)') routineP//": Not yet implemented"
    END IF

  END SUBROUTINE arpack_op_tddfpt

! *****************************************************************************

!!****f* cp2k/qs_arpack_tddfpt/setup_op_tddfpt *
!!
!!   NAME
!!     setup_op_tddfpt
!!
!!   FUNCTION
!!     -
!!
!!   ARGUMENTS
!!     -
!!
!!   NOTES
!!     -
!!
!!   MODIFICATION HISTORY
!!     , TCH, created
!!
!!   SEE ALSO
!!     
!!
!!****
  SUBROUTINE setup_op_tddfpt(op, p_env, mode, glob_env)
  
    IMPLICIT NONE

  	! arguments
    TYPE(op_tddfpt_type), INTENT(inout) :: op
    TYPE(qs_p_env_type), TARGET, INTENT(inout) :: p_env
    CHARACTER, INTENT(in)               :: mode
    TYPE(global_environment_type), INTENT(in) :: glob_env

  	! locals
    CHARACTER(len=*), PARAMETER               :: routineN='setup_op_tddfpt',&
                                                 routineP=moduleN//':'//routineN
  
  	! code
    IF (glob_env%ionode) THEN
       WRITE (glob_env%scr, '(A)') routineP//": Not yet implemented"
    END IF
  
    op%p_env => p_env
    op%mode = mode

  END SUBROUTINE setup_op_tddfpt

! *****************************************************************************

END MODULE qs_arpack_tddfpt
