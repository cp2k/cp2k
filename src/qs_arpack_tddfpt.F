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

#include "cp_prep_globals.h"

MODULE qs_arpack_tddfpt
  USE cp_error_handling,               ONLY: cp_simple_assert
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_full_matrix_type
  USE cp_fm_utils,                     ONLY: copy_fm_to_rep_vec,&
                                             copy_rep_vec_to_fm
  USE cp_log_handling,                 ONLY: cp_fatal_level,&
                                             cp_to_string
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl
  USE qs_blacs,                        ONLY: cp_fm_release,&
                                             replicate_blacs_matrix
  USE qs_environment_types,            ONLY: qs_environment_type
  USE qs_p_env_methods,                ONLY: p_op_l1,&
                                             p_op_l1_spin,&
                                             p_op_l2
  USE qs_p_env_types,                  ONLY: qs_p_env_type
  USE sparse_matrix_types,             ONLY: copy_sparse_to_full_matrix,&
                                             real_matrix_type
  IMPLICIT NONE
  
  TYPE op_tddfpt_type
     TYPE(qs_p_env_type), POINTER       :: p_env
     TYPE(qs_environment_type), POINTER :: qs_env
     CHARACTER                          :: mode  ! operator for x or y
     INTEGER                            :: spin  ! which spin to handle
  END TYPE op_tddfpt_type

  ! some constants
  CHARACTER(len=*), PARAMETER :: moduleN="qs_arpack_tddfpt"
  LOGICAL, PARAMETER          :: DEBUG_THIS_MODULE = .TRUE.

  PRIVATE
  PUBLIC :: arpack_op_tddfpt, op_tddfpt_type, setup_op_tddfpt
  
CONTAINS

! *****************************************************************************

!!****f* cp2k/qs_arpack_tddfpt/arpack_op_tddfpt *
!!
!!   NAME
!!     arpack_op_tddfpt
!!
!!   SYNOPSIS
!!     Subroutine arpack_op_tddfpt(op, vin, vout, glob_env)
!!       Implicit None
!!       Type(op_tddfpt_type), Intent (INOUT):: op
!!       Real(Kind=dbl), Dimension(:), Pointer:: vin, vout
!!       Type(global_environment_type), Intent (IN):: glob_env
!!     End Subroutine arpack_op_tddfpt
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
    LOGICAL                                   :: condition
    INTEGER                                   :: stat, row_global, col_global
    TYPE(cp_full_matrix_type), POINTER        :: v, res
    REAL(dbl), DIMENSION(:), POINTER          :: dummy_vec
    CHARACTER(len=*), PARAMETER               :: routineN='arpack_op_tddfpt',&
                                                 routineP=moduleN//'/'//routineN

    INTEGER :: m, n, lda

    ! code
    CPAAssert((0<op%spin).AND.(op%spin<=op%qs_env%dft_control%nspins))
    CPAAssert(SIZE(vin)==op%p_env%n_mo(op%spin)*op%p_env%n_ao(op%spin))

    ! allocate v and res
    CALL replicate_blacs_matrix(op%p_env%psi0d(op%spin)%matrix, v, &
         routineP//": v matrix")
    CALL replicate_blacs_matrix(op%p_env%psi0d(op%spin)%matrix, res, &
         routineP//": res matrix")

    CALL copy_rep_vec_to_fm(vin, v)
    
    CALL p_op_l1_spin(op%p_env, op%qs_env, op%spin, v, res)
!!!!!!!!! l2 part missing

    CALL copy_fm_to_rep_vec(res, vout)

    ! free v and res
    CALL cp_fm_release(v); NULLIFY(v)
    CALL cp_fm_release(res); NULLIFY(res)

  END SUBROUTINE arpack_op_tddfpt

! *****************************************************************************

!!****f* cp2k/qs_arpack_tddfpt/setup_op_tddfpt *
!!
!!   NAME
!!     setup_op_tddfpt
!!
!!   SYNOPSIS
!!     Subroutine setup_op_tddfpt(op, p_env, qs_env, mode, spin, glob_env)
!!       Implicit None
!!       Type(op_tddfpt_type), Intent (INOUT):: op
!!       Type(qs_p_env_type), Pointer:: p_env
!!       Type(qs_environment_type), Intent (INOUT), Target:: qs_env
!!       Character, Intent (IN):: mode
!!       Integer, Intent (IN):: spin
!!       Type(global_environment_type), Intent (IN):: glob_env
!!     End Subroutine setup_op_tddfpt
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
  SUBROUTINE setup_op_tddfpt(op, p_env, qs_env, mode, spin, glob_env)
  
    IMPLICIT NONE

  ! arguments
    TYPE(op_tddfpt_type), INTENT(inout)        :: op
    TYPE(qs_p_env_type), POINTER               :: p_env
    TYPE(qs_environment_type), INTENT(inout), TARGET :: qs_env
    CHARACTER, INTENT(in)                      :: mode
    INTEGER, INTENT(IN)                        :: spin
    TYPE(global_environment_type), INTENT(in)  :: glob_env

  ! locals
    CHARACTER(len=*), PARAMETER               :: routineN='setup_op_tddfpt',&
                                                 routineP=moduleN//':'//routineN
  
  ! code
    op%p_env => p_env
    op%qs_env => qs_env
    op%mode = mode
    op%spin = spin

  END SUBROUTINE setup_op_tddfpt

! *****************************************************************************

END MODULE qs_arpack_tddfpt
