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
module qs_arpack_tddfpt

  use kinds,      only : dbl
  use qs_p_types, only : qs_p_env_type
  use qs_p_utils, only : p_op_l1, p_op_l2

  implicit none

  private
  public :: arpack_op_tddfpt
  
contains

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
!!     - p_env : environment that contains necessary quantities
!!     - vin   : vector to operate on
!!     - vout  : result
!!
!!   NOTES
!!     Doesn't work yet.
!!
!!   MODIFICATION HISTORY
!!     09-2002, TCH, created
!!
!!****
  subroutine arpack_op_tddfpt(p_env, vin, vout)
  
    implicit none

  	! arguments
    type(qs_p_env_type), intent(inout) :: p_env
    real(dbl), dimension(:), pointer    :: vin, vout
  
  	! locals
  
  
  	! code
    ! calls to p_op_l1, p_op_l2 and cholesky factorization of 
    ! S matrix and so on
  
  end subroutine arpack_op_tddfpt

! *****************************************************************************

end module qs_arpack_tddfpt
