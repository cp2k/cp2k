!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****m* cp2k/qs_arpack_full *
!!
!!   NAME
!!     qs_arpack_full
!!
!!   FUNCTION
!!     
!!   NOTES
!!
!!   SYNOPSIS
!!     
!!
!!   MODIFICATION HISTORY
!!
!!
!!****

! seems to be efficient for neig << n as the work is +- neig^2*n
! the prefactor is such that, it seems not to be highly efficient for 
! finding all the occupied orbitals if a modest (normal) basis set is used.
! with small modifications should also be able to compute eigenvectors/eigenvalues
! of different parts of the spectrum / without cholesky decomposition
!
! access of the amatrix is only through multiplications of a vector with amatrix

module qs_arpack_full

  use kinds, only: wp => dp, wp_size => dp_size
  use termination,         only: stop_memory, stop_program
  use qs_blacs,            only: cp_full_matrix_type, &
                                 blacs_vector_symv, &
                                 cp_fm_get_info

  implicit none
  
  type full_m_op_type
     type(cp_full_matrix_type), pointer :: a
     integer                            :: a_size
     logical                            :: general
     type(cp_full_matrix_type), pointer :: b
     type(cp_full_matrix_type), pointer :: buffer
  end type full_m_op_type

  character(len=*), parameter :: moduleN = "qs_arpack_full"

  private
  public :: full_m_op_type, arpack_full_m_op, setup_full_m_op

contains

! *****************************************************************************

!
! notice that bmatrix can be used to solve the general eigenvalue problem.
! bmatrix has to be the cholesky decomp.
! L^X A L^X  should be calculate as described in the ARPACK manual
! right now, it is not implemented
!
subroutine arpack_full_m_op(full_m_op, vin, vout)
  implicit none
  type(full_m_op_type), intent(in)          :: full_m_op
  real(wp), dimension(:), pointer           :: vin, vout

  call blacs_vector_symv(full_m_op%a, vin, vout, full_m_op%buffer)
  
end subroutine arpack_full_m_op

! *****************************************************************************

subroutine setup_full_m_op(op, a, buffer, general, b)
  implicit none
  ! arguments
  type(full_m_op_type), intent(out)            :: op
  type(cp_full_matrix_type), pointer           :: a, buffer
  logical, intent(in)                          :: general
  type(cp_full_matrix_type), pointer, optional :: b

  ! locals
  character(len=*), parameter :: routineN = "setup_full_m_op", &
                                 routineP = moduleN//"/"//routineN

  ! code
  op%a => a
  call cp_fm_get_info(op%a, nrow_global=op%a_size)
  op%buffer => buffer
  op%general = general
  if (op%general .and. .not.present(b)) then
     call stop_program(routineP, "The general problem requires a b matrix")
  end if
  if (present(b)) op%b => b
end subroutine setup_full_m_op

! *****************************************************************************

end module qs_arpack_full

