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

  use global_types, only : global_environment_type
  use kinds, only: wp => dp, wp_size => dp_size
  use termination,         only: stop_memory, stop_program
  USE cp_fm_types,         ONLY: cp_fm_type, &
                                 cp_fm_get_info
  use qs_blacs,            only: blacs_vector_symv
  USE cp_para_types, ONLY: cp_para_env_type, cp_blacs_env_type
  USE cp_error_handling, ONLY: cp_error_type

  implicit none
  
  type op_symv_type
     type(cp_fm_type), pointer :: a
     integer                            :: a_size
     logical                            :: general
     type(cp_fm_type), pointer :: b
     type(cp_fm_type), pointer :: buffer
  end type op_symv_type

  character(len=*), parameter :: moduleN = "qs_arpack_full"

  private
  public :: arpack_op_symv, op_symv_type, setup_op_symv

contains

! *****************************************************************************

SUBROUTINE arpack_op_symv(op, vin, vout, blacs_env, para_env, error)
  implicit none
  ! arguments
  type(op_symv_type), intent(in)            :: op
  real(wp), dimension(:), pointer           :: vin, vout
  TYPE(cp_blacs_env_type), pointer          :: blacs_env
  TYPE(cp_para_env_type), pointer           :: para_env
  TYPE(cp_error_type), INTENT(inout), optional :: error
  ! code
  call blacs_vector_symv(op%a, vin, vout, op%buffer)
end subroutine arpack_op_symv

! *****************************************************************************

!!****f* cp2k/qs_arpack_full/setup_op_symv *
!!
!!   NAME
!!     setup_op_symv
!!
!!   FUNCTION
!!     Sets up the op_symv_type object that stores information on
!!     how to compute the arpack_op_symv operation.
!!
!!   ARGUMENTS
!!     - op     : the op_symv_type that will store the information
!!     - a      : the matrix A
!!     - buffer : a buffer for the calculations
!!     - general: wether the problem is standard or general
!!     - b      : the matrix B (OPTIONAL)
!!
!!   NOTES
!!     - b must be symmetric positive definite.
!!
!!   MODIFICATION HISTORY
!!     09-2002, TCH, created
!!
!!****
subroutine setup_op_symv(op, a, buffer, general, b)
  implicit none
  ! arguments
  type(op_symv_type), intent(out)              :: op
  type(cp_fm_type), pointer           :: a, buffer
  logical, intent(in)                          :: general
  type(cp_fm_type), pointer, optional :: b

  ! locals
  integer                     :: buffer_cols, buffer_rows
  character(len=*), parameter :: routineN = "setup_op_symv", &
                                 routineP = moduleN//"/"//routineN

  ! code
  op%a => a
  call cp_fm_get_info(op%a, nrow_global=op%a_size)

  op%buffer => buffer
  call cp_fm_get_info(op%buffer, ncol_global=buffer_cols, nrow_global=buffer_rows)
  if (buffer_cols < 2 .or. buffer_rows < op%a_size) then
     call stop_program(routineP, "buffer too small")
  end if
  
  op%general = general
  if (op%general .and. .not.present(b)) then
     call stop_program(routineP, "The general problem requires a b matrix")
  end if

  if (present(b)) op%b => b

end subroutine setup_op_symv

! *****************************************************************************

end module qs_arpack_full

