!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****m* cp2k/m *
!!
!!   NAME
!!     m
!!
!!   FUNCTION
!!     A module to cure the fucked up naming conventions in the whole
!!     blacs / full_matrix mess.
!!
!!
!!   MODIFICATION HISTORY
!!     10-2002, TCH, created
!!
!!****
module cp_fm_utils

  use cp_error_handling, only : cp_simple_assert
  use cp_full_matrix
  use kinds,             only : dbl, dp_size
  use message_passing,   only : mp_sum
  use qs_blacs,          only : blacs_replicated_copy, blacs_get_element
  use termination,       only : stop_memory

  implicit none

  character(len=*), parameter :: moduleN='cp_fm_utils'
  logical, parameter :: DEBUG_THIS_MODULE=.FALSE.

  private
  public :: copy_fm_col_to_rep_vec, &
            copy_rep_vec_to_fm_col, &
            copy_fm_to_rep_vec, &
            copy_rep_vec_to_fm, &
            copy_fm_col_to_fm

contains

! *****************************************************************************

  subroutine copy_fm_col_to_rep_vec(mat_a, col_a, vec_b)
  
    implicit none

  	! arguments
    type(cp_full_matrix_type), pointer :: mat_a
    integer, intent(in)                :: col_a
    real(dbl), dimension(:), pointer   :: vec_b
  
  	! locals
    integer :: stat
    integer :: a_row_global, a_col_global, &
               a_row_block, a_col_block, &
               a_row_local, a_col_local
    integer :: row, col
    real(dbl) :: element
    integer, dimension(:), pointer :: a_row_indices, a_col_indices
    character(len=*), parameter :: routineN='copy_fm_col_to_rep_vec', &
                                   routineP=moduleN//'/'//routineN

    ! code
    call cp_fm_get_info(matrix=mat_a, &
                        nrow_global=a_row_global, ncol_global=a_col_global, &
                        nrow_block=a_row_block,   ncol_block=a_col_block, &
                        nrow_local=a_row_local,   ncol_local=a_col_local, &
                        row_indices=a_row_indices, col_indices=a_col_indices)

    CPAAssert(a_row_global == size(vec_b))
    CPAAssert(col_a >=1 .and. a_col_global>=col_a)

    vec_b(:) = 0.0_dbl

    do col=1, a_col_local
       if (a_col_indices(col)==col_a) then
          do row=1, a_row_local
             vec_b(a_row_indices(row)) = mat_a%local_data(row, col)
          end do
       end if
    end do

#if defined(__SCALAPACK)
    call mp_sum(vec_b, mat_a%matrix_struct%para_env%group)
#endif

    if (DEBUG_THIS_MODULE) then
       do row=1, a_row_global
          call blacs_get_element(mat_a, row, col_a, element)
#if defined(__SCALAPACK)
          call mp_sum(element, mat_a%matrix_struct%para_env%group)
#endif
          if (element/=vec_b(row)) then
             write (*,*) routineP//': inconsistent data'
             stop 1
          end if
       end do
    end if

  end subroutine copy_fm_col_to_rep_vec

! *****************************************************************************

  subroutine copy_rep_vec_to_fm_col(vec_a, mat_b, col_b)

    implicit none

    ! arguments
    real(dbl), dimension(:), pointer   :: vec_a
    type(cp_full_matrix_type), pointer :: mat_b
    integer, intent(in)                :: col_b

    ! locals
    integer :: b_row_global, b_col_global, &
               b_row_block, b_col_block, &
               b_row_local, b_col_local
    integer :: row, col
    integer, dimension(:), pointer :: b_row_indices, b_col_indices
    real(DBL) :: element
    character(len=*), parameter :: routineN='copy_rep_vec_to_fm', &
                                   routineP=moduleN//'/'//routineN

    ! code
    call cp_fm_get_info(matrix=mat_b, &
                        nrow_global=b_row_global, ncol_global=b_col_global, &
                        nrow_block=b_row_block,   ncol_block=b_col_block, &
                        nrow_local=b_row_local,   ncol_local=b_col_local, &
                        row_indices=b_row_indices, col_indices=b_col_indices)

    CPAAssert(size(vec_a) == b_row_global)
    CPAAssert(col_b >= 1 .and. col_b <= b_col_global)

    ! very simple code not very efficient
    do col=1, b_col_local
          if (b_col_indices(col) == col_b) then
             do row=1, b_row_local
                mat_b%local_data(row,col) = vec_a(b_row_indices(row))
             end do
          end if
    end do

    if (DEBUG_THIS_MODULE) then
          do row=1, b_row_global 
             call blacs_get_element(mat_b, row, col_b, element)
#if defined(__SCALAPACK)
             call mp_sum(element, mat_b%matrix_struct%para_env%group)
#endif
             if (element/=vec_a(row)) then
                write (*,*) routineP//': Inconsistent data'
                stop 1
             end if
          end do
    end if

  end subroutine copy_rep_vec_to_fm_col

! *****************************************************************************

  subroutine copy_fm_to_rep_vec(mat_a, vec_b)

    implicit none

    ! arguments
    type(cp_full_matrix_type), pointer :: mat_a
    real(dbl), dimension(:), pointer   :: vec_b

    ! locals
    integer :: a_row_global, a_col_global, &
               a_row_block, a_col_block, &
               a_row_local, a_col_local
    integer :: row, col, vec_b_pos, i, j
    integer, dimension(:), pointer :: a_row_indices, a_col_indices
    real(DBL) :: element
    character(len=*), parameter :: routineN='copy_fm_to_rep_vec', &
                                   routineP=moduleN//'/'//routineN

    ! code
    call cp_fm_get_info(matrix=mat_a, &
                        nrow_global=a_row_global, ncol_global=a_col_global, &
                        nrow_block=a_row_block,   ncol_block=a_col_block, &
                        nrow_local=a_row_local,   ncol_local=a_col_local, &
                        row_indices=a_row_indices, col_indices=a_col_indices)

    CPAAssert(a_row_global*a_col_global==size(vec_b)) ! the vector must fit
    
       do row=1, a_row_local
    do col=1, a_col_local
          i = a_row_indices(row)
          j = a_col_indices(col)
          vec_b_pos = (a_row_global*(j-1))+i
          CPAAssert(vec_b_pos>=1.and.vec_b_pos<=size(vec_b))
          vec_b(vec_b_pos) = mat_a%local_data(row,col)
       end do
    end do

#if defined(__SCALAPACK)
    call mp_sum(vec_b, mat_a%matrix_struct%para_env%group)
#endif
    
    if (DEBUG_THIS_MODULE) then
       do col=1, a_col_global
          do row=1, a_row_global
             call blacs_get_element(mat_a, row, col, element)
#if defined(__SCALAPACK)
             call mp_sum(element, mat_a%matrix_struct%para_env%group)
#endif
             if (element/=vec_b(((col-1)*a_row_global)+row)) then
                write (*,*) routineP//': Inconsistent data'
                stop 1
             end if
          end do
       end do
    end if

  end subroutine copy_fm_to_rep_vec

! *****************************************************************************

  subroutine copy_rep_vec_to_fm(vec_a, mat_b)

    implicit none

    ! arguments
    real(dbl), dimension(:), pointer   :: vec_a
    type(cp_full_matrix_type), pointer :: mat_b

    ! locals
    integer :: b_row_global, b_col_global, &
               b_row_block, b_col_block, &
               b_row_local, b_col_local
    integer :: row, col, vec_a_pos, i, j
    integer, dimension(:), pointer :: b_row_indices, b_col_indices
    real(DBL) :: element
    character(len=*), parameter :: routineN='copy_rep_vec_to_fm', &
                                   routineP=moduleN//'/'//routineN

    ! code
    call cp_fm_get_info(matrix=mat_b, &
                        nrow_global=b_row_global, ncol_global=b_col_global, &
                        nrow_block=b_row_block,   ncol_block=b_col_block, &
                        nrow_local=b_row_local,   ncol_local=b_col_local, &
                        row_indices=b_row_indices, col_indices=b_col_indices)

    CPAAssert(b_row_global*b_col_global==size(vec_a)) ! the matrix must fit

    ! very simple code not very efficient
       do row=1, b_row_local
    do col=1, b_col_local
          i = b_row_indices(row)
          j = b_col_indices(col)
          vec_a_pos = (b_row_global*(j-1)) + i
          CPAAssert(vec_a_pos>=1.and.vec_a_pos<=size(vec_a))
          mat_b%local_data(row,col) = vec_a(vec_a_pos)
       end do
    end do

    if (DEBUG_THIS_MODULE) then
       do col=1, b_col_global
          do row=1, b_row_global 
             call blacs_get_element(mat_b, row, col, element)
#if defined(__SCALAPACK)
             call mp_sum(element, mat_b%matrix_struct%para_env%group)
#endif
             if (element/=vec_a(((col-1)*b_row_global)+row)) then
                write (*,*) routineP//': Inconsistent data'
                stop 1
             end if
          end do
       end do
    end if

  end subroutine copy_rep_vec_to_fm

! *****************************************************************************

  subroutine copy_fm_col_to_fm(mat_a, col_a, mat_b)
    
    implicit none

    ! arguments
    type(cp_full_matrix_type), pointer :: mat_a
    integer, intent(in)                :: col_a
    type(cp_full_matrix_type), pointer :: mat_b

    ! locals
    integer :: a_row_global, stat
    real(dbl), dimension(:), pointer :: temp_vec
    character(len=*), parameter :: routineN='copy_fm_col_to_fm', &
                                   routineP=moduleN//'/'//routineN
    

    call cp_fm_get_info(matrix=mat_a, nrow_global=a_row_global)

    allocate(temp_vec(a_row_global), stat=stat)
    if (stat/=0) call stop_memory(routineP, 'temp_vec', a_row_global*dp_size)

    call copy_fm_col_to_rep_vec(mat_a, col_a, temp_vec)
    call copy_rep_vec_to_fm(temp_vec, mat_b)

    deallocate(temp_vec, stat=stat)
    if (stat/=0) call stop_memory(routineP, 'temp_vec')

  end subroutine copy_fm_col_to_fm

! *****************************************************************************

end module cp_fm_utils

