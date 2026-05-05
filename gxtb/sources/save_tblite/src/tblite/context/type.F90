! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later
!
! tblite is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! tblite is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with tblite.  If not, see <https://www.gnu.org/licenses/>.

!> @file tblite/context/type.f90
!> Provides a context manager for storing persistent environment settings

#ifndef WITH_MKL
#define WITH_MKL 0
#endif

!> Calculation context for storing and communicating with the environment
module tblite_context_type
   use iso_fortran_env, only : output_unit
   use mctc_env, only : wp, error_type
   use tblite_context_logger, only : context_logger
   use tblite_context_solver, only : context_solver
   use tblite_context_terminal, only : context_terminal
   use tblite_scf_solver, only : solver_type
   implicit none
   private


   !> Calculation context type for error handling and output messages
   type, public :: context_type
      !> Default output unit for this context
      integer :: unit = output_unit
      !> Default verbosity for procedures using this context
      integer :: verbosity = 1
      !> Stack containing the error messages of this context
      type(error_type), allocatable :: error_log(:)
      !> Optional logger to be used for writing messages
      class(context_logger), allocatable :: io
      !> Optional factory for creating electronic solvers
      class(context_solver), allocatable :: solver
      !> Color support for output
      type(context_terminal) :: terminal = context_terminal()
   contains
      !> Write a message to the output
      procedure :: message
      !> Push an error message to the context
      procedure :: set_error
      !> Pop an error message from the context
      procedure :: get_error
      !> Query the context for errors
      procedure :: failed
      !> Create electronic solver instance
      procedure :: new_solver
      !> Delete an electronic solver instance
      procedure :: delete_solver
   end type context_type


contains


!> Add an error message to the context
subroutine set_error(self, error)
   !> Instance of the calculation context
   class(context_type), intent(inout) :: self
   !> Error handling
   type(error_type), intent(in), optional :: error

   if (present(error)) then
      if (.not.allocated(self%error_log)) allocate(self%error_log(0))
      self%error_log = [self%error_log, error]
   end if
end subroutine set_error


!> Pop an error message from the context
subroutine get_error(self, error)
   !> Instance of the calculation context
   class(context_type), intent(inout) :: self
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   if (.not.allocated(self%error_log)) allocate(self%error_log(0))

   if (size(self%error_log) > 0) then
      error = self%error_log(size(self%error_log))
      self%error_log = self%error_log(:size(self%error_log)-1)
   end if
end subroutine get_error


!> Write a message to the output
subroutine message(self, msg)
   !> Instance of the calculation context
   class(context_type), intent(inout) :: self
   !> Message to write
   character(len=*), intent(in) :: msg
   type(error_type), allocatable :: error

   if (allocated(self%io)) then
      call self%io%message(msg, error)
      if (allocated(error)) then
         call self%set_error(error)
      end if
   else
      write(self%unit, '(a)') msg
   end if
end subroutine message


!> Query the context for errors
pure function failed(self)
   !> Instance of the calculation context
   class(context_type), intent(in) :: self
   !> Error status of the context
   logical :: failed

   failed = .false.
   if (allocated(self%error_log)) then
      failed = size(self%error_log) > 0
   end if
end function failed


!> Create new electronic solver
subroutine new_solver(self, solver, overlap, nel, kt)
   use tblite_lapack_solver, only : lapack_solver
   !> Instance of the calculation context
   class(context_type), intent(inout) :: self
   !> New electronic solver
   class(solver_type), allocatable, intent(out) :: solver
   !> Overlap matrix
   real(wp), intent(in) :: overlap(:, :)
   !> Number of electrons per spin channel
   real(wp), intent(in) :: nel(:)
   !> Electronic temperature
   real(wp), intent(in) :: kt

   if (.not.allocated(self%solver)) then
      self%solver = lapack_solver()
   end if

   call self%solver%new(solver, overlap, nel, kt)
end subroutine new_solver


!> Delete electronic solver instance
subroutine delete_solver(self, solver)
   !> Instance of the calculation context
   class(context_type), intent(inout) :: self
   !> Electronic solver instance
   class(solver_type), allocatable, intent(inout) :: solver

   if (allocated(self%solver)) then
      call self%solver%delete(solver)
   end if
#if WITH_MKL 
   call mkl_free_buffers()
#endif
   if (allocated(solver)) deallocate(solver)
end subroutine delete_solver


end module tblite_context_type
