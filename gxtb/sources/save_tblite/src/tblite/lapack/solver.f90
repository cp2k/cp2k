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

!> @file tblite_lapack/solver.f90
!> Provides a wrapper for the eigenvalue solvers provided by LAPACK

!> LAPACK based eigenvalue solvers
module tblite_lapack_solver
   use mctc_env, only : wp
   use tblite_context_solver, only : context_solver
   use tblite_lapack_sygvd, only : sygvd_solver, new_sygvd
   use tblite_lapack_sygvr, only : sygvr_solver, new_sygvr
   use tblite_scf_solver, only : solver_type
   implicit none
   private

   public :: solver_type, lapack_algorithm


   !> Possible solvers provided by LAPACK
   type :: enum_lapack
      !> Divide-and-conquer solver
      integer :: gvd = 1
      !> Relatively robust solver
      integer :: gvr = 2
   end type enum_lapack

   !> Actual enumerator of possible solvers
   type(enum_lapack), parameter :: lapack_algorithm = enum_lapack()


   !> Generator for LAPACK based electronic solvers
   type, public, extends(context_solver) :: lapack_solver
      !> Selected electronic solver algorithm
      integer :: algorithm = lapack_algorithm%gvd
   contains
      !> Create new instance of electronic solver
      procedure :: new
      !> Delete an electronic solver instance
      procedure :: delete
   end type lapack_solver


contains


!> Create new electronic solver
subroutine new(self, solver, overlap, nel, kt)
   !> Instance of the solver factory
   class(lapack_solver), intent(inout) :: self
   !> New electronic solver
   class(solver_type), allocatable, intent(out) :: solver
   !> Overlap matrix
   real(wp), intent(in) :: overlap(:, :)
   !> Number of electrons per spin channel
   real(wp), intent(in) :: nel(:)
   !> Electronic temperature
   real(wp), intent(in) :: kt

   select case(self%algorithm)
   case(lapack_algorithm%gvd)
      block
         type(sygvd_solver), allocatable :: tmp
         allocate(tmp)
         call new_sygvd(tmp, overlap, nel, kt)
         call move_alloc(tmp, solver)
      end block
   case(lapack_algorithm%gvr)
      block
         type(sygvr_solver), allocatable :: tmp
         allocate(tmp)
         call new_sygvr(tmp, overlap, nel, kt)
         call move_alloc(tmp, solver)
      end block
   end select
end subroutine new


!> Delete electronic solver instance
subroutine delete(self, solver)
   !> Instance of the solver factory
   class(lapack_solver), intent(inout) :: self
   !> Electronic solver instance
   class(solver_type), allocatable, intent(inout) :: solver

   if (allocated(solver)) then
      call solver%delete()
      deallocate(solver)
   end if
end subroutine delete


end module tblite_lapack_solver
