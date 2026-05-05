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

!> @file tblite/scf/solver.f90
!> Provides a base class for defining electronic solvers

!> Declaration of the abstract base class for electronic solvers
module tblite_scf_solver
   use mctc_env, only : sp, dp, wp, error_type
   use tblite_blas, only : gemm
   implicit none
   private

   !> Abstract base class for electronic solvers
   type, public, abstract :: solver_type
      !> Electronic temperature
      real(wp) :: kt
      !> Number of electrons per spin channel
      real(wp), allocatable :: nel(:)
   contains
      procedure(get_density), deferred :: get_density
      procedure(get_density), deferred :: get_wdensity
      procedure(delete), deferred :: delete
   end type solver_type

   abstract interface
      subroutine get_density(self, hmat, smat, eval, focc, density, error)
         import :: wp, solver_type, error_type
         !> Solver for the general eigenvalue problem
         class(solver_type), intent(inout) :: self
         !> Overlap matrix
         real(wp), contiguous, intent(in) :: smat(:, :)
         !> Hamiltonian matrix, can contains eigenvectors on output
         real(wp), contiguous, intent(inout) :: hmat(:, :, :)
         !> Eigenvalues
         real(wp), contiguous, intent(inout) :: eval(:, :)
         !> Occupation numbers
         real(wp), contiguous, intent(inout) :: focc(:, :)
         !> Density matrix
         real(wp), contiguous, intent(inout) :: density(:, :, :)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine get_density

      subroutine delete(self)
         import :: solver_type
         !> Solver for the general eigenvalue problem
         class(solver_type), intent(inout) :: self
      end subroutine delete
   end interface

end module tblite_scf_solver
