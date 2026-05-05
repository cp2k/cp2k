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

!> @file tblite/scf/diag.f90
!> Provides a base class for defining diagonalization based electronic solvers

!> Declaration of the abstract base class for electronic solvers based on diagonalization
module tblite_scf_diag
   use mctc_env, only : sp, dp, wp, error_type
   use tblite_blas, only : gemm
   use tblite_scf_solver, only : solver_type
   use tblite_wavefunction_fermi, only : get_fermi_filling
   implicit none
   private

   !> Abstract base class for electronic solvers
   type, public, abstract, extends(solver_type) :: diag_solver_type
   contains
      generic :: solve => solve_sp, solve_dp
      procedure(solve_sp), deferred :: solve_sp
      procedure(solve_dp), deferred :: solve_dp
      procedure :: get_density
      procedure :: get_wdensity
      procedure :: delete
   end type diag_solver_type

   abstract interface
      subroutine solve_sp(self, hmat, smat, eval, error)
         import :: diag_solver_type, error_type, sp
         class(diag_solver_type), intent(inout) :: self
         real(sp), contiguous, intent(inout) :: hmat(:, :)
         real(sp), contiguous, intent(in) :: smat(:, :)
         real(sp), contiguous, intent(inout) :: eval(:)
         type(error_type), allocatable, intent(out) :: error
      end subroutine solve_sp
      subroutine solve_dp(self, hmat, smat, eval, error)
         import :: diag_solver_type, error_type, dp
         class(diag_solver_type), intent(inout) :: self
         real(dp), contiguous, intent(inout) :: hmat(:, :)
         real(dp), contiguous, intent(in) :: smat(:, :)
         real(dp), contiguous, intent(inout) :: eval(:)
         type(error_type), allocatable, intent(out) :: error
      end subroutine solve_dp
   end interface


contains

subroutine get_density(self, hmat, smat, eval, focc, density, error)
   !> Solver for the general eigenvalue problem
   class(diag_solver_type), intent(inout) :: self
   !> Overlap matrix
   real(wp), contiguous, intent(in) :: smat(:, :)
   !> Hamiltonian matrix, contains eigenvectors on output
   real(wp), contiguous, intent(inout) :: hmat(:, :, :)
   !> Eigenvalues
   real(wp), contiguous, intent(inout) :: eval(:, :)
   !> Occupation numbers
   real(wp), contiguous, intent(inout) :: focc(:, :)
   !> Density matrix
   real(wp), contiguous, intent(inout) :: density(:, :, :)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp) :: e_fermi
   integer :: nspin, spin, homo

   nspin = min(size(eval, 2), size(hmat, 3), size(density, 3))

   density(:, :, :) = 0.0_wp
   select case(nspin)
   case default
      ! Restricted case diagonalize the Fock matrix for doubly occupied orbitals
      ! and occupy based on the electronic temperature
      call self%solve(hmat(:, :, 1), smat, eval(:, 1), error)
      if (allocated(error)) return

      focc(:, :) = 0.0_wp
      do spin = 1, 2
         call get_fermi_filling(self%nel(spin), self%kt, eval(:, 1), &
            & homo, focc(:, spin), e_fermi)
      end do

      focc(:, 1) = focc(:, 1) + focc(:, 2)
      call get_density_matrix(focc(:, 1), hmat(:, :, 1), density(:, :, 1))
      focc(:, 1) = focc(:, 1) - focc(:, 2)
   case(2)
      ! Unrestricted case diagonalize the alpha and beta Fock matrix
      ! and occupy both independently based on the electronic temperature
      do spin = 1, 2
         call self%solve(hmat(:, :, spin), smat, eval(:, spin), error)
         if (allocated(error)) return

         call get_fermi_filling(self%nel(spin), self%kt, eval(:, spin), &
            & homo, focc(:, spin), e_fermi)
         call get_density_matrix(focc(:, spin), hmat(:, :, spin), density(:, :, spin))
      end do
   end select
end subroutine get_density

subroutine get_wdensity(self, hmat, smat, eval, focc, density, error)
   !> Solver for the general eigenvalue problem
   class(diag_solver_type), intent(inout) :: self
   !> Overlap matrix
   real(wp), contiguous, intent(in) :: smat(:, :)
   !> Hamiltonian matrix containing eigenvectors from SCF
   real(wp), contiguous, intent(inout) :: hmat(:, :, :)
   !> Eigenvalues
   real(wp), contiguous, intent(inout) :: eval(:, :)
   !> Occupation numbers
   real(wp), contiguous, intent(inout) :: focc(:, :)
   !> Energy weighted density matrix
   real(wp), contiguous, intent(inout) :: density(:, :, :)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: spin, nspin
   real(wp), allocatable :: tmp(:)

   nspin = min(size(eval, 2), size(hmat, 3), size(density, 3))

   if (nspin == 1 .and. size(focc, 2) == 2) then
      focc(:, 1) = focc(:, 1) + focc(:, 2)
   end if

   do spin = 1, nspin
      tmp = focc(:, spin) * eval(:, spin)
      call get_density_matrix(tmp, hmat(:, :, spin), density(:, :, spin))
   end do

   if (nspin == 1 .and. size(focc, 2) == 2) then
      focc(:, 1) = focc(:, 1) - focc(:, 2)
   end if
end subroutine get_wdensity

!> Get the density matrix from the coefficients and occupation numbers
subroutine get_density_matrix(focc, coeff, pmat)
   !> Occupation numbers
   real(wp), intent(in) :: focc(:)
   !> Coefficients of the wavefunction
   real(wp), contiguous, intent(in) :: coeff(:, :)
   !> Density matrix to be computed
   real(wp), contiguous, intent(out) :: pmat(:, :)

   real(wp), allocatable :: scratch(:, :)
   integer :: iao, jao

   allocate(scratch(size(pmat, 1), size(pmat, 2)))
   !$omp parallel do collapse(2) default(none) schedule(runtime) &
   !$omp shared(scratch, coeff, focc, pmat) private(iao, jao)
   do iao = 1, size(pmat, 1)
      do jao = 1, size(pmat, 2)
         scratch(jao, iao) = coeff(jao, iao) * focc(iao)
      end do
   end do
   call gemm(scratch, coeff, pmat, transb='t', beta=1.0_wp)
end subroutine get_density_matrix

!> Delete the solver instance
subroutine delete(self)
   !> Solver for the general eigenvalue problem
   class(diag_solver_type), intent(inout) :: self

   ! No specific resources to free in this base class
end subroutine delete

end module tblite_scf_diag