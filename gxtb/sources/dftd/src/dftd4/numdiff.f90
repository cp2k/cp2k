! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! Lesser GNU General Public License for more details.
!
! You should have received a copy of the Lesser GNU General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

!> Numerical differentation of DFT-D4 model
module dftd4_numdiff
   use dftd4_cutoff, only : realspace_cutoff
   use dftd4_damping, only : damping_type
   use dftd4_disp, only : get_dispersion
   use dftd4_model, only : dispersion_model
   use dftd4_param, only : param_type
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   implicit none
   private

   public :: get_dispersion_hessian


contains


!> Evaluate hessian matrix by numerical differention
subroutine get_dispersion_hessian(mol, disp, damp, param, cutoff, hessian)
   !DEC$ ATTRIBUTES DLLEXPORT :: get_dispersion_hessian

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Dispersion model
   class(dispersion_model), intent(in) :: disp

   !> Damping function
   type(damping_type), intent(in) :: damp

   !> Damping parameters
   type(param_type), intent(in) :: param

   !> Realspace cutoffs
   type(realspace_cutoff), intent(in) :: cutoff

   !> Dispersion hessian
   real(wp), intent(out) :: hessian(:, :, :, :)

   integer :: iat, ix
   real(wp), parameter :: step = 1.0e-4_wp
   type(structure_type) :: displ
   real(wp) :: el, er
   real(wp), allocatable :: gl(:, :), gr(:, :), sl(:, :), sr(:, :)

   hessian(:, :, :, :) = 0.0_wp
   !$omp parallel default(none) &
   !$omp private(iat, ix, displ, er, el, gr, gl, sr, sl) &
   !$omp shared(mol, disp, damp, param, cutoff, hessian)
   displ = mol
   allocate(gl(3, mol%nat), gr(3, mol%nat), sl(3, 3), sr(3, 3))
   !$omp do schedule(dynamic) collapse(2)
   do iat = 1, mol%nat
      do ix = 1, 3
         displ%xyz(ix, iat) = mol%xyz(ix, iat) + step
         call get_dispersion(displ, disp, damp, param, cutoff, el, gl, sl)

         displ%xyz(ix, iat) = mol%xyz(ix, iat) - step
         call get_dispersion(displ, disp, damp, param, cutoff, er, gr, sr)

         displ%xyz(ix, iat) = mol%xyz(ix, iat)
         hessian(:, :, ix, iat) = (gl - gr) / (2 * step)
      end do
   end do
   !$omp end parallel
end subroutine get_dispersion_hessian

end module dftd4_numdiff
