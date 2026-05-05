! This file is part of s-dftd3.
! SPDX-Identifier: LGPL-3.0-or-later
!
! s-dftd3 is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! s-dftd3 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with s-dftd3.  If not, see <https://www.gnu.org/licenses/>.

module dftd3_damping
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   implicit none

   public :: damping_param, dispersion_interface


   type, abstract :: damping_param
   contains
      procedure(dispersion_interface), deferred :: get_dispersion2
      procedure(dispersion_interface), deferred :: get_dispersion3
      procedure(pairwise_dispersion_interface), deferred :: get_pairwise_dispersion2
      procedure(pairwise_dispersion_interface), deferred :: get_pairwise_dispersion3
   end type damping_param


   abstract interface
      !> Evaluation of the dispersion energy expression
      subroutine dispersion_interface(self, mol, trans, cutoff, rvdw, r4r2, c6, dc6dcn, &
            & energy, dEdcn, gradient, sigma)
         import :: structure_type, damping_param, wp

         !> Damping parameters
         class(damping_param), intent(in) :: self

         !> Molecular structure data
         class(structure_type), intent(in) :: mol

         !> Lattice points
         real(wp), intent(in) :: trans(:, :)

         !> Real space cutoff
         real(wp), intent(in) :: cutoff

         !> Van-der-Waals radii for damping function
         real(wp), intent(in) :: rvdw(:, :)

         !> Expectation values for C8 extrapolation
         real(wp), intent(in) :: r4r2(:)

         !> C6 coefficients for all atom pairs.
         real(wp), intent(in) :: c6(:, :)

         !> Derivative of the C6 w.r.t. the coordination number
         real(wp), intent(in), optional :: dc6dcn(:, :)

         !> Dispersion energy
         real(wp), intent(inout) :: energy(:)

         !> Derivative of the energy w.r.t. the coordination number
         real(wp), intent(inout), optional :: dEdcn(:)

         !> Dispersion gradient
         real(wp), intent(inout), optional :: gradient(:, :)

         !> Dispersion virial
         real(wp), intent(inout), optional :: sigma(:, :)
      end subroutine dispersion_interface


      !> Evaluation of the pairwise representation of the dispersion energy
      subroutine pairwise_dispersion_interface(self, mol, trans, cutoff, rvdw, r4r2, c6, &
            & energy)
         import :: structure_type, damping_param, wp

         !> Damping parameters
         class(damping_param), intent(in) :: self

         !> Molecular structure data
         class(structure_type), intent(in) :: mol

         !> Lattice points
         real(wp), intent(in) :: trans(:, :)

         !> Real space cutoff
         real(wp), intent(in) :: cutoff

         !> Van-der-Waals radii for damping function
         real(wp), intent(in) :: rvdw(:, :)

         !> Expectation values for r4 over r2 operator
         real(wp), intent(in) :: r4r2(:)

         !> C6 coefficients for all atom pairs.
         real(wp), intent(in) :: c6(:, :)

         !> Pairwise representation of the dispersion energy
         real(wp), intent(inout) :: energy(:, :)
      end subroutine pairwise_dispersion_interface
   end interface


end module dftd3_damping
