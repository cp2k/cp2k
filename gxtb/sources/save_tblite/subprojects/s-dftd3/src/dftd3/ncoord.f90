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

module dftd3_ncoord
   use, intrinsic :: iso_fortran_env, only : error_unit
   use mctc_env, only : error_type, wp
   use mctc_io, only : structure_type
   use mctc_ncoord, only : ncoord_type, new_ncoord, cn_count
   implicit none
   private

   public :: get_coordination_number, add_coordination_number_derivs

   !> Steepness of counting function
   real(wp), parameter :: default_kcn = 16.0_wp

contains


!> Wrapper for geometric fractional coordination number 
!> with standard exponential counting function.
subroutine get_coordination_number(mol, trans, cutoff, rcov, cn, dcndr, dcndL)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Covalent radius
   real(wp), intent(in) :: rcov(:)

   !> Error function coordination number.
   real(wp), intent(out) :: cn(:)

   !> Derivative of the CN with respect to the Cartesian coordinates.
   real(wp), intent(out), optional :: dcndr(:, :, :)

   !> Derivative of the CN with respect to strain deformations.
   real(wp), intent(out), optional :: dcndL(:, :, :)

   class(ncoord_type), allocatable :: ncoord
   type(error_type), allocatable :: error

   call new_ncoord(ncoord, mol, cn_count%exp, error, &
      & kcn=default_kcn, cutoff=cutoff, rcov=rcov)
   if(allocated(error)) then
      write(error_unit, '("[Error]:", 1x, a)') error%message
      error stop
   end if

   call ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)

end subroutine get_coordination_number


subroutine add_coordination_number_derivs(mol, trans, cutoff, rcov, dEdcn, gradient, sigma)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Covalent radius
   real(wp), intent(in) :: rcov(:)

   !> Derivative of expression with respect to the coordination number
   real(wp), intent(in) :: dEdcn(:)

   !> Derivative of the CN with respect to the Cartesian coordinates
   real(wp), intent(inout) :: gradient(:, :)

   !> Derivative of the CN with respect to strain deformations
   real(wp), intent(inout) :: sigma(:, :)

   class(ncoord_type), allocatable :: ncoord
   type(error_type), allocatable :: error

   call new_ncoord(ncoord, mol, cn_count%exp, error, &
      & kcn=default_kcn, cutoff=cutoff, rcov=rcov)
   if(allocated(error)) then
      write(error_unit, '("[Error]:", 1x, a)') error%message
      error stop
   end if

   call ncoord%add_coordination_number_derivs(mol, trans, dEdcn, gradient, sigma)

end subroutine add_coordination_number_derivs


end module dftd3_ncoord