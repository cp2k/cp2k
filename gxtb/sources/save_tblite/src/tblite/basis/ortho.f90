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

!> @file tblite/basis/ortho.f90
!> Provides utilities for orthonormalizing basis functions

!> Gram-Schmidt orthonormalization routines for contracted Gaussian basis functions
module tblite_basis_ortho
   use mctc_env, only : wp
   use mctc_io_constants, only : pi
   use tblite_basis_type, only : cgto_type
   implicit none
   private

   public :: orthogonalize


contains

!> Orthogonalize a contracted Gaussian basis function to an existing basis function
pure subroutine orthogonalize(cgtoi, cgtoj)
   !> Existing basis function
   class(cgto_type), intent(in) :: cgtoi
   !> Basis function to orthogonalize
   class(cgto_type), intent(inout) :: cgtoj

   integer :: ipr, jpr
   real(wp) :: eab, oab, kab, overlap

   if (cgtoi%ang /= cgtoj%ang) return

   overlap = 0.0_wp
   do ipr = 1, cgtoi%nprim
      do jpr = 1, cgtoj%nprim
         eab = cgtoi%alpha(ipr) + cgtoj%alpha(jpr)
         oab = 1.0_wp / eab
         kab = sqrt(pi*oab)**3

         overlap = overlap + cgtoi%coeff(ipr) * cgtoj%coeff(jpr) * kab
      end do
   end do

   cgtoj%alpha(cgtoj%nprim+1:cgtoj%nprim+cgtoi%nprim) = cgtoi%alpha(:cgtoi%nprim)
   cgtoj%coeff(cgtoj%nprim+1:cgtoj%nprim+cgtoi%nprim) = -overlap*cgtoi%coeff(:cgtoi%nprim)
   cgtoj%nprim = cgtoj%nprim + cgtoi%nprim

   overlap = 0.0_wp
   do ipr = 1, cgtoj%nprim
      do jpr = 1, cgtoj%nprim
         eab = cgtoj%alpha(ipr) + cgtoj%alpha(jpr)
         oab = 1.0_wp / eab
         kab = sqrt(pi*oab)**3

         overlap = overlap + cgtoj%coeff(ipr) * cgtoj%coeff(jpr) * kab
      end do
   end do

   cgtoj%coeff(:cgtoj%nprim) = cgtoj%coeff(:cgtoj%nprim) / sqrt(overlap)

end subroutine orthogonalize

end module tblite_basis_ortho
