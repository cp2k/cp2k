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

module dftd3_utils
   use mctc_env, only : wp
   use mctc_io_math, only : matinv_3x3
   implicit none

   public :: wrap_to_central_cell


contains


subroutine wrap_to_central_cell(xyz, lattice, periodic)
   real(wp), intent(inout) :: xyz(:, :)
   real(wp), intent(in) :: lattice(:, :)
   logical, intent(in) :: periodic(:)
   real(wp) :: invlat(3, 3), vec(3)
   integer :: iat, idir

   if (.not.any(periodic)) return

   invlat = matinv_3x3(lattice)
   do iat = 1, size(xyz, 2)
      vec(:) = matmul(invlat, xyz(:, iat))
      vec(:) = shift_back_abc(vec)
      xyz(:, iat) = matmul(lattice, vec)
   end do

end subroutine wrap_to_central_cell


elemental function shift_back_abc(in) result(out)
   !> fractional coordinate in (-∞,+∞)
   real(wp),intent(in) :: in
   !> fractional coordinate in [0,1)
   real(wp) :: out
   real(wp),parameter :: p_pbc_eps = 1.0e-14_wp
   out = in
   if(in < (0.0_wp - p_pbc_eps)) &
      out = in + real(ceiling(-in),wp)
   if(in > (1.0_wp + p_pbc_eps)) &
      out = in - real(floor  ( in),wp)
   if (abs(in - 1.0_wp) < p_pbc_eps) &
      out = in - 1.0_wp
end function shift_back_abc


end module dftd3_utils
