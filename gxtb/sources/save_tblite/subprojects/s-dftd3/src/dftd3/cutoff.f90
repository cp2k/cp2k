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

module dftd3_cutoff
   use mctc_env, only : wp
   implicit none

   public :: realspace_cutoff, get_lattice_points


   !> Coordination number cutoff
   real(wp), parameter :: cn_default = 40.0_wp

   !> Two-body interaction cutoff
   real(wp), parameter :: disp2_default = 60.0_wp

   !> Three-body interaction cutoff
   real(wp), parameter :: disp3_default = 40.0_wp

   !> Counter-poise correction cutoff
   real(wp), parameter :: gcp_default = 60.0_wp

   !> Short-range bond correction cutoff
   real(wp), parameter :: srb_default = 60.0_wp


   !> Collection of real space cutoffs
   type :: realspace_cutoff
      sequence

      !> Coordination number cutoff
      real(wp) :: cn = cn_default

      !> Two-body interaction cutoff
      real(wp) :: disp2 = disp2_default

      !> Three-body interaction cutoff
      real(wp) :: disp3 = disp3_default

      !> Counter-poise correction cutoff
      real(wp) :: gcp = gcp_default

      !> Short-range bond correction cutoff
      real(wp) :: srb = srb_default

   end type realspace_cutoff


contains


subroutine get_lattice_points(periodic, lat, rthr, trans)
   logical, intent(in) :: periodic(:)
   real(wp), intent(in) :: rthr
   real(wp), intent(in) :: lat(:, :)
   real(wp), allocatable, intent(out) :: trans(:, :)
   real(wp) :: vec(size(lat, 2)), norms(size(lat, 1), size(lat, 2))
   integer :: rep(3)
   integer :: itr, ix, iy, iz

   if (.not.any(periodic)) then
      allocate(trans(3, 1))
      trans(:, :) = 0.0_wp
      return
   end if

   if (all(periodic)) then
      call get_translations(lat, rthr, rep)
   else if (count(periodic) == 1) then
      vec = sqrt(sum(lat**2, 1))
      where(periodic)
         rep = ceiling(rthr / vec)
      elsewhere
         rep = 0
      endwhere
   else if (count(periodic) == 2) then
      call get_normals(lat, norms)
      where(spread(periodic, 2, 3))
         norms = lat
      endwhere
      call get_translations(norms, rthr, rep)
      where(.not.periodic)
         rep = 0
      endwhere
   end if

   allocate(trans(3, product(2*rep+1)))
   itr = 0
   do ix = -rep(1), rep(1)
      do iy = -rep(2), rep(2)
         do iz = -rep(3), rep(3)
            itr = itr + 1
            trans(:, itr) = lat(:, 1)*ix + lat(:, 2)*iy + lat(:, 3)*iz
         end do
      end do
   end do

end subroutine get_lattice_points


subroutine get_normals(lattice, normal)
   real(wp), intent(in) :: lattice(:, :)
   real(wp), intent(out) :: normal(:, :)

   integer :: itr
   do itr = 1, 3
      call crossproduct(lattice(:, mod(itr, 3) + 1), lattice(:, mod(itr + 1, 3) + 1), &
         & normal(:, itr))
   end do
end subroutine get_normals


!> generate a supercell based on a realspace cutoff, this subroutine
!  doesn't know anything about the convergence behaviour of the
!  associated property.
pure subroutine get_translations(lat, rthr, rep)
   real(wp), intent(in)  :: rthr
   real(wp), intent(in)  :: lat(3, 3)
   integer, intent(out) :: rep(3)
   real(wp) :: norm(3, 3), normy(3), normz(3)
   real(wp) :: cos10, cos21, cos32

   ! find normal to the plane...
   call crossproduct(lat(:, 2), lat(:, 3), norm(:, 1))
   call crossproduct(lat(:, 3), lat(:, 1), norm(:, 2))
   call crossproduct(lat(:, 1), lat(:, 2), norm(:, 3))
   ! ...normalize it...
   norm(:, 1) = norm(:, 1)/norm2(norm(:, 1))
   norm(:, 2) = norm(:, 2)/norm2(norm(:, 2))
   norm(:, 3) = norm(:, 3)/norm2(norm(:, 3))
   ! cos angles between normals and lattice vectors
   cos10 = sum(norm(:, 1)*lat(:, 1))
   cos21 = sum(norm(:, 2)*lat(:, 2))
   cos32 = sum(norm(:, 3)*lat(:, 3))
   rep(1) = ceiling(abs(rthr/cos10))
   rep(2) = ceiling(abs(rthr/cos21))
   rep(3) = ceiling(abs(rthr/cos32))

end subroutine get_translations

pure subroutine crossproduct(a, b, c)
   real(wp), intent(in)  :: a(3), b(3)
   real(wp), intent(out) :: c(3)
   real(wp) :: x, y, z

   x = a(2)*b(3)-b(2)*a(3)
   y = a(3)*b(1)-b(3)*a(1)
   z = a(1)*b(2)-b(1)*a(2)
   c = [x, y, z]
end subroutine crossproduct

end module dftd3_cutoff
