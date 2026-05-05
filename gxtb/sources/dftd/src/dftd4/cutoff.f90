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

!> Realspace cutoff and lattice point generator utilities
module dftd4_cutoff
   use mctc_env, only : wp
   implicit none
   private

   public :: realspace_cutoff, get_lattice_points


   !> Coordination number cutoff
   real(wp), parameter :: cn_default = 30.0_wp

   !> Two-body interaction cutoff
   real(wp), parameter :: disp2_default = 60.0_wp

   !> Three-body interaction cutoff
   real(wp), parameter :: disp3_default = 40.0_wp


   !> Collection of real space cutoffs
   type :: realspace_cutoff
      sequence

      !> Coordination number cutoff
      real(wp) :: cn = cn_default

      !> Two-body interaction cutoff
      real(wp) :: disp2 = disp2_default

      !> Three-body interaction cutoff
      real(wp) :: disp3 = disp3_default

   end type realspace_cutoff


   interface get_lattice_points
      module procedure :: get_lattice_points_cutoff
      module procedure :: get_lattice_points_rep_3d
   end interface get_lattice_points


contains


!> Generate lattice points from repeatitions
subroutine get_lattice_points_rep_3d(lat, rep, origin, trans)

   !> Lattice vectors
   real(wp), intent(in) :: lat(:, :)

   !> Repeatitions of lattice points to generate
   integer, intent(in) :: rep(:)

   !> Include the origin in the generated lattice points
   logical, intent(in) :: origin

   !> Generated lattice points
   real(wp), allocatable, intent(out) :: trans(:, :)

   integer :: itr, ix, iy, iz, jx, jy, jz

   itr = 0
   if (origin) then
      allocate(trans(3, product(2*rep+1)))
      do ix = 0, rep(1)
         do iy = 0, rep(2)
            do iz = 0, rep(3)
               do jx = 1, merge(-1, 1, ix > 0), -2
                  do jy = 1, merge(-1, 1, iy > 0), -2
                     do jz = 1, merge(-1, 1, iz > 0), -2
                        itr = itr + 1
                        trans(:, itr) = lat(:, 1)*ix*jx &
                           & + lat(:, 2)*iy*jy + lat(:, 3)*iz*jz
                     end do
                  end do
               end do
            end do
         end do
      end do
   else
      allocate(trans(3, product(2*rep+1)-1))
      do ix = 0, rep(1)
         do iy = 0, rep(2)
            do iz = 0, rep(3)
               if (ix == 0 .and. iy == 0 .and. iz == 0) cycle
               do jx = 1, merge(-1, 1, ix > 0), -2
                  do jy = 1, merge(-1, 1, iy > 0), -2
                     do jz = 1, merge(-1, 1, iz > 0), -2
                        itr = itr + 1
                        trans(:, itr) = lat(:, 1)*ix*jx &
                           & + lat(:, 2)*iy*jy + lat(:, 3)*iz*jz
                     end do
                  end do
               end do
            end do
         end do
      end do
   end if

end subroutine get_lattice_points_rep_3d


!> Create lattice points within a given cutoff
subroutine get_lattice_points_cutoff(periodic, lat, rthr, trans)
   !DEC$ ATTRIBUTES DLLEXPORT :: get_lattice_points_cutoff

   !> Periodic dimensions
   logical, intent(in) :: periodic(:)

   !> Real space cutoff
   real(wp), intent(in) :: rthr

   !> Lattice parameters
   real(wp), intent(in) :: lat(:, :)

   !> Generated lattice points
   real(wp), allocatable, intent(out) :: trans(:, :)

   integer :: rep(3)

   if (.not.any(periodic)) then
      allocate(trans(3, 1))
      trans(:, :) = 0.0_wp
   else
      call get_translations(lat, rthr, rep)
      call get_lattice_points(lat, rep, .true., trans)
   end if

end subroutine get_lattice_points_cutoff


!> Generate a supercell based on a realspace cutoff, this subroutine
!> doesn't know anything about the convergence behaviour of the
!> associated property.
pure subroutine get_translations(lat, rthr, rep)
   real(wp), intent(in) :: rthr
   real(wp), intent(in) :: lat(3, 3)
   integer, intent(out) :: rep(3)
   real(wp) :: normx(3), normy(3), normz(3)
   real(wp) :: cos10, cos21, cos32

   ! find normal to the plane...
   call crossproduct(lat(:, 2), lat(:, 3), normx)
   call crossproduct(lat(:, 3), lat(:, 1), normy)
   call crossproduct(lat(:, 1), lat(:, 2), normz)
   ! ...normalize it...
   normx = normx/norm2(normx)
   normy = normy/norm2(normy)
   normz = normz/norm2(normz)
   ! cos angles between normals and lattice vectors
   cos10 = sum(normx*lat(:, 1))
   cos21 = sum(normy*lat(:, 2))
   cos32 = sum(normz*lat(:, 3))
   rep(1) = ceiling(abs(rthr/cos10))
   rep(2) = ceiling(abs(rthr/cos21))
   rep(3) = ceiling(abs(rthr/cos32))

contains

   pure subroutine crossproduct(a, b, c)
      real(wp), intent(in) :: a(3)
      real(wp), intent(in) :: b(3)
      real(wp), intent(out) :: c(3)
      c(1)=a(2)*b(3)-b(2)*a(3)
      c(2)=a(3)*b(1)-b(3)*a(1)
      c(3)=a(1)*b(2)-b(1)*a(2)
   end subroutine crossproduct

end subroutine get_translations


end module dftd4_cutoff
