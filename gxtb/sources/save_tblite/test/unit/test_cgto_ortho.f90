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

module test_cgto_ortho
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use tblite_basis_cache, only : cgto_cache
   use tblite_basis_ortho, only : orthogonalize
   use tblite_basis_type, only : cgto_type
   use tblite_basis_slater, only : slater_to_gauss
   use tblite_integral_overlap, only : overlap_cgto
   implicit none
   private

   public :: collect_cgto_ortho

   real(wp), parameter :: thr = 5e+6_wp*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

contains


!> Collect all exported unit tests
subroutine collect_cgto_ortho(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("ortho-1s-2s", test_ortho_1s_2s) &
      ]

end subroutine collect_cgto_ortho


subroutine test_ortho_1s_2s(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: stat
   type(cgto_type) :: cgtoi, cgtoj
   type(cgto_cache) :: icache, jcache
   real(wp), parameter :: vec(3) = 0.0_wp, r2 = 0.0_wp
   real(wp) :: overlap(2, 2:2)

   call slater_to_gauss(6, 1, 0, 1.207940_wp, cgtoi, .true., stat)
   call check(error, stat)
   if (allocated(error)) return
   call slater_to_gauss(2, 2, 0, 1.993207_wp, cgtoj, .true., stat)
   call check(error, stat)
   if (allocated(error)) return

   call orthogonalize(cgtoi, cgtoj)

   call cgtoi%update(icache, .false.)
   call cgtoj%update(jcache, .false.)

   call overlap_cgto(cgtoj, cgtoj, jcache, jcache, r2, vec, 100.0_wp, overlap(2:2, 2))
   call check(error, overlap(2, 2), 1.0_wp, thr=thr)
   if (allocated(error)) then
      print*, (overlap(2, 2)) - 1.0_wp
      return
   end if

   call overlap_cgto(cgtoi, cgtoj, icache, jcache, r2, vec, 100.0_wp, overlap(1:1, 2))
   call check(error, overlap(1, 2), 0.0_wp, thr=thr)
   if (allocated(error)) then
      print*, (overlap(1, 2))
      return
   end if

end subroutine test_ortho_1s_2s


end module test_cgto_ortho
