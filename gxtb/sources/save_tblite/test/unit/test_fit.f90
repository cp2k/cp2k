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

module test_fit
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use tblite_fit_newuoa, only : newuoa
   use tblite_fit_settings, only : fit_settings
   use tblite_toml, only : toml_table
   implicit none
   private

   public :: collect_fit

   real(wp), parameter :: thr = 1e-6_wp


contains


!> Collect all exported unit tests
subroutine collect_fit(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("settings", test_settings), &
      new_unittest("newuoa", test_newuoa) &
      ]

end subroutine collect_fit


subroutine test_settings(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table) :: table
   type(fit_settings) :: set

   table = toml_table()

   call set%load(table, error)
   if (allocated(error)) return

   table = toml_table()

   call set%dump(table, error)
   if (allocated(error)) return
end subroutine test_settings


subroutine test_newuoa(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp) :: x(10)
   integer :: maxfun, n, npt, i
   real(wp) :: rhoend, rhobeg
   class(*), allocatable :: h

   real(wp), parameter :: ref(8, 4) = reshape([&
      & 2.113249e-1_wp, 7.886751e-1_wp, 0.000000e+0_wp, 0.000000e+0_wp, &
      & 0.000000e+0_wp, 0.000000e+0_wp, 0.000000e+0_wp, 0.000000e+0_wp, &
      & 1.026728e-1_wp, 4.062038e-1_wp, 5.937962e-1_wp, 8.973272e-1_wp, &
      & 0.000000e+0_wp, 0.000000e+0_wp, 0.000000e+0_wp, 0.000000e+0_wp, &
      & 6.687652e-2_wp, 2.887405e-1_wp, 3.666823e-1_wp, 6.333176e-1_wp, &
      & 7.112593e-1_wp, 9.331234e-1_wp, 0.000000e+0_wp, 0.000000e+0_wp, &
      & 4.315278e-2_wp, 1.930909e-1_wp, 2.663289e-1_wp, 5.000001e-1_wp, &
      & 5.000002e-1_wp, 7.336713e-1_wp, 8.069092e-1_wp, 9.568473e-1_wp],&
      & shape(ref))
   allocate(integer :: h)

   maxfun = 5000
   rhoend = 1.0e-6_wp
   do n = 2, 8, 2
      npt = 2 * n + 1
      do i = 1, n
         x(i) = real(i, wp) / real(n+1, wp)
      end do
      rhobeg = 0.2_wp * x(1)
      call newuoa(n, npt, x, rhobeg, rhoend, 0, maxfun, calfun_chebyquad, h, error)
      if (allocated(error)) exit
      do i = 1, n
         call check(error, x(i), ref(i, n/2), thr=thr)
         if (allocated(error)) exit
      end do
      if (allocated(error)) exit
   end do
   if (allocated(error)) return

end subroutine test_newuoa


!> The Chebyquad test problem (Fletcher, 1965)
function calfun_chebyquad(n, x, h, error) result(f)
   integer, intent(in) :: n
   real(wp), intent(in) :: x(*)
   class(*), intent(in) :: h
   type(error_type), allocatable, intent(out) :: error
   real(wp) :: f

   real(wp) :: y(10, 10)
   integer :: i,j,np,iw
   real(wp) :: sum

   do j = 1, n
      y(1, j) = 1.0_wp
      y(2, j) = 2.0_wp * x(j) - 1.0_wp
   end do
   do i = 2, n
      do j = 1, n
         y(i+1, j) = 2.0_wp * y(2, j) * y(i, j) - y(i-1, j)
      end do
   end do
   f = 0.0_wp
   np = n + 1
   iw = 1
   do i = 1, np
      sum = 0.0_wp
      do j = 1, n
         sum = sum + y(i, j)
      end do
      sum = sum / real(n, wp)
      if (iw > 0) sum = sum + 1.0_wp / real(i*i-2*i, wp)
      iw = - iw
      f = f + sum * sum
   end do
end function calfun_chebyquad

end module test_fit
