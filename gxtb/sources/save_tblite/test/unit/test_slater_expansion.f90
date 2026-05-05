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

module test_slater_expansion
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use tblite_basis_cache, only : cgto_cache
   use tblite_basis_type, only : cgto_type
   use tblite_basis_slater, only : slater_to_gauss
   use tblite_integral_overlap, only : overlap_cgto
   implicit none
   private

   public :: collect_slater_expansion

   real(wp), parameter :: thr = 5e+6_wp*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

contains


!> Collect all exported unit tests
subroutine collect_slater_expansion(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("norm-1s-sto1g", test_norm_1s_sto1g), &
      new_unittest("norm-2s-sto1g", test_norm_2s_sto1g), &
      new_unittest("norm-3s-sto1g", test_norm_3s_sto1g), &
      new_unittest("norm-4s-sto1g", test_norm_4s_sto1g), &
      new_unittest("norm-5s-sto1g", test_norm_5s_sto1g), &
      new_unittest("norm-2p-sto1g", test_norm_2p_sto1g), &
      new_unittest("norm-3p-sto1g", test_norm_3p_sto1g), &
      new_unittest("norm-4p-sto1g", test_norm_4p_sto1g), &
      new_unittest("norm-5p-sto1g", test_norm_5p_sto1g), &
      new_unittest("norm-3d-sto1g", test_norm_3d_sto1g), &
      new_unittest("norm-4d-sto1g", test_norm_4d_sto1g), &
      new_unittest("norm-5d-sto1g", test_norm_5d_sto1g), &
      new_unittest("norm-4f-sto1g", test_norm_4f_sto1g), &
      new_unittest("norm-5f-sto1g", test_norm_5f_sto1g), &
      new_unittest("norm-5g-sto1g", test_norm_5g_sto1g), &
      new_unittest("norm-1s-sto2g", test_norm_1s_sto2g), &
      new_unittest("norm-2s-sto2g", test_norm_2s_sto2g), &
      new_unittest("norm-3s-sto2g", test_norm_3s_sto2g), &
      new_unittest("norm-4s-sto2g", test_norm_4s_sto2g), &
      new_unittest("norm-5s-sto2g", test_norm_5s_sto2g), &
      new_unittest("norm-2p-sto2g", test_norm_2p_sto2g), &
      new_unittest("norm-3p-sto2g", test_norm_3p_sto2g), &
      new_unittest("norm-4p-sto2g", test_norm_4p_sto2g), &
      new_unittest("norm-5p-sto2g", test_norm_5p_sto2g), &
      new_unittest("norm-3d-sto2g", test_norm_3d_sto2g), &
      new_unittest("norm-4d-sto2g", test_norm_4d_sto2g), &
      new_unittest("norm-5d-sto2g", test_norm_5d_sto2g), &
      new_unittest("norm-4f-sto2g", test_norm_4f_sto2g), &
      new_unittest("norm-5f-sto2g", test_norm_5f_sto2g), &
      new_unittest("norm-5g-sto2g", test_norm_5g_sto2g), &
      new_unittest("norm-1s-sto3g", test_norm_1s_sto3g), &
      new_unittest("norm-2s-sto3g", test_norm_2s_sto3g), &
      new_unittest("norm-3s-sto3g", test_norm_3s_sto3g), &
      new_unittest("norm-4s-sto3g", test_norm_4s_sto3g), &
      new_unittest("norm-5s-sto3g", test_norm_5s_sto3g), &
      new_unittest("norm-2p-sto3g", test_norm_2p_sto3g), &
      new_unittest("norm-3p-sto3g", test_norm_3p_sto3g), &
      new_unittest("norm-4p-sto3g", test_norm_4p_sto3g), &
      new_unittest("norm-5p-sto3g", test_norm_5p_sto3g), &
      new_unittest("norm-3d-sto3g", test_norm_3d_sto3g), &
      new_unittest("norm-4d-sto3g", test_norm_4d_sto3g), &
      new_unittest("norm-5d-sto3g", test_norm_5d_sto3g), &
      new_unittest("norm-4f-sto3g", test_norm_4f_sto3g), &
      new_unittest("norm-5f-sto3g", test_norm_5f_sto3g), &
      new_unittest("norm-5g-sto3g", test_norm_5g_sto3g), &
      new_unittest("norm-1s-sto4g", test_norm_1s_sto4g), &
      new_unittest("norm-2s-sto4g", test_norm_2s_sto4g), &
      new_unittest("norm-3s-sto4g", test_norm_3s_sto4g), &
      new_unittest("norm-4s-sto4g", test_norm_4s_sto4g), &
      new_unittest("norm-5s-sto4g", test_norm_5s_sto4g), &
      new_unittest("norm-2p-sto4g", test_norm_2p_sto4g), &
      new_unittest("norm-3p-sto4g", test_norm_3p_sto4g), &
      new_unittest("norm-4p-sto4g", test_norm_4p_sto4g), &
      new_unittest("norm-5p-sto4g", test_norm_5p_sto4g), &
      new_unittest("norm-3d-sto4g", test_norm_3d_sto4g), &
      new_unittest("norm-4d-sto4g", test_norm_4d_sto4g), &
      new_unittest("norm-5d-sto4g", test_norm_5d_sto4g), &
      new_unittest("norm-4f-sto4g", test_norm_4f_sto4g), &
      new_unittest("norm-5f-sto4g", test_norm_5f_sto4g), &
      new_unittest("norm-5g-sto4g", test_norm_5g_sto4g), &
      new_unittest("norm-1s-sto5g", test_norm_1s_sto5g), &
      new_unittest("norm-2s-sto5g", test_norm_2s_sto5g), &
      new_unittest("norm-3s-sto5g", test_norm_3s_sto5g), &
      new_unittest("norm-4s-sto5g", test_norm_4s_sto5g), &
      new_unittest("norm-5s-sto5g", test_norm_5s_sto5g), &
      new_unittest("norm-2p-sto5g", test_norm_2p_sto5g), &
      new_unittest("norm-3p-sto5g", test_norm_3p_sto5g), &
      new_unittest("norm-4p-sto5g", test_norm_4p_sto5g), &
      new_unittest("norm-5p-sto5g", test_norm_5p_sto5g), &
      new_unittest("norm-3d-sto5g", test_norm_3d_sto5g), &
      new_unittest("norm-4d-sto5g", test_norm_4d_sto5g), &
      new_unittest("norm-5d-sto5g", test_norm_5d_sto5g), &
      new_unittest("norm-4f-sto5g", test_norm_4f_sto5g), &
      new_unittest("norm-5f-sto5g", test_norm_5f_sto5g), &
      new_unittest("norm-5g-sto5g", test_norm_5g_sto5g), &
      new_unittest("norm-1s-sto6g", test_norm_1s_sto6g), &
      new_unittest("norm-2s-sto6g", test_norm_2s_sto6g), &
      new_unittest("norm-3s-sto6g", test_norm_3s_sto6g), &
      new_unittest("norm-4s-sto6g", test_norm_4s_sto6g), &
      new_unittest("norm-5s-sto6g", test_norm_5s_sto6g), &
      new_unittest("norm-6s-sto6g", test_norm_6s_sto6g), &
      new_unittest("norm-2p-sto6g", test_norm_2p_sto6g), &
      new_unittest("norm-3p-sto6g", test_norm_3p_sto6g), &
      new_unittest("norm-4p-sto6g", test_norm_4p_sto6g), &
      new_unittest("norm-5p-sto6g", test_norm_5p_sto6g), &
      new_unittest("norm-6p-sto6g", test_norm_6p_sto6g), &
      new_unittest("norm-3d-sto6g", test_norm_3d_sto6g), &
      new_unittest("norm-4d-sto6g", test_norm_4d_sto6g), &
      new_unittest("norm-5d-sto6g", test_norm_5d_sto6g), &
      new_unittest("norm-4f-sto6g", test_norm_4f_sto6g), &
      new_unittest("norm-5f-sto6g", test_norm_5f_sto6g), &
      new_unittest("norm-5g-sto6g", test_norm_5g_sto6g) &
      ]

end subroutine collect_slater_expansion


subroutine test_norm_s(error, n, ng)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, intent(in) :: ng, n

   integer :: stat
   type(cgto_type) :: cgto
   type(cgto_cache) :: cache
   real(wp), parameter :: vec(3) = 0.0_wp, r2 = 0.0_wp
   real(wp) :: overlap(1, 1)

   call slater_to_gauss(ng, n, 0, 1.0_wp, cgto, .true., stat)
   call check(error, stat)
   if (allocated(error)) return

   call cgto%update(cache, .false.)

   call overlap_cgto(cgto, cgto, cache, cache, r2, vec, 100.0_wp, overlap)

   call check(error, overlap(1, 1), 1.0_wp, thr=thr)
   if (allocated(error)) then
      print*, (overlap(1, 1) - 1.0_wp)
   end if

end subroutine test_norm_s


subroutine test_norm_1s_sto1g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_s(error, 1, 1)
end subroutine test_norm_1s_sto1g

subroutine test_norm_2s_sto1g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_s(error, 2, 1)
end subroutine test_norm_2s_sto1g

subroutine test_norm_3s_sto1g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_s(error, 3, 1)
end subroutine test_norm_3s_sto1g

subroutine test_norm_4s_sto1g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_s(error, 4, 1)
end subroutine test_norm_4s_sto1g

subroutine test_norm_5s_sto1g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_s(error, 5, 1)
end subroutine test_norm_5s_sto1g


subroutine test_norm_1s_sto2g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_s(error, 1, 2)
end subroutine test_norm_1s_sto2g

subroutine test_norm_2s_sto2g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_s(error, 2, 2)
end subroutine test_norm_2s_sto2g

subroutine test_norm_3s_sto2g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_s(error, 3, 2)
end subroutine test_norm_3s_sto2g

subroutine test_norm_4s_sto2g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_s(error, 4, 2)
end subroutine test_norm_4s_sto2g

subroutine test_norm_5s_sto2g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_s(error, 5, 2)
end subroutine test_norm_5s_sto2g


subroutine test_norm_1s_sto3g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_s(error, 1, 3)
end subroutine test_norm_1s_sto3g

subroutine test_norm_2s_sto3g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_s(error, 2, 3)
end subroutine test_norm_2s_sto3g

subroutine test_norm_3s_sto3g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_s(error, 3, 3)
end subroutine test_norm_3s_sto3g

subroutine test_norm_4s_sto3g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_s(error, 4, 3)
end subroutine test_norm_4s_sto3g

subroutine test_norm_5s_sto3g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_s(error, 5, 3)
end subroutine test_norm_5s_sto3g


subroutine test_norm_1s_sto4g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_s(error, 1, 4)
end subroutine test_norm_1s_sto4g

subroutine test_norm_2s_sto4g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_s(error, 2, 4)
end subroutine test_norm_2s_sto4g

subroutine test_norm_3s_sto4g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_s(error, 3, 4)
end subroutine test_norm_3s_sto4g

subroutine test_norm_4s_sto4g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_s(error, 4, 4)
end subroutine test_norm_4s_sto4g

subroutine test_norm_5s_sto4g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_s(error, 5, 4)
end subroutine test_norm_5s_sto4g


subroutine test_norm_1s_sto5g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_s(error, 1, 5)
end subroutine test_norm_1s_sto5g

subroutine test_norm_2s_sto5g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_s(error, 2, 5)
end subroutine test_norm_2s_sto5g

subroutine test_norm_3s_sto5g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_s(error, 3, 5)
end subroutine test_norm_3s_sto5g

subroutine test_norm_4s_sto5g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_s(error, 4, 5)
end subroutine test_norm_4s_sto5g

subroutine test_norm_5s_sto5g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_s(error, 5, 5)
end subroutine test_norm_5s_sto5g


subroutine test_norm_1s_sto6g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_s(error, 1, 6)
end subroutine test_norm_1s_sto6g

subroutine test_norm_2s_sto6g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_s(error, 2, 6)
end subroutine test_norm_2s_sto6g

subroutine test_norm_3s_sto6g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_s(error, 3, 6)
end subroutine test_norm_3s_sto6g

subroutine test_norm_4s_sto6g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_s(error, 4, 6)
end subroutine test_norm_4s_sto6g

subroutine test_norm_5s_sto6g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_s(error, 5, 6)
end subroutine test_norm_5s_sto6g

subroutine test_norm_6s_sto6g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_s(error, 6, 6)
end subroutine test_norm_6s_sto6g


subroutine test_norm_p(error, n, ng)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, intent(in) :: ng, n

   integer :: stat, i, j
   type(cgto_type) :: cgto
   type(cgto_cache) :: cache
   real(wp), parameter :: vec(3) = 0.0_wp, r2 = 0.0_wp
   real(wp) :: overlap(3, 3)

   call slater_to_gauss(ng, n, 1, 1.0_wp, cgto, .true., stat)
   call check(error, stat)
   if (allocated(error)) return

   call cgto%update(cache, .false.)

   call overlap_cgto(cgto, cgto, cache, cache, r2, vec, 100.0_wp, overlap)

   lp: do i = 1, 3
      do j = 1, 3
         call check(error, overlap(j, i), merge(1.0_wp, 0.0_wp, j == i), thr=thr)
         if (allocated(error)) exit lp
      end do
   end do lp
   if (allocated(error)) then
      print '(3es20.14e1)', overlap - reshape([1,0,0, 0,1,0, 0,0,1], [3,3])
   end if

end subroutine test_norm_p


subroutine test_norm_2p_sto1g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_p(error, 2, 1)
end subroutine test_norm_2p_sto1g

subroutine test_norm_3p_sto1g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_p(error, 3, 1)
end subroutine test_norm_3p_sto1g

subroutine test_norm_4p_sto1g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_p(error, 4, 1)
end subroutine test_norm_4p_sto1g

subroutine test_norm_5p_sto1g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_p(error, 5, 1)
end subroutine test_norm_5p_sto1g


subroutine test_norm_2p_sto2g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_p(error, 2, 2)
end subroutine test_norm_2p_sto2g

subroutine test_norm_3p_sto2g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_p(error, 3, 2)
end subroutine test_norm_3p_sto2g

subroutine test_norm_4p_sto2g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_p(error, 4, 2)
end subroutine test_norm_4p_sto2g

subroutine test_norm_5p_sto2g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_p(error, 5, 2)
end subroutine test_norm_5p_sto2g


subroutine test_norm_2p_sto3g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_p(error, 2, 3)
end subroutine test_norm_2p_sto3g

subroutine test_norm_3p_sto3g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_p(error, 3, 3)
end subroutine test_norm_3p_sto3g

subroutine test_norm_4p_sto3g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_p(error, 4, 3)
end subroutine test_norm_4p_sto3g

subroutine test_norm_5p_sto3g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_p(error, 5, 3)
end subroutine test_norm_5p_sto3g


subroutine test_norm_2p_sto4g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_p(error, 2, 4)
end subroutine test_norm_2p_sto4g

subroutine test_norm_3p_sto4g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_p(error, 3, 4)
end subroutine test_norm_3p_sto4g

subroutine test_norm_4p_sto4g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_p(error, 4, 4)
end subroutine test_norm_4p_sto4g

subroutine test_norm_5p_sto4g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_p(error, 5, 4)
end subroutine test_norm_5p_sto4g


subroutine test_norm_2p_sto5g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_p(error, 2, 5)
end subroutine test_norm_2p_sto5g

subroutine test_norm_3p_sto5g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_p(error, 3, 5)
end subroutine test_norm_3p_sto5g

subroutine test_norm_4p_sto5g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_p(error, 4, 5)
end subroutine test_norm_4p_sto5g

subroutine test_norm_5p_sto5g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_p(error, 5, 5)
end subroutine test_norm_5p_sto5g


subroutine test_norm_2p_sto6g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_p(error, 2, 6)
end subroutine test_norm_2p_sto6g

subroutine test_norm_3p_sto6g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_p(error, 3, 6)
end subroutine test_norm_3p_sto6g

subroutine test_norm_4p_sto6g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_p(error, 4, 6)
end subroutine test_norm_4p_sto6g

subroutine test_norm_5p_sto6g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_p(error, 5, 6)
end subroutine test_norm_5p_sto6g

subroutine test_norm_6p_sto6g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_p(error, 6, 6)
end subroutine test_norm_6p_sto6g


subroutine test_norm_d(error, n, ng)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, intent(in) :: ng, n

   integer :: stat, i, j
   type(cgto_type) :: cgto
   type(cgto_cache) :: cache
   real(wp), parameter :: vec(3) = 0.0_wp, r2 = 0.0_wp
   real(wp) :: overlap(5, 5)

   call slater_to_gauss(ng, n, 2, 1.0_wp, cgto, .true., stat)
   call check(error, stat)
   if (allocated(error)) return

   call cgto%update(cache, .false.)

   call overlap_cgto(cgto, cgto, cache, cache, r2, vec, 100.0_wp, overlap)

   lp: do i = 1, 5
      do j = 1, 5
         call check(error, overlap(j, i), merge(1.0_wp, 0.0_wp, j == i), thr=thr)
         if (allocated(error)) exit lp
      end do
   end do lp
   if (allocated(error)) then
      print '(5es20.14e1)', overlap &
         & - reshape([1,0,0,0,0, 0,1,0,0,0, 0,0,1,0,0, &
         &            0,0,0,1,0, 0,0,0,0,1], [5,5])
   end if

end subroutine test_norm_d


subroutine test_norm_3d_sto1g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_d(error, 3, 1)
end subroutine test_norm_3d_sto1g

subroutine test_norm_4d_sto1g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_d(error, 4, 1)
end subroutine test_norm_4d_sto1g

subroutine test_norm_5d_sto1g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_d(error, 5, 1)
end subroutine test_norm_5d_sto1g


subroutine test_norm_3d_sto2g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_d(error, 3, 2)
end subroutine test_norm_3d_sto2g

subroutine test_norm_4d_sto2g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_d(error, 4, 2)
end subroutine test_norm_4d_sto2g

subroutine test_norm_5d_sto2g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_d(error, 5, 2)
end subroutine test_norm_5d_sto2g


subroutine test_norm_3d_sto3g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_d(error, 3, 3)
end subroutine test_norm_3d_sto3g

subroutine test_norm_4d_sto3g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_d(error, 4, 3)
end subroutine test_norm_4d_sto3g

subroutine test_norm_5d_sto3g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_d(error, 5, 3)
end subroutine test_norm_5d_sto3g


subroutine test_norm_3d_sto4g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_d(error, 3, 4)
end subroutine test_norm_3d_sto4g

subroutine test_norm_4d_sto4g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_d(error, 4, 4)
end subroutine test_norm_4d_sto4g

subroutine test_norm_5d_sto4g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_d(error, 5, 4)
end subroutine test_norm_5d_sto4g


subroutine test_norm_3d_sto5g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_d(error, 3, 5)
end subroutine test_norm_3d_sto5g

subroutine test_norm_4d_sto5g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_d(error, 4, 5)
end subroutine test_norm_4d_sto5g

subroutine test_norm_5d_sto5g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_d(error, 5, 5)
end subroutine test_norm_5d_sto5g


subroutine test_norm_3d_sto6g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_d(error, 3, 6)
end subroutine test_norm_3d_sto6g

subroutine test_norm_4d_sto6g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_d(error, 4, 6)
end subroutine test_norm_4d_sto6g

subroutine test_norm_5d_sto6g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_d(error, 5, 6)
end subroutine test_norm_5d_sto6g


subroutine test_norm_f(error, n, ng)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, intent(in) :: ng, n

   integer :: stat, i, j
   type(cgto_type) :: cgto
   type(cgto_cache) :: cache
   real(wp), parameter :: vec(3) = 0.0_wp, r2 = 0.0_wp
   real(wp) :: overlap(7, 7)

   call slater_to_gauss(ng, n, 3, 1.0_wp, cgto, .true., stat)
   call check(error, stat)
   if (allocated(error)) return

   call cgto%update(cache, .false.)

   call overlap_cgto(cgto, cgto, cache, cache, r2, vec, 100.0_wp, overlap)

   lp: do i = 1, 7
      do j = 1, 7
         call check(error, overlap(j, i), merge(1.0_wp, 0.0_wp, j == i), thr=thr)
         if (allocated(error)) exit lp
      end do
   end do lp
   if (allocated(error)) then
      print '(7es20.14e1)', overlap &
         & - reshape([1,0,0,0,0,0,0, 0,1,0,0,0,0,0, 0,0,1,0,0,0,0, &
         &            0,0,0,1,0,0,0, 0,0,0,0,1,0,0, 0,0,0,0,0,1,0, &
         &            0,0,0,0,0,0,1], [7,7])
   end if

end subroutine test_norm_f


subroutine test_norm_4f_sto1g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_f(error, 4, 1)
end subroutine test_norm_4f_sto1g

subroutine test_norm_5f_sto1g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_f(error, 5, 1)
end subroutine test_norm_5f_sto1g


subroutine test_norm_4f_sto2g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_f(error, 4, 2)
end subroutine test_norm_4f_sto2g

subroutine test_norm_5f_sto2g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_f(error, 5, 2)
end subroutine test_norm_5f_sto2g


subroutine test_norm_4f_sto3g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_f(error, 4, 3)
end subroutine test_norm_4f_sto3g

subroutine test_norm_5f_sto3g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_f(error, 5, 3)
end subroutine test_norm_5f_sto3g


subroutine test_norm_4f_sto4g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_f(error, 4, 4)
end subroutine test_norm_4f_sto4g

subroutine test_norm_5f_sto4g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_f(error, 5, 4)
end subroutine test_norm_5f_sto4g


subroutine test_norm_4f_sto5g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_f(error, 4, 5)
end subroutine test_norm_4f_sto5g

subroutine test_norm_5f_sto5g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_f(error, 5, 5)
end subroutine test_norm_5f_sto5g


subroutine test_norm_4f_sto6g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_f(error, 4, 6)
end subroutine test_norm_4f_sto6g

subroutine test_norm_5f_sto6g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_f(error, 5, 6)
end subroutine test_norm_5f_sto6g


subroutine test_norm_g(error, n, ng)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, intent(in) :: ng, n

   integer :: stat, i, j
   type(cgto_type) :: cgto
   type(cgto_cache) :: cache
   real(wp), parameter :: vec(3) = 0.0_wp, r2 = 0.0_wp
   real(wp) :: overlap(9, 9)

   call slater_to_gauss(ng, n, 4, 1.0_wp, cgto, .true., stat)
   call check(error, stat)
   if (allocated(error)) return

   call cgto%update(cache, .false.)

   call overlap_cgto(cgto, cgto, cache, cache, r2, vec, 100.0_wp, overlap)

   lp: do i = 1, 9
      do j = 1, 9
         call check(error, overlap(j, i), merge(1.0_wp, 0.0_wp, j == i), thr=thr)
         if (allocated(error)) exit lp
      end do
   end do lp
   if (allocated(error)) then
      print '(9es20.14e1)', overlap
   end if

end subroutine test_norm_g


subroutine test_norm_5g_sto1g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_f(error, 5, 1)
end subroutine test_norm_5g_sto1g


subroutine test_norm_5g_sto2g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_g(error, 5, 2)
end subroutine test_norm_5g_sto2g


subroutine test_norm_5g_sto3g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_g(error, 5, 3)
end subroutine test_norm_5g_sto3g


subroutine test_norm_5g_sto4g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_g(error, 5, 4)
end subroutine test_norm_5g_sto4g


subroutine test_norm_5g_sto5g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_g(error, 5, 5)
end subroutine test_norm_5g_sto5g


subroutine test_norm_5g_sto6g(error)
   type(error_type), allocatable, intent(out) :: error
   call test_norm_g(error, 5, 6)
end subroutine test_norm_5g_sto6g


end module test_slater_expansion
