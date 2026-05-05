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

module test_tagged_io
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use tblite_io_tag, only : tagged_data, tagged_entry, write_tagged
   implicit none
   private

   public :: collect_tagged_io

   real(wp), parameter :: thr = 10*epsilon(1.0_wp)

contains


!> Collect all exported unit tests
subroutine collect_tagged_io(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("write-read", test_write_read) &
      ]

end subroutine collect_tagged_io


subroutine test_write_read(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: io, stat

   type(tagged_data) :: tdata
   real(wp) :: gradient1(3, 8), sigma1(3, 3), energy1, charges1(8)
   !real(wp) :: gradient2(3, 8), sigma2(3, 3), energy2, charges2(8)

   call random_number(energy1)
   call random_number(gradient1)
   call random_number(sigma1)
   call random_number(charges1)

   open(newunit=io, status="scratch")

   call write_tagged(io, "energy", energy1, stat)
   call check(error, stat)
   if (allocated(error)) return

   call write_tagged(io, "gradient", gradient1, stat)
   call check(error, stat)
   if (allocated(error)) return

   call write_tagged(io, "sigma", sigma1, stat)
   call check(error, stat)
   if (allocated(error)) return

   call write_tagged(io, "charges", charges1, stat)
   call check(error, stat)
   if (allocated(error)) return

   rewind io

   call tdata%load(io, stat)
   call check(error, stat)
   if (allocated(error)) return

   call check(error, size(tdata%val) == 4)
   if (allocated(error)) return

   close(io)

end subroutine test_write_read


end module test_tagged_io
