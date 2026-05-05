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

module test_npy
   use mctc_env, only : dp, error_type, i4
   use mctc_env_testing, only : new_unittest, unittest_type, check, &
      & test_failed
   use tblite_io_numpy, only : save_npy, load_npy
   use tblite_output_format, only : format_string
   implicit none
   private

   public :: collect_npy

contains

subroutine collect_npy(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("write-npy-i4-r1", test_write_npy_i4_rank1), &
      new_unittest("write-npy-rdp-r1", test_write_npy_rdp_rank1), &
      new_unittest("write-npy-rdp-r2", test_write_npy_rdp_rank2), &
      new_unittest("write-npy-rdp-r3", test_write_npy_rdp_rank3), &
      new_unittest("read-npy-i4-r1", test_read_npy_i4_rank1), &
      new_unittest("read-npy-rdp-r1", test_read_npy_rdp_rank1), &
      new_unittest("read-npy-rdp-r2", test_read_npy_rdp_rank2), &
      new_unittest("read-npy-rdp-r3", test_read_npy_rdp_rank3), &
      new_unittest("missing-npy-i4-r1", test_missing_npy_i4_rank1, should_fail=.true.), &
      new_unittest("missing-npy-rdp-r1", test_missing_npy_rdp_rank1, should_fail=.true.), &
      new_unittest("missing-npy-rdp-r2", test_missing_npy_rdp_rank2, should_fail=.true.), &
      new_unittest("missing-npy-rdp-r3", test_missing_npy_rdp_rank3, should_fail=.true.), &
      new_unittest("type-mismatch", test_type_mismatch), &
      new_unittest("rank-mismatch-rdp", test_rank_mismatch_rdp), &
      new_unittest("rank-mismatch-i4", test_rank_mismatch_i4), &
      new_unittest("invalid-major-version", test_invalid_major_version), &
      new_unittest("invalid-minor-version", test_invalid_minor_version), &
      new_unittest("invalid-header-len", test_invalid_header_len), &
      new_unittest("invalid-nul-byte", test_invalid_nul_byte), &
      new_unittest("invalid-key", test_invalid_key), &
      new_unittest("invalid-comma", test_invalid_comma), &
      new_unittest("invalid-string", test_invalid_string), &
      new_unittest("invalid-tuple-literal", test_invalid_tuple_literal), &
      new_unittest("invalid-tuple-comma", test_invalid_tuple_comma), &
      new_unittest("invalid-tuple", test_invalid_tuple), &
      new_unittest("invalid-character", test_invalid_character), &
      new_unittest("invalid-true", test_invalid_true), &
      new_unittest("invalid-false", test_invalid_false), &
      new_unittest("invalid-magic-number", test_invalid_magic_number), &
      new_unittest("invalid-magic-string", test_invalid_magic_string), &
      new_unittest("duplicate-descr", test_duplicate_descr), &
      new_unittest("missing-descr", test_missing_descr), &
      new_unittest("missing-fortran_order", test_missing_fortran_order), &
      new_unittest("missing-shape", test_missing_shape), &
      new_unittest("iomsg-deallocated", test_iomsg_deallocated) &
   ]
end subroutine collect_npy

subroutine test_write_npy_i4_rank1(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: stat, i
   character(len=*), parameter :: filename = ".write-i4-r1.npy"
   integer(i4), allocatable :: input(:), output(:)

   allocate(input(64))
   input = [(i*(i-1), i=size(input), 1, -1)]
   call save_npy(filename, input, stat)

   call check(error, stat, "Writing of npy file failed")
   if (allocated(error)) return

   call load_npy(filename, output, stat)
   call delete_file(filename)

   call check(error, stat, "Reading of npy file failed")
   if (allocated(error)) return

   call check(error, size(output), size(input))
   if (allocated(error)) return

   call check(error, all((output - input) == 0), &
      "Precision loss when rereading array")
end subroutine test_write_npy_i4_rank1

subroutine test_write_npy_rdp_rank1(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: stat
   character(len=*), parameter :: filename = ".write-rdp-r1.npy"
   real(dp), allocatable :: input(:), output(:)

   allocate(input(65535))
   call random_number(input)
   call save_npy(filename, input, stat)

   call check(error, stat, "Writing of npy file failed")
   if (allocated(error)) return

   call load_npy(filename, output, stat)
   call delete_file(filename)

   call check(error, stat, "Reading of npy file failed")
   if (allocated(error)) return

   call check(error, size(output), size(input))
   if (allocated(error)) return

   call check(error, all(abs(output - input) <= epsilon(1.0_dp)), &
      "Precision loss when rereading array")
end subroutine test_write_npy_rdp_rank1

subroutine test_write_npy_rdp_rank2(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: stat
   character(len=*), parameter :: filename = ".write-rdp-r2.npy"
   real(dp), allocatable :: input(:, :), output(:, :)

   allocate(input(10, 4))
   call random_number(input)
   call save_npy(filename, input, stat)

   call check(error, stat, "Writing of npy file failed")
   if (allocated(error)) return

   call load_npy(filename, output, stat)
   call delete_file(filename)

   call check(error, stat, "Reading of npy file failed")
   if (allocated(error)) return

   call check(error, size(output), size(input))
   if (allocated(error)) return

   call check(error, all(abs(output - input) <= epsilon(1.0_dp)), &
      "Precision loss when rereading array")
end subroutine test_write_npy_rdp_rank2

subroutine test_write_npy_rdp_rank3(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: stat
   character(len=*), parameter :: filename = ".write-rdp-r3.npy"
   real(dp), allocatable :: input(:, :, :), output(:, :, :)

   allocate(input(10, 4, 5))
   call random_number(input)
   call save_npy(filename, input, stat)

   call check(error, stat, "Writing of npy file failed")
   if (allocated(error)) return

   call load_npy(filename, output, stat)
   call delete_file(filename)

   call check(error, stat, "Reading of npy file failed")
   if (allocated(error)) return

   call check(error, size(output), size(input))
   if (allocated(error)) return

   call check(error, all(abs(output - input) <= epsilon(1.0_dp)), &
      "Precision loss when rereading array")
end subroutine test_write_npy_rdp_rank3

subroutine test_read_npy_i4_rank1(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: dict = &
       "{'descr':'<i4', 'shape' : (88,   ) ,'fortran_order' :False }         " //  &
       char(10)
   character(len=*), parameter :: header = &
       char(int(z"93")) // "NUMPY" // char(1) // char(0) // &
       char(len(dict)) // char(0) // dict

   integer :: io, stat, i
   character(len=*), parameter :: filename = ".read-i4-r1.npy"
   integer(i4), allocatable :: input(:), output(:)
   character(len=:), allocatable :: msg

   allocate(input(88))
   input = [(i*(i+1), i=size(input), 1, -1)]

   open(newunit=io, file=filename, form="unformatted", access="stream")
   write(io) header
   write(io) input
   close(io)

   call load_npy(filename, output, stat, msg)
   call delete_file(filename)

   call check(error, stat, "Reading of npy file failed", msg)
   if (allocated(error)) return

   call check(error, size(output), size(input))
   if (allocated(error)) return

   call check(error, all((output - input) == 0), &
      "Precision loss when rereading array")
end subroutine test_read_npy_i4_rank1

subroutine test_read_npy_rdp_rank1(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: dict = &
       "{'descr':'<f8', 'shape' : (37 ,  ) , 'fortran_order': False}         " //  &
       char(10)
   character(len=*), parameter :: header = &
       char(int(z"93")) // "NUMPY" // char(1) // char(0) // &
       char(len(dict)) // char(0) // dict

   integer :: io, stat
   character(len=*), parameter :: filename = ".read-rdp-r1.npy"
   real(dp), allocatable :: input(:), output(:)

   allocate(input(37))
   call random_number(input)

   open(newunit=io, file=filename, form="unformatted", access="stream")
   write(io) header
   write(io) input
   close(io)

   call load_npy(filename, output, stat)
   call delete_file(filename)

   call check(error, stat, "Reading of npy file failed")
   if (allocated(error)) return

   call check(error, size(output), 37)
   if (allocated(error)) return

   call check(error, all(abs(output - input) <= epsilon(1.0_dp)), &
      "Precision loss when rereading array")
end subroutine test_read_npy_rdp_rank1

subroutine test_read_npy_rdp_rank2(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: dict = &
       "{'descr': '<f8', 'fortran_order': True, 'shape': (10, 4, ), }        " //  &
       char(10)
   character(len=*), parameter :: header = &
       char(int(z"93")) // "NUMPY" // char(1) // char(0) // &
       char(len(dict)) // char(0) // dict

   integer :: io, stat
   character(len=*), parameter :: filename = ".read-rdp-r2.npy"
   real(dp), allocatable :: input(:, :), output(:, :)

   allocate(input(10, 4))
   call random_number(input)

   open(newunit=io, file=filename, form="unformatted", access="stream")
   write(io) header
   write(io) input
   close(io)

   call load_npy(filename, output, stat)
   call delete_file(filename)

   call check(error, stat, "Reading of npy file failed")
   if (allocated(error)) return

   call check(error, size(output), size(input))
   if (allocated(error)) return

   call check(error, all(abs(output - input) <= epsilon(1.0_dp)), &
      "Precision loss when rereading array")
end subroutine test_read_npy_rdp_rank2

subroutine test_read_npy_rdp_rank3(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: dict = &
       "{'descr': '<f8', 'fortran_order': True, 'shape': (10, 2, 2), }       " //  &
       char(10)
   character(len=*), parameter :: header = &
       char(int(z"93")) // "NUMPY" // char(2) // char(0) // &
       char(len(dict)) // char(0) // char(0) // char(0) // dict

   integer :: io, stat
   character(len=*), parameter :: filename = ".read-rdp-r3.npy"
   real(dp), allocatable :: input(:, :, :), output(:, :, :)

   allocate(input(10, 2, 2))
   call random_number(input)

   open(newunit=io, file=filename, form="unformatted", access="stream")
   write(io) header
   write(io) input
   close(io)

   call load_npy(filename, output, stat)
   call delete_file(filename)

   call check(error, stat, "Reading of npy file failed")
   if (allocated(error)) return

   call check(error, size(output), size(input))
   if (allocated(error)) return

   call check(error, all(abs(output - input) <= epsilon(1.0_dp)), &
      "Precision loss when rereading array")
end subroutine test_read_npy_rdp_rank3

subroutine test_missing_npy_i4_rank1(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer(i4), allocatable :: output(:)
   integer :: stat
   character(len=:), allocatable :: msg

   call load_npy(".missing-i4-r1.npy", output, stat, msg)
   call check(error, stat, "Loading npy file failed", msg)
end subroutine test_missing_npy_i4_rank1

subroutine test_missing_npy_rdp_rank1(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(dp), allocatable :: output(:)
   integer :: stat
   character(len=:), allocatable :: msg

   call load_npy(".missing-rdp-r1.npy", output, stat, msg)
   call check(error, stat, "Loading npy file failed", msg)
end subroutine test_missing_npy_rdp_rank1

subroutine test_missing_npy_rdp_rank2(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(dp), allocatable :: output(:, :)
   integer :: stat
   character(len=:), allocatable :: msg

   call load_npy(".missing-rdp-r2.npy", output, stat, msg)
   call check(error, stat, "Loading npy file failed", msg)
end subroutine test_missing_npy_rdp_rank2

subroutine test_missing_npy_rdp_rank3(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(dp), allocatable :: output(:, :, :)
   integer :: stat
   character(len=:), allocatable :: msg

   call load_npy(".missing-rdp-r3.npy", output, stat, msg)
   call check(error, stat, "Loading npy file failed", msg)
end subroutine test_missing_npy_rdp_rank3

subroutine test_invalid_character(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: dict = &
       "{'descr': '<f8', %'fortran_order': True, 'shape': (10, 4, ), }       " //  &
       char(10)
   character(len=*), parameter :: header = &
       char(int(z"93")) // "NUMPY" // char(1) // char(0) // &
       char(len(dict)) // char(0) // dict

   integer :: io, stat
   character(len=*), parameter :: filename = ".test-invalid-character.npy"
   real(dp), allocatable :: array(:, :)

   open(newunit=io, file=filename, form="unformatted", access="stream")
   write(io) header
   write(io) spread(0.0_dp, 1, 40)
   close(io)

   call load_npy(filename, array, stat)
   call delete_file(filename)

   call check(error, stat, 305, "Expected error code 305 got "//format_string(stat, '(i0)'))
end subroutine test_invalid_character

subroutine test_invalid_true(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: dict = &
       "{'descr': '<f8',  'fortran_order': Tru , 'shape': (10, 4, ), }       " //  &
       char(10)
   character(len=*), parameter :: header = &
       char(int(z"93")) // "NUMPY" // char(1) // char(0) // &
       char(len(dict)) // char(0) // dict

   integer :: io, stat
   character(len=*), parameter :: filename = ".test-invalid-character.npy"
   real(dp), allocatable :: array(:, :)

   open(newunit=io, file=filename, form="unformatted", access="stream")
   write(io) header
   write(io) spread(0.0_dp, 1, 40)
   close(io)

   call load_npy(filename, array, stat)
   call delete_file(filename)

   call check(error, stat, 321, "Expected error code 305 got "//format_string(stat, '(i0)'))
end subroutine test_invalid_true

subroutine test_invalid_false(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: dict = &
       "{'descr': '<f8',  'fortran_order': Fls , 'shape': (10, 4, ), }       " //  &
       char(10)
   character(len=*), parameter :: header = &
       char(int(z"93")) // "NUMPY" // char(1) // char(0) // &
       char(len(dict)) // char(0) // dict

   integer :: io, stat
   character(len=*), parameter :: filename = ".test-invalid-character.npy"
   real(dp), allocatable :: array(:, :)

   open(newunit=io, file=filename, form="unformatted", access="stream")
   write(io) header
   write(io) spread(0.0_dp, 1, 40)
   close(io)

   call load_npy(filename, array, stat)
   call delete_file(filename)

   call check(error, stat, 321, "Expected error code 305 got "//format_string(stat, '(i0)'))
end subroutine test_invalid_false

subroutine test_invalid_magic_number(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: dict = &
       "{'descr': '<f8', 'fortran_order': True, 'shape': (10, 4, ), }        " //  &
       char(10)
   character(len=*), parameter :: header = &
       char(50) // "NUMPY" // char(1) // char(0) // &
       char(len(dict)) // char(0) // dict

   integer :: io, stat
   character(len=*), parameter :: filename = ".test-invalid-magic-num.npy"
   real(dp), allocatable :: array(:, :)

   open(newunit=io, file=filename, form="unformatted", access="stream")
   write(io) header
   write(io) spread(0.0_dp, 1, 40)
   close(io)

   call load_npy(filename, array, stat)
   call delete_file(filename)

   call check(error, stat, 201, "Expected error code 201 got "//format_string(stat, '(i0)'))
end subroutine test_invalid_magic_number

subroutine test_invalid_magic_string(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: dict = &
       "{'descr': '<f8', 'fortran_order': True, 'shape': (10, 4, ), }        " //  &
       char(10)
   character(len=*), parameter :: header = &
       char(int(z"93")) // "numpy" // char(1) // char(0) // &
       char(len(dict)) // char(0) // dict

   integer :: io, stat
   character(len=*), parameter :: filename = ".test-invalid-magic-str.npy"
   real(dp), allocatable :: array(:, :)

   open(newunit=io, file=filename, form="unformatted", access="stream")
   write(io) header
   write(io) spread(0.0_dp, 1, 40)
   close(io)

   call load_npy(filename, array, stat)
   call delete_file(filename)

   call check(error, stat, 202, "Expected error code 202 got "//format_string(stat, '(i0)'))
end subroutine test_invalid_magic_string

subroutine test_invalid_major_version(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: dict = &
       "{'descr': '<f8', 'fortran_order': True, 'shape': (10, 4, ), }        " //  &
       char(10)
   character(len=*), parameter :: header = &
       char(int(z"93")) // "NUMPY" // char(0) // char(0) // &
       char(len(dict)) // char(0) // dict

   integer :: io, stat
   character(len=*), parameter :: filename = ".test-invalid-major-version.npy"
   real(dp), allocatable :: array(:, :)

   open(newunit=io, file=filename, form="unformatted", access="stream")
   write(io) header
   write(io) spread(0.0_dp, 1, 40)
   close(io)

   call load_npy(filename, array, stat)
   call delete_file(filename)

   call check(error, stat, 203,"Expected error code 203 got "//format_string(stat, '(i0)'))
end subroutine test_invalid_major_version

subroutine test_invalid_minor_version(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: dict = &
       "{'descr': '<f8', 'fortran_order': True, 'shape': (10, 4, ), }        " //  &
       char(10)
   character(len=*), parameter :: header = &
       char(int(z"93")) // "NUMPY" // char(1) // char(9) // &
       char(len(dict)) // char(0) // dict

   integer :: io, stat
   character(len=*), parameter :: filename = ".test-invalid-minor-version.npy"
   real(dp), allocatable :: array(:, :)

   open(newunit=io, file=filename, form="unformatted", access="stream")
   write(io) header
   write(io) spread(0.0_dp, 1, 40)
   close(io)

   call load_npy(filename, array, stat)
   call delete_file(filename)

   call check(error, stat, 204, "Expected error code 204 got "//format_string(stat, '(i0)'))
end subroutine test_invalid_minor_version

subroutine test_invalid_header_len(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: dict = &
       "{'descr': '<f8', 'fortran_order': True, 'shape': (10, 4, ), }        " //  &
       char(10)
   character(len=*), parameter :: header = &
       char(int(z"93")) // "NUMPY" // char(1) // char(0) // &
       char(len(dict)-1) // char(0) // dict

   integer :: io, stat
   character(len=*), parameter :: filename = ".test-invalid-header-len.npy"
   real(dp), allocatable :: array(:, :)

   open(newunit=io, file=filename, form="unformatted", access="stream")
   write(io) header
   write(io) spread(0.0_dp, 1, 40)
   close(io)

   call load_npy(filename, array, stat)
   call delete_file(filename)

   call check(error, stat, 101, "Expected error code 101 got "//format_string(stat, '(i0)'))
end subroutine test_invalid_header_len

subroutine test_invalid_nul_byte(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: dict = &
       "{'descr': '<f8', 'fortran_order': True, 'shape': (10, 4, ), }       " //  &
       char(0) // char(10)
   character(len=*), parameter :: header = &
       char(int(z"93")) // "NUMPY" // char(1) // char(0) // &
       char(len(dict)) // char(0) // dict

   integer :: io, stat
   character(len=*), parameter :: filename = ".test-invalid-nul-byte.npy"
   real(dp), allocatable :: array(:, :)

   open(newunit=io, file=filename, form="unformatted", access="stream")
   write(io) header
   write(io) spread(0.0_dp, 1, 40)
   close(io)

   call load_npy(filename, array, stat)
   call delete_file(filename)

   call check(error, stat, 102, "Expected error code 102 got "//format_string(stat, '(i0)'))
end subroutine test_invalid_nul_byte

subroutine test_invalid_key(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: dict = &
       "{'fortran_order': True, 'shape': (10, 4, ), 'descr': '<f8', 'x': 1, }" //  &
       char(10)
   character(len=*), parameter :: header = &
       char(int(z"93")) // "NUMPY" // char(1) // char(0) // &
       char(len(dict)) // char(0) // dict

   integer :: io, stat
   character(len=*), parameter :: filename = ".test-invalid-key.npy"
   real(dp), allocatable :: array(:, :)

   open(newunit=io, file=filename, form="unformatted", access="stream")
   write(io) header
   write(io) spread(0.0_dp, 1, 40)
   close(io)

   call load_npy(filename, array, stat)
   call delete_file(filename)

   call check(error, stat, 304, "Expected error code 304 got "//format_string(stat, '(i0)'))
end subroutine test_invalid_key

subroutine test_invalid_comma(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: dict = &
       "{'fortran_order': True,, 'shape': (10, 4, ), 'descr': '<f8', }       " //  &
       "   " //  &
       char(10)
   character(len=*), parameter :: header = &
       char(int(z"93")) // "NUMPY" // char(1) // char(0) // &
       char(len(dict)) // char(0) // dict

   integer :: io, stat
   character(len=*), parameter :: filename = ".test-invalid-comma.npy"
   real(dp), allocatable :: array(:)

   open(newunit=io, file=filename, form="unformatted", access="stream")
   write(io) header
   write(io) spread(0.0_dp, 1, 40)
   close(io)

   call load_npy(filename, array, stat)
   call delete_file(filename)

   call check(error, stat, 301, "Expected error code 301 got "//format_string(stat, '(i0)'))
end subroutine test_invalid_comma

subroutine test_invalid_string(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: dict = &
       "{'fortran_order': True, 'shape': (10, 4, ), 'descr': '<f8' '<f4', }  " //  &
       "   " //  &
       char(10)
   character(len=*), parameter :: header = &
       char(int(z"93")) // "NUMPY" // char(1) // char(0) // &
       char(len(dict)) // char(0) // dict

   integer :: io, stat
   character(len=:), allocatable :: msg
   character(len=*), parameter :: filename = ".test-invalid-string.npy"
   real(dp), allocatable :: array(:)

   open(newunit=io, file=filename, form="unformatted", access="stream")
   write(io) header
   write(io) spread(0.0_dp, 1, 40)
   close(io)

   call load_npy(filename, array, stat, msg)
   call delete_file(filename)

   call check(error, stat, 302, "Expected error code 302 got "//format_string(stat, '(i0)'))
end subroutine test_invalid_string

subroutine test_invalid_tuple_literal(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: dict = &
       "{'fortran_order': True, 'shape': (10  4, ), 'descr': '<f8'        }  " //  &
       "   " //  &
       char(10)
   character(len=*), parameter :: header = &
       char(int(z"93")) // "NUMPY" // char(1) // char(0) // &
       char(len(dict)) // char(0) // dict

   integer :: io, stat
   character(len=:), allocatable :: msg
   character(len=*), parameter :: filename = ".test-invalid-tuple-literal.npy"
   real(dp), allocatable :: array(:)

   open(newunit=io, file=filename, form="unformatted", access="stream")
   write(io) header
   write(io) spread(0.0_dp, 1, 40)
   close(io)

   call load_npy(filename, array, stat, msg)
   call delete_file(filename)

   call check(error, stat, 311, "Expected error code 311 got "//format_string(stat, '(i0)'))
end subroutine test_invalid_tuple_literal

subroutine test_invalid_tuple_comma(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: dict = &
       "{'fortran_order': True, 'shape': (10,,4, ), 'descr': '<f8'        }  " //  &
       "   " //  &
       char(10)
   character(len=*), parameter :: header = &
       char(int(z"93")) // "NUMPY" // char(1) // char(0) // &
       char(len(dict)) // char(0) // dict

   integer :: io, stat
   character(len=:), allocatable :: msg
   character(len=*), parameter :: filename = ".test-invalid-tuple-comma.npy"
   real(dp), allocatable :: array(:)

   open(newunit=io, file=filename, form="unformatted", access="stream")
   write(io) header
   write(io) spread(0.0_dp, 1, 40)
   close(io)

   call load_npy(filename, array, stat, msg)
   call delete_file(filename)

   call check(error, stat, 312, "Expected error code 312 got "//format_string(stat, '(i0)'))
end subroutine test_invalid_tuple_comma

subroutine test_invalid_tuple(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: dict = &
       "{'fortran_order': True, 'shape': (10,%4, ), 'descr': '<f8'        }  " //  &
       "   " //  &
       char(10)
   character(len=*), parameter :: header = &
       char(int(z"93")) // "NUMPY" // char(1) // char(0) // &
       char(len(dict)) // char(0) // dict

   integer :: io, stat
   character(len=:), allocatable :: msg
   character(len=*), parameter :: filename = ".test-invalid-tuple.npy"
   real(dp), allocatable :: array(:)

   open(newunit=io, file=filename, form="unformatted", access="stream")
   write(io) header
   write(io) spread(0.0_dp, 1, 40)
   close(io)

   call load_npy(filename, array, stat, msg)
   call delete_file(filename)

   call check(error, stat, 313, "Expected error code 313 got "//format_string(stat, '(i0)'))
end subroutine test_invalid_tuple

subroutine test_type_mismatch(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: dict = &
       "{'descr': '<f4',                 'fortran_order': True, 'shape': (40, ), }"//&
       "   " //  &
       char(10)
   character(len=*), parameter :: header = &
       char(int(z"93")) // "NUMPY" // char(1) // char(0) // &
       char(len(dict)) // char(0) // dict

   integer :: io, stat
   character(len=:), allocatable :: msg
   character(len=*), parameter :: filename = ".test-type-mismatch.npy"
   integer, allocatable :: iarray1(:)
   real(dp), allocatable :: rarray1(:), rarray2(:, :), rarray3(:, :, :)

   open(newunit=io, file=filename, form="unformatted", access="stream")
   write(io) header
   write(io) spread(0.0_dp, 1, 40)
   close(io)

   call load_npy(filename, iarray1, stat, msg)
   call check(error, stat, 11, "Expected error code 11")

   if (.not.allocated(error)) then
      call load_npy(filename, rarray1, stat, msg)
      call check(error, stat, 11, "Expected error code 11 got "//format_string(stat, '(i0)'))
   end if

   if (.not.allocated(error)) then
      call load_npy(filename, rarray2, stat, msg)
      call check(error, stat, 11, "Expected error code 11 got "//format_string(stat, '(i0)'))
   end if

   if (.not.allocated(error)) then
      call load_npy(filename, rarray3, stat, msg)
      call check(error, stat, 11, "Expected error code 11 got "//format_string(stat, '(i0)'))
   end if
   call delete_file(filename)
end subroutine test_type_mismatch

subroutine test_rank_mismatch_i4(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: dict = &
       "{'descr': '<i4', 'fortran_order': True, 'shape': (40, 1, )               }"//&
       "   " //  &
       char(10)
   character(len=*), parameter :: header = &
       char(int(z"93")) // "NUMPY" // char(1) // char(0) // &
       char(len(dict)) // char(0) // dict

   integer :: io, stat
   character(len=:), allocatable :: msg
   character(len=*), parameter :: filename = ".test-rank-mismatch-rdp.npy"
   integer(i4), allocatable :: iarray1(:)

   open(newunit=io, file=filename, form="unformatted", access="stream")
   write(io) header
   write(io) spread(0.0_dp, 1, 40)
   close(io)

   call load_npy(filename, iarray1, stat, msg)
   call check(error, stat, 12, "Expected error code 12 got "//format_string(stat, '(i0)'))

   call delete_file(filename)
end subroutine test_rank_mismatch_i4

subroutine test_rank_mismatch_rdp(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: dict = &
       "{'descr': '<f8', 'fortran_order': True, 'shape': (40, 1, 1, 1, )         }"//&
       "   " //  &
       char(10)
   character(len=*), parameter :: header = &
       char(int(z"93")) // "NUMPY" // char(1) // char(0) // &
       char(len(dict)) // char(0) // dict

   integer :: io, stat
   character(len=:), allocatable :: msg
   character(len=*), parameter :: filename = ".test-rank-mismatch-rdp.npy"
   real(dp), allocatable :: rarray1(:), rarray2(:, :), rarray3(:, :, :)

   open(newunit=io, file=filename, form="unformatted", access="stream")
   write(io) header
   write(io) spread(0.0_dp, 1, 40)
   close(io)

   call load_npy(filename, rarray1, stat, msg)
   call check(error, stat, 12, "Expected error code 12 got "//format_string(stat, '(i0)'))

   if (.not.allocated(error)) then
      call load_npy(filename, rarray2, stat, msg)
      call check(error, stat, 12, "Expected error code 12 got "//format_string(stat, '(i0)'))
   end if

   if (.not.allocated(error)) then
      call load_npy(filename, rarray3, stat, msg)
      call check(error, stat, 12, "Expected error code 12 got "//format_string(stat, '(i0)'))
   end if

   call delete_file(filename)
end subroutine test_rank_mismatch_rdp

subroutine test_duplicate_descr(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: dict = &
       "{'descr': '<f8', 'descr': '<f8', 'fortran_order': True, 'shape': (40, ), }"//&
       "   " //  &
       char(10)
   character(len=*), parameter :: header = &
       char(int(z"93")) // "NUMPY" // char(1) // char(0) // &
       char(len(dict)) // char(0) // dict

   integer :: io, stat
   character(len=*), parameter :: filename = ".test-invalid-descr.npy"
   real(dp), allocatable :: array(:)

   open(newunit=io, file=filename, form="unformatted", access="stream")
   write(io) header
   write(io) spread(0.0_dp, 1, 40)
   close(io)

   call load_npy(filename, array, stat)
   call delete_file(filename)

   call check(error, stat, 303, "Expected error code 303 got "//format_string(stat, '(i0)'))
end subroutine test_duplicate_descr

subroutine test_missing_descr(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: dict = &
       "{'fortran_order': True, 'shape': (10, 4, ), }                        " //  &
       char(10)
   character(len=*), parameter :: header = &
       char(int(z"93")) // "NUMPY" // char(1) // char(0) // &
       char(len(dict)) // char(0) // dict

   integer :: io, stat
   character(len=*), parameter :: filename = ".test-missing-descr.npy"
   real(dp), allocatable :: array(:, :)

   open(newunit=io, file=filename, form="unformatted", access="stream")
   write(io) header
   write(io) spread(0.0_dp, 1, 40)
   close(io)

   call load_npy(filename, array, stat)
   call delete_file(filename)

   call check(error, stat, 306, "Expected error code 306 got "//format_string(stat, '(i0)'))
end subroutine test_missing_descr

subroutine test_missing_fortran_order(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: dict = &
       "{'descr': '<f8', 'shape': (10, 4, ), }                               " //  &
       char(10)
   character(len=*), parameter :: header = &
       char(int(z"93")) // "NUMPY" // char(1) // char(0) // &
       char(len(dict)) // char(0) // dict

   integer :: io, stat
   character(len=*), parameter :: filename = ".test-missing-fortran_order.npy"
   real(dp), allocatable :: array(:, :)

   open(newunit=io, file=filename, form="unformatted", access="stream")
   write(io) header
   write(io) spread(0.0_dp, 1, 40)
   close(io)

   call load_npy(filename, array, stat)
   call delete_file(filename)

   call check(error, stat, 308, "Expected error code 308 got "//format_string(stat, '(i0)'))
end subroutine test_missing_fortran_order

subroutine test_missing_shape(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: dict = &
       "{'fortran_order': True, 'descr': '<f8'}                              " //  &
       char(10)
   character(len=*), parameter :: header = &
       char(int(z"93")) // "NUMPY" // char(1) // char(0) // &
       char(len(dict)) // char(0) // dict

   integer :: io, stat
   character(len=*), parameter :: filename = ".test-missing-shape.npy"
   real(dp), allocatable :: array(:, :)

   open(newunit=io, file=filename, form="unformatted", access="stream")
   write(io) header
   write(io) spread(0.0_dp, 1, 40)
   close(io)

   call load_npy(filename, array, stat)
   call delete_file(filename)

   call check(error, stat, 307, "Expected error code 307 got "//format_string(stat, '(i0)'))
end subroutine test_missing_shape

subroutine test_iomsg_deallocated(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: stat
   character(len=:), allocatable :: msg

   character(len=*), parameter :: filename = ".test-iomsg-deallocated.npy"
   real(dp), allocatable :: input(:, :)

   msg = "This message should be deallocated."

   allocate(input(12, 5))
   call random_number(input)
   call save_npy(filename, input, stat, msg)
   call delete_file(filename)

   call check(error, .not.allocated(msg), "Message wrongly allocated.")

end subroutine test_iomsg_deallocated

subroutine delete_file(filename)
   character(len=*), intent(in) :: filename

   integer :: io

   open(newunit=io, file=filename)
   close(io, status="delete")
end subroutine delete_file

end module test_npy