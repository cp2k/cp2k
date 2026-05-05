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

module test_npz
   use mctc_env, only : dp, error_type, i2, i4
   use mctc_env_testing, only : new_unittest, unittest_type, check, &
      & test_failed
   use tblite_io_numpy, only : save_npz, load_npz
   use tblite_io_numpy_crc32, only : crc32_hash
   use tblite_output_format, only : format_string
   implicit none
   private

   public :: collect_npz

   ! zip file headers
   character(len=*), parameter :: local_header_start(*) = &
       &[char(int(z"50")), char(int(z"4b")), char(int(z"03")), char(int(z"04")), &
       & char(int(z"14")), char(0), char(0), char(0), char(0), char(0), &
       & char(0), char(0), char(0), char(0)]
   character(len=*), parameter :: global_header_start(*) = &
       &[char(int(z"50")), char(int(z"4b")), char(int(z"01")), char(int(z"02")), &
       & char(int(z"14")), char(0), local_header_start(5:)]
   character(len=*), parameter :: global_header_final(*) = &
      &[char(0), char(0), char(0), char(0), char(0), char(0), char(0), char(0), &
      & char(0), char(0), char(0), char(0), char(0), char(0)]
   character(len=*), parameter :: footer_start(*) = &
      &[char(int(z"50")), char(int(z"4b")), char(int(z"05")), char(int(z"06")), &
      & char(0), char(0), char(0), char(0), char(1), char(0), char(1), char(0)]

contains

subroutine collect_npz(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("write-npz-i4-r1", test_write_npz_i4_rank1), &
      new_unittest("write-npz-rdp-r1", test_write_npz_rdp_rank1), &
      new_unittest("write-npz-rdp-r2", test_write_npz_rdp_rank2), &
      new_unittest("write-npz-rdp-r3", test_write_npz_rdp_rank3), &
      new_unittest("write-npz-i4-rank1-multi", test_write_npz_i4_rank1_multi), &
      new_unittest("write-npz-rdp-rank1-multi", test_write_npz_rdp_rank1_multi), &
      new_unittest("write-npz-rdp-rank2-multi", test_write_npz_rdp_rank2_multi), &
      new_unittest("write-npz-rdp-rank3-multi", test_write_npz_rdp_rank3_multi), &
      new_unittest("read-npz-rdp-r1", test_read_npz_rdp_rank1), &
      new_unittest("read-npz-rdp-r2", test_read_npz_rdp_rank2), &
      new_unittest("read-npz-rdp-r3", test_read_npz_rdp_rank3), &
      new_unittest("read-npz-file-not-found", test_read_file_not_found, should_fail=.true.), &
      new_unittest("read-npz-checksum-mismatch", test_read_checksum_mismatch, should_fail=.true.), &
      new_unittest("write-npz-missing-global1", test_write_missing_global1, should_fail=.true.), &
      new_unittest("write-npz-missing-global2", test_write_missing_global2, should_fail=.true.), &
      new_unittest("write-npz-missing-footer1", test_write_missing_footer1, should_fail=.true.), &
      new_unittest("write-npz-missing-footer2", test_write_missing_footer2, should_fail=.true.), &
      new_unittest("write-npz-mismatch-header", test_write_mismatch_header, should_fail=.true.), &
      new_unittest("write-npz-mismatch-footer", test_write_mismatch_footer, should_fail=.true.), &
      new_unittest("crc32", test_crc32) &
   ]
end subroutine collect_npz

subroutine test_write_npz_i4_rank1(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: stat, i
   character(len=*), parameter :: filename = ".write-i4-r1.npz", varname = "test"
   integer(i4), allocatable :: input(:), output(:)
   character(len=:), allocatable :: msg

   allocate(input(72))
   input = [(2*i, i=1, size(input))]
   call save_npz(filename, varname, input, stat, msg)

   call check(error, stat, "Writing of npz file failed", msg)
   if (allocated(error)) return

   call load_npz(filename, varname, output, stat, msg)
   call delete_file(filename)

   call check(error, stat, "Reading of npz file failed", msg)
   if (allocated(error)) return

   call check(error, size(output), size(input))
   if (allocated(error)) return

   call check(error, all((output - input) == 0), &
      "Precision loss when rereading array")
end subroutine test_write_npz_i4_rank1

subroutine test_write_npz_rdp_rank1(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: stat
   character(len=*), parameter :: filename = ".write-rdp-r1.npz", varname = "test"
   real(dp), allocatable :: input(:), output(:)
   character(len=:), allocatable :: msg

   allocate(input(51))
   call random_number(input)
   call save_npz(filename, varname, input, stat)

   call check(error, stat, "Writing of npz file failed")
   if (allocated(error)) return

   call load_npz(filename, varname, output, stat, msg)
   call delete_file(filename)

   call check(error, stat, "Reading of npy file failed", msg)
   if (allocated(error)) return

   call check(error, size(output), size(input))
   if (allocated(error)) return

   call check(error, all(abs(output - input) <= epsilon(1.0_dp)), &
      "Precision loss when rereading array")
end subroutine test_write_npz_rdp_rank1

subroutine test_write_npz_rdp_rank2(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: stat
   character(len=*), parameter :: filename = ".write-rdp-r2.npz", varname = "test"
   real(dp), allocatable :: input(:, :), output(:, :)
   character(len=:), allocatable :: msg

   allocate(input(13, 7))
   call random_number(input)
   call save_npz(filename, varname, input, stat, msg)

   call check(error, stat, "Writing of npz file failed", msg)
   if (allocated(error)) return

   call load_npz(filename, varname, output, stat, msg)
   call delete_file(filename)

   call check(error, stat, "Reading of npy file failed", msg)
   if (allocated(error)) return

   call check(error, size(output), size(input))
   if (allocated(error)) return

   call check(error, all(abs(output - input) <= epsilon(1.0_dp)), &
      "Precision loss when rereading array")
end subroutine test_write_npz_rdp_rank2

subroutine test_write_npz_rdp_rank3(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: stat
   character(len=*), parameter :: filename = ".write-rdp-r1.npz", varname = "test"
   real(dp), allocatable :: input(:, :, :), output(:, :, :)
   character(len=:), allocatable :: msg

   allocate(input(9, 2, 11))
   call random_number(input)
   call save_npz(filename, varname, input, stat, msg)

   call check(error, stat, "Writing of npz file failed", msg)
   if (allocated(error)) return

   call load_npz(filename, varname, output, stat, msg)
   call delete_file(filename)

   call check(error, stat, "Reading of npy file failed", msg)
   if (allocated(error)) return

   call check(error, size(output), size(input))
   if (allocated(error)) return

   call check(error, all(abs(output - input) <= epsilon(1.0_dp)), &
      "Precision loss when rereading array")
end subroutine test_write_npz_rdp_rank3

subroutine test_write_npz_i4_rank1_multi(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: stat, i
   character(len=*), parameter :: filename = ".mwrite-i4-r1.npz"
   integer(i4), allocatable :: input(:), output1(:), output2(:)
   character(len=:), allocatable :: msg

   allocate(input(93))
   input = [(i*(i-1), i=size(input), 1, -1)]
   call save_npz(filename, "array1", input, stat, msg)

   if (stat == 0) &
      call save_npz(filename, "arr2", input, stat, msg)

   call check(error, stat, "Writing of npz file failed", msg)
   if (allocated(error)) return

   call load_npz(filename, "array1", output1, stat, msg)
   if (stat == 0) &
      call load_npz(filename, "arr2", output2, stat, msg)
   call delete_file(filename)

   call check(error, stat, "Reading of npz file failed", msg)
   if (allocated(error)) return

   call check(error, size(output1), size(input))
   if (allocated(error)) return

   call check(error, size(output2), size(input))
   if (allocated(error)) return

   call check(error, all((output1 - input) == 0), &
      "Precision loss when rereading array")

   call check(error, all((output2 - input) == 0), &
      "Precision loss when rereading array")
end subroutine test_write_npz_i4_rank1_multi

subroutine test_write_npz_rdp_rank1_multi(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: stat, i
   character(len=*), parameter :: filename = ".mwrite-rdp-r1.npz"
   real(dp), allocatable :: input(:), output1(:), output2(:)
   character(len=:), allocatable :: msg

   allocate(input(66))
   call random_number(input)
   call save_npz(filename, "array1", input, stat, msg)

   if (stat == 0) &
      call save_npz(filename, "arr2", input, stat, msg)

   call check(error, stat, "Writing of npz file failed", msg)
   if (allocated(error)) return

   call load_npz(filename, "array1", output1, stat, msg)
   if (stat == 0) &
      call load_npz(filename, "arr2", output2, stat, msg)
   call delete_file(filename)

   call check(error, stat, "Reading of npz file failed", msg)
   if (allocated(error)) return

   call check(error, size(output1), size(input))
   if (allocated(error)) return

   call check(error, size(output2), size(input))
   if (allocated(error)) return

   call check(error, all(abs(output1 - input) <= epsilon(1.0_dp)), &
      "Precision loss when rereading array")

   call check(error, all(abs(output2 - input) <= epsilon(1.0_dp)), &
      "Precision loss when rereading array")
end subroutine test_write_npz_rdp_rank1_multi

subroutine test_write_npz_rdp_rank2_multi(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: stat, i
   character(len=*), parameter :: filename = ".mwrite-rdp-r2.npz"
   real(dp), allocatable :: input(:, :), output1(:, :), output2(:, :)
   character(len=:), allocatable :: msg

   allocate(input(19, 23))
   call random_number(input)
   call save_npz(filename, "array1", input, stat, msg)

   if (stat == 0) &
      call save_npz(filename, "arr2", input, stat, msg)

   call check(error, stat, "Writing of npz file failed", msg)
   if (allocated(error)) return

   call load_npz(filename, "array1", output1, stat, msg)
   if (stat == 0) &
      call load_npz(filename, "arr2", output2, stat, msg)
   call delete_file(filename)

   call check(error, stat, "Reading of npz file failed", msg)
   if (allocated(error)) return

   call check(error, size(output1), size(input))
   if (allocated(error)) return

   call check(error, size(output2), size(input))
   if (allocated(error)) return

   call check(error, all(abs(output1 - input) <= epsilon(1.0_dp)), &
      "Precision loss when rereading array")

   call check(error, all(abs(output2 - input) <= epsilon(1.0_dp)), &
      "Precision loss when rereading array")
end subroutine test_write_npz_rdp_rank2_multi

subroutine test_write_npz_rdp_rank3_multi(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: stat, i
   character(len=*), parameter :: filename = ".mwrite-rdp-r3.npz"
   real(dp), allocatable :: input(:, :, :), output1(:, :, :), output2(:, :, :)
   character(len=:), allocatable :: msg

   allocate(input(9, 10, 3))
   call random_number(input)
   call save_npz(filename, "array1", input, stat, msg)

   if (stat == 0) &
      call save_npz(filename, "arr2", input, stat, msg)

   call check(error, stat, "Writing of npz file failed", msg)
   if (allocated(error)) return

   call load_npz(filename, "array1", output1, stat, msg)
   if (stat == 0) &
      call load_npz(filename, "arr2", output2, stat, msg)
   call delete_file(filename)

   call check(error, stat, "Reading of npz file failed", msg)
   if (allocated(error)) return

   call check(error, size(output1), size(input))
   if (allocated(error)) return

   call check(error, size(output2), size(input))
   if (allocated(error)) return

   call check(error, all(abs(output1 - input) <= epsilon(1.0_dp)), &
      "Precision loss when rereading array")

   call check(error, all(abs(output2 - input) <= epsilon(1.0_dp)), &
      "Precision loss when rereading array")
end subroutine test_write_npz_rdp_rank3_multi

subroutine test_read_npz_rdp_rank1(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: dict = &
       "{'descr':'<f8', 'shape' : ( 61 ,  ) , 'fortran_order': False}         " //  &
       char(10)
   character(len=*), parameter :: header = &
       char(int(z"93")) // "NUMPY" // char(1) // char(0) // &
       char(len(dict)) // char(0) // dict

   integer, parameter :: vsize(*) = [61]

   character(len=*), parameter :: varname = "test", pathname = varname // ".npy"
   integer(i2), parameter :: path_size = len(pathname)
   integer(i4), parameter :: nbytes = len(header) + product(vsize) * 8
   integer(i4), parameter :: global_header_size = 46 + path_size
   integer(i4), parameter :: global_header_offset = 30 + path_size + nbytes

   integer :: io, stat
   character(len=*), parameter :: filename = ".read-rdp-r1.npz"
   integer(i4) :: crc_expected
   real(dp), allocatable :: input(:), output(:)
   character(len=:), allocatable :: msg

   allocate(input(vsize(1)))
   call random_number(input)

   crc_expected = crc32_hash(header)
   crc_expected = crc32_hash(input, crc_expected)

   open(newunit=io, file=filename, form="unformatted", access="stream")
   write(io) local_header_start
   write(io) crc_expected, nbytes, nbytes, path_size, 0_i2
   write(io) pathname
   write(io) header
   write(io) input
   write(io) global_header_start
   write(io) crc_expected, nbytes, nbytes, path_size, 0_i2
   write(io) global_header_final
   write(io) pathname
   write(io) footer_start
   write(io) global_header_size , global_header_offset, 0_i2
   close(io)

   call load_npz(filename, varname, output, stat, msg)
   call delete_file(filename)

   call check(error, stat, "Reading of npz file failed", msg)
   if (allocated(error)) return

   call check(error, size(output), product(vsize))
   if (allocated(error)) return

   call check(error, all(abs(output - input) <= epsilon(1.0_dp)), &
      "Precision loss when rereading array")
end subroutine test_read_npz_rdp_rank1

subroutine test_read_npz_rdp_rank2(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: dict = &
       "{'shape' : (21 , 5 ) , 'fortran_order':True,  'descr':'<f8', }    " //  &
       char(10)
   character(len=*), parameter :: header = &
       char(int(z"93")) // "NUMPY" // char(1) // char(0) // &
       char(len(dict)) // char(0) // dict

   integer, parameter :: vsize(*) = [21, 5]

   character(len=*), parameter :: varname = "test", pathname = varname // ".npy"
   integer(i2), parameter :: path_size = len(pathname)
   integer(i4), parameter :: nbytes = len(header) + product(vsize) * 8
   integer(i4), parameter :: global_header_size = 46 + path_size
   integer(i4), parameter :: global_header_offset = 30 + path_size + nbytes

   integer :: io, stat
   character(len=*), parameter :: filename = ".read-rdp-r2.npz"
   integer(i4) :: crc_expected
   real(dp), allocatable :: input(:, :), output(:, :)
   character(len=:), allocatable :: msg

   allocate(input(vsize(1), vsize(2)))
   call random_number(input)

   crc_expected = crc32_hash(header)
   crc_expected = crc32_hash(reshape(input, [size(input)]), crc_expected)

   open(newunit=io, file=filename, form="unformatted", access="stream")
   write(io) local_header_start
   write(io) crc_expected, nbytes, nbytes, path_size, 0_i2
   write(io) pathname
   write(io) header
   write(io) input
   write(io) global_header_start
   write(io) crc_expected, nbytes, nbytes, path_size, 0_i2
   write(io) global_header_final
   write(io) pathname
   write(io) footer_start
   write(io) global_header_size , global_header_offset, 0_i2
   close(io)

   call load_npz(filename, varname, output, stat, msg)
   call delete_file(filename)

   call check(error, stat, "Reading of npz file failed", msg)
   if (allocated(error)) return

   call check(error, size(output), product(vsize))
   if (allocated(error)) return

   call check(error, all(abs(output - input) <= epsilon(1.0_dp)), &
      "Precision loss when rereading array")
end subroutine test_read_npz_rdp_rank2

subroutine test_read_npz_rdp_rank3(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: dict = &
       "{'descr': '<f8', 'fortran_order': True, 'shape': (3,7,4), }" // &
       repeat(char(int(z'20')), 59)// char(10)
   character(len=*), parameter :: header = &
       char(int(z"93")) // "NUMPY" // char(1) // char(0) // &
       char(len(dict)) // char(0) // dict

   integer, parameter :: vsize(*) = [3, 7, 4]

   character(len=*), parameter :: varname = "test", pathname = varname // ".npy"
   integer(i2), parameter :: path_size = len(pathname)
   integer(i4), parameter :: nbytes = len(header) + product(vsize) * 8
   integer(i4), parameter :: global_header_size = 46 + path_size
   integer(i4), parameter :: global_header_offset = 30 + path_size + nbytes

   integer :: io, stat
   character(len=*), parameter :: filename = ".read-rdp-r3.npz"
   integer(i4) :: crc_expected
   real(dp), allocatable :: input(:, :, :), output(:, :, :)
   character(len=:), allocatable :: msg

   allocate(input(vsize(1), vsize(2), vsize(3)))
   call random_number(input)

   crc_expected = crc32_hash(header)
   crc_expected = crc32_hash(reshape(input, [size(input)]), crc_expected)

   open(newunit=io, file=filename, form="unformatted", access="stream")
   write(io) local_header_start
   write(io) crc_expected, nbytes, nbytes, path_size, 0_i2
   write(io) pathname
   write(io) header
   write(io) input
   write(io) global_header_start
   write(io) crc_expected, nbytes, nbytes, path_size, 0_i2
   write(io) global_header_final
   write(io) pathname
   write(io) footer_start
   write(io) global_header_size , global_header_offset, 0_i2
   close(io)

   call load_npz(filename, varname, output, stat, msg)
   call delete_file(filename)

   call check(error, stat, "Reading of npz file failed", msg)
   if (allocated(error)) return

   call check(error, size(output), product(vsize))
   if (allocated(error)) return

   call check(error, all(abs(output - input) <= epsilon(1.0_dp)), &
      "Precision loss when rereading array")
end subroutine test_read_npz_rdp_rank3

subroutine test_crc32(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call check_hash(error, crc32_hash("a", 0), int(z'e8b7be43'))
   if (allocated(error)) return
   call check_hash(error, crc32_hash(["a"], 0), int(z'e8b7be43'))
   if (allocated(error)) return
   call check_hash(error, crc32_hash([0.0_dp], 0), int(z'6522df69'))
   if (allocated(error)) return
   call check_hash(error, crc32_hash([1.1_dp], 0), int(z'ce0696c5'))
   if (allocated(error)) return
   call check_hash(error, crc32_hash([1.1_dp, 1.4_dp, -9.7_dp, 1.8e-23_dp, 3.5e12_dp], 0), int(z'e36ce8e0'))

contains
   subroutine check_hash(error, computed, expected)
      type(error_type), allocatable, intent(out) :: error
      integer(i4), intent(in) :: computed
      integer(i4), intent(in) :: expected

      if (computed /= expected) then
         call test_failed(error, "CRC32 hash mismatch, expected "//format_string(expected, '(z08)')//&
            & " got "//format_string(computed, '(z08)'))
      end if
   end subroutine check_hash
end subroutine test_crc32

subroutine test_read_file_not_found(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: dict = &
       "{'descr':'<f8', 'shape' : ( 61 ,  ) , 'fortran_order': False}         " //  &
       char(10)
   character(len=*), parameter :: header = &
       char(int(z"93")) // "NUMPY" // char(1) // char(0) // &
       char(len(dict)) // char(0) // dict

   integer, parameter :: vsize(*) = [61]

   character(len=*), parameter :: varname = "test", pathname = varname // ".npy"
   integer(i2), parameter :: path_size = len(pathname)
   integer(i4), parameter :: nbytes = len(header) + product(vsize) * 8
   integer(i4), parameter :: global_header_size = 46 + path_size
   integer(i4), parameter :: global_header_offset = 30 + path_size + nbytes

   integer :: io, stat
   character(len=*), parameter :: filename = ".read-file-not-found.npz"
   integer(i4) :: crc_expected
   real(dp), allocatable :: input(:), output(:)
   character(len=:), allocatable :: msg

   allocate(input(vsize(1)))
   call random_number(input)

   crc_expected = crc32_hash(header)
   crc_expected = crc32_hash(input, crc_expected)

   open(newunit=io, file=filename, form="unformatted", access="stream")
   write(io) local_header_start
   write(io) crc_expected, nbytes, nbytes, path_size, 0_i2
   write(io) pathname
   write(io) header
   write(io) input
   write(io) global_header_start
   write(io) crc_expected, nbytes, nbytes, path_size, 0_i2
   write(io) global_header_final
   write(io) pathname
   write(io) footer_start
   write(io) global_header_size , global_header_offset, 0_i2
   close(io)

   call load_npz(filename, "other_name", output, stat, msg)
   call delete_file(filename)

   if (stat /= 0) call test_failed(error, msg)
end subroutine test_read_file_not_found

subroutine test_read_checksum_mismatch(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: dict = &
       "{'descr':'<f8', 'shape' : ( 61 ,  ) , 'fortran_order': False}         " //  &
       char(10)
   character(len=*), parameter :: header = &
       char(int(z"93")) // "NUMPY" // char(1) // char(0) // &
       char(len(dict)) // char(0) // dict

   integer, parameter :: vsize(*) = [61]

   character(len=*), parameter :: varname = "test", pathname = varname // ".npy"
   integer(i2), parameter :: path_size = len(pathname)
   integer(i4), parameter :: nbytes = len(header) + product(vsize) * 8
   integer(i4), parameter :: global_header_size = 46 + path_size
   integer(i4), parameter :: global_header_offset = 30 + path_size + nbytes

   integer :: io, stat
   character(len=*), parameter :: filename = ".read-checksum-mismatch.npz"
   integer(i4) :: crc_expected
   real(dp), allocatable :: input(:), output(:)
   character(len=:), allocatable :: msg

   allocate(input(vsize(1)))
   call random_number(input)

   crc_expected = -1

   open(newunit=io, file=filename, form="unformatted", access="stream")
   write(io) local_header_start
   write(io) crc_expected, nbytes, nbytes, path_size, 0_i2
   write(io) pathname
   write(io) header
   write(io) input
   write(io) global_header_start
   write(io) crc_expected, nbytes, nbytes, path_size, 0_i2
   write(io) global_header_final
   write(io) pathname
   write(io) footer_start
   write(io) global_header_size , global_header_offset, 0_i2
   close(io)

   call load_npz(filename, varname, output, stat, msg)
   call delete_file(filename)

   if (stat /= 0) call test_failed(error, msg)
end subroutine test_read_checksum_mismatch

subroutine test_write_missing_global1(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: dict = &
       "{'descr':'<f8', 'shape' : ( 61 ,  ) , 'fortran_order': False}         " //  &
       char(10)
   character(len=*), parameter :: header = &
       char(int(z"93")) // "NUMPY" // char(1) // char(0) // &
       char(len(dict)) // char(0) // dict

   integer, parameter :: vsize(*) = [61]

   character(len=*), parameter :: varname = "test", pathname = varname // ".npy"
   integer(i2), parameter :: path_size = len(pathname)
   integer(i4), parameter :: nbytes = len(header) + product(vsize) * 8
   integer(i4), parameter :: global_header_size = 46 + path_size
   integer(i4), parameter :: global_header_offset = 30 + path_size + nbytes

   integer :: io, stat
   character(len=*), parameter :: filename = ".write-missing-global1.npz"
   integer(i4) :: crc_expected
   real(dp), allocatable :: input(:), output(:)
   character(len=:), allocatable :: msg

   allocate(input(vsize(1)))
   call random_number(input)

   crc_expected = crc32_hash(header)
   crc_expected = crc32_hash(input, crc_expected)

   open(newunit=io, file=filename, form="unformatted", access="stream")
   write(io) local_header_start
   write(io) crc_expected, nbytes, nbytes, path_size, 0_i2
   write(io) pathname
   write(io) header
   write(io) input
   close(io)

   call save_npz(filename, "other", input, stat, msg)
   call delete_file(filename)

   if (stat /= 0) call test_failed(error, msg)
end subroutine test_write_missing_global1

subroutine test_write_missing_global2(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: dict = &
       "{'descr':'<f8', 'shape' : ( 61 ,  ) , 'fortran_order': False}         " //  &
       char(10)
   character(len=*), parameter :: header = &
       char(int(z"93")) // "NUMPY" // char(1) // char(0) // &
       char(len(dict)) // char(0) // dict

   integer, parameter :: vsize(*) = [61]

   character(len=*), parameter :: varname = "test", pathname = varname // ".npy"
   integer(i2), parameter :: path_size = len(pathname)
   integer(i4), parameter :: nbytes = len(header) + product(vsize) * 8
   integer(i4), parameter :: global_header_size = 46 + path_size
   integer(i4), parameter :: global_header_offset = 30 + path_size + nbytes

   integer :: io, stat
   character(len=*), parameter :: filename = ".write-missing-global2.npz"
   integer(i4) :: crc_expected
   real(dp), allocatable :: input(:), output(:)
   character(len=:), allocatable :: msg

   allocate(input(vsize(1)))
   call random_number(input)

   crc_expected = crc32_hash(header)
   crc_expected = crc32_hash(input, crc_expected)

   open(newunit=io, file=filename, form="unformatted", access="stream")
   write(io) local_header_start
   write(io) crc_expected, nbytes, nbytes, path_size, 0_i2
   write(io) pathname
   write(io) header
   write(io) input
   write(io) global_header_start(:4)
   close(io)

   call save_npz(filename, "other", input, stat, msg)
   call delete_file(filename)

   if (stat /= 0) call test_failed(error, msg)
end subroutine test_write_missing_global2

subroutine test_write_missing_footer1(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: dict = &
       "{'descr':'<f8', 'shape' : ( 61 ,  ) , 'fortran_order': False}         " //  &
       char(10)
   character(len=*), parameter :: header = &
       char(int(z"93")) // "NUMPY" // char(1) // char(0) // &
       char(len(dict)) // char(0) // dict

   integer, parameter :: vsize(*) = [61]

   character(len=*), parameter :: varname = "test", pathname = varname // ".npy"
   integer(i2), parameter :: path_size = len(pathname)
   integer(i4), parameter :: nbytes = len(header) + product(vsize) * 8
   integer(i4), parameter :: global_header_size = 46 + path_size
   integer(i4), parameter :: global_header_offset = 30 + path_size + nbytes

   integer :: io, stat
   character(len=*), parameter :: filename = ".write-missing-footer1.npz"
   integer(i4) :: crc_expected
   real(dp), allocatable :: input(:), output(:)
   character(len=:), allocatable :: msg

   allocate(input(vsize(1)))
   call random_number(input)

   crc_expected = crc32_hash(header)
   crc_expected = crc32_hash(input, crc_expected)

   open(newunit=io, file=filename, form="unformatted", access="stream")
   write(io) local_header_start
   write(io) crc_expected, nbytes, nbytes, path_size, 0_i2
   write(io) pathname
   write(io) header
   write(io) input
   write(io) global_header_start
   write(io) crc_expected, nbytes, nbytes, path_size, 0_i2
   write(io) global_header_final
   write(io) pathname
   close(io)

   call save_npz(filename, "other", input, stat, msg)
   call delete_file(filename)

   if (stat /= 0) call test_failed(error, msg)
end subroutine test_write_missing_footer1

subroutine test_write_missing_footer2(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: dict = &
       "{'descr':'<f8', 'shape' : ( 61 ,  ) , 'fortran_order': False}         " //  &
       char(10)
   character(len=*), parameter :: header = &
       char(int(z"93")) // "NUMPY" // char(1) // char(0) // &
       char(len(dict)) // char(0) // dict

   integer, parameter :: vsize(*) = [61]

   character(len=*), parameter :: varname = "test", pathname = varname // ".npy"
   integer(i2), parameter :: path_size = len(pathname)
   integer(i4), parameter :: nbytes = len(header) + product(vsize) * 8
   integer(i4), parameter :: global_header_size = 46 + path_size
   integer(i4), parameter :: global_header_offset = 30 + path_size + nbytes

   integer :: io, stat
   character(len=*), parameter :: filename = ".write-missing-footer2.npz"
   integer(i4) :: crc_expected
   real(dp), allocatable :: input(:), output(:)
   character(len=:), allocatable :: msg

   allocate(input(vsize(1)))
   call random_number(input)

   crc_expected = crc32_hash(header)
   crc_expected = crc32_hash(input, crc_expected)

   open(newunit=io, file=filename, form="unformatted", access="stream")
   write(io) local_header_start
   write(io) crc_expected, nbytes, nbytes, path_size, 0_i2
   write(io) pathname
   write(io) header
   write(io) input
   write(io) global_header_start
   write(io) crc_expected, nbytes, nbytes, path_size, 0_i2
   write(io) global_header_final
   write(io) pathname
   write(io) footer_start(:4)
   close(io)

   call save_npz(filename, "other", input, stat, msg)
   call delete_file(filename)

   if (stat /= 0) call test_failed(error, msg)
end subroutine test_write_missing_footer2

subroutine test_write_mismatch_header(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: dict = &
       "{'descr':'<f8', 'shape' : ( 61 ,  ) , 'fortran_order': False}         " //  &
       char(10)
   character(len=*), parameter :: header = &
       char(int(z"93")) // "NUMPY" // char(1) // char(0) // &
       char(len(dict)) // char(0) // dict

   integer, parameter :: vsize(*) = [61]

   character(len=*), parameter :: varname = "test", pathname = varname // ".npy"
   integer(i2), parameter :: path_size = len(pathname)
   integer(i4), parameter :: nbytes = len(header) + product(vsize) * 8
   integer(i4), parameter :: global_header_size = 46 + path_size
   integer(i4), parameter :: global_header_offset = 30 + path_size + nbytes

   integer :: io, stat
   character(len=*), parameter :: filename = ".write-mismatch-header.npz"
   integer(i4) :: crc_expected
   real(dp), allocatable :: input(:), output(:)
   character(len=:), allocatable :: msg

   allocate(input(vsize(1)))
   call random_number(input)

   crc_expected = crc32_hash(header)
   crc_expected = crc32_hash(input, crc_expected)

   open(newunit=io, file=filename, form="unformatted", access="stream")
   write(io) local_header_start
   write(io) crc_expected, nbytes, nbytes, path_size, 0_i2
   write(io) pathname
   write(io) header
   write(io) input
   write(io) footer_start
   close(io)

   call save_npz(filename, "other", input, stat, msg)
   call delete_file(filename)

   if (stat /= 0) call test_failed(error, msg)
end subroutine test_write_mismatch_header

subroutine test_write_mismatch_footer(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: dict = &
       "{'descr':'<f8', 'shape' : ( 61 ,  ) , 'fortran_order': False}         " //  &
       char(10)
   character(len=*), parameter :: header = &
       char(int(z"93")) // "NUMPY" // char(1) // char(0) // &
       char(len(dict)) // char(0) // dict

   integer, parameter :: vsize(*) = [61]

   character(len=*), parameter :: varname = "test", pathname = varname // ".npy"
   integer(i2), parameter :: path_size = len(pathname)
   integer(i4), parameter :: nbytes = len(header) + product(vsize) * 8
   integer(i4), parameter :: global_header_size = 46 + path_size
   integer(i4), parameter :: global_header_offset = 30 + path_size + nbytes

   integer :: io, stat
   character(len=*), parameter :: filename = ".write-mismatch-footer.npz"
   integer(i4) :: crc_expected
   real(dp), allocatable :: input(:), output(:)
   character(len=:), allocatable :: msg

   allocate(input(vsize(1)))
   call random_number(input)

   crc_expected = crc32_hash(header)
   crc_expected = crc32_hash(input, crc_expected)

   open(newunit=io, file=filename, form="unformatted", access="stream")
   write(io) local_header_start
   write(io) crc_expected, nbytes, nbytes, path_size, 0_i2
   write(io) pathname
   write(io) header
   write(io) input
   write(io) global_header_start
   write(io) crc_expected, nbytes, nbytes, path_size, 0_i2
   write(io) global_header_final
   write(io) pathname
   write(io) local_header_start(:4)
   close(io)

   call save_npz(filename, "other", input, stat, msg)
   call delete_file(filename)

   if (stat /= 0) call test_failed(error, msg)
end subroutine test_write_mismatch_footer

subroutine delete_file(filename)
   character(len=*), intent(in) :: filename

   integer :: io

   open(newunit=io, file=filename)
   close(io, status="delete")
end subroutine delete_file

subroutine hexdump(filename)
   character(len=*), intent(in) :: filename

   integer :: io, stat, line, i
   logical :: skip
   character(len=1) :: chunk(16)
   character(len=*), parameter :: hex_format(16) = [&
      '(z0.8, 1x, 1(1x, z0.2),                  t61, "|",  1a, "|")', &
      '(z0.8, 1x, 2(1x, z0.2),                  t61, "|",  2a, "|")', &
      '(z0.8, 1x, 3(1x, z0.2),                  t61, "|",  3a, "|")', &
      '(z0.8, 1x, 4(1x, z0.2),                  t61, "|",  4a, "|")', &
      '(z0.8, 1x, 5(1x, z0.2),                  t61, "|",  5a, "|")', &
      '(z0.8, 1x, 6(1x, z0.2),                  t61, "|",  6a, "|")', &
      '(z0.8, 1x, 7(1x, z0.2),                  t61, "|",  7a, "|")', &
      '(z0.8, 1x, 8(1x, z0.2),                  t61, "|",  8a, "|")', &
      '(z0.8, 1x, 8(1x, z0.2), 1x, 1(1x, z0.2), t61, "|",  9a, "|")', &
      '(z0.8, 1x, 8(1x, z0.2), 1x, 2(1x, z0.2), t61, "|", 10a, "|")', &
      '(z0.8, 1x, 8(1x, z0.2), 1x, 3(1x, z0.2), t61, "|", 11a, "|")', &
      '(z0.8, 1x, 8(1x, z0.2), 1x, 4(1x, z0.2), t61, "|", 12a, "|")', &
      '(z0.8, 1x, 8(1x, z0.2), 1x, 5(1x, z0.2), t61, "|", 13a, "|")', &
      '(z0.8, 1x, 8(1x, z0.2), 1x, 6(1x, z0.2), t61, "|", 14a, "|")', &
      '(z0.8, 1x, 8(1x, z0.2), 1x, 7(1x, z0.2), t61, "|", 15a, "|")', &
      '(z0.8, 1x, 8(1x, z0.2), 1x, 8(1x, z0.2), t61, "|", 16a, "|")']

   print '(3a)', "Hexdump of file '", filename, "'"
   open(newunit=io, file=filename, form="unformatted", access="stream", status="old")
   line = 0
   skip = .false.
   read(io, iostat=stat) chunk
   do while(stat == 0)
      if (all(chunk == ' ')) then
         if (.not.skip) print '(A)', "*"
         skip = .true.
      else
         print hex_format(16), &
            & line, chunk, merge(chunk, '.', is_printable(chunk))
         skip = .false.
      end if
      line = line + 16
      read(io, pos=line + 1, iostat=stat) chunk
   end do
   do i = 1, 16
      read(io, pos=line + i, iostat=stat) chunk(i)
      if (stat /= 0) then
         print hex_format(i-1), &
            & line, chunk(:i-1), merge(chunk(:i-1), '.', is_printable(chunk(:i-1)))
         exit
      end if
   end do
   close(io)
end subroutine hexdump

elemental logical function is_printable(c)
   character(len=1), intent(in) :: c !! The character to test.
   integer :: ic
   ic = iachar(c)
   !The character is printable if it's between ' ' and '~' in the ASCII table
   is_printable = ic >= iachar(' ') .and. ic <= int(z'7E')
end function is_printable

end module test_npz