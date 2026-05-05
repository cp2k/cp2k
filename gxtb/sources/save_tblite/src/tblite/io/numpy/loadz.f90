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

!> @file tblite/io/numpy/load.f90
!> Provides npy input routines

!> Implementation of npy input routines
module tblite_io_numpy_loadz
   use mctc_env, only : dp, i2, i4
   use tblite_io_numpy_constants, only : type_rdp, type_i4, &
      & zip_global_sig, zip_local_sig, zip_footer_sig, zip_min_version
   use tblite_io_numpy_load, only : load_npy, get_npy_descriptor
   use tblite_io_numpy_utils, only : reader_type, new_reader, delete_reader
   use tblite_io_numpy_crc32, only : crc32_hash
   use tblite_output_format, only : format_string
   implicit none
   private

   public :: load_npz
   public :: get_npz_descriptor

   !> Interface for loading npz files
   interface load_npz
      module procedure load_npz_i4_r1
      module procedure load_npz_rdp_r1
      module procedure load_npz_rdp_r2
      module procedure load_npz_rdp_r3
   end interface load_npz

   character(len=*), parameter :: nl = achar(10)

contains

!> Get numpy array descriptor from a npz file
subroutine get_npz_descriptor(filename, varname, vtype, vshape, iostat, iomsg)
   !> Filename of the npz file
   character(len=*), intent(in) :: filename
   !> Name of the variable to load
   character(len=*), intent(in) :: varname
   !> Variable type of the numpy array
   character(len=:), allocatable, intent(out) :: vtype
   !> Variable shape of the numpy array
   integer, allocatable, intent(out) :: vshape(:)
   !> Status of the read operation
   integer, intent(out), optional :: iostat
   !> Error message in case of failure
   character(len=:), allocatable, intent(out), optional :: iomsg

   integer :: stat
   character(len=:), allocatable :: msg
   character(len=1), allocatable :: buffer(:)
   type(reader_type) :: reader

   call load_npz_buffer(filename, varname, buffer, stat, msg)
   if (stat == 0) then
      call new_reader(reader, buffer, filename//"::"//varname//".npy")
      call get_npy_descriptor(reader, vtype, vshape, stat, msg)
      call delete_reader(reader)
   end if

   call handle_iostat(stat, msg, filename//"::"//varname//".npy", iostat)
   if (present(iomsg).and.allocated(msg).and.stat /= 0) call move_alloc(msg, iomsg)
end subroutine get_npz_descriptor

!> Load a numpy array from a npz file
subroutine load_npz_i4_r1(filename, varname, array, iostat, iomsg)
   !> Filename of the npz file
   character(len=*), intent(in) :: filename
   !> Name of the variable to load
   character(len=*), intent(in) :: varname
   !> Array to load the data into
   integer(i4), allocatable, intent(out) :: array(:)
   !> Status of the read operation
   integer, intent(out), optional :: iostat
   !> Error message in case of failure
   character(len=:), allocatable, intent(out), optional :: iomsg

   integer :: stat
   character(len=:), allocatable :: msg
   character(len=1), allocatable :: buffer(:)
   type(reader_type) :: reader

   call load_npz_buffer(filename, varname, buffer, stat, msg)
   if (stat == 0) then
      call new_reader(reader, buffer, filename//"::"//varname//".npy")
      call load_npy(reader, array, stat, msg)
      call delete_reader(reader)
   end if

   call handle_iostat(stat, msg, filename//"::"//varname//".npy", iostat)
   if (present(iomsg).and.allocated(msg).and.stat /= 0) call move_alloc(msg, iomsg)
end subroutine load_npz_i4_r1

!> Load a numpy array from a npz file
subroutine load_npz_rdp_r1(filename, varname, array, iostat, iomsg)
   !> Filename of the npz file
   character(len=*), intent(in) :: filename
   !> Name of the variable to load
   character(len=*), intent(in) :: varname
   !> Array to load the data into
   real(dp), allocatable, intent(out) :: array(:)
   !> Status of the read operation
   integer, intent(out), optional :: iostat
   !> Error message in case of failure
   character(len=:), allocatable, intent(out), optional :: iomsg

   integer :: stat
   character(len=:), allocatable :: msg
   character(len=1), allocatable :: buffer(:)
   type(reader_type) :: reader

   call load_npz_buffer(filename, varname, buffer, stat, msg)
   if (stat == 0) then
      call new_reader(reader, buffer, filename//"::"//varname//".npy")
      call load_npy(reader, array, stat, msg)
      call delete_reader(reader)
   end if

   call handle_iostat(stat, msg, filename//"::"//varname//".npy", iostat)
   if (present(iomsg).and.allocated(msg).and.stat /= 0) call move_alloc(msg, iomsg)
end subroutine load_npz_rdp_r1

!> Load a numpy array from a npz file
subroutine load_npz_rdp_r2(filename, varname, array, iostat, iomsg)
   !> Filename of the npz file
   character(len=*), intent(in) :: filename
   !> Name of the variable to load
   character(len=*), intent(in) :: varname
   !> Array to load the data into
   real(dp), allocatable, intent(out) :: array(:, :)
   !> Status of the read operation
   integer, intent(out), optional :: iostat
   !> Error message in case of failure
   character(len=:), allocatable, intent(out), optional :: iomsg

   integer :: stat
   character(len=:), allocatable :: msg
   character(len=1), allocatable :: buffer(:)
   type(reader_type) :: reader

   call load_npz_buffer(filename, varname, buffer, stat, msg)
   if (stat == 0) then
      call new_reader(reader, buffer, filename//"::"//varname//".npy")
      call load_npy(reader, array, stat, msg)
      call delete_reader(reader)
   end if

   call handle_iostat(stat, msg, filename//"::"//varname//".npy", iostat)
   if (present(iomsg).and.allocated(msg).and.stat /= 0) call move_alloc(msg, iomsg)
end subroutine load_npz_rdp_r2

!> Load a numpy array from a npz file
subroutine load_npz_rdp_r3(filename, varname, array, iostat, iomsg)
   !> Filename of the npz file
   character(len=*), intent(in) :: filename
   !> Name of the variable to load
   character(len=*), intent(in) :: varname
   !> Array to load the data into
   real(dp), allocatable, intent(out) :: array(:, :, :)
   !> Status of the read operation
   integer, intent(out), optional :: iostat
   !> Error message in case of failure
   character(len=:), allocatable, intent(out), optional :: iomsg

   integer :: stat
   character(len=:), allocatable :: msg
   character(len=1), allocatable :: buffer(:)
   type(reader_type) :: reader

   call load_npz_buffer(filename, varname, buffer, stat, msg)
   if (stat == 0) then
      call new_reader(reader, buffer, filename//"::"//varname//".npy")
      call load_npy(reader, array, stat, msg)
      call delete_reader(reader)
   end if

   call handle_iostat(stat, msg, filename//"::"//varname//".npy", iostat)
   if (present(iomsg).and.allocated(msg).and.stat /= 0) call move_alloc(msg, iomsg)
end subroutine load_npz_rdp_r3


!> Load a numpy array from a npz file
subroutine load_npz_buffer(filename, varname, buffer, stat, msg)
   !> Filename of the npz file
   character(len=*), intent(in) :: filename
   !> Name of the variable to load
   character(len=*), intent(in) :: varname
   !> Buffer to load the data into
   character(len=1), allocatable, intent(out) :: buffer(:)
   !> Status of the read operation
   integer, intent(out) :: stat
   !> Error message in case of failure
   character(len=:), allocatable, intent(out) :: msg

   character(len=30) :: local_header
   integer(i4) :: nbytes, nbytes_compressed, crc_expected
   integer(i2) :: compression_method
   character(len=:), allocatable :: path
   character(len=512) :: errmsg
   integer :: io
   logical :: exist

   exist = .false.
   open(newunit=io, file=filename, form="unformatted", access="stream", iostat=stat, iomsg=errmsg)
   if (stat /= 0) msg = "Failed to open file '"//filename//"': "//trim(errmsg)
   do while (stat == 0)
      read(io, iostat=stat) local_header
      if (stat /= 0) then
         msg = "Failed to read local header from '"//filename//"'"
         exit
      end if
      if (.not.is_local_header(local_header)) exit

      call read_data(io, local_header, path, nbytes, nbytes_compressed, compression_method, &
         & stat, msg)

      if (stat == 0) then
         if (allocated(buffer)) deallocate(buffer)
         allocate(buffer(nbytes_compressed), stat=stat)
         msg = "Failed to allocate buffer for data from '"//filename//"::"//path//"'"
      end if
      if (stat == 0) then
         read(io, iostat=stat) buffer
         msg = "Failed to read data from '"//filename//"::"//path//"'"
      end if

      if (stat == 0 .and. varname//".npy" == path) then
         exist = .true.
         crc_expected = transfer(local_header(15:18), crc_expected)
         call check_crc32(buffer, crc_expected, filename//"::"//path, stat, msg)
         exit
      end if
   end do
   close(io)

   if (stat == 0 .and. .not.exist) then
      stat = 501
      msg = file_not_found_error(filename, varname)
   end if
end subroutine load_npz_buffer

!> Check if the header is a local header
pure function is_local_header(header) result(is_local)
   !> Header to check
   character(len=30), intent(in) :: header
   !> Status of the header
   logical :: is_local

   integer(i4) :: header_sig

   header_sig = transfer(header(1:4), header_sig)
   is_local = header_sig == zip_local_sig
end function is_local_header

!> Read the data from the local header
subroutine read_data(io, local_header, path, nbytes, nbytes_compressed, compression_method, stat, msg)
   !> Unformatted I/O unit
   integer, intent(in) :: io
   !> Local header
   character(len=30), intent(in) :: local_header
   !> Path to the data
   character(len=:), allocatable, intent(out) :: path
   !> Number of bytes in the data
   integer(i4), intent(out) :: nbytes, nbytes_compressed
   !> Compression method
   integer(i2), intent(out) :: compression_method
   !> Status of the read operation
   integer, intent(out) :: stat
   !> Error message in case of failure
   character(len=:), allocatable, intent(out) :: msg

   integer(i2) :: path_size, extra_field_size
   character(len=:), allocatable :: buffer

   path_size = transfer(local_header(27:28), path_size)
   allocate(character(len=path_size) :: path)
   read(io) path

   extra_field_size = transfer(local_header(29:30), extra_field_size)
   if (extra_field_size > 0_i2) then
      allocate(character(len=extra_field_size) :: buffer)
      read(io) buffer
      deallocate(buffer)
   end if

   compression_method = transfer(local_header(9:10), compression_method)
   nbytes_compressed = transfer(local_header(19:22), nbytes)
   nbytes = transfer(local_header(23:26), nbytes_compressed)
end subroutine read_data

!> Check the CRC32 of the data
subroutine check_crc32(buffer, expected, filename, stat, msg)
   !> Buffer to check
   character(len=1), intent(in) :: buffer(:)
   !> Expected CRC32 value
   integer(i4), intent(in) :: expected
   !> Filename of the data for error message
   character(len=*), intent(in) :: filename
   !> Status of the check
   integer, intent(out) :: stat
   !> Error message in case of failure
   character(len=:), allocatable, intent(out) :: msg

   integer(i4) :: actual

   actual = crc32_hash(buffer)
   if (actual /= expected) then
      stat = 502
      msg = "CRC mismatch for " // filename // &
         & " (expected: " // format_string(expected, '(z08)') // &
         & " actual: " // format_string(actual, '(z08)') // ")"
   end if
end subroutine check_crc32

!> Handle the iostat of the read operation
subroutine handle_iostat(stat, msg, filename, iostat)
   !> Status of the read operation
   integer, intent(in) :: stat
   !> Error message in case of failure
   character(len=:), allocatable, intent(in) :: msg
   !> Filename of the data for error message
   character(len=*), intent(in) :: filename
   !> Status of the read operation
   integer, intent(out), optional :: iostat

   if (present(iostat)) then
       iostat = stat
   else if (stat /= 0) then
       if (allocated(msg)) then
           error stop "Failed to read array from file '"//filename//"'"//nl//&
               & msg
       else
           error stop "Failed to read array from file '"//filename//"'"
       end if
   end if
end subroutine handle_iostat

!> Generate an error message for file not found
pure function file_not_found_error(filename, varname) result(msg)
   !> Filename of the npz file
   character(len=*), intent(in) :: filename
   !> Name of the variable to load
   character(len=*), intent(in) :: varname
   !> Error message in case of failure
   character(len=:), allocatable :: msg

   msg = "File '"//filename//"::"//varname//".npy' not found in npz archive"
end function file_not_found_error

end module tblite_io_numpy_loadz
