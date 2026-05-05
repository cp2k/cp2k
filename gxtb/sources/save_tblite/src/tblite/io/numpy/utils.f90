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

!> @file tblite/io/numpy/utils.f90
!> Provide helper routines for numpy file handling

!> This module provides helper routines for numpy file handling
module tblite_io_numpy_utils
   use iso_fortran_env, only : iostat_end
   use mctc_env, only : sp, dp, i2, i4
   implicit none
   private

   public :: reader_type, new_reader, delete_reader, read

   type :: reader_type
      integer :: io = -1
      character(len=:), allocatable :: filename
      character(len=1), allocatable :: buffer(:)
      integer :: pos = 0
   end type reader_type

   interface new_reader
      module procedure new_reader_from_file
      module procedure new_reader_from_buffer
   end interface new_reader

   interface read
      module procedure read_char_r0
      module procedure read_char_r1
      module procedure read_i4_r1
      module procedure read_rdp_r1
      module procedure read_rdp_r2
      module procedure read_rdp_r3
   end interface read

contains

subroutine new_reader_from_file(reader, filename, stat)
   type(reader_type), intent(out) :: reader
   character(len=*), intent(in) :: filename
   integer, intent(out) :: stat

   logical :: exist

   inquire(file=filename, exist=exist)
   if (.not. exist) then
      stat = 10
      return
   end if
   open(newunit=reader%io, file=filename, form="unformatted", access="stream", iostat=stat)
   reader%filename = filename
end subroutine new_reader_from_file

subroutine new_reader_from_buffer(reader, buffer, filename)
   type(reader_type), intent(out) :: reader
   character(len=1), allocatable :: buffer(:)
   character(len=*), intent(in) :: filename

   reader%filename = filename
   reader%buffer = buffer
   reader%pos = 0
end subroutine new_reader_from_buffer

subroutine read_char_r0(reader, var, stat)
   type(reader_type), intent(inout) :: reader
   character(len=*), intent(out) :: var
   integer, intent(out) :: stat

   if (allocated(reader%buffer)) then
      if (reader%pos + len(var) > size(reader%buffer)) then
         stat = iostat_end
         return
      end if
      var = transfer(reader%buffer(reader%pos+1:reader%pos+len(var)), var)
      reader%pos = reader%pos + len(var)
   else
      read(reader%io, iostat=stat) var
   end if
end subroutine read_char_r0

subroutine delete_reader(reader)
   type(reader_type), intent(inout) :: reader

   reader%pos = 0
   if (allocated(reader%buffer)) then
      deallocate(reader%buffer)
   end if
   if (reader%io /= -1) then
      close(reader%io)
      reader%io = -1
   end if
end subroutine delete_reader

subroutine read_char_r1(reader, var, stat)
   type(reader_type), intent(inout) :: reader
   character(len=1), intent(out) :: var(:)
   integer, intent(out) :: stat

   if (allocated(reader%buffer)) then
      if (reader%pos + size(var) > size(reader%buffer)) then
         stat = iostat_end
         return
      end if
      var = reader%buffer(reader%pos+1:reader%pos+size(var))
      reader%pos = reader%pos + size(var)
   else
      read(reader%io, iostat=stat) var
   end if
end subroutine read_char_r1

subroutine read_i4_r1(reader, var, stat)
   type(reader_type), intent(inout) :: reader
   integer(i4), intent(out) :: var(:)
   integer, intent(out) :: stat

   if (allocated(reader%buffer)) then
      if (reader%pos + size(var) * 4 > size(reader%buffer)) then
         stat = iostat_end
         return
      end if
      var = transfer(reader%buffer(reader%pos+1:reader%pos+size(var)*4), var)
      reader%pos = reader%pos + size(var) * 4
   else
      read(reader%io, iostat=stat) var
   end if
end subroutine read_i4_r1

subroutine read_rdp_r1(reader, var, stat)
   type(reader_type), intent(inout) :: reader
   real(dp), intent(out) :: var(:)
   integer, intent(out) :: stat

   if (allocated(reader%buffer)) then
      if (reader%pos + size(var) * 8 > size(reader%buffer)) then
         stat = iostat_end
         return
      end if
      var = transfer(reader%buffer(reader%pos+1:reader%pos+size(var)*8), var)
      reader%pos = reader%pos + size(var) * 8
   else
      read(reader%io, iostat=stat) var
   end if
end subroutine read_rdp_r1

subroutine read_rdp_r2(reader, var, stat)
   type(reader_type), intent(inout) :: reader
   real(dp), intent(out) :: var(:, :)
   integer, intent(out) :: stat

   if (allocated(reader%buffer)) then
      if (reader%pos + size(var) * 8 > size(reader%buffer)) then
         stat = iostat_end
         return
      end if
      var = reshape(transfer(reader%buffer(reader%pos+1:reader%pos+size(var)*8), var), shape(var))
      reader%pos = reader%pos + size(var) * 8
   else
      read(reader%io, iostat=stat) var
   end if
end subroutine read_rdp_r2

subroutine read_rdp_r3(reader, var, stat)
   type(reader_type), intent(inout) :: reader
   real(dp), intent(out) :: var(:, :, :)
   integer, intent(out) :: stat

   if (allocated(reader%buffer)) then
      if (reader%pos + size(var) * 8 > size(reader%buffer)) then
         stat = iostat_end
         return
      end if
      var = reshape(transfer(reader%buffer(reader%pos+1:reader%pos+size(var)*8), var), shape(var))
      reader%pos = reader%pos + size(var) * 8
   else
      read(reader%io, iostat=stat) var
   end if
end subroutine read_rdp_r3

end module tblite_io_numpy_utils