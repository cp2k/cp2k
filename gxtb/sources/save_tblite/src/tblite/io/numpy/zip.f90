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

!> @file tblite/io/numpy/zip.f90
!> Provides zip implementation for numpy input/output.

!> Provide routines for handling zip files
module tblite_io_numpy_zip
   use mctc_env, only : i2, i4
   use tblite_io_numpy_constants, only : &
      & zip_global_sig, zip_local_sig, zip_footer_sig, zip_min_version
   use tblite_output_format, only : format_string
   implicit none
   private

   public :: list_zip_file, zip_file

   type :: zip_record
      character(len=:), allocatable :: path
   end type zip_record

   type :: zip_file
      type(zip_record), allocatable :: records(:)
      character(len=:), allocatable :: global_header
      integer(i4) :: global_header_offset = 0_i4
      integer(i2) :: nrecs = 0_i2
   end type zip_file


contains

subroutine list_zip_file(io, filename, zip, stat, msg)
   integer, intent(in) :: io
   character(len=*), intent(in) :: filename
   type(zip_file), intent(out) :: zip
   integer, intent(out) :: stat
   character(len=:), allocatable :: msg

   integer :: irec
   integer(i2) :: path_size, extra_field_size, comment_size
   integer(i2) :: disk_no, disk_start, nrecs_on_disk
   integer(i4) :: nbytes_compressed, global_header_size
   character(len=512) :: errmsg
   integer :: res, length, pos
   integer(i4) :: header_sig
   character(len=:), allocatable :: path

   stat = 0
   irec = 0
   pos = 1
   read(io, pos=pos, iostat=stat, iomsg=errmsg) header_sig
   do while(stat == 0 .and. is_local_header(header_sig))
      irec = irec + 1

      if (stat == 0) &
         read(io, pos=pos+18, iostat=stat, iomsg=errmsg) nbytes_compressed
      if (stat == 0) &
         read(io, pos=pos+26, iostat=stat, iomsg=errmsg) path_size
      if (stat == 0) &
         read(io, pos=pos+28, iostat=stat, iomsg=errmsg) extra_field_size

      if (stat == 0) then
         if (allocated(path)) deallocate(path)
         allocate(character(len=path_size) :: path, stat=stat)
      end if
      if (stat == 0) &
         read(io, pos=pos+30, iostat=stat, iomsg=errmsg) path

      pos = pos + 30 + path_size + extra_field_size + nbytes_compressed
      read(io, pos=pos, iostat=stat, iomsg=errmsg) header_sig
   end do
   if (stat /= 0) then
      msg = "Failed to read local header block for '"//filename//"'"
      if (len_trim(errmsg) > 0) &
         msg = msg // " ("//trim(errmsg)//")"
      return
   end if

   if (.not.is_global_header(header_sig)) then
      stat = 400
      msg = "Expected global header signature for '"//filename//"' got "// &
         & format_string(header_sig, '(z0.8)')
      return
   end if

   allocate(zip%records(irec))

   irec = 0
   ! global_header_offset = pos - 1
   do while(stat == 0 .and. is_global_header(header_sig))
      irec = irec + 1

      if (stat == 0) &
         read(io, pos=pos+28, iostat=stat, iomsg=errmsg) path_size
      if (stat == 0) &
         read(io, pos=pos+30, iostat=stat, iomsg=errmsg) extra_field_size
      if (stat == 0) &
         read(io, pos=pos+32, iostat=stat, iomsg=errmsg) comment_size

      if (stat == 0) then
         if (allocated(path)) deallocate(path)
         allocate(character(len=path_size) :: path, stat=stat)
      end if
      if (stat == 0) &
         read(io, pos=pos+46, iostat=stat, iomsg=errmsg) path

      zip%records(irec)%path = path

      pos = pos + 46 + path_size + extra_field_size + comment_size
      read(io, pos=pos, iostat=stat, iomsg=errmsg) header_sig
   end do
   if (stat /= 0) then
      msg = "Failed to read global header block for '"//filename//"'"
      if (len_trim(errmsg) > 0) &
         msg = msg // " ("//trim(errmsg)//")"
      return
   end if
   if (.not.is_footer_header(header_sig)) then
      stat = 401
      msg = "Expected footer signature for '"//filename//"' got "// &
         & format_string(header_sig, '(z0.8)')
      return
   end if
   ! global_header_size = pos - global_header_offset + 1

   read(io, pos=pos+4, iostat=stat, iomsg=errmsg) &
      & disk_no, disk_start, nrecs_on_disk, zip%nrecs, &
      & global_header_size, zip%global_header_offset, comment_size

   if (stat == 0) &
      allocate(character(len=global_header_size) :: zip%global_header, stat=stat)
   if (stat == 0) &
      read(io, iostat=stat, pos=zip%global_header_offset+1) zip%global_header

   if (stat /= 0) then
      msg = "Failed to read footer for '"//filename//"'"
      if (len_trim(errmsg) > 0) &
         msg = msg // " ("//trim(errmsg)//")"
      return
   end if

   if (disk_no /= 0) then
      stat = 402
      msg = "Cannot read zip file with disk_no != 0"
   end if

   if (disk_start /= 0) then
      stat = 403
      msg = "Cannot read zip file with disk_start != 0"
   end if

   if (nrecs_on_disk /= zip%nrecs) then
      stat = 404
      msg = "Cannot read zip file with nrecs_on_disk != nrecs"
   end if
end subroutine list_zip_file

pure function is_local_header(header_sig) result(is_local)
   integer(i4), intent(in) :: header_sig
   logical :: is_local

   is_local = header_sig == zip_local_sig
end function is_local_header

pure function is_global_header(header_sig) result(is_global)
   integer(i4), intent(in) :: header_sig
   logical :: is_global

   is_global = header_sig == zip_global_sig
end function is_global_header

pure function is_footer_header(header_sig) result(is_footer)
   integer(i4), intent(in) :: header_sig
   logical :: is_footer

   is_footer = header_sig == zip_footer_sig
end function is_footer_header
end module tblite_io_numpy_zip