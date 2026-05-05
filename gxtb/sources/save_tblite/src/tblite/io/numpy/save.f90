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

!> @file tblite/io/numpy/save.f90
!> Provides npy output routines

!> Implementation of saving multidimensional arrays to npy files
module tblite_io_numpy_save
   use mctc_env, only : dp, i4
   use tblite_io_numpy_constants, only : magic_number, magic_string, type_rdp, type_i4
   use tblite_output_format, only : format_string
   implicit none
   private

   public :: save_npy, npy_header

   interface save_npy
      module procedure save_npy_i4_r1
      module procedure save_npy_rdp_r1
      module procedure save_npy_rdp_r2
      module procedure save_npy_rdp_r3
   end interface save_npy

   character(len=*), parameter :: nl = achar(10)

contains

!> Generate magic header string for npy format
pure function magic_header(major, minor) result(str)
   !> Major version of npy format
   integer, intent(in) :: major
   !> Minor version of npy format
   integer, intent(in) :: minor
   !> Magic string for npy format
   character(len=8) :: str

   str = magic_number // magic_string // achar(major) // achar(minor)
end function magic_header


!> Generate header for npy format
pure function npy_header(vtype, vshape) result(str)
   !> Type of variable
   character(len=*), intent(in) :: vtype
   !> Shape of variable
   integer, intent(in) :: vshape(:)
   !> Header string for npy format
   character(len=:), allocatable :: str

   integer, parameter :: len_v10 = 8 + 2, len_v20 = 8 + 4, block_size = 64

   str = &
      "{'descr': '"//vtype//&
      "', 'fortran_order': True, 'shape': "//&
      shape_str(vshape)//", }"

    if (len(str) + len_v10 >= 65535) then
      str = str // &
         & repeat(" ", block_size - mod(len(str) + len_v20 + 1, block_size)) // nl
      str = magic_header(2, 0) // to_bytes_i4(int(len(str))) // str
    else
      str = str // &
         & repeat(" ", block_size - mod(len(str) + len_v10 + 1, block_size)) // nl
      str = magic_header(1, 0) // to_bytes_i2(int(len(str))) // str
    end if
end function npy_header

!> Write integer as byte string in little endian encoding
pure function to_bytes_i4(val) result(str)
   !> Integer value to convert to bytes
   integer, intent(in) :: val
   !> String of bytes
   character(len=4) :: str

   str = achar(mod(val, 256**1)) // &
      & achar(mod(val, 256**2) / 256**1) // &
      & achar(mod(val, 256**3) / 256**2) // &
      & achar(val / 256**3)
end function to_bytes_i4


!> Write integer as byte string in little endian encoding, 2-byte truncated version
pure function to_bytes_i2(val) result(str)
   !> Integer value to convert to bytes
   integer, intent(in) :: val
   !> String of bytes
   character(len=2) :: str

   str = achar(mod(val, 2**8)) // &
      & achar(mod(val, 2**16) / 2**8)
end function to_bytes_i2


!> Print array shape as tuple of int
pure function shape_str(vshape) result(str)
   !> Shape of variable
   integer, intent(in) :: vshape(:)
   !> Shape string for npy format
   character(len=:), allocatable :: str

   integer :: i

   str = "("
   do i = 1, size(vshape)
      str = str//format_string(vshape(i), '(i0)')//", "
   enddo
   str = str//")"
end function shape_str


!> Save 1-dimensional array in npy format
subroutine save_npy_i4_r1(filename, array, iostat, iomsg)
   !> Name of the npy file to load from
   character(len=*), intent(in) :: filename
   !> Array to be loaded from the npy file
   integer(i4), intent(in) :: array(:)
   !> Error status of loading, zero on success
   integer, intent(out), optional :: iostat
   !> Associated error message in case of non-zero status code
   character(len=:), allocatable, intent(out), optional :: iomsg

   character(len=*), parameter :: vtype = type_i4
   integer :: io, stat

   open(newunit=io, file=filename, form="unformatted", access="stream", iostat=stat)
   if (stat == 0) then
      write(io, iostat=stat) npy_header(vtype, shape(array))
   end if
   if (stat == 0) then
      write(io, iostat=stat) array
   end if
   close(io, iostat=stat)

   call handle_iostat(stat, filename, iostat, iomsg)
end subroutine save_npy_i4_r1

!> Save 1-dimensional array in npy format
subroutine save_npy_rdp_r1(filename, array, iostat, iomsg)
   !> Name of the npy file to load from
   character(len=*), intent(in) :: filename
   !> Array to be loaded from the npy file
   real(dp), intent(in) :: array(:)
   !> Error status of loading, zero on success
   integer, intent(out), optional :: iostat
   !> Associated error message in case of non-zero status code
   character(len=:), allocatable, intent(out), optional :: iomsg

   character(len=*), parameter :: vtype = type_rdp
   integer :: io, stat

   open(newunit=io, file=filename, form="unformatted", access="stream", iostat=stat)
   if (stat == 0) then
      write(io, iostat=stat) npy_header(vtype, shape(array))
   end if
   if (stat == 0) then
      write(io, iostat=stat) array
   end if
   close(io, iostat=stat)

   call handle_iostat(stat, filename, iostat, iomsg)
end subroutine save_npy_rdp_r1

!> Save 2-dimensional array in npy format
subroutine save_npy_rdp_r2(filename, array, iostat, iomsg)
   !> Name of the npy file to load from
   character(len=*), intent(in) :: filename
   !> Array to be loaded from the npy file
   real(dp), intent(in) :: array(:, :)
   !> Error status of loading, zero on success
   integer, intent(out), optional :: iostat
   !> Associated error message in case of non-zero status code
   character(len=:), allocatable, intent(out), optional :: iomsg

   character(len=*), parameter :: vtype = type_rdp
   integer :: io, stat

   open(newunit=io, file=filename, form="unformatted", access="stream", iostat=stat)
   if (stat == 0) then
      write(io, iostat=stat) npy_header(vtype, shape(array))
   end if
   if (stat == 0) then
      write(io, iostat=stat) array
   end if
   close(io, iostat=stat)

   call handle_iostat(stat, filename, iostat, iomsg)
end subroutine save_npy_rdp_r2

!> Save 3-dimensional array in npy format
subroutine save_npy_rdp_r3(filename, array, iostat, iomsg)
   !> Name of the npy file to load from
   character(len=*), intent(in) :: filename
   !> Array to be loaded from the npy file
   real(dp), intent(in) :: array(:, :, :)
   !> Error status of loading, zero on success
   integer, intent(out), optional :: iostat
   !> Associated error message in case of non-zero status code
   character(len=:), allocatable, intent(out), optional :: iomsg

   character(len=*), parameter :: vtype = type_rdp
   integer :: io, stat

   open(newunit=io, file=filename, form="unformatted", access="stream", iostat=stat)
   if (stat == 0) then
      write(io, iostat=stat) npy_header(vtype, shape(array))
   end if
   if (stat == 0) then
      write(io, iostat=stat) array
   end if
   close(io, iostat=stat)

   call handle_iostat(stat, filename, iostat, iomsg)
end subroutine save_npy_rdp_r3

!> Handle iostat and iomsg of write operation
subroutine handle_iostat(stat, filename, iostat, iomsg)
   !> Error status of loading, zero on success
   integer, intent(in) :: stat
   !> Name of the npy file to load from
   character(len=*), intent(in) :: filename
   !> Error status of loading, zero on success
   integer, intent(out), optional :: iostat
   !> Associated error message in case of non-zero status code
   character(len=:), allocatable, intent(out), optional :: iomsg

   if (present(iostat)) then
      iostat = stat
   else if (stat /= 0) then
      error stop "Failed to write array to file '"//filename//"'"
   end if

   if (present(iomsg)) then
      if (stat /= 0) then
         iomsg = "Failed to write array to file '"//filename//"'"
      end if
   end if
end subroutine handle_iostat

end module tblite_io_numpy_save