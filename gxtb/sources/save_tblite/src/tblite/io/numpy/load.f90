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
module tblite_io_numpy_load
   use mctc_env, only : dp, i4
   use tblite_io_numpy_constants, only : magic_number, magic_string, type_rdp, type_i4
   use tblite_io_numpy_utils, only : reader_type, new_reader, read, delete_reader
   use tblite_output_format, only : format_string
   implicit none
   private

   public :: load_npy
   public :: get_npy_descriptor

   interface load_npy
      module procedure load_npy_i4_r1
      module procedure load_npy_i4_r1_file
      module procedure load_npy_rdp_r1
      module procedure load_npy_rdp_r1_file
      module procedure load_npy_rdp_r2
      module procedure load_npy_rdp_r2_file
      module procedure load_npy_rdp_r3
      module procedure load_npy_rdp_r3_file
   end interface load_npy

   character(len=*), parameter :: nl = achar(10)

   enum, bind(c)
       enumerator :: invalid, string, lbrace, rbrace, comma, colon, &
           lparen, rparen, bool, literal, space
   end enum

   type :: token_type
       integer :: first, last, kind
   end type token_type

contains

subroutine load_npy_i4_r1_file(filename, array, iostat, iomsg)
   !> Name of the npy file to load from
   character(len=*), intent(in) :: filename
   !> Array to be loaded from the npy file
   integer(i4), allocatable, intent(out) :: array(:)
   !> Error status of loading, zero on success
   integer, intent(out), optional :: iostat
   !> Associated error message in case of non-zero status code
   character(len=:), allocatable, intent(out), optional :: iomsg

   type(reader_type) :: reader
   integer :: stat
   character(len=:), allocatable :: msg

   call new_reader(reader, filename, stat)
   if (stat /= 0) then
      call handle_iostat(stat, msg, filename, iostat)
      return
   end if
   call load_npy(reader, array, iostat, iomsg)
   call delete_reader(reader)
end subroutine load_npy_i4_r1_file

subroutine load_npy_i4_r1(reader, array, iostat, iomsg)
   !> Reader for the binary file
   type(reader_type), intent(inout) :: reader
   !> Array to be loaded from the npy file
   integer(i4), allocatable, intent(out) :: array(:)
   !> Error status of loading, zero on success
   integer, intent(out), optional :: iostat
   !> Associated error message in case of non-zero status code
   character(len=:), allocatable, intent(out), optional :: iomsg

   character(len=*), parameter :: vtype = type_i4
   integer, parameter :: rank = 1

   integer :: stat
   character(len=:), allocatable :: msg

   catch: block
       character(len=:), allocatable :: this_type
       integer, allocatable :: vshape(:)

       call get_npy_descriptor(reader, this_type, vshape, stat, msg)
       if (stat /= 0) exit catch

       if (this_type /= vtype) then
           stat = 11
           msg = type_error(reader%filename, this_type, vtype)
           exit catch
       end if

       if (size(vshape) /= rank) then
           stat = 12
           msg = rank_error(reader%filename, vshape, rank)
           exit catch
       end if

       allocate(array(vshape(1)), stat=stat)
       if (stat /= 0) then
           msg = allocation_error(vtype, vshape)
           exit catch
       end if

       call read(reader, array, stat)
   end block catch
   call handle_iostat(stat, msg, reader%filename, iostat)

   if (present(iomsg).and.allocated(msg)) call move_alloc(msg, iomsg)
end subroutine load_npy_i4_r1

subroutine load_npy_rdp_r1_file(filename, array, iostat, iomsg)
   !> Name of the npy file to load from
   character(len=*), intent(in) :: filename
   !> Array to be loaded from the npy file
   real(dp), allocatable, intent(out) :: array(:)
   !> Error status of loading, zero on success
   integer, intent(out), optional :: iostat
   !> Associated error message in case of non-zero status code
   character(len=:), allocatable, intent(out), optional :: iomsg

   type(reader_type) :: reader
   integer :: stat
   character(len=:), allocatable :: msg

   call new_reader(reader, filename, stat)
   if (stat /= 0) then
      call handle_iostat(stat, msg, filename, iostat)
      return
   end if
   call load_npy(reader, array, iostat, iomsg)
   call delete_reader(reader)
end subroutine load_npy_rdp_r1_file

subroutine load_npy_rdp_r1(reader, array, iostat, iomsg)
   !> Reader for the binary file
   type(reader_type), intent(inout) :: reader
   !> Array to be loaded from the npy file
   real(dp), allocatable, intent(out) :: array(:)
   !> Error status of loading, zero on success
   integer, intent(out), optional :: iostat
   !> Associated error message in case of non-zero status code
   character(len=:), allocatable, intent(out), optional :: iomsg

   character(len=*), parameter :: vtype = type_rdp
   integer, parameter :: rank = 1

   integer :: stat
   character(len=:), allocatable :: msg

   catch: block
       character(len=:), allocatable :: this_type
       integer, allocatable :: vshape(:)

       call get_npy_descriptor(reader, this_type, vshape, stat, msg)
       if (stat /= 0) exit catch

       if (this_type /= vtype) then
           stat = 11
           msg = type_error(reader%filename, this_type, vtype)
           exit catch
       end if

       if (size(vshape) /= rank) then
           stat = 12
           msg = rank_error(reader%filename, vshape, rank)
           exit catch
       end if

       allocate(array(vshape(1)), stat=stat)
       if (stat /= 0) then
           msg = allocation_error(vtype, vshape)
           exit catch
       end if

       call read(reader, array, stat)
   end block catch
   call handle_iostat(stat, msg, reader%filename, iostat)

   if (present(iomsg).and.allocated(msg)) call move_alloc(msg, iomsg)
end subroutine load_npy_rdp_r1

subroutine load_npy_rdp_r2_file(filename, array, iostat, iomsg)
   !> Name of the npy file to load from
   character(len=*), intent(in) :: filename
   !> Array to be loaded from the npy file
   real(dp), allocatable, intent(out) :: array(:, :)
   !> Error status of loading, zero on success
   integer, intent(out), optional :: iostat
   !> Associated error message in case of non-zero status code
   character(len=:), allocatable, intent(out), optional :: iomsg

   type(reader_type) :: reader
   integer :: stat
   character(len=:), allocatable :: msg

    call new_reader(reader, filename, stat)
   if (stat /= 0) then
      call handle_iostat(stat, msg, filename, iostat)
      return
   end if
    call load_npy(reader, array, iostat, iomsg)
    call delete_reader(reader)
end subroutine load_npy_rdp_r2_file

subroutine load_npy_rdp_r2(reader, array, iostat, iomsg)
   !> Reader for the binary file
   type(reader_type), intent(inout) :: reader
   !> Array to be loaded from the npy file
   real(dp), allocatable, intent(out) :: array(:, :)
   !> Error status of loading, zero on success
   integer, intent(out), optional :: iostat
   !> Associated error message in case of non-zero status code
   character(len=:), allocatable, intent(out), optional :: iomsg

   character(len=*), parameter :: vtype = type_rdp
   integer, parameter :: rank = 2

   integer :: stat
   character(len=:), allocatable :: msg

   catch: block
       character(len=:), allocatable :: this_type
       integer, allocatable :: vshape(:)

       call get_npy_descriptor(reader, this_type, vshape, stat, msg)
       if (stat /= 0) exit catch

       if (this_type /= vtype) then
           stat = 11
           msg = type_error(reader%filename, this_type, vtype)
           exit catch
       end if

       if (size(vshape) /= rank) then
           stat = 12
           msg = rank_error(reader%filename, vshape, rank)
           exit catch
       end if

       allocate(array(vshape(1), vshape(2)), stat=stat)
       if (stat /= 0) then
           msg = allocation_error(vtype, vshape)
           exit catch
       end if

       call read(reader, array, stat)
   end block catch
   call handle_iostat(stat, msg, reader%filename, iostat)

   if (present(iomsg).and.allocated(msg)) call move_alloc(msg, iomsg)
end subroutine load_npy_rdp_r2

subroutine load_npy_rdp_r3_file(filename, array, iostat, iomsg)
   !> Name of the npy file to load from
   character(len=*), intent(in) :: filename
   !> Array to be loaded from the npy file
   real(dp), allocatable, intent(out) :: array(:, :, :)
   !> Error status of loading, zero on success
   integer, intent(out), optional :: iostat
   !> Associated error message in case of non-zero status code
   character(len=:), allocatable, intent(out), optional :: iomsg

   type(reader_type) :: reader
   integer :: stat
   character(len=:), allocatable :: msg

   call new_reader(reader, filename, stat)
   if (stat /= 0) then
      call handle_iostat(stat, msg, filename, iostat)
      return
   end if
   call load_npy(reader, array, iostat, iomsg)
   call delete_reader(reader)
end subroutine load_npy_rdp_r3_file

subroutine load_npy_rdp_r3(reader, array, iostat, iomsg)
   !> Reader for the binary file
   type(reader_type), intent(inout) :: reader
   !> Array to be loaded from the npy file
   real(dp), allocatable, intent(out) :: array(:, :, :)
   !> Error status of loading, zero on success
   integer, intent(out), optional :: iostat
   !> Associated error message in case of non-zero status code
   character(len=:), allocatable, intent(out), optional :: iomsg

   character(len=*), parameter :: vtype = type_rdp
   integer, parameter :: rank = 3

   integer :: stat
   character(len=:), allocatable :: msg

   catch: block
       character(len=:), allocatable :: this_type
       integer, allocatable :: vshape(:)

       call get_npy_descriptor(reader , this_type, vshape, stat, msg)
       if (stat /= 0) exit catch

       if (this_type /= vtype) then
           stat = 11
           msg = type_error(reader%filename, this_type, vtype)
           exit catch
       end if

       if (size(vshape) /= rank) then
           stat = 12
           msg = rank_error(reader%filename, vshape, rank)
           exit catch
       end if

       allocate(array(vshape(1), vshape(2), vshape(3)), stat=stat)
       if (stat /= 0) then
           msg = allocation_error(vtype, vshape)
           exit catch
       end if

       call read(reader, array, stat)
   end block catch
   call handle_iostat(stat, msg, reader%filename, iostat)

   if (present(iomsg).and.allocated(msg)) call move_alloc(msg, iomsg)
end subroutine load_npy_rdp_r3

!> Read the npy header from a binary file and retrieve the descriptor string.
subroutine get_npy_descriptor(reader, vtype, vshape, stat, msg)
   !> Reader for the binary file
   type(reader_type), intent(inout) :: reader
   !> Type of data saved in npy file
   character(len=:), allocatable, intent(out) :: vtype
   !> Shape descriptor of the
   integer, allocatable, intent(out) :: vshape(:)
   !> Status of operation
   integer, intent(out) :: stat
   !> Associated error message in case of non-zero status
   character(len=:), allocatable, intent(out) :: msg

   integer :: major, header_len, i
   character(len=:), allocatable :: dict
   character(len=8) :: header
   character :: buf4(4), buf2(2)
   logical :: fortran_order
   
   ! stat should be zero if no error occurred
   stat = 0
   
   call read(reader, header, stat)
   if (stat /= 0) return

   call parse_header(header, major, stat, msg)
   if (stat /= 0) return

   if (major > 1) then
      call read(reader, buf4, stat)
      if (stat /= 0) return
      header_len = ichar(buf4(1)) &
         &       + ichar(buf4(2)) * 256**1 &
         &       + ichar(buf4(3)) * 256**2 &
         &       + ichar(buf4(4)) * 256**3
   else
      call read(reader, buf2, stat)
      if (stat /= 0) return
      header_len = ichar(buf2(1)) &
         &       + ichar(buf2(2)) * 256**1
   end if
   allocate(character(header_len) :: dict, stat=stat)
   if (stat /= 0) return

   call read(reader, dict, stat)
   if (stat /= 0) return

   if (dict(header_len:header_len) /= nl) then
       stat = 101
       msg = "Descriptor length does not match"
       return
   end if

   if (scan(dict, achar(0)) > 0) then
       stat = 102
       msg = "Nul byte not allowed in descriptor string"
       return
   end if

   call parse_descriptor(trim(dict(:len(dict)-1)), reader%filename, &
       & vtype, fortran_order, vshape, stat, msg)
   if (stat /= 0) return

   if (.not.fortran_order) then
       vshape = [(vshape(i), i = size(vshape), 1, -1)]
   end if
end subroutine get_npy_descriptor


!> Parse the first eight bytes of the npy header to verify the data
subroutine parse_header(header, major, stat, msg)
   !> Header of the binary file
   character(len=*), intent(in) :: header
   !> Major version of the npy format
   integer, intent(out) :: major
   !> Status of operation
   integer, intent(out) :: stat
   !> Associated error message in case of non-zero status
   character(len=:), allocatable, intent(out) :: msg

   integer :: minor

   ! stat should be zero if no error occurred
   stat = 0

   if (header(1:1) /= magic_number) then
       stat = 201
       msg = "Expected z'93' but got z'"//format_string(ichar(header(1:1)), '(i0)')//"' "//&
           & "as first byte"
       return
   end if

   if (header(2:6) /= magic_string) then
       stat = 202
       msg = "Expected identifier '"//magic_string//"'"
       return
   end if

   major = ichar(header(7:7))
   if (.not.any(major == [1, 2, 3])) then
       stat = 203
       msg = "Unsupported format major version number '"//format_string(major, '(i0)')//"'"
       return
   end if

   minor = ichar(header(8:8))
   if (minor /= 0) then
       stat = 204
       msg = "Unsupported format version "// &
           & "'"//format_string(major, '(i0)')//"."//format_string(minor, '(i0)')//"'"
       return
   end if
end subroutine parse_header

!> Parse the descriptor in the npy header. This routine implements a minimal
!> non-recursive parser for serialized Python dictionaries.
subroutine parse_descriptor(input, filename, vtype, fortran_order, vshape, stat, msg)
   !> Input string to parse as descriptor
   character(len=*), intent(in) :: input
   !> Filename for error reporting
   character(len=*), intent(in) :: filename
   !> Type of the data stored, retrieved from field `descr`
   character(len=:), allocatable, intent(out) :: vtype
   !> Whether the data is in left layout, retrieved from field `fortran_order`
   logical, intent(out) :: fortran_order
   !> Shape of the stored data, retrieved from field `shape`
   integer, allocatable, intent(out) :: vshape(:)
   !> Status of operation
   integer, intent(out) :: stat
   !> Associated error message in case of non-zero status
   character(len=:), allocatable, intent(out) :: msg

   integer :: pos
   character(len=:), allocatable :: key
   type(token_type) :: token, last
   logical :: has_descr, has_shape, has_fortran_order

   has_descr = .false.
   has_shape = .false.
   has_fortran_order = .false.
   pos = 0
   call next_token(input, filename, pos, token, [lbrace], stat, msg)
   if (stat /= 0) return

   last = token_type(pos, pos, comma)
   do while (pos < len(input))
       call get_token(input, pos, token)
       select case(token%kind)
       case(space)
           continue
       case(comma)
           if (token%kind == last%kind) then
               stat = 301
               msg = make_message(filename, input, token%first, token%last, &
                   & "Comma cannot appear at this point")
               return
           end if
           last = token
       case(rbrace)
           exit
       case(string)
           if (token%kind == last%kind) then
               stat = 302
               msg = make_message(filename, input, token%first, token%last, &
                   & "String cannot appear at this point")
               return
           end if
           last = token

           key = input(token%first+1:token%last-1)
           call next_token(input, filename, pos, token, [colon], stat, msg)
           if (stat /= 0) return

           if (key == "descr" .and. has_descr &
               & .or. key == "fortran_order" .and. has_fortran_order &
               & .or. key == "shape" .and. has_shape) then
               stat = 303
               msg = make_message(filename, input, last%first, last%last, &
                   & "Duplicate entry for '"//key//"' found")
               return
           end if

           select case(key)
           case("descr")
               call next_token(input, filename, pos, token, [string], stat, msg)
               if (stat /= 0) return

               vtype = input(token%first+1:token%last-1)
               has_descr = .true.

           case("fortran_order")
               call next_token(input, filename, pos, token, [bool], stat, msg)
               if (stat /= 0) return

               fortran_order = input(token%first:token%last) == "True"
               has_fortran_order = .true.

           case("shape")
               call parse_tuple(input, filename, pos, vshape, stat, msg)
               if (stat /= 0) return

               has_shape = .true.

           case default
               stat = 304
               msg = make_message(filename, input, last%first, last%last, &
                   & "Invalid entry '"//key//"' in dictionary encountered")
               return
           end select
       case default
           stat = 305
           msg = make_message(filename, input, token%first, token%last, &
               & "Invalid token encountered")
           return
       end select
   end do

   if (.not.has_descr) then
       stat = 306
       msg = make_message(filename, input, 1, pos, &
           & "Dictionary does not contain required entry 'descr'")
   end if

   if (.not.has_shape) then
       stat = 307
       msg = make_message(filename, input, 1, pos, &
           & "Dictionary does not contain required entry 'shape'")
   end if

   if (.not.has_fortran_order) then
       stat = 308
       msg = make_message(filename, input, 1, pos, &
           & "Dictionary does not contain required entry 'fortran_order'")
   end if

end subroutine parse_descriptor

pure function make_message(filename, input, first, last, message) result(str)
   !> Filename for context
   character(len=*), intent(in) :: filename
   !> Input string to parse
   character(len=*), intent(in) :: input
   !> Offset in the input
   integer, intent(in) :: first, last
   !> Error message
   character(len=*), intent(in) :: message
   !> Final output message
   character(len=:), allocatable :: str

   character(len=*), parameter :: nl = new_line('a')

   str = message // nl // &
       & " --> " // filename // ":1:" // format_string(first, '(i0)') // "-" // format_string(last, '(i0)') // nl // &
       & "  |" // nl // &
       & "1 | " // input // nl // &
       & "  |" // repeat(" ", first) // repeat("^", last - first + 1) // nl // &
       & "  |"
end function make_message

!> Parse a tuple of integers into an array of integers
pure subroutine parse_tuple(input, filename, pos, tuple, stat, msg)
   !> Input string to parse
   character(len=*), intent(in) :: input
   !> Filename for error reporting
    character(len=*), intent(in) :: filename
   !> Offset in the input, will be advanced after reading
   integer, intent(inout) :: pos
   !> Array representing tuple of integers
   integer, allocatable, intent(out) :: tuple(:)
   !> Status of operation
   integer, intent(out) :: stat
   !> Associated error message in case of non-zero status
   character(len=:), allocatable, intent(out) :: msg

   type(token_type) :: token
   integer :: last, itmp

   allocate(tuple(0), stat=stat)
   if (stat /= 0) return

   call next_token(input, filename, pos, token, [lparen], stat, msg)
   if (stat /= 0) return

   last = comma
   do while (pos < len(input))
       call get_token(input, pos, token)
       select case(token%kind)
       case(space)
           continue
       case(literal)
           if (token%kind == last) then
               stat = 311
               msg = make_message(filename, input, token%first, token%last, &
                   & "Invalid token encountered")
               return
           end if
           last = token%kind
           read(input(token%first:token%last), *, iostat=stat) itmp
           if (stat /= 0) then
               return
           end if
           tuple = [tuple, itmp]
       case(comma)
           if (token%kind == last) then
               stat = 312
               msg = make_message(filename, input, token%first, token%last, &
                   & "Invalid token encountered")
               return
           end if
           last = token%kind
       case(rparen)
           exit
       case default
           stat = 313
           msg = make_message(filename, input, token%first, token%last, &
               & "Invalid token encountered")
           return
       end select
   end do
end subroutine parse_tuple

!> Get the next allowed token
pure subroutine next_token(input, filename, pos, token, allowed_token, stat, msg)
   !> Input string to parse
   character(len=*), intent(in) :: input
   !> Filename for error reporting
   character(len=*), intent(in) :: filename
   !> Current offset in the input string
   integer, intent(inout) :: pos
   !> Last token parsed
   type(token_type), intent(out) :: token
   !> Tokens allowed in the current context
   integer, intent(in) :: allowed_token(:)
   !> Status of operation
   integer, intent(out) :: stat
   !> Associated error message in case of non-zero status
   character(len=:), allocatable, intent(out) :: msg

   stat = pos
   do while (pos < len(input))
       call get_token(input, pos, token)
       if (token%kind == space) then
           continue
       else if (any(token%kind == allowed_token)) then
           stat = 0
           exit
       else
           stat = 321
           msg = make_message(filename, input, token%first, token%last, &
               & "Invalid token encountered")
           exit
       end if
   end do
end subroutine next_token

!> Tokenize input string
pure subroutine get_token(input, pos, token)
   !> Input strin to tokenize
   character(len=*), intent(in) :: input
   !> Offset in input string, will be advanced
   integer, intent(inout) :: pos
   !> Returned token from the next position
   type(token_type), intent(out) :: token

   character :: quote

   pos = pos + 1
   select case(input(pos:pos))
   case("""", "'")
       quote = input(pos:pos)
       token%first = pos
       pos = pos + 1
       do while (pos <= len(input))
           if (input(pos:pos) == quote) then
               token%last = pos
               exit
           else
               pos = pos + 1
           end if
       end do
       token%kind = string
   case("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")
       token%first = pos
       do while (pos <= len(input))
           if (.not.any(input(pos:pos) == ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"])) then
               pos = pos - 1
               token%last = pos
               exit
           else
               pos = pos + 1
           end if
       end do
       token%kind = literal
   case("T")
       if (starts_with(input(pos:), "True")) then
           token = token_type(pos, pos+3, bool)
           pos = pos + 3
       else
           token = token_type(pos, pos, invalid)
       end if
   case("F")
       if (starts_with(input(pos:), "False")) then
           token = token_type(pos, pos+4, bool)
           pos = pos + 4
       else
           token = token_type(pos, pos, invalid)
       end if
   case("{")
       token = token_type(pos, pos, lbrace)
   case("}")
       token = token_type(pos, pos, rbrace)
   case(",")
       token = token_type(pos, pos, comma)
   case(":")
       token = token_type(pos, pos, colon)
   case("(")
       token = token_type(pos, pos, lparen)
   case(")")
       token = token_type(pos, pos, rparen)
   case(" ", nl)
       token = token_type(pos, pos, space)
   case default
       token = token_type(pos, pos, invalid)
   end select

end subroutine get_token

!> Create an error message for a type mismatch
pure function type_error(filename, this_type, vtype) result(msg)
   !> Filename for error reporting
   character(len=*), intent(in) :: filename
   !> Type of the data stored, retrieved from field `descr`
   character(len=*), intent(in) :: this_type
   !> Type of the data expected
   character(len=*), intent(in) :: vtype
   !> Error message
   character(len=:), allocatable :: msg

   msg = "File '"//filename//"' contains data of type '"//this_type//"', "//&
      & "but expected '"//vtype//"'"
end function type_error

!> Create an error message for a rank mismatch
pure function rank_error(filename, vshape, rank) result(msg)
   !> Filename for error reporting
   character(len=*), intent(in) :: filename
   !> Shape of the data stored, retrieved from field `shape`
   integer, intent(in) :: vshape(:)
   !> Rank of the data expected
   integer, intent(in) :: rank
   !> Error message
   character(len=:), allocatable :: msg

   msg = "File '"//filename//"' contains data of rank "//&
      & format_string(size(vshape), '(i0)')//", but expected "//&
      & format_string(rank, '(i0)')
end function rank_error

!> Create an error message for a failed allocation
pure function allocation_error(vtype, vshape) result(msg)
   !> Type of the data expected
   character(len=*), intent(in) :: vtype
   !> Shape of the data expected
   integer, intent(in) :: vshape(:)
   !> Error message
   character(len=:), allocatable :: msg

   msg = "Failed to allocate array of type '"//vtype//"' "//&
      & "with total size of "//format_string(product(vshape), '(i0)')
end function allocation_error

!> Handle the iostat of the read operation
subroutine handle_iostat(stat, msg, filename, iostat)
   integer, intent(in) :: stat
   character(len=:), allocatable, intent(in) :: msg
   character(len=*), intent(in) :: filename
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

!> Check whether a string starts with substring or not
pure function starts_with(string, substring) result(match)
   character(len=*), intent(in) :: string
   character(len=*), intent(in) :: substring
   logical :: match
   integer :: nsub

   nsub = len(substring)
   if (len(string) < nsub) then
      match = .false.
      return
   end if
   match = string(1:nsub) == substring
end function starts_with

end module tblite_io_numpy_load
