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

!> @dir tblite/api/double_dictionary
!> Contains API access to a double dictionary instance i.e. to retrieve teh post processing dict.

!> @file tblite/api/double_dictionery.f90
!> Implements API calls to double dictionary.
module tblite_api_double_dictionary
   use mctc_env, only : wp, error_type, fatal_error
   use iso_c_binding
   use tblite_double_dictionary
   use tblite_api_version, only : namespace
   use tblite_api_utils, only : f_c_character, c_f_character
   use tblite_api_error, only : vp_error
   implicit none
   private

   public :: vp_double_dictionary
   public :: get_array_entry_index_api, get_n_entries_dict_api, get_label_entry_index_api, &
      & get_array_size_index_api, get_array_entry_label_api, get_array_size_label_api

   type :: vp_double_dictionary
      type(double_dictionary_type) :: ptr
   end type

   logical, parameter :: debug = .false.
contains

function get_n_entries_dict_api(verror, vdict) result(n) &
      & bind(C, name = namespace//"get_n_entries_dict")
   type(vp_double_dictionary), pointer :: dict
   type(c_ptr), value :: vdict
   type(c_ptr), value :: verror 
   type(vp_error), pointer :: error
   integer(kind=c_int) :: n
   logical :: ok
   if (debug) print '("[Info]", 1x, a)', "get n_entries form dict"
   n = 0
   call check_dict(verror, vdict, error, dict, ok)
   if (.not.(ok)) return
   n = dict%ptr%get_n_entries()
end function

subroutine check_dict(verror, vdict, error, dict, ok)
   type(vp_double_dictionary), pointer :: dict
   type(c_ptr), value :: vdict
   type(c_ptr), value :: verror 
   type(vp_error), pointer :: error
   logical :: ok
   ok = .false.

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vdict)) then
      call fatal_error(error%ptr, "Dictionary is missing.")
      return
   end if
   call c_f_pointer(vdict, dict)
   ok = .true.
end subroutine

subroutine get_array_size_label_api(verror, vdict, label, dim1, dim2, dim3) &
   & bind(C, name = namespace//"get_array_size_label")

   type(vp_double_dictionary), pointer :: dict
   type(c_ptr), value :: vdict
   type(c_ptr), value :: verror 
   type(vp_error), pointer :: error
   character(kind=c_char), intent(in) :: label(*)
   character(len=:), allocatable :: f_char
   integer(kind=c_int), intent(out) :: dim1, dim2, dim3
   logical :: ok
   real(kind=wp), allocatable :: array1(:), array2(:, :), array3(:, :, :)
   integer :: index
   
   if (debug) print '("[Info]", 1x, a)', "get array size by label"

   dim1 = 0
   dim2 = 0
   dim3 = 0

   call check_dict(verror, vdict, error, dict, ok)
   if (.not.(ok)) return
   call c_f_character(label, f_char)
   index = dict%ptr%get_index(f_char)
   if (index == 0) then 
      call fatal_error(error%ptr, "Label is not a key in the dictionary.")
      return
   end if

   call dict%ptr%get_entry(f_char, array1)
   if (allocated(array1)) then
      dim1 = size(array1, dim=1)
      return
   end if

   call dict%ptr%get_entry(f_char, array2)
   if (allocated(array2)) then
      dim1 = size(array2, dim=1)
      dim2 = size(array2, dim=2)
      return
   end if

   call dict%ptr%get_entry(f_char, array3)
   if (allocated(array3)) then
      dim1 = size(array3, dim=1)
      dim2 = size(array3, dim=2)
      dim3 = size(array3, dim=3)
      return
   end if

end subroutine

subroutine get_array_entry_label_api(verror, vdict, label, array_out) &
   & bind(C, name = namespace//"get_array_entry_label")
type(vp_double_dictionary), pointer :: dict
type(c_ptr), value :: vdict
type(c_ptr), value :: verror 
type(vp_error), pointer :: error
character(kind=c_char), intent(in) :: label(*)
character(len=:), allocatable :: f_char
real(c_double), intent(out) :: array_out(*)
real(kind=wp), allocatable :: array1(:), array2(:, :), array3(:, :, :)
integer :: index

logical :: ok
if (debug) print '("[Info]", 1x, a)', "get array values by label"

call check_dict(verror, vdict, error, dict, ok)
if (.not.(ok)) return

call c_f_character(label, f_char)
index = dict%ptr%get_index(f_char)
if (index == 0) then 
   call fatal_error(error%ptr, "Label is not a key in the dictionary.")
   return
end if

call dict%ptr%get_entry(index, array1)
if (allocated(array1)) then
   array_out(:size(array1)) = &
      & reshape(array1, [size(array1)])
   return
end if

call dict%ptr%get_entry(index, array2)
if (allocated(array2)) then
   array_out(:size(array2)) = &
      & reshape(array2, [size(array2)])
   return
end if

call dict%ptr%get_entry(index, array3)
if (allocated(array3)) then
   array_out(:size(array3)) = &
      & reshape(array3, [size(array3)])
   return
end if

end subroutine

subroutine get_array_size_index_api(verror, vdict, index, dim1, dim2, dim3) &
      & bind(C, name = namespace//"get_array_size_index")
   type(vp_double_dictionary), pointer :: dict
   type(c_ptr), value :: vdict
   type(c_ptr), value :: verror 
   type(vp_error), pointer :: error
   integer(kind=c_int), intent(in) :: index
   integer(kind=c_int), intent(out) :: dim1, dim2, dim3
   logical :: ok
   real(kind=wp), allocatable :: array1(:), array2(:, :), array3(:, :, :)
   if (debug) print '("[Info]", 1x, a)', "get array size by index"
   dim1 = 0
   dim2 = 0
   dim3 = 0

   call check_dict(verror, vdict, error, dict, ok)
   if (.not.(ok)) return


   call dict%ptr%get_entry(index, array1)
   if (allocated(array1)) then
      dim1 = size(array1, dim=1)
      return
   end if

   call dict%ptr%get_entry(index, array2)
   if (allocated(array2)) then
      dim1 = size(array2, dim=1)
      dim2 = size(array2, dim=2)
      return
   end if

   call dict%ptr%get_entry(index, array3)
   if (allocated(array3)) then
      dim1 = size(array3, dim=1)
      dim2 = size(array3, dim=2)
      dim3 = size(array3, dim=3)
      return
   end if

   call fatal_error(error%ptr, "Index is not a key in the dictionary.")
   return

end subroutine

subroutine get_array_entry_index_api(verror, vdict, index, array_out) &
      & bind(C, name = namespace//"get_array_entry_index")
   type(vp_double_dictionary), pointer :: dict
   type(c_ptr), value :: vdict
   type(c_ptr), value :: verror 
   type(vp_error), pointer :: error
   integer(kind=c_int), intent(in) :: index
   real(c_double), intent(out) :: array_out(*)
   real(kind=wp), allocatable :: array1(:), array2(:, :), array3(:, :, :)
   logical :: ok
   if (debug) print '("[Info]", 1x, a)', "get array value by index"
   call check_dict(verror, vdict, error, dict, ok)
   if (.not.(ok)) return

   call dict%ptr%get_entry(index, array1)
   if (allocated(array1)) then
      array_out(:size(array1)) = &
         & reshape(array1, [size(array1)])
      return
   end if

   call dict%ptr%get_entry(index, array2)
   if (allocated(array2)) then
      array_out(:size(array2)) = &
         & reshape(array2, [size(array2)])
      return
   end if

   call dict%ptr%get_entry(index, array3)
   if (allocated(array3)) then
      array_out(:size(array3)) = &
         & reshape(array3, [size(array3)])
      return
   end if

   call fatal_error(error%ptr, "Index is not a key in the dictionary.")
   return

end subroutine

subroutine get_label_entry_index_api(verror, vdict, index, charptr, buffersize) &
      & bind(C, name= namespace//"get_label_entry_index")
   type(vp_double_dictionary), pointer :: dict
   type(c_ptr), value :: vdict
   type(c_ptr), value :: verror 
   type(vp_error), pointer :: error
   integer(kind=c_int), intent(in) :: index
   character(kind=c_char), intent(out) :: charptr(*)
   integer(c_int), intent(in), optional :: buffersize
   character(len=:), allocatable :: fchar
   integer :: max_length
   logical :: ok
   if (debug) print '("[Info]", 1x, a)', "get dict label from index"
   call check_dict(verror, vdict, error, dict, ok)
   if (.not.(ok)) return

   associate(dict => dict%ptr)
      call dict%get_label(index, fchar)
   end associate
   
   if (present(buffersize)) then
      max_length = buffersize
   else
      max_length = huge(max_length) - 2
   end if

   call f_c_character(fchar, charptr, max_length)

end subroutine

subroutine delete_post_processing_api(vdict) &
      & bind(C, name=namespace//"delete_double_dictionary")
   type(c_ptr), intent(inout) :: vdict
   type(vp_double_dictionary), pointer :: dict
if (debug) print '("[Info]", 1x, a)', "get n_entries form dict"
   if (debug) print '("[Info]", 1x, a)', "delete_double_dictionary"

   if (c_associated(vdict)) then
      call c_f_pointer(vdict, dict)

      deallocate(dict)
      vdict = c_null_ptr
   end if
end subroutine delete_post_processing_api

end module
