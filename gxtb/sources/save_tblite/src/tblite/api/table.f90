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

!> @file tblite/api/table.f90
!> Provides API exports for the #tblite_table handle.

!> API export for managing data tables
module tblite_api_table
   use, intrinsic :: iso_c_binding
   use mctc_env, only : fatal_error
   use tblite_api_error, only : vp_error
   use tblite_api_version, only : namespace
   use tblite_api_utils, only : c_f_character, strlen
   use tblite_toml, only : toml_table, toml_array, toml_value, add_array, set_value, get_value
   implicit none
   private

   public :: vp_table
   public :: new_table_api, delete_table_api
   public :: table_set_double_api, table_set_int64_t_api, table_set_bool_api, &
      & table_set_char_api, table_add_table_api

   !> Void pointer to manage general data tables
   type :: vp_table
      !> Actual payload
      type(toml_table), pointer :: ptr
      !> Data is owned by the object
      logical :: owned
   end type vp_table

   logical, parameter :: debug = .false.

contains

!> Create data table reference object
function new_table_api(vtable) &
      & result(vval) &
      & bind(C, name=namespace//"new_table")
   type(c_ptr), value :: vtable
   type(vp_table), pointer :: table
   type(vp_table), pointer :: val
   type(c_ptr) :: vval
   type(toml_table), pointer :: dat

   if (debug) print '("[Info]", 1x, a)', "new_table"

   allocate(val)
   if (c_associated(vtable)) then
      call c_f_pointer(vtable, table)
      val%ptr => table%ptr
      val%owned = .false.
   else
      allocate(dat)
      dat = toml_table()
      val%ptr => dat
      val%owned = .true.
   end if
   vval = c_loc(val)
end function new_table_api

!> Delete data table object
subroutine delete_table_api(vtable) &
      & bind(C, name=namespace//"delete_table")
   type(c_ptr), intent(inout) :: vtable
   type(vp_table), pointer :: table

   if (debug) print '("[Info]", 1x, a)', "delete_table"

   if (c_associated(vtable)) then
      call c_f_pointer(vtable, table)

      if (table%owned) deallocate(table%ptr)
      deallocate(table)
      vtable = c_null_ptr
   end if

end subroutine delete_table_api

subroutine table_set_double_api(verror, vtable, ckey, val, n) &
      & bind(C, name=namespace//"table_set_double")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vtable
   type(vp_table), pointer :: table
   character(kind=c_char), intent(in) :: ckey(*)
   character(len=:), allocatable :: key
   real(c_double), intent(in) :: val(*)
   integer(c_int), value :: n
   type(toml_array), pointer :: array
   integer :: i, stat

   if (debug) print '("[Info]", 1x, a)', "table_set_double"

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vtable)) then
      call fatal_error(error%ptr, "Data table object is missing")
      return
   end if

   call c_f_pointer(vtable, table)
   call c_f_character(ckey, key)

   if (table%ptr%has_key(key)) call table%ptr%delete(key)

   if (n > 0) then
      call add_array(table%ptr, key, array)
      do i = 1, n
         call set_value(array, i, val(i), stat=stat)
         if (stat /= 0) exit
      end do
   else
      call set_value(table%ptr, key, val(1), stat=stat)
   end if

   if (stat /= 0) then
      call fatal_error(error%ptr, "Failed to push back double value(s) to data table")
   end if
end subroutine table_set_double_api

subroutine table_set_int64_t_api(verror, vtable, ckey, val, n) &
      & bind(C, name=namespace//"table_set_int64_t")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vtable
   type(vp_table), pointer :: table
   character(kind=c_char), intent(in) :: ckey(*)
   character(len=:), allocatable :: key
   integer(c_int64_t), intent(in) :: val(*)
   integer(c_int), value :: n
   type(toml_array), pointer :: array
   integer :: i, stat

   if (debug) print '("[Info]", 1x, a)', "table_set_int64_t"

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vtable)) then
      call fatal_error(error%ptr, "Data table object is missing")
      return
   end if

   call c_f_pointer(vtable, table)
   call c_f_character(ckey, key)

   if (table%ptr%has_key(key)) call table%ptr%delete(key)

   if (n > 0) then
      call add_array(table%ptr, key, array)
      do i = 1, n
         call set_value(array, i, val(i), stat=stat)
         if (stat /= 0) exit
      end do
   else
      call set_value(table%ptr, key, val(1), stat=stat)
   end if

   if (stat /= 0) then
      call fatal_error(error%ptr, "Failed to push back integer value(s) to data table")
   end if
end subroutine table_set_int64_t_api

subroutine table_set_bool_api(verror, vtable, ckey, val, n) &
      & bind(C, name=namespace//"table_set_bool")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vtable
   type(vp_table), pointer :: table
   character(kind=c_char), intent(in) :: ckey(*)
   character(len=:), allocatable :: key
   logical(c_bool), intent(in) :: val(*)
   integer(c_int), value :: n
   type(toml_array), pointer :: array
   integer :: i, stat

   if (debug) print '("[Info]", 1x, a)', "table_set_bool"

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vtable)) then
      call fatal_error(error%ptr, "Data table object is missing")
      return
   end if

   call c_f_pointer(vtable, table)
   call c_f_character(ckey, key)

   if (table%ptr%has_key(key)) call table%ptr%delete(key)

   if (n > 0) then
      call add_array(table%ptr, key, array)
      do i = 1, n
         call set_value(array, i, merge(.true., .false., val(i)), stat=stat)
         if (stat /= 0) exit
      end do
   else
      call set_value(table%ptr, key, merge(.true., .false., val(1)), stat=stat)
   end if

   if (stat /= 0) then
      call fatal_error(error%ptr, "Failed to push back boolean value(s) to data table")
   end if
end subroutine table_set_bool_api

subroutine table_set_char_api(verror, vtable, ckey, cval, n) &
      & bind(C, name=namespace//"table_set_char")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vtable
   type(vp_table), pointer :: table
   character(kind=c_char), intent(in) :: ckey(*)
   character(len=:), allocatable :: key, val
   type(c_ptr), value :: cval
   character(kind=c_char), pointer :: carr(:, :)
   integer(c_int), value :: n
   type(toml_array), pointer :: array
   integer :: i, stat

   if (debug) print '("[Info]", 1x, a)', "table_set_char"

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vtable)) then
      call fatal_error(error%ptr, "Data table object is missing")
      return
   end if

   call c_f_pointer(vtable, table)
   call c_f_character(ckey, key)
   call c_f_pointer(cval, carr, [strlen(cval)+1, min(n, 1)])

   if (table%ptr%has_key(key)) call table%ptr%delete(key)

   if (n > 0) then
      call add_array(table%ptr, key, array)
      do i = 1, n
         call c_f_character(carr(:, i), val)
         call set_value(array, i, val, stat=stat)
         if (stat /= 0) exit
      end do
   else
      call c_f_character(carr(:, 1), val)
      call set_value(table%ptr, key, val, stat=stat)
   end if

   if (stat /= 0) then
      call fatal_error(error%ptr, "Failed to push back character string to data table")
   end if
end subroutine table_set_char_api

function table_add_table_api(verror, vtable, ckey) &
      & result(vval) &
      & bind(C, name=namespace//"table_add_table")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vtable
   type(vp_table), pointer :: table
   character(kind=c_char), intent(in) :: ckey(*)
   character(len=:), allocatable :: key
   type(c_ptr) :: vval
   type(vp_table), pointer :: val
   type(toml_table), pointer :: tmp
   integer :: stat

   if (debug) print '("[Info]", 1x, a)', "table_add_table"

   vval = c_null_ptr
   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vtable)) then
      call fatal_error(error%ptr, "Data table object is missing")
      return
   end if

   vval = new_table_api(vtable)
   call c_f_pointer(vtable, table)
   call c_f_pointer(vval, val)
   call c_f_character(ckey, key)

   call get_value(table%ptr, key, tmp, stat=stat)
   val%ptr => tmp

   if (stat /= 0) then
      call fatal_error(error%ptr, "Failed to push back subtable to data table")
      call delete_table_api(vval)
   end if
end function table_add_table_api


end module tblite_api_table
