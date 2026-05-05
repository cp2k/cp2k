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

!> @file tblite/api/param.f90
!> Provides API exports for the #tblite_param handle.

!> API export for managing parametrization records.
module tblite_api_param
   use, intrinsic :: iso_c_binding
   use mctc_env, only : fatal_error
   use tblite_api_error, only : vp_error
   use tblite_api_table, only : vp_table
   use tblite_api_version, only : namespace
   use tblite_api_utils, only : f_c_character
   use tblite_param, only : param_record
   use tblite_toml, only : toml_table
   use tblite_xtb_gfn1, only : export_gfn1_param
   use tblite_xtb_gfn2, only : export_gfn2_param
   use tblite_xtb_ipea1, only : export_ipea1_param
   implicit none
   private

   public :: vp_param
   public :: new_param_api, delete_param_api
   public :: load_param_api, dump_param_api
   public :: export_gfn1_param_api, export_gfn2_param_api, export_ipea1_param_api

   !> Void pointer to manage parametrization records
   type :: vp_param
      !> Actual payload
      type(param_record) :: ptr
   end type vp_param

   logical, parameter :: debug = .false.

contains

!> Create new parametrization records object
function new_param_api() &
      & result(vparam) &
      & bind(C, name=namespace//"new_param")
   type(vp_param), pointer :: param
   type(c_ptr) :: vparam

   if (debug) print '("[Info]", 1x, a)', "new_param"

   allocate(param)
   vparam = c_loc(param)
end function new_param_api

!> Delete paramrization records object
subroutine delete_param_api(vparam) &
      & bind(C, name=namespace//"delete_param")
   type(c_ptr), intent(inout) :: vparam
   type(vp_param), pointer :: param

   if (debug) print '("[Info]", 1x, a)', "delete_param"

   if (c_associated(vparam)) then
      call c_f_pointer(vparam, param)

      deallocate(param)
      vparam = c_null_ptr
   end if

end subroutine delete_param_api

!> Load parametrization record from data table
subroutine load_param_api(verror, vparam, vtable) &
      & bind(C, name=namespace//"load_param")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vparam
   type(vp_param), pointer :: param
   type(c_ptr), value :: vtable
   type(vp_table), pointer :: table

   if (debug) print '("[Info]", 1x, a)', "load_param"

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vparam)) then
      call fatal_error(error%ptr, "Parametrization record is missing")
      return
   end if
   call c_f_pointer(vparam, param)

   if (.not.c_associated(vtable)) then
      call fatal_error(error%ptr, "Data table object is missing")
      return
   end if
   call c_f_pointer(vtable, table)

   call param%ptr%load(table%ptr, error%ptr)

end subroutine load_param_api

!> Dump parametrization record to data table
subroutine dump_param_api(verror, vparam, vtable) &
      & bind(C, name=namespace//"dump_param")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vparam
   type(vp_param), pointer :: param
   type(c_ptr), value :: vtable
   type(vp_table), pointer :: table

   if (debug) print '("[Info]", 1x, a)', "dump_param"

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vparam)) then
      call fatal_error(error%ptr, "Parametrization record is missing")
      return
   end if
   call c_f_pointer(vparam, param)

   if (.not.c_associated(vtable)) then
      call fatal_error(error%ptr, "Data table object is missing")
      return
   end if
   call c_f_pointer(vtable, table)

   call param%ptr%dump(table%ptr, error%ptr)

end subroutine dump_param_api


subroutine export_gfn1_param_api(verror, vparam) &
      & bind(C, name=namespace//"export_gfn1_param")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vparam
   type(vp_param), pointer :: param

   if (debug) print '("[Info]", 1x, a)', "export_gfn1_param"

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vparam)) then
      call fatal_error(error%ptr, "Parametrization record is missing")
      return
   end if
   call c_f_pointer(vparam, param)

   call export_gfn1_param(param%ptr)
end subroutine export_gfn1_param_api


subroutine export_gfn2_param_api(verror, vparam) &
      & bind(C, name=namespace//"export_gfn2_param")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vparam
   type(vp_param), pointer :: param

   if (debug) print '("[Info]", 1x, a)', "export_gfn2_param"

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vparam)) then
      call fatal_error(error%ptr, "Parametrization record is missing")
      return
   end if
   call c_f_pointer(vparam, param)

   call export_gfn2_param(param%ptr)
end subroutine export_gfn2_param_api


subroutine export_ipea1_param_api(verror, vparam) &
      & bind(C, name=namespace//"export_ipea1_param")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vparam
   type(vp_param), pointer :: param

   if (debug) print '("[Info]", 1x, a)', "ipea1_param"

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vparam)) then
      call fatal_error(error%ptr, "Parametrization record is missing")
      return
   end if
   call c_f_pointer(vparam, param)

   call export_ipea1_param(param%ptr)
end subroutine export_ipea1_param_api

end module tblite_api_param
