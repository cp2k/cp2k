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

!> Implementation of the driver entry point for handling parametrization data
module tblite_driver_param
   use, intrinsic :: iso_fortran_env, only : output_unit
   use mctc_env, only : error_type, fatal_error, wp
   use tblite_cli, only : param_config
   use tblite_output_ascii
   use tblite_param, only : param_record
   use tblite_toml, only : toml_table, merge_table
   use tblite_xtb_gfn2, only : export_gfn2_param
   use tblite_xtb_gfn1, only : export_gfn1_param
   use tblite_xtb_gxtb, only : export_gxtb_param
   use tblite_xtb_ipea1, only : export_ipea1_param
   implicit none
   private

   public :: main

   interface main
      module procedure :: param_main
   end interface

contains


subroutine param_main(config, error)
   type(param_config), intent(in) :: config
   type(error_type), allocatable, intent(out) :: error

   type(param_record), allocatable :: param, base
   type(toml_table) :: table1, table2

   if (allocated(config%method)) then
      allocate(base)
      select case(config%method)
      case default
         call fatal_error(error, "Unknown method '"//config%method//"' requested")
      case("gfn2")
         call export_gfn2_param(base)
      case("gfn1")
         call export_gfn1_param(base)
      case("gxtb")
         call export_gxtb_param(base)
      case("ipea1")
         call export_ipea1_param(base)
      end select
   end if
   if (allocated(error)) return

   if (allocated(config%input)) then
      allocate(param)
      call param%load(config%input, error)

      if (.not.allocated(error) .and. allocated(base)) then
         table1 = toml_table()
         call param%dump(table1, error)
         if (.not. allocated(error)) then
            table2 = toml_table()
            call base%dump(table2, error)
         end if
         if (.not. allocated(error)) then
            call merge_table(table1, table2)
            call table2%destroy
            deallocate(param, base)
            allocate(param)
            call param%load(table1, error)
            call table1%destroy
         end if
      end if
   else
      call move_alloc(base, param)
   end if
   if (allocated(error)) return

   if (config%verbosity > 1) then
      if (allocated(param%name)) &
         write(output_unit, '(a)') param%name
      if (allocated(param%reference)) &
         write(output_unit, '(a)') param%reference
   end if

   if (allocated(config%output)) then
      call param%dump(config%output, error)
      if (.not.allocated(error)) then
         if (config%verbosity > 0) write(output_unit, '(a)') &
            "[Info] Parameter file written to '"//config%output//"'"
      end if
   end if
   if (allocated(error)) return

end subroutine param_main


end module tblite_driver_param
