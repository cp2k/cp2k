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

!> @dir tblite/fit
!> Contains routines for parameter optimization

!> @file tblite/fit/settings.f90
!> Provides setting for the optimization driver

!> Declaration of options for the optimization driver
module tblite_fit_settings
   use mctc_env, only : wp, error_type, fatal_error
   use tblite_param_serde, only : serde_record
   use tblite_param, only : param_record, param_mask
   use tblite_toml, only : toml_table, get_value, set_value, add_table
   implicit none
   private

   !> Optimization configuration
   type, public, extends(serde_record) :: fit_settings
      character(len=:), allocatable :: method
      character(len=:), allocatable :: script
      character(len=:), allocatable :: output
      character(len=:), allocatable :: fitpar
      logical :: relative
      integer :: max_iter
      real(wp) :: trustr
      type(param_record), allocatable :: base
      type(param_mask) :: mask
   contains
      !> Read parametrization data from TOML data structure
      procedure :: load_from_toml
      !> Write parametrization data to TOML data structure
      procedure :: dump_to_toml
   end type fit_settings

contains

!> Read parametrization data from TOML data structure
subroutine load_from_toml(self, table, error)
   !> Instance of the parametrization data
   class(fit_settings), intent(inout) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child

   call get_value(table, "max-iter", self%max_iter, 100)
   if (self%max_iter <= 0) then
      call fatal_error(error, "Positive number of iterations required")
      return
   end if

   call get_value(table, "method", self%method, "newuoa")
   if (self%method /= "newuoa") then
      call fatal_error(error, "Unknown optimization method '"//self%method//"'")
      return
   end if

   call get_value(table, "relative", self%relative, .true.)

   call get_value(table, "trust-rad", self%trustr, 0.1_wp)
   if (self%trustr <= 0.0_wp) then
      call fatal_error(error, "Positive value for trust radius required")
      return
   end if

   call get_value(table, "script", self%script, "./run.sh")
   call get_value(table, "data-file", self%output, ".data")
   call get_value(table, "param-file", self%fitpar, "fitpar.toml")

   call set_base_param(self%mask, self%base)
   call get_value(table, "mask", child)
   call self%mask%load(child, error)
   if (allocated(error)) return
end subroutine load_from_toml

subroutine set_base_param(mask, base)
   type(param_mask), intent(out) :: mask
   type(param_record), intent(in), target, optional :: base
   if (present(base)) mask%ref => base%record
end subroutine set_base_param

!> Write parametrization data to TOML datastructure
subroutine dump_to_toml(self, table, error)
   !> Instance of the parametrization data
   class(fit_settings), intent(in) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child

   call set_value(table, "max-iter", self%max_iter)
   call set_value(table, "method", self%method)
   call set_value(table, "relative", self%relative)
   call set_value(table, "trust-rad", self%trustr)
   call set_value(table, "script", self%script)
   call set_value(table, "data-file", self%output)
   call set_value(table, "param-file", self%fitpar)

   call add_table(table, "mask", child)
   call self%mask%dump(child, error)
   if (allocated(error)) return
end subroutine dump_to_toml

end module tblite_fit_settings
