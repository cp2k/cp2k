! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! Lesser GNU General Public License for more details.
!
! You should have received a copy of the Lesser GNU General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

!> @dir dftd4/param
!> Contains the parameterizations for the DFT-D models

!> @file dftd4/param.f90
!> Reexports of the parameter types and parameterizations.

!> Proxy module for handling the parametrizations
module dftd4_param
   use dftd4_model_type, only : dftd_models
   use dftd4_param_d4, only : get_damping_d4
   use dftd4_param_functionals, only : functional_group, get_functionals, &
      & get_functional_id, p_invalid, p_r2scan_3c, p_default
   use dftd4_param_type, only : param_type
   use dftd4_utils, only : lowercase
   use mctc_env, only : error_type, fatal_error, wp
   implicit none
   private

   public :: param_type, get_damping_params
   public :: functional_group, get_functionals
   public :: get_functional_id, p_invalid, p_r2scan_3c, p_default

   !> Retrieve damping parameters from functional name or ID
   interface get_damping_params
      module procedure :: get_damping_params_name
      module procedure :: get_damping_params_id
   end interface get_damping_params
   !DEC$ ATTRIBUTES DLLEXPORT :: get_damping_params

contains

!> Retrieve damping parameters from functional name
subroutine get_damping_params_name(error, functional, model, damping_2b, &
   & damping_3b, param)
   !DEC$ ATTRIBUTES DLLEXPORT :: get_damping_params_name

   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Functional name for which to retrieve the damping parameters
   character(len=*), intent(in) :: functional
   !> Type of dispersion model to use
   integer, intent(in) :: model
   !> Type of two-body damping function to use
   integer, intent(in) :: damping_2b
   !> Type of three-body damping function to use
   integer, intent(in) :: damping_3b
   !> Damping parameters for the functional
   type(param_type), intent(inout) :: param

   character(len=:), allocatable :: fname
   integer :: is, id

   is = index(functional, '/')
   if (is == 0) is = len_trim(functional) + 1
   fname = lowercase(functional(:is-1))

   id = get_functional_id(fname)

   call get_damping_params_id(error, id, model, damping_2b, damping_3b, param)

end subroutine get_damping_params_name

subroutine get_damping_params_id(error, id, model, damping_2b, damping_3b, param)
   !DEC$ ATTRIBUTES DLLEXPORT :: get_damping_params_id

   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Functional ID for which to retrieve the damping parameters
   integer, intent(in) :: id
   !> Type of dispersion model to use
   integer, intent(in) :: model
   !> Type of two-body damping function to use
   integer, intent(in) :: damping_2b
   !> Type of three-body damping function to use
   integer, intent(in) :: damping_3b
   !> Damping parameters for the functional
   type(param_type), intent(inout) :: param

   select case(model)
   case(dftd_models%d4, dftd_models%d4s)
      call get_damping_d4(error, id, damping_2b, damping_3b, param)
   case default
      call fatal_error(error, "No damping parameters available for this dispersion model.")
      return
   end select

end subroutine get_damping_params_id


end module dftd4_param
