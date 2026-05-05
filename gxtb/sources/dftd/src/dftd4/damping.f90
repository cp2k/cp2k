! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
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

!> @dir dftd4/damping
!> Contains the damping functions for the DFT-D models

!> @file dftd4/damping.f90
!> Reexports of the damping functions. 

!> Proxy module for handling the damping functions
module dftd4_damping
   use dftd4_damping_cso, only : cso_damping_twobody, new_cso_damping_twobody
   use dftd4_damping_koide, only : koide_damping_twobody, new_koide_damping_twobody
   use dftd4_damping_mzero, only : mzero_damping_twobody, new_mzero_damping_twobody
   use dftd4_damping_optpower, only : optpower_damping_twobody, new_optpower_damping_twobody
   use dftd4_damping_rational, only : rational_damping_twobody, rational_damping_threebody, &
      & new_rational_damping_twobody, new_rational_damping_threebody
   use dftd4_damping_screened, only : screened_damping_twobody, screened_damping_threebody, &
      & new_screened_damping_twobody, new_screened_damping_threebody
   use dftd4_damping_type, only : damping_type, damping_twobody, damping_threebody, &
      & twobody_damping_function, threebody_damping_function
   use dftd4_damping_zero, only : zero_damping_twobody, zero_damping_threebody, &
      & zero_avg_damping_threebody, new_zero_damping_twobody, new_zero_damping_threebody, &
      & new_zero_avg_damping_threebody
   use dftd4_utils, only : lowercase
   use mctc_env, only : error_type, fatal_error, wp
   implicit none
   private

   public :: damping_type, damping_twobody, damping_threebody, new_damping
   public :: twobody_damping_function, threebody_damping_function
   public :: cso_damping_twobody
   public :: koide_damping_twobody
   public :: mzero_damping_twobody
   public :: optpower_damping_twobody
   public :: rational_damping_twobody, rational_damping_threebody
   public :: screened_damping_twobody, screened_damping_threebody
   public :: zero_damping_twobody, zero_damping_threebody, zero_avg_damping_threebody
   public :: get_damping_function_id

   !> Damping function identifiers from the function names
   interface get_damping_function_id
      module procedure :: get_damping_function_id_model
      module procedure :: get_damping_function_id_name
   end interface get_damping_function_id

contains

subroutine new_damping(error, damp, damping_2b, damping_3b)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Damping object
   class(damping_type), intent(out) :: damp
   !> Type of two-body damping function to use
   integer, intent(in) :: damping_2b
   !> Type of three-body damping function to use
   integer, intent(in) :: damping_3b

   ! Select the two-body damping function
   select case(damping_2b)
   case(twobody_damping_function%rational)
      block
         type(rational_damping_twobody), allocatable :: tmp
         allocate(tmp)
         call new_rational_damping_twobody(tmp)
         call move_alloc(tmp, damp%damping_2b)
      end block
   case(twobody_damping_function%screened)
      block
         type(screened_damping_twobody), allocatable :: tmp
         allocate(tmp)
         call new_screened_damping_twobody(tmp)
         call move_alloc(tmp, damp%damping_2b)
      end block
   case(twobody_damping_function%zero)
      block
         type(zero_damping_twobody), allocatable :: tmp
         allocate(tmp)
         call new_zero_damping_twobody(tmp)
         call move_alloc(tmp, damp%damping_2b)
      end block
   case(twobody_damping_function%mzero)
      block
         type(mzero_damping_twobody), allocatable :: tmp
         allocate(tmp)
         call new_mzero_damping_twobody(tmp)
         call move_alloc(tmp, damp%damping_2b)
      end block
   case(twobody_damping_function%optpower)
      block
         type(optpower_damping_twobody), allocatable :: tmp
         allocate(tmp)
         call new_optpower_damping_twobody(tmp)
         call move_alloc(tmp, damp%damping_2b)
      end block
   case(twobody_damping_function%cso)
      block
         type(cso_damping_twobody), allocatable :: tmp
         allocate(tmp)
         call new_cso_damping_twobody(tmp)
         call move_alloc(tmp, damp%damping_2b)
      end block
   case(twobody_damping_function%koide)
      block
         type(koide_damping_twobody), allocatable :: tmp
         allocate(tmp)
         call new_koide_damping_twobody(tmp)
         call move_alloc(tmp, damp%damping_2b)
      end block
   case default
      call fatal_error(error, "Unsupported option for the two-body damping function.")
      return
   end select

   ! Select the optional three-body damping function
   select case(damping_3b)
   case(threebody_damping_function%rational)
      block
         type(rational_damping_threebody), allocatable :: tmp
         allocate(tmp)
         call new_rational_damping_threebody(tmp)
         call move_alloc(tmp, damp%damping_3b)
      end block
   case(threebody_damping_function%screened)
      block
         type(screened_damping_threebody), allocatable :: tmp
         allocate(tmp)
         call new_screened_damping_threebody(tmp)
         call move_alloc(tmp, damp%damping_3b)
      end block
   case(threebody_damping_function%zero)
      block
         type(zero_damping_threebody), allocatable :: tmp
         allocate(tmp)
         call new_zero_damping_threebody(tmp)
         call move_alloc(tmp, damp%damping_3b)
      end block
   case(threebody_damping_function%zero_avg)
      block
         type(zero_avg_damping_threebody), allocatable :: tmp
         allocate(tmp)
         call new_zero_avg_damping_threebody(tmp)
         call move_alloc(tmp, damp%damping_3b)
      end block
   end select

end subroutine new_damping


!> Return the damping function identifiers for a given damping function object
subroutine get_damping_function_id_model(error, damp, damping_2b_id, damping_3b_id)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Damping function
   type(damping_type), intent(in) :: damp

   !> Identifier of the two-body damping function
   integer, intent(out), optional :: damping_2b_id

   !> Identifier of the three-body damping function
   integer, intent(out), optional :: damping_3b_id

   if (allocated(damp%damping_2b) .and. present(damping_2b_id)) then
      associate(damping_2b => damp%damping_2b)
         select type(damping_2b)
         type is (rational_damping_twobody)
            damping_2b_id = twobody_damping_function%rational
         type is (screened_damping_twobody)
            damping_2b_id = twobody_damping_function%screened
         type is (zero_damping_twobody)
            damping_2b_id = twobody_damping_function%zero
         type is (mzero_damping_twobody)
            damping_2b_id = twobody_damping_function%mzero
         type is (optpower_damping_twobody)
            damping_2b_id = twobody_damping_function%optpower
         type is (cso_damping_twobody)
            damping_2b_id = twobody_damping_function%cso
         type is (koide_damping_twobody)
            damping_2b_id = twobody_damping_function%koide
         class default
            call fatal_error(error, "Unknown two-body damping function type")
         end select
      end associate
   end if

   if (allocated(damp%damping_3b) .and. present(damping_3b_id)) then
      associate(damping_3b => damp%damping_3b)
         select type(damping_3b)
         type is (rational_damping_threebody)
            damping_3b_id = threebody_damping_function%rational
         type is (screened_damping_threebody)
            damping_3b_id = threebody_damping_function%screened
         type is (zero_damping_threebody)
            damping_3b_id = threebody_damping_function%zero
         type is (zero_avg_damping_threebody)
            damping_3b_id = threebody_damping_function%zero_avg
         class default
            call fatal_error(error, "Unknown three-body damping function type")
         end select
      end associate
   end if

end subroutine get_damping_function_id_model


!> Return the damping function identifiers for a case insensitive function names
subroutine get_damping_function_id_name(error, damping_2b_name, damping_3b_name, &
   & damping_2b_id, damping_3b_id)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Name of the two-body damping function
   character(len=*), intent(in), optional :: damping_2b_name

   !> Name of the three-body damping function
   character(len=*), intent(in), optional :: damping_3b_name

   !> Identifier of the two-body damping function
   integer, intent(out), optional :: damping_2b_id

   !> Identifier of the three-body damping function
   integer, intent(out), optional :: damping_3b_id

   if (present(damping_2b_name) .and. present(damping_2b_id)) then
      select case(lowercase(damping_2b_name))
      case("rational")
         damping_2b_id = twobody_damping_function%rational
      case("screened")
         damping_2b_id = twobody_damping_function%screened
      case("zero")
         damping_2b_id = twobody_damping_function%zero
      case("mzero")
         damping_2b_id = twobody_damping_function%mzero
      case("optpower")
         damping_2b_id = twobody_damping_function%optpower
      case("cso")
         damping_2b_id = twobody_damping_function%cso
      case("koide")
         damping_2b_id = twobody_damping_function%koide
      case default
         call fatal_error(error, "Unknown two-body damping function name")
      end select
   end if

   if (present(damping_3b_name) .and. present(damping_3b_id)) then
      select case(lowercase(damping_3b_name))
      case("rational")
         damping_3b_id = threebody_damping_function%rational
      case("screened")
         damping_3b_id = threebody_damping_function%screened
      case("zero")
         damping_3b_id = threebody_damping_function%zero
      case("zero_avg")
         damping_3b_id = threebody_damping_function%zero_avg
      case default
         call fatal_error(error, "Unknown three-body damping function name")
      end select
   end if

end subroutine get_damping_function_id_name

end module dftd4_damping
