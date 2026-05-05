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

!> @file tblite/param/firstorder.f90
!> Provides a model for the isotropic first-order electrostatic

!> Definition of the isotropic first-order electrostatic contributions.
module tblite_param_firstorder
   use mctc_env, only : wp, error_type, fatal_error
   use tblite_param_serde, only : serde_record
   use tblite_toml, only : toml_table, get_value, set_value, add_table
   implicit none
   private

   public :: count


   character(len=*), parameter :: k_split = "split", k_split_exp = "split_exp", &
      & k_split_slope = "split_slope", k_split_offset = "split_offset"

   !> Parametrization record for first-order electrostatic contributions
   type, public, extends(serde_record) :: firstorder_record
      !> Exponent of IP and EA spliting function
      real(wp) :: split_exp
      !> Slope of IP and EA spliting function
      real(wp) :: split_slope
      !> Offset of IP and EA spliting function
      real(wp) :: split_offset
      !> Maximum angular momentum for shell-specific CN-dependence of IP/EA
      integer :: lmax
   contains
      generic :: load => load_from_array
      generic :: dump => dump_to_array
      !> Read parametrization data from TOML data structure
      procedure :: load_from_toml
      !> Write parametrization data to TOML data structure
      procedure :: dump_to_toml
      !> Read parametrization data from parameter array
      procedure, private :: load_from_array
      !> Write parametrization data to parameter array
      procedure, private :: dump_to_array
   end type


   !> Masking for the first order electrostatics
   type, public :: firstorder_mask
   end type firstorder_mask


   interface count
      module procedure :: count_mask
   end interface count


contains


!> Read parametrization data from TOML data structure
subroutine load_from_toml(self, table, error)
   !> Instance of the parametrization data
   class(firstorder_record), intent(inout) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child, child_shell
   real(wp), allocatable :: last
   integer :: l, stat

   call get_value(table, k_split, child, requested=.false.)
   if (.not.associated(child)) then
      call fatal_error(error, "No entry for first-order IP/EA splitting found")
      return
   end if

   call get_value(child, k_split_exp, self%split_exp, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for first-order IP/EA splitting exponent")
      return
   end if
   call get_value(child, k_split_slope, self%split_slope, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for first-order IP/EA splitting slope")
      return
   end if
   call get_value(child, k_split_offset, self%split_offset, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for first-order IP/EA splitting offset")
      return
   end if

end subroutine load_from_toml


!> Write parametrization data to TOML datastructure
subroutine dump_to_toml(self, table, error)
   !> Instance of the parametrization data
   class(firstorder_record), intent(in) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child, child_shell
   integer :: l

   call add_table(table, k_split, child)
   call set_value(child, k_split_exp, self%split_exp)
   call set_value(child, k_split_slope, self%split_slope)
   call set_value(child, k_split_offset, self%split_offset)
   
end subroutine dump_to_toml


!> Read parametrization data from parameter array
subroutine load_from_array(self, array, offset, base, mask, error)
   class(firstorder_record), intent(inout) :: self
   real(wp), intent(in) :: array(:)
   integer, intent(inout) :: offset
   type(firstorder_record), intent(in) :: base
   type(firstorder_mask), intent(in) :: mask
   type(error_type), allocatable, intent(out) :: error

   select type(self)
   type is (firstorder_record)
      self = base
   end select

end subroutine load_from_array

!> Write parametrization data to parameter array
subroutine dump_to_array(self, array, offset, mask, error)
   class(firstorder_record), intent(in) :: self
   real(wp), intent(inout) :: array(:)
   integer, intent(inout) :: offset
   type(firstorder_mask), intent(in) :: mask
   type(error_type), allocatable, intent(out) :: error

end subroutine dump_to_array

elemental function count_mask(mask) result(ncount)
   type(firstorder_mask), intent(in) :: mask
   integer :: ncount
   ncount = 0
end function count_mask


end module tblite_param_firstorder
