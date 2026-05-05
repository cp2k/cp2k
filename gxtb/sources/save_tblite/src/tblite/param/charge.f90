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

!> @file tblite/param/charge.f90
!> Provides model for the isotropic second-order electrostatic

!> Definition of the isotropic second-order electrostatic model
module tblite_param_charge
   use mctc_env, only : wp, error_type, fatal_error
   use tblite_coulomb_charge, only : coulomb_kernel
   use tblite_param_serde, only : serde_record
   use tblite_toml, only : toml_table, get_value, set_value, add_table
   implicit none
   private

   public :: count


   character(len=*), parameter :: k_effective = "effective", k_gexp = "gexp", &
      & k_average = "average", k_gamma = "gamma", k_hubbard_exp = "hubbard_exp"
   real(wp), parameter :: default_hubbard_exp = 0.0_wp

   !> Parametrization record for the isotropic second-order electrostatics
   type, public, extends(serde_record) :: charge_record
      !> Coulomb interaction kernel
      integer :: kernel
      !> Averaging scheme for the chemical hardness / Hubbard parameters
      character(len=:), allocatable :: average
      !> Exponent manipulating the long range behaviour of the Coulombic kernel
      real(wp) :: gexp
      !> Exponent of radius dependent hubbard scaling
      real(wp) :: hubbard_exp
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


   !> Masking for the isotropic electrostatic model
   type, public :: charge_mask
   end type charge_mask


   interface count
      module procedure :: count_mask
   end interface count


contains


!> Read parametrization data from TOML data structure
subroutine load_from_toml(self, table, error)
   !> Instance of the parametrization data
   class(charge_record), intent(inout) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child
   integer :: stat

   self%kernel = coulomb_kernel%dftb_gamma
   call get_value(table, k_gamma, child, requested=.false.)
   if (.not.associated(child)) then
      self%kernel = coulomb_kernel%effective
      call get_value(table, k_effective, child, requested=.false.)
   end if
   if (.not.associated(child)) then
      call fatal_error(error, "No entry for Coulomb electrostatic found")
      return
   end if

   if (self%kernel == coulomb_kernel%effective) then
      call get_value(child, k_gexp, self%gexp, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Invalid entry for effective Coulomb exponent")
         return
      end if

      call get_value(child, k_average, self%average, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Invalid entry for hardness averaging")
         return
      end if
      select case(self%average)
      case default
         call fatal_error(error, "Invalid '"//self%average//"' averaging for hardness")
         return
      case("harmonic", "geometric", "arithmetic", "general")
      end select

      call get_value(child, k_hubbard_exp, self%hubbard_exp, default_hubbard_exp, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Invalid entry for radius dependent hubbard scaling")
         return
      end if
   end if
end subroutine load_from_toml


!> Write parametrization data to TOML datastructure
subroutine dump_to_toml(self, table, error)
   !> Instance of the parametrization data
   class(charge_record), intent(in) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child

   if (self%kernel == coulomb_kernel%effective) then
      call add_table(table, k_effective, child)
      call set_value(child, k_gexp, self%gexp)
      call set_value(child, k_average, self%average)
      call set_value(child, k_hubbard_exp, self%hubbard_exp)
   else
      call add_table(table, k_gamma, child)
   end if
end subroutine dump_to_toml


!> Read parametrization data from parameter array
subroutine load_from_array(self, array, offset, base, mask, error)
   class(charge_record), intent(inout) :: self
   real(wp), intent(in) :: array(:)
   integer, intent(inout) :: offset
   type(charge_record), intent(in) :: base
   type(charge_mask), intent(in) :: mask
   type(error_type), allocatable, intent(out) :: error

   select type(self)
   type is (charge_record)
      self = base
   end select

end subroutine load_from_array

!> Write parametrization data to parameter array
subroutine dump_to_array(self, array, offset, mask, error)
   class(charge_record), intent(in) :: self
   real(wp), intent(inout) :: array(:)
   integer, intent(inout) :: offset
   type(charge_mask), intent(in) :: mask
   type(error_type), allocatable, intent(out) :: error

end subroutine dump_to_array

elemental function count_mask(mask) result(ncount)
   type(charge_mask), intent(in) :: mask
   integer :: ncount
   ncount = 0
end function count_mask

end module tblite_param_charge
