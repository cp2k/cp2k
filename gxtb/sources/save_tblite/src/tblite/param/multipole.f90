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

!> @file tblite/param/multipole.f90
!> Provides a model for the anisotropic second-order electrostatic

!> Definition of the anisotropic second-order electrostatic contributions
module tblite_param_multipole
   use mctc_env, only : wp, error_type, fatal_error
   use tblite_coulomb_multipole, only : multipole_damping
   use tblite_param_serde, only : serde_record
   use tblite_toml, only : toml_table, get_value, set_value, add_table
   implicit none
   private

   public :: count


   character(len=*), parameter :: k_gfn2 = "gfn2", k_gxtb = "gxtb", &
      & k_kdmp3 = "kdmp3", k_kdmp5 = "kdmp5", k_kdmp7 = "kdmp7", k_kdmp9 = "kdmp9", &
      & k_kexp3 = "kexp3", k_kexp5 = "kexp5", k_kexp7 = "kexp7", k_kexp9 = "kexp9", &
      & k_kradexp = "kradexp", k_shift = "shift", k_rmax = "rmax", k_average = "average"

   !> Representation of the multipolar electrostatics
   type, public, extends(serde_record) :: multipole_record
      !> Damping type for the multipole electrostatics
      integer :: damping
      !> Damping exponent for quadratic terms
      real(wp) :: kdmp3
      !> Damping exponent for cubic terms
      real(wp) :: kdmp5
      !> Damping exponent for quartic terms
      real(wp) :: kdmp7
      !> Damping exponent for quintic terms
      real(wp) :: kdmp9
      !> Damping function exponent for quadratic terms
      real(wp) :: kexp3
      !> Damping function exponent for cubic terms
      real(wp) :: kexp5
      !> Damping function exponent for quartic terms
      real(wp) :: kexp7
      !> Damping function exponent for quintic terms
      real(wp) :: kexp9
      !> Exponent for multipole radii
      real(wp) :: kradexp
      !> Shift for valence CN
      real(wp) :: shift
      !> Maximum multipole radius
      real(wp) :: rmax
      !> Averaging scheme for the multipole radii
      character(len=:), allocatable :: average
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


   !> Masking for the anisotropic electrostatic parametrization
   type, public :: multipole_mask
   end type multipole_mask


   interface count
      module procedure :: count_mask
   end interface count


contains


!> Read parametrization data from TOML data structure
subroutine load_from_toml(self, table, error)
   !> Instance of the parametrization data
   class(multipole_record), intent(inout) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child
   integer :: stat

   self%damping = multipole_damping%gfn2
   call get_value(table, k_gfn2, child, requested=.false.)
   if (.not.associated(child)) then
      self%damping = multipole_damping%gxtb
      call get_value(table, k_gxtb, child, requested=.false.)
   end if
   if (.not.associated(child)) then
      call fatal_error(error, "No entry for damped multipole electrostatics found")
      return
   end if

   ! Read all general damping parameters
   call get_value(child, k_kdmp3, self%kdmp3, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for quadratic multipole damping prefactor")
      return
   end if
   call get_value(child, k_kdmp5, self%kdmp5, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for cubic multipole damping prefactor")
      return
   end if
   call get_value(child, k_kdmp7, self%kdmp7, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for quartic multipole damping prefactor")
      return
   end if
   call get_value(child, k_kdmp9, self%kdmp9, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for quintic multipole damping prefactor")
      return
   end if

   call get_value(child, k_kexp3, self%kexp3, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for quadratic multipole damping function exponent")
      return
   end if
   call get_value(child, k_kexp5, self%kexp5, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for cubic multipole damping function exponent")
      return
   end if
   call get_value(child, k_kexp7, self%kexp7, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for quartic multipole damping function exponent")
      return
   end if
   call get_value(child, k_kexp9, self%kexp9, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for quintic multipole damping function exponent")
      return
   end if

   ! Read GFN2 specific damping parameters
   if (self%damping == multipole_damping%gfn2) then
      call get_value(child, k_kradexp, self%kradexp, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Invalid entry for multipole damping function exponent")
         return
      end if
      call get_value(child, k_shift, self%shift, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Invalid entry for CN multipole shift")
         return
      end if
      call get_value(child, k_rmax, self%rmax, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Invalid entry for maximum multipole radius")
         return
      end if
   end if

   call get_value(child, k_average, self%average, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for multipole radii averaging")
      return
   end if
   select case(self%average)
   case default
      call fatal_error(error, "Invalid '"//self%average//"' averaging for multipole radii")
      return
   case("harmonic", "geometric", "arithmetic", "general")
   end select

end subroutine load_from_toml


!> Write parametrization data to TOML datastructure
subroutine dump_to_toml(self, table, error)
   !> Instance of the parametrization data
   class(multipole_record), intent(in) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child

   if (self%damping == multipole_damping%gfn2) then
      call add_table(table, k_gfn2, child)
      call set_value(child, k_kradexp, self%kradexp)
      call set_value(child, k_shift, self%shift)
      call set_value(child, k_rmax, self%rmax)
   else if (self%damping == multipole_damping%gxtb) then
      call add_table(table, k_gxtb, child)
   end if

   call set_value(child, k_kdmp3, self%kdmp3)
   call set_value(child, k_kdmp5, self%kdmp5)
   call set_value(child, k_kdmp7, self%kdmp7)
   call set_value(child, k_kdmp9, self%kdmp9)
   call set_value(child, k_kexp3, self%kexp3)
   call set_value(child, k_kexp5, self%kexp5)
   call set_value(child, k_kexp7, self%kexp7)
   call set_value(child, k_kexp9, self%kexp9)
   call set_value(child, k_average, self%average)

end subroutine dump_to_toml


!> Read parametrization data from parameter array
subroutine load_from_array(self, array, offset, base, mask, error)
   class(multipole_record), intent(inout) :: self
   real(wp), intent(in) :: array(:)
   integer, intent(inout) :: offset
   type(multipole_record), intent(in) :: base
   type(multipole_mask), intent(in) :: mask
   type(error_type), allocatable, intent(out) :: error

   select type(self)
   type is (multipole_record)
      self = base
   end select

end subroutine load_from_array

!> Write parametrization data to parameter array
subroutine dump_to_array(self, array, offset, mask, error)
   class(multipole_record), intent(in) :: self
   real(wp), intent(inout) :: array(:)
   integer, intent(inout) :: offset
   type(multipole_mask), intent(in) :: mask
   type(error_type), allocatable, intent(out) :: error

end subroutine dump_to_array

elemental function count_mask(mask) result(ncount)
   type(multipole_mask), intent(in) :: mask
   integer :: ncount
   ncount = 0
end function count_mask


end module tblite_param_multipole
