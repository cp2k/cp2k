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

!> @file tblite/param/repulsion.f90
!> Provides a module for the repulsion interactions

!> Definition of the repulsion interactions
module tblite_param_repulsion
   use mctc_env, only : wp, error_type, fatal_error
   use tblite_param_serde, only : serde_record
   use tblite_repulsion_type, only : repulsion_kernel
   use tblite_toml, only : toml_table, get_value, set_value, add_table
   implicit none
   private

   public :: count


   character(len=*), parameter :: k_gfn = "gfn", k_gxtb = "gxtb", k_kexp = "kexp", &
      & k_klight = "klight", k_rexp = "rexp", k_k2="k2", k_k2light="k2light", &
      & k_k3="k3", k_k4="k4", k_short="short", k_short_alpha="short_alpha", &
      & k_short_exp="short_exp"

   !> Parametrization records describing the repulsion interactions
   type, public, extends(serde_record) :: repulsion_record
      !> Repulsion kernel
      integer :: kernel
      !> Exponent of the repulsion polynomial
      real(wp) :: rexp
      !> Distance exponent for repulsion damping
      real(wp) :: kexp
      !> Distance exponent for repulsion damping of light atom pairs
      real(wp) :: klight
      !> Second-order expansion coefficient for the repulsion
      !> for elements other than H or He
      real(wp) :: k2
      !> Second-order expansion coefficient for the repulsion 
      !> for light elements (H or He)
      real(wp) :: k2light
      !> Third-order expansion coefficient for the repulsion
      real(wp) :: k3
      !> Fourth-order expansion coefficient for the repulsion
      real(wp) :: k4
      !> Short-range repulsion correction prefactor
      real(wp) :: short
      !> Short-range repulsion correction damping exponent scaling
      real(wp) :: short_alpha
      !> Short-range repulsion correction distance damping exponent 
      real(wp) :: short_exp
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


   !> Provides a mask for the repulsion model
   type, public :: repulsion_mask
   end type repulsion_mask


   interface count
      module procedure :: count_mask
   end interface count


contains


!> Read parametrization data from TOML data structure
subroutine load_from_toml(self, table, error)
   !> Instance of the parametrization data
   class(repulsion_record), intent(inout) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child
   integer :: stat

   self%kernel = repulsion_kernel%gfn
   call get_value(table, k_gfn, child, requested=.false.)
   if (.not.associated(child)) then
      self%kernel = repulsion_kernel%gxtb
      call get_value(table, k_gxtb, child, requested=.false.)
   end if
   if (.not.associated(child)) then
      call fatal_error(error, "No entry for repulsion model found")
      return
   end if

   ! Read the GFN type repulsion model 
   if (self%kernel == repulsion_kernel%gfn) then

      call get_value(child, k_kexp, self%kexp, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Invalid entry for repulsion exponent")
         return
      end if

      call get_value(child, k_klight, self%klight, self%kexp, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Invalid entry for light-atom repulsion exponent")
         return
      end if
      
      call get_value(child, k_rexp, self%rexp, 1.0_wp, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Invalid entry for the repulsion polynomial exponent")
         return
      end if

   ! Read the g-xTB type repulsion model
   else

      call get_value(child, k_kexp, self%kexp, 1.0_wp, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Invalid entry for the repulsion exponent")
         return
      end if

      call get_value(child, k_k2, self%k2, 1.0_wp, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Invalid entry for the second order repulsion expansion coefficient")
         return
      end if

      call get_value(child, k_k2light, self%k2light, 1.0_wp, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Invalid entry for the light-atom second order repulsion expansion coefficient")
         return
      end if

      call get_value(child, k_k3, self%k3, 1.0_wp, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Invalid entry for the third order repulsion expansion coefficient")
         return
      end if

      call get_value(child, k_k4, self%k4, 1.0_wp, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Invalid entry for the fourth order repulsion expansion coefficient")
         return
      end if

      call get_value(child, k_short, self%short, 0.0_wp, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Invalid entry for the short-range repulsion correction prefactor")
         return
      end if

      call get_value(child, k_short_alpha, self%short_alpha, 1.0_wp, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Invalid entry for the short-range repulsion correction damping exponent scaling")
         return
      end if

      call get_value(child, k_short_exp, self%short_exp, 2.0_wp, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Invalid entry for the short-range repulsion correction distance damping exponent")
         return
      end if

   end if

end subroutine load_from_toml


!> Write parametrization data to TOML datastructure
subroutine dump_to_toml(self, table, error)
   !> Instance of the parametrization data
   class(repulsion_record), intent(in) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child

   if (self%kernel == repulsion_kernel%gfn) then
      call add_table(table, k_gfn, child)
      call set_value(child, k_kexp, self%kexp)
      if (abs(self%kexp - self%klight) > epsilon(self%kexp)) then
         call set_value(child, k_klight, self%klight)
      end if
      call set_value(child, k_rexp, self%rexp)
   else
      call add_table(table, k_gxtb, child)
      call set_value(child, k_kexp, self%kexp)
      call set_value(child, k_k2, self%k2)
      call set_value(child, k_k2light, self%k2light)
      call set_value(child, k_k3, self%k3)
      call set_value(child, k_k4, self%k4)
      call set_value(child, k_short, self%short)
      call set_value(child, k_short_alpha, self%short_alpha)
      call set_value(child, k_short_exp, self%short_exp)
   end if

end subroutine dump_to_toml


!> Read parametrization data from parameter array
subroutine load_from_array(self, array, offset, base, mask, error)
   class(repulsion_record), intent(inout) :: self
   real(wp), intent(in) :: array(:)
   integer, intent(inout) :: offset
   type(repulsion_record), intent(in) :: base
   type(repulsion_mask), intent(in) :: mask
   type(error_type), allocatable, intent(out) :: error

   select type(self)
   type is (repulsion_record)
      self = base
   end select

end subroutine load_from_array

!> Write parametrization data to parameter array
subroutine dump_to_array(self, array, offset, mask, error)
   class(repulsion_record), intent(in) :: self
   real(wp), intent(inout) :: array(:)
   integer, intent(inout) :: offset
   type(repulsion_mask), intent(in) :: mask
   type(error_type), allocatable, intent(out) :: error

end subroutine dump_to_array

elemental function count_mask(mask) result(ncount)
   type(repulsion_mask), intent(in) :: mask
   integer :: ncount
   ncount = 0
end function count_mask


end module tblite_param_repulsion
