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

!> @file tblite/param/thirdorder.f90
!> Provides a model for the isotropic third-order tight-binding

!> Definition of the isotropic third-order tight-binding contributions.
module tblite_param_thirdorder
   use mctc_env, only : wp, error_type, fatal_error
   use tblite_coulomb_thirdorder, only : thirdorder_kernel
   use tblite_param_serde, only : serde_record
   use tblite_toml, only : toml_table, get_value, set_value, add_table
   implicit none
   private

   public :: count


   character(len=*), parameter :: k_onsite = "onsite", k_twobody="twobody", &
      & k_shell = "shell", k_ang(0:4) = ["s", "p", "d", "f", "g"], &
      & k_texp = "texp", k_average = "average", k_onsite_scale = "onsite_scale", &
      & k_offsite_scale = "offsite_scale"

   !> Parametrization record for third-order tight-binding contributions
   type, public, extends(serde_record) :: thirdorder_record
      !> Third-order interaction kernel
      integer :: kernel   
      !> Maximum angular momentum for shell-specific Hubbard derivatives
      integer :: lmax
      !> Whether the third-order contribution is shell-resolved
      logical :: shell
      !> Scaling factor for shell-specific Hubbard derivatives
      real(wp) :: ksh(0:4)
      !> Averaging scheme for the chemical hardness / Hubbard parameters
      character(len=:), allocatable :: average
      !> Exponent for two-body third-order interaction
      real(wp) :: texp
      !> Onsite scaling factor for two-body interaction
      real(wp) :: onsite_scale
      !> Offsite scaling factor for two-body interaction
      real(wp) :: offsite_scale
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


   !> Masking for the third order tight-binding
   type, public :: thirdorder_mask
   end type thirdorder_mask


   interface count
      module procedure :: count_mask
   end interface count


contains


!> Read parametrization data from TOML data structure
subroutine load_from_toml(self, table, error)
   !> Instance of the parametrization data
   class(thirdorder_record), intent(inout) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child, child_shell
   real(wp), allocatable :: last
   integer :: l, stat

   ! Read onsite third-order kernel parameters
   self%kernel = thirdorder_kernel%onsite
   call get_value(table, k_onsite, child, requested=.false.)
   if (.not.associated(child)) then
      self%kernel = thirdorder_kernel%twobody
      call get_value(table, k_twobody, child, requested=.false.)
   end if
   if (.not.associated(child)) then
      call fatal_error(error, "No entry for third-order tight-binding found")
      return
   end if

   ! Read two-body third-order kernel parameters
   if (self%kernel == thirdorder_kernel%twobody) then
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

      call get_value(child, k_texp, self%texp, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Invalid entry for two-body third-order exponent")
         return
      end if
      call get_value(child, k_onsite_scale, self%onsite_scale, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Invalid entry for onsite scaling of two-body third-order term")
         return
      end if
      call get_value(child, k_offsite_scale, self%offsite_scale, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Invalid entry for offsite scaling of two-body third-order term")
         return
      end if
   end if

   ! Read general third-order shell-dependence parameters
   call get_value(child, k_shell, child_shell, requested=.false.)
   self%shell = associated(child_shell)
   if (self%shell) then
      do l = 0, 4
         if (.not.child_shell%has_key(k_ang(l))) then
            if (allocated(last)) then
               self%ksh(l) = last
               cycle
            end if
            call fatal_error(error, "No entry for "//k_ang(l)//"-shell Hubbard derivative")
            exit
         end if
         call get_value(child_shell, k_ang(l), self%ksh(l), stat=stat)
         if (stat /= 0) then
            call fatal_error(error, "Cannot read "//k_ang(l)//"-shell Hubbard derivative")
            exit
         end if
         if (stat == 0) then
            last = self%ksh(l)
            self%lmax = l
         end if
      end do
      if (allocated(error)) return
   end if
end subroutine load_from_toml


!> Write parametrization data to TOML datastructure
subroutine dump_to_toml(self, table, error)
   !> Instance of the parametrization data
   class(thirdorder_record), intent(in) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child, child_shell
   integer :: l

   if (self%kernel == thirdorder_kernel%onsite) then
      call add_table(table, k_onsite, child)
      if (self%shell) then
         call add_table(child, k_shell, child_shell)
         do l = 0, self%lmax
            call set_value(child_shell, k_ang(l), self%ksh(l))
         end do
      else
         call set_value(child, k_shell, .false.)
      end if
   else
      call add_table(table, k_twobody, child)
      if (self%shell) then
         call add_table(child, k_shell, child_shell)
         do l = 0, self%lmax
            call set_value(child_shell, k_ang(l), self%ksh(l))
         end do
      else
         call set_value(child, k_shell, .false.)
      end if
      call set_value(child, k_average, self%average)
      call set_value(child, k_texp, self%texp)
      call set_value(child, k_onsite_scale, self%onsite_scale)
      call set_value(child, k_offsite_scale, self%offsite_scale)
   end if
end subroutine dump_to_toml


!> Read parametrization data from parameter array
subroutine load_from_array(self, array, offset, base, mask, error)
   class(thirdorder_record), intent(inout) :: self
   real(wp), intent(in) :: array(:)
   integer, intent(inout) :: offset
   type(thirdorder_record), intent(in) :: base
   type(thirdorder_mask), intent(in) :: mask
   type(error_type), allocatable, intent(out) :: error

   select type(self)
   type is (thirdorder_record)
      self = base
   end select

end subroutine load_from_array

!> Write parametrization data to parameter array
subroutine dump_to_array(self, array, offset, mask, error)
   class(thirdorder_record), intent(in) :: self
   real(wp), intent(inout) :: array(:)
   integer, intent(inout) :: offset
   type(thirdorder_mask), intent(in) :: mask
   type(error_type), allocatable, intent(out) :: error

end subroutine dump_to_array

elemental function count_mask(mask) result(ncount)
   type(thirdorder_mask), intent(in) :: mask
   integer :: ncount
   ncount = 0
end function count_mask


end module tblite_param_thirdorder
