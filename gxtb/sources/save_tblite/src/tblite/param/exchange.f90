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

!> @file tblite/param/exchange.f90
!> Provides parameter record for the exchange

!> Defines model for the exchange parameters
module tblite_param_exchange
   use mctc_env, only : wp, error_type, fatal_error
   use tblite_toml, only : toml_table, get_value, set_value, add_table
   use tblite_param_serde, only : serde_record
   implicit none
   private

   public :: count

   character(len=*), parameter :: k_mulliken = "mulliken", k_frscale = "frscale", &
      & k_omega = "omega", k_lrscale = "lrscale", &
      & k_hubbard_average = "hubbard_average", k_gexp = "gexp", k_hubbard_exp = "hubbard_exp", &
      & k_hubbard_exp_r0 = "hubbard_exp_r0", k_ondiag = "ondiag", k_offdiag = "offdiag", &
      & k_ang(0:4) = ["s", "p", "d", "f", "g"], k_offdiag_average="offdiag_average", &
      & k_onsite = "onsite", k_kq = "kq", k_bond_order = "bond_order", &
      & k_cscale_average = "cscale_average", k_cexp = "cexp", k_crad_average = "crad_average"

   !> Parametrization record for the exchange interaction
   type, public, extends(serde_record) :: exchange_record
      !> Full-range scaling for the Fock exchange
      real(wp) :: frscale 
      !> Range separation parameter
      real(wp) :: omega
      !> Scaling of the long range exchange
      real(wp) :: lrscale
      !> Averaging function for the chemical hardness / Hubbard parameters
      character(len=:), allocatable :: hubbard_average
      !> Smoothening exponent: (1) Mataga-Nishimoto, (2) Klopman-Ohno
      real(wp) :: gexp
      !> Diagonal Fock exchange scaling 
      real(wp) :: ondiag
      !> Shell dependence of off-diagonal Fock exchange scaling
      real(wp) :: offdiag_l(0:4)
      !> Averaging function for the off-diagonal Fock exchange scaling
      character(len=:), allocatable :: offdiag_average
      !> maximum angular momentum supported
      integer :: lmax
      !> Exponent of radius dependent hubbard scaling
      real(wp) :: hubbard_exp
      !> Radius prefactor of radius dependent hubbard scaling
      real(wp) :: hubbard_exp_r0
      !> Shell charge-dependence of the onsite Fock exchange
      real(wp) :: kq(0:4)
      !> Averaging function for correlation scaling
      character(len=:), allocatable :: cscale_average
      !> Bond-order correlation damping exponent
      real(wp) :: cexp
      !> Averaging function for correlation radius
      character(len=:), allocatable :: crad_average
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

   !> Masking for the exchange model
   type, public :: exchange_mask
   end type exchange_mask

   interface count
      module procedure :: count_mask
   end interface count

contains


!> Read parametrization data from TOML data structure
subroutine load_from_toml(self, table, error)
   !> Instance of the parametrization data
   class(exchange_record), intent(inout) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child, child_shell
   real(wp), allocatable :: last
   integer :: l, stat

   ! Read global parameters
   call get_value(table, k_frscale, self%frscale, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for fullrange scale of exchange")
      return
   end if

   if (table%has_key(k_omega)) then
      call get_value(table, k_omega, self%omega, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Invalid entry for range separation parameter")
         return
      end if

      call get_value(table, k_lrscale, self%lrscale, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Invalid entry for long range scaling")
         return
      end if
   else
      self%omega = 0.0_wp
      self%lrscale = 0.0_wp
   endif

   ! Read parameters for Mulliken-approximated exchange
   call get_value(table, k_mulliken, child, requested=.false.)
   if (.not.associated(child)) then
      call fatal_error(error, "No entry for Mulliken-approximated Fock exchange found")
      return
   end if

   call get_value(child, k_hubbard_average, self%hubbard_average, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for hardness averaging")
      return
   end if
   select case(self%hubbard_average)
   case default
      call fatal_error(error, "Invalid '"//self%hubbard_average//"' averaging for hardness")
      return
   case("harmonic", "geometric", "arithmetic", "general")
   end select

   call get_value(child, k_gexp, self%gexp, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for effective interaction exponent")
      return
   end if

   call get_value(child, k_ondiag, self%ondiag, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for diagonal Fock exchange scaling")
      return
   end if

   call get_value(child, k_offdiag, child_shell, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for off-diagonal scaling")
      return
   else
      do l = 0, 4
         if (.not. child_shell%has_key(k_ang(l))) then
            if (allocated(last)) then
               self%offdiag_l(l) = last
               cycle
            end if
            call fatal_error(error, "No entry for "//k_ang(l)//"-shell off-diagonal scaling")
            exit
         end if
         call get_value(child_shell, k_ang(l), self%offdiag_l(l), stat=stat)
         if (stat /= 0) then
            call fatal_error(error, "Cannot read "//k_ang(l)//"-shell off-diagonal scaling")
            exit
         end if
         if (stat == 0) then
            last = self%offdiag_l(l)
            self%lmax = l
         end if
      end do
      if (allocated(error)) return
   end if

   call get_value(child, k_offdiag_average, self%offdiag_average, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for off-diagonal scaling averaging")
      return
   end if
   select case(self%offdiag_average)
   case default
      call fatal_error(error, "Invalid '"//self%offdiag_average//"' averaging for off-diagonal scaling")
      return
   case("harmonic", "geometric", "arithmetic", "general")
   end select

   call get_value(child, k_hubbard_exp, self%hubbard_exp, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for radius dependent hubbard scaling")
      return
   end if

   call get_value(child, k_hubbard_exp_r0, self%hubbard_exp_r0, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for radius prefactor of radius dependent hubbard scaling")
      return
   end if

   ! Read parameters for onsite exchange correction
   call get_value(table, k_onsite, child, requested=.false.)
   if (.not.associated(child)) then
      call fatal_error(error, "No entry for onsite Fock exchange found")
      return
   end if

   call get_value(child, k_kq, child_shell, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for charge-dependence of onsite Fock exchange")
      return
   else
      do l = 0, 4
         if (.not. child_shell%has_key(k_ang(l))) then
            if (allocated(last)) then
               self%kq(l) = last
               cycle
            end if
            call fatal_error(error, "No entry for "//k_ang(l)//"-shell charge-dependence of onsite Fock exchange")
            exit
         end if
         call get_value(child_shell, k_ang(l), self%kq(l), stat=stat)
         if (stat /= 0) then
            call fatal_error(error, "Cannot read "//k_ang(l)//"-shell charge-dependence of onsite Fock exchange")
            exit
         end if
         if (stat == 0) then
            last = self%kq(l)
            self%lmax = l
         end if
      end do
      if (allocated(error)) return
   end if

   ! Read parameters for bond-order correlation correction
   call get_value(table, k_bond_order, child, requested=.false.)
   if (.not.associated(child)) then
      call fatal_error(error, "No entry for bond-order correlation correction found")
      return
   end if
   call get_value(child, k_cscale_average, self%cscale_average, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for correlation scaling averaging")
      return
   end if
   select case(self%cscale_average)
   case default
      call fatal_error(error, "Invalid '"//self%cscale_average//"' averaging for correlation scaling")
      return
   case("harmonic", "geometric", "arithmetic", "general")
   end select

   call get_value(child, k_cexp, self%cexp, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for bond-order correlation damping exponent")
      return
   end if

   call get_value(child, k_crad_average, self%crad_average, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for correlation radius averaging")
      return
   end if
   select case(self%crad_average)
   case default
      call fatal_error(error, "Invalid '"//self%crad_average//"' averaging for correlation radius")
      return
   case("harmonic", "geometric", "arithmetic", "general")
   end select

end subroutine load_from_toml


!> Write parametrization data to TOML datastructure
subroutine dump_to_toml(self, table, error)
   !> Instance of the parametrization data
   class(exchange_record), intent(in) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child, child_shell
   integer :: l

   call set_value(table, k_frscale, self%frscale)
   if (self%omega .ne. 0.0_wp) then
      call set_value(table, k_omega, self%omega)
      call set_value(table, k_lrscale, self%lrscale)
   end if

   call add_table(table, k_mulliken, child)
   call set_value(child, k_hubbard_average, self%hubbard_average)
   call set_value(child, k_gexp, self%gexp)
   
   call set_value(child, k_ondiag, self%ondiag)
   call add_table(child, k_offdiag, child_shell)
   do l = 0, self%lmax
      call set_value(child_shell, k_ang(l), self%offdiag_l(l))
   end do
   call set_value(child, k_offdiag_average, self%offdiag_average)

   call set_value(child, k_hubbard_exp, self%hubbard_exp)
   call set_value(child, k_hubbard_exp_r0, self%hubbard_exp_r0)

   call add_table(table, k_onsite, child)
   call add_table(child, k_kq, child_shell)
   do l = 0, self%lmax
      call set_value(child_shell, k_ang(l), self%kq(l))
   end do

   call add_table(table, k_bond_order, child)
   call set_value(child, k_cscale_average, self%cscale_average)
   call set_value(child, k_cexp, self%cexp)
   call set_value(child, k_crad_average, self%crad_average)

end subroutine dump_to_toml


!> Read parametrization data from parameter array
subroutine load_from_array(self, array, offset, base, mask, error)
   class(exchange_record), intent(inout) :: self
   real(wp), intent(in) :: array(:)
   integer, intent(inout) :: offset
   type(exchange_record), intent(in) :: base
   type(exchange_mask), intent(in) :: mask
   type(error_type), allocatable, intent(out) :: error

   select type(self)
   type is (exchange_record)
      self = base
   end select

end subroutine load_from_array

!> Write parametrization data to parameter array
subroutine dump_to_array(self, array, offset, mask, error)
   class(exchange_record), intent(in) :: self
   real(wp), intent(inout) :: array(:)
   integer, intent(inout) :: offset
   type(exchange_mask), intent(in) :: mask
   type(error_type), allocatable, intent(out) :: error

end subroutine dump_to_array

elemental function count_mask(mask) result(ncount)
   type(exchange_mask), intent(in) :: mask
   integer :: ncount
   ncount = 0
end function count_mask

end module tblite_param_exchange

