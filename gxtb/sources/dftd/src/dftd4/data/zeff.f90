! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

module dftd4_data_zeff
   use mctc_env, only : wp
   use mctc_io_symbols, only : to_number
   implicit none
   private

   public :: get_effective_charge


   interface get_effective_charge
      module procedure :: get_effective_charge_num
      module procedure :: get_effective_charge_sym
   end interface get_effective_charge


   integer, parameter :: max_elem = 118


  !> Effective nuclear charges from the def2-ECPs used for calculating the 
  !> reference polarizabilities for DFT-D4.
  real(wp), parameter :: effective_nuclear_charge(max_elem) = [ &
    &   1,                                                 2,  & ! H-He
    &   3, 4,                               5, 6, 7, 8, 9,10,  & ! Li-Ne
    &  11,12,                              13,14,15,16,17,18,  & ! Na-Ar
    &  19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,  & ! K-Kr
    &   9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,  & ! Rb-Xe
    &   9,10,11,30,31,32,33,34,35,36,37,38,39,40,41,42,43,  & ! Cs-Lu
    &  12,13,14,15,16,17,18,19,20,21,22,23,24,25,26, & ! Hf-Rn
    !  just copy & paste from above
    &   9,10,11,30,31,32,33,34,35,36,37,38,39,40,41,42,43,  & ! Fr-Lr
    &  12,13,14,15,16,17,18,19,20,21,22,23,24,25,26] ! Rf-Og

contains


!> Get effective nuclear charge for a given element symbol
elemental function get_effective_charge_sym(sym) result(zeff)

   !> Element symbol
   character(len=*), intent(in) :: sym

   !> Effective nuclear charge
   real(wp) :: zeff

   zeff = get_effective_charge(to_number(sym))

end function get_effective_charge_sym


!> Get effective nuclear charge for a given atomic number
elemental function get_effective_charge_num(num) result(zeff)

   !> Atomic number
   integer, intent(in) :: num

   !> Effective nuclear charge
   real(wp) :: zeff

   if (num > 0 .and. num <= size(effective_nuclear_charge)) then
      zeff = effective_nuclear_charge(num)
   else
      zeff = 0.0_wp
   end if

end function get_effective_charge_num


end module dftd4_data_zeff
