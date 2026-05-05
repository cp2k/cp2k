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

module dftd4_data_en
   use mctc_env, only : wp
   use mctc_io_symbols, only : to_number
   implicit none
   private

   public :: get_electronegativity


   interface get_electronegativity
      module procedure :: get_electronegativity_num
      module procedure :: get_electronegativity_sym
   end interface get_electronegativity


   integer, parameter :: max_elem = 118

   !> Pauling electronegativities, used for the covalent coordination number.
   real(wp), parameter :: pauling_en(max_elem) = [ &
      & 2.20_wp,3.00_wp, & ! H,He
      & 0.98_wp,1.57_wp,2.04_wp,2.55_wp,3.04_wp,3.44_wp,3.98_wp,4.50_wp, & ! Li-Ne
      & 0.93_wp,1.31_wp,1.61_wp,1.90_wp,2.19_wp,2.58_wp,3.16_wp,3.50_wp, & ! Na-Ar
      & 0.82_wp,1.00_wp, & ! K,Ca
      &                 1.36_wp,1.54_wp,1.63_wp,1.66_wp,1.55_wp, & ! Sc-
      &                 1.83_wp,1.88_wp,1.91_wp,1.90_wp,1.65_wp, & ! -Zn
      &                 1.81_wp,2.01_wp,2.18_wp,2.55_wp,2.96_wp,3.00_wp, & ! Ga-Kr
      & 0.82_wp,0.95_wp, & ! Rb,Sr
      &                 1.22_wp,1.33_wp,1.60_wp,2.16_wp,1.90_wp, & ! Y-
      &                 2.20_wp,2.28_wp,2.20_wp,1.93_wp,1.69_wp, & ! -Cd
      &                 1.78_wp,1.96_wp,2.05_wp,2.10_wp,2.66_wp,2.60_wp, & ! In-Xe
      & 0.79_wp,0.89_wp, & ! Cs,Ba
      &         1.10_wp,1.12_wp,1.13_wp,1.14_wp,1.15_wp,1.17_wp,1.18_wp, & ! La-Eu
      &         1.20_wp,1.21_wp,1.22_wp,1.23_wp,1.24_wp,1.25_wp,1.26_wp, & ! Gd-Yb
      &                 1.27_wp,1.30_wp,1.50_wp,2.36_wp,1.90_wp, & ! Lu-
      &                 2.20_wp,2.20_wp,2.28_wp,2.54_wp,2.00_wp, & ! -Hg
      &                 1.62_wp,2.33_wp,2.02_wp,2.00_wp,2.20_wp,2.20_wp, & ! Tl-Rn

      & 0.79_wp,0.90_wp, & ! Fr,Ra
      &         1.10_wp,1.30_wp,1.50_wp,1.38_wp,1.36_wp,1.28_wp,1.30_wp, & ! Ac-Am
      &         1.30_wp,1.30_wp,1.30_wp,1.30_wp,1.30_wp,1.30_wp,1.30_wp, & ! Cm-No
      &         1.30_wp, & ! Lr
      ! only dummies below
      &                 1.50_wp,1.50_wp,1.50_wp,1.50_wp, & ! Rf-Bh
      &                 1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp, & ! Hs-Cn
      &                 1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp,1.50_wp ] ! Nh-Og

contains


!> Get electronegativity for a given element symbol
elemental function get_electronegativity_sym(sym) result(en)
   !DEC$ ATTRIBUTES DLLEXPORT :: get_electronegativity_sym

   !> Element symbol
   character(len=*), intent(in) :: sym

   !> Electronegativity
   real(wp) :: en

   en = get_electronegativity(to_number(sym))

end function get_electronegativity_sym


!> Get electronegativity for a given atomic number
elemental function get_electronegativity_num(num) result(en)
   !DEC$ ATTRIBUTES DLLEXPORT :: get_electronegativity_num

   !> Atomic number
   integer, intent(in) :: num

   !> Electronegativity
   real(wp) :: en

   if (num > 0 .and. num <= size(pauling_en)) then
      en = pauling_en(num)
   else
      en = 0.0_wp
   end if

end function get_electronegativity_num


end module dftd4_data_en
