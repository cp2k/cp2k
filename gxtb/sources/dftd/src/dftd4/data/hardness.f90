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

module dftd4_data_hardness
   use mctc_env, only : wp
   use mctc_io_symbols, only : to_number
   implicit none
   private

   public :: get_hardness


   interface get_hardness
      module procedure :: get_hardness_num
      module procedure :: get_hardness_sym
   end interface get_hardness


   integer, parameter :: max_elem = 118

  !> Element-specific chemical hardnesses for the charge scaling function used
  !> to extrapolate the C6 coefficients in DFT-D4.
  real(wp), parameter :: chemical_hardness(max_elem) = [ &
    & 0.47259288_wp, 0.92203391_wp, 0.17452888_wp, 0.25700733_wp, 0.33949086_wp, &
    & 0.42195412_wp, 0.50438193_wp, 0.58691863_wp, 0.66931351_wp, 0.75191607_wp, &
    & 0.17964105_wp, 0.22157276_wp, 0.26348578_wp, 0.30539645_wp, 0.34734014_wp, &
    & 0.38924725_wp, 0.43115670_wp, 0.47308269_wp, 0.17105469_wp, 0.20276244_wp, &
    & 0.21007322_wp, 0.21739647_wp, 0.22471039_wp, 0.23201501_wp, 0.23933969_wp, &
    & 0.24665638_wp, 0.25398255_wp, 0.26128863_wp, 0.26859476_wp, 0.27592565_wp, &
    & 0.30762999_wp, 0.33931580_wp, 0.37235985_wp, 0.40273549_wp, 0.43445776_wp, &
    & 0.46611708_wp, 0.15585079_wp, 0.18649324_wp, 0.19356210_wp, 0.20063311_wp, &
    & 0.20770522_wp, 0.21477254_wp, 0.22184614_wp, 0.22891872_wp, 0.23598621_wp, &
    & 0.24305612_wp, 0.25013018_wp, 0.25719937_wp, 0.28784780_wp, 0.31848673_wp, &
    & 0.34912431_wp, 0.37976593_wp, 0.41040808_wp, 0.44105777_wp, 0.05019332_wp, &
    & 0.06762570_wp, 0.08504445_wp, 0.10247736_wp, 0.11991105_wp, 0.13732772_wp, &
    & 0.15476297_wp, 0.17218265_wp, 0.18961288_wp, 0.20704760_wp, 0.22446752_wp, &
    & 0.24189645_wp, 0.25932503_wp, 0.27676094_wp, 0.29418231_wp, 0.31159587_wp, &
    & 0.32902274_wp, 0.34592298_wp, 0.36388048_wp, 0.38130586_wp, 0.39877476_wp, &
    & 0.41614298_wp, 0.43364510_wp, 0.45104014_wp, 0.46848986_wp, 0.48584550_wp, &
    & 0.12526730_wp, 0.14268677_wp, 0.16011615_wp, 0.17755889_wp, 0.19497557_wp, & ! Tl-At
    & 0.21240778_wp, 0.07263525_wp, 0.09422158_wp, 0.09920295_wp, 0.10418621_wp, & ! Rn-Th
    & 0.14235633_wp, 0.16394294_wp, 0.18551941_wp, 0.22370139_wp, 0.25110000_wp, & ! Pa-Am 
    & 0.25030000_wp, 0.28840000_wp, 0.31000000_wp, 0.33160000_wp, 0.35320000_wp, & ! Cm-Fm
    & 0.36820000_wp, 0.39630000_wp, 0.40140000_wp, 0.00000000_wp, 0.00000000_wp, & ! Md-Db
    & 0.00000000_wp, 0.00000000_wp, 0.00000000_wp, 0.00000000_wp, 0.00000000_wp, & ! Sg-Ds
    & 0.00000000_wp, 0.00000000_wp, 0.00000000_wp, 0.00000000_wp, 0.00000000_wp, & ! Rg-Mc
    & 0.00000000_wp, 0.00000000_wp, 0.00000000_wp] ! Lv,Ts,Og 

contains


!> Get chemical hardness for a given element symbol
elemental function get_hardness_sym(sym) result(eta)

   !> Element symbol
   character(len=*), intent(in) :: sym

   !> Chemical hardness
   real(wp) :: eta

   eta = get_hardness(to_number(sym))

end function get_hardness_sym


!> Get chemical hardness for a given atomic number
elemental function get_hardness_num(num) result(eta)

   !> Atomic number
   integer, intent(in) :: num

   !> Chemical hardness
   real(wp) :: eta

   if (num > 0 .and. num <= size(chemical_hardness)) then
      eta = chemical_hardness(num)
   else
      eta = 0.0_wp
   end if

end function get_hardness_num


end module dftd4_data_hardness
