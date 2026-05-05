! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
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

module dftd4_data_r4r2
   use mctc_env, only : wp
   use mctc_io_convert, only : aatoau
   use mctc_io_symbols, only : to_number
   implicit none
   private

   public :: get_r4r2_val


   !> Get r4/r2 expectation value
   interface get_r4r2_val
      module procedure :: get_r4r2_val_num
      module procedure :: get_r4r2_val_sym
   end interface get_r4r2_val


   integer, parameter :: max_elem = 118

   !  PBE0/def2-QZVP atomic values calculated by S. Grimme in Gaussian (2010)
   !  rare gases recalculated by J. Mewes with PBE0/aug-cc-pVQZ in Dirac (2018)
   !  He: 3.4698 -> 3.5544, Ne: 3.1036 -> 3.7943, Ar: 5.6004 -> 5.6638,
   !  Kr: 6.1971 -> 6.2312, Xe: 7.5152 -> 8.8367
   !  not replaced but recalculated (PBE0/cc-pVQZ) were
   !   H: 8.0589 ->10.9359, Li:29.0974 ->39.7226, Be:14.8517 ->17.7460
   !  also new super heavies Cn,Nh,Fl,Lv,Og
   !  Am-Rg calculated at 4c-PBE/Dyall-AE4Z (Dirac 2022)
   real(wp), parameter :: r4_over_r2(max_elem) = [  &
      &  8.0589_wp, 3.4698_wp, & ! H,He
      & 29.0974_wp,14.8517_wp,11.8799_wp, 7.8715_wp, 5.5588_wp, 4.7566_wp, 3.8025_wp, 3.1036_wp, & ! Li-Ne
      & 26.1552_wp,17.2304_wp,17.7210_wp,12.7442_wp, 9.5361_wp, 8.1652_wp, 6.7463_wp, 5.6004_wp, & ! Na-Ar
      & 29.2012_wp,22.3934_wp, & ! K,Ca
      &           19.0598_wp,16.8590_wp,15.4023_wp,12.5589_wp,13.4788_wp, & ! Sc-
      &           12.2309_wp,11.2809_wp,10.5569_wp,10.1428_wp, 9.4907_wp, & ! -Zn
      &                      13.4606_wp,10.8544_wp, 8.9386_wp, 8.1350_wp, 7.1251_wp, 6.1971_wp, & ! Ga-Kr
      & 30.0162_wp,24.4103_wp, & ! Rb,Sr
      &            20.3537_wp,17.4780_wp,13.5528_wp,11.8451_wp,11.0355_wp, & ! Y-
      &            10.1997_wp, 9.5414_wp, 9.0061_wp, 8.6417_wp, 8.9975_wp, & ! -Cd
      &                       14.0834_wp,11.8333_wp,10.0179_wp, 9.3844_wp, 8.4110_wp, 7.5152_wp, & ! In-Xe
      & 32.7622_wp,27.5708_wp, & ! Cs,Ba
      &            23.1671_wp,21.6003_wp,20.9615_wp,20.4562_wp,20.1010_wp,19.7475_wp,19.4828_wp, & ! La-Eu
      &            15.6013_wp,19.2362_wp,17.4717_wp,17.8321_wp,17.4237_wp,17.1954_wp,17.1631_wp, & ! Gd-Yb
      &            14.5716_wp,15.8758_wp,13.8989_wp,12.4834_wp,11.4421_wp, & ! Lu-
      &            10.2671_wp, 8.3549_wp, 7.8496_wp, 7.3278_wp, 7.4820_wp, & ! -Hg
      &                       13.5124_wp,11.6554_wp,10.0959_wp, 9.7340_wp, 8.8584_wp, 8.0125_wp, & ! Tl-Rn
      & 29.8135_wp,26.3157_wp, & ! Fr,Ra
      &            19.1885_wp,15.8542_wp,16.1305_wp,15.6161_wp,15.1226_wp,16.1576_wp,14.6510_wp, & ! Ac-Am
      &            14.7178_wp,13.9108_wp,13.5623_wp,13.2326_wp,12.9189_wp,12.6133_wp,12.3142_wp, & ! Cm-No
      &            14.8326_wp,12.3771_wp,10.6378_wp, 9.3638_wp, 8.2297_wp, & ! Lr-
      &             7.5667_wp, 6.9456_wp, 6.3946_wp, 5.9159_wp, 5.4929_wp, & ! -Cn
      &                        6.7286_wp, 6.5144_wp,10.9169_wp,10.3600_wp, 9.4723_wp, 8.6641_wp ] ! Nh-Og

   integer :: idum
   real(wp), parameter :: sqrt_z_r4_over_r2(max_elem) = &
      &  sqrt(0.5_wp*(r4_over_r2*[(sqrt(real(idum,wp)),idum=1,max_elem)]))


contains


!> Get r4/r2 expectation value for a given element symbol
elemental function get_r4r2_val_sym(sym) result(rad)

   !> Element symbol
   character(len=*), intent(in) :: sym

   !> r4/r2 expectation value
   real(wp) :: rad

   rad = get_r4r2_val(to_number(sym))

end function get_r4r2_val_sym


!> Get r4/r2 expectation value for a given atomic number
elemental function get_r4r2_val_num(num) result(rad)

   !> Atomic number
   integer, intent(in) :: num

   !> r4/r2 expectation value
   real(wp) :: rad

   if (num > 0 .and. num <= size(sqrt_z_r4_over_r2)) then
      rad = sqrt_z_r4_over_r2(num)
   else
      rad = 0.0_wp
   end if

end function get_r4r2_val_num


end module dftd4_data_r4r2
