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

!> Element specific data needed for the DFT-D4 model
module dftd4_data
   use dftd4_data_covrad, only : get_covalent_rad
   use dftd4_data_en, only : get_electronegativity
   use dftd4_data_hardness, only : get_hardness
   use dftd4_data_r4r2, only : get_r4r2_val
   use dftd4_data_wfpair, only : get_wfpair_val
   use dftd4_data_zeff, only : get_effective_charge
   implicit none
   public

end module dftd4_data
