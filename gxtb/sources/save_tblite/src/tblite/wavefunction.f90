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

!> @dir tblite/wavefunction
!> Contains wavefunction related information

!> @file tblite/wavefunction.f90
!> Provides reexports of wavefunction related procedures and types

!> Proxy module for wavefunction related types and procedures
module tblite_wavefunction
   use tblite_wavefunction_guess, only : sad_guess, eeq_guess, eeqbc_guess, &
      & shell_partition
   use tblite_wavefunction_mulliken, only : get_molecular_dipole_moment, &
      & get_molecular_quadrupole_moment, get_mayer_bond_orders
   use tblite_wavefunction_restart, only : load_wavefunction, save_wavefunction
   use tblite_wavefunction_spin, only : magnet_to_updown, updown_to_magnet, &
      & pot_magnet_to_updown, pot_updown_to_magnet
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction, &
      & get_density_matrix, get_alpha_beta_occupation
   implicit none
   public

end module tblite_wavefunction
