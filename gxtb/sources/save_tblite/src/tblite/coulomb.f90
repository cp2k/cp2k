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

!> @dir tblite/coulomb
!> Contains the Coulomb related interactions.

!> @file tblite/coulomb.f90
!> Reexports of the Coulombic interaction containers.

!> Proxy module for handling Coulomb-type interactions
module tblite_coulomb
   use tblite_coulomb_cache, only : coulomb_cache
   use tblite_coulomb_charge, only : coulomb_charge_type, coulomb_kernel, &
      & effective_coulomb, new_effective_coulomb, gamma_coulomb, new_gamma_coulomb
   use tblite_coulomb_firstorder, only : onsite_firstorder, new_onsite_firstorder
   use tblite_coulomb_fourthorder, only : onsite_fourthorder, new_onsite_fourthorder
   use tblite_coulomb_multipole, only : damped_multipole, multipole_damping, &
      & gfn2_multipole, new_gfn2_multipole, gxtb_multipole, new_gxtb_multipole
   use tblite_coulomb_thirdorder, only : thirdorder_type, thirdorder_kernel, & 
      & onsite_thirdorder, new_onsite_thirdorder, twobody_thirdorder, &
      & new_twobody_thirdorder
   use tblite_coulomb_type, only : coulomb_type
   use tblite_utils, only : average_type, new_average, average_id
   implicit none

   public :: coulomb_type
   public :: coulomb_cache
   public :: onsite_firstorder, new_onsite_firstorder
   public :: coulomb_charge_type, coulomb_kernel
   public :: effective_coulomb, new_effective_coulomb
   public :: average_type, new_average, average_id
   public :: gamma_coulomb, new_gamma_coulomb
   public :: damped_multipole, multipole_damping
   public :: gfn2_multipole, new_gfn2_multipole
   public :: gxtb_multipole, new_gxtb_multipole
   public :: thirdorder_type, thirdorder_kernel
   public :: onsite_thirdorder, new_onsite_thirdorder
   public :: twobody_thirdorder, new_twobody_thirdorder
   public :: onsite_fourthorder, new_onsite_fourthorder
end module tblite_coulomb
