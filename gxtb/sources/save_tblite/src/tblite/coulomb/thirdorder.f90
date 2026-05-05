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

!> @dir tblite/coulomb/thirdorder
!> Contains the implementation isotropic third-order tight-binding

!> @file tblite/coulomb/thirdoder.f90
!> Provides a reexport of the isotropic third-order tigh-binding implementations

!> Proxy module for defining isotropic third-order tight-binding
module tblite_coulomb_thirdorder
   use tblite_coulomb_thirdorder_type, only : thirdorder_type
   use tblite_coulomb_thirdorder_onsite, only : onsite_thirdorder, new_onsite_thirdorder
   use tblite_coulomb_thirdorder_twobody, only : twobody_thirdorder, new_twobody_thirdorder
   implicit none
   private

   public :: thirdorder_type, thirdorder_kernel
   public :: onsite_thirdorder, new_onsite_thirdorder
   public :: twobody_thirdorder, new_twobody_thirdorder

   !> Possible third-order interaction kernels
   type :: enum_thirdorder_kernel
      !> Purely onsite third-order tight-binding interaction kernel
      integer :: onsite = 1
      !> Two-body third-order tight-binding interaction kernel
      integer :: twobody = 2
   end type enum_thirdorder_kernel

   !> Actual enumerator for possible third-order interaction kernels
   type(enum_thirdorder_kernel), parameter :: thirdorder_kernel = enum_thirdorder_kernel()

end module tblite_coulomb_thirdorder
