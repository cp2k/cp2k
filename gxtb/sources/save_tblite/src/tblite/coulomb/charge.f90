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

!> @dir tblite/coulomb/charge
!> Contains the implementation of the isotropic second-order electrostatics

!> @file tblite/coulomb/charge.f90
!> Provides a reexport of the isotropic second-order electrostatic implementations

!> Proxy module for defining isotropic electrostatic interactions
module tblite_coulomb_charge
   use tblite_coulomb_charge_effective, only : effective_coulomb, new_effective_coulomb
   use tblite_coulomb_charge_gamma, only : gamma_coulomb, new_gamma_coulomb
   use tblite_coulomb_charge_type, only : coulomb_charge_type
   use tblite_utils_average, only : average_type, new_average, average_id
   implicit none
   private

   public :: coulomb_charge_type, coulomb_kernel
   public :: effective_coulomb, new_effective_coulomb
   public :: average_type, new_average, average_id
   public :: gamma_coulomb, new_gamma_coulomb


   !> Possible interaction kernels
   type :: enum_coulomb_kernel
      !> Effective Klopman-Ohno interaction kernel
      integer :: effective = 1
      !> DFTB γ-functional interaction kernel
      integer :: dftb_gamma = 2
   end type enum_coulomb_kernel

   !> Actual enumerator for possible interaction kernels
   type(enum_coulomb_kernel), parameter :: coulomb_kernel = enum_coulomb_kernel()

end module tblite_coulomb_charge
