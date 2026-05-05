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

!> @dir tblite/basis
!> Contains basis set related implementations

!> @file tblite/basis.f90
!> Reexports the basis set related data classes

!> Proxy module for the basis set related procedures and types
module tblite_basis
   use tblite_basis_ortho, only : orthogonalize
   use tblite_basis_slater, only : slater_to_gauss
   use tblite_basis_type, only : basis_type, new_basis, get_cutoff, basis_set, &
      & cgto_type, new_cgto, cgto_container
   use tblite_basis_qvszp, only : new_qvszp_basis, new_qvszp_cgto, &
      & qvszp_cgto_type, qvszp_basis_type
   implicit none
   private

   public :: basis_type, new_basis, get_cutoff, basis_set
   public :: cgto_type, new_cgto, cgto_container
   public :: orthogonalize
   public :: slater_to_gauss
   public :: new_qvszp_basis, new_qvszp_cgto, qvszp_cgto_type, qvszp_basis_type
end module tblite_basis
