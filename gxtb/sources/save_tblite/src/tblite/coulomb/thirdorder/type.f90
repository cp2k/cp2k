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

!> @file tblite/coulomb/thirdoder/type.f90
!> Provides a general base class for third-order tight-binding

!> Proxy module for defining isotropic third-order tight-binding
module tblite_coulomb_thirdorder_type
   use mctc_env, only : wp
   use tblite_coulomb_type, only : coulomb_type
   implicit none
   private

   !> Abstract base class for third-order tight-binding
   type, public, extends(coulomb_type), abstract :: thirdorder_type
      !> Whether the third-order contribution is shell-dependent
      logical :: shell_resolved
      !> Number of shell for each atom
      integer, allocatable :: nsh_at(:)
      !> Shell offset for each atom
      integer, allocatable :: ish_at(:)
      !> Hubbard derivatives for each species
      real(wp), allocatable :: hubbard_derivs(:, :)
      !> Logical indicating whether there is a CN dependence 
      logical :: cn_dep
   contains 
      !> Get information about density dependent quantities used in the energy
      procedure :: variable_info
   end type thirdorder_type

contains

!> Get information about density dependent quantities used in the energy
pure function variable_info(self) result(info)
   use tblite_scf_info, only : scf_info, atom_resolved, shell_resolved
   !> Instance of the third-order tight-binding container
   class(thirdorder_type), intent(in) :: self
   !> Information on the required potential data
   type(scf_info) :: info

   info = scf_info(charge=merge(shell_resolved, atom_resolved, self%shell_resolved))
end function variable_info

end module tblite_coulomb_thirdorder_type
