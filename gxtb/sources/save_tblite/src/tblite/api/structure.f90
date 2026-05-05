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

!> @file tblite/api/structure.f90
!> Provides API exports for the #tblite_structure handle.

!> API export for working with molecular structure data objects
module tblite_api_structure
   use, intrinsic :: iso_c_binding
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type, new
   use tblite_api_error, only : vp_error
   use tblite_api_version, only : namespace
   use tblite_cutoff, only : wrap_to_central_cell
   implicit none
   private

   public :: vp_structure
   public :: new_structure_api, delete_structure_api, update_structure_geometry_api
   public :: update_structure_charge_api, update_structure_uhf_api


   !> Void pointer to molecular structure data
   type :: vp_structure
      !> Actual payload
      type(structure_type) :: ptr
   end type vp_structure


   logical, parameter :: debug = .false.


contains


!> Create new molecular structure data (quantities in Bohr)
function new_structure_api(verror, natoms, numbers, positions, c_charge, c_uhf, &
      & c_lattice, c_periodic) result(vmol) &
      & bind(C, name=namespace//"new_structure")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   integer(c_int), value, intent(in) :: natoms
   integer(c_int), intent(in) :: numbers(natoms)
   real(c_double), intent(in) :: positions(3, natoms)
   real(c_double), intent(in), optional :: c_charge
   real(wp), allocatable :: charge
   integer(c_int), intent(in), optional :: c_uhf
   integer, allocatable :: uhf
   real(c_double), intent(in), optional :: c_lattice(3, 3)
   real(wp), allocatable :: lattice(:, :)
   logical(c_bool), intent(in), optional :: c_periodic(3)
   logical, allocatable :: periodic(:)
   type(vp_structure), pointer :: mol
   type(c_ptr) :: vmol

   if (debug) print '("[Info]",1x, a)', "new_structure"

   vmol = c_null_ptr

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (present(c_lattice)) then
      allocate(lattice(3, 3))
      lattice(:, :) = c_lattice
   end if
   if (present(c_periodic)) then
      allocate(periodic(3))
      periodic(:) = c_periodic
   end if
   if (present(c_charge)) then
      charge = c_charge
   end if
   if (present(c_uhf)) then
      uhf = c_uhf
   end if

   allocate(mol)
   call new(mol%ptr, numbers, positions, lattice=lattice, periodic=periodic, &
      charge=charge, uhf=uhf)
   vmol = c_loc(mol)

   call wrap_to_central_cell(mol%ptr%xyz, mol%ptr%lattice, mol%ptr%periodic)

   call verify_structure(error%ptr, mol%ptr)

end function new_structure_api


!> Delete molecular structure data
subroutine delete_structure_api(vmol) &
      & bind(C, name=namespace//"delete_structure")
   type(c_ptr), intent(inout) :: vmol
   type(vp_structure), pointer :: mol

   if (debug) print '("[Info]",1x, a)', "delete_structure"

   if (c_associated(vmol)) then
      call c_f_pointer(vmol, mol)

      deallocate(mol)
      vmol = c_null_ptr
   end if

end subroutine delete_structure_api


!> Update coordinates and lattice parameters (quantities in Bohr)
subroutine update_structure_geometry_api(verror, vmol, positions, lattice) &
      & bind(C, name=namespace//"update_structure_geometry")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vmol
   type(vp_structure), pointer :: mol
   real(c_double), intent(in) :: positions(3, *)
   real(c_double), intent(in), optional :: lattice(3, 3)

   if (debug) print '("[Info]",1x, a)', "update_structure"

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vmol)) then
      call fatal_error(error%ptr, "Molecular structure data is missing")
      return
   end if
   call c_f_pointer(vmol, mol)

   if (mol%ptr%nat <= 0 .or. mol%ptr%nid <= 0 .or. .not.allocated(mol%ptr%num) &
      & .or. .not.allocated(mol%ptr%id) .or. .not.allocated(mol%ptr%xyz)) then
      call fatal_error(error%ptr, "Invalid molecular structure data provided")
      return
   end if

   mol%ptr%xyz(:, :) = positions(:3, :mol%ptr%nat)
   if (present(lattice)) then
      mol%ptr%lattice(:, :) = lattice(:3, :3)
   end if

   call wrap_to_central_cell(mol%ptr%xyz, mol%ptr%lattice, mol%ptr%periodic)

   call verify_structure(error%ptr, mol%ptr)

end subroutine update_structure_geometry_api


!> Update coordinates and lattice parameters (quantities in Bohr)
subroutine update_structure_charge_api(verror, vmol, c_charge) &
      & bind(C, name=namespace//"update_structure_charge")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vmol
   type(vp_structure), pointer :: mol
   real(c_double), intent(in), optional :: c_charge

   if (debug) print '("[Info]",1x, a)', "update_structure_charge"

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vmol)) then
      call fatal_error(error%ptr, "Molecular structure data is missing")
      return
   end if
   call c_f_pointer(vmol, mol)

   if (mol%ptr%nat <= 0 .or. mol%ptr%nid <= 0 .or. .not.allocated(mol%ptr%num) &
      & .or. .not.allocated(mol%ptr%id) .or. .not.allocated(mol%ptr%xyz)) then
      call fatal_error(error%ptr, "Invalid molecular structure data provided")
      return
   end if

   if (present(c_charge)) then
      mol%ptr%charge = c_charge
   end if

end subroutine update_structure_charge_api


!> Update coordinates and lattice parameters (quantities in Bohr)
subroutine update_structure_uhf_api(verror, vmol, c_uhf) &
      & bind(C, name=namespace//"update_structure_uhf")
   type(c_ptr), value :: verror
   type(vp_error), pointer :: error
   type(c_ptr), value :: vmol
   type(vp_structure), pointer :: mol
   integer(c_int), intent(in), optional :: c_uhf

   if (debug) print '("[Info]",1x, a)', "update_structure_uhf"

   if (.not.c_associated(verror)) return
   call c_f_pointer(verror, error)

   if (.not.c_associated(vmol)) then
      call fatal_error(error%ptr, "Molecular structure data is missing")
      return
   end if
   call c_f_pointer(vmol, mol)

   if (mol%ptr%nat <= 0 .or. mol%ptr%nid <= 0 .or. .not.allocated(mol%ptr%num) &
      & .or. .not.allocated(mol%ptr%id) .or. .not.allocated(mol%ptr%xyz)) then
      call fatal_error(error%ptr, "Invalid molecular structure data provided")
      return
   end if

   if (present(c_uhf)) then
      mol%ptr%uhf = c_uhf
   end if

end subroutine update_structure_uhf_api


!> Cold fusion check
subroutine verify_structure(error, mol)
   type(error_type), allocatable, intent(out) :: error
   type(structure_type), intent(in) :: mol
   integer :: iat, jat, stat
   stat = 0
   do iat = 1, mol%nat
      do jat = 1, iat - 1
         if (norm2(mol%xyz(:, jat) - mol%xyz(:, iat)) < 1.0e-9_wp) stat = stat + 1
      end do
   end do
   if (stat > 0) then
      call fatal_error(error, "Too close interatomic distances found")
   end if
end subroutine verify_structure


end module tblite_api_structure
