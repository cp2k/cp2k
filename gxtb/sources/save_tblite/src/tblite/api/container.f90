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

!> @file tblite/api/container.f90
!> Provides API exports for the #tblite_container handle.

!> API export for managing interaction containers
module tblite_api_container
   use, intrinsic :: iso_c_binding
   use mctc_env, only : error_type, fatal_error, wp
   use mctc_io, only : structure_type
   use tblite_api_version, only : namespace
   use tblite_api_calculator, only : vp_calculator
   use tblite_api_context, only : vp_context
   use tblite_api_structure, only : vp_structure
   use tblite_basis, only : basis_type
   use tblite_container, only : container_type
   use tblite_data_spin, only : get_spin_constant
   use tblite_external_field, only : electric_field
   use tblite_spin, only : spin_polarization, new_spin_polarization
   use tblite_api_utils, only: c_f_character
   implicit none
   private

   public :: vp_container, delete_container_api


   !> Void pointer to a container instance
   type :: vp_container
      !> Actual container
      class(container_type), allocatable :: ptr
   end type vp_container

   logical, parameter :: debug = .false.

contains


subroutine push_back_api(vctx, vcalc, vcont) &
      & bind(C, name=namespace//"calculator_push_back")
   type(c_ptr), value :: vctx
   type(vp_context), pointer :: ctx
   type(c_ptr), value :: vcalc
   type(vp_calculator), pointer :: calc
   type(c_ptr), intent(inout) :: vcont
   type(vp_container), pointer :: cont

   type(error_type), allocatable :: error

   if (debug) print '("[Info]", 1x, a)', "calculator_push_back"

   if (.not.c_associated(vctx)) return
   call c_f_pointer(vctx, ctx)

   if (.not.c_associated(vcalc)) then
      call fatal_error(error, "Calculator object is missing")
      call ctx%ptr%set_error(error)
      return
   end if
   call c_f_pointer(vcalc, calc)

   if (.not.c_associated(vcont)) then
      call fatal_error(error, "Container object is missing")
      call ctx%ptr%set_error(error)
      return
   end if
   call c_f_pointer(vcont, cont)

   ! Hacky way to propagate spin channel information to calculator
   select type(tcont => cont%ptr)
   type is(spin_polarization)
      calc%nspin = 2
   end select

   call calc%ptr%push_back(cont%ptr)

   deallocate(cont)
   vcont = c_null_ptr
end subroutine push_back_api


function new_electric_field_api(efield) result(vcont) &
      & bind(C, name=namespace//"new_electric_field")
   real(c_double), intent(in) :: efield(3)
   type(c_ptr) :: vcont
   type(vp_container), pointer :: cont

   allocate(cont)
   cont%ptr = electric_field(efield)
   vcont = c_loc(cont)
end function new_electric_field_api


function new_spin_polarization_api(vctx, vmol, vcalc, wscale) result(vcont) &
      & bind(C, name=namespace//"new_spin_polarization")
   type(c_ptr), value :: vctx
   type(vp_context), pointer :: ctx
   type(c_ptr), value :: vmol
   type(vp_structure), pointer :: mol
   type(c_ptr), value :: vcalc
   type(vp_calculator), pointer :: calc
   type(c_ptr) :: vcont
   type(vp_container), pointer :: cont
   real(c_double), value :: wscale

   type(error_type), allocatable :: error
   type(spin_polarization), allocatable :: spin
   real(wp), allocatable :: wll(:, :, :)

   if (debug) print '("[Info]", 1x, a)', "new_spin_polarization"
   vcont = c_null_ptr

   if (.not.c_associated(vctx)) return
   call c_f_pointer(vctx, ctx)

   if (.not.c_associated(vmol)) then
      call fatal_error(error, "Molecular structure data is missing")
      call ctx%ptr%set_error(error)
      return
   end if
   call c_f_pointer(vmol, mol)

   if (.not.c_associated(vcalc)) then
      call fatal_error(error, "Calculator object is missing")
      call ctx%ptr%set_error(error)
      return
   end if
   call c_f_pointer(vcalc, calc)

   allocate(spin)
   call get_spin_constants(wll, mol%ptr, calc%ptr%bas)
   wll(:, :, :) = wscale * wll
   call new_spin_polarization(spin, mol%ptr, wll, calc%ptr%bas%nsh_id)

   allocate(cont)
   call move_alloc(spin, cont%ptr)
   vcont = c_loc(cont)
end function new_spin_polarization_api


subroutine get_spin_constants(wll, mol, bas)
   real(wp), allocatable, intent(out) :: wll(:, :, :)
   type(structure_type), intent(in) :: mol
   type(basis_type), intent(in) :: bas

   integer :: izp, ish, jsh, il, jl

   allocate(wll(bas%nsh, bas%nsh, mol%nid), source=0.0_wp)

   do izp = 1, mol%nid
      do ish = 1, bas%nsh_id(izp)
         il = bas%cgto(ish, izp)%raw%ang
         do jsh = 1, bas%nsh_id(izp)
            jl = bas%cgto(jsh, izp)%raw%ang
            wll(jsh, ish, izp) = get_spin_constant(jl, il, mol%num(izp))
         end do
      end do
   end do
end subroutine get_spin_constants


subroutine delete_container_api(vcont) &
      & bind(C, name=namespace//"delete_container")
   type(c_ptr), intent(inout) :: vcont
   type(vp_container), pointer :: cont

   if (debug) print '("[Info]", 1x, a)', "delete_container"

   if (c_associated(vcont)) then
      call c_f_pointer(vcont, cont)

      deallocate(cont)
      vcont = c_null_ptr
   end if
end subroutine delete_container_api


end module tblite_api_container
