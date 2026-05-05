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
module tblite_api_solvation
   use, intrinsic :: iso_c_binding
   use mctc_env, only : error_type, fatal_error, wp
   use mctc_io, only : structure_type
   use tblite_api_version, only : namespace
   use tblite_api_container, only : vp_container
   use tblite_api_error, only : vp_error
   use tblite_api_structure, only : vp_structure
   use tblite_container, only : container_type, container_list
   use tblite_solvation, only : solvation_input, cpcm_input, alpb_input, &
      & solvent_data, get_solvent_data, solvation_type, new_solvation, solution_state, &
      & new_solvation_cds, new_solvation_shift, cds_input, shift_input, born_kernel
   use tblite_api_utils, only: c_f_character
   implicit none
   private

   public :: new_gb_solvation_epsilon_api, new_alpb_solvation_solvent_api, &
      & new_cpcm_solvation_epsilon_api

   enum, bind(c)
      enumerator :: &
         solvation_gbe = 10, &
         solvation_alpb_gfn1 = 11, &
         solvation_alpb_gfn2 = 12, &
         solvation_gb = 20, &
         solvation_gbsa_gfn1 = 21, &
         solvation_gbsa_gfn2 = 22
   end enum


   logical, parameter :: debug = .false.

contains


function new_cpcm_solvation_epsilon_api(verr, vmol, eps) result(vcont) &
   & bind(C, name=namespace//"new_cpcm_solvation_epsilon")
   type(c_ptr), value :: verr
   type(vp_error), pointer :: err
   type(c_ptr), value :: vmol
   type(vp_structure), pointer :: mol
   real(kind=c_double), value :: eps
   type(c_ptr) :: vcont
   type(vp_container), pointer :: cont

   type(solvation_input) :: solvmodel
   class(solvation_type), allocatable :: solv

   if (debug) print '("[Info]", 1x, a)', "new_cpcm_solvation_epsilon"
   vcont = c_null_ptr

   if (.not.c_associated(verr)) return
   call c_f_pointer(verr, err)

   if (.not.c_associated(vmol)) then
      call fatal_error(err%ptr, "Molecular structure data is missing")
      return
   end if
   call c_f_pointer(vmol, mol)

   solvmodel%cpcm = cpcm_input(eps)
   call new_solvation(solv, mol%ptr, solvmodel, err%ptr)
   if (allocated(err%ptr)) return
   
   allocate(cont)
   call move_alloc(solv, cont%ptr)
   
   vcont = c_loc(cont)
end function new_cpcm_solvation_epsilon_api

function new_gb_solvation_epsilon_api(verr, vmol, eps, version, born_type) result(vcont) &
   & bind(C, name=namespace//"new_gb_solvation_epsilon")
   type(c_ptr), value :: verr
   type(vp_error), pointer :: err
   type(c_ptr), value :: vmol
   type(vp_structure), pointer :: mol
   real(c_double), value :: eps
   integer(c_int), value :: version
   integer(c_int), value :: born_type
   type(c_ptr) :: vcont
   type(vp_container), pointer :: cont

   type(solvation_input) :: solvmodel
   class(solvation_type), allocatable :: solv
   logical :: alpb

   if (debug) print '("[Info]", 1x, a)', "new_alpb_solvation_epsilon"
   vcont = c_null_ptr

   if (.not.c_associated(verr)) return
   call c_f_pointer(verr, err)

   if (.not.c_associated(vmol)) then
      call fatal_error(err%ptr, "Molecular structure data is missing")
      return
   end if
   call c_f_pointer(vmol, mol)

   select case(version)
   case default
      call fatal_error(err%ptr, "Unknown solvation model requested")
      return
   case(solvation_gbe)
      alpb = .true.
   case(solvation_gb)
      alpb = .false.
   end select

   solvmodel%alpb = alpb_input(eps, alpb=alpb, kernel=born_type)
   call new_solvation(solv, mol%ptr, solvmodel, err%ptr)
   if (allocated(err%ptr)) return
   
   allocate(cont)
   call move_alloc(solv, cont%ptr)
   
   vcont = c_loc(cont)
end function new_gb_solvation_epsilon_api

function new_alpb_solvation_solvent_api(verr, vmol, csolvstr, version, crefstate) result(vcont) &
   & bind(C, name=namespace//"new_alpb_solvation_solvent")
   type(c_ptr), value :: verr
      type(vp_error), pointer :: err
   type(c_ptr), value :: vmol
   type(vp_structure), pointer :: mol
   type(c_ptr) :: vcont
   character(kind=c_char), intent(in) :: csolvstr(*)
   integer(c_int), value :: version
   integer(c_int), value :: crefstate
   type(vp_container), pointer :: cont

   type(solvation_input) :: solvmodel
   type(solvent_data) :: solvent
   character(len=:), allocatable :: solvstr, method
   integer :: kernel, sol_state
   logical :: alpb

   class(container_list), allocatable :: list
   class(container_type), allocatable :: tmp_cont
   class(solvation_type), allocatable :: solv, cds, shift

   if (debug) print '("[Info]", 1x, a)', "new_alpb_solvation_solvent"
   vcont = c_null_ptr

   if (.not.c_associated(verr)) return
   call c_f_pointer(verr, err)

   if (.not.c_associated(vmol)) then
      call fatal_error(err%ptr, "Molecular structure data is missing")
      return
   end if
   call c_f_pointer(vmol, mol)

   call c_f_character(csolvstr, solvstr)
   solvent = get_solvent_data(solvstr)
   if (solvent%eps <= 0.0_wp) then
      call fatal_error(err%ptr, "String value for epsilon was not found among database of solvents")
      return
   end if

   select case(version)
   case(solvation_alpb_gfn1, solvation_alpb_gfn2)
      method = merge("gfn2", "gfn1", version == solvation_alpb_gfn2)
      kernel = born_kernel%p16
      alpb = .true.
   case(solvation_gbsa_gfn1, solvation_gbsa_gfn2)
      method = merge("gfn2", "gfn1", version == solvation_gbsa_gfn2)
      kernel = born_kernel%still
      alpb = .false.
   case default
      call fatal_error(err%ptr, "Unknown value for model version")
      return
   end select

   sol_state = int(crefstate)

   solvmodel%alpb = alpb_input(solvent%eps, solvent=solvent%solvent, kernel=kernel, alpb=alpb)
   solvmodel%cds = cds_input(alpb=alpb, solvent=solvent%solvent)
   solvmodel%shift = shift_input(alpb=alpb, solvent=solvent%solvent, state=sol_state)

   allocate(list)
   call new_solvation(solv, mol%ptr, solvmodel, err%ptr, method)
   if (allocated(err%ptr)) return
   call move_alloc(solv, tmp_cont)
   call list%push_back(tmp_cont)
   
   call new_solvation_cds(cds, mol%ptr, solvmodel, err%ptr, method)
   if (allocated(err%ptr)) return
   call move_alloc(cds, tmp_cont)
   call list%push_back(tmp_cont)

   call new_solvation_shift(shift, solvmodel, err%ptr, method)
   if (allocated(err%ptr)) return
   call move_alloc(shift, tmp_cont)
   call list%push_back(tmp_cont)

   allocate(cont)
   call move_alloc(list, cont%ptr)

   vcont = c_loc(cont)
end function new_alpb_solvation_solvent_api

end module tblite_api_solvation