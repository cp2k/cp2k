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

!> @dir tblite/solvation
!> Contains the implementation of the implicit solvation models.

!> @file tblite/solvation.f90
!> Provides reexports of the implict solvation model related implementations.

!> Proxy module for implicit solvation models.
module tblite_solvation
   use mctc_env, only : error_type, fatal_error
   use mctc_io, only : structure_type
   use tblite_solvation_alpb, only : alpb_solvation, new_alpb, alpb_input, born_kernel
   use tblite_solvation_cpcm, only : cpcm_solvation, new_cpcm, cpcm_input
   use tblite_solvation_cds, only : cds_solvation, new_cds, cds_input
   use tblite_solvation_shift, only : shift_solvation, new_shift, shift_input, solution_state
   use tblite_solvation_data, only : solvent_data, get_solvent_data
   use tblite_solvation_input, only : solvation_input
   use tblite_solvation_type, only : solvation_type
   use tblite_solvation_data_alpb, only: get_alpb_param
   use tblite_solvation_data_cds, only: get_cds_param
   use tblite_solvation_data_shift, only: get_shift_param
   implicit none
   private

   public :: alpb_solvation, new_alpb, alpb_input, born_kernel
   public :: cpcm_solvation, new_cpcm, cpcm_input
   public :: cds_solvation, new_cds, cds_input
   public :: shift_solvation, new_shift, shift_input
   public :: solvent_data, get_solvent_data
   public :: solvation_input, new_solvation, solvation_type
   public :: new_solvation_cds
   public :: new_solvation_shift, solution_state

contains

!> Create new solvation model from input data
subroutine new_solvation(solv, mol, input, error, method)
   !> Instance of the solvation model
   class(solvation_type), allocatable, intent(out) :: solv
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Input data
   type(solvation_input), intent(in) :: input
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Method for parameter selection
   character(len=*), optional, intent(in) :: method
   !> scratch input
   type(alpb_input), allocatable :: scratch_input

   if (allocated(input%alpb)) then
      scratch_input = input%alpb
      ! ALPB/GBSA with empirical parameters
      if (allocated(input%alpb%solvent) .and. present(method)) then
         call get_alpb_param(scratch_input, mol, method, error)
         if(allocated(error)) then
            call fatal_error(error, "No ALPB/GBSA parameters found for the method/solvent")
            return
         end if
      end if
      solv = alpb_solvation(mol, scratch_input, method)
      return
   end if

   if (allocated(input%cpcm)) then
      solv = cpcm_solvation(mol, input%cpcm)
      return
   end if

   call fatal_error(error, "Unknown solvation model")
end subroutine new_solvation

!> Create new cds solvation model from input data
subroutine new_solvation_cds(solv, mol, input, error, method)
   !> Instance of the solvation model
   class(solvation_type), allocatable, intent(out) :: solv
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Input data
   type(solvation_input), intent(in) :: input
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Method for parameter selection
   character(len=*), optional, intent(in) :: method
   !> scratch input
   type(cds_input), allocatable :: scratch_input

   if (allocated(input%cds) .and. &
      & allocated(input%cds%solvent) .and. present(method)) then
      scratch_input = input%cds
      call get_cds_param(scratch_input, mol, method, error)
      if(allocated(error)) then
         call fatal_error(error, "No CDS parameters found for the method/solvent")
         return
      end if
      solv = cds_solvation(mol, scratch_input, method)
      return
   end if

   call fatal_error(error, "Unknown CDS model")
end subroutine new_solvation_cds

!> Create new solvation shift from input data
subroutine new_solvation_shift(solv, input, error, method)
   !> Instance of the solvation model
   class(solvation_type), allocatable, intent(out) :: solv
   !> Input data
   type(solvation_input), intent(in) :: input
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Method for parameter selection
   character(len=*), optional, intent(in) :: method
   !> scratch input
   type(shift_input), allocatable :: scratch_input

   if (allocated(input%shift) .and. &
      & allocated(input%shift%solvent) .and. present(method)) then
      scratch_input = input%shift
      call get_shift_param(scratch_input, method, error)
      if(allocated(error)) then
         call fatal_error(error, "No shift parameters found for the method/solvent")
         return
      end if
      solv = shift_solvation(scratch_input)
      return
   end if

   call fatal_error(error, "Unknown solvation shift")
end subroutine new_solvation_shift

end module tblite_solvation
