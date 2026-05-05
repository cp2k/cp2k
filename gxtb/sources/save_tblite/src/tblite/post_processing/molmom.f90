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

!> @file tblite/post_processing/molmom.f90
!> Implements the calculation of molecular moments as post processing methods.
module tblite_post_processing_molecular_moments
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_basis_type, only : basis_type
   use tblite_container, only : container_cache
   use tblite_context, only : context_type
   use tblite_double_dictionary, only : double_dictionary_type
   use tblite_integral_type, only : integral_type
   use tblite_output_format, only : format_string
   use tblite_param_molecular_moments, only : molecular_multipole_record
   use tblite_post_processing_type, only : post_processing_type
   use tblite_results, only : results_type
   use tblite_timer, only : timer_type, format_time
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_wavefunction_mulliken, only : get_molecular_dipole_moment, get_molecular_quadrupole_moment
   implicit none
   private

   public :: molecular_moments, new_molecular_moments

   type, extends(post_processing_type) :: molecular_moments
      logical :: comp_dipm = .true. , comp_qm = .true.
   contains
      procedure :: compute
      procedure :: print_timer
   end type molecular_moments

   character(len=27), parameter :: label = "Molecular Multipole Moments"
   ! This is not thread-safe
   !type(timer_type) :: timer

contains

subroutine new_molecular_moments(new_molmom_type, param)
   type(molecular_moments), intent(inout) :: new_molmom_type
   type(molecular_multipole_record), intent(in) :: param
   
   new_molmom_type%label = label

   new_molmom_type%comp_dipm = param%moldipm
   new_molmom_type%comp_qm = param%molqp

end subroutine new_molecular_moments

subroutine compute(self, mol, wfn, integrals, calc, cache_list, ctx, prlevel, dict)
   class(molecular_moments),intent(inout) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction strcuture data
   type(wavefunction_type), intent(in) :: wfn
   !> integral container
   type(integral_type), intent(in) :: integrals
   !> calculator instance
   type(xtb_calculator), intent(in) :: calc
   !> Cache list for storing caches of various interactions
   type(container_cache), intent(inout) :: cache_list(:)
   !> Context container for writing to stdout
   type(context_type), intent(inout) :: ctx
   !> Print level
   integer, intent(in) :: prlevel
   !> Dictionary for storing results
   type(double_dictionary_type), intent(inout) :: dict
   real(wp) :: dipm(3), qp(6)

   !call timer%push("total")
   if (self%comp_dipm) then
      !call timer%push("dipole")
      call get_molecular_dipole_moment(mol, wfn%qat(:, 1), wfn%dpat(:, :, 1), dipm)
      call dict%add_entry("molecular-dipole", dipm)
      !call timer%pop()
   end if
   if (self%comp_qm) then
      !call timer%push("quadrupole")
      call get_molecular_quadrupole_moment(mol, wfn%qat(:, 1), wfn%dpat(:, :, 1), wfn%qpat(:, :, 1), qp)
      call dict%add_entry("molecular-quadrupole", qp)
      !call timer%pop()
   end if
   !call timer%pop()
end subroutine compute

subroutine print_timer(self, prlevel, ctx)
   !> Instance of the interaction container
   class(molecular_moments), intent(in) :: self
   integer :: prlevel
   type(context_type) :: ctx
   real(wp) :: ttime, stime
   integer :: it
   character(len=*), parameter :: labels(*) = [character(len=20):: &
      "dipole", "quadrupole"  ]



   ! if (prlevel > 2) then
   !    call ctx%message(label//" timing details:")
   !    ttime = timer%get("total")
   !    call ctx%message(" total:"//repeat(" ", 16)//format_time(ttime))
   !    do it = 1, size(labels)
   !       stime = timer%get(labels(it))
   !       if (stime <= epsilon(0.0_wp)) cycle
   !       call ctx%message(" - "//labels(it)//format_time(stime) &
   !          & //" ("//format_string(int(stime/ttime*100), '(i3)')//"%)")
   !    end do
   !    call ctx%message("")
   ! end if
end subroutine print_timer

end module tblite_post_processing_molecular_moments
