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

!> @file tblite/post_processing/type.f90
!> Implements post processing container abstract class, and the collection of computing caches.
module tblite_post_processing_type
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_basis_type, only : basis_type
   use tblite_container_cache, only : container_cache
   use tblite_context, only : context_type
   use tblite_double_dictionary, only : double_dictionary_type
   use tblite_integral_type, only : integral_type
   use tblite_results, only : results_type
   use tblite_timer, only : timer_type, format_time
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_xtb_calculator, only : xtb_calculator
   implicit none
   private
   public :: post_processing_type, collect_containers_caches

   type, abstract :: post_processing_type
      character(len=:), allocatable :: label

   contains
      !> Setup container
      procedure(compute), deferred :: compute
      procedure :: info
      procedure :: print_timer
   end type  post_processing_type

   type(timer_type) :: timer
   abstract interface
      subroutine compute(self, mol, wfn, integrals, calc, cache_list, ctx, prlevel, dict)
      import :: post_processing_type, structure_type, wavefunction_type, integral_type, xtb_calculator, &
         & context_type, container_cache, double_dictionary_type
      class(post_processing_type),intent(inout) :: self
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
      end subroutine compute
   end interface
contains


pure function info(self, verbosity, indent) result(str)
   !> Instance of the interaction container
   class(post_processing_type), intent(in) :: self
   !> Verbosity level
   integer, intent(in) :: verbosity
   !> Indentation level
   character(len=*), intent(in) :: indent
   !> Information on the container
   character(len=:), allocatable :: str

   if (allocated(self%label)) then
      str = indent // self%label
   else
      str = "Unknown"
   end if
end function info

subroutine print_timer(self, prlevel, ctx)
   !> Instance of the interaction container
   class(post_processing_type), intent(in) :: self
   !> Print level
   integer :: prlevel
   !> Context container for writing to stdout
   type(context_type) :: ctx
   real(wp) :: ttime

   if (prlevel > 1) then
      ttime = timer%get("total")
      call ctx%message(" total:"//repeat(" ", 16)//format_time(ttime))
      call ctx%message("")
   end if
end subroutine print_timer

subroutine collect_containers_caches(rcache, ccache, hcache, dcache, icache, calc, cache_list)
   type(container_cache), allocatable, intent(inout) :: rcache, ccache, hcache, dcache, icache
   type(container_cache), allocatable, intent(inout) :: cache_list(:)
   type(xtb_calculator), intent(in) :: calc

   allocate(cache_list(5))
   if (allocated(calc%repulsion)) call move_alloc(rcache%raw, cache_list(1)%raw)
   if (allocated(calc%coulomb)) call move_alloc(ccache%raw, cache_list(2)%raw)
   if (allocated(calc%halogen)) call move_alloc(hcache%raw, cache_list(3)%raw)
   if (allocated(calc%dispersion)) call move_alloc(dcache%raw, cache_list(4)%raw)
   if (allocated(calc%interactions)) call move_alloc(icache%raw, cache_list(5)%raw)
end subroutine collect_containers_caches

end module tblite_post_processing_type
