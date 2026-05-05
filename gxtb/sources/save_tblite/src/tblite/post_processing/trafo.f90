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

!> @file tblite/post_processing/trafo.f90
!> Implements the spherical to cartesian transformation of MO coeffs as post processing method.
module tblite_post_processing_trafo
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_basis_type, only : basis_type
   use tblite_container_cache, only : container_cache
   use tblite_context, only : context_type
   use tblite_double_dictionary, only : double_dictionary_type
   use tblite_integral_type, only : integral_type
   use tblite_integral_trafo, only : adjoint_transform0
   use tblite_output_format, only : format_string
   use tblite_post_processing_type, only : post_processing_type
   use tblite_timer, only : timer_type, format_time
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_xtb_calculator, only : xtb_calculator
   implicit none
   private

   public :: new_sph_cart_trafo, sph_cart_trafo

   !> Spherical to cartesian transformation of MO coeffs as post-processing method
   type, extends(post_processing_type) :: sph_cart_trafo
   contains
      !> Calculate transformed MO coeffs in the Cartesian basis
      procedure :: compute
      !> Print timings
      procedure :: print_timer
   end type sph_cart_trafo

   character(len=*), parameter :: label = "Spherical-to-cartesian MO transformation"

   !type(timer_type) :: timer

contains

subroutine new_sph_cart_trafo(self)
   !> Instance of the spherical to cartesian transformation post-processing
   type(sph_cart_trafo), intent(out) :: self

   self%label = label
end subroutine new_sph_cart_trafo

subroutine compute(self, mol, wfn, integrals, calc, cache_list, ctx, prlevel, dict)
   !> Instance of the spherical to cartesian transformation post-processing
   class(sph_cart_trafo), intent(inout) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction strcuture data
   type(wavefunction_type), intent(in) :: wfn
   !> Integral container
   type(integral_type), intent(in) :: integrals
   !> Calculator instance
   type(xtb_calculator), intent(in) :: calc
   !> Cache list for storing caches of various interactions
   type(container_cache), intent(inout) :: cache_list(:)
   !> Context container for writing to stdout
   type(context_type), intent(inout) :: ctx
   !> Print level
   integer, intent(in) :: prlevel
   !> Dictionary for storing results
   type(double_dictionary_type), intent(inout) :: dict

   integer :: nspin, ispin, iat, jat, izp, jzp, ish, jsh, is, js, li, lj
   integer :: ii, jj, ni, nj, iicart, jjcart, nicart, njcart
   real(wp), allocatable :: coeff_cart(:, :, :)

   !call timer%push("total")
   nspin = size(wfn%density, dim=3)

   allocate(coeff_cart(calc%bas%nao_cart, calc%bas%nao, nspin), source=0.0_wp)

   !$omp parallel do default(none) schedule(runtime) collapse(2) &
   !$omp shared(mol, calc, wfn, coeff_cart, nspin) &
   !$omp private(ispin, iat, izp, is, ish, li, ii, ni, iicart, nicart)
   do ispin = 1, nspin
      do iat = 1, mol%nat
         izp = mol%id(iat)
         is = calc%bas%ish_at(iat)

         do ish = 1, calc%bas%nsh_at(iat)
            li = calc%bas%cgto(ish, izp)%raw%ang

            ii = calc%bas%iao_sh(is+ish)
            ni = calc%bas%nao_sh(is+ish)
            iicart = calc%bas%iao_cart_sh(is+ish)
            nicart = calc%bas%nao_cart_sh(is+ish)

            ! Transform all MO coefficients for the current spherical AO shell
            call adjoint_transform0(li, 0, wfn%coeff(ii+1:ii+ni, :, ispin), &
               & coeff_cart(iicart+1:iicart+nicart, :, ispin), &
               & bra=.true., ket=.false.)
         end do
      end do
   end do

   call dict%add_entry("cartesian-mos", coeff_cart)

   !call timer%pop()

end subroutine compute


subroutine print_timer(self, prlevel, ctx)
   !> Instance of the spherical to cartesian transformation post-processing
   class(sph_cart_trafo), intent(in) :: self
   !> Print level
   integer :: prlevel
   !> Context container for writing to stdout
   type(context_type) :: ctx

   real(wp) :: ttime, stime
   integer :: it
   character(len=*), parameter :: labels(*) = [character(len=20):: &
      & ]

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

end module tblite_post_processing_trafo
