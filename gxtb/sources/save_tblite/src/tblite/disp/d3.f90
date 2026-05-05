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

!> @file tblite/disp/d3.f90
!> Provides a proxy for the [DFT-D3 dispersion correction](https://dftd3.readthedocs.io)

!> Semiclassical DFT-D3 dispersion correction
module tblite_disp_d3
   use mctc_env, only : error_type, wp
   use mctc_io, only : structure_type
   use mctc_ncoord, only : new_ncoord, ncoord_type, cn_count
   use tblite_container_cache, only : container_cache
   use tblite_cutoff, only : get_lattice_points
   use tblite_disp_cache, only : dispersion_cache
   use tblite_disp_type, only : dispersion_type
   use dftd3, only : d3_model, new_d3_model, rational_damping_param, realspace_cutoff
   implicit none
   private

   public :: new_d3_dispersion


   !> Container for DFT-D3 type dispersion correction
   type, public, extends(dispersion_type) :: d3_dispersion
      !> Instance of the actual D3 dispersion model
      type(d3_model) :: model
      !> Rational damping parameters
      type(rational_damping_param) :: param
      !> Selected real space cutoffs for this instance
      type(realspace_cutoff) :: cutoff
       !> Coordination number instance
      class(ncoord_type), allocatable :: ncoord
   contains
      procedure :: get_engrad
   end type d3_dispersion

   character(len=*), parameter :: label = "DFT-D3(BJ) dispersion correction"

contains

subroutine new_d3_dispersion(self, mol, s6, s8, a1, a2, s9, error)
   !> Instance of the dispersion correction
   type(d3_dispersion), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Damping parameters
   real(wp), intent(in) :: s6, s8, a1, a2, s9
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   self%label = label
   call new_d3_model(self%model, mol)
   self%param = rational_damping_param(s6=s6, s8=s8, s9=s9, a1=a1, a2=a2, alp=14.0_wp)
   self%cutoff = realspace_cutoff(cn=25.0_wp, disp3=25.0_wp, disp2=50.0_wp)
   call new_ncoord(self%ncoord, mol, cn_count%exp, error, &
      & cutoff=self%cutoff%cn, rcov=self%model%rcov)
end subroutine new_d3_dispersion

!> Evaluate non-selfconsistent part of the dispersion correction
subroutine get_engrad(self, mol, cache, energies, gradient, sigma)
   !> Instance of the dispersion correction
   class(d3_dispersion), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different dispersion runs
   type(container_cache), intent(inout) :: cache
   !> Dispersion energy
   real(wp), intent(inout) :: energies(:)
   !> Dispersion gradient
   real(wp), contiguous, intent(inout), optional :: gradient(:, :)
   !> Dispersion virial
   real(wp), contiguous, intent(inout), optional :: sigma(:, :)

   logical :: grad
   integer :: mref
   real(wp), allocatable :: cn(:)
   real(wp), allocatable :: gwvec(:, :), gwdcn(:, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :)
   real(wp), allocatable :: dEdcn(:), lattr(:, :)
   class(ncoord_type), allocatable :: ncoord

   mref = maxval(self%model%ref)
   grad = present(gradient).and.present(sigma)

   allocate(cn(mol%nat))
   call get_lattice_points(mol%periodic, mol%lattice, self%cutoff%cn, lattr)
   call self%ncoord%get_coordination_number(mol, lattr, cn)

   allocate(gwvec(mref, mol%nat))
   if (grad) allocate(gwdcn(mref, mol%nat))
   call self%model%weight_references(mol, cn, gwvec, gwdcn)

   allocate(c6(mol%nat, mol%nat))
   if (grad) allocate(dc6dcn(mol%nat, mol%nat))
   call self%model%get_atomic_c6(mol, gwvec, gwdcn, c6, dc6dcn)

   if (grad) then
      allocate(dEdcn(mol%nat))
      dEdcn(:) = 0.0_wp
   end if
   call get_lattice_points(mol%periodic, mol%lattice, self%cutoff%disp2, lattr)
   call self%param%get_dispersion2(mol, lattr, self%cutoff%disp2, self%model%rvdw, &
      & self%model%r4r2, c6, dc6dcn, energies, dEdcn, gradient, sigma)

   call get_lattice_points(mol%periodic, mol%lattice, self%cutoff%disp3, lattr)
   call self%param%get_dispersion3(mol, lattr, self%cutoff%disp3, self%model%rvdw, &
      & self%model%r4r2, c6, dc6dcn, energies, dEdcn, gradient, sigma)
   if (grad) then
      call self%ncoord%add_coordination_number_derivs(mol, lattr, dEdcn, gradient, sigma)
   end if
end subroutine get_engrad


end module tblite_disp_d3
