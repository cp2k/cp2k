! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! Lesser GNU General Public License for more details.
!
! You should have received a copy of the Lesser GNU General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

!> High-level wrapper to obtain the dispersion energy for a DFT-D4 calculation
module dftd4_disp
   use, intrinsic :: iso_fortran_env, only : error_unit
   use dftd4_blas, only : d4_gemv
   use dftd4_cache, only : dispersion_cache
   use dftd4_cutoff, only : realspace_cutoff, get_lattice_points
   use dftd4_damping, only : damping_type
   use dftd4_data, only : get_covalent_rad
   use dftd4_model, only : dispersion_model
   use dftd4_ncoord, only : get_coordination_number, add_coordination_number_derivs
   use dftd4_param, only : param_type
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type
   use mctc_io_convert, only : autoaa
   use multicharge, only : get_charges
   implicit none
   private

   public :: get_dispersion, get_properties, get_pairwise_dispersion


contains


!> Wrapper to handle the evaluation of dispersion energy and derivatives
subroutine get_dispersion(mol, disp, damp, param, cutoff, energy, gradient, sigma)
   !DEC$ ATTRIBUTES DLLEXPORT :: get_dispersion

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Dispersion model
   class(dispersion_model), intent(in) :: disp

   !> Damping function
   type(damping_type), intent(in) :: damp

   !> Damping parameters
   type(param_type), intent(in) :: param

   !> Realspace cutoffs
   type(realspace_cutoff), intent(in) :: cutoff

   !> Dispersion energy
   real(wp), intent(out) :: energy

   !> Dispersion gradient
   real(wp), intent(out), contiguous, optional :: gradient(:, :)

   !> Dispersion virial
   real(wp), intent(out), contiguous, optional :: sigma(:, :)
 
   logical :: grad
   integer :: mref
   real(wp), allocatable :: cn(:)
   real(wp), allocatable :: q(:), dqdr(:, :, :), dqdL(:, :, :)
   real(wp), allocatable :: dEdcn(:), dEdq(:), energies(:)
   real(wp), allocatable :: lattr(:, :)
   type(dispersion_cache) :: cache
   type(error_type), allocatable :: error

   mref = maxval(disp%ref)
   grad = present(gradient).and.present(sigma)

   if (.not. allocated(disp%mchrg)) then
      write(error_unit, '("[Error]:", 1x, a)') "Not supported for self-consistent DFT-D"
      error stop
   end if

   allocate(cn(mol%nat))
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
   call get_coordination_number(mol, lattr, cutoff%cn, disp%rcov, disp%en, cn)

   allocate(q(mol%nat))
   if (grad) allocate(dqdr(3, mol%nat, mol%nat), dqdL(3, 3, mol%nat))
   call get_charges(disp%mchrg, mol, error, q, dqdr, dqdL)
   if(allocated(error)) then
      write(error_unit, '("[Error]:", 1x, a)') error%message
      error stop
   end if

   call disp%update(mol, cache, cn, q, grad=grad)

   allocate(energies(mol%nat))
   energies(:) = 0.0_wp
   if (grad) then
      allocate(dEdcn(mol%nat), dEdq(mol%nat))
      dEdcn(:) = 0.0_wp
      dEdq(:) = 0.0_wp
      gradient(:, :) = 0.0_wp
      sigma(:, :) = 0.0_wp
   end if

   call get_lattice_points(mol%periodic, mol%lattice, cutoff%disp2, lattr)
   call disp%get_dispersion2(mol, cache, damp, param, lattr, cutoff%disp2, &
      & energies, dEdcn, dEdq, gradient, sigma)
   
   if (.not.disp%charge_dep_3b) then
      q(:) = 0.0_wp
      if (grad) then
         call d4_gemv(dqdr, dEdq, gradient, beta=1.0_wp)
         call d4_gemv(dqdL, dEdq, sigma, beta=1.0_wp)
      end if
   end if

   call disp%update(mol, cache, cn, q, grad=grad)

   call get_lattice_points(mol%periodic, mol%lattice, cutoff%disp3, lattr)
   call disp%get_dispersion3(mol, cache, damp, param, lattr, cutoff%disp3, &
      & energies, dEdcn, dEdq, gradient, sigma)
   if (grad) then
      call add_coordination_number_derivs(mol, lattr, cutoff%cn, &
         & disp%rcov, disp%en, dEdcn, gradient, sigma)
   end if

   if (disp%charge_dep_3b) then
      if (grad) then
         call d4_gemv(dqdr, dEdq, gradient, beta=1.0_wp)
         call d4_gemv(dqdL, dEdq, sigma, beta=1.0_wp)
      end if
   end if

   energy = sum(energies)

end subroutine get_dispersion


!> Wrapper to handle the evaluation of properties related to this dispersion model
subroutine get_properties(mol, disp, cutoff, cn, q, c6, alpha, alphaqq)
   !DEC$ ATTRIBUTES DLLEXPORT :: get_properties

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Dispersion model
   class(dispersion_model), intent(in) :: disp

   !> Realspace cutoffs
   type(realspace_cutoff), intent(in) :: cutoff

   !> Coordination number
   real(wp), intent(out) :: cn(:)

   !> Atomic partial charges
   real(wp), intent(out), contiguous :: q(:)

   !> C6 coefficients
   real(wp), intent(out) :: c6(:, :)

   !> Static dipole-dipole polarizabilities
   real(wp), intent(out) :: alpha(:)

   !> Static quadrupole-quadrupole polarizabilities
   real(wp), intent(out) :: alphaqq(:)

   integer :: mref
   real(wp), allocatable :: lattr(:, :)
   type(dispersion_cache) :: cache
   type(error_type), allocatable :: error

   if (.not. allocated(disp%mchrg)) then
      write(error_unit, '("[Error]:", 1x, a)') "Not supported for self-consistent DFT-D"
      error stop
   end if

   mref = maxval(disp%ref)

   call get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
   call get_coordination_number(mol, lattr, cutoff%cn, disp%rcov, disp%en, cn)

   call get_charges(disp%mchrg, mol, error, q)
   if(allocated(error)) then
      write(error_unit, '("[Error]:", 1x, a)') error%message
      error stop
   end if

   call disp%update(mol, cache, cn, q, grad=.false.)
   c6(:, :) = cache%c6

   call disp%get_polarizabilities(cache, alpha, alphaqq)

end subroutine get_properties


!> Wrapper to handle the evaluation of pairwise representation of the dispersion energy
subroutine get_pairwise_dispersion(mol, disp, damp, param, cutoff, energy2, energy3)
   !DEC$ ATTRIBUTES DLLEXPORT :: get_pairwise_dispersion

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Dispersion model
   class(dispersion_model), intent(in) :: disp

   !> Damping function
   type(damping_type), intent(in) :: damp

   !> Damping parameters
   type(param_type), intent(in) :: param

   !> Realspace cutoffs
   type(realspace_cutoff), intent(in) :: cutoff

   !> Pairwise representation of additive dispersion energy
   real(wp), intent(out) :: energy2(:, :)

   !> Pairwise representation of non-additive dispersion energy
   real(wp), intent(out) :: energy3(:, :)

   integer :: mref
   real(wp), allocatable :: cn(:), q(:), lattr(:, :)
   type(dispersion_cache) :: cache
   type(error_type), allocatable :: error

   if (.not. allocated(disp%mchrg)) then
      write(error_unit, '("[Error]:", 1x, a)') "Not supported for self-consistent DFT-D"
      error stop
   end if

   mref = maxval(disp%ref)

   allocate(cn(mol%nat))
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
   call get_coordination_number(mol, lattr, cutoff%cn, disp%rcov, disp%en, cn)

   allocate(q(mol%nat))
   call get_charges(disp%mchrg, mol, error, q)
   if(allocated(error)) then
      write(error_unit, '("[Error]:", 1x, a)') error%message
      error stop
   end if

   call disp%update(mol, cache, cn, q, grad=.false.)

   energy2(:, :) = 0.0_wp
   energy3(:, :) = 0.0_wp
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%disp2, lattr)
   call disp%get_pairwise_dispersion2(mol, cache, damp%damping_2b, param, &
      & lattr, cutoff%disp2, energy2)

   q(:) = 0.0_wp
   call disp%update(mol, cache, cn, q, grad=.false.)

   call get_lattice_points(mol%periodic, mol%lattice, cutoff%disp3, lattr)
   call disp%get_pairwise_dispersion3(mol, cache, damp%damping_3b, param, &
      & lattr, cutoff%disp3, energy3)

end subroutine get_pairwise_dispersion


end module dftd4_disp
