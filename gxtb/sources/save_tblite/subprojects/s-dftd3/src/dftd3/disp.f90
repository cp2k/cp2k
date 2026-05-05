! This file is part of s-dftd3.
! SPDX-Identifier: LGPL-3.0-or-later
!
! s-dftd3 is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! s-dftd3 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with s-dftd3.  If not, see <https://www.gnu.org/licenses/>.

module dftd3_disp
   use dftd3_cutoff, only : realspace_cutoff, get_lattice_points
   use dftd3_damping, only : damping_param
   use dftd3_model, only : d3_model
   use dftd3_ncoord, only : get_coordination_number, add_coordination_number_derivs
   use mctc_data, only : get_covalent_rad
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_convert, only : autoaa
   implicit none
   private

   public :: get_dispersion, get_pairwise_dispersion


   !> Calculate dispersion energy
   interface get_dispersion
      module procedure :: get_dispersion_atomic
      module procedure :: get_dispersion_scalar
   end interface get_dispersion

contains


!> Calculate atom-resolved dispersion energies.
subroutine get_dispersion_atomic(mol, disp, param, cutoff, energies, gradient, sigma)

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Dispersion model
   class(d3_model), intent(in) :: disp

   !> Damping parameters
   class(damping_param), intent(in) :: param

   !> Realspace cutoffs
   type(realspace_cutoff), intent(in) :: cutoff

   !> Dispersion energy
   real(wp), intent(out) :: energies(:)

   !> Dispersion gradient
   real(wp), intent(out), contiguous, optional :: gradient(:, :)

   !> Dispersion virial
   real(wp), intent(out), contiguous, optional :: sigma(:, :)

   logical :: grad
   integer :: mref
   real(wp), allocatable :: cn(:)
   real(wp), allocatable :: gwvec(:, :), gwdcn(:, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :)
   real(wp), allocatable :: dEdcn(:)
   real(wp), allocatable :: lattr(:, :)

   mref = maxval(disp%ref)
   grad = present(gradient).and.present(sigma)

   allocate(cn(mol%nat))
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
   call get_coordination_number(mol, lattr, cutoff%cn, disp%rcov, cn)

   allocate(gwvec(mref, mol%nat))
   if (grad) allocate(gwdcn(mref, mol%nat))
   call disp%weight_references(mol, cn, gwvec, gwdcn)

   allocate(c6(mol%nat, mol%nat))
   if (grad) allocate(dc6dcn(mol%nat, mol%nat))
   call disp%get_atomic_c6(mol, gwvec, gwdcn, c6, dc6dcn)

   energies(:) = 0.0_wp
   if (grad) then
      allocate(dEdcn(mol%nat))
      dEdcn(:) = 0.0_wp
      gradient(:, :) = 0.0_wp
      sigma(:, :) = 0.0_wp
   end if
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%disp2, lattr)
   call param%get_dispersion2(mol, lattr, cutoff%disp2, disp%rvdw, disp%r4r2, c6, dc6dcn, &
      & energies, dEdcn, gradient, sigma)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%disp3, lattr)
   call param%get_dispersion3(mol, lattr, cutoff%disp3, disp%rvdw, disp%r4r2, c6, dc6dcn, &
      & energies, dEdcn, gradient, sigma)
   if (grad) then
      call get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
      call add_coordination_number_derivs(mol, lattr, cutoff%cn, disp%rcov, dEdcn, &
         & gradient, sigma)
   end if

end subroutine get_dispersion_atomic


!> Calculate scalar dispersion energy.
subroutine get_dispersion_scalar(mol, disp, param, cutoff, energy, gradient, sigma)

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Dispersion model
   class(d3_model), intent(in) :: disp

   !> Damping parameters
   class(damping_param), intent(in) :: param

   !> Realspace cutoffs
   type(realspace_cutoff), intent(in) :: cutoff

   !> Dispersion energy
   real(wp), intent(out) :: energy

   !> Dispersion gradient
   real(wp), intent(out), contiguous, optional :: gradient(:, :)

   !> Dispersion virial
   real(wp), intent(out), contiguous, optional :: sigma(:, :)

   real(wp), allocatable :: energies(:)

   allocate(energies(mol%nat))

   call get_dispersion_atomic(mol, disp, param, cutoff, energies, gradient, sigma)

   energy = sum(energies)

end subroutine get_dispersion_scalar


!> Wrapper to handle the evaluation of pairwise representation of the dispersion energy
subroutine get_pairwise_dispersion(mol, disp, param, cutoff, energy2, energy3)

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Dispersion model
   class(d3_model), intent(in) :: disp

   !> Damping parameters
   class(damping_param), intent(in) :: param

   !> Realspace cutoffs
   type(realspace_cutoff), intent(in) :: cutoff

   !> Pairwise representation of additive dispersion energy
   real(wp), intent(out) :: energy2(:, :)

   !> Pairwise representation of non-additive dispersion energy
   real(wp), intent(out) :: energy3(:, :)

   integer :: mref
   real(wp), allocatable :: cn(:), gwvec(:, :), c6(:, :), lattr(:, :)

   mref = maxval(disp%ref)

   allocate(cn(mol%nat))
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
   call get_coordination_number(mol, lattr, cutoff%cn, disp%rcov, cn)

   allocate(gwvec(mref, mol%nat))
   call disp%weight_references(mol, cn, gwvec)

   allocate(c6(mol%nat, mol%nat))
   call disp%get_atomic_c6(mol, gwvec, c6=c6)

   energy2(:, :) = 0.0_wp
   energy3(:, :) = 0.0_wp
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%disp2, lattr)
   call param%get_pairwise_dispersion2(mol, lattr, cutoff%disp2, disp%rvdw, disp%r4r2, &
      & c6, energy2)

   call get_lattice_points(mol%periodic, mol%lattice, cutoff%disp3, lattr)
   call param%get_pairwise_dispersion3(mol, lattr, cutoff%disp3, disp%rvdw, disp%r4r2, &
      & c6, energy3)

end subroutine get_pairwise_dispersion


end module dftd3_disp
