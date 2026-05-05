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

!> @file tblite/repulsion/gfn.f90
!> Provides a screened Coulomb repulsion interaction

!> Classical repulsion interaction as used with the GFNn-xTB Hamiltonian
module tblite_repulsion_gfn
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_container, only : container_cache
   use tblite_cutoff, only : get_lattice_points
   use tblite_repulsion_cache, only : repulsion_cache
   use tblite_repulsion_type, only : repulsion_type
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   public :: new_repulsion_gfn


   !> Container to evaluate classical repulsion interactions for the GFN Hamiltonian
   type, public, extends(repulsion_type) :: gfn_repulsion
      !> Exponent for the damping of the repulsion
      real(wp), allocatable :: alpha(:, :)
      !> Effective nuclear charge
      real(wp), allocatable :: zeff(:, :)
      !> Exponent of the repulsion polynomial
      real(wp), allocatable :: rexp(:, :)
      !> Distance exponent for repulsion damping
      real(wp), allocatable :: kexp(:, :)
   contains
      !> Evaluate non-selfconsistent repulsion energy and gradient
      procedure :: get_engrad
   end type gfn_repulsion

   character(len=*), parameter :: label = "screened Coulomb repulsion"

contains


subroutine new_repulsion_gfn(self, mol, alpha, zeff, kexp, kexp_light, rexp, cutoff)
   !> Instance of the repulsion container
   type(gfn_repulsion), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Exponent for the damping of the repulsion
   real(wp), intent(in) :: alpha(:)
   !> Effective nuclear charge
   real(wp), intent(in) :: zeff(:)
   !> Distance exponent for repulsion damping
   real(wp), intent(in) :: kexp
   !> Distance exponent for repulsion damping of light atom pairs
   real(wp), intent(in) :: kexp_light
   !> Exponent of the repulsion polynomial
   real(wp), intent(in) :: rexp
   !> Real-space cutoff
   real(wp), intent(in), optional :: cutoff

   integer :: isp, izp, jsp, jzp

   self%label = label
   allocate(self%alpha(mol%nid, mol%nid))
   do isp = 1, mol%nid
      do jsp = 1, mol%nid
         self%alpha(jsp, isp) = sqrt(alpha(isp)*alpha(jsp))
      end do
   end do

   allocate(self%zeff(mol%nid, mol%nid))
   do isp = 1, mol%nid
      do jsp = 1, mol%nid
         self%zeff(jsp, isp) = zeff(isp)*zeff(jsp)
      end do
   end do

   allocate(self%kexp(mol%nid, mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do jsp = 1, mol%nid
         jzp = mol%num(jsp)
         self%kexp(jsp, isp) = merge(kexp, kexp_light, izp > 2 .or. jzp > 2)
      end do
   end do

   allocate(self%rexp(mol%nid, mol%nid))
   self%rexp(:, :) = rexp
   if (present(cutoff)) self%cutoff = cutoff

end subroutine new_repulsion_gfn


!> Evaluate classical interaction for energy and derivatives
subroutine get_engrad(self, mol, cache, energies, gradient, sigma)
   !> Instance of the repulsion container
   class(gfn_repulsion), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different runs
   type(container_cache), intent(inout) :: cache
   !> Repulsion energy
   real(wp), intent(inout) :: energies(:)
   !> Molecular gradient of the repulsion energy
   real(wp), contiguous, intent(inout), optional :: gradient(:, :)
   !> Strain derivatives of the repulsion energy
   real(wp), contiguous, intent(inout), optional :: sigma(:, :)

   real(wp), allocatable :: trans(:, :)

   call get_lattice_points(mol%periodic, mol%lattice, self%cutoff, trans)

   if (present(gradient) .and. present(sigma)) then
      call get_repulsion_derivs(mol, trans, self%cutoff, self%alpha, self%zeff, &
         & self%kexp, self%rexp, energies, gradient, sigma)
   else
      call get_repulsion_energy(mol, trans, self%cutoff, self%alpha, self%zeff, &
         & self%kexp, self%rexp, energies)
   end if

end subroutine get_engrad


subroutine get_repulsion_energy(mol, trans, cutoff, alpha, zeff, kexp, rexp, &
   & energies)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points
   real(wp), intent(in) :: trans(:, :)
   !> Real space cutoff
   real(wp), intent(in) :: cutoff
   !> Exponents for all element pairs
   real(wp), intent(in) :: alpha(:, :)
   !> Effective nuclear charges for all element pairs
   real(wp), intent(in) :: zeff(:, :)
   !> Pairwise parameters for all element pairs
   real(wp), intent(in) :: kexp(:, :)
   !> Pairwise parameters for all element pairs
   real(wp), intent(in) :: rexp(:, :)
   !> Repulsion energy
   real(wp), intent(inout) :: energies(:)

   integer :: iat, jat, izp, jzp, itr
   real(wp) :: r1, r2, rij(3), r1k, r1r, exa, cutoff2, dE

   cutoff2 = cutoff**2

   !$omp parallel do default(none) schedule(runtime) reduction(+:energies) &
   !$omp shared(mol, trans, cutoff2, kexp, rexp, alpha, zeff) &
   !$omp private(iat, jat, izp, jzp, itr, r1, r2, rij, r1k, r1r, exa, dE)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)
         do itr = 1, size(trans, dim=2)
            rij = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
            r2 = sum(rij**2)
            if (r2 > cutoff2 .or. r2 < 1.0e-12_wp) cycle
            r1 = sqrt(r2)

            r1k = r1**kexp(jzp, izp)
            exa = exp(-alpha(jzp, izp)*r1k)
            r1r = r1**rexp(jzp, izp)
            dE = zeff(jzp, izp) * exa/r1r
            energies(iat) = energies(iat) + 0.5_wp * dE
            if (iat /= jat) then
               energies(jat) = energies(jat) + 0.5_wp * dE
            end if
         end do
      end do
   end do

end subroutine get_repulsion_energy


subroutine get_repulsion_derivs(mol, trans, cutoff, alpha, zeff, kexp, rexp, &
   & energies, gradient, sigma)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points
   real(wp), intent(in) :: trans(:, :)
   !> Real space cutoff
   real(wp), intent(in) :: cutoff
   !> Exponents for all element pairs
   real(wp), intent(in) :: alpha(:, :)
   !> Effective nuclear charges for all element pairs
   real(wp), intent(in) :: zeff(:, :)
   !> Pairwise parameters for all element pairs
   real(wp), intent(in) :: kexp(:, :)
   !> Pairwise parameters for all element pairs
   real(wp), intent(in) :: rexp(:, :)
   !> Repulsion energy
   real(wp), intent(inout) :: energies(:)
   !> Molecular gradient of the repulsion energy
   real(wp), intent(inout) :: gradient(:, :)
   !> Strain derivatives of the repulsion energy
   real(wp), intent(inout) :: sigma(:, :)

   integer :: iat, jat, izp, jzp, itr
   real(wp) :: r1, r2, rij(3), r1k, r1r, exa, cutoff2, dE, dG(3), dS(3, 3)

   cutoff2 = cutoff**2

   !$omp parallel do default(none) schedule(runtime) reduction(+:energies, gradient, sigma) &
   !$omp shared(mol, trans, cutoff2, kexp, rexp, alpha, zeff) &
   !$omp private(iat, jat, izp, jzp, itr, r1, r2, rij, r1k, r1r, exa, dE, dG, dS)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)
         do itr = 1, size(trans, dim=2)
            rij = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
            r2 = sum(rij**2)
            if (r2 > cutoff2 .or. r2 < 1.0e-12_wp) cycle
            r1 = sqrt(r2)

            r1k = r1**kexp(jzp, izp)
            exa = exp(-alpha(jzp, izp)*r1k)
            r1r = r1**rexp(jzp, izp)
            dE = zeff(jzp, izp) * exa/r1r
            dG = -(alpha(jzp, izp)*r1k*kexp(jzp, izp) + rexp(jzp, izp)) * dE * rij/r2
            dS = spread(dG, 1, 3) * spread(rij, 2, 3)
            energies(iat) = energies(iat) + 0.5_wp * dE
            if (iat /= jat) then
               energies(jat) = energies(jat) + 0.5_wp * dE
               gradient(:, iat) = gradient(:, iat) + dG
               gradient(:, jat) = gradient(:, jat) - dG
               sigma(:, :) = sigma + dS
            else
               sigma(:, :) = sigma + 0.5_wp * dS
            end if
         end do
      end do
   end do

end subroutine get_repulsion_derivs


end module tblite_repulsion_gfn
