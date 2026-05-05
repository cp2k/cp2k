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

!> @file tblite/spin.f90
!> Provides spin-interactions for spin-polarized tight-binding calculations.

!> Implementation of spin-interactions for spin-polarized tight-binding.
!>
!> The spin-components are represented as difference between alpha and beta population
!> resulting in a magnetization density.
module tblite_spin
   use mctc_env, only : wp
   use mctc_io, only : structure_type, new
   use tblite_blas, only : dot
   use tblite_container_cache, only : container_cache
   use tblite_container_type, only : container_type
   use tblite_scf_info, only : scf_info, shell_resolved
   use tblite_scf_potential, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   public :: spin_polarization, new_spin_polarization


   !> On-site spin-interaction
   type, extends(container_type) :: spin_polarization
      !> Number of shells for each atom
      integer, allocatable :: nsh_at(:)
      !> Index offset for each shell block
      integer, allocatable :: ish_at(:)
      !> Spin constants for each element
      real(wp), allocatable :: wll(:, :, :)
   contains
      !> Return dependency on density
      procedure :: variable_info
      !> Get spin-polarization energy
      procedure :: get_energy
      !> Get spin-polarization potential
      procedure :: get_potential
   end type spin_polarization

   !> Identifier for container
   character(len=*), parameter :: label = "spin polarization"

contains


!> Create new spin-polarization container
pure subroutine new_spin_polarization(self, mol, wll, nshell)
   !> Instance of the spin-polarization container
   type(spin_polarization), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Spin constants for each species
   real(wp), intent(in) :: wll(:, :, :)
   !> Number of shells for each species
   integer, intent(in), optional :: nshell(:)

   integer :: ind, iat

   self%label = label
   self%wll = wll
   self%nsh_at = nshell(mol%id)

   allocate(self%ish_at(mol%nat))
   ind = 0
   do iat = 1, mol%nat
      self%ish_at(iat) = ind
      ind = ind + self%nsh_at(iat)
   end do
end subroutine new_spin_polarization


!> Get spin-polarization energy
subroutine get_energy(self, mol, cache, wfn, energies)
   !> Instance of the spin-polarization container
   class(spin_polarization), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> spin-polarization energy
   real(wp), intent(inout) :: energies(:)

   integer :: iat, izp, ish, jsh, ii, spin

   do spin = 2, wfn%nspin
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ii = self%ish_at(iat)
         do ish = 1, self%nsh_at(iat)
            do jsh = 1, self%nsh_at(iat)
               energies(iat) = energies(iat) + &
                  & 0.5_wp*wfn%qsh(ii+ish, spin)*wfn%qsh(ii+jsh, spin)*self%wll(jsh, ish, izp)
            end do
         end do
      end do
   end do
end subroutine get_energy


!> Get spin-polarization potential
subroutine get_potential(self, mol, cache, wfn, pot)
   !> Instance of the spin-polarization container
   class(spin_polarization), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential
   type(potential_type), intent(inout) :: pot

   integer :: iat, izp, ish, jsh, ii, spin

   do spin = 2, wfn%nspin
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ii = self%ish_at(iat)
         do ish = 1, self%nsh_at(iat)
            do jsh = 1, self%nsh_at(iat)
               pot%vsh(ii+ish, spin) = pot%vsh(ii+ish, spin) + &
                  & wfn%qsh(ii+jsh, spin) * self%wll(jsh, ish, izp)
            end do
         end do
      end do
   end do
end subroutine get_potential


!> Return dependency on density
pure function variable_info(self) result(info)
   !> Instance of the spin-polarization container
   class(spin_polarization), intent(in) :: self
   !> Information on the required potential data
   type(scf_info) :: info

   info = scf_info(charge=shell_resolved)
end function variable_info


end module tblite_spin
