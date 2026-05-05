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

!> @dir tblite/external
!> Contains external potentials

!> @file tblite/external/field.f90
!> Provides interaction of charge density with external fields

!> Interaction of Hamiltonian with external fields
module tblite_external_field
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_blas, only : dot
   use tblite_container_cache, only : container_cache
   use tblite_container_type, only : container_type
   use tblite_scf_info, only : scf_info, atom_resolved
   use tblite_scf_potential, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   public :: new_electric_field

   !> Container for contribution from an instantaneous electric field
   type, public, extends(container_type) :: electric_field
      !> Instantaneous electric field
      real(wp) :: efield(3) = 0.0_wp
   contains
      !> Update cache from container
      procedure :: update
      !> Return dependency on density
      procedure :: variable_info
      !> Get electric field energy
      procedure :: get_energy
      !> Get electric field potential
      procedure :: get_potential
      !> Get electric field potential
      procedure :: get_potential_gradient
      !> Get electric field gradient
      procedure :: get_gradient
   end type electric_field


   !> Create electric field container from field vector
   interface electric_field
      module procedure :: create_electric_field
   end interface electric_field


   !> Identifier for container
   character(len=*), parameter :: label = "electric field"


contains


!> Create new electric field container
pure subroutine new_electric_field(self, efield)
   !> Instance of the electric field container
   type(electric_field), intent(out) :: self
   !> Instantaneous electric field
   real(wp), intent(in) :: efield(:)

   self%efield(:) = efield(:3)
   self%label = label
end subroutine new_electric_field


!> Create new electric field container
pure function create_electric_field(efield) result(new)
   !> Instantaneous electric field
   real(wp), intent(in) :: efield(:)
   !> Instance of the electric field container
   type(electric_field) :: new

   call new_electric_field(new, efield)
end function create_electric_field


!> Update cache from container
subroutine update(self, mol, cache)
   !> Instance of the electric field container
   class(electric_field), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
end subroutine update


!> Get electric field energy
subroutine get_energy(self, mol, cache, wfn, energies)
   !> Instance of the electric field container
   class(electric_field), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Electric field energy
   real(wp), intent(inout) :: energies(:)

   real(wp), allocatable :: vdp(:, :), vat(:)

   vat = matmul(self%efield, mol%xyz)
   vdp = spread(self%efield, 2, mol%nat)
   energies(:) = energies - vat * wfn%qat(:, 1) - sum(vdp * wfn%dpat(:, :, 1), 1)
end subroutine get_energy


!> Get electric field potential
subroutine get_potential(self, mol, cache, wfn, pot)
   !> Instance of the electric field container
   class(electric_field), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential
   type(potential_type), intent(inout) :: pot

   pot%vat(:, 1) = pot%vat(:, 1) - matmul(self%efield, mol%xyz)
   pot%vdp(:, :, 1) = pot%vdp(:, :, 1) - spread(self%efield, 2, mol%nat)
end subroutine get_potential


!> Get electric field potential gradient (empty)
subroutine get_potential_gradient(self, mol, cache, wfn, pot)
   !> Instance of the electric field container
   class(electric_field), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential
   type(potential_type), intent(inout) :: pot

end subroutine get_potential_gradient


!> Get electric field gradient
subroutine get_gradient(self, mol, cache, wfn, gradient, sigma)
   !> Instance of the electric field container
   class(electric_field), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Molecular gradient of the repulsion energy
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Strain derivatives of the repulsion energy
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   real(wp), allocatable :: vdp(:, :), stmp(:, :)

   vdp = spread(self%efield, 2, mol%nat)
   stmp = matmul(vdp, transpose(mol%xyz))
   gradient(:, :) = gradient - vdp
   sigma(:, :) = sigma - 0.5_wp * (stmp + transpose(stmp))
end subroutine get_gradient


!> Return dependency on density
pure function variable_info(self) result(info)
   !> Instance of the electric field container
   class(electric_field), intent(in) :: self
   !> Information on the required potential data
   type(scf_info) :: info

   info = scf_info(charge=atom_resolved, dipole=atom_resolved)
end function variable_info


end module tblite_external_field
