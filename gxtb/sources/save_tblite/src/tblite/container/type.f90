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

!> @file tblite/container/type.f90
!> Provides an abstract base class for interaction containers

!> Definition of abstract base class for general interactions
module tblite_container_type
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_container_cache, only : container_cache
   use tblite_scf_info, only : scf_info
   use tblite_scf_potential, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   !> Abstract base class for interactions containers
   type, public, abstract :: container_type
      !> Label identifying this contribution
      character(len=:), allocatable :: label
   contains
      !> Update container cache
      procedure :: update
      !> Get information about density dependent quantities used in the energy
      procedure :: variable_info
      !> Evaluate non-selfconsistent part of the interaction
      procedure :: get_engrad
      !> Evaluate selfconsistent energy of the interaction
      procedure :: get_energy
      !> Evaluate selfconsistent and overlap-dependent energy of the interaction
      procedure :: get_energy_w_overlap
      !> Evaluate charge dependent potential shift from the interaction
      procedure :: get_potential
      !> Evaluate charge and overlap dependent potential shift from the interaction
      procedure :: get_potential_w_overlap
      !> Evaluate gradient of charge dependent potential shift from the interaction
      procedure :: get_potential_gradient
      !> Evaluate gradient contributions from the selfconsistent interaction
      procedure :: get_gradient
      !> Evaluate gradient contributions with overlap dependence from the selfconsistent interaction
      procedure :: get_gradient_w_overlap
      !> Information on container
      procedure :: info
   end type container_type


contains


!> Update container cache
subroutine update(self, mol, cache)
   !> Instance of the interaction container
   class(container_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different runs
   type(container_cache), intent(inout) :: cache
end subroutine update


!> Evaluate non-selfconsistent part of the interaction
subroutine get_engrad(self, mol, cache, energies, gradient, sigma)
   !> Instance of the interaction container
   class(container_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different runs
   type(container_cache), intent(inout) :: cache
   !> Interaction energy
   real(wp), intent(inout) :: energies(:)
   !> Interaction gradient
   real(wp), contiguous, intent(inout), optional :: gradient(:, :)
   !> Interaction virial
   real(wp), contiguous, intent(inout), optional :: sigma(:, :)
end subroutine get_engrad


!> Evaluate selfconsistent energy of the interaction
subroutine get_energy(self, mol, cache, wfn, energies)
   !> Instance of the interaction container
   class(container_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different runs
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Interaction energy
   real(wp), intent(inout) :: energies(:)
end subroutine get_energy


!> Evaluate selfconsistent and overlap-dependent energy of the interaction
subroutine get_energy_w_overlap(self, mol, cache, wfn, overlap, energies)
   !> Instance of the interaction container
   class(container_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different runs
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Overlap integral matrix
   real(wp), contiguous, intent(in) :: overlap(:, :)
   !> Interaction energy
   real(wp), intent(inout) :: energies(:)
end subroutine get_energy_w_overlap


!> Evaluate charge dependent potential shift from the interaction
subroutine get_potential(self, mol, cache, wfn, pot)
   !> Instance of the interaction container
   class(container_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different runs
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential
   type(potential_type), intent(inout) :: pot
end subroutine get_potential


!> Evaluate charge and overlap dependent potential shift from the interaction
subroutine get_potential_w_overlap(self, mol, cache, wfn, overlap, pot)
   !> Instance of the interaction container
   class(container_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different runs
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Overlap integral matrix
   real(wp), contiguous, intent(in) :: overlap(:, :)
   !> Density dependent potential
   type(potential_type), intent(inout) :: pot
end subroutine get_potential_w_overlap


!> Evaluate charge dependent potential shift from the interaction
subroutine get_potential_gradient(self, mol, cache, wfn, pot)
   !> Instance of the interaction container
   class(container_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different runs
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential
   type(potential_type), intent(inout) :: pot
end subroutine get_potential_gradient


!> Evaluate gradient contributions from the selfconsistent interaction
subroutine get_gradient(self, mol, cache, wfn, gradient, sigma)
   !> Instance of the interaction container
   class(container_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different runs
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Interaction gradient
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Interaction virial
   real(wp), contiguous, intent(inout) :: sigma(:, :)
end subroutine get_gradient


!> Evaluate gradient contributions with overlap dependence from the selfconsistent interaction
subroutine get_gradient_w_overlap(self, mol, cache, wfn, overlap, &
   & ao_grad, gradient, sigma)
   !> Instance of the interaction container
   class(container_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different runs
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Overlap integral matrix
   real(wp), contiguous, intent(in) :: overlap(:, :)
   !> Orbital gradient contribution to the energy-weighted density matrix
   real(wp), contiguous, intent(inout) :: ao_grad(:, :)
   !> Molecular gradient of the interaction energy
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Strain derivatives of the interaction energy
   real(wp), contiguous, intent(inout) :: sigma(:, :)
end subroutine get_gradient_w_overlap


!> Get information about density dependent quantities used in the energy
pure function variable_info(self) result(info)
   !> Instance of the interaction container
   class(container_type), intent(in) :: self
   !> Information on the required potential data
   type(scf_info) :: info

   info = scf_info()
end function variable_info


!> Information on container
pure function info(self, verbosity, indent) result(str)
   !> Instance of the interaction container
   class(container_type), intent(in) :: self
   !> Verbosity level
   integer, intent(in) :: verbosity
   !> Indentation level
   character(len=*), intent(in) :: indent
   !> Information on the container
   character(len=:), allocatable :: str

   if (allocated(self%label)) then
      str = self%label
   else
      str = "Unknown"
   end if
end function info


end module tblite_container_type
