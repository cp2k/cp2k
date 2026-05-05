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

!> @file tblite/container/list.f90
!> Provides a collection for several interaction containers

!> Definition of list of general interaction contaienrs
module tblite_container_list
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_container_cache, only : container_cache, resize
   use tblite_container_type, only : container_type
   use tblite_scf_info, only : scf_info
   use tblite_scf_potential, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_output_format, only : format_string
   implicit none
   private


   !> Wrapped container type for creation of lists
   type :: container_node
      !> Actual interaction container
      class(container_type), allocatable :: raw
   end type container_node

   !> List of interaction containers
   type, public, extends(container_type) :: container_list
      private
      !> Number of stored containers
      integer :: nc = 0
      !> Raw list of interaction containers
      type(container_node), allocatable :: list(:)
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
      !> Add a container
      procedure :: push_back
      !> Remove a container
      procedure :: pop
      !> Information about the container list
      procedure :: info
   end type container_list


   !> Reallocate list of containers
   interface resize
      module procedure :: resize_node
   end interface


   !> List of container cache instances
   type, public :: cache_list
      !> Actual cache instances
      type(container_cache), allocatable :: list(:)
   end type cache_list


contains


!> Update container cache
subroutine update(self, mol, cache)
   !> Instance of the interaction container
   class(container_list), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different runs
   type(container_cache), intent(inout) :: cache

   integer :: ic
   type(cache_list), pointer :: ptr

   call taint(cache, ptr)
   call resize(ptr%list, self%nc)

   do ic = 1, self%nc
      if (allocated(self%list(ic)%raw)) then
         associate(cont => self%list(ic)%raw)
            call cont%update(mol, ptr%list(ic))
         end associate
      end if
   end do
end subroutine update


!> Evaluate non-selfconsistent part of the interaction
subroutine get_engrad(self, mol, cache, energies, gradient, sigma)
   !> Instance of the interaction container
   class(container_list), intent(in) :: self
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

   integer :: ic
   type(cache_list), pointer :: ptr

   call view(cache, ptr)

   do ic = 1, self%nc
      if (allocated(self%list(ic)%raw)) then
         associate(cont => self%list(ic)%raw)
            call cont%get_engrad(mol, ptr%list(ic), energies, gradient, sigma)
         end associate
      end if
   end do
end subroutine get_engrad


!> Evaluate selfconsistent energy of the interaction
subroutine get_energy(self, mol, cache, wfn, energies)
   !> Instance of the interaction container
   class(container_list), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different runs
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Interaction energy
   real(wp), intent(inout) :: energies(:)

   integer :: ic
   type(cache_list), pointer :: ptr

   call view(cache, ptr)

   do ic = 1, self%nc
      if (allocated(self%list(ic)%raw)) then
         associate(cont => self%list(ic)%raw)
            call cont%get_energy(mol, ptr%list(ic), wfn, energies)
         end associate
      end if
   end do
end subroutine get_energy


!> Evaluate selfconsistent and overlap-dependent energy of the interaction
subroutine get_energy_w_overlap(self, mol, cache, wfn, overlap, energies)
   !> Instance of the interaction container
   class(container_list), intent(in) :: self
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

   integer :: ic
   type(cache_list), pointer :: ptr

   call view(cache, ptr)

   do ic = 1, self%nc
      if (allocated(self%list(ic)%raw)) then
         associate(cont => self%list(ic)%raw)
            call cont%get_energy_w_overlap(mol, ptr%list(ic), wfn, overlap, energies)
         end associate
      end if
   end do
end subroutine get_energy_w_overlap


!> Evaluate charge dependent potential shift from the interaction
subroutine get_potential(self, mol, cache, wfn, pot)
   !> Instance of the interaction container
   class(container_list), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different runs
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential
   type(potential_type), intent(inout) :: pot

   integer :: ic
   type(cache_list), pointer :: ptr

   call view(cache, ptr)

   do ic = 1, self%nc
      if (allocated(self%list(ic)%raw)) then
         associate(cont => self%list(ic)%raw)
            call cont%get_potential(mol, ptr%list(ic), wfn, pot)
         end associate
      end if
   end do
end subroutine get_potential


!> Evaluate charge and overlap dependent potential shift from the interaction
subroutine get_potential_w_overlap(self, mol, cache, wfn, overlap, pot)
   !> Instance of the interaction container
   class(container_list), intent(in) :: self
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

   integer :: ic
   type(cache_list), pointer :: ptr

   call view(cache, ptr)

   do ic = 1, self%nc
      if (allocated(self%list(ic)%raw)) then
         associate(cont => self%list(ic)%raw)
            call cont%get_potential_w_overlap(mol, ptr%list(ic), wfn, overlap, pot)
         end associate
      end if
   end do
end subroutine get_potential_w_overlap


!> Evaluate charge dependent potential shift from the interaction
subroutine get_potential_gradient(self, mol, cache, wfn, pot)
   !> Instance of the interaction container
   class(container_list), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different runs
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential
   type(potential_type), intent(inout) :: pot

   integer :: ic
   type(cache_list), pointer :: ptr

   call view(cache, ptr)

   do ic = 1, self%nc
      if (allocated(self%list(ic)%raw)) then
         associate(cont => self%list(ic)%raw)
            call cont%get_potential_gradient(mol, ptr%list(ic), wfn, pot)
         end associate
      end if
   end do
end subroutine get_potential_gradient


!> Evaluate gradient contributions from the selfconsistent interaction
subroutine get_gradient(self, mol, cache, wfn, gradient, sigma)
   !> Instance of the interaction container
   class(container_list), intent(in) :: self
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

   integer :: ic
   type(cache_list), pointer :: ptr

   call view(cache, ptr)

   do ic = 1, self%nc
      if (allocated(self%list(ic)%raw)) then
         associate(cont => self%list(ic)%raw)
            call cont%get_gradient(mol, ptr%list(ic), wfn, gradient, sigma)
         end associate
      end if
   end do
end subroutine get_gradient



!> Evaluate gradient contributions with overlap dependence from the selfconsistent interaction
subroutine get_gradient_w_overlap(self, mol, cache, wfn, overlap, &
   & ao_grad, gradient, sigma)
   !> Instance of the interaction container
   class(container_list), intent(in) :: self
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
   !> Interaction gradient
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Interaction virial
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   integer :: ic
   type(cache_list), pointer :: ptr

   call view(cache, ptr)

   do ic = 1, self%nc
      if (allocated(self%list(ic)%raw)) then
         associate(cont => self%list(ic)%raw)
            call cont%get_gradient_w_overlap(mol, ptr%list(ic), wfn, overlap, &
               & ao_grad, gradient, sigma)
         end associate
      end if
   end do
end subroutine get_gradient_w_overlap


!> Get information about density dependent quantities used in the energy
pure function variable_info(self) result(info)
   use tblite_scf_info, only : scf_info, max
   !> Instance of the interaction container
   class(container_list), intent(in) :: self
   !> Information on the required potential data
   type(scf_info) :: info

   integer :: ic

   info = scf_info()

   do ic = 1, self%nc
      if (allocated(self%list(ic)%raw)) then
         info = max(info, self%list(ic)%raw%variable_info())
      end if
   end do
end function variable_info


!> Information on container
pure function info(self, verbosity, indent) result(str)
   !> Instance of the interaction container
   class(container_list), intent(in) :: self
   !> Verbosity level
   integer, intent(in) :: verbosity
   !> Indentation level
   character(len=*), intent(in) :: indent
   !> Information on the container
   character(len=:), allocatable :: str

   integer :: ic
   character(len=*), parameter :: nl = new_line('a')

   if (allocated(self%label)) then
      str = self%label
   else
      str = "Interactions"
   end if

   str = str // format_string(self%nc, '(1x,"(",i0,")")')

   do ic = 1, self%nc
      str = str // nl // indent // format_string(ic, '("[",i0,"]")') // " " // &
         & self%list(ic)%raw%info(verbosity, indent//"* ")
   end do
end function info


!> Add a container
subroutine push_back(self, cont)
   !> Instance of the container list
   class(container_list), intent(inout) :: self
   !> Container to be added
   class(container_type), allocatable, intent(inout) :: cont

   if (.not.allocated(self%list)) call resize(self%list)

   if (allocated(cont)) then
      if (self%nc >= size(self%list)) call resize(self%list)
      self%nc = self%nc + 1
      call move_alloc(cont, self%list(self%nc)%raw)
   end if
end subroutine push_back


!> Add a container
subroutine pop(self, cont, idx)
   !> Instance of the container list
   class(container_list), intent(inout) :: self
   !> Container to be removed
   class(container_type), allocatable, intent(out) :: cont
   !> Index to remove container from
   integer, intent(in), optional :: idx

   integer :: ic, pos

   if (present(idx)) then
      pos = max(1, min(self%nc, idx))
   else
      pos = self%nc
   end if

   if (pos <= 0) return

   call move_alloc(self%list(pos)%raw, cont)

   self%nc = self%nc - 1
   do ic = pos, self%nc
      call move_alloc(self%list(ic+1)%raw, self%list(ic)%raw)
   end do
end subroutine pop


!> Reallocate list of containers
subroutine resize_node(list, n)
   !> Instance of the array to be resized
   type(container_node), allocatable, intent(inout) :: list(:)
   !> Dimension of the final array size
   integer, intent(in), optional :: n

   type(container_node), allocatable :: tmp(:)
   integer, parameter :: initial_size = 20
   integer :: this_size, new_size, item

   if (allocated(list)) then
      this_size = size(list, 1)
      call move_alloc(list, tmp)
   else
      this_size = initial_size
   end if

   if (present(n)) then
      new_size = n
   else
      new_size = this_size + this_size/2 + 1
   end if

   allocate(list(new_size))

   if (allocated(tmp)) then
      this_size = min(size(tmp, 1), size(list, 1))
      do item = 1, this_size
         call move_alloc(tmp(item)%raw, list(item)%raw)
      end do
      deallocate(tmp)
   end if

end subroutine resize_node


!> Inspect container cache and reallocate it in case of type mismatch
subroutine taint(cache, ptr)
   !> Instance of the container cache
   type(container_cache), target, intent(inout) :: cache
   !> Reference to the container cache
   type(cache_list), pointer, intent(out) :: ptr

   if (allocated(cache%raw)) then
      call view(cache, ptr)
      if (associated(ptr)) return
      deallocate(cache%raw)
   end if

   if (.not.allocated(cache%raw)) then
      block
         type(cache_list), allocatable :: tmp
         allocate(tmp)
         call move_alloc(tmp, cache%raw)
      end block
   end if

   call view(cache, ptr)
end subroutine taint

!> Return reference to container cache after resolving its type
subroutine view(cache, ptr)
   !> Instance of the container cache
   type(container_cache), target, intent(inout) :: cache
   !> Reference to the container cache
   type(cache_list), pointer, intent(out) :: ptr
   nullify(ptr)
   select type(target => cache%raw)
   type is(cache_list)
      ptr => target
   end select
end subroutine view


end module tblite_container_list
