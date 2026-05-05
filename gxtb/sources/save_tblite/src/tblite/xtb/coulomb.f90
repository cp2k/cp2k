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

!> @file tblite/xtb/coulomb.f90
!> Provides a collection for all Coulomb related interactions

!> Collects all Coulomb related interaction in as single interaction container.
module tblite_xtb_coulomb
   use mctc_data_covrad, only : get_covalent_rad
   use mctc_env, only : error_type, wp
   use mctc_io, only : structure_type
   use mctc_ncoord, only : new_ncoord, ncoord_type
   use tblite_container, only : container_cache, container_type
   use tblite_coulomb, only : onsite_firstorder, coulomb_charge_type, &
      & damped_multipole, thirdorder_type, onsite_fourthorder, coulomb_cache
   use tblite_scf_potential, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   public :: new_coulomb

   !> Collection of Coulombic/tight-binding interactions
   type, public, extends(container_type) :: tb_coulomb
      !> Onsite first-order tight-binding
      type(onsite_firstorder), allocatable :: es1
      !> Isotroptic second-order electrostatics
      class(coulomb_charge_type), allocatable :: es2
      !> Anisotropic second-order electrostatic
      class(damped_multipole), allocatable :: aes2
      !> Isotropic third-order tight-binding
      class(thirdorder_type), allocatable :: es3
      !> Onsite fourth-order tight-binding
      type(onsite_fourthorder), allocatable :: es4
      !> Coordination number for Coulomb/tight-binding
      class(ncoord_type), allocatable :: ncoord
      !> Valence coordination number
      real(wp), allocatable :: valence_cn(:)
   contains
      procedure :: update
      procedure :: variable_info
      procedure :: get_energy
      procedure :: get_potential
      procedure :: get_potential_gradient
      procedure :: get_gradient
      procedure :: info
   end type tb_coulomb

   character(len=*), parameter :: label = "Coulomb electrostatics and tight-binding"

contains


!> Create new collection of coulombic/tight-binding interaction
subroutine new_coulomb(self, mol, error, cn_count_type, cn_rcov, cn_exp, valence_cn)
   !> Instance of coulomb and tight-binding collection
   type(tb_coulomb), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Optional type of coordination number counting function
   integer, intent(in), optional :: cn_count_type
   !> Optional covalent atomic radii for Coulomb/tight-binding coordination number
   real(wp), intent(in), optional :: cn_rcov(:)
   !> Optional exponent of the Coulomb/tight-binding coordination number
   real(wp), intent(in), optional :: cn_exp
   !> Valence coordination number
   real(wp), intent(in), optional :: valence_cn(:)

   real(wp), allocatable :: rad(:)

   self%label = label

   ! Obtain covalent radii for Coulomb/tight-binding coordination number
   allocate(rad(mol%nid))
   if (present(cn_rcov)) then
      rad(:) = cn_rcov
   else
      rad(:) = get_covalent_rad(mol%num)
   end if

   ! Create the coordination number for Coulomb/tight-binding interactions
   ! if required by any of the interactions
   if (present(cn_count_type)) then 
      call new_ncoord(self%ncoord, mol, cn_count_type=cn_count_type, &
         & error=error, kcn=cn_exp, rcov=rad)
   end if

   ! Check for alternative valence coordination number
   allocate(self%valence_cn(mol%nid))
   if (present(valence_cn)) then
      self%valence_cn = valence_cn
   else
      self%valence_cn(:) = 0.0_wp
   end if

end subroutine new_coulomb


subroutine update(self, mol, cache)
   !> Instance of the tight-binding and tight-binding container
   class(tb_coulomb), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache

   type(coulomb_cache), pointer :: ptr

   call taint(cache, ptr)

   ! Obtain optional coordination number for all Coulomb/tight-binding interactions
   if (.not.allocated(ptr%cn)) allocate(ptr%cn(mol%nat))
   if (.not.allocated(ptr%dcndr)) allocate(ptr%dcndr(3, mol%nat, mol%nat))
   if (.not.allocated(ptr%dcndL)) allocate(ptr%dcndL(3, 3, mol%nat))
   if (allocated(self%ncoord)) then
      call self%ncoord%get_cn(mol, ptr%cn, ptr%dcndr, ptr%dcndL)
   else 
      ! In case of no coordination number, we use the valence CN as default
      ptr%cn(:) = self%valence_cn(mol%id)
      ptr%dcndr(:, :, :) = 0.0_wp
      ptr%dcndL(:, :, :) = 0.0_wp
   end if

   if (allocated(self%es1)) then
      call self%es1%update(mol, cache)
   end if

   if (allocated(self%es2)) then
      call self%es2%update(mol, cache)
   end if

   if (allocated(self%aes2)) then
      call self%aes2%update(mol, cache)
   end if

   if (allocated(self%es3)) then
      call self%es3%update(mol, cache)
   end if

   if (allocated(self%es4)) then
      call self%es4%update(mol, cache)
   end if
end subroutine update


subroutine get_energy(self, mol, cache, wfn, energies)
   !> Instance of the electrostatic and tight-binding container
   class(tb_coulomb), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Electrostatic and tight-binding energy
   real(wp), intent(inout) :: energies(:)

   if (allocated(self%es1)) then
      call self%es1%get_energy(mol, cache, wfn, energies)
   end if
   
   if (allocated(self%es2)) then
      call self%es2%get_energy(mol, cache, wfn, energies)
   end if

   if (allocated(self%aes2)) then
      call self%aes2%get_energy(mol, cache, wfn, energies)
   end if

   if (allocated(self%es3)) then
      call self%es3%get_energy(mol, cache, wfn, energies)
   end if

   if (allocated(self%es4)) then
      call self%es4%get_energy(mol, cache, wfn, energies)
   end if
end subroutine get_energy


subroutine get_potential(self, mol, cache, wfn, pot)
   !> Instance of the electrostatic and tight-binding container
   class(tb_coulomb), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential
   type(potential_type), intent(inout) :: pot

   if (allocated(self%es1)) then
      call self%es1%get_potential(mol, cache, wfn, pot)
   end if

   if (allocated(self%es2)) then
      call self%es2%get_potential(mol, cache, wfn, pot)
   end if

   if (allocated(self%aes2)) then
      call self%aes2%get_potential(mol, cache, wfn, pot)
   end if

   if (allocated(self%es3)) then
      call self%es3%get_potential(mol, cache, wfn, pot)
   end if

   if (allocated(self%es4)) then
      call self%es4%get_potential(mol, cache, wfn, pot)
   end if
end subroutine get_potential


subroutine get_potential_gradient(self, mol, cache, wfn, pot)
   !> Instance of the electrostatic and tight-binding container
   class(tb_coulomb), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential
   type(potential_type), intent(inout) :: pot

   if (allocated(self%es1)) then
      call self%es1%get_potential_gradient(mol, cache, wfn, pot)
   end if

   if (allocated(self%es2)) then
      call self%es2%get_potential_gradient(mol, cache, wfn, pot)
   end if

   if (allocated(self%aes2)) then
      call self%aes2%get_potential_gradient(mol, cache, wfn, pot)
   end if

   if (allocated(self%es3)) then
      call self%es3%get_potential_gradient(mol, cache, wfn, pot)
   end if

   if (allocated(self%es4)) then
      call self%es4%get_potential_gradient(mol, cache, wfn, pot)
   end if
end subroutine get_potential_gradient


subroutine get_gradient(self, mol, cache, wfn, gradient, sigma)
   !> Instance of the electrostatic and tight-binding container
   class(tb_coulomb), intent(in) :: self
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

   if (allocated(self%es1)) then
      call self%es1%get_gradient(mol, cache, wfn, gradient, sigma)
   end if

   if (allocated(self%es2)) then
      call self%es2%get_gradient(mol, cache, wfn, gradient, sigma)
   end if

   if (allocated(self%aes2)) then
      call self%aes2%get_gradient(mol, cache, wfn, gradient, sigma)
   end if

   if (allocated(self%es3)) then
      call self%es3%get_gradient(mol, cache, wfn, gradient, sigma)
   end if

   if (allocated(self%es4)) then
      call self%es4%get_gradient(mol, cache, wfn, gradient, sigma)
   end if
end subroutine get_gradient


pure function variable_info(self) result(info)
   use tblite_scf_info, only : scf_info, max
   !> Instance of the electrostatic and tight-binding container
   class(tb_coulomb), intent(in) :: self
   !> Information on the required potential data
   type(scf_info) :: info

   info = scf_info()

   if (allocated(self%es1)) then
      info = max(info, self%es1%variable_info())
   end if

   if (allocated(self%es2)) then
      info = max(info, self%es2%variable_info())
   end if

   if (allocated(self%aes2)) then
      info = max(info, self%aes2%variable_info())
   end if

   if (allocated(self%es3)) then
      info = max(info, self%es3%variable_info())
   end if

   if (allocated(self%es4)) then
      info = max(info, self%es4%variable_info())
   end if
end function variable_info


!> Information on container
pure function info(self, verbosity, indent) result(str)
   !> Instance of the electrostatic and tight-binding container
   class(tb_coulomb), intent(in) :: self
   !> Verbosity level
   integer, intent(in) :: verbosity
   !> Indentation level
   character(len=*), intent(in) :: indent
   !> Information on the container
   character(len=:), allocatable :: str

   character(len=*), parameter :: nl = new_line('a'), marker = " * "

   if (allocated(self%label)) then
      str = self%label
   else
      str = "Coulomb electrostatics and tight-binding"
   end if

   if (allocated(self%es1)) then
      str = str // nl // indent // marker // &
         & self%es1%info(verbosity, indent//marker)
   end if

   if (allocated(self%es2)) then
      str = str // nl // indent // marker // &
         & self%es2%info(verbosity, indent//marker)
   end if

   if (allocated(self%aes2)) then
      str = str // nl // indent // marker // &
         & self%aes2%info(verbosity, indent//marker)
   end if

   if (allocated(self%es3)) then
      str = str // nl // indent // marker // &
         & self%es3%info(verbosity, indent//marker)
   end if

   if (allocated(self%es4)) then
      str = str // nl // indent // marker // &
         & self%es4%info(verbosity, indent//marker)
   end if
end function info


!> Inspect container cache and reallocate it in case of type mismatch
subroutine taint(cache, ptr)
   !> Instance of the container cache
   type(container_cache), target, intent(inout) :: cache
   !> Reference to the coulomb cache
   type(coulomb_cache), pointer, intent(out) :: ptr

   if (allocated(cache%raw)) then
      call view(cache, ptr)
      if (associated(ptr)) return
      deallocate(cache%raw)
   end if

   if (.not.allocated(cache%raw)) then
      block
         type(coulomb_cache), allocatable :: tmp
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
   type(coulomb_cache), pointer, intent(out) :: ptr
   nullify(ptr)
   select type(target => cache%raw)
   type is(coulomb_cache)
      ptr => target
   end select
end subroutine view

end module tblite_xtb_coulomb
