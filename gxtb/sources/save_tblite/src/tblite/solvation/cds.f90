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

!> @file tblite/solvation/cds.f90
!> Provides non-polar contributions to the solvation free energy

!> Cavity, dispersion and surface contribution to the solvation free energy
module tblite_solvation_cds
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use mctc_io_math, only : matdet_3x3
   use tblite_blas, only : dot, gemv, symv
   use tblite_container_cache, only : container_cache
   use tblite_mesh_lebedev, only : grid_size, get_angular_grid, list_bisection
   use tblite_scf_info, only : scf_info, atom_resolved, not_used
   use tblite_scf_potential, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_solvation_surface, only : surface_integrator, new_surface_integrator
   use tblite_solvation_data, only : get_vdw_rad_cosmo
   use tblite_solvation_type, only : solvation_type
   use tblite_solvation_cm5, only : get_cm5_charges
   use mctc_io_convert, only : aatoau
   implicit none
   private

   public :: cds_solvation, new_cds, cds_input

   !> Input for CDS model
   type :: cds_input
      !> Probe radius of the solvent
      real(wp) :: probe = 1.00000_wp * aatoau
      !> Number of angular grid points for integration
      integer :: nang = 230
      !> Van-der-Waals radii for each species
      real(wp), allocatable :: rad(:)
      !> Surface tension for each species
      real(wp), allocatable :: tension(:)
      !> Hydrogen bonding strength for each species
      real(wp), allocatable :: hbond(:)
      !> Linearized Poisson-Boltzmann model for parameter selection
      logical :: alpb = .true.
      !> Solvent for parameter selection
      character(len=:), allocatable :: solvent
   end type cds_input

   !> Provide constructor for CDS input
   interface cds_input
      module procedure :: create_cds_input
   end interface cds_input

   !> Definition of the CDS model
   type, extends(solvation_type) :: cds_solvation
      !> Integrator for solvent accessible surface area
      type(surface_integrator) :: sasa
      !> Surface tension for each species
      real(wp), allocatable :: tension(:)
      !> Hydrogen bonding correction
      real(wp), allocatable :: hbond(:)
      !> Use CM5 charges (for GFN1-xTB compatibility)
      logical :: useCM5 = .false.
   contains
      !> Update cache from container
      procedure :: update
      !> Return dependency on density
      procedure :: variable_info
      !> Get solvation energy and gradient
      procedure :: get_engrad
      !> Get solvation energy
      procedure :: get_energy
      !> Get solvation potential
      procedure :: get_potential
      !> Get solvation gradient
      procedure :: get_gradient
   end type cds_solvation

   !> Provide constructor for CDS solvation
   interface cds_solvation
      module procedure :: create_cds
   end interface cds_solvation

   !> Restart data for CDS calculation
   type :: cds_cache
      !> Surface tension
      real(wp), allocatable :: tension(:)
      !> Derivative of surface tension w.r.t. cartesian displacements
      real(wp), allocatable :: dtdr(:, :, :)
      !> Solvent accessible surface area for each atom
      real(wp), allocatable :: surface(:)
      !> Derivative of the solvent accessible surface area w.r.t. cartesian displacements
      real(wp), allocatable :: dsdr(:, :, :)
      !> Scratch array to store potential shift
      real(wp), allocatable :: scratch(:)
      !> Hydrogen bonding correction
      real(wp), allocatable :: hbond(:)
      !> Scratch for atomic charges
      real(wp), allocatable :: qscratch(:)
      !> CM5 charges (only required for GFN1 compatibility)
      real(wp), allocatable :: cm5(:)
      !> CM5 charge derivatives
      real(wp), allocatable :: dcm5dr(:,:,:)
   end type cds_cache

   !> Identifier for container
   character(len=*), parameter :: label = "cds contribution to solvation"

contains


!> Consturctor for CDS input to properly assign allocatable strings
function create_cds_input(solvent, alpb, probe, nang, rad, tension, hbond) result(self)
   !> Solvent for parameter selection
   character(len=*), intent(in), optional :: solvent
   !> Use analytical linearized Poisson-Boltzmann model
   logical, intent(in), optional :: alpb
   !> Probe radius of the solvent
   real(wp), intent(in), optional :: probe
   !> Number of angular grid points for integration
   integer, intent(in), optional :: nang
   !> Van-der-Waals radii for each species
   real(wp), allocatable, intent(in), optional :: rad(:)
   !> Surface tension for each species
   real(wp), allocatable, intent(in), optional :: tension(:)
   !> Hydrogen bonding strength for each species
   real(wp), allocatable, intent(in), optional :: hbond(:)

   type(cds_input) :: self

   if (present(solvent)) then 
      self%solvent = solvent
   end if

   if (present(alpb)) then 
      self%alpb = alpb
   end if

   if (present(probe)) then
      self%probe = probe
   end if

   if (present(nang)) then
      self%nang = nang
   end if

   if (present(rad)) then
      self%rad = rad
   end if

   if (present(tension)) then
      self%tension = tension
   end if

   if (present(hbond)) then
      self%hbond = hbond
   end if

end function create_cds_input


!> Create new CDS solvation model
subroutine new_cds(self, mol, input, method)
   !> Instance of the solvation model
   type(cds_solvation), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Input for CDS solvation
   type(cds_input), intent(in) :: input
   !> Method for parameter selection
   character(len=*), intent(in), optional :: method

   real(wp), allocatable :: rad(:)

   self%label = label
   if (allocated(input%rad)) then
      rad = input%rad
   else
      rad = get_vdw_rad_cosmo(mol%num)
   end if

   self%tension = input%tension
  
   if (allocated(input%hbond)) then
      self%hbond = input%hbond/(4*pi*(rad+input%probe)**2)
   end if

   if (allocated(input%solvent) .and. present(method))then
      self%useCM5 = trim(method)=='gfn1'
   endif

   call new_surface_integrator(self%sasa, mol%id, rad, input%probe, input%nang)
end subroutine new_cds

!> Type constructor for CDS solvation
function create_cds(mol, input, method) result(self)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Input for CDS solvation
   type(cds_input), intent(in) :: input
   !> Method for parameter selection
   character(len=*), intent(in), optional :: method
   !> Instance of the solvation model
   type(cds_solvation) :: self

   call new_cds(self, mol, input, method)
end function create_cds

!> Update cache from container
subroutine update(self, mol, cache)
   !> Instance of the solvation model
   class(cds_solvation), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache

   type(cds_cache), pointer :: ptr

   call taint(cache, ptr)

   if (.not.allocated(ptr%surface)) then
      allocate(ptr%surface(mol%nat))
   end if
   if (.not.allocated(ptr%dsdr)) then
      allocate(ptr%dsdr(3, mol%nat, mol%nat))
   end if
   if (.not.allocated(ptr%scratch)) then
      allocate(ptr%scratch(mol%nat))
   end if
   if (.not.allocated(ptr%qscratch))then
      allocate(ptr%qscratch(mol%nat))
   endif
   if (self%useCM5)then
      if (.not.allocated(ptr%cm5))then
         allocate(ptr%cm5(mol%nat))
      endif
      if (.not.allocated(ptr%dcm5dr))then
         allocate(ptr%dcm5dr(3, mol%nat, mol%nat))
      endif
      call get_cm5_charges(mol, ptr%cm5, ptr%dcm5dr)
   endif

   call self%sasa%get_surface(mol, ptr%surface, ptr%dsdr)

   ptr%tension = self%tension(mol%id)
   if (allocated(self%hbond)) then
      ptr%hbond = self%hbond(mol%id)
   end if
end subroutine update


!> Evaluate non-selfconsistent part of the interaction
subroutine get_engrad(self, mol, cache, energies, gradient, sigma)
   !> Instance of the solvation model
   class(cds_solvation), intent(in) :: self
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

   type(cds_cache), pointer :: ptr

   call view(cache, ptr)

   energies(:) = energies + ptr%surface * ptr%tension

   if (present(gradient)) then
      if (allocated(ptr%dsdr)) &
         call gemv(ptr%dsdr, ptr%tension, gradient, beta=1.0_wp)
      if (allocated(ptr%dtdr)) &
         call gemv(ptr%dtdr, ptr%surface, gradient, beta=1.0_wp)
   end if
end subroutine get_engrad


!> Get solvation energy
subroutine get_energy(self, mol, cache, wfn, energies)
   !> Instance of the solvation model
   class(cds_solvation), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Solvation free energy
   real(wp), intent(inout) :: energies(:)

   type(cds_cache), pointer :: ptr

   call view(cache, ptr)

   if (allocated(self%hbond)) then
      if(self%useCM5)then
         ptr%qscratch(:) = wfn%qat(:, 1) + ptr%cm5(:)
      else
         ptr%qscratch(:) = wfn%qat(:, 1)
      endif
      ptr%scratch(:) = ptr%hbond * ptr%surface * ptr%qscratch(:)**2
      energies(:) = energies + ptr%scratch
   end if
end subroutine get_energy


!> Get solvation potential
subroutine get_potential(self, mol, cache, wfn, pot)
   !> Instance of the solvation model
   class(cds_solvation), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential
   type(potential_type), intent(inout) :: pot

   type(cds_cache), pointer :: ptr

   call view(cache, ptr)

   if (allocated(self%hbond)) then
      if(self%useCM5)then
         ptr%qscratch(:) = wfn%qat(:, 1) + ptr%cm5(:)
      else
         ptr%qscratch(:) = wfn%qat(:, 1)
      endif
      ptr%scratch(:) = 2*ptr%hbond * ptr%surface * ptr%qscratch(:)
      pot%vat(:, 1) = pot%vat(:, 1) + ptr%scratch(:)
   end if
end subroutine get_potential


!> Get solvation gradient
subroutine get_gradient(self, mol, cache, wfn, gradient, sigma)
   !> Instance of the solvation model
   class(cds_solvation), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Molecular gradient of the solvation free energy
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Strain derivatives of the solvation free energy
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   type(cds_cache), pointer :: ptr

   call view(cache, ptr)

   if (allocated(self%hbond)) then
      if(self%useCM5)then
         ptr%qscratch(:) = wfn%qat(:, 1) + ptr%cm5(:)
      else
         ptr%qscratch(:) = wfn%qat(:, 1)
      endif
      ptr%scratch(:) = ptr%hbond * ptr%qscratch(:)**2
      if (allocated(ptr%dsdr)) &
         call gemv(ptr%dsdr, ptr%scratch, gradient, beta=1.0_wp)
      if(self%useCM5)then 
        ptr%scratch(:) = 2*ptr%hbond * ptr%surface * ptr%qscratch(:)  
        call gemv(ptr%dcm5dr, ptr%scratch, gradient, beta=1.0_wp) 
      endif  
   end if
end subroutine get_gradient


!> Return dependency on density
pure function variable_info(self) result(info)
   !> Instance of the solvation model
   class(cds_solvation), intent(in) :: self
   !> Information on the required potential data
   type(scf_info) :: info

   info = scf_info(charge=merge(atom_resolved, not_used, allocated(self%hbond)))
end function variable_info


subroutine taint(cache, ptr)
   type(container_cache), target, intent(inout) :: cache
   type(cds_cache), pointer, intent(out) :: ptr

   if (allocated(cache%raw)) then
      call view(cache, ptr)
      if (associated(ptr)) return
      deallocate(cache%raw)
   end if

   if (.not.allocated(cache%raw)) then
      block
         type(cds_cache), allocatable :: tmp
         allocate(tmp)
         call move_alloc(tmp, cache%raw)
      end block
   end if

   call view(cache, ptr)
end subroutine taint

subroutine view(cache, ptr)
   type(container_cache), target, intent(inout) :: cache
   type(cds_cache), pointer, intent(out) :: ptr
   nullify(ptr)
   select type(target => cache%raw)
   type is(cds_cache)
      ptr => target
   end select
end subroutine view

end module tblite_solvation_cds
