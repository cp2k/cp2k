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

!> @file tblite/scf/iterator.f90
!> Provides the implementation of the actual self-consistent field iteractions

!> Iterator for evaluating the Hamiltonian self-consistently
module tblite_scf_iterator
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type
   use tblite_basis_type, only : basis_type
   use tblite_container, only : container_cache, container_list
   use tblite_disp, only : dispersion_type
   use tblite_exchange, only : exchange_type
   use tblite_integral_type, only : integral_type
   use tblite_repulsion_type, only : repulsion_type
   use tblite_scf_cache, only : iterator_cache
   use tblite_scf_info, only : scf_info, atom_resolved, shell_resolved, &
      & orbital_resolved
   use tblite_scf_mixer, only : new_mixer
   use tblite_scf_mixer_broyden, only : broyden_mixer
   use tblite_scf_mixer_diis, only : diis_mixer
   use tblite_scf_mixer_input, only : mixer_input_container, mixer_mode
   use tblite_scf_mixer_simple, only : simple_mixer
   use tblite_scf_mixer_type, only : mixer_container
   use tblite_scf_potential, only : potential_type, add_pot_to_h1
   use tblite_scf_solver, only : solver_type
   use tblite_spin, only : spin_polarization
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_wavefunction_fermi, only : get_fermi_filling
   use tblite_wavefunction_mulliken, only : get_mulliken_shell_charges, &
      & get_mulliken_atomic_multipoles
   use tblite_xtb_coulomb, only : tb_coulomb
   implicit none
   private

   public :: next_scf, get_electronic_energy, reduce
   public :: next_density, get_qat_from_qsh, new_iterator

   !> Default maximum number of self-consistent iterations
   integer, parameter :: max_iter_default = 250

   !> Iterator data type controlling the SCF iterations
   type, public :: iterator_type
      !> Energy convergence threshold
      real(wp) :: econv
      !> Density convergence threshold
      real(wp) :: pconv
      !> Maximum number of self-consistent iteractions
      integer :: max_iter
      !> List of mixers to be used during the SCF
      type(mixer_container), allocatable :: mixers(:)
   contains
      !> Update cache
      procedure :: update
      !> Peform the next SCF iteration
      procedure :: next_scf
      !> Check convergence
      procedure :: check_convergence
      !> Check maximum number of iterations
      procedure :: check_max_iter
      !> Modify the starting iteration for the mixers
      procedure :: set_mixer_start
      !> Modify the damping for the mixers
      procedure :: set_mixer_damping
      !> Set density/potential in the mixer
      procedure :: set_mixer
      !> Set difference of density/potential in the mixer
      procedure :: diff_mixer
      !> Apply the mixer for the next density/potential matrix
      procedure :: next_mixer
      !> Get density/potential matrix from the mixer
      procedure :: get_mixer
      !> Collect history for the mixer
      procedure :: collect_mixer
      !> Cleanup mixer data
      procedure :: cleanup
   end type iterator_type


contains


subroutine new_iterator(self, mol, bas, mixer_inputs, max_iter, accuracy)
   !> Instance of the iterator
   type(iterator_type), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Optional list of input parameters for all used electronic mixers
   type(mixer_input_container), intent(in), optional :: mixer_inputs(:)
   !> Optional maximum number of SCF iterations
   integer, intent(in), optional :: max_iter
   !> Optional accuracy specification
   real(wp), intent(in), optional :: accuracy

   integer :: imix
   real(wp) :: acc

   ! Only add a mixer if inputs are provided
   if (present(mixer_inputs)) then
      allocate(self%mixers(size(mixer_inputs)))
      do imix = 1, size(mixer_inputs)
         call new_mixer(self%mixers(imix)%raw, mixer_inputs(imix)%raw)
      end do
   end if

   if (present(max_iter)) then
      self%max_iter = max_iter
   else
      self%max_iter = max_iter_default
   end if

   if (present(accuracy)) then
      acc = accuracy
   else
      acc = 1.0_wp
   end if

   self%econv = 1.e-6_wp * acc
   self%pconv = 2.e-5_wp * acc

end subroutine new_iterator


subroutine update(self, cache, mol, bas, wfn, info, energy)
   !> Instance of the iterator
   class(iterator_type), intent(in) :: self
   !> Cache for the iterator
   type(iterator_cache), intent(inout) :: cache
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Tight-binding wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Information about the used SCF contributions
   type(scf_info), intent(in) :: info
   !> Current total energy
   real(wp), intent(in) :: energy

   integer :: imix, ndim_density, ndim_potential, ndim
   real(wp) :: de
   logical :: new_current_mixer, new_precollect_mixer

   ! Increment iteration count
   cache%iscf = cache%iscf + 1

   ! Save the SCF information
   cache%info = info

   ! Update mixer and mixer cache
   if (allocated(self%mixers)) then

      ! Allocate list of caches to be used in the scf   
      if (.not. allocated(cache%mcache)) allocate(cache%mcache(size(self%mixers)))

      new_current_mixer = .false.
      new_precollect_mixer = .false.
      de = huge(1.0_wp)
      if (cache%iscf > 1) then
         de = abs(energy - cache%elast)
         ! Check if we should switch to the next mixer or even further
         do imix = cache%current_mixer + 1, size(self%mixers)
            associate(mixer => self%mixers(imix)%raw)
               if (cache%iscf >= mixer%start .or. de < mixer%conv_start) then
                  ! Remove old mixer data
                  call self%mixers(cache%current_mixer)%raw%cleanup(cache%mcache(cache%current_mixer))
                  cache%current_mixer = imix
                  new_current_mixer = .true.
               end if
            end associate
         end do
      else
         cache%current_mixer = 1
         new_current_mixer = .true.
      end if

      ! Check if the next mixer should already precollect
      cache%precollect_mixer = -1
      do imix = cache%current_mixer + 1, size(self%mixers)
         associate(mixer => self%mixers(imix)%raw)
            ! Setup precollect one iteration before to allow a seeding iteration
            if (cache%iscf >= mixer%start - mixer%precollect - 1 &
               & .or. de < mixer%conv_start) then
               cache%precollect_mixer = imix
               new_precollect_mixer = .true.
               exit
            end if
         end associate
      end do
      
      if (new_current_mixer .or. new_precollect_mixer) then
         ! Collect the information about the mixer dimensions
         call get_mixer_dimensions(mol, bas, wfn, cache%info, ndim_density, ndim_potential)
         if (new_current_mixer) then
            ! Select actual mixer dimension based on the mode of the current mixer
            ndim = merge(ndim_density, ndim_potential, &
               & self%mixers(cache%current_mixer)%raw%mode == mixer_mode%density)
            ! Update and initialize the current mixer cache
            call self%mixers(cache%current_mixer)%raw%update(cache%mcache(cache%current_mixer), ndim)
         end if
         if (new_precollect_mixer) then
            ! Select actual mixer dimension based on the mode of the precollecting mixer
            ndim = merge(ndim_density, ndim_potential, &
               & self%mixers(cache%precollect_mixer)%raw%mode == mixer_mode%density)
            ! Update and initialize the precollecting mixer cache
            call self%mixers(cache%precollect_mixer)%raw%update(cache%mcache(cache%precollect_mixer), ndim)
         end if
      end if
   end if

   ! Save current energy
   cache%elast = energy

end subroutine update


subroutine get_mixer_dimensions(mol, bas, wfn, info, ndim_density, ndim_potential)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Tight-binding wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Information about the used SCF contributions
   type(scf_info), intent(in) :: info
   !> Size of the extrapolated quantity for density mixing
   integer, intent(out) :: ndim_density
   !> Size of the extrapolated quantity for potential mixing
   integer, intent(out) :: ndim_potential

   ndim_density = 0
   ndim_potential = 0

   select case(info%charge)
   case(atom_resolved)
      ndim_density = ndim_density + mol%nat
   case(shell_resolved)
      ndim_density = ndim_density + bas%nsh
   end select

   select case(info%dipole)
   case(atom_resolved)
      ndim_density = ndim_density + 3*mol%nat
   end select

   select case(info%quadrupole)
   case(atom_resolved)
      ndim_density = ndim_density + 6*mol%nat
   end select

   select case(info%density)
   case(orbital_resolved)
      ndim_density = ndim_density + bas%nao * bas%nao
   end select

   ndim_density = ndim_density * wfn%nspin

   ! For the potential so far only based on the Fock matrix
   ndim_potential = ndim_potential + bas%nao * bas%nao
   ndim_potential = ndim_potential * wfn%nspin

end subroutine get_mixer_dimensions


!> Evaluate self-consistent iteration for the density-dependent Hamiltonian
subroutine next_scf(self, mol, bas, wfn, solver, iter_cache, coulomb, &
   & dispersion, exchange, repulsion, spinpol, interactions, ints, pot, &
   & ccache, dcache, ecache, rcache, scache, icache, energies, error)
   !> Instance of the iterator
   class(iterator_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Tight-binding wavefunction data
   type(wavefunction_type), intent(inout) :: wfn
   !> Solver for the general eigenvalue problem
   class(solver_type), intent(inout) :: solver
   !> Cache for the iterator
   type(iterator_cache), intent(inout) :: iter_cache
   !> Container for coulombic and tight-binding interactions
   type(tb_coulomb), intent(in), optional :: coulomb
   !> Container for dispersion interactions
   class(dispersion_type), intent(in), optional :: dispersion
   !> Container for exchange interactions
   class(exchange_type), intent(in), optional :: exchange
   !> Container for repulsion interactions
   class(repulsion_type), intent(in), optional :: repulsion
   !> Container for spin-polarization interactions
   type(spin_polarization), intent(in), optional :: spinpol
   !> Container for general interactions
   type(container_list), intent(in), optional :: interactions
   !> Integral container
   type(integral_type), intent(in) :: ints
   !> Density dependent potential shifts
   type(potential_type), intent(inout) :: pot
   !> Restart data for coulombic and tight-binding interactions
   type(container_cache), intent(inout), optional :: ccache
   !> Restart data for dispersion interactions
   type(container_cache), intent(inout), optional :: dcache
   !> Restart data for exchange containers
   type(container_cache), intent(inout), optional :: ecache
   !> Restart data for repulsion containers
   type(container_cache), intent(inout), optional :: rcache
   !> Restart data for spin-polarization containers
   type(container_cache), intent(inout), optional :: scache
   !> Restart data for interaction containers
   type(container_cache), intent(inout), optional :: icache

   !> Self-consistent energy
   real(wp), intent(inout) :: energies(:)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), allocatable :: eao(:)
   real(wp) :: ts

   if (iter_cache%iscf > 1) then
      ! Mix density-related quantities
      call self%next_mixer(iter_cache, mixer_mode%density, error)
      if (allocated(error)) return
      ! Retrieve mixed density-related quantities
      call self%get_mixer(iter_cache, bas, wfn, pot, mixer_mode%density)
   end if

   ! Construct the density-dependent potentials and Fock matrix contributions
   call pot%reset
   if (present(coulomb) .and. present(ccache)) then
      call coulomb%get_potential(mol, ccache, wfn, pot)
   end if
   if (present(dispersion) .and. present(dcache)) then
      call dispersion%get_potential(mol, dcache, wfn, pot)
   end if
   if (present(exchange) .and. present(ecache)) then
      call exchange%get_potential_w_overlap(mol, ecache, wfn, ints%overlap, pot)
   end if
   if (present(repulsion) .and. present(rcache)) then
      call repulsion%get_potential(mol, rcache, wfn, pot)
   end if
   if (present(spinpol) .and. present(scache)) then
      call spinpol%get_potential(mol, scache, wfn, pot)
   end if
   if (present(interactions) .and. present(icache)) then
      call interactions%get_potential(mol, icache, wfn, pot)
   end if
   ! Add the potential contributions to the Fock matrix
   call add_pot_to_h1(bas, ints, pot, wfn%coeff)

   if (iter_cache%iscf > 1) then
      ! Calculate potential difference relativ to last iteration
      call self%diff_mixer(iter_cache, wfn, pot, mixer_mode%potential)
      ! Collect potential-related quantities in the mixer cache
      call self%collect_mixer(iter_cache, mixer_mode%potential)
      ! Mix potential-related quantities
      call self%next_mixer(iter_cache, mixer_mode%potential, error)
      if (allocated(error)) return
      ! Retrieve mixed potential-related quantities
      call self%get_mixer(iter_cache, bas, wfn, pot, mixer_mode%potential)
   end if

   ! Save current extrapolated density-related quantities before diagonalization
   call self%set_mixer(iter_cache, wfn, pot, mixer_mode%density)
   ! Save current extrapolated potential-related quantities before diagonalization
   call self%set_mixer(iter_cache, wfn, pot, mixer_mode%potential)

   ! Solve the Roothaan-Hall equations and construct a new density matrix
   call next_density(wfn, solver, ints, ts, error)
   if (allocated(error)) return

   ! Partition the density matrix into Mulliken charges and multipoles
   call get_mulliken_shell_charges(bas, ints%overlap, wfn%density, wfn%n0sh, &
      & wfn%qsh)
   call get_qat_from_qsh(bas, wfn%qsh, wfn%qat)

   call get_mulliken_atomic_multipoles(bas, ints%dipole, wfn%density, &
      & wfn%dpat)
   call get_mulliken_atomic_multipoles(bas, ints%quadrupole, wfn%density, &
      & wfn%qpat)

   ! Calculate the difference for density-related quantities after diagonalization
   call self%diff_mixer(iter_cache, wfn, pot, mixer_mode%density)
   ! Collect density-related quantities in the mixer cache
   call self%collect_mixer(iter_cache, mixer_mode%density)

   ! Compute the electronic energy as Tr(P x F)
   allocate(eao(bas%nao), source=0.0_wp)
   call get_electronic_energy(ints%hamiltonian, wfn%density, eao)
   energies(:) = ts / size(energies)
   call reduce(energies, eao, bas%ao2at)

   ! Add energy from density-dependent interactions
   if (present(coulomb) .and. present(ccache)) then
      call coulomb%get_energy(mol, ccache, wfn, energies)
   end if
   if (present(dispersion) .and. present(dcache)) then
      call dispersion%get_energy(mol, dcache, wfn, energies)
   end if
   if (present(exchange) .and. present(ecache)) then
      call exchange%get_energy_w_overlap(mol, ecache, wfn, ints%overlap, energies)
   end if
   if (present(repulsion) .and. present(rcache)) then
      call repulsion%get_energy(mol, rcache, wfn, energies)
   end if
   if (present(spinpol) .and. present(scache)) then
      call spinpol%get_energy(mol, scache, wfn, energies)
   end if
   if (present(interactions) .and. present(icache)) then
      call interactions%get_energy(mol, icache, wfn, energies)
   end if

end subroutine next_scf


subroutine get_electronic_energy(h0, density, energies)
   real(wp), intent(in) :: h0(:, :)
   real(wp), intent(in) :: density(:, :, :)
   real(wp), intent(inout) :: energies(:)

   integer :: iao, jao, spin

   !$omp parallel do collapse(3) schedule(runtime) default(none) &
   !$omp reduction(+:energies) shared(h0, density) private(spin, iao, jao)
   do spin = 1, size(density, 3)
      do iao = 1, size(density, 2)
         do jao = 1, size(density, 1)
            energies(iao) = energies(iao) + h0(jao, iao) * density(jao, iao, spin)
         end do
      end do
   end do
end subroutine get_electronic_energy


subroutine reduce(reduced, full, map)
   real(wp), intent(inout) :: reduced(:)
   real(wp), intent(in) :: full(:)
   integer, intent(in) :: map(:)

   integer :: ix

   do ix = 1, size(map)
      reduced(map(ix)) = reduced(map(ix)) + full(ix)
   end do
end subroutine reduce


subroutine get_qat_from_qsh(bas, qsh, qat)
   type(basis_type), intent(in) :: bas
   real(wp), intent(in) :: qsh(:, :)
   real(wp), intent(out) :: qat(:, :)

   integer :: ish, ispin

   qat(:, :) = 0.0_wp
   !$omp parallel do schedule(runtime) collapse(2) default(none) &
   !$omp reduction(+:qat) shared(bas, qsh) private(ish)
   do ispin = 1, size(qsh, 2)
      do ish = 1, size(qsh, 1)
         qat(bas%sh2at(ish), ispin) = qat(bas%sh2at(ish), ispin) + qsh(ish, ispin)
      end do
   end do
end subroutine get_qat_from_qsh


subroutine next_density(wfn, solver, ints, ts, error)
   !> Tight-binding wavefunction data
   type(wavefunction_type), intent(inout) :: wfn
   !> Solver for the general eigenvalue problem
   class(solver_type), intent(inout) :: solver
   !> Integral container
   type(integral_type), intent(in) :: ints
   !> Electronic entropy
   real(wp), intent(out) :: ts
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp) :: e_fermi, stmp(2)
   real(wp), allocatable :: focc(:)
   integer :: spin

   call solver%get_density(wfn%coeff, ints%overlap, wfn%emo, wfn%focc, &
      & wfn%density, error)
   do spin = 1, 2
      call get_electronic_entropy(wfn%focc(:, spin), wfn%kt, stmp(spin))
   end do
   ts = sum(stmp)
end subroutine next_density


subroutine get_electronic_entropy(occ, kt, s)
   real(wp), intent(in) :: occ(:)
   real(wp), intent(in) :: kt
   real(wp), intent(out) :: s

   s = sum(log(occ ** occ * (1 - occ) ** (1 - occ))) * kt
end subroutine get_electronic_entropy


subroutine check_convergence(self, cache, energy, energy_error, density_error, &
   & econverged, pconverged, converged)
   !> Instance of the iterator
   class(iterator_type), intent(in) :: self
   !> Cache for the electronic mixer
   type(iterator_cache), intent(inout) :: cache
   !> Current total energy
   real(wp), intent(in) :: energy
   !> Energy error
   real(wp), intent(out) :: energy_error
   !> Density error
   real(wp), intent(out) :: density_error
   !> Energy convergence flag
   logical, intent(out) :: econverged
   !> Density convergence flag
   logical, intent(out) :: pconverged
   !> Convergence flag
   logical, intent(out) :: converged

   energy_error = energy - cache%elast
   if (allocated(self%mixers)) then
      density_error = self%mixers(cache%current_mixer)%raw%get_error( &
         & cache%mcache(cache%current_mixer))
   else
      density_error = 0.0_wp
   end if

   if (abs(energy_error) < self%econv) then 
      econverged = .true. 
   end if
   if (density_error < self%pconv .and. cache%iscf > 1) then
      pconverged = .true.
   end if

   converged = econverged .and. pconverged

end subroutine check_convergence


pure function check_max_iter(self, cache) result(below)
   !> Instance of the iterator
   class(iterator_type), intent(in) :: self
   !> Cache for the electronic mixer
   type(iterator_cache), intent(in) :: cache
   !> Result flag indicating if maximum iterations is not exceeded
   logical :: below

   below = cache%iscf < self%max_iter

end function check_max_iter


subroutine set_mixer_start(self, broyden_start, diis_start)
   !> Instance of the iterator
   class(iterator_type), intent(inout) :: self
   !> Starting iterator for Broyden mixer
   integer, intent(in), optional :: broyden_start
   !> Starting iterator for DIIS mixer
   integer, intent(in), optional :: diis_start

   integer :: imix

   if (.not. allocated(self%mixers)) return

   do imix = 1, size(self%mixers)
      select type(mixer => self%mixers(imix)%raw)
      type is (broyden_mixer)
         if (present(broyden_start)) then
            mixer%start = broyden_start
         end if
      type is (diis_mixer)
         if (present(diis_start)) then
            mixer%start = diis_start
         end if
      end select
   end do 

end subroutine set_mixer_start


subroutine set_mixer_damping(self, simple_damping, broyden_damping)
   !> Instance of the iterator
   class(iterator_type), intent(inout) :: self
   !> Damping for simple mixer
   real(wp), intent(in), optional :: simple_damping
   !> Damping for Broyden mixer
   real(wp), intent(in), optional :: broyden_damping

   integer :: imix

   if (.not. allocated(self%mixers)) return

   do imix = 1, size(self%mixers)
      select type(mixer => self%mixers(imix)%raw)
      type is (simple_mixer)
         if (present(simple_damping)) then
            mixer%damp = simple_damping
         end if
      type is (broyden_mixer)
         if (present(broyden_damping)) then
            mixer%damp = broyden_damping
         end if
      end select
   end do 

end subroutine set_mixer_damping


subroutine set_mixer(self, cache, wfn, pot, mode)
   !> Instance of the iterator
   class(iterator_type), intent(in) :: self
   !> Cache for the iterator
   type(iterator_cache), intent(inout) :: cache
   !> Tight-binding wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential shifts
   type(potential_type), intent(in) :: pot
   !> Mode of the mixer
   integer, intent(in) :: mode

   integer :: mixer_index, targets(2), it

   if (.not. allocated(self%mixers)) return

   targets = [cache%current_mixer, cache%precollect_mixer]

   do it = 1, 2
      mixer_index = targets(it)
      
      ! Skip invalid mixer indices
      if (mixer_index <= 0) cycle
      if (mixer_index > size(self%mixers)) cycle
      if (.not. allocated(cache%mcache(mixer_index)%raw)) cycle

      ! Check if the mixer mode matches
      if (self%mixers(mixer_index)%raw%mode /= mode) cycle

      ! Reset chunk index
      cache%mcache(mixer_index)%raw%iset = 0
      
      associate (mixer => self%mixers(mixer_index)%raw)
         if (mode == mixer_mode%density) then
            select case(cache%info%charge)
            case(atom_resolved)
               call mixer%set(cache%mcache(mixer_index), wfn%qat)
            case(shell_resolved)
               call mixer%set(cache%mcache(mixer_index), wfn%qsh)
            end select

            select case(cache%info%dipole)
            case(atom_resolved)
               call mixer%set(cache%mcache(mixer_index), wfn%dpat)
            end select

            select case(cache%info%quadrupole)
            case(atom_resolved)
               call mixer%set(cache%mcache(mixer_index), wfn%qpat)
            end select

            select case(cache%info%density)
            case(orbital_resolved)
               call mixer%set(cache%mcache(mixer_index), wfn%density)
            end select
         else if (mode == mixer_mode%potential) then
            ! select case(cache%info%charge)
            ! case(atom_resolved)
            !    call mixer%set(cache%mcache(mixer_index), pot%vat)
            ! case(shell_resolved)
            !    call mixer%set(cache%mcache(mixer_index), pot%vsh)
            ! end select

            ! select case(cache%info%dipole)
            ! case(atom_resolved)
            !    call mixer%set(cache%mcache(mixer_index), pot%vdp)
            ! end select

            ! select case(cache%info%quadrupole)
            ! case(atom_resolved)
            !    call mixer%set(cache%mcache(mixer_index), pot%vqp)
            ! end select

            ! select case(cache%info%density)
            ! case(orbital_resolved)
            !    call mixer%set(cache%mcache(mixer_index), wfn%coeff)  
            ! end select

            ! For the potential so far only based on the Fock matrix
            call mixer%set(cache%mcache(mixer_index), wfn%coeff)
         end if
      end associate
      cache%mcache(mixer_index)%raw%initialized = .true.
   end do

end subroutine set_mixer


subroutine diff_mixer(self, cache, wfn, pot, mode)
   !> Instance of the iterator
   class(iterator_type), intent(in) :: self
   !> Cache for the iterator
   type(iterator_cache), intent(inout) :: cache
   !> Tight-binding wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential shifts
   type(potential_type), intent(in) :: pot
   !> Mode of the mixer
   integer, intent(in) :: mode

   integer :: mixer_index, targets(2), it

   if (.not. allocated(self%mixers)) return

   targets = [cache%current_mixer, cache%precollect_mixer]

   do it = 1, 2
      mixer_index = targets(it)
      
      ! Skip invalid mixer indices
      if (mixer_index <= 0) cycle
      if (mixer_index > size(self%mixers)) cycle
      if (.not. allocated(cache%mcache(mixer_index)%raw)) cycle
      if (.not. cache%mcache(mixer_index)%raw%initialized) cycle

      ! Check if the mixer mode matches
      if (self%mixers(mixer_index)%raw%mode /= mode) cycle

      ! Reset chunk index
      cache%mcache(mixer_index)%raw%idif = 0

      associate(mixer => self%mixers(mixer_index)%raw)
         if (mode == mixer_mode%density) then
            select case(cache%info%charge)
            case(atom_resolved)
               call mixer%diff(cache%mcache(mixer_index), wfn%qat)
            case(shell_resolved)
               call mixer%diff(cache%mcache(mixer_index), wfn%qsh)
            end select

            select case(cache%info%dipole)
            case(atom_resolved)
               call mixer%diff(cache%mcache(mixer_index), wfn%dpat)
            end select

            select case(cache%info%quadrupole)
            case(atom_resolved)
               call mixer%diff(cache%mcache(mixer_index), wfn%qpat)
            end select

            select case(cache%info%density)
            case(orbital_resolved)
               call mixer%diff(cache%mcache(mixer_index), wfn%density)
            end select
         else if (mode == mixer_mode%potential) then
            ! select case(cache%info%charge)
            ! case(atom_resolved)
            !    call mixer%diff(cache%mcache(mixer_index), pot%vat)
            ! case(shell_resolved)
            !    call mixer%diff(cache%mcache(mixer_index), pot%vsh)
            ! end select

            ! select case(cache%info%dipole)
            ! case(atom_resolved)
            !    call mixer%diff(cache%mcache(mixer_index), pot%vdp)
            ! end select

            ! select case(cache%info%quadrupole)
            ! case(atom_resolved)
            !    call mixer%diff(cache%mcache(mixer_index), pot%vqp)
            ! end select

            ! select case(cache%info%density)
            ! case(orbital_resolved)
            !    call mixer%diff(cache%mcache(mixer_index), wfn%coeff)
            ! end select

            ! For the potential so far only based on the Fock matrix
            call mixer%diff(cache%mcache(mixer_index), wfn%coeff)
         end if
      end associate
   end do

end subroutine diff_mixer


subroutine next_mixer(self, cache, mode, error)
   !> Instance of iterator
   class(iterator_type), intent(in) :: self
   !> Cache for the iterator
   type(iterator_cache), intent(inout) :: cache
   !> Mode of the mixer
   integer, intent(in) :: mode
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   if (.not. allocated(self%mixers)) return
   if (.not. cache%mcache(cache%current_mixer)%raw%initialized) return

   if (self%mixers(cache%current_mixer)%raw%mode == mode) then
      call self%mixers(cache%current_mixer)%raw%next( &
         & cache%mcache(cache%current_mixer), error)
   end if

end subroutine next_mixer


subroutine get_mixer(self, cache, bas, wfn, pot, mode)
   !> Instance of the iterator
   class(iterator_type), intent(in) :: self
   !> Cache for the iterator
   type(iterator_cache), intent(inout) :: cache
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Tight-binding wavefunction data
   type(wavefunction_type), intent(inout) :: wfn
   !> Density dependent potential shifts
   type(potential_type), intent(inout) :: pot
   !> Mode of the mixer
   integer, intent(in) :: mode

   if (.not. allocated(self%mixers)) return
   if (.not. cache%mcache(cache%current_mixer)%raw%initialized) return

   ! Reset chunk index 
   cache%mcache(cache%current_mixer)%raw%iget = 0

   if (self%mixers(cache%current_mixer)%raw%mode == mode) then
      associate(mixer => self%mixers(cache%current_mixer)%raw)

         if (mode == mixer_mode%density) then
            select case(cache%info%charge)
            case(atom_resolved)
               call mixer%get(cache%mcache(cache%current_mixer), wfn%qat)
            case(shell_resolved)
               call mixer%get(cache%mcache(cache%current_mixer), wfn%qsh)
               call get_qat_from_qsh(bas, wfn%qsh, wfn%qat)
            end select

            select case(cache%info%dipole)
            case(atom_resolved)
               call mixer%get(cache%mcache(cache%current_mixer), wfn%dpat)
            end select

            select case(cache%info%quadrupole)
            case(atom_resolved)
               call mixer%get(cache%mcache(cache%current_mixer), wfn%qpat)
            end select

            select case(cache%info%density)
            case(orbital_resolved)
               call mixer%get(cache%mcache(cache%current_mixer), wfn%density)
            end select
         else if (mode == mixer_mode%potential) then
            ! select case(cache%info%charge)
            ! case(atom_resolved)
            !    call mixer%get(cache%mcache(cache%current_mixer), pot%vat)
            ! case(shell_resolved)
            !    call mixer%get(cache%mcache(cache%current_mixer), pot%vsh)
            ! end select

            ! select case(cache%info%dipole)
            ! case(atom_resolved)
            !    call mixer%get(cache%mcache(cache%current_mixer), pot%vdp)
            ! end select

            ! select case(cache%info%quadrupole)
            ! case(atom_resolved)
            !    call mixer%get(cache%mcache(cache%current_mixer), pot%vqp)
            ! end select

            ! select case(cache%info%density)
            ! case(orbital_resolved)
            !    call mixer%get(cache%mcache(cache%current_mixer), wfn%coeff)
            ! end select

            ! For the potential so far only based on the Fock matrix
            call mixer%get(cache%mcache(cache%current_mixer), wfn%coeff)
         end if
      end associate
   end if

end subroutine get_mixer


subroutine collect_mixer(self, cache, mode)
   !> Instance of the iterator
   class(iterator_type), intent(in) :: self
   !> Cache for the iterator
   type(iterator_cache), intent(inout) :: cache
   !> Mode of the mixer
   integer, intent(in) :: mode

   if (.not. allocated(self%mixers)) return

   if (cache%current_mixer > 0) then
      if (cache%mcache(cache%current_mixer)%raw%initialized) then
         if(self%mixers(cache%current_mixer)%raw%mode == mode) then
            call self%mixers(cache%current_mixer)%raw%collect( &
               & cache%mcache(cache%current_mixer))
         end if
      end if
   end if
   if (cache%precollect_mixer > 0) then
      if (cache%mcache(cache%precollect_mixer)%raw%initialized) then
         if (self%mixers(cache%precollect_mixer)%raw%mode == mode) then
            call self%mixers(cache%precollect_mixer)%raw%collect( &
               & cache%mcache(cache%precollect_mixer))
         end if
      end if
   end if
end subroutine collect_mixer


!> Remove cached mixer data
subroutine cleanup(self, cache)
   !> Instance of the iterator
   class(iterator_type), intent(in) :: self
   !> Cache for the iterator
   type(iterator_cache), intent(inout) :: cache

   integer :: imix

   if (.not. allocated(cache%mcache)) return

   do imix = 1, size(cache%mcache)
      if (.not. allocated(cache%mcache(imix)%raw)) cycle
      call self%mixers(imix)%raw%cleanup(cache%mcache(imix))
      deallocate(cache%mcache(imix)%raw)
   end do

end subroutine cleanup

end module tblite_scf_iterator
