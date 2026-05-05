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

!> @file tblite/solvation/shift.f90
!> Provides solvent and state specific shifts to the free energy.

module tblite_solvation_shift
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use tblite_container_cache, only : container_cache
   use tblite_scf_potential, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_solvation_type, only : solvation_type
   use mctc_io_codata2018, only : molar_gas_constant, atomic_unit_of_energy
   use mctc_io_constants, only : codata
   use mctc_io_convert, only : kjtoau
   implicit none
   private
 
   public :: shift_solvation, new_shift, shift_input, solution_state

   type :: shift_input
      !> Reference state
      integer :: state = 1
      !> Solvent specific free energy shift
      real(wp) :: gshift = 0.0_wp
      !> Molar mass of solvent
      real(wp), allocatable :: molar_mass 
      !> Density of solvent
      real(wp), allocatable :: rho
      !> Temperature in Kelvin
      real(wp) :: temperature = 298.15_wp
      !> Linearized Poisson-Boltzmann model for parameter selection
      logical :: alpb = .true.
      !> Solvent for parameter selection
      character(len=:), allocatable :: solvent
   end type shift_input

   !> Provide constructor for shift input
   interface shift_input
      module procedure :: create_shift_input
   end interface shift_input

   !> Possible reference states for the solution
   type :: enum_solution_state
      !> 1 l of ideal gas in 1 l of liquid solution
      integer :: gsolv = 1
      !> 1 bar of ideal gas and 1 mol/L of liquid solution
      integer :: bar1mol = 2
      !> 1 bar of ideal gas to 1 mol/L of liquid solution at infinite dilution
      integer :: reference = 3
   end type enum_solution_state

   !> Actual solvation state enumerator
   type(enum_solution_state), parameter :: solution_state = enum_solution_state()

   real(wp), parameter :: refDensity = 1.0e-3_wp         ! kg__au/(1.0e10_wp*AA__Bohr)**3
   real(wp), parameter :: refPressure = 1.0e2_wp         ! kPa__au
   real(wp), parameter :: refMolecularMass = 1.0_wp      ! amu__au
   real(wp), parameter :: ambientTemperature = 298.15_wp

   type, extends(solvation_type) :: shift_solvation
      !> total shift
      real(wp) :: total_shift 
   contains
      !> Update cache from container
      procedure :: update
      !> Get shift energy
      procedure :: get_engrad
   end type shift_solvation

   !> Provide constructor for shift
   interface shift_solvation
      module procedure :: create_shift
   end interface shift_solvation

   type :: shift_cache
      !> Total free energy shift
      real(wp) :: total_shift
   end type shift_cache

   !> Identifier for container 
   character(len=*), parameter :: label = "empirical free energy shift and standard state/solvent correction"

contains


!> Consturctor for shift input to properly assign allocatable strings
function create_shift_input(state, solvent, alpb, gshift, molar_mass, rho, temperature) result(self)
   !> Reference state
   integer, intent(in), optional :: state
   !> Solvent for parameter selection
   character(len=*), intent(in), optional :: solvent
   !> Use analytical linearized Poisson-Boltzmann model
   logical, intent(in), optional :: alpb
   !> Solvent specific free energy shift
   real(wp), intent(in), optional :: gshift
   !> Molar mass of solvent
   real(wp), allocatable, intent(in), optional :: molar_mass 
   !> Density of solvent
   real(wp), allocatable, intent(in), optional :: rho
   !> Temperature in Kelvin
   real(wp), intent(in), optional :: temperature

   type(shift_input) :: self

   if (present(state)) then
      self%state = state
   end if

   if (present(solvent)) then 
      self%solvent = solvent
   end if

   if (present(alpb)) then 
      self%alpb = alpb
   end if

   if (present(gshift)) then 
      self%gshift = gshift
   end if

   if (present(molar_mass)) then 
      self%molar_mass = molar_mass
   end if

   if (present(rho)) then
      self%rho = rho
   end if

   if (present(temperature)) then
      self%temperature = temperature
   end if

end function create_shift_input

!> Calculate the solvent and state shift
subroutine new_shift(self, input)
   !> Instance of the solvation model
   type(shift_solvation) :: self
   !> Input for shift solvation
   type(shift_input), intent(in) :: input

   !> State shift
   real(wp) :: stateshift

   self%label = label
   stateshift = get_stateshift(input%state, input%temperature, input%rho, input%molar_mass)
   self%total_shift = input%gshift + stateshift 
end subroutine new_shift

!> Type constructor for shift
function create_shift(input) result(self)
   !> Instance of the solvation model
   type(shift_solvation) :: self
   !> Input for solvation shift
   type(shift_input), intent(in) :: input

   call new_shift(self, input)
end function create_shift

!> Update cache from container
subroutine update(self, mol, cache)
   !> Instance of the solvation model
   class(shift_solvation), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache

   type(shift_cache), pointer :: ptr

   call taint(cache, ptr)
   ptr%total_shift = self%total_shift
end subroutine update

!> Evaluate non-selfconsistent part of the interaction
subroutine get_engrad(self, mol, cache, energies, gradient, sigma)
   !> Instance of the solvation model
   class(shift_solvation), intent(in) :: self
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

   type(shift_cache), pointer :: ptr

   call view(cache, ptr)

   energies(:) = energies + ptr%total_shift/real(mol%nat)
end subroutine get_engrad

!> Calculate the solvation state shift
function get_stateshift(state, temperature, density, molecularMass) &
   & result(stateshift)
   !> Reference state
   integer, intent(in) :: state
   !> Temperature of the solution
   real(wp), intent(in) :: temperature
   !> Mass density of the solvent
   real(wp), intent(in) :: density
   !> Molecular mass of the solvent
   real(wp), intent(in) :: molecularMass
   !> Resulting shift to the solvation free energy
   real(wp) :: stateshift

   ! Boltzmann constant in au/K
   real(wp) :: boltzmann = codata%kb / atomic_unit_of_energy
   ! Ideal gas molar volume in m^3/mol 
   real(wp) :: idealGasMolVolume = molar_gas_constant * ambientTemperature / refPressure

   select case(state)
   case default
      stateshift = 0.0_wp
   case(solution_state%bar1mol)
      stateshift = temperature * boltzmann &
         & * log(idealGasMolVolume * temperature / ambientTemperature)
   case(solution_state%reference)
      stateshift = temperature * boltzmann &
         & * (log(idealGasMolVolume * temperature / ambientTemperature) &
         & + log(density/refDensity * refMolecularMass/molecularMass))
   end select

end function get_stateshift

subroutine taint(cache, ptr)
   type(container_cache), target, intent(inout) :: cache
   type(shift_cache), pointer, intent(out) :: ptr

   if (allocated(cache%raw)) then
      call view(cache, ptr)
      if (associated(ptr)) return
      deallocate(cache%raw)
   end if

   if (.not.allocated(cache%raw)) then
      block
         type(shift_cache), allocatable :: tmp
         allocate(tmp)
         call move_alloc(tmp, cache%raw)
      end block
   end if

   call view(cache, ptr)
end subroutine taint

subroutine view(cache, ptr)
   type(container_cache), target, intent(inout) :: cache
   type(shift_cache), pointer, intent(out) :: ptr
   nullify(ptr)
   select type(target => cache%raw)
   type is(shift_cache)
      ptr => target
   end select
end subroutine view

end module tblite_solvation_shift
