! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! Lesser GNU General Public License for more details.
!
! You should have received a copy of the Lesser GNU General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

!> Definition of the abstract base dispersion model and generic implementations
module dftd4_model_type
   use dftd4_cache, only : dispersion_cache
   use dftd4_damping_type, only : damping_type, damping_twobody, damping_threebody
   use dftd4_integrator, only : integrator_type
   use dftd4_param_type, only : param_type
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type
   use multicharge, only : mchrg_model_type
   implicit none
   private

   public :: dispersion_model, dftd_models, d4_qmod


   !> Abstract base dispersion model to evaluate C6 coefficients
   type, abstract :: dispersion_model

      !> Number of atoms coupled to by pairwise parameters
      integer :: ncoup

      !> Number of frequency grid points for dynamic polarizabilities
      integer :: ngrid

      !> Charge scaling height
      real(wp) :: ga

      !> Charge scaling steepness
      real(wp) :: gc

      !> Effective nuclear charges
      real(wp), allocatable :: zeff(:)

      !> Chemical hardness
      real(wp), allocatable :: eta(:)

      !> Electronegativity
      real(wp), allocatable :: en(:)

      !> Covalent radii for coordination number
      real(wp), allocatable :: rcov(:)

      !> Expectation values for C8 extrapolation
      real(wp), allocatable :: r4r2(:)

      !> Number of reference systems
      integer, allocatable :: ref(:)

      !> Number of Gaussian weights for each reference
      integer, allocatable :: ngw(:, :)

      !> Reference coordination numbers
      real(wp), allocatable :: cn(:, :)

      !> Reference partial charges
      real(wp), allocatable :: q(:, :)

      !> Reference dynamic polarizabilities
      real(wp), allocatable :: aiw(:, :, :)

      !> Reference C6 coefficients
      real(wp), allocatable :: c6(:, :, :, :)

      !> Indicates whether three-body dispersion term is charge dependent
      logical :: charge_dep_3b = .false.

      !> Multicharge model
      class(mchrg_model_type), allocatable :: mchrg 

      !> Integrator for polarizabilities
      class(integrator_type), allocatable :: integrator

      !> Label identifying the dispersion model
      character(len=:), allocatable :: label

      !> Default two-body damping function
      integer :: default_damping_2b

      !> Default three-body damping function
      integer :: default_damping_3b

   contains

      !> Update cache with dispersion coefficients and properties
      procedure(update), deferred :: update

      !> Evaluate atomic polarizabilities from cache
      procedure(get_polarizabilities), deferred :: get_polarizabilities

      !> Get two-body dispersion coefficients for an atom pair
      procedure(get_2b_coeffs), deferred :: get_2b_coeffs

      !> Get two-body dispersion coefficients and derivatives for an atom pair
      procedure(get_2b_derivs), deferred :: get_2b_derivs

      !> Get three-body dispersion coefficients for an atom triple
      procedure(get_3b_coeffs), deferred :: get_3b_coeffs

      !> Get three-body dispersion coefficients and derivatives for an atom triple
      procedure(get_3b_derivs), deferred :: get_3b_derivs

      !> Calculate two-body dispersion energy and derivatives
      procedure :: get_dispersion2
      
      !> Calculate three-body dispersion energy and derivatives
      procedure :: get_dispersion3
      
      !> Calculate pairwise two-body dispersion energies and derivatives
      procedure :: get_pairwise_dispersion2
      
      !> Calculate pairwise three-body dispersion energies and derivatives
      procedure :: get_pairwise_dispersion3

      !> Calculate damping radius for two-body interactions
      procedure(get_2b_rdamp), deferred :: get_2b_rdamp

      !> Calculate damping radius for three-body interactions
      procedure(get_3b_rdamp), deferred :: get_3b_rdamp

   end type dispersion_model

   abstract interface

      !> Update dispersion cache
      subroutine update(self, mol, cache, cn, q, grad, only_c6)
         import dispersion_model, dispersion_cache, structure_type, wp
         !> Dispersion model
         class(dispersion_model), intent(in) :: self
         !> Molecular structure data
         type(structure_type), intent(in) :: mol
         !> Dispersion cache
         type(dispersion_cache), intent(inout) :: cache
         !> Coordination numbers
         real(wp), intent(in) :: cn(:)
         !> Atomic partial charges
         real(wp), intent(in) :: q(:)
         !> Whether to compute derivatives
         logical, intent(in), optional :: grad
         !> Whether to compute only C6 coefficients
         logical, intent(in), optional :: only_c6
      end subroutine update

      !> Calculate atomic polarizabilities
      subroutine get_polarizabilities(self, cache, alpha, alphaqq, &
         & dadcn, dadq, daqqdcn, daqqdq)
         import dispersion_model, dispersion_cache, wp
         !> Dispersion model
         class(dispersion_model), intent(in) :: self
         !> Dispersion cache
         type(dispersion_cache), intent(in) :: cache
         !> Static dipole-dipole polarizabilities
         real(wp), intent(out) :: alpha(:)
         !> Static quadrupole-quadrupole polarizabilities
         real(wp), intent(out) :: alphaqq(:)
         !> Derivative of dip-dip polarizibility w.r.t. the coordination numbers
         real(wp), intent(out), optional :: dadcn(:)
         !> Derivative of dip-dip polarizibility w.r.t. the atomic partial charges
         real(wp), intent(out), optional :: dadq(:)
         !> Derivative of quad-quad polarizibility w.r.t. the coordination numbers
         real(wp), intent(out), optional :: daqqdcn(:)
         !> Derivative of quad-quad polarizibility w.r.t. the atomic partial charges
         real(wp), intent(out), optional :: daqqdq(:)
      end subroutine get_polarizabilities

      !> Get two-body dispersion coefficients for atom pair ij
      subroutine get_2b_coeffs(self, cache, iat, jat, izp, jzp, c6, c8)
         import dispersion_model, dispersion_cache, wp
         !> Dispersion model
         class(dispersion_model), intent(in) :: self
         !> Dispersion cache
         type(dispersion_cache), intent(in) :: cache
         !> First atom index
         integer, intent(in) :: iat
         !> Second atom index
         integer, intent(in) :: jat
         !> Atomic number of first atom
         integer, intent(in) :: izp
         !> Atomic number of second atom
         integer, intent(in) :: jzp
         !> C6 dispersion coefficient
         real(wp), intent(out) :: c6
         !> C8 dispersion coefficient
         real(wp), intent(out) :: c8
      end subroutine

      !> Get two-body dispersion coefficients for atom pair ij with derivatives
      subroutine get_2b_derivs(self, cache, iat, jat, izp, jzp, c6, c8, &
         & dc6dcni, dc6dqi, dc6dcnj, dc6dqj, dc8dcni, dc8dqi, dc8dcnj, dc8dqj)
         import dispersion_model, dispersion_cache, wp
         !> Dispersion model
         class(dispersion_model), intent(in) :: self
         !> Dispersion cache
         type(dispersion_cache), intent(in) :: cache
         !> First atom index
         integer, intent(in) :: iat
         !> Second atom index
         integer, intent(in) :: jat
         !> Atomic number of first atom
         integer, intent(in) :: izp
         !> Atomic number of second atom
         integer, intent(in) :: jzp
         !> C6 dispersion coefficient
         real(wp), intent(out) :: c6
         !> C8 dispersion coefficient
         real(wp), intent(out) :: c8
         !> Derivative of C6 w.r.t the coordination number of atom i
         real(wp), intent(out) :: dc6dcni
         !> Derivative of C6 w.r.t the partial charge of atom i
         real(wp), intent(out) :: dc6dqi
         !> Derivative of C6 w.r.t the coordination number of atom i
         real(wp), intent(out) :: dc6dcnj
         !> Derivative of C6 w.r.t the partial charge of atom i
         real(wp), intent(out) :: dc6dqj
         !> Derivative of C8 w.r.t the coordination number of atom i
         real(wp), intent(out) :: dc8dcni
         !> Derivative of C8 w.r.t the partial charge of atom i
         real(wp), intent(out) :: dc8dqi
         !> Derivative of C8 w.r.t the coordination number of atom i
         real(wp), intent(out) :: dc8dcnj
         !> Derivative of C8 w.r.t the partial charge of atom i
         real(wp), intent(out) :: dc8dqj
      end subroutine

      !> Get three-body dispersion coefficient for atom triple ijk
      subroutine get_3b_coeffs(self, cache, iat, jat, kat, c9)
         import dispersion_model, dispersion_cache, wp
         !> Dispersion model
         class(dispersion_model), intent(in) :: self
         !> Dispersion cache
         type(dispersion_cache), intent(in) :: cache
         !> First atom index
         integer, intent(in) :: iat
         !> Second atom index
         integer, intent(in) :: jat
         !> Third atom index
         integer, intent(in) :: kat
         !> C9 dispersion coefficient
         real(wp), intent(out) :: c9
      end subroutine

      !> Get three-body dispersion coefficient for atom triple ijk with derivatives
      subroutine get_3b_derivs(self, cache, iat, jat, kat, c9, &
         & dc9dcni, dc9dqi, dc9dcnj, dc9dqj, dc9dcnk, dc9dqk)
         import dispersion_model, dispersion_cache, wp
         !> Dispersion model
         class(dispersion_model), intent(in) :: self
         !> Dispersion cache
         type(dispersion_cache), intent(in) :: cache
         !> First atom index
         integer, intent(in) :: iat
         !> Second atom index
         integer, intent(in) :: jat
         !> Third atom index
         integer, intent(in) :: kat
         !> C9 dispersion coefficient
         real(wp), intent(out) :: c9
         !> Derivative of C9 w.r.t the coordination number of atom i
         real(wp), intent(out) :: dc9dcni
         !> Derivative of C9 w.r.t the partial charge of atom i
         real(wp), intent(out) :: dc9dqi
         !> Derivative of C9 w.r.t the coordination number of atom j
         real(wp), intent(out) :: dc9dcnj
         !> Derivative of C9 w.r.t the partial charge of atom j
         real(wp), intent(out) :: dc9dqj
         !> Derivative of C9 w.r.t the coordination number of atom k
         real(wp), intent(out) :: dc9dcnk
         !> Derivative of C9 w.r.t the partial charge of atom k
         real(wp), intent(out) :: dc9dqk
      end subroutine

      !> Calculate damping radius for two-body interactions
      pure function get_2b_rdamp(self, izp, jzp) result(rdamp)
         import dispersion_model, wp
         !> Dispersion model
         class(dispersion_model), intent(in) :: self
         !> Atomic number of first atom
         integer, intent(in) :: izp
         !> Atomic number of second atom
         integer, intent(in) :: jzp
         !> Damping radius
         real(wp) :: rdamp
      end function get_2b_rdamp

      !> Calculate damping radius for three-body interactions
      pure function get_3b_rdamp(self, izp, jzp) result(rdamp)
         import dispersion_model, wp
         !> Dispersion model
         class(dispersion_model), intent(in) :: self
         !> Atomic number of first atom
         integer, intent(in) :: izp
         !> Atomic number of second atom
         integer, intent(in) :: jzp
         !> Damping radius
         real(wp) :: rdamp
      end function get_3b_rdamp

   end interface

   !> Possible DFT-D models
   type :: enum_dftd_models
      !> DFT-D4 model
      integer :: d4 = 1
      !> DFT-D4S model
      integer :: d4s = 2
      !> Reviset DFT-D4Srev model
      integer :: d4srev = 3
   end type enum_dftd_models

   !> Actual enumerator for possible DFT-D models
   type(enum_dftd_models), parameter :: dftd_models = enum_dftd_models()
   !DEC$ ATTRIBUTES DLLEXPORT :: dftd_models

   !> Possible reference charges for D4
   type :: enum_qmod

      !> Electronegativity equilibration charges
      integer :: eeq = 1

      !> GFN2-xTB Mulliken partial charges
      integer :: gfn2 = 2

      !> Bond-Capcity Electronegativity equilibration charges
      integer :: eeqbc = 3

      !> g-xTB Mulliken partial charges
      integer :: gxtb = 4

   end type enum_qmod

   !> Actual enumerator for D4 reference charges
   type(enum_qmod), parameter :: d4_qmod = enum_qmod()
   !DEC$ ATTRIBUTES DLLEXPORT :: d4_qmod

contains

!> Wrapper for the evaluation of two-body dispersion energy and derivatives
subroutine get_dispersion2(self, mol, cache, damp, param, trans, cutoff, energy, &
   & dEdcn, dEdq, gradient, sigma)
   !DEC$ ATTRIBUTES DLLEXPORT :: get_dispersion2
   !> Dispersion model
   class(dispersion_model), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Dispersion cache
   type(dispersion_cache), intent(in) :: cache
   !> Damping function
   type(damping_type), intent(in) :: damp
   !> Damping parameters
   type(param_type), intent(in) :: param
   !> Lattice translations
   real(wp), intent(in) :: trans(:, :)
   !> Real space cutoff
   real(wp), intent(in) :: cutoff
   !> Dispersion energy
   real(wp), intent(inout) :: energy(:)
   !> Derivative of the energy w.r.t. the coordination number
   real(wp), intent(inout), optional :: dEdcn(:)
   !> Derivative of the energy w.r.t. the partial charges
   real(wp), intent(inout), optional :: dEdq(:)
   !> Dispersion gradient
   real(wp), intent(inout), optional :: gradient(:, :)
   !> Dispersion virial
   real(wp), intent(inout), optional :: sigma(:, :)

   logical :: grad

   if (.not. allocated(damp%damping_2b) .or. .not. allocated(param%s6) .or. &
      & .not. allocated(param%s8)) return
   if ((abs(param%s6) < epsilon(1.0_wp) .and. abs(param%s8) < epsilon(1.0_wp))) return
   grad = present(gradient) .and. present(sigma) .and. present(dEdcn) .and. present(dEdq)

   if (grad) then
      call get_dispersion2_derivs(self, mol, cache, damp%damping_2b, param, &
         & trans, cutoff, energy, dEdcn, dEdq, gradient, sigma)
   else
      call get_dispersion2_energy(self, mol, cache, damp%damping_2b, param, &
         & trans, cutoff, energy)
   end if
end subroutine get_dispersion2


!> Evaluation of the two-body dispersion energy expression
subroutine get_dispersion2_energy(self, mol, cache, damp, param, trans, cutoff, &
   & energy)
   !> Dispersion model
   class(dispersion_model), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Dispersion cache
   type(dispersion_cache), intent(in) :: cache
   !> Two-body damping function
   class(damping_twobody), intent(in) :: damp
   !> Damping parameters
   type(param_type), intent(in) :: param
   !> Lattice translations
   real(wp), intent(in) :: trans(:, :)
   !> Real space cutoff
   real(wp), intent(in) :: cutoff
   !> Dispersion energy
   real(wp), intent(inout) :: energy(:)

   integer :: iat, jat, jtr, izp, jzp
   real(wp) :: vec(3), r2, cutoff2, c6, c8, rdamp, d6, d8, dE

   ! Thread-private array for reduction
   ! Set to 0 explicitly as the shared variants are potentially non-zero (inout)
   real(wp), allocatable :: energy_local(:)

   cutoff2 = cutoff*cutoff
   
   !$omp parallel default(none) &
   !$omp shared(mol, self, cache, damp, param, trans, cutoff2) &
   !$omp private(iat, jat, jtr, izp, jzp, vec, r2, c6, c8, rdamp, d6, d8, dE) &
   !$omp shared(energy) &
   !$omp private(energy_local)
   allocate(energy_local(size(energy, 1)), source=0.0_wp)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)
         call self%get_2b_coeffs(cache, iat, jat, izp, jzp, c6, c8)
         rdamp = self%get_2b_rdamp(izp, jzp)
         do jtr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, jtr))
            r2 = vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3)
            if (r2 > cutoff2 .or. r2 < epsilon(1.0_wp)) cycle

            call damp%get_2b_damp(param, r2, rdamp, d6, d8)

            dE = -0.5_wp * (c6 * d6 + c8 * d8)

            energy_local(iat) = energy_local(iat) + dE
            if (iat /= jat) then
               energy_local(jat) = energy_local(jat) + dE
            end if
         end do
      end do
   end do
   !$omp end do
   !$omp critical (get_dispersion2_energy_)
   energy(:) = energy + energy_local
   !$omp end critical (get_dispersion2_energy_)
   deallocate(energy_local)
   !$omp end parallel

end subroutine get_dispersion2_energy


!> Evaluation of the two-body dispersion energy and gradient expression
subroutine get_dispersion2_derivs(self, mol, cache, damp, param, trans, cutoff, &
   & energy, dEdcn, dEdq, gradient, sigma)
   !> Dispersion model
   class(dispersion_model), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Dispersion cache
   type(dispersion_cache), intent(in) :: cache
   !> Two-body damping function
   class(damping_twobody), intent(in) :: damp
   !> Damping parameters
   type(param_type), intent(in) :: param
   !> Lattice translations
   real(wp), intent(in) :: trans(:, :)
   !> Real space cutoff
   real(wp), intent(in) :: cutoff
   !> Dispersion energy
   real(wp), intent(inout) :: energy(:)
   !> Derivative of the energy w.r.t. the coordination number
   real(wp), intent(inout) :: dEdcn(:)
   !> Derivative of the energy w.r.t. the partial charges
   real(wp), intent(inout) :: dEdq(:)
   !> Dispersion gradient
   real(wp), intent(inout) :: gradient(:, :)
   !> Dispersion virial
   real(wp), intent(inout) :: sigma(:, :)

   integer :: iat, jat, jtr, izp, jzp
   real(wp) :: vec(3), r2, cutoff2, r, c6, c8, rdamp
   real(wp) :: d6, d8, d6dr, d8dr, gdisp
   real(wp) :: dE, dG(3), dS(3, 3)
   real(wp) :: dc6dcni, dc6dqi, dc6dcnj, dc6dqj
   real(wp) :: dc8dcni, dc8dqi, dc8dcnj, dc8dqj

   ! Thread-private arrays for reduction
   ! Set to 0 explicitly as the shared variants are potentially non-zero (inout)
   real(wp), allocatable :: energy_local(:)
   real(wp), allocatable :: dEdcn_local(:)
   real(wp), allocatable :: dEdq_local(:)
   real(wp), allocatable :: gradient_local(:, :)
   real(wp), allocatable :: sigma_local(:, :)

   cutoff2 = cutoff*cutoff
   
   !$omp parallel default(none) &
   !$omp shared(mol, self, cache, damp, param, trans, cutoff2) &
   !$omp private(iat, jat, jtr, izp, jzp, vec, r2, c6, c8, rdamp, &
   !$omp& d6, d8, d6dr, d8dr, dc6dcni, dc6dqi, dc6dcnj, dc6dqj, &
   !$omp& dc8dcni, dc8dqi, dc8dcnj, dc8dqj, gdisp, dE, dG, dS) &
   !$omp shared(energy, gradient, sigma, dEdcn, dEdq) &
   !$omp private(energy_local, gradient_local, sigma_local, dEdcn_local, dEdq_local)
   allocate(energy_local(size(energy, 1)), source=0.0_wp)
   allocate(dEdcn_local(size(energy, 1)), source=0.0_wp)
   allocate(dEdq_local(size(energy, 1)), source=0.0_wp)
   allocate(gradient_local(size(gradient, 1), size(gradient, 2)), source=0.0_wp)
   allocate(sigma_local(size(sigma, 1), size(sigma, 2)), source=0.0_wp)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)
         call self%get_2b_derivs(cache, iat, jat, izp, jzp, c6, c8, &
            & dc6dcni, dc6dqi, dc6dcnj, dc6dqj, dc8dcni, dc8dqi, dc8dcnj, dc8dqj)
         rdamp = self%get_2b_rdamp(izp, jzp)
         do jtr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, jtr))
            r2 = vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3)
            if (r2 > cutoff2 .or. r2 < epsilon(1.0_wp)) cycle

            call damp%get_2b_derivs(param, r2, rdamp, d6, d8, d6dr, d8dr)
            
            dE = -0.5_wp * (c6 * d6 + c8 * d8)
            
            gdisp = -(c6 * d6dr + c8 * d8dr)
            dG(:) = gdisp * vec(:)
            dS(:, :) = spread(dG, 1, 3) * spread(vec, 2, 3) * 0.5_wp

            energy_local(iat) = energy_local(iat) + dE
            dEdcn_local(iat) = dEdcn_local(iat) - (dc6dcni*d6 + dc8dcni*d8)
            dEdq_local(iat) = dEdq_local(iat) - (dc6dqi*d6 + dc8dqi*d8)
            sigma_local(:, :) = sigma_local + dS
            if (iat /= jat) then
               energy_local(jat) = energy_local(jat) + dE
               dEdcn_local(jat) = dEdcn_local(jat) - (dc6dcnj*d6 + dc8dcnj*d8)
               dEdq_local(jat) = dEdq_local(jat) - (dc6dqj*d6 + dc8dqj*d8)
               gradient_local(:, iat) = gradient_local(:, iat) + dG
               gradient_local(:, jat) = gradient_local(:, jat) - dG
               sigma_local(:, :) = sigma_local + dS
            end if
         end do
      end do
   end do
   !$omp end do
   !$omp critical (get_dispersion2_derivs_)
   energy(:) = energy + energy_local
   dEdcn(:) = dEdcn + dEdcn_local
   dEdq(:) = dEdq + dEdq_local
   gradient(:, :) = gradient + gradient_local
   sigma(:, :) = sigma + sigma_local
   !$omp end critical (get_dispersion2_derivs_)
   deallocate(energy_local)
   deallocate(dEdcn_local)
   deallocate(dEdq_local)
   deallocate(gradient_local)
   deallocate(sigma_local)
   !$omp end parallel

end subroutine get_dispersion2_derivs


!> Wrapper to handle the evaluation of three-body dispersion energy and derivatives
subroutine get_dispersion3(self, mol, cache, damp, param, trans, cutoff, energy, &
   & dEdcn, dEdq, gradient, sigma)
   !DEC$ ATTRIBUTES DLLEXPORT :: get_dispersion3
   !> Dispersion model
   class(dispersion_model), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Dispersion cache
   type(dispersion_cache), intent(in) :: cache
   !> Damping function
   type(damping_type), intent(in) :: damp
   !> Damping parameters
   type(param_type), intent(in) :: param
   !> Lattice translations
   real(wp), intent(in) :: trans(:, :)
   !> Real space cutoff
   real(wp), intent(in) :: cutoff
   !> Dispersion energy
   real(wp), intent(inout) :: energy(:)
   !> Derivative of the energy w.r.t. the coordination number
   real(wp), intent(inout), optional :: dEdcn(:)
   !> Derivative of the energy w.r.t. the partial charges
   real(wp), intent(inout), optional :: dEdq(:)
   !> Dispersion gradient
   real(wp), intent(inout), optional :: gradient(:, :)
   !> Dispersion virial
   real(wp), intent(inout), optional :: sigma(:, :)

   logical :: grad

   if (.not. allocated(damp%damping_3b) .or. .not.allocated(param%s9)) return
   if (abs(param%s9) < epsilon(1.0_wp)) return
   grad = present(gradient) .and. present(sigma) .and. present(dEdcn) .and. present(dEdq)

   if (grad) then
      call get_dispersion3_derivs(self, mol, cache, damp%damping_3b, param, &
         & trans, cutoff, energy, dEdcn, dEdq, gradient, sigma)
   else
      call get_dispersion3_energy(self, mol, cache, damp%damping_3b, param, &
         & trans, cutoff, energy)
   end if
end subroutine get_dispersion3


!> Evaluation of the three-body dispersion energy
subroutine get_dispersion3_energy(self, mol, cache, damp, param, trans, cutoff, energy)
   !> Dispersion model
   class(dispersion_model), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Dispersion cache
   type(dispersion_cache), intent(in) :: cache
   !> Three-body damping function
   class(damping_threebody), intent(in) :: damp
   !> Damping parameters
   type(param_type), intent(in) :: param
   !> Lattice translations
   real(wp), intent(in) :: trans(:, :)
   !> Real space cutoff
   real(wp), intent(in) :: cutoff
   !> Dispersion energy
   real(wp), intent(inout) :: energy(:)

   integer :: iat, jat, kat, izp, jzp, kzp, jtr, ktr
   real(wp) :: vij(3), vjk(3), vik(3), r2ij, r2jk, r2ik, r1, r2, r3, r5
   real(wp) :: rdamp, rdampij, rdampik, rdampjk, c9, d9, triple, ang, dE
   real(wp) :: cutoff2

   ! Thread-private arrays for reduction
   ! Set to 0 explicitly as the shared variants are potentially non-zero (inout)
   real(wp), allocatable :: energy_local(:)

   cutoff2 = cutoff*cutoff
   
   !$omp parallel default(none) &
   !$omp shared(mol, self, cache, damp, param, trans, cutoff2) &
   !$omp private(iat, jat, kat, izp, jzp, kzp, jtr, ktr, vij, vjk, vik, &
   !$omp& r1, r2, r3, r5, r2ij, r2jk, r2ik, rdamp, rdampij, rdampik, rdampjk, &
   !$omp& d9, ang, triple, c9, dE) &
   !$omp shared(energy) &
   !$omp private(energy_local)
   allocate(energy_local(size(energy, 1)), source=0.0_wp)
   !$omp do schedule(dynamic)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)
         rdampij = self%get_3b_rdamp(izp, jzp)
         do jtr = 1, size(trans, 2)
            vij(:) = mol%xyz(:, jat) + trans(:, jtr) - mol%xyz(:, iat)
            r2ij = vij(1)*vij(1) + vij(2)*vij(2) + vij(3)*vij(3)
            if (r2ij > cutoff2 .or. r2ij < epsilon(1.0_wp)) cycle
            do kat = 1, jat
               kzp = mol%id(kat)
               call self%get_3b_coeffs(cache, iat, jat, kat, c9)
               rdampik = self%get_3b_rdamp(izp, kzp)
               rdampjk = self%get_3b_rdamp(jzp, kzp)
               rdamp = rdampij * rdampik * rdampjk
               triple = triple_scale(iat, jat, kat)
               do ktr = 1, size(trans, 2)
                  vik(:) = mol%xyz(:, kat) + trans(:, ktr) - mol%xyz(:, iat)
                  r2ik = vik(1)*vik(1) + vik(2)*vik(2) + vik(3)*vik(3)
                  if (r2ik > cutoff2 .or. r2ik < epsilon(1.0_wp)) cycle
                  vjk(:) = mol%xyz(:, kat) + trans(:, ktr) - mol%xyz(:, jat) &
                     & - trans(:, jtr)
                  r2jk = vjk(1)*vjk(1) + vjk(2)*vjk(2) + vjk(3)*vjk(3)
                  if (r2jk > cutoff2 .or. r2jk < epsilon(1.0_wp)) cycle
                  r2 = r2ij*r2ik*r2jk
                  r1 = sqrt(r2)
                  r3 = r2 * r1
                  r5 = r3 * r2

                  call damp%get_3b_damp(param, r1, r2ij, r2ik, r2jk, &
                     & rdamp, rdampij, rdampik, rdampjk, d9)

                  ang = 0.375_wp*(r2ij + r2jk - r2ik)*(r2ij - r2jk + r2ik)&
                     & *(-r2ij + r2jk + r2ik) / r5 + 1.0_wp / r3

                  dE = c9 * triple * d9 * ang / 3.0_wp 
                  energy_local(iat) = energy_local(iat) - dE
                  energy_local(jat) = energy_local(jat) - dE
                  energy_local(kat) = energy_local(kat) - dE
               end do
            end do
         end do
      end do
   end do
   !$omp end do
   !$omp critical (get_atm_dispersion_energy_)
   energy(:) = energy + energy_local
   !$omp end critical (get_atm_dispersion_energy_)
   deallocate(energy_local)
   !$omp end parallel

end subroutine get_dispersion3_energy


!> Evaluation of the three-body dispersion energy and gradient expression
subroutine get_dispersion3_derivs(self, mol, cache, damp, param, trans, cutoff, &
   & energy, dEdcn, dEdq, gradient, sigma)
   !> Dispersion model
   class(dispersion_model), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Dispersion cache
   type(dispersion_cache), intent(in) :: cache
   !> Three-body damping function
   class(damping_threebody), intent(in) :: damp
   !> Damping parameters
   type(param_type), intent(in) :: param
   !> Lattice translations
   real(wp), intent(in) :: trans(:, :)
   !> Real space cutoff
   real(wp), intent(in) :: cutoff
   !> Dispersion energy
   real(wp), intent(inout) :: energy(:)
   !> Derivative of the energy w.r.t. the coordination number
   real(wp), intent(inout) :: dEdcn(:)
   !> Derivative of the energy w.r.t. the partial charges
   real(wp), intent(inout) :: dEdq(:)
   !> Dispersion gradient
   real(wp), intent(inout) :: gradient(:, :)
   !> Dispersion virial
   real(wp), intent(inout) :: sigma(:, :)

   integer :: iat, jat, kat, izp, jzp, kzp, jtr, ktr
   real(wp) :: vij(3), vjk(3), vik(3), r2ij, r2jk, r2ik, r1, r2, r3, r5
   real(wp) :: rdamp, rdampij, rdampik, rdampjk
   real(wp) :: c9, dc9dcni, dc9dqi, dc9dcnj, dc9dqj, dc9dcnk, dc9dqk
   real(wp) :: d9, d9drij, d9drik, d9drjk, ang, dang, triple
   real(wp) :: cutoff2
   real(wp) :: dE, dGij(3), dGik(3), dGjk(3), dS(3,3)
   
   ! Thread-private arrays for reduction
   ! Set to 0 explicitly as the shared variants are potentially non-zero (inout)
   real(wp), allocatable :: energy_local(:), dEdcn_local(:), dEdq_local(:)
   real(wp), allocatable :: gradient_local(:, :), sigma_local(:, :)

   cutoff2 = cutoff*cutoff

   !$omp parallel default(none) &
   !$omp shared(mol, self, cache, damp, param, trans, cutoff2) &
   !$omp private(iat, jat, kat, izp, jzp, kzp, jtr, ktr, vij, vjk, vik, &
   !$omp& r1, r2, r3, r5, r2ij, r2jk, r2ik, rdamp, rdampij, rdampik, rdampjk, &
   !$omp& d9, d9drij, d9drik, d9drjk, ang, dang, triple, &
   !$omp& c9, dc9dcni, dc9dqi, dc9dcnj, dc9dqj, dc9dcnk, dc9dqk, &
   !$omp& dE, dGij, dGik, dGjk, dS) &
   !$omp shared(energy, gradient, sigma, dEdcn, dEdq) &
   !$omp private(energy_local, dEdcn_local, dEdq_local, gradient_local, sigma_local)
   allocate(energy_local(size(energy, 1)), source=0.0_wp)
   allocate(dEdcn_local(size(energy, 1)), source=0.0_wp)
   allocate(dEdq_local(size(energy, 1)), source=0.0_wp)
   allocate(gradient_local(3, size(energy, 1)), source=0.0_wp)
   allocate(sigma_local(3, 3), source=0.0_wp)
   !$omp do schedule(dynamic)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)
         rdampij = self%get_3b_rdamp(izp, jzp)
         do jtr = 1, size(trans, 2)
            vij(:) = mol%xyz(:, jat) + trans(:, jtr) - mol%xyz(:, iat)
            r2ij = vij(1)*vij(1) + vij(2)*vij(2) + vij(3)*vij(3)
            if (r2ij > cutoff2 .or. r2ij < epsilon(1.0_wp)) cycle
            do kat = 1, jat
               kzp = mol%id(kat)
               call self%get_3b_derivs(cache, iat, jat, kat, c9, &
                  & dc9dcni, dc9dqi, dc9dcnj, dc9dqj, dc9dcnk, dc9dqk)
               rdampik = self%get_3b_rdamp(izp, kzp)
               rdampjk = self%get_3b_rdamp(jzp, kzp)
               rdamp = rdampij * rdampik * rdampjk
               triple = triple_scale(iat, jat, kat)
               do ktr = 1, size(trans, 2)
                  vik(:) = mol%xyz(:, kat) + trans(:, ktr) - mol%xyz(:, iat)
                  r2ik = vik(1)*vik(1) + vik(2)*vik(2) + vik(3)*vik(3)
                  if (r2ik > cutoff2 .or. r2ik < epsilon(1.0_wp)) cycle                  
                  vjk(:) = mol%xyz(:, kat) + trans(:, ktr) - mol%xyz(:, jat) &
                     & - trans(:, jtr)
                  r2jk = vjk(1)*vjk(1) + vjk(2)*vjk(2) + vjk(3)*vjk(3)
                  if (r2jk > cutoff2 .or. r2jk < epsilon(1.0_wp)) cycle
                  r2 = r2ij*r2ik*r2jk
                  r1 = sqrt(r2)
                  r3 = r2 * r1
                  r5 = r3 * r2

                  call damp%get_3b_derivs(param, r1, r2ij, r2ik, r2jk, rdamp, &
                        & rdampij, rdampik, rdampjk, d9, d9drij, d9drik, d9drjk)
                  ang = 0.375_wp*(r2ij + r2jk - r2ik)*(r2ij - r2jk + r2ik)&
                     & *(-r2ij + r2jk + r2ik) / r5 + 1.0_wp / r3
                  
                  ! d/drij
                  dang = -0.375_wp * (r2ij**3 + r2ij**2 * (r2jk + r2ik) &
                     & + r2ij * (3.0_wp * r2jk**2 + 2.0_wp * r2jk*r2ik &
                     & + 3.0_wp * r2ik**2)&
                     & - 5.0_wp * (r2jk - r2ik)**2 * (r2jk + r2ik)) / r5
                  dGij(:) = c9 * (-dang * d9 / r2ij - ang * d9drij) * vij

                  ! d/drik
                  dang = -0.375_wp * (r2ik**3 + r2ik**2 * (r2jk + r2ij) &
                     & + r2ik * (3.0_wp * r2jk**2 + 2.0_wp * r2jk * r2ij &
                     & + 3.0_wp * r2ij**2) &
                     & - 5.0_wp * (r2jk - r2ij)**2 * (r2jk + r2ij)) / r5
                  dGik(:) = c9 * (-dang * d9 / r2ik - ang * d9drik) * vik

                  ! d/drjk
                  dang = -0.375_wp * (r2jk**3 + r2jk**2*(r2ik + r2ij) &
                     & + r2jk * (3.0_wp * r2ik**2 + 2.0_wp * r2ik * r2ij &
                     & + 3.0_wp * r2ij**2) &
                     & - 5.0_wp * (r2ik - r2ij)**2 * (r2ik + r2ij)) / r5
                  dGjk(:) = c9 * (-dang * d9 / r2jk - ang * d9drjk) * vjk

                  dE = triple * d9 * ang / 3.0_wp
                  energy_local(iat) = energy_local(iat) - c9 * dE
                  energy_local(jat) = energy_local(jat) - c9 * dE
                  energy_local(kat) = energy_local(kat) - c9 * dE

                  gradient_local(:, iat) = gradient_local(:, iat) - dGij - dGik
                  gradient_local(:, jat) = gradient_local(:, jat) + dGij - dGjk
                  gradient_local(:, kat) = gradient_local(:, kat) + dGik + dGjk
                  
                  dS(:, :) = spread(dGij, 1, 3) * spread(vij, 2, 3)&
                        & + spread(dGik, 1, 3) * spread(vik, 2, 3)&
                        & + spread(dGjk, 1, 3) * spread(vjk, 2, 3)
                  sigma_local(:, :) = sigma_local + dS * triple

                  dE = dE * 3.0_wp 
                  dEdcn_local(iat) = dEdcn_local(iat) - dE * dc9dcni
                  dEdcn_local(jat) = dEdcn_local(jat) - dE * dc9dcnj
                  dEdcn_local(kat) = dEdcn_local(kat) - dE * dc9dcnk

                  dEdq_local(iat) = dEdq_local(iat) - dE * dc9dqi
                  dEdq_local(jat) = dEdq_local(jat) - dE * dc9dqj
                  dEdq_local(kat) = dEdq_local(kat) - dE * dc9dqk
               end do
            end do
         end do
      end do
   end do
   !$omp end do
   !$omp critical (get_atm_dispersion_derivs_)
   energy(:) = energy + energy_local
   dEdcn(:) = dEdcn + dEdcn_local
   dEdq(:) = dEdq + dEdq_local
   gradient(:, :) = gradient + gradient_local
   sigma(:, :) = sigma + sigma_local
   !$omp end critical (get_atm_dispersion_derivs_)
   deallocate(energy_local)
   deallocate(dEdcn_local)
   deallocate(dEdq_local)
   deallocate(gradient_local)
   deallocate(sigma_local)
   !$omp end parallel

end subroutine get_dispersion3_derivs


!> Calculate pairwise two-body dispersion energy
subroutine get_pairwise_dispersion2(self, mol, cache, damp, param, trans, cutoff, energy)
   !> Dispersion model
   class(dispersion_model), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Dispersion cache
   type(dispersion_cache), intent(in) :: cache
   !> Two-body damping function
   class(damping_twobody), intent(in) :: damp
   !> Damping parameters
   type(param_type), intent(in) :: param
   !> Lattice translations
   real(wp), intent(in) :: trans(:, :)
   !> Real space cutoff
   real(wp), intent(in) :: cutoff
   !> Pairwise dispersion energy
   real(wp), intent(inout) :: energy(:, :)
   
   integer :: iat, jat, izp, jzp, jtr
   real(wp) :: vec(3), r2, cutoff2, c6, c8, rdamp, d6, d8, dE

   ! Thread-private array for reduction
   ! Set to 0 explicitly as the shared variants are potentially non-zero (inout)
   real(wp), allocatable :: energy_local(:, :)
   
   if (.not. allocated(param%s6) .or. .not. allocated(param%s8)) return
   if ((abs(param%s6) < epsilon(1.0_wp) .and. abs(param%s8) < epsilon(1.0_wp))) return

   cutoff2 = cutoff*cutoff
   !$omp parallel default(none) &
   !$omp shared(mol, self, cache, damp, param, trans, cutoff2) &
   !$omp private(iat, jat, jtr, izp, jzp, vec, r2, c6, c8, rdamp, d6, d8, dE) &
   !$omp shared(energy) &
   !$omp private(energy_local)
   allocate(energy_local(size(energy, 1), size(energy, 2)), source=0.0_wp)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)
         call self%get_2b_coeffs(cache, iat, jat, izp, jzp, c6, c8)
         rdamp = self%get_2b_rdamp(izp, jzp)
         do jtr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, jtr))
            r2 = vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3)
            if (r2 > cutoff2 .or. r2 < epsilon(1.0_wp)) cycle
            
            call damp%get_2b_damp(param, r2, rdamp, d6, d8)

            dE = -0.5_wp * (c6 * d6 + c8 * d8)

            energy_local(jat, iat) = energy_local(jat, iat) + dE
            if (iat /= jat) then
               energy_local(iat, jat) = energy_local(iat, jat) + dE
            end if
         end do
      end do
   end do
   !$omp end do
   !$omp critical (get_pairwise_dispersion2_)
   energy(:, :) = energy + energy_local
   !$omp end critical (get_pairwise_dispersion2_)
   deallocate(energy_local)
   !$omp end parallel

end subroutine get_pairwise_dispersion2

!> Calculate pairwise three-body dispersion energy
subroutine get_pairwise_dispersion3(self, mol, cache, damp, param, trans, cutoff, energy)
   !> Dispersion model
   class(dispersion_model), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Dispersion cache
   type(dispersion_cache), intent(in) :: cache
   !> Three-body damping function
   class(damping_threebody), intent(in) :: damp
   !> Damping parameters
   type(param_type), intent(in) :: param
   !> Lattice translations
   real(wp), intent(in) :: trans(:, :)
   !> Real space cutoff
   real(wp), intent(in) :: cutoff
   !> Pairwise dispersion energy
   real(wp), intent(inout) :: energy(:, :)
   
   integer :: iat, jat, kat, izp, jzp, kzp, jtr, ktr
   real(wp) :: vij(3), vjk(3), vik(3), r2ij, r2jk, r2ik, r1, r2, r3, r5
   real(wp) :: rdamp, rdampij, rdampik, rdampjk, c9, d9, triple, ang, dE
   real(wp) :: cutoff2

   ! Thread-private arrays for reduction
   ! Set to 0 explicitly as the shared variants are potentially non-zero (inout)
   real(wp), allocatable :: energy_local(:, :)

   if (.not.allocated(param%s9)) return
   if (abs(param%s9) < epsilon(1.0_wp)) return

   cutoff2 = cutoff*cutoff
   !$omp parallel default(none) &
   !$omp shared(mol, self, cache, damp, param, trans, cutoff2) &
   !$omp private(iat, jat, kat, izp, jzp, kzp, jtr, ktr, vij, vjk, vik, &
   !$omp& r1, r2, r3, r5, r2ij, r2jk, r2ik, rdamp, rdampij, rdampik, rdampjk, &
   !$omp& d9, ang, triple, c9, dE) &
   !$omp shared(energy) &
   !$omp private(energy_local)
   allocate(energy_local(size(energy, 1), size(energy, 2)), source=0.0_wp)
   !$omp do schedule(dynamic)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)
         rdampij = self%get_3b_rdamp(izp, jzp)
         do jtr = 1, size(trans, 2)
            vij(:) = mol%xyz(:, jat) + trans(:, jtr) - mol%xyz(:, iat)
            r2ij = vij(1)*vij(1) + vij(2)*vij(2) + vij(3)*vij(3)
            if (r2ij > cutoff2 .or. r2ij < epsilon(1.0_wp)) cycle
            do kat = 1, jat
               kzp = mol%id(kat)
               call self%get_3b_coeffs(cache, iat, jat, kat, c9)
               rdampik = self%get_3b_rdamp(izp, kzp)
               rdampjk = self%get_3b_rdamp(jzp, kzp)
               rdamp = rdampij * rdampik * rdampjk
               triple = triple_scale(iat, jat, kat)
               do ktr = 1, size(trans, 2)
                  vik(:) = mol%xyz(:, kat) + trans(:, ktr) - mol%xyz(:, iat)
                  r2ik = vik(1)*vik(1) + vik(2)*vik(2) + vik(3)*vik(3)
                  if (r2ik > cutoff2 .or. r2ik < epsilon(1.0_wp)) cycle
                  vjk(:) = mol%xyz(:, kat) + trans(:, ktr) - mol%xyz(:, jat) &
                     & - trans(:, jtr)
                  r2jk = vjk(1)*vjk(1) + vjk(2)*vjk(2) + vjk(3)*vjk(3)
                  if (r2jk > cutoff2 .or. r2jk < epsilon(1.0_wp)) cycle
                  r2 = r2ij*r2ik*r2jk
                  r1 = sqrt(r2)
                  r3 = r2 * r1
                  r5 = r3 * r2

                  call damp%get_3b_damp(param, r1, r2ij, r2ik, r2jk, &
                     & rdamp, rdampij, rdampik, rdampjk, d9)

                  ang = 0.375_wp*(r2ij + r2jk - r2ik)*(r2ij - r2jk + r2ik)&
                     & *(-r2ij + r2jk + r2ik) / r5 + 1.0_wp / r3

                  dE = c9 * triple * d9 * ang / 6.0_wp 
                  energy_local(jat, iat) = energy_local(jat, iat) - dE
                  energy_local(kat, iat) = energy_local(kat, iat) - dE
                  energy_local(iat, jat) = energy_local(iat, jat) - dE
                  energy_local(kat, jat) = energy_local(kat, jat) - dE
                  energy_local(iat, kat) = energy_local(iat, kat) - dE
                  energy_local(jat, kat) = energy_local(jat, kat) - dE
               end do
            end do
         end do
      end do
   end do
   !$omp end do
   !$omp critical (get_pairwise_dispersion3_)
   energy(:, :) = energy + energy_local
   !$omp end critical (get_pairwise_dispersion3_)
   deallocate(energy_local)
   !$omp end parallel

end subroutine get_pairwise_dispersion3

!> Logic exercise to distribute a triple energy to atomwise energies.
elemental function triple_scale(ii, jj, kk) result(triple)

   !> Atom indices
   integer, intent(in) :: ii, jj, kk

   !> Fraction of energy
   real(wp) :: triple

   if (ii == jj) then
      if (ii == kk) then
         ! ii'i" -> 1/6
         triple = 1.0_wp/6.0_wp
      else
         ! ii'j -> 1/2
         triple = 0.5_wp
      end if
   else
      if (ii /= kk .and. jj /= kk) then
         ! ijk -> 1 (full)
         triple = 1.0_wp
      else
         ! ijj' and iji' -> 1/2
         triple = 0.5_wp
      end if
   end if

end function triple_scale

end module dftd4_model_type
