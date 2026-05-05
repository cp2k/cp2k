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

subroutine get_disp2_switch(cutoff, inner, active)

   real(wp), intent(in) :: cutoff
   real(wp), intent(out) :: inner
   logical, intent(out) :: active

   character(len=64) :: env
   integer :: stat, io
   real(wp) :: width

   inner = cutoff
   active = .false.

   call get_environment_variable("DFTD4_DISP2_SMOOTH_WIDTH", env, status=stat)
   if (stat /= 0) then
      call get_environment_variable("TBLITE_D4_DISP2_SMOOTH_WIDTH", env, status=stat)
   end if
   if (stat == 0) then
      read(env, *, iostat=io) width
      active = io == 0 .and. width > 0.0_wp .and. width < cutoff
      if (active) inner = cutoff - width
   end if

end subroutine get_disp2_switch


subroutine get_disp3_switch(cutoff, inner, active)

   real(wp), intent(in) :: cutoff
   real(wp), intent(out) :: inner
   logical, intent(out) :: active

   character(len=64) :: env
   integer :: stat, io
   real(wp) :: width

   inner = cutoff
   active = .false.

   call get_environment_variable("DFTD4_DISP3_SMOOTH_WIDTH", env, status=stat)
   if (stat /= 0) then
      call get_environment_variable("TBLITE_D4_DISP3_SMOOTH_WIDTH", env, status=stat)
   end if
   if (stat == 0) then
      read(env, *, iostat=io) width
      active = io == 0 .and. width > 0.0_wp .and. width < cutoff
      if (active) inner = cutoff - width
   end if

end subroutine get_disp3_switch


pure subroutine smooth_cutoff(r, cutoff, inner, active, sw, dswdr)

   real(wp), intent(in) :: r, cutoff, inner
   logical, intent(in) :: active
   real(wp), intent(out) :: sw, dswdr

   real(wp) :: width, x

   if (.not. active .or. r <= inner) then
      sw = 1.0_wp
      dswdr = 0.0_wp
   else if (r >= cutoff) then
      sw = 0.0_wp
      dswdr = 0.0_wp
   else
      width = cutoff - inner
      x = (cutoff - r) / width
      sw = x**3 * (10.0_wp + x*(-15.0_wp + 6.0_wp*x))
      dswdr = -30.0_wp * x**2 * (1.0_wp - x)**2 / width
   end if

end subroutine smooth_cutoff


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
   logical :: use_switch
   real(wp) :: vec(3), r2, r, cutoff2, cutoff_inner, c6, c8, rdamp, d6, d8
   real(wp) :: dE, sw, dswdr

   cutoff2 = cutoff*cutoff
   call get_disp2_switch(cutoff, cutoff_inner, use_switch)
   
   !$omp parallel do schedule(runtime) default(none) reduction(+:energy) &
   !$omp shared(mol, self, cache, damp, param, trans, cutoff, cutoff2, cutoff_inner, use_switch) &
   !$omp private(iat, jat, jtr, izp, jzp, vec, r2, r, c6, c8, rdamp, d6, d8, dE, sw, dswdr)
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
            if (use_switch) then
               r = sqrt(r2)
               call smooth_cutoff(r, cutoff, cutoff_inner, use_switch, sw, dswdr)
               if (sw <= 0.0_wp) cycle
            else
               sw = 1.0_wp
            end if

            call damp%get_2b_damp(param, r2, rdamp, d6, d8)

            dE = -0.5_wp * sw * (c6 * d6 + c8 * d8)

            energy(iat) = energy(iat) + dE
            if (iat /= jat) then
               energy(jat) = energy(jat) + dE
            end if
         end do
      end do
   end do

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
   logical :: use_switch
   real(wp) :: vec(3), r2, r, cutoff2, cutoff_inner, c6, c8, rdamp
   real(wp) :: d6, d8, d6dr, d8dr, edisp0, gdisp0, gdisp
   real(wp) :: dE, dG(3), dS(3, 3), sw, dswdr
   real(wp) :: dc6dcni, dc6dqi, dc6dcnj, dc6dqj
   real(wp) :: dc8dcni, dc8dqi, dc8dcnj, dc8dqj

   cutoff2 = cutoff*cutoff
   call get_disp2_switch(cutoff, cutoff_inner, use_switch)
   
   !$omp parallel do schedule(runtime) default(none) &
   !$omp reduction(+:energy, gradient, sigma, dEdcn, dEdq) &
   !$omp shared(mol, self, cache, damp, param, trans, cutoff, cutoff2, cutoff_inner, use_switch) &
   !$omp private(iat, jat, jtr, izp, jzp, vec, r2, c6, c8, rdamp, &
   !$omp& d6, d8, d6dr, d8dr, dc6dcni, dc6dqi, dc6dcnj, dc6dqj, &
   !$omp& dc8dcni, dc8dqi, dc8dcnj, dc8dqj, edisp0, gdisp0, gdisp, dE, dG, dS, &
   !$omp& r, sw, dswdr)
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
            if (use_switch) then
               r = sqrt(r2)
               call smooth_cutoff(r, cutoff, cutoff_inner, use_switch, sw, dswdr)
               if (sw <= 0.0_wp) cycle
            else
               r = 1.0_wp
               sw = 1.0_wp
               dswdr = 0.0_wp
            end if

            call damp%get_2b_derivs(param, r2, rdamp, d6, d8, d6dr, d8dr)
            
            edisp0 = c6 * d6 + c8 * d8
            dE = -0.5_wp * sw * edisp0
            
            gdisp0 = -(c6 * d6dr + c8 * d8dr)
            gdisp = sw * gdisp0 - edisp0 * dswdr / r
            dG(:) = gdisp * vec(:)
            dS(:, :) = spread(dG, 1, 3) * spread(vec, 2, 3) * 0.5_wp

            energy(iat) = energy(iat) + dE
            dEdcn(iat) = dEdcn(iat) - sw * (dc6dcni*d6 + dc8dcni*d8)
            dEdq(iat) = dEdq(iat) - sw * (dc6dqi*d6 + dc8dqi*d8)
            sigma(:, :) = sigma(:, :) + dS
            if (iat /= jat) then
               energy(jat) = energy(jat) + dE
               dEdcn(jat) = dEdcn(jat) - sw * (dc6dcnj*d6 + dc8dcnj*d8)
               dEdq(jat) = dEdq(jat) - sw * (dc6dqj*d6 + dc8dqj*d8)
               gradient(:, iat) = gradient(:, iat) + dG
               gradient(:, jat) = gradient(:, jat) - dG
               sigma(:, :) = sigma(:, :) + dS
            end if
         end do
      end do
   end do

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
   logical :: use_switch
   real(wp) :: vij(3), vjk(3), vik(3), r2ij, r2jk, r2ik, rij, rjk, rik
   real(wp) :: r1, r2, r3, r5
   real(wp) :: rdamp, rdampij, rdampik, rdampjk, c9, d9, triple, ang, dE
   real(wp) :: cutoff2, cutoff_inner, swij, swjk, swik, dswdr, sw

   cutoff2 = cutoff*cutoff
   call get_disp3_switch(cutoff, cutoff_inner, use_switch)
   
   !$omp parallel do schedule(dynamic) default(none) reduction(+:energy) &
   !$omp shared(mol, self, cache, damp, param, trans, cutoff, cutoff2, cutoff_inner, use_switch) &
   !$omp private(iat, jat, kat, izp, jzp, kzp, jtr, ktr, vij, vjk, vik, &
   !$omp& r1, r2, r3, r5, r2ij, r2jk, r2ik, rij, rjk, rik, rdamp, rdampij, &
   !$omp& rdampik, rdampjk, d9, ang, triple, c9, dE, swij, swjk, swik, dswdr, sw)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)
         rdampij = self%get_3b_rdamp(izp, jzp)
         do jtr = 1, size(trans, 2)
            vij(:) = mol%xyz(:, jat) + trans(:, jtr) - mol%xyz(:, iat)
            r2ij = vij(1)*vij(1) + vij(2)*vij(2) + vij(3)*vij(3)
            if (r2ij > cutoff2 .or. r2ij < epsilon(1.0_wp)) cycle
            if (use_switch) then
               rij = sqrt(r2ij)
               call smooth_cutoff(rij, cutoff, cutoff_inner, use_switch, swij, dswdr)
            else
               swij = 1.0_wp
            end if
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
                  if (use_switch) then
                     rik = sqrt(r2ik)
                     call smooth_cutoff(rik, cutoff, cutoff_inner, use_switch, swik, dswdr)
                  else
                     swik = 1.0_wp
                  end if
                  vjk(:) = mol%xyz(:, kat) + trans(:, ktr) - mol%xyz(:, jat) &
                     & - trans(:, jtr)
                  r2jk = vjk(1)*vjk(1) + vjk(2)*vjk(2) + vjk(3)*vjk(3)
                  if (r2jk > cutoff2 .or. r2jk < epsilon(1.0_wp)) cycle
                  if (use_switch) then
                     rjk = sqrt(r2jk)
                     call smooth_cutoff(rjk, cutoff, cutoff_inner, use_switch, swjk, dswdr)
                  else
                     swjk = 1.0_wp
                  end if
                  sw = swij * swik * swjk
                  if (sw <= 0.0_wp) cycle
                  r2 = r2ij*r2ik*r2jk
                  r1 = sqrt(r2)
                  r3 = r2 * r1
                  r5 = r3 * r2

                  call damp%get_3b_damp(param, r1, r2ij, r2ik, r2jk, &
                     & rdamp, rdampij, rdampik, rdampjk, d9)

                  ang = 0.375_wp*(r2ij + r2jk - r2ik)*(r2ij - r2jk + r2ik)&
                     & *(-r2ij + r2jk + r2ik) / r5 + 1.0_wp / r3

                  dE = c9 * triple * d9 * ang / 3.0_wp * sw
                  energy(iat) = energy(iat) - dE
                  energy(jat) = energy(jat) - dE
                  energy(kat) = energy(kat) - dE
               end do
            end do
         end do
      end do
   end do

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
   logical :: use_switch
   real(wp) :: vij(3), vjk(3), vik(3), r2ij, r2jk, r2ik, rij, rjk, rik
   real(wp) :: r1, r2, r3, r5
   real(wp) :: rdamp, rdampij, rdampik, rdampjk
   real(wp) :: c9, dc9dcni, dc9dqi, dc9dcnj, dc9dqj, dc9dcnk, dc9dqk
   real(wp) :: d9, d9drij, d9drik, d9drjk, ang, dang, triple
   real(wp) :: cutoff2, cutoff_inner
   real(wp) :: dE, dGij(3), dGik(3), dGjk(3), dS(3,3)
   real(wp) :: dE0, swij, swjk, swik, dswijdr, dswjkdr, dswikdr, sw

   cutoff2 = cutoff*cutoff
   call get_disp3_switch(cutoff, cutoff_inner, use_switch)

   !$omp parallel do schedule(dynamic) default(none) &
   !$omp reduction(+:energy, gradient, sigma, dEdcn, dEdq) &
   !$omp shared(mol, self, cache, damp, param, trans, cutoff, cutoff2, cutoff_inner, use_switch) &
   !$omp private(iat, jat, kat, izp, jzp, kzp, jtr, ktr, vij, vjk, vik, &
   !$omp& r1, r2, r3, r5, r2ij, r2jk, r2ik, rij, rjk, rik, rdamp, rdampij, rdampik, rdampjk, &
   !$omp& d9, d9drij, d9drik, d9drjk, ang, dang, triple, &
   !$omp& c9, dc9dcni, dc9dqi, dc9dcnj, dc9dqj, dc9dcnk, dc9dqk, &
   !$omp& dE, dE0, dGij, dGik, dGjk, dS, swij, swjk, swik, dswijdr, dswjkdr, &
   !$omp& dswikdr, sw)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)
         rdampij = self%get_3b_rdamp(izp, jzp)
         do jtr = 1, size(trans, 2)
            vij(:) = mol%xyz(:, jat) + trans(:, jtr) - mol%xyz(:, iat)
            r2ij = vij(1)*vij(1) + vij(2)*vij(2) + vij(3)*vij(3)
            if (r2ij > cutoff2 .or. r2ij < epsilon(1.0_wp)) cycle
            if (use_switch) then
               rij = sqrt(r2ij)
               call smooth_cutoff(rij, cutoff, cutoff_inner, use_switch, swij, dswijdr)
            else
               rij = 1.0_wp
               swij = 1.0_wp
               dswijdr = 0.0_wp
            end if
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
                  if (use_switch) then
                     rik = sqrt(r2ik)
                     call smooth_cutoff(rik, cutoff, cutoff_inner, use_switch, swik, dswikdr)
                  else
                     rik = 1.0_wp
                     swik = 1.0_wp
                     dswikdr = 0.0_wp
                  end if
                  vjk(:) = mol%xyz(:, kat) + trans(:, ktr) - mol%xyz(:, jat) &
                     & - trans(:, jtr)
                  r2jk = vjk(1)*vjk(1) + vjk(2)*vjk(2) + vjk(3)*vjk(3)
                  if (r2jk > cutoff2 .or. r2jk < epsilon(1.0_wp)) cycle
                  if (use_switch) then
                     rjk = sqrt(r2jk)
                     call smooth_cutoff(rjk, cutoff, cutoff_inner, use_switch, swjk, dswjkdr)
                  else
                     rjk = 1.0_wp
                     swjk = 1.0_wp
                     dswjkdr = 0.0_wp
                  end if
                  sw = swij * swik * swjk
                  if (sw <= 0.0_wp) cycle
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
                  dE0 = c9 * d9 * ang
                  dGij(:) = sw * c9 * (-dang * d9 / r2ij - ang * d9drij) * vij &
                     & - dE0 * dswijdr / rij * swik * swjk * vij

                  ! d/drik
                  dang = -0.375_wp * (r2ik**3 + r2ik**2 * (r2jk + r2ij) &
                     & + r2ik * (3.0_wp * r2jk**2 + 2.0_wp * r2jk * r2ij &
                     & + 3.0_wp * r2ij**2) &
                     & - 5.0_wp * (r2jk - r2ij)**2 * (r2jk + r2ij)) / r5
                  dGik(:) = sw * c9 * (-dang * d9 / r2ik - ang * d9drik) * vik &
                     & - dE0 * dswikdr / rik * swij * swjk * vik

                  ! d/drjk
                  dang = -0.375_wp * (r2jk**3 + r2jk**2*(r2ik + r2ij) &
                     & + r2jk * (3.0_wp * r2ik**2 + 2.0_wp * r2ik * r2ij &
                     & + 3.0_wp * r2ij**2) &
                     & - 5.0_wp * (r2ik - r2ij)**2 * (r2ik + r2ij)) / r5
                  dGjk(:) = sw * c9 * (-dang * d9 / r2jk - ang * d9drjk) * vjk &
                     & - dE0 * dswjkdr / rjk * swij * swik * vjk

                  dE = dE0 * triple * sw / 3.0_wp
                  energy(iat) = energy(iat) - dE
                  energy(jat) = energy(jat) - dE
                  energy(kat) = energy(kat) - dE

                  gradient(:, iat) = gradient(:, iat) - (dGij + dGik) * triple
                  gradient(:, jat) = gradient(:, jat) + (dGij - dGjk) * triple
                  gradient(:, kat) = gradient(:, kat) + (dGik + dGjk) * triple
                  
                  dS(:, :) = spread(dGij, 1, 3) * spread(vij, 2, 3)&
                        & + spread(dGik, 1, 3) * spread(vik, 2, 3)&
                        & + spread(dGjk, 1, 3) * spread(vjk, 2, 3)
                  sigma(:, :) = sigma(:, :) + dS * triple

                  dE = d9 * ang * triple * sw
                  dEdcn(iat) = dEdcn(iat) - dE * dc9dcni
                  dEdcn(jat) = dEdcn(jat) - dE * dc9dcnj
                  dEdcn(kat) = dEdcn(kat) - dE * dc9dcnk

                  dEdq(iat) = dEdq(iat) - dE * dc9dqi
                  dEdq(jat) = dEdq(jat) - dE * dc9dqj
                  dEdq(kat) = dEdq(kat) - dE * dc9dqk
               end do
            end do
         end do
      end do
   end do

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
   logical :: use_switch
   real(wp) :: vec(3), r2, r, cutoff2, cutoff_inner, c6, c8, rdamp, d6, d8
   real(wp) :: dE, sw, dswdr
   
   if (.not. allocated(param%s6) .or. .not. allocated(param%s8)) return
   if ((abs(param%s6) < epsilon(1.0_wp) .and. abs(param%s8) < epsilon(1.0_wp))) return

   cutoff2 = cutoff*cutoff
   call get_disp2_switch(cutoff, cutoff_inner, use_switch)
   !$omp parallel do schedule(runtime) default(none) reduction(+:energy) &
   !$omp shared(mol, self, cache, damp, param, trans, cutoff, cutoff2, cutoff_inner, use_switch) &
   !$omp private(iat, jat, jtr, izp, jzp, vec, r2, r, c6, c8, rdamp, d6, d8, dE, sw, dswdr)
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
            if (use_switch) then
               r = sqrt(r2)
               call smooth_cutoff(r, cutoff, cutoff_inner, use_switch, sw, dswdr)
               if (sw <= 0.0_wp) cycle
            else
               sw = 1.0_wp
            end if
            
            call damp%get_2b_damp(param, r2, rdamp, d6, d8)

            dE = -0.5_wp * sw * (c6 * d6 + c8 * d8)

            energy(jat, iat) = energy(jat, iat) + dE
            if (iat /= jat) then
               energy(iat, jat) = energy(iat, jat) + dE
            end if
         end do
      end do
   end do

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
   logical :: use_switch
   real(wp) :: vij(3), vjk(3), vik(3), r2ij, r2jk, r2ik, rij, rjk, rik
   real(wp) :: r1, r2, r3, r5
   real(wp) :: rdamp, rdampij, rdampik, rdampjk, c9, d9, triple, ang, dE
   real(wp) :: cutoff2, cutoff_inner, swij, swjk, swik, dswdr, sw

   if (.not.allocated(param%s9)) return
   if (abs(param%s9) < epsilon(1.0_wp)) return

   cutoff2 = cutoff*cutoff
   call get_disp3_switch(cutoff, cutoff_inner, use_switch)
   !$omp parallel do schedule(dynamic) default(none) reduction(+:energy) &
   !$omp shared(mol, self, cache, damp, param, trans, cutoff, cutoff2, cutoff_inner, use_switch) &
   !$omp private(iat, jat, kat, izp, jzp, kzp, jtr, ktr, vij, vjk, vik, &
   !$omp& r1, r2, r3, r5, r2ij, r2jk, r2ik, rij, rjk, rik, rdamp, rdampij, &
   !$omp& rdampik, rdampjk, d9, ang, triple, c9, dE, swij, swjk, swik, dswdr, sw)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)
         rdampij = self%get_3b_rdamp(izp, jzp)
         do jtr = 1, size(trans, 2)
            vij(:) = mol%xyz(:, jat) + trans(:, jtr) - mol%xyz(:, iat)
            r2ij = vij(1)*vij(1) + vij(2)*vij(2) + vij(3)*vij(3)
            if (r2ij > cutoff2 .or. r2ij < epsilon(1.0_wp)) cycle
            if (use_switch) then
               rij = sqrt(r2ij)
               call smooth_cutoff(rij, cutoff, cutoff_inner, use_switch, swij, dswdr)
            else
               swij = 1.0_wp
            end if
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
                  if (use_switch) then
                     rik = sqrt(r2ik)
                     call smooth_cutoff(rik, cutoff, cutoff_inner, use_switch, swik, dswdr)
                  else
                     swik = 1.0_wp
                  end if
                  vjk(:) = mol%xyz(:, kat) + trans(:, ktr) - mol%xyz(:, jat) &
                     & - trans(:, jtr)
                  r2jk = vjk(1)*vjk(1) + vjk(2)*vjk(2) + vjk(3)*vjk(3)
                  if (r2jk > cutoff2 .or. r2jk < epsilon(1.0_wp)) cycle
                  if (use_switch) then
                     rjk = sqrt(r2jk)
                     call smooth_cutoff(rjk, cutoff, cutoff_inner, use_switch, swjk, dswdr)
                  else
                     swjk = 1.0_wp
                  end if
                  sw = swij * swik * swjk
                  if (sw <= 0.0_wp) cycle
                  r2 = r2ij*r2ik*r2jk
                  r1 = sqrt(r2)
                  r3 = r2 * r1
                  r5 = r3 * r2

                  call damp%get_3b_damp(param, r1, r2ij, r2ik, r2jk, &
                     & rdamp, rdampij, rdampik, rdampjk, d9)

                  ang = 0.375_wp*(r2ij + r2jk - r2ik)*(r2ij - r2jk + r2ik)&
                     & *(-r2ij + r2jk + r2ik) / r5 + 1.0_wp / r3

                  dE = c9 * triple * d9 * ang / 6.0_wp * sw
                  energy(jat, iat) = energy(jat, iat) - dE
                  energy(kat, iat) = energy(kat, iat) - dE
                  energy(iat, jat) = energy(iat, jat) - dE
                  energy(kat, jat) = energy(kat, jat) - dE
                  energy(iat, kat) = energy(iat, kat) - dE
                  energy(jat, kat) = energy(jat, kat) - dE
               end do
            end do
         end do
      end do
   end do

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
