! This file is part of s-dftd3.
! SPDX-Identifier: LGPL-3.0-or-later
!
! s-dftd3 is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! s-dftd3 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with s-dftd3.  If not, see <https://www.gnu.org/licenses/>.

module dftd3_damping_rational
   use dftd3_damping, only : damping_param
   use dftd3_damping_atm, only : get_atm_dispersion, get_atm_pairwise_dispersion
   use dftd3_param, only : d3_param
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   implicit none

   public :: rational_damping_param, new_rational_damping


   !> Rational (Becke-Johnson) damping model
   type, extends(damping_param) :: rational_damping_param
      real(wp) :: s6
      real(wp) :: s8
      real(wp) :: s9
      real(wp) :: a1
      real(wp) :: a2
      real(wp) :: alp
   contains

      !> Evaluate pairwise dispersion energy expression
      procedure :: get_dispersion2

      !> Evaluate ATM three-body dispersion energy expression
      procedure :: get_dispersion3

      !> Evaluate pairwise representation of additive dispersion energy
      procedure :: get_pairwise_dispersion2

      !> Evaluate pairwise representation of non-additive dispersion energy
      procedure :: get_pairwise_dispersion3

   end type rational_damping_param


   real(wp), parameter :: rs9 = 4.0_wp/3.0_wp


contains


subroutine get_disp2_switch(cutoff, inner, active)

   real(wp), intent(in) :: cutoff
   real(wp), intent(out) :: inner
   logical, intent(out) :: active

   character(len=64) :: env
   integer :: stat, io
   real(wp) :: width

   inner = cutoff
   width = 0.05_wp

   call get_environment_variable("SDFTD3_DISP2_SMOOTH_WIDTH", env, status=stat)
   if (stat /= 0 .or. len_trim(env) == 0) then
      call get_environment_variable("DFTD3_DISP2_SMOOTH_WIDTH", env, status=stat)
   end if
   if (stat /= 0 .or. len_trim(env) == 0) then
      call get_environment_variable("TBLITE_D3_DISP2_SMOOTH_WIDTH", env, status=stat)
   end if
   if (stat == 0 .and. len_trim(env) > 0) then
      read(env, *, iostat=io) width
   end if
   active = width > 0.0_wp .and. width < cutoff
   if (active) inner = cutoff - width

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

   call get_environment_variable("SDFTD3_DISP3_SMOOTH_WIDTH", env, status=stat)
   if (stat /= 0 .or. len_trim(env) == 0) then
      call get_environment_variable("DFTD3_DISP3_SMOOTH_WIDTH", env, status=stat)
   end if
   if (stat /= 0 .or. len_trim(env) == 0) then
      call get_environment_variable("TBLITE_D3_DISP3_SMOOTH_WIDTH", env, status=stat)
   end if
   if (stat == 0 .and. len_trim(env) > 0) then
      read(env, *, iostat=io) width
      if (io == 0 .and. width > 0.0_wp .and. width < cutoff) then
         inner = cutoff - width
         active = .true.
      end if
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


!> Create new rational damping model
subroutine new_rational_damping(self, param)

   !> Rational damping parameters
   type(rational_damping_param), intent(out) :: self

   !> Parameters
   type(d3_param), intent(in) :: param

   self%s6 = param%s6
   self%s8 = param%s8
   self%s9 = param%s9
   self%a1 = param%a1
   self%a2 = param%a2
   self%alp = param%alp

end subroutine new_rational_damping


!> Evaluation of the dispersion energy expression
subroutine get_dispersion2(self, mol, trans, cutoff, rvdw, r4r2, c6, dc6dcn, &
      & energy, dEdcn, gradient, sigma)

   !> Damping parameters
   class(rational_damping_param), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Van-der-Waals radii for damping function
   real(wp), intent(in) :: rvdw(:, :)

   !> Expectation values for C8 extrapolation
   real(wp), intent(in) :: r4r2(:)

   !> C6 coefficients for all atom pairs.
   real(wp), intent(in) :: c6(:, :)

   !> Derivative of the C6 w.r.t. the coordination number
   real(wp), intent(in), optional :: dc6dcn(:, :)

   !> Dispersion energy
   real(wp), intent(inout) :: energy(:)

   !> Derivative of the energy w.r.t. the coordination number
   real(wp), intent(inout), optional :: dEdcn(:)

   !> Dispersion gradient
   real(wp), intent(inout), optional :: gradient(:, :)

   !> Dispersion virial
   real(wp), intent(inout), optional :: sigma(:, :)

   logical :: grad

   if (abs(self%s6) < epsilon(1.0_wp) .and. abs(self%s8) < epsilon(1.0_wp)) return
   grad = present(dc6dcn) .and. present(dEdcn) .and. present(gradient) &
      & .and. present(sigma)

   if (grad) then
      call get_dispersion_derivs(self, mol, trans, cutoff, rvdw, r4r2, c6, dc6dcn, &
         & energy, dEdcn, gradient, sigma)
   else
      call get_dispersion_energy(self, mol, trans, cutoff, rvdw, r4r2, c6, energy)
   end if

end subroutine get_dispersion2


!> Evaluation of the dispersion energy expression
subroutine get_dispersion_energy(self, mol, trans, cutoff, rvdw, r4r2, c6, energy)

   !> Damping parameters
   class(rational_damping_param), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Van-der-Waals radii for damping function
   real(wp), intent(in) :: rvdw(:, :)

   !> Expectation values for C8 extrapolation
   real(wp), intent(in) :: r4r2(:)

   !> C6 coefficients for all atom pairs.
   real(wp), intent(in) :: c6(:, :)

   !> Dispersion energy
   real(wp), intent(inout) :: energy(:)

   integer :: iat, jat, izp, jzp, jtr
   logical :: use_switch
   real(wp) :: vec(3), r2, r, cutoff2, cutoff_inner, r0ij, rrij, c6ij, t6, t8
   real(wp) :: edisp, dE, sw, dswdr
   real(wp) :: r0ij2, r0ij6, r0ij8, r4

   cutoff2 = cutoff*cutoff
   call get_disp2_switch(cutoff, cutoff_inner, use_switch)

   !$omp parallel do schedule(runtime) default(none) reduction(+:energy) &
   !$omp shared(mol, self, c6, trans, cutoff, cutoff2, cutoff_inner, use_switch, r4r2) &
   !$omp private(iat, jat, izp, jzp, jtr, vec, r2, r0ij, rrij, c6ij, t6, &
   !$omp& t8, edisp, dE, r0ij2, r0ij6, r0ij8, r4, r, sw, dswdr)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)
         rrij = 3*r4r2(izp)*r4r2(jzp)
         r0ij = self%a1 * sqrt(rrij) + self%a2
         c6ij = c6(jat, iat)
         r0ij2 = r0ij * r0ij
         r0ij6 = r0ij2 * r0ij2 * r0ij2
         r0ij8 = r0ij6 * r0ij2
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

            r4 = r2 * r2
            t6 = 1.0_wp/(r4 * r2 + r0ij6)
            t8 = 1.0_wp/(r4 * r4 + r0ij8)

            edisp = sw * (self%s6*t6 + self%s8*rrij*t8)

            dE = -c6ij*edisp * 0.5_wp

            energy(iat) = energy(iat) + dE
            if (iat /= jat) then
               energy(jat) = energy(jat) + dE
            end if
         end do
      end do
   end do

end subroutine get_dispersion_energy


!> Evaluation of the dispersion energy expression
subroutine get_dispersion_derivs(self, mol, trans, cutoff, rvdw, r4r2, c6, dc6dcn, &
      & energy, dEdcn, gradient, sigma)

   !> Damping parameters
   class(rational_damping_param), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Van-der-Waals radii for damping function
   real(wp), intent(in) :: rvdw(:, :)

   !> Expectation values for C8 extrapolation
   real(wp), intent(in) :: r4r2(:)

   !> C6 coefficients for all atom pairs.
   real(wp), intent(in) :: c6(:, :)

   !> Derivative of the C6 w.r.t. the coordination number
   real(wp), intent(in) :: dc6dcn(:, :)

   !> Dispersion energy
   real(wp), intent(inout) :: energy(:)

   !> Derivative of the energy w.r.t. the coordination number
   real(wp), intent(inout) :: dEdcn(:)

   !> Dispersion gradient
   real(wp), intent(inout) :: gradient(:, :)

   !> Dispersion virial
   real(wp), intent(inout) :: sigma(:, :)

   integer :: iat, jat, izp, jzp, jtr, ic, jc
   logical :: use_switch
   real(wp) :: vec(3), r2, r, cutoff2, cutoff_inner, r0ij, rrij, c6ij, t6, t8
   real(wp) :: d6, d8, edisp0, gdisp0, edisp, gdisp, sw, dswdr
   real(wp) :: dE, dG(3), dS(3, 3)
   real(wp) :: r0ij2, r0ij6, r0ij8, r4

   cutoff2 = cutoff*cutoff
   call get_disp2_switch(cutoff, cutoff_inner, use_switch)

   !$omp parallel do schedule(runtime) default(none) reduction(+:energy, gradient, sigma, dEdcn) &
   !$omp shared(mol, self, c6, dc6dcn, trans, cutoff, cutoff2, cutoff_inner, use_switch, r4r2) &
   !$omp private(iat, jat, izp, jzp, jtr, ic, jc, vec, r2, r0ij, rrij, c6ij, t6, t8, &
   !$omp& d6, d8, edisp0, gdisp0, edisp, gdisp, dE, dG, dS, r0ij2, r0ij6, &
   !$omp& r0ij8, r4, r, sw, dswdr)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)
         rrij = 3*r4r2(izp)*r4r2(jzp)
         r0ij = self%a1 * sqrt(rrij) + self%a2
         c6ij = c6(jat, iat)
         r0ij2 = r0ij * r0ij
         r0ij6 = r0ij2 * r0ij2 * r0ij2
         r0ij8 = r0ij6 * r0ij2
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

            r4 = r2 * r2
            t6 = 1.0_wp/(r4 * r2 + r0ij6)
            t8 = 1.0_wp/(r4 * r4 + r0ij8)

            d6 = -6*r4*t6*t6
            d8 = -8*r4*r2*t8*t8

            edisp0 = self%s6*t6 + self%s8*rrij*t8
            gdisp0 = self%s6*d6 + self%s8*rrij*d8
            edisp = sw * edisp0
            gdisp = sw * gdisp0 + dswdr * edisp0 / r

            dE = -c6ij*edisp * 0.5_wp
            dG(:) = -c6ij*gdisp*vec
            do ic = 1, 3
               do jc = 1, 3
                  dS(ic, jc) = dG(ic) * vec(jc) * 0.5_wp
               end do
            end do

            energy(iat) = energy(iat) + dE
            dEdcn(iat) = dEdcn(iat) - dc6dcn(iat, jat) * edisp
            sigma(:, :) = sigma(:, :) + dS
            if (iat /= jat) then
               energy(jat) = energy(jat) + dE
               dEdcn(jat) = dEdcn(jat) - dc6dcn(jat, iat) * edisp
               gradient(:, iat) = gradient(:, iat) + dG
               gradient(:, jat) = gradient(:, jat) - dG
               sigma(:, :) = sigma(:, :) + dS
            end if
         end do
      end do
   end do

end subroutine get_dispersion_derivs


!> Evaluation of the dispersion energy expression
subroutine get_dispersion3(self, mol, trans, cutoff, rvdw, r4r2, c6, dc6dcn, &
      & energy, dEdcn, gradient, sigma)

   !> Damping parameters
   class(rational_damping_param), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Van-der-Waals radii for damping function
   real(wp), intent(in) :: rvdw(:, :)

   !> Expectation values for C8 extrapolation
   real(wp), intent(in) :: r4r2(:)

   !> C6 coefficients for all atom pairs.
   real(wp), intent(in) :: c6(:, :)

   !> Derivative of the C6 w.r.t. the coordination number
   real(wp), intent(in), optional :: dc6dcn(:, :)

   !> Dispersion energy
   real(wp), intent(inout) :: energy(:)

   !> Derivative of the energy w.r.t. the coordination number
   real(wp), intent(inout), optional :: dEdcn(:)

   !> Dispersion gradient
   real(wp), intent(inout), optional :: gradient(:, :)

   !> Dispersion virial
   real(wp), intent(inout), optional :: sigma(:, :)

   call get_atm_dispersion(mol, trans, cutoff, self%s9, rs9, self%alp+2, &
      & rvdw, c6, dc6dcn, energy, dEdcn, gradient, sigma)

end subroutine get_dispersion3


!> Evaluation of the dispersion energy expression projected on atomic pairs
subroutine get_pairwise_dispersion2(self, mol, trans, cutoff, rvdw, r4r2, c6, energy)

   !> Damping parameters
   class(rational_damping_param), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Van-der-Waals radii for damping function
   real(wp), intent(in) :: rvdw(:, :)

   !> Expectation values for r4 over r2 operator
   real(wp), intent(in) :: r4r2(:)

   !> C6 coefficients for all atom pairs.
   real(wp), intent(in) :: c6(:, :)

   !> Dispersion energy
   real(wp), intent(inout) :: energy(:, :)

   integer :: iat, jat, izp, jzp, jtr
   logical :: use_switch
   real(wp) :: vec(3), r2, r, cutoff2, cutoff_inner, r0ij, rrij, c6ij, t6, t8
   real(wp) :: edisp, dE, sw, dswdr

   cutoff2 = cutoff*cutoff
   call get_disp2_switch(cutoff, cutoff_inner, use_switch)

   !$omp parallel do schedule(runtime) default(none) reduction(+:energy) &
   !$omp shared(mol, self, c6, trans, cutoff, cutoff2, cutoff_inner, use_switch, r4r2) &
   !$omp private(iat, jat, izp, jzp, jtr, vec, r2, r0ij, rrij, c6ij, t6, &
   !$omp& t8, edisp, dE, r, sw, dswdr)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)
         rrij = 3*r4r2(izp)*r4r2(jzp)
         r0ij = self%a1 * sqrt(rrij) + self%a2
         c6ij = c6(jat, iat)
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

            t6 = 1.0_wp/(r2**3 + r0ij**6)
            t8 = 1.0_wp/(r2**4 + r0ij**8)

            edisp = sw * (self%s6*t6 + self%s8*rrij*t8)

            dE = -c6ij*edisp * 0.5_wp

            energy(jat, iat) = energy(jat, iat) + dE
            if (iat /= jat) then
               energy(iat, jat) = energy(iat, jat) + dE
            end if
         end do
      end do
   end do

end subroutine get_pairwise_dispersion2


!> Evaluation of the dispersion energy expression
subroutine get_pairwise_dispersion3(self, mol, trans, cutoff, rvdw, r4r2, c6, energy)

   !> Damping parameters
   class(rational_damping_param), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Van-der-Waals radii for damping function
   real(wp), intent(in) :: rvdw(:, :)

   !> Expectation values for r4 over r2 operator
   real(wp), intent(in) :: r4r2(:)

   !> C6 coefficients for all atom pairs.
   real(wp), intent(in) :: c6(:, :)

   !> Dispersion energy
   real(wp), intent(inout) :: energy(:, :)

   call get_atm_pairwise_dispersion(mol, trans, cutoff, self%s9, rs9, self%alp+2, &
      & rvdw, c6, energy)

end subroutine get_pairwise_dispersion3


end module dftd3_damping_rational
