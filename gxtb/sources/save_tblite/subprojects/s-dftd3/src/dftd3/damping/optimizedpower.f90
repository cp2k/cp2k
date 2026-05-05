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

module dftd3_damping_optimizedpower
   use dftd3_damping, only : damping_param
   use dftd3_damping_atm, only : get_atm_dispersion, get_atm_pairwise_dispersion
   use dftd3_param, only : d3_param
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   implicit none

   public :: optimizedpower_damping_param, new_optimizedpower_damping


   !> Optimized power damping model
   type, extends(damping_param) :: optimizedpower_damping_param
      real(wp) :: s6
      real(wp) :: s8
      real(wp) :: s9
      real(wp) :: a1
      real(wp) :: a2
      real(wp) :: alp
      real(wp) :: bet
   contains

      !> Evaluate pairwise dispersion energy expression
      procedure :: get_dispersion2

      !> Evaluate ATM three-body dispersion energy expression
      procedure :: get_dispersion3

      !> Evaluate pairwise representation of additive dispersion energy
      procedure :: get_pairwise_dispersion2

      !> Evaluate pairwise representation of non-additive dispersion energy
      procedure :: get_pairwise_dispersion3

   end type optimizedpower_damping_param


   real(wp), parameter :: rs9 = 4.0_wp/3.0_wp


contains


!> Create new optimized power damping model
subroutine new_optimizedpower_damping(self, param)

   !> Optimized power damping parameters
   type(optimizedpower_damping_param), intent(out) :: self

   !> Parameters
   type(d3_param), intent(in) :: param

   self%s6 = param%s6
   self%s8 = param%s8
   self%s9 = param%s9
   self%a1 = param%a1
   self%a2 = param%a2
   self%alp = param%alp
   self%bet = param%bet

end subroutine new_optimizedpower_damping


!> Evaluation of the dispersion energy expression
subroutine get_dispersion2(self, mol, trans, cutoff, rvdw, r4r2, c6, dc6dcn, &
      & energy, dEdcn, gradient, sigma)

   !> Damping parameters
   class(optimizedpower_damping_param), intent(in) :: self

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
   class(optimizedpower_damping_param), intent(in) :: self

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
   real(wp) :: vec(3), r2, cutoff2, r0ij, rrij, c6ij, t6, t8, edisp, dE, rb, ab
   real(wp) :: r0ij2, r0ij6, r0ij8, abr0ij6, abr0ij8, r4


   cutoff2 = cutoff*cutoff

   !$omp parallel do schedule(runtime) default(none) reduction(+:energy) &
   !$omp shared(mol, self, c6, trans, cutoff2, r4r2) &
   !$omp private(iat, jat, izp, jzp, jtr, vec, r2, r0ij, rrij, c6ij, t6, &
   !$omp& t8, edisp, dE, rb, ab, r0ij2, r0ij6, r0ij8, abr0ij6, abr0ij8, r4)
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
         ab = r0ij**self%bet
         abr0ij6 = ab * r0ij6
         abr0ij8 = ab * r0ij8
         do jtr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, jtr))
            r2 = vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3)
            if (r2 > cutoff2 .or. r2 < epsilon(1.0_wp)) cycle

            rb = r2**(self%bet*0.5_wp)
            r4 = r2 * r2
            t6 = rb/(rb*r4*r2 + abr0ij6)
            t8 = rb/(rb*r4*r4 + abr0ij8)

            edisp = self%s6*t6 + self%s8*rrij*t8

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
   class(optimizedpower_damping_param), intent(in) :: self

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
   real(wp) :: vec(3), r2, cutoff2, r0ij, rrij, c6ij, t6, t8, d6, d8, edisp, gdisp
   real(wp) :: dE, dG(3), dS(3, 3), rb, ab
   real(wp) :: r0ij2, r0ij6, r0ij8, abr0ij6, abr0ij8, r4


   cutoff2 = cutoff*cutoff

   !$omp parallel do schedule(runtime) default(none) reduction(+:energy, gradient, sigma, dEdcn) &
   !$omp shared(mol, self, c6, dc6dcn, trans, cutoff2, r4r2) &
   !$omp private(iat, jat, izp, jzp, jtr, ic, jc, vec, r2, r0ij, rrij, c6ij, t6, t8, &
   !$omp& d6, d8, edisp, gdisp, dE, dG, dS, rb, ab, &
   !$omp& r0ij2, r0ij6, r0ij8, abr0ij6, abr0ij8, r4)
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
         ab = r0ij**self%bet
         abr0ij6 = ab * r0ij6
         abr0ij8 = ab * r0ij8
         do jtr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, jtr))
            r2 = vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3)
            if (r2 > cutoff2 .or. r2 < epsilon(1.0_wp)) cycle

            rb = r2**(self%bet*0.5_wp)
            r4 = r2 * r2
            t6 = rb/(rb*r4*r2 + abr0ij6)
            t8 = rb/(rb*r4*r4 + abr0ij8)

            d6 = -6*r4*t6*t6 + self%bet*abr0ij6 * t6*t6 / (rb*r2)
            d8 = -8*r4*r2*t8*t8 + self%bet*abr0ij8 * t8*t8 / (rb*r2)

            edisp = self%s6*t6 + self%s8*rrij*t8
            gdisp = self%s6*d6 + self%s8*rrij*d8

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
   class(optimizedpower_damping_param), intent(in) :: self

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
   class(optimizedpower_damping_param), intent(in) :: self

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
   real(wp) :: vec(3), r2, cutoff2, r0ij, rrij, c6ij, t6, t8, edisp, dE, rb, ab


   cutoff2 = cutoff*cutoff

   !$omp parallel do schedule(runtime) default(none) reduction(+:energy) &
   !$omp shared(mol, self, c6, trans, cutoff2, r4r2) &
   !$omp private(iat, jat, izp, jzp, jtr, vec, r2, r0ij, rrij, c6ij, t6, &
   !$omp& t8, edisp, dE, rb, ab)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)
         rrij = 3*r4r2(izp)*r4r2(jzp)
         r0ij = self%a1 * sqrt(rrij) + self%a2
         c6ij = c6(jat, iat)
         ab = r0ij**self%bet
         do jtr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, jtr))
            r2 = vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3)
            if (r2 > cutoff2 .or. r2 < epsilon(1.0_wp)) cycle

            rb = r2**(self%bet*0.5_wp)
            t6 = rb/(rb*r2**3 + ab*r0ij**6)
            t8 = rb/(rb*r2**4 + ab*r0ij**8)

            edisp = self%s6*t6 + self%s8*rrij*t8

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
   class(optimizedpower_damping_param), intent(in) :: self

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


end module dftd3_damping_optimizedpower
