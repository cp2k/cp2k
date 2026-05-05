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

!> Implementation of the rational (Becke--Johnson) damping function.
module dftd4_damping_rational
   use dftd4_damping_type, only : damping_twobody, damping_threebody
   use dftd4_param_type, only : param_type
   use mctc_env, only : error_type, fatal_error, wp
   implicit none
   private

   public :: rational_damping_twobody, new_rational_damping_twobody
   public :: rational_damping_threebody, new_rational_damping_threebody

   !> Rational (Becke-Johnson) damping for two-body terms
   type, extends(damping_twobody) :: rational_damping_twobody
   contains
      !> Evaluate two-body damping factor
      procedure :: get_2b_damp
      !> Evaluate two-body damping factor with derivatives
      procedure :: get_2b_derivs
      !> Check the availability of the parameters for two-body damping
      procedure :: check_2b_params
   end type rational_damping_twobody

   !> Rational (Becke-Johnson) damping for three-body terms
   type, extends(damping_threebody) :: rational_damping_threebody
   contains
      !> Evaluate three-body damping factor
      procedure :: get_3b_damp
      !> Evaluate three-body damping factor with derivatives
      procedure :: get_3b_derivs
      !> Check the availability of the parameters for three-body damping
      procedure :: check_3b_params
   end type rational_damping_threebody

   character(len=*), parameter :: rational_label = "Rational (Becke-Johnson) damping"
   character(len=*), parameter :: rational_label_short = "BJ"

contains


!> Create a new instance of the rational damping for two-body terms
subroutine new_rational_damping_twobody(self)
   !> Instance of the two-body rational damping function
   type(rational_damping_twobody), intent(out) :: self

   self%label = rational_label
   self%label_short = rational_label_short
end subroutine new_rational_damping_twobody


!> Create a new instance of the rational damping for three-body terms
subroutine new_rational_damping_threebody(self)
   !> Instance of the three-body rational damping function
   type(rational_damping_threebody), intent(out) :: self

   self%label = rational_label
   self%label_short = rational_label_short
end subroutine new_rational_damping_threebody


!> Evaluate two-body damping factor
pure subroutine get_2b_damp(self, param, r2, rdamp, d6, d8)
   !> Instance of the two-body damping function
   class(rational_damping_twobody), intent(in) :: self
   !> Damping parameters
   type(param_type), intent(in) :: param
   !> Square of interatomic distance
   real(wp), intent(in) :: r2
   !> Damping radius
   real(wp), intent(in) :: rdamp
   !> Damping factor for C6/R^6 term
   real(wp), intent(out) :: d6
   !> Damping factor for C8/R^8 term
   real(wp), intent(out) :: d8

   real(wp) :: r6, r8, rpoly, rpoly2, rpoly6, rpoly8

   r6 = r2 * r2 * r2
   r8 = r6 * r2

   ! First-order polynomial scaling of the damping radius
   rpoly = param%a1 * rdamp + param%a2

   rpoly2 = rpoly * rpoly
   rpoly6 = rpoly2 * rpoly2 * rpoly2
   rpoly8 = rpoly6 * rpoly2

   d6 = param%s6 / (r6 + rpoly6)
   d8 = param%s8 / (r8 + rpoly8)


end subroutine get_2b_damp


!> Evaluate two-body damping factor with derivarives
pure subroutine get_2b_derivs(self, param, r2, rdamp, d6, d8, d6dr, d8dr)
   !> Instance of the two-body damping function
   class(rational_damping_twobody), intent(in) :: self
   !> Damping parameters
   type(param_type), intent(in) :: param
   !> Square of interatomic distance
   real(wp), intent(in) :: r2
   !> Damping radius
   real(wp), intent(in) :: rdamp
   !> Damping factor for C6/R^6 term
   real(wp), intent(out) :: d6
   !> Damping factor for C8/R^8 term
   real(wp), intent(out) :: d8
   !> Derivative of damping factor for C6/R^6 w.r.t the interatomic distance
   real(wp), intent(out) :: d6dr
   !> Derivative of damping factor for C8/R^8 w.r.t the interatomic distance
   real(wp), intent(out) :: d8dr

   real(wp) :: r4, r6, r8, rpoly, rpoly2, rpoly6, rpoly8

   r4 = r2 * r2
   r6 = r4 * r2
   r8 = r6 * r2

   ! First-order polynomial scaling of the damping radius
   rpoly = param%a1 * rdamp + param%a2

   rpoly2 = rpoly * rpoly
   rpoly6 = rpoly2 * rpoly2 * rpoly2
   rpoly8 = rpoly6 * rpoly2

   d6 = param%s6 / (r6 + rpoly6)
   d8 = param%s8 / (r8 + rpoly8)

   d6dr = -6.0_wp * r4 * d6**2 / param%s6
   d8dr = -8.0_wp * r6 * d8**2 / param%s8

end subroutine get_2b_derivs


!> Evaluate three-body damping factor
pure subroutine get_3b_damp(self, param, r, r2ij, r2ik, r2jk, &
   & rdamp, rdampij, rdampik, rdampjk, d9)
   !> Instance of the three-body damping function
   class(rational_damping_threebody), intent(in) :: self
   !> Damping parameters
   type(param_type), intent(in) :: param
   !> Product of pairwise interatomic distances
   real(wp), intent(in) :: r
   !> Square of interatomic distance between atoms i and j
   real(wp), intent(in) :: r2ij
   !> Square of interatomic distance between atoms i and k
   real(wp), intent(in) :: r2ik
   !> Square of interatomic distance between atoms j and k
   real(wp), intent(in) :: r2jk
   !> Product of pairwise damping radii
   real(wp), intent(in) :: rdamp
   !> Pairwise damping radius of atoms i and j
   real(wp), intent(in) :: rdampij
   !> Pairwise damping radius of atoms i and k
   real(wp), intent(in) :: rdampik
   !> Pairwise damping radius of atoms j and k
   real(wp), intent(in) :: rdampjk
   !> Damping factor for C9/R^9 term
   real(wp), intent(out) :: d9

   real(wp) :: r3ij, r3ik, r3jk, fij, fik, fjk, rpolyij, rpolyik, rpolyjk

   ! First-order polynomial scaling of the damping radius
   rpolyij = param%a1 * rdampij + param%a2
   rpolyik = param%a1 * rdampik + param%a2
   rpolyjk = param%a1 * rdampjk + param%a2

   r3ij = r2ij * sqrt(r2ij)
   r3ik = r2ik * sqrt(r2ik)
   r3jk = r2jk * sqrt(r2jk)

   fij = r3ij / (r3ij + rpolyij**3)
   fik = r3ik / (r3ik + rpolyik**3)
   fjk = r3jk / (r3jk + rpolyjk**3)

   d9 = param%s9 * fij * fik * fjk

end subroutine get_3b_damp


!> Evaluate three-body damping factor
pure subroutine get_3b_derivs(self, param, r, r2ij, r2ik, r2jk, &
   & rdamp, rdampij, rdampik, rdampjk, d9, d9drij, d9drik, d9drjk)
   !> Instance of the three-body damping function
   class(rational_damping_threebody), intent(in) :: self
   !> Damping parameters
   type(param_type), intent(in) :: param
   !> Product of pairwise interatomic distances
   real(wp), intent(in) :: r
   !> Square of interatomic distance between atoms i and j
   real(wp), intent(in) :: r2ij
   !> Square of interatomic distance between atoms i and k
   real(wp), intent(in) :: r2ik
   !> Square of interatomic distance between atoms j and k
   real(wp), intent(in) :: r2jk
   !> Product of pairwise damping radii
   real(wp), intent(in) :: rdamp
   !> Pairwise damping radius of atoms i and j
   real(wp), intent(in) :: rdampij
   !> Pairwise damping radius of atoms i and k
   real(wp), intent(in) :: rdampik
   !> Pairwise damping radius of atoms j and k
   real(wp), intent(in) :: rdampjk
   !> Damping factor for C9/R^9 term
   real(wp), intent(out) :: d9
   !> Derivative of damping factor for C9/R^9 
   !> w.r.t the interatomic distance between atoms i and j
   real(wp), intent(out) :: d9drij
   !> Derivative of damping factor for C9/R^9 
   !> w.r.t the interatomic distance between atoms i and k
   real(wp), intent(out) :: d9drik
   !> Derivative of damping factor for C9/R^9 
   !> w.r.t the interatomic distance between atoms j and k
   real(wp), intent(out) :: d9drjk

   real(wp) :: r3ij, r3ik, r3jk, fij, fik, fjk, pre, rpolyij, rpolyik, rpolyjk

   ! First-order polynomial scaling of the damping radius
   rpolyij = param%a1 * rdampij + param%a2
   rpolyik = param%a1 * rdampik + param%a2
   rpolyjk = param%a1 * rdampjk + param%a2

   r3ij = r2ij * sqrt(r2ij)
   r3ik = r2ik * sqrt(r2ik)
   r3jk = r2jk * sqrt(r2jk)

   fij = r3ij / (r3ij + rpolyij**3)
   fik = r3ik / (r3ik + rpolyik**3)
   fjk = r3jk / (r3jk + rpolyjk**3)

   d9 = param%s9 * fij * fik * fjk

   pre = 3.0_wp * d9
   d9drij = pre * (1.0_wp - fij) / r2ij
   d9drik = pre * (1.0_wp - fik) / r2ik
   d9drjk = pre * (1.0_wp - fjk) / r2jk

end subroutine get_3b_derivs


!> Check the availability of the parameters for two-body damping
subroutine check_2b_params(self, error, param)
   !> Instance of the two-body damping function
   class(rational_damping_twobody), intent(in) :: self
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Damping parameters
   type(param_type), intent(in) :: param

   if (.not. allocated(param%s6)) then
      call fatal_error(error, "Missing parameter in rational damping: s6")
   end if
   if (.not. allocated(param%s8)) then
      call fatal_error(error, "Missing parameter in rational damping: s8")
   end if

end subroutine check_2b_params


!> Check the availability of the parameters for three-body damping
subroutine check_3b_params(self, error, param)
   !> Instance of the three-body damping function
   class(rational_damping_threebody), intent(in) :: self
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Damping parameters
   type(param_type), intent(in) :: param

   if (.not. allocated(param%s9)) then
      call fatal_error(error, "Missing parameter in rational damping: s9")
   end if

end subroutine check_3b_params

end module dftd4_damping_rational
