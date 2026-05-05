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

!> Implementation of the screened rational (Becke--Johnson) damping function.
module dftd4_damping_screened
   use dftd4_damping_type, only : damping_twobody, damping_threebody
   use dftd4_param_type, only : param_type
   use mctc_env, only : error_type, fatal_error, wp
   use mctc_io_constants, only : pi
   implicit none
   private

   public :: screened_damping_twobody, new_screened_damping_twobody
   public :: screened_damping_threebody, new_screened_damping_threebody

   !> Screened rational (Becke-Johnson) damping for two-body terms
   type, extends(damping_twobody) :: screened_damping_twobody
   contains
      !> Evaluate two-body damping factor
      procedure :: get_2b_damp
      !> Evaluate two-body damping factor with derivarives
      procedure :: get_2b_derivs
      !> Check the availability of the parameters for two-body damping
      procedure :: check_2b_params
   end type screened_damping_twobody

   !> Screened rational (Becke-Johnson) damping for three-body terms
   type, extends(damping_threebody) :: screened_damping_threebody
   contains
      !> Evaluate three-body damping factor
      procedure :: get_3b_damp
      !> Evaluate three-body damping factor with derivarives
      procedure :: get_3b_derivs
      !> Check the availability of the parameters for three-body damping
      procedure :: check_3b_params
   end type screened_damping_threebody

   real(wp), parameter :: sqrtpi = sqrt(pi)

   character(len=*), parameter :: screened_label = "Screened rational (Becke-Johnson) damping"
   character(len=*), parameter :: screened_label_short = "scBJ"

contains


!> Create a new instance of the screened damping for two-body terms
subroutine new_screened_damping_twobody(self)
   !> Instance of the two-body screened damping function
   type(screened_damping_twobody), intent(out) :: self

   self%label = screened_label
   self%label_short = screened_label_short
end subroutine new_screened_damping_twobody


!> Create a new instance of the screened damping for three-body terms
subroutine new_screened_damping_threebody(self)
   !> Instance of the three-body screened damping function
   type(screened_damping_threebody), intent(out) :: self

   self%label = screened_label
   self%label_short = screened_label_short
end subroutine new_screened_damping_threebody


!> Evaluate two-body damping factor
pure subroutine get_2b_damp(self, param, r2, rdamp, d6, d8)
   !> Instance of the two-body damping function
   class(screened_damping_twobody), intent(in) :: self
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

   real(wp) :: r1, rpoly, arg, erf_val, reff, inv, inv2, inv6

   r1 = sqrt(r2)

   ! First-order polynomial scaling of the damping radius
   rpoly = param%a1 * rdamp + param%a2

   arg = -param%a3 * (r1 - param%a4 * rdamp)
   erf_val = 0.5_wp * (1.0_wp + erf(arg))
   reff = rpoly * erf_val

   inv = 1.0_wp / (r1 + reff)   
   inv2 = inv * inv
   inv6 = inv2 * inv2 * inv2

   d6 = param%s6 * inv6
   d8 = param%s8 * inv6 * inv2

end subroutine get_2b_damp


!> Evaluate two-body damping factor with derivarives
pure subroutine get_2b_derivs(self, param, r2, rdamp, d6, d8, d6dr, d8dr)
   !> Instance of the two-body damping function
   class(screened_damping_twobody), intent(in) :: self
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

   real(wp) :: r1, rpoly, arg, erf_val, reff, inv, inv2, inv6,  dreff

   r1 = sqrt(r2)

   ! First-order polynomial scaling of the damping radius
   rpoly = param%a1 * rdamp + param%a2

   arg = -param%a3 * (r1 - param%a4 * rdamp)
   erf_val = 0.5_wp * (1.0_wp + erf(arg))
   reff = rpoly * erf_val

   inv = 1.0_wp / (r1 + reff)
   inv2 = inv * inv
   inv6 = inv2 * inv2 * inv2

   d6 = param%s6 * inv6
   d8 = param%s8 * inv6 * inv2

   dreff = -rpoly * param%a3 * exp(-arg * arg) / sqrtpi

   d6dr = -6.0_wp * d6 * inv * (1.0_wp + dreff) / r1
   d8dr = -8.0_wp * d8 * inv * (1.0_wp + dreff) / r1

end subroutine get_2b_derivs


!> Evaluate three-body damping factor
pure subroutine get_3b_damp(self, param, r, r2ij, r2ik, r2jk, &
   & rdamp, rdampij, rdampik, rdampjk, d9)
   !> Instance of the three-body damping function
   class(screened_damping_threebody), intent(in) :: self
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

   real(wp) :: rij, rik, rjk, rpolyij, rpolyik, rpolyjk
   real(wp) :: arg_ij, arg_ik, arg_jk, reff_ij, reff_ik, reff_jk
   real(wp) :: fij, fik, fjk

   rij = sqrt(r2ij)
   rik = sqrt(r2ik)
   rjk = sqrt(r2jk)

   ! First-order polynomial scaling of the damping radius
   rpolyij = param%a1 * rdampij + param%a2
   rpolyik = param%a1 * rdampik + param%a2
   rpolyjk = param%a1 * rdampjk + param%a2

   arg_ij = -param%a3 * (rij - param%a4 * rdampij)
   arg_ik = -param%a3 * (rik - param%a4 * rdampik)
   arg_jk = -param%a3 * (rjk - param%a4 * rdampjk)

   reff_ij = rpolyij * 0.5_wp * (1.0_wp + erf(arg_ij))
   reff_ik = rpolyik * 0.5_wp * (1.0_wp + erf(arg_ik))
   reff_jk = rpolyjk * 0.5_wp * (1.0_wp + erf(arg_jk))

   fij = (rij / (rij + reff_ij))**3
   fik = (rik / (rik + reff_ik))**3
   fjk = (rjk / (rjk + reff_jk))**3

   d9 = param%s9 * fij * fik * fjk

end subroutine get_3b_damp


!> Evaluate three-body damping factor
pure subroutine get_3b_derivs(self, param, r, r2ij, r2ik, r2jk, &
   & rdamp, rdampij, rdampik, rdampjk, d9, d9drij, d9drik, d9drjk)
   !> Instance of the three-body damping function
   class(screened_damping_threebody), intent(in) :: self
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

   real(wp) :: rij, rik, rjk, rpolyij, rpolyik, rpolyjk
   real(wp) :: arg_ij, arg_ik, arg_jk, reff_ij, reff_ik, reff_jk
   real(wp) :: fij, fik, fjk, dreff_ij, dreff_ik, dreff_jk
   real(wp) :: pre_ij, pre_ik, pre_jk

   rij = sqrt(r2ij)
   rik = sqrt(r2ik)
   rjk = sqrt(r2jk)

   ! First-order polynomial scaling of the damping radius
   rpolyij = param%a1 * rdampij + param%a2
   rpolyik = param%a1 * rdampik + param%a2
   rpolyjk = param%a1 * rdampjk + param%a2

   arg_ij = -param%a3 * (rij - param%a4 * rdampij)
   arg_ik = -param%a3 * (rik - param%a4 * rdampik)
   arg_jk = -param%a3 * (rjk - param%a4 * rdampjk)

   reff_ij = rpolyij * 0.5_wp * (1.0_wp + erf(arg_ij))
   reff_ik = rpolyik * 0.5_wp * (1.0_wp + erf(arg_ik))
   reff_jk = rpolyjk * 0.5_wp * (1.0_wp + erf(arg_jk))

   fij = (rij / (rij + reff_ij))**3
   fik = (rik / (rik + reff_ik))**3
   fjk = (rjk / (rjk + reff_jk))**3

   d9 = param%s9 * fij * fik * fjk

   dreff_ij = -rpolyij * param%a3 * exp(-arg_ij * arg_ij) / sqrtpi
   dreff_ik = -rpolyik * param%a3 * exp(-arg_ik * arg_ik) / sqrtpi
   dreff_jk = -rpolyjk * param%a3 * exp(-arg_jk * arg_jk) / sqrtpi

   pre_ij = 3.0_wp * d9 / r2ij
   pre_ik = 3.0_wp * d9 / r2ik
   pre_jk = 3.0_wp * d9 / r2jk

   d9drij = pre_ij * (reff_ij - rij * dreff_ij) / (rij + reff_ij)
   d9drik = pre_ik * (reff_ik - rik * dreff_ik) / (rik + reff_ik)
   d9drjk = pre_jk * (reff_jk - rjk * dreff_jk) / (rjk + reff_jk)

end subroutine get_3b_derivs


!> Check the availability of the parameters for two-body damping
subroutine check_2b_params(self, error, param)
   !> Instance of the two-body damping function
   class(screened_damping_twobody), intent(in) :: self
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Damping parameters
   type(param_type), intent(in) :: param

   if (.not. allocated(param%s6)) then
      call fatal_error(error, "Missing parameter in screened damping: s6")
   end if
   if (.not. allocated(param%s8)) then
      call fatal_error(error, "Missing parameter in screened damping: s8")
   end if
   if (.not. allocated(param%a3)) then
      call fatal_error(error, "Missing parameter in screened damping: a3")
   end if
   if (.not. allocated(param%a4)) then
      call fatal_error(error, "Missing parameter in screened damping: a4")
   end if

end subroutine check_2b_params


!> Check the availability of the parameters for three-body damping
subroutine check_3b_params(self, error, param)
   !> Instance of the three-body damping function
   class(screened_damping_threebody), intent(in) :: self
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Damping parameters
   type(param_type), intent(in) :: param

   if (.not. allocated(param%s9)) then
      call fatal_error(error, "Missing parameter in screened damping: s9")
   end if
   if (.not. allocated(param%a3)) then
      call fatal_error(error, "Missing parameter in screened damping: a3")
   end if
   if (.not. allocated(param%a4)) then
      call fatal_error(error, "Missing parameter in screened damping: a4")
   end if

end subroutine check_3b_params

end module dftd4_damping_screened
