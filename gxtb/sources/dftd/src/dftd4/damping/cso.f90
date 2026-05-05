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

!> Implementation of the C-Six-Only (CSO) damping function.
module dftd4_damping_cso
   use dftd4_damping_type, only : damping_twobody
   use dftd4_param_type, only : param_type
   use mctc_env, only : error_type, fatal_error, wp
   implicit none
   private

   public :: cso_damping_twobody, new_cso_damping_twobody

   !> C-Six-Only (CSO) damping for two-body terms
   type, extends(damping_twobody) :: cso_damping_twobody
   contains
      !> Evaluate two-body damping factor
      procedure :: get_2b_damp
      !> Evaluate two-body damping factor with derivatives
      procedure :: get_2b_derivs
      !> Check the availability of the parameters for two-body damping
      procedure :: check_2b_params
   end type cso_damping_twobody

   character(len=*), parameter :: cso_label = "C-Six-Only damping"
   character(len=*), parameter :: cso_label_short = "CSO"

contains


!> Create a new instance of the C-Six-Only damping for two-body terms
subroutine new_cso_damping_twobody(self)
   !> Instance of the two-body C-Six-Only damping function
   type(cso_damping_twobody), intent(out) :: self

   self%label = cso_label
   self%label_short = cso_label_short
end subroutine new_cso_damping_twobody


!> Evaluate two-body damping factor
pure subroutine get_2b_damp(self, param, r2, rdamp, d6, d8)
   !> Instance of the two-body damping function
   class(cso_damping_twobody), intent(in) :: self
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

   real(wp) :: r1, r6, a2_2, a2_6, ef, sf, sig

   r1 = sqrt(r2)
   r6 = r2 * r2 * r2

   ! In CSO the rational damping radii are replaced by only the constant a2
   a2_2 = param%a2 * param%a2
   a2_6 = a2_2 * a2_2 * a2_2

   ! Use the damping radius without linear scaling for the exponential function
   ef = exp(r1 - param%a4 * rdamp)
   sf = 1.0_wp / (1.0_wp + ef)
   sig = param%s6 + param%a3 * sf 

   d6 = sig / (r6 + a2_6)
   d8 = 0.0_wp

end subroutine get_2b_damp


!> Evaluate two-body damping factor with derivarives
pure subroutine get_2b_derivs(self, param, r2, rdamp, d6, d8, d6dr, d8dr)
   !> Instance of the two-body damping function
   class(cso_damping_twobody), intent(in) :: self
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

   real(wp) :: r1, r4, r6, a2_2, a2_6, ef, sf, sig, t6, dsig

   r1 = sqrt(r2)
   r4 = r2 * r2
   r6 = r4 * r2

   ! In CSO the rational damping radii are replaced by only the constant a2
   a2_2 = param%a2 * param%a2
   a2_6 = a2_2 * a2_2 * a2_2

   ! Use the damping radius without linear scaling for the exponential function
   ef = exp(r1 - param%a4 * rdamp)
   sf = 1.0_wp / (1.0_wp + ef)
   sig = param%s6 + param%a3 * sf

   t6 = 1.0_wp / (r6 + a2_6)

   d6 = sig * t6
   d8 = 0.0_wp

   dsig = -param%a3 * sf * (1.0_wp - sf) / r1

   d6dr = dsig * t6 - 6.0_wp * r4 * t6**2 * sig
   d8dr = 0.0_wp

end subroutine get_2b_derivs


!> Check the availability of the parameters for two-body damping
subroutine check_2b_params(self, error, param)
   !> Instance of the two-body damping function
   class(cso_damping_twobody), intent(in) :: self
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Damping parameters
   type(param_type), intent(in) :: param

   if (.not. allocated(param%s6)) then
      call fatal_error(error, "Missing parameter in C-Six-Only damping: s6")
   end if
   if (.not. allocated(param%a3)) then
      call fatal_error(error, "Missing parameter in C-Six-Only damping: a3")
   end if
   if (.not. allocated(param%a4)) then
      call fatal_error(error, "Missing parameter in C-Six-Only damping: a4")
   end if
   if (abs(param%a1) < epsilon(param%a1)) then
      call fatal_error(error, "a1 cannot be zero in C-Six-Only damping")
   end if

end subroutine check_2b_params

end module dftd4_damping_cso