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

!> Implementation of the optimized power (Head-Gordon) damping function.
module dftd4_damping_optpower
   use dftd4_damping_type, only : damping_twobody
   use dftd4_param_type, only : param_type
   use mctc_env, only : error_type, fatal_error, wp
   implicit none
   private

   public :: optpower_damping_twobody, new_optpower_damping_twobody

   !> Optimized power (Head-Gordon) damping for two-body terms
   type, extends(damping_twobody) :: optpower_damping_twobody
   contains
      !> Evaluate two-body damping factor
      procedure :: get_2b_damp
      !> Evaluate two-body damping factor with derivatives
      procedure :: get_2b_derivs
      !> Check the availability of the parameters for two-body damping
      procedure :: check_2b_params
   end type optpower_damping_twobody

   character(len=*), parameter :: optpower_label = "Optimized power damping"
   character(len=*), parameter :: optpower_label_short = "op"

contains


!> Create a new instance of the optimized power damping for two-body terms
subroutine new_optpower_damping_twobody(self)
   !> Instance of the two-body screened damping function
   type(optpower_damping_twobody), intent(out) :: self

   self%label = optpower_label
   self%label_short = optpower_label_short
end subroutine new_optpower_damping_twobody


!> Evaluate two-body damping factor
pure subroutine get_2b_damp(self, param, r2, rdamp, d6, d8)
   !> Instance of the two-body damping function
   class(optpower_damping_twobody), intent(in) :: self
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

   real(wp) :: r6, r8, rpoly, rpoly2, rpoly6, rpoly8, rb, ab, t6, t8

   r6 = r2 * r2 * r2
   r8 = r6 * r2

   ! First-order polynomial scaling of the damping radius
   rpoly = param%a1 * rdamp + param%a2

   rpoly2 = rpoly * rpoly
   rpoly6 = rpoly2 * rpoly2 * rpoly2
   rpoly8 = rpoly6 * rpoly2

   rb = r2**(param%bet * 0.5_wp)
   ab = rpoly**param%bet

   t6 = rb / (rb * r6 + ab * rpoly6)
   t8 = rb / (rb * r8 + ab * rpoly8)

   d6 = param%s6 * t6
   d8 = param%s8 * t8

end subroutine get_2b_damp


!> Evaluate two-body damping factor with derivarives
pure subroutine get_2b_derivs(self, param, r2, rdamp, d6, d8, d6dr, d8dr)
   !> Instance of the two-body damping function
   class(optpower_damping_twobody), intent(in) :: self
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

   real(wp) :: r4, r6, r8, rpoly, rpoly2, rpoly6, rpoly8, rb, ab, t6, t8, pre

   r4 = r2 * r2
   r6 = r4 * r2
   r8 = r6 * r2

   ! First-order polynomial scaling of the damping radius
   rpoly = param%a1 * rdamp + param%a2

   rpoly2 = rpoly * rpoly
   rpoly6 = rpoly2 * rpoly2 * rpoly2
   rpoly8 = rpoly6 * rpoly2

   rb = r2**(param%bet * 0.5_wp)
   ab = rpoly**param%bet

   t6 = rb / (rb * r6 + ab * rpoly6)
   t8 = rb / (rb * r8 + ab * rpoly8)

   d6 = param%s6 * t6
   d8 = param%s8 * t8

   pre = param%bet * ab / (rb * r2)
   d6dr = param%s6 * t6**2 * (-6.0_wp * r4 + pre * rpoly6)
   d8dr = param%s8 * t8**2 * (-8.0_wp * r6 + pre * rpoly8)

end subroutine get_2b_derivs


!> Check the availability of the parameters for two-body damping
subroutine check_2b_params(self, error, param)
   !> Instance of the two-body damping function
   class(optpower_damping_twobody), intent(in) :: self
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Damping parameters
   type(param_type), intent(in) :: param

   if (.not. allocated(param%s6)) then
      call fatal_error(error, "Missing parameter in optimized power damping: s6")
   end if
   if (.not. allocated(param%s8)) then
      call fatal_error(error, "Missing parameter in optimized power damping: s8")
   end if
   if (.not. allocated(param%bet)) then
      call fatal_error(error, "Missing parameter in optimized power damping: bet")
   end if

end subroutine check_2b_params

end module dftd4_damping_optpower