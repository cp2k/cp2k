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

!> Implementation of the modified zero (Sherrill) damping function.
module dftd4_damping_mzero
   use dftd4_damping_type, only : damping_twobody
   use dftd4_param_type, only : param_type
   use mctc_env, only : error_type, fatal_error, wp
   implicit none
   private

   public :: mzero_damping_twobody, new_mzero_damping_twobody

   !> Modified zero (Sherrill) damping for two-body terms
   type, extends(damping_twobody) :: mzero_damping_twobody
   contains
      !> Evaluate two-body damping factor
      procedure :: get_2b_damp
      !> Evaluate two-body damping factor with derivatives
      procedure :: get_2b_derivs
      !> Check the availability of the parameters for two-body damping
      procedure :: check_2b_params
   end type mzero_damping_twobody

   character(len=*), parameter :: mzero_label = "Modified zero damping"
   character(len=*), parameter :: mzero_label_short = "0M"

contains


!> Create a new instance of the modified zero damping for two-body terms
subroutine new_mzero_damping_twobody(self)
   !> Instance of the two-body modified zero damping function
   type(mzero_damping_twobody), intent(out) :: self

   self%label = mzero_label
   self%label_short = mzero_label_short
end subroutine new_mzero_damping_twobody


!> Evaluate two-body damping factor
pure subroutine get_2b_damp(self, param, r2, rdamp, d6, d8)
   !> Instance of the two-body damping function
   class(mzero_damping_twobody), intent(in) :: self
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

   real(wp) :: rpoly, r1, r6, r8, x6, x8, t6, t8, inv6, inv8, alp6, alp8

   ! First-order polynomial scaling of the damping radius
   rpoly = param%a1 * rdamp + param%a2

   alp6 = param%alp
   alp8 = param%alp + 2.0_wp

   r1 = sqrt(r2)
   r6 = r2*r2*r2
   r8 = r6*r2

   x6 = r1 / (param%rs6 * rpoly) + param%bet * rpoly
   x8 = r1 / (param%rs8 * rpoly) + param%bet * rpoly

   t6 = x6**(-alp6)
   t8 = x8**(-alp8)

   inv6 = 1.0_wp / (1.0_wp + 6.0_wp*t6)
   inv8 = 1.0_wp / (1.0_wp + 6.0_wp*t8)

   d6 = param%s6 * inv6 / r6
   d8 = param%s8 * inv8 / r8

end subroutine get_2b_damp


!> Evaluate two-body damping factor with derivarives
pure subroutine get_2b_derivs(self, param, r2, rdamp, d6, d8, d6dr, d8dr)
   !> Instance of the two-body damping function
   class(mzero_damping_twobody), intent(in) :: self
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

   real(wp) :: rpoly, r1, r6, r8, t6, t8, x6, x8, dt6, dt8, inv6, inv8, alp6, alp8

   ! First-order polynomial scaling of the damping radius
   rpoly = param%a1 * rdamp + param%a2

   alp6 = param%alp
   alp8 = param%alp + 2.0_wp

   r1 = sqrt(r2)
   r6 = r2*r2*r2
   r8 = r6*r2

   x6 = r1 / (param%rs6 * rpoly) + param%bet * rpoly
   x8 = r1 / (param%rs8 * rpoly) + param%bet * rpoly

   t6 = x6**(-alp6)
   t8 = x8**(-alp8)

   inv6 = 1.0_wp / (1.0_wp + 6.0_wp*t6)
   inv8 = 1.0_wp / (1.0_wp + 6.0_wp*t8)

   d6 = param%s6 * inv6 / r6
   d8 = param%s8 * inv8 / r8

   dt6 = r1 / (x6 * param%rs6 * rpoly)
   dt8 = r1 / (x8 * param%rs8 * rpoly)

   d6dr = (d6 / r2) * (-6.0_wp + 6.0_wp * alp6 * t6 * inv6 * dt6)
   d8dr = (d8 / r2) * (-8.0_wp + 6.0_wp * alp8 * t8 * inv8 * dt8)

end subroutine get_2b_derivs


!> Check the availability of the parameters for two-body damping
subroutine check_2b_params(self, error, param)
   !> Instance of the two-body damping function
   class(mzero_damping_twobody), intent(in) :: self
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Damping parameters
   type(param_type), intent(in) :: param

   if (.not. allocated(param%s6)) then
      call fatal_error(error, "Missing parameter in modified zero damping: s6")
   end if
   if (.not. allocated(param%s8)) then
      call fatal_error(error, "Missing parameter in modified zero damping: s8")
   end if
   if (.not. allocated(param%alp)) then
      call fatal_error(error, "Missing parameter in modified zero damping: alp")
   end if
   if (.not. allocated(param%bet)) then
      call fatal_error(error, "Missing parameter in modified zero damping: bet")
   end if
   if (.not. allocated(param%rs6)) then
      call fatal_error(error, "Missing parameter in modified zero damping: rs6")
   end if
   if (.not. allocated(param%rs8)) then
      call fatal_error(error, "Missing parameter in modified zero damping: rs8")
   end if

end subroutine check_2b_params

end module dftd4_damping_mzero