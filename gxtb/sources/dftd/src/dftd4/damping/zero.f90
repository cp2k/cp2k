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

!> Implementation of the zero (Chai--Head-Gordon) damping function.
module dftd4_damping_zero
   use dftd4_damping_type, only : damping_twobody, damping_threebody
   use dftd4_param_type, only : param_type
   use mctc_env, only : error_type, fatal_error, wp
   implicit none
   private

   public :: zero_damping_twobody, new_zero_damping_twobody
   public :: zero_damping_threebody, new_zero_damping_threebody
   public :: zero_avg_damping_threebody, new_zero_avg_damping_threebody

   !> Zero (Chai-Head-Gordon) damping for two-body terms
   type, extends(damping_twobody) :: zero_damping_twobody
   contains
      !> Evaluate two-body damping factor
      procedure :: get_2b_damp
      !> Evaluate two-body damping factor with derivarives
      procedure :: get_2b_derivs
      !> Check the availability of the parameters for two-body damping
      procedure :: check_2b_params
   end type zero_damping_twobody

   !> Zero (Chai-Head-Gordon) damping for three-body terms
   type, extends(damping_threebody) :: zero_damping_threebody
   contains
      !> Evaluate three-body damping factor
      procedure :: get_3b_damp
      !> Evaluate three-body damping factor with derivarives
      procedure :: get_3b_derivs
      !> Check the availability of the parameters for three-body damping
      procedure :: check_3b_params
   end type zero_damping_threebody

   !> Zero (Chai-Head-Gordon) damping based on average distances for three-body terms
   type, extends(damping_threebody) :: zero_avg_damping_threebody
   contains
      !> Evaluate three-body damping factor
      procedure :: get_3b_damp => get_3b_damp_avg
      !> Evaluate three-body damping factor with derivarives
      procedure :: get_3b_derivs => get_3b_derivs_avg
      !> Check the availability of the parameters for three-body damping
      procedure :: check_3b_params => check_3b_params_avg
   end type zero_avg_damping_threebody

   character(len=*), parameter :: zero_label = "Zero (Chai-Head-Gordon) damping"
   character(len=*), parameter :: zero_label_short = "0"
   character(len=*), parameter :: zero_avg_label = "Averaged distance zero (Chai-Head-Gordon) damping"
   character(len=*), parameter :: zero_avg_label_short = "0avg"

contains


!> Create a new instance of the zero damping for two-body terms
subroutine new_zero_damping_twobody(self)
   !> Instance of the two-body zero damping function
   type(zero_damping_twobody), intent(out) :: self

   self%label = zero_label
   self%label_short = zero_label_short
end subroutine new_zero_damping_twobody


!> Create a new instance of the zero damping for three-body terms
subroutine new_zero_damping_threebody(self)
   !> Instance of the three-body zero damping function
   type(zero_damping_threebody), intent(out) :: self

   self%label = zero_label
   self%label_short = zero_label_short
end subroutine new_zero_damping_threebody


!> Create a new instance of the zero damping for three-body terms
subroutine new_zero_avg_damping_threebody(self)
   !> Instance of the three-body zero damping function
   type(zero_avg_damping_threebody), intent(out) :: self

   self%label = zero_avg_label
   self%label_short = zero_avg_label_short
end subroutine new_zero_avg_damping_threebody

!> Evaluate two-body damping factor
pure subroutine get_2b_damp(self, param, r2, rdamp, d6, d8)
   !> Instance of the two-body damping function
   class(zero_damping_twobody), intent(in) :: self
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

   real(wp) :: rpoly, ratior, r6, r8, t6, t8, inv6, inv8, alp6, alp8
   
   alp6 = param%alp
   alp8 = param%alp + 2.0_wp

   ! First-order polynomial scaling of the damping radius
   rpoly = param%a1 * rdamp + param%a2

   ratior = rpoly / sqrt(r2)
   r6 = r2*r2*r2
   r8 = r6*r2

   t6  = (param%rs6 * ratior)**alp6
   t8  = (param%rs8 * ratior)**alp8

   inv6 = 1.0_wp / (1.0_wp + 6.0_wp*t6)
   inv8 = 1.0_wp / (1.0_wp + 6.0_wp*t8)

   d6 = param%s6 * inv6 / r6
   d8 = param%s8 * inv8 / r8

end subroutine get_2b_damp


!> Evaluate two-body damping factor with derivarives
pure subroutine get_2b_derivs(self, param, r2, rdamp, d6, d8, d6dr, d8dr)
   !> Instance of the two-body damping function
   class(zero_damping_twobody), intent(in) :: self
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

   real(wp) :: rpoly, ratior, r6, r8, t6, t8, inv6, inv8, alp6, alp8

   alp6 = param%alp
   alp8 = param%alp + 2.0_wp

   ! First-order polynomial scaling of the damping radius
   rpoly = param%a1 * rdamp + param%a2

   ratior = rpoly / sqrt(r2)
   r6 = r2*r2*r2
   r8 = r6*r2

   t6  = (param%rs6 * ratior)**alp6
   t8  = (param%rs8 * ratior)**alp8
   inv6 = 1.0_wp / (1.0_wp + 6.0_wp*t6)
   inv8 = 1.0_wp / (1.0_wp + 6.0_wp*t8)

   d6 = param%s6 * inv6 / r6
   d8 = param%s8 * inv8 / r8

   d6dr = d6 * ( -6.0_wp/r2 + 6.0_wp * alp6 * t6 * inv6 / r2 )
   d8dr = d8 * ( -8.0_wp/r2 + 6.0_wp * alp8 * t8 * inv8 / r2 )

end subroutine get_2b_derivs


!> Evaluate three-body damping factor
pure subroutine get_3b_damp(self, param, r, r2ij, r2ik, r2jk, &
   & rdamp, rdampij, rdampik, rdampjk, d9)
   !> Instance of the three-body damping function
   class(zero_damping_threebody), intent(in) :: self
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
   real(wp) :: t3ij, t3ik, t3jk, inv3ij, inv3ik, inv3jk

   rij = sqrt(r2ij)
   rik = sqrt(r2ik)
   rjk = sqrt(r2jk)

   ! First-order polynomial scaling of the damping radius
   rpolyij = param%a1 * rdampij + param%a2
   rpolyik = param%a1 * rdampik + param%a2
   rpolyjk = param%a1 * rdampjk + param%a2

   t3ij = (param%rs9 * rpolyij / rij)**param%alp
   t3ik = (param%rs9 * rpolyik / rik)**param%alp
   t3jk = (param%rs9 * rpolyjk / rjk)**param%alp

   inv3ij = 1.0_wp / (1.0_wp + 6.0_wp * t3ij)
   inv3ik = 1.0_wp / (1.0_wp + 6.0_wp * t3ik)
   inv3jk = 1.0_wp / (1.0_wp + 6.0_wp * t3jk)

   d9 = param%s9 * inv3ij * inv3ik * inv3jk

end subroutine get_3b_damp


!> Evaluate three-body damping factor
pure subroutine get_3b_derivs(self, param, r, r2ij, r2ik, r2jk, &
   & rdamp, rdampij, rdampik, rdampjk, d9, d9drij, d9drik, d9drjk)
   !> Instance of the three-body damping function
   class(zero_damping_threebody), intent(in) :: self
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
   real(wp) :: t3ij, t3ik, t3jk, inv3ij, inv3ik, inv3jk, pre

   rij = sqrt(r2ij)
   rik = sqrt(r2ik)
   rjk = sqrt(r2jk)

   ! First-order polynomial scaling of the damping radius
   rpolyij = param%a1 * rdampij + param%a2
   rpolyik = param%a1 * rdampik + param%a2
   rpolyjk = param%a1 * rdampjk + param%a2


   t3ij = (param%rs9 * rpolyij / rij)**param%alp
   t3ik = (param%rs9 * rpolyik / rik)**param%alp
   t3jk = (param%rs9 * rpolyjk / rjk)**param%alp

   inv3ij = 1.0_wp / (1.0_wp + 6.0_wp * t3ij)
   inv3ik = 1.0_wp / (1.0_wp + 6.0_wp * t3ik)
   inv3jk = 1.0_wp / (1.0_wp + 6.0_wp * t3jk)

   d9 = param%s9 * inv3ij * inv3ik * inv3jk

   pre = 6.0_wp * param%alp * d9
   
   d9drij = pre * t3ij * inv3ij / r2ij
   d9drik = pre * t3ik * inv3ik / r2ik
   d9drjk = pre * t3jk * inv3jk / r2jk

end subroutine get_3b_derivs


!> Evaluate three-body damping factor based on distance product
pure subroutine get_3b_damp_avg(self, param, r, r2ij, r2ik, r2jk, &
   & rdamp, rdampij, rdampik, rdampjk, d9)
   !> Instance of the three-body damping function
   class(zero_avg_damping_threebody), intent(in) :: self
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

   real(wp) :: rpoly, rpolyij, rpolyik, rpolyjk

   ! First-order polynomial scaling of the damping radius
   rpolyij = param%a1 * rdampij + param%a2
   rpolyik = param%a1 * rdampik + param%a2
   rpolyjk = param%a1 * rdampjk + param%a2

   rpoly = rpolyij * rpolyik * rpolyjk

   d9 = param%s9 / (1.0_wp + 6.0_wp * (rpoly / r)**(param%alp / 3.0_wp))

end subroutine get_3b_damp_avg


!> Evaluate three-body damping factor based on distance product
pure subroutine get_3b_derivs_avg(self, param, r, r2ij, r2ik, r2jk, &
   & rdamp, rdampij, rdampik, rdampjk, d9, d9drij, d9drik, d9drjk)
   !> Instance of the three-body damping function
   class(zero_avg_damping_threebody), intent(in) :: self
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

   real(wp) :: rpolyij, rpolyik, rpolyjk, rpoly, dftmp, t9, inv9

   ! First-order polynomial scaling of the damping radius
   rpolyij = param%a1 * rdampij + param%a2
   rpolyik = param%a1 * rdampik + param%a2
   rpolyjk = param%a1 * rdampjk + param%a2

   rpoly = rpolyij * rpolyik * rpolyjk

   t9 = (rpoly / r)**(param%alp / 3.0_wp)

   inv9 = 1.0_wp / (1.0_wp + 6.0_wp * t9)

   d9 = param%s9 * inv9

   dftmp = 2.0_wp * param%alp * d9 * t9 * inv9
   d9drij = dftmp / r2ij
   d9drik = dftmp / r2ik
   d9drjk = dftmp / r2jk

end subroutine get_3b_derivs_avg


!> Check the availability of the parameters for two-body damping
subroutine check_2b_params(self, error, param)
   !> Instance of the two-body damping function
   class(zero_damping_twobody), intent(in) :: self
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Damping parameters
   type(param_type), intent(in) :: param

   if (.not. allocated(param%s6)) then
      call fatal_error(error, "Missing parameter in zero damping: s6")
   end if
   if (.not. allocated(param%s8)) then
      call fatal_error(error, "Missing parameter in zero damping: s8")
   end if
   if (.not. allocated(param%alp)) then
      call fatal_error(error, "Missing parameter in zero damping: alp")
   end if
   if (.not. allocated(param%rs6)) then
      call fatal_error(error, "Missing parameter in zero damping: rs6")
   end if
   if (.not. allocated(param%rs8)) then
      call fatal_error(error, "Missing parameter in zero damping: rs8")
   end if

end subroutine check_2b_params


!> Check the availability of the parameters for three-body damping
subroutine check_3b_params(self, error, param)
   !> Instance of the three-body damping function
   class(zero_damping_threebody), intent(in) :: self
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Damping parameters
   type(param_type), intent(in) :: param

   if (.not. allocated(param%s9)) then
      call fatal_error(error, "Missing parameter in zero damping: s9")
   end if
   if (.not. allocated(param%alp)) then
      call fatal_error(error, "Missing parameter in zero damping: alp")
   end if
   if (.not. allocated(param%rs9)) then
      call fatal_error(error, "Missing parameter in zero damping: rs9")
   end if

end subroutine check_3b_params


!> Check the availability of the parameters for three-body damping
subroutine check_3b_params_avg(self, error, param)
   !> Instance of the three-body damping function
   class(zero_avg_damping_threebody), intent(in) :: self
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Damping parameters
   type(param_type), intent(in) :: param

   if (.not. allocated(param%s9)) then
      call fatal_error(error, "Missing parameter in averaged zero damping: s9")
   end if
   if (.not. allocated(param%alp)) then
      call fatal_error(error, "Missing parameter in averaged zero damping: alp")
   end if
   if (.not. allocated(param%rs9)) then
      call fatal_error(error, "Missing parameter in averaged zero damping: rs9")
   end if

end subroutine check_3b_params_avg

end module dftd4_damping_zero