! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
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

!> Generic interface to define damping functions for the DFT-D models
module dftd4_damping_type
   use dftd4_param_type, only : param_type
   use mctc_env, only : error_type, wp
   implicit none
   private

   public :: damping_type, damping_twobody, damping_threebody
   public :: twobody_damping_function, threebody_damping_function

   !> Abstract base type for damping parameterizations
   type :: damping_type
      !> Damping function for two-body terms
      class(damping_twobody), allocatable :: damping_2b
      !> Damping function for three-body terms
      class(damping_threebody), allocatable :: damping_3b
   contains
      !> Check the availability of parameters for the damping functions
      procedure :: check_params
   end type damping_type

   type, abstract :: damping_twobody
      !> Label identifying this two-body damping function
      character(len=:), allocatable :: label
      !> Short form label identifying this two-body damping function
      character(len=:), allocatable :: label_short
   contains
      !> Evaluate two-body damping factor
      procedure(get_2b_damp), deferred :: get_2b_damp
      !> Evaluate two-body damping factor with derivarives
      procedure(get_2b_derivs), deferred :: get_2b_derivs
      !> Check the availability of the parameters for two-body damping
      procedure(check_2b_params), deferred :: check_2b_params
   end type damping_twobody

   type, abstract :: damping_threebody
      !> Label identifying this three-body damping function
      character(len=:), allocatable :: label
      !> Short form label identifying this three-body damping function
      character(len=:), allocatable :: label_short
   contains
      !> Evaluate three-body damping factor
      procedure(get_3b_damp), deferred :: get_3b_damp
      !> Evaluate three-body damping factor with derivarives
      procedure(get_3b_derivs), deferred :: get_3b_derivs
      !> Check the availability of the parameters for three-body damping
      procedure(check_3b_params), deferred :: check_3b_params
   end type damping_threebody

   abstract interface
      !> Evaluate two-body damping factor
      pure subroutine get_2b_damp(self, param, r2, rdamp, d6, d8)
         import :: damping_twobody, param_type, wp
         !> Instance of the two-body damping function
         class(damping_twobody), intent(in) :: self
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
      end subroutine get_2b_damp

      !> Evaluate two-body damping factor with derivarives
      pure subroutine get_2b_derivs(self, param, r2, rdamp, d6, d8, d6dr, d8dr)
         import :: damping_twobody, param_type, wp
         !> Instance of the two-body damping function
         class(damping_twobody), intent(in) :: self
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
      end subroutine get_2b_derivs

      !> Check the availability of the parameters for two-body damping
      subroutine check_2b_params(self, error, param)
         import :: damping_twobody, param_type, error_type
         !> Instance of the two-body damping function
         class(damping_twobody), intent(in) :: self
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
         !> Damping parameters
         type(param_type), intent(in) :: param
      end subroutine check_2b_params

      !> Evaluate three-body damping factor
      pure subroutine get_3b_damp(self, param, r, r2ij, r2ik, r2jk, &
         & rdamp, rdampij, rdampik, rdampjk, d9)
         import :: damping_threebody, param_type, wp
         !> Instance of the three-body damping function
         class(damping_threebody), intent(in) :: self
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
      end subroutine get_3b_damp

      !> Evaluate three-body damping factor
      pure subroutine get_3b_derivs(self, param, r, r2ij, r2ik, r2jk, &
         & rdamp, rdampij, rdampik, rdampjk, d9, d9drij, d9drik, d9drjk)
         import :: damping_threebody, param_type, wp
         !> Instance of the three-body damping function
         class(damping_threebody), intent(in) :: self
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
      end subroutine get_3b_derivs

      !> Check the availability of the parameters for three-body damping
      subroutine check_3b_params(self, error, param)
         import :: damping_threebody, param_type, error_type
         !> Instance of the three-body damping function
         class(damping_threebody), intent(in) :: self
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
         !> Damping parameters
         type(param_type), intent(in) :: param
      end subroutine check_3b_params
   end interface


   !> Possible two-body damping functions
   type :: enum_twobody_damping_function
      !> Rational Becke-Johnson damping
      integer :: rational = 1
      !> Screened rational Becke-Johnson damping
      integer :: screened = 2
      !> Chai-Head-Gordon zero-damping
      integer :: zero = 3
      !> Modified Chai-Head-Gordon zero-damping
      integer :: mzero = 4
      !> Optimized power damping
      integer :: optpower = 5
      !> C-Six-Only damping
      integer :: cso = 6
      !> Spherical wave based dispersion damping by Koide
      integer :: koide = 7
   end type enum_twobody_damping_function

   !> Actual enumerator for possible two-body damping functions
   type(enum_twobody_damping_function), parameter :: twobody_damping_function = &
      & enum_twobody_damping_function()


   !> Possible two-body damping functions
   type :: enum_threebody_damping_function
      !> Rational Becke-Johnson damping
      integer :: rational = 1
      !> Screened rational Becke-Johnson damping
      integer :: screened = 2
      !> Chai-Head-Gordon zero-damping
      integer :: zero = 3
      !> Average of the Chai-Head-Gordon zero-damping for three-body terms
      integer :: zero_avg = 4
   end type enum_threebody_damping_function

   !> Actual enumerator for possible two-body damping functions
   type(enum_threebody_damping_function), parameter :: threebody_damping_function = &
      & enum_threebody_damping_function()

contains

!> Check the availability of parameters for the damping functions
subroutine check_params(self, error, param)
   !> Damping object
   class(damping_type), intent(in) :: self
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Damping parameters
   type(param_type), intent(in) :: param

   call self%damping_2b%check_2b_params(error, param)
   if (allocated(error)) return

   if (allocated(self%damping_3b)) then
      call self%damping_3b%check_3b_params(error, param)
      if (allocated(error)) return
   end if

end subroutine check_params

end module dftd4_damping_type
