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

module dftd4_param_type
   use mctc_env, only : wp
   implicit none
   private

   public :: param_type

   type :: param_type
      !> Scaling factor for C6/R^6 term
      real(wp), allocatable :: s6
      !> Scaling factor for C8/R^8 term
      real(wp), allocatable :: s8
      !> Scaling factor for C9/R^9 term
      real(wp), allocatable :: s9

      !> Linear damping radius dependence (always present for radii adaptation)
      !> used in: rational, cso
      real(wp) :: a1 = 1.0_wp
      !> Constant damping radius shift (always present for radii adaptation)
      !> used in: rational, cso
      real(wp) :: a2 = 0.0_wp

      !> Damping radius scaling in screening function 
      !> Short-range s6 modification (a1 in the original publication)
      !> used in: screened, cso
      real(wp), allocatable :: a3
      !> Screening function exponent
      !> Radius scaling for short-range s6 modification (a2 in the original publication)
      !> used in: screened, cso
      real(wp), allocatable :: a4
      !> Distance-radius-fraction scaling factor for C6/R^6 term
      !> used in: zero, mzero, koide
      real(wp), allocatable :: rs6
      !> Distance-radius-fraction scaling factor for C8/R^8 term 
      !> used in: zero, mzero, koide
      real(wp), allocatable :: rs8
      !> Distance-radius-fraction scaling factor for C9/R^9 term 
      !> used in: zero
      real(wp), allocatable :: rs9
      !> Zero-damping exponent parameter
      !> used in: zero, mzero
      real(wp), allocatable :: alp
      !> Additional linear damping radius dependence in modified zero damping
      !> Exponent for of rational damping function for optimized power damping
      !> used in: mzero, optpower
      real(wp), allocatable :: bet
   contains
      !> Reset all damping parameters
      procedure :: reset
   end type param_type

   !> Provide constructor for parameter type
   interface param_type
      module procedure :: create_param
   end interface param_type


contains

!> Consturctor for parameter opject
pure function create_param(s6, s8, s9, a1, a2, a3, a4, rs6, rs8, rs9, alp, bet) result(self)
   !> Scaling factor for C6/R^6 term
   real(wp), intent(in), optional :: s6
   !> Scaling factor for C8/R^8 term
   real(wp), intent(in), optional :: s8
   !> Scaling factor for C9/R^9 term
   real(wp), intent(in), optional :: s9
   !> Linear damping radius dependence
   real(wp), intent(in), optional :: a1
   !> Constant damping radius shift
   real(wp), intent(in), optional :: a2
   !> Damping radius scaling in screening function or short-range s6 modification
   real(wp), intent(in), optional :: a3
   !> Screening function exponent or radius scaling for short-range s6 modification
   real(wp), intent(in), optional :: a4
   !> Radius fraction scaling factor for C6/R^6 term
   real(wp), intent(in), optional :: rs6
   !> Radius fraction scaling factor for C8/R^8 term
   real(wp), intent(in), optional :: rs8
   !> Radius fraction scaling factor for C9/R^9 term 
   real(wp), intent(in), optional :: rs9
   !> Zero-damping exponent parameter
   real(wp), intent(in), optional :: alp
   !> Additional linear damping radius dependence in modified zero damping
   !> Exponent for of rational damping function for optimized power damping
   real(wp), intent(in), optional :: bet
   !> Parameter object
   type(param_type) :: self

   if (present(s6)) self%s6 = s6
   if (present(s8)) self%s8 = s8
   if (present(s9)) self%s9 = s9
   if (present(a1)) self%a1 = a1
   if (present(a2)) self%a2 = a2
   if (present(a3)) self%a3 = a3
   if (present(a4)) self%a4 = a4
   if (present(rs6)) self%rs6 = rs6
   if (present(rs8)) self%rs8 = rs8
   if (present(rs9)) self%rs9 = rs9
   if (present(alp)) self%alp = alp
   if (present(bet)) self%bet = bet

end function create_param

!> Reset all damping parameters
subroutine reset(self)
   !> Damping parameters
   class(param_type), intent(inout) :: self

   if (allocated(self%s6)) deallocate(self%s6)
   if (allocated(self%s8)) deallocate(self%s8)
   if (allocated(self%s9)) deallocate(self%s9)
   if (allocated(self%a3)) deallocate(self%a3)
   if (allocated(self%a4)) deallocate(self%a4)
   if (allocated(self%rs6)) deallocate(self%rs6)
   if (allocated(self%rs8)) deallocate(self%rs8)
   if (allocated(self%rs9)) deallocate(self%rs9)
   if (allocated(self%alp)) deallocate(self%alp)
   if (allocated(self%bet)) deallocate(self%bet)
   self%a1 = 1.0_wp
   self%a2 = 0.0_wp

end subroutine reset

end module dftd4_param_type