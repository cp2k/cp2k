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

!> Integrator for Casimir-Polder dispersion coefficients
module dftd4_integrator_type
   use mctc_env, only : wp
   implicit none
   private

   public :: integrator_type

   !> Abstract integrator type
   type, abstract :: integrator_type
      !> Number of grid points
      integer :: ngrid
      !> Frequencies
      real(wp), allocatable :: freq(:)
      !> Weights
      real(wp), allocatable :: weights(:)
   contains
      !> Integration interface
      procedure(integrate), deferred :: integrate
   end type integrator_type

   abstract interface
      !> Interface for integration
      function integrate(self, pol) result(val)
         import :: integrator_type, wp
         !> Integrator instance
         class(integrator_type), intent(in) :: self
         !> Polarizabilities at imaginary frequencies
         real(wp), intent(in) :: pol(:)
         !> Integrated value
         real(wp) :: val
      end function integrate
   end interface

end module dftd4_integrator_type
