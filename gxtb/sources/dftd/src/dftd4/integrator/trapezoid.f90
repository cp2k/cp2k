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

!> Trapezoid integrator for Casimir-Polder dispersion coefficients
module dftd4_integrator_trapezoid
   use dftd4_integrator_type , only : integrator_type
   use mctc_env, only : wp
   implicit none
   private

   public :: trapezoid_integrator, new_trapezoid_integrator

   !> Trapezoidal integrator
   type, extends(integrator_type) :: trapezoid_integrator
   contains
      procedure :: integrate => trapzoid_integrate
   end type trapezoid_integrator

contains

   !> Factory for numerical trapezoid integrator
   subroutine new_trapezoid_integrator(self, ngrid, freq)
      !> New integrator instance
      type(trapezoid_integrator), intent(out) :: self
      !> Number of grid points
      integer, intent(in) :: ngrid
      !> Frequencies
      real(wp), intent(in) :: freq(:)

      integer :: igrid

      self%ngrid = ngrid
      self%freq = freq

      allocate(self%weights(ngrid))
      self%weights(1) = 0.5_wp * (freq(2) - freq(1))
      do igrid = 2, ngrid-1
         self%weights(igrid) = 0.5_wp * (freq(igrid+1) - freq(igrid-1))
      end do
      self%weights(ngrid) = 0.5_wp * (freq(ngrid) - freq(ngrid-1))

   end subroutine new_trapezoid_integrator

   !> Numerical Casimir--Polder integration using trapezoid rule
   function trapzoid_integrate(self, pol) result(val)
      !> Integrator instance
      class(trapezoid_integrator), intent(in) :: self
      !> Polarizabilities at imaginary frequencies
      real(wp), intent(in) :: pol(:)
      !> Integrated value
      real(wp) :: val
      
      val = sum(pol(1:self%ngrid)*self%weights)

   end function trapzoid_integrate

end module dftd4_integrator_trapezoid
