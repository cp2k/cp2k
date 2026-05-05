! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later
!
! tblite is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! tblite is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with tblite.  If not, see <https://www.gnu.org/licenses/>.

!> @file tblite/utils/average.f90
!> Contains averaging utilities

module tblite_utils_average
   use mctc_env, only : wp
   implicit none
   private

   public :: average_type, new_average, average_id

   !> Averager including derivatives
   type :: average_type
      !> Average of two parameters
      procedure(average_interface), pointer, nopass :: value => null()
      !> Derivative of the average w.r.t. the first parameter
      procedure(average_deriv_interface), pointer, nopass :: deriv => null()
   end type average_type

   abstract interface
      !> Average parameters
      pure function average_interface(gi, gj, interpolate) result(gij)
         import :: wp
         !> Parameter i
         real(wp), intent(in) :: gi
         !> Parameter j
         real(wp), intent(in) :: gj
         !> Parameter to smoothely interpolate between different averages 
         real(wp), intent(in), optional :: interpolate
         !> Averaged parameter
         real(wp) :: gij
      end function average_interface
      !> Parameter derivative of the average
      pure function average_deriv_interface(gi, gj, interpolate) result(dgij)
         import :: wp
         !> Parameter i
         real(wp), intent(in) :: gi
         !> Parameter j
         real(wp), intent(in) :: gj
         !> Parameter to smoothely interpolate between different averages 
         real(wp), intent(in), optional :: interpolate
         !> Derivative of the averaged parameter with respect to gi
         real(wp) :: dgij
      end function average_deriv_interface
   end interface

   !> Possible averaging functions
   type :: enum_average_id
      !> Arithmetic average
      integer :: arithmetic = 0
      !> Geometric average
      integer :: geometric = 1
      !> Harmonic average
      integer :: harmonic = 2
      !> General average, interpolated between 
      !> arithmetic (0), geometric (1), and harmonic (2)
      integer :: general = 3
   end type enum_average_id

   !> Actual enumerator for possible averaging functions
   type(enum_average_id), parameter :: average_id = enum_average_id()

contains

!> Factory to create a new average 
subroutine new_average(self, kind)
   !> Averager object
   type(average_type), intent(out) :: self
   !> Kind of averaging function
   integer, intent(in)  :: kind

   select case(kind)
   case default
      nullify(self%value)
      nullify(self%deriv)
   case (average_id%arithmetic)
      self%value => arithmetic_average
      self%deriv => arithmetic_average_deriv
   case (average_id%geometric)
      self%value => geometric_average
      self%deriv => geometric_average_deriv
   case (average_id%harmonic)
      self%value => harmonic_average
      self%deriv => harmonic_average_deriv
   case (average_id%general)
      self%value => general_average
      self%deriv => general_average_deriv
   end select
end subroutine new_average


!> Arithmetic averaging function
pure function arithmetic_average(gi, gj, interpolate) result(gij)
   !> Parameter i
   real(wp), intent(in) :: gi
   !> Parameter j
   real(wp), intent(in) :: gj
   !> Parameter to smoothely interpolate between different averages (not used)
   real(wp), intent(in), optional :: interpolate
   !> Averaged parameter
   real(wp) :: gij

   gij = 0.5_wp*(gi+gj)

end function arithmetic_average

!> Arithmetic averaging function derivative
pure function arithmetic_average_deriv(gi, gj, interpolate) result(dgij)
   !> Parameter i
   real(wp), intent(in) :: gi
   !> Parameter j
   real(wp), intent(in) :: gj
   !> Parameter to smoothely interpolate between different averages 
   real(wp), intent(in), optional :: interpolate
   !> Derivative of the averaged parameter with respect to gi
   real(wp) :: dgij

   dgij = 0.5_wp

end function arithmetic_average_deriv


!> Geometric averaging function
pure function geometric_average(gi, gj, interpolate) result(gij)
   !> Parameter i
   real(wp), intent(in) :: gi
   !> Parameter j
   real(wp), intent(in) :: gj
   !> Parameter to smoothely interpolate between different averages (not used)
   real(wp), intent(in), optional :: interpolate
   !> Averaged parameter
   real(wp) :: gij

   gij = sqrt(gi*gj)

end function geometric_average

!> Geometric averaging function derivative
pure function geometric_average_deriv(gi, gj, interpolate) result(dgij)
   !> Parameter i
   real(wp), intent(in) :: gi
   !> Parameter j
   real(wp), intent(in) :: gj
   !> Parameter to smoothely interpolate between different averages 
   real(wp), intent(in), optional :: interpolate
   !> Derivative of the averaged parameter with respect to gi
   real(wp) :: dgij

   dgij = 0.5_wp*sqrt(gj/gi)

end function geometric_average_deriv


!> Harmonic averaging function
pure function harmonic_average(gi, gj, interpolate) result(gij)
   !> Parameter i
   real(wp), intent(in) :: gi
   !> Parameter j
   real(wp), intent(in) :: gj
   !> Parameter to smoothely interpolate between different averages (not used)
   real(wp), intent(in), optional :: interpolate
   !> Averaged parameter
   real(wp) :: gij

   gij = 2.0_wp/(1.0_wp/gi+1.0_wp/gj)

end function harmonic_average

!> Harmonic averaging function derivative
pure function harmonic_average_deriv(gi, gj, interpolate) result(dgij)
   !> Parameter i
   real(wp), intent(in) :: gi
   !> Parameter j
   real(wp), intent(in) :: gj
   !> Parameter to smoothely interpolate between different averages 
   real(wp), intent(in), optional :: interpolate
   !> Derivative of the averaged parameter with respect to gi
   real(wp) :: dgij

   dgij = 2.0_wp*gj**2/(gi+gj)**2

end function harmonic_average_deriv


!> Interpolated general averaging function between
!> (0) arithmetic (default), (1) geometric, (2) harmonic
pure function general_average(gi, gj, interpolate) result(gij)
   !> Parameter i
   real(wp), intent(in) :: gi
   !> Parameter j
   real(wp), intent(in) :: gj
   !> Parameter to smoothely interpolate between different averages
   real(wp), intent(in), optional :: interpolate
   !> Averaged parameter
   real(wp) :: gij

   real(wp) :: tmp, alpha

   alpha = merge(interpolate, 0.0_wp, present(interpolate))
   tmp = 2.0_wp/(gi+gj)
   gij = tmp**(alpha-1.0_wp)*(gi*gj)**(alpha/2.0_wp)

end function general_average

!> Interpolated general averaging function derivative 
pure function general_average_deriv(gi, gj, interpolate) result(dgij)
   !> Parameter i
   real(wp), intent(in) :: gi
   !> Parameter j
   real(wp), intent(in) :: gj
   !> Parameter to smoothely interpolate between different averages 
   real(wp), intent(in), optional :: interpolate
   !> Derivative of the averaged parameter with respect to gi
   real(wp) :: dgij 

   real(wp) :: alpha, gij, tmp

   alpha = merge(interpolate, 0.0_wp, present(interpolate))
   tmp = 2.0_wp/(gi+gj)
   gij = tmp**(alpha-1.0_wp)*(gi*gj)**(alpha/2.0_wp)
   dgij = gij * ((alpha/2.0_wp)/gi - (alpha-1.0_wp)/(gi+gj))

end function general_average_deriv

end module tblite_utils_average