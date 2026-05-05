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

!> @file tblite/coulomb/ewald.f90
!> Provides an utilities for implementing Ewald summations

!> Helper tools for dealing with Ewald summation related calculations
module tblite_coulomb_ewald
   use mctc_env, only : wp
   use mctc_io_constants, only : pi
   use mctc_io_math, only : matdet_3x3, matinv_3x3
   implicit none
   private

   public :: get_alpha, get_dir_cutoff, get_rec_cutoff

   real(wp), parameter :: twopi = 2 * pi
   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))

   !> Evaluator for interaction term
   type, abstract :: term_type
   contains
      !> Interaction at specified distance
      procedure(get_value), deferred :: get_value
   end type term_type

   abstract interface
      !> Interaction at specified distance
      pure function get_value(self, dist, alpha, vol) result(val)
         import :: term_type, wp
         !> Instance of interaction
         class(term_type), intent(in) :: self
         !> Distance between the two atoms
         real(wp), intent(in) :: dist
         !> Parameter of the Ewald summation
         real(wp), intent(in) :: alpha
         !> Volume of the real space unit cell
         real(wp), intent(in) :: vol
         !> Value of the interaction
         real(wp) :: val
      end function get_value
   end interface

   !> Real space Coulombic interaction 1/R
   type, extends(term_type) :: dir_term
   contains
      procedure :: get_value => dir_value
   end type dir_term

   !> Reciprocal space Coulombic interaction 1/R in 3D
   type, extends(term_type) :: rec_3d_term
   contains
      procedure :: get_value => rec_3d_value
   end type rec_3d_term

   !> Real space Coulombic interaction 1/R^2
   type, extends(term_type) :: dir_mp_term
   contains
      procedure :: get_value => dir_mp_value
   end type dir_mp_term

   !> Reciprocal space Coulombic interaction 1/R^2 in 3D
   type, extends(term_type) :: rec_3d_mp_term
   contains
      procedure :: get_value => rec_3d_mp_value
   end type rec_3d_mp_term

contains


!> Convenience interface to determine Ewald splitting parameter
subroutine get_alpha(lattice, alpha, multipole)
   !> Lattice vectors
   real(wp), intent(in) :: lattice(:, :)
   !> Estimated Ewald splitting parameter
   real(wp), intent(out) :: alpha
   !> Multipole expansion is used
   logical, intent(in) :: multipole

   real(wp) :: vol, rec_lat(3, 3)
   class(term_type), allocatable :: dirv, recv

   vol = abs(matdet_3x3(lattice))
   rec_lat = twopi*transpose(matinv_3x3(lattice))
   if (multipole) then
      dirv = dir_mp_term()
      recv = rec_3d_mp_term()
   else
      dirv = dir_term()
      recv = rec_3d_term()
   end if

   call search_alpha(lattice, rec_lat, vol, eps, dirv, recv, alpha)
end subroutine get_alpha


!> Get optimal alpha-parameter for the Ewald summation by finding alpha, where
!> decline of real and reciprocal part of Ewald are equal.
subroutine search_alpha(lattice, rec_lat, volume, tolerance, dirv, recv, alpha)
   !> Lattice vectors
   real(wp), intent(in) :: lattice(:,:)
   !> Reciprocal vectors
   real(wp), intent(in) :: rec_lat(:,:)
   !> Volume of the unit cell
   real(wp), intent(in) :: volume
   !> Tolerance for difference in real and rec. part
   real(wp), intent(in) :: tolerance
   !> Real-space interaction term
   class(term_type), intent(in) :: dirv
   !> Reciprocal space interaction term
   class(term_type), intent(in) :: recv
   !> Optimal alpha
   real(wp), intent(out) :: alpha

   real(wp) :: alpl, alpr, rlen, dlen, diff
   real(wp), parameter :: alpha0 = sqrt(epsilon(0.0_wp))
   integer, parameter :: niter = 30
   integer :: ibs, stat

   rlen = sqrt(minval(sum(rec_lat(:,:)**2, dim=1)))
   dlen = sqrt(minval(sum(lattice(:,:)**2, dim=1)))

   stat = 0
   alpha = alpha0
   diff = rec_dir_diff(alpha, dirv, recv, rlen, dlen, volume)
   do while (diff < -tolerance .and. alpha <= huge(1.0_wp))
      alpha = 2.0_wp * alpha
      diff = rec_dir_diff(alpha, dirv, recv, rlen, dlen, volume)
   end do
   if (alpha > huge(1.0_wp)) then
      stat = 1
   elseif (alpha == alpha0) then
      stat = 2
   end if

   if (stat == 0) then
      alpl = 0.5_wp * alpha
      do while (diff < tolerance .and. alpha <= huge(1.0_wp))
         alpha = 2.0_wp * alpha
         diff = rec_dir_diff(alpha, dirv, recv, rlen, dlen, volume)
      end do
      if (alpha > huge(1.0_wp)) then
         stat = 3
      end if
   end if

   if (stat == 0) then
      alpr = alpha
      alpha = (alpl + alpr) * 0.5_wp
      ibs = 0
      diff = rec_dir_diff(alpha, dirv, recv, rlen, dlen, volume)
      do while (abs(diff) > tolerance .and. ibs <= niter)
         if (diff < 0) then
            alpl = alpha
         else
            alpr = alpha
         end if
         alpha = (alpl + alpr) * 0.5_wp
         diff = rec_dir_diff(alpha, dirv, recv, rlen, dlen, volume)
         ibs = ibs + 1
      end do
      if (ibs > niter) then
         stat = 4
      end if
   end if

   if (stat /= 0) then
      alpha = 0.25_wp
   end if

end subroutine search_alpha


!> Return cutoff for reciprocal contributions
function get_rec_cutoff(alpha, volume, conv) result(x)
   !> Parameter of Ewald summation
   real(wp), intent(in) :: alpha
   !> Volume of the unit cell
   real(wp), intent(in) :: volume
   !> Tolerance value
   real(wp), intent(in) :: conv
   !> Magnitude of reciprocal vector
   real(wp) :: x

   class(term_type), allocatable :: term

   term = rec_3d_term()
   x = search_cutoff(term, alpha, volume, conv)

end function get_rec_cutoff


!> Return cutoff for real-space contributions
function get_dir_cutoff(alpha, conv) result(x)
   !> Parameter of Ewald summation
   real(wp), intent(in) :: alpha
   !> Tolerance for Ewald summation
   real(wp), intent(in) :: conv
   !> Magnitude of real-space vector
   real(wp) :: x

   real(wp), parameter :: volume = 0.0_wp
   class(term_type), allocatable :: term

   term = dir_term()
   x = search_cutoff(term, alpha, volume, conv)

end function get_dir_cutoff


!> Search for cutoff value of interaction term
function search_cutoff(term, alpha, volume, conv) result(x)
   !> Interaction term
   class(term_type), intent(in) :: term
   !> Parameter of Ewald summation
   real(wp), intent(in) :: alpha
   !> Volume of the unit cell
   real(wp), intent(in) :: volume
   !> Tolerance value
   real(wp), intent(in) :: conv
   !> Magnitude of reciprocal vector
   real(wp) :: x

   real(wp), parameter :: init = sqrt(epsilon(0.0_wp))
   integer, parameter :: miter = 30
   real(wp) :: xl, xr, yl, yr, y
   integer :: iter

   x = init
   y = term%get_value(x, alpha, volume)
   do while (y > conv .and. x <= huge(1.0_wp))
      x = 2.0_wp * x
      y = term%get_value(x, alpha, volume)
   end do

   xl = 0.5_wp * x
   yl = term%get_value(xl, alpha, volume)
   xr = x
   yr = y

   do iter = 1, miter
      if (yl - yr <= conv) exit
      x = 0.5_wp * (xl + xr)
      y = term%get_value(x, alpha, volume)
      if (y >= conv) then
         xl = x
         yl = y
      else
         xr = x
         yr = y
      end if
   end do

end function search_cutoff


!> Returns the difference in the decrease of the real and reciprocal parts of the
!> Ewald sum. In order to make the real space part shorter than the reciprocal
!> space part, the values are taken at different distances for the real and the
!> reciprocal space parts.
pure function rec_dir_diff(alpha, dirv, recv, rlen, dlen, volume) result(diff)
   !> Parameter for the Ewald summation
   real(wp), intent(in) :: alpha
   !> Procedure pointer to real-space routine
   class(term_type), intent(in) :: dirv
   !> Procedure pointer to reciprocal routine
   class(term_type), intent(in) :: recv
   !> Length of the shortest reciprocal space vector in the sum
   real(wp), intent(in) :: rlen
   !> Length of the shortest real space vector in the sum
   real(wp), intent(in) :: dlen
   !> Volume of the real space unit cell
   real(wp), intent(in) :: volume
   !> Difference between changes in the two terms
   real(wp) :: diff

   diff = ((recv%get_value(4*rlen, alpha, volume) - recv%get_value(5*rlen, alpha, volume))) &
      & - (dirv%get_value(2*dlen, alpha, volume) - dirv%get_value(3*dlen, alpha, volume))

end function rec_dir_diff


!> Returns the max. value of a term in the reciprocal space part of the Ewald
!> summation for a given vector length.
pure function rec_3d_value(self, dist, alpha, vol) result(rval)
   !> Instance of interaction
   class(rec_3d_term), intent(in) :: self
   !> Length of the reciprocal space vector
   real(wp), intent(in) :: dist
   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha
   !> Volume of the real space unit cell
   real(wp), intent(in) :: vol
   !> Reciprocal term
   real(wp) :: rval

   rval = 4.0_wp*pi*(exp(-0.25_wp*dist*dist/(alpha**2))/(vol*dist*dist))

end function rec_3d_value


!> Returns the max. value of a term in the reciprocal space part of the Ewald
!> summation for a given vector length.
pure function rec_3d_mp_value(self, dist, alpha, vol) result(rval)
   !> Instance of interaction
   class(rec_3d_mp_term), intent(in) :: self
   !> Length of the reciprocal space vector
   real(wp), intent(in) :: dist
   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha
   !> Volume of the real space unit cell
   real(wp), intent(in) :: vol
   !> Reciprocal term
   real(wp) :: rval

   real(wp) :: g2

   g2 = dist*dist

   rval = 4.0_wp*pi*(exp(-0.25_wp*g2/(alpha**2))/vol)

end function rec_3d_mp_value


!> Direct space interaction at specified distance
pure function dir_value(self, dist, alpha, vol) result(val)
   !> Instance of interaction
   class(dir_term), intent(in) :: self
   !> Distance between the two atoms
   real(wp), intent(in) :: dist
   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha
   !> Volume of the real space unit cell
   real(wp), intent(in) :: vol
   !> Value of the interaction
   real(wp) :: val

   val = erfc(alpha*dist)/dist
end function dir_value

!> Direct space interaction at specified distance
pure function dir_mp_value(self, dist, alpha, vol) result(val)
   !> Instance of interaction
   class(dir_mp_term), intent(in) :: self
   !> Distance between the two atoms
   real(wp), intent(in) :: dist
   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha
   !> Volume of the real space unit cell
   real(wp), intent(in) :: vol
   !> Value of the interaction
   real(wp) :: val

   real(wp) :: arg

   arg = alpha * dist

   val = (2/sqrtpi*arg*exp(-arg*arg) + erfc(arg))/dist**3
end function dir_mp_value

end module tblite_coulomb_ewald
