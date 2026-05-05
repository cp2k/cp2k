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

!> @file tblite/blas/level1.f90
!> Provides interfactes to level 1 BLAS routines

!> High-level interface to level 1 basic linear algebra subprogram operations
module tblite_blas_level1
   use mctc_env, only : sp, dp
   implicit none
   private

   public :: wrap_dot, wrap_axpy


   !> Forms the dot product of two vectors.
   interface wrap_dot
      module procedure :: wrap_sdot
      module procedure :: wrap_ddot
      module procedure :: wrap_sdot12
      module procedure :: wrap_sdot21
      module procedure :: wrap_sdot22
      module procedure :: wrap_ddot12
      module procedure :: wrap_ddot21
      module procedure :: wrap_ddot22
   end interface wrap_dot

   !> Forms constant times a vector plus a vector
   interface wrap_axpy
      module procedure :: wrap_saxpy
      module procedure :: wrap_daxpy
   end interface wrap_axpy

   !> Forms the dot product of two vectors.
   !> Uses unrolled loops for increments equal to one.
   interface blas_dot
      pure function sdot(n, x, incx, y, incy)
         import :: sp
         real(sp) :: sdot
         real(sp), intent(in) :: x(*)
         real(sp), intent(in) :: y(*)
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end function sdot
      pure function ddot(n, x, incx, y, incy)
         import :: dp
         real(dp) :: ddot
         real(dp), intent(in) :: x(*)
         real(dp), intent(in) :: y(*)
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end function ddot
   end interface blas_dot

   !> Forms constant times a vector plus a vector.
   !> Uses unrolled loops for increments equal to one.
   interface blas_axpy
      pure subroutine saxpy(n, a, x, incx, y, incy)
         import :: sp
         integer, intent(in) :: n
         real(sp), intent(in) :: a
         real(sp), intent(in) :: x(*)
         integer, intent(in) :: incx
         real(sp), intent(inout) :: y(*)
         integer, intent(in) :: incy
      end subroutine saxpy
      pure subroutine daxpy(n, a, x, incx, y, incy)
         import :: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: a
         real(dp), intent(in) :: x(*)
         integer, intent(in) :: incx
         real(dp), intent(inout) :: y(*)
         integer, intent(in) :: incy
      end subroutine daxpy
   end interface blas_axpy

contains


function wrap_sdot(xvec, yvec) result(dot)
   real(sp) :: dot
   real(sp), intent(in) :: xvec(:)
   real(sp), intent(in) :: yvec(:)
   integer :: incx, incy, n
   incx = 1
   incy = 1
   n = size(xvec)
   dot = blas_dot(n, xvec, incx, yvec, incy)
end function wrap_sdot


function wrap_ddot(xvec, yvec) result(dot)
   real(dp) :: dot
   real(dp), intent(in) :: xvec(:)
   real(dp), intent(in) :: yvec(:)
   integer :: incx, incy, n
   incx = 1
   incy = 1
   n = size(xvec)
   dot = blas_dot(n, xvec, incx, yvec, incy)
end function wrap_ddot


function wrap_sdot12(xvec, yvec) result(dot)
   real(sp) :: dot
   real(sp), intent(in) :: xvec(:)
   real(sp), intent(in), contiguous, target :: yvec(:, :)
   real(sp), pointer :: yptr(:)
   yptr(1:size(yvec)) => yvec
   dot = wrap_dot(xvec, yptr)
end function wrap_sdot12


function wrap_sdot21(xvec, yvec) result(dot)
   real(sp) :: dot
   real(sp), intent(in), contiguous, target :: xvec(:, :)
   real(sp), intent(in) :: yvec(:)
   real(sp), pointer :: xptr(:)
   xptr(1:size(xvec)) => xvec
   dot = wrap_dot(xptr, yvec)
end function wrap_sdot21


function wrap_sdot22(xvec, yvec) result(dot)
   real(sp) :: dot
   real(sp), intent(in), contiguous, target :: xvec(:, :)
   real(sp), intent(in), contiguous, target :: yvec(:, :)
   real(sp), pointer :: xptr(:), yptr(:)
   xptr(1:size(xvec)) => xvec
   yptr(1:size(yvec)) => yvec
   dot = wrap_dot(xptr, yptr)
end function wrap_sdot22


function wrap_ddot12(xvec, yvec) result(dot)
   real(dp) :: dot
   real(dp), intent(in) :: xvec(:)
   real(dp), intent(in), contiguous, target :: yvec(:, :)
   real(dp), pointer :: yptr(:)
   yptr(1:size(yvec)) => yvec
   dot = wrap_dot(xvec, yptr)
end function wrap_ddot12


function wrap_ddot21(xvec, yvec) result(dot)
   real(dp) :: dot
   real(dp), intent(in), contiguous, target :: xvec(:, :)
   real(dp), intent(in) :: yvec(:)
   real(dp), pointer :: xptr(:)
   xptr(1:size(xvec)) => xvec
   dot = wrap_dot(xptr, yvec)
end function wrap_ddot21


function wrap_ddot22(xvec, yvec) result(dot)
   real(dp) :: dot
   real(dp), intent(in), contiguous, target :: xvec(:, :)
   real(dp), intent(in), contiguous, target :: yvec(:, :)
   real(dp), pointer :: xptr(:), yptr(:)
   xptr(1:size(xvec)) => xvec
   yptr(1:size(yvec)) => yvec
   dot = wrap_dot(xptr, yptr)
end function wrap_ddot22


subroutine wrap_saxpy(xvec, yvec, alpha)
   real(sp), intent(in) :: xvec(:)
   real(sp), intent(inout) :: yvec(:)
   real(sp), intent(in), optional :: alpha
   
   integer :: incx, incy, n
   real(sp) :: a

   n = size(xvec)
   a = 1.0_sp
   incx = 1
   incy = 1
   if (present(alpha)) a = alpha
   call blas_axpy(n, a, xvec, incx, yvec, incy)
end subroutine wrap_saxpy


subroutine wrap_daxpy(xvec, yvec, alpha)
   real(dp), intent(in) :: xvec(:)
   real(dp), intent(inout) :: yvec(:)
   real(dp), intent(in), optional :: alpha

   integer :: incx, incy, n
   real(dp) :: a

   n = size(xvec)
   a = 1.0_dp
   incx = 1
   incy = 1
   if (present(alpha)) a = alpha
   call blas_axpy(n, a, xvec, incx, yvec, incy)
end subroutine wrap_daxpy


end module tblite_blas_level1
