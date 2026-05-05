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

!> @file tblite/lapack/getrf.f90
!> Provides wrappers for LU factorization routines

!> Wrapper rountines for LU factorization
module tblite_lapack_getrf
   use mctc_env, only : sp, dp
   implicit none
   private

   public :: wrap_getrf


   !> Computes an LU factorization of a general M-by-N matrix A
   !> using partial pivoting with row interchanges.
   !>
   !> The factorization has the form
   !>    A = P * L * U
   !> where P is a permutation matrix, L is lower triangular with unit
   !> diagonal elements (lower trapezoidal if m > n), and U is upper
   !> triangular (upper trapezoidal if m < n).
   interface wrap_getrf
      module procedure :: wrap_sgetrf
      module procedure :: wrap_dgetrf
   end interface wrap_getrf


   !> Computes an LU factorization of a general M-by-N matrix A
   !> using partial pivoting with row interchanges.
   !>
   !> The factorization has the form
   !>    A = P * L * U
   !> where P is a permutation matrix, L is lower triangular with unit
   !> diagonal elements (lower trapezoidal if m > n), and U is upper
   !> triangular (upper trapezoidal if m < n).
   interface lapack_getrf
      pure subroutine sgetrf(m, n, a, lda, ipiv, info)
         import :: sp
         real(sp), intent(inout) :: a(lda, *)
         integer, intent(out) :: ipiv(*)
         integer, intent(out) :: info
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine sgetrf
      pure subroutine dgetrf(m, n, a, lda, ipiv, info)
         import :: dp
         real(dp), intent(inout) :: a(lda, *)
         integer, intent(out) :: ipiv(*)
         integer, intent(out) :: info
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine dgetrf
   end interface lapack_getrf

contains

subroutine wrap_sgetrf(amat, ipiv, info)
   real(sp), intent(inout) :: amat(:, :)
   integer, intent(out) :: ipiv(:)
   integer, intent(out) :: info
   integer :: m, n, lda
   lda = max(1, size(amat, 1))
   m = size(amat, 1)
   n = size(amat, 2)
   call lapack_getrf(m, n, amat, lda, ipiv, info)
end subroutine wrap_sgetrf


subroutine wrap_dgetrf(amat, ipiv, info)
   real(dp), intent(inout) :: amat(:, :)
   integer, intent(out) :: ipiv(:)
   integer, intent(out) :: info
   integer :: m, n, lda
   lda = max(1, size(amat, 1))
   m = size(amat, 1)
   n = size(amat, 2)
   call lapack_getrf(m, n, amat, lda, ipiv, info)
end subroutine wrap_dgetrf

end module tblite_lapack_getrf
