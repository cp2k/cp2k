! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later
!
! tblite is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! tblite is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! Lesser GNU General Public License for more details.
!
! You should have received a copy of the Lesser GNU General Public License
! along with tblite.  If not, see <https://www.gnu.org/licenses/>.

!> @file tblite/blas/level2.f90
!> Provides interfactes to level 2 BLAS routines

!> High-level interface to level 2 basic linear algebra subprogram operations
module tblite_blas_level2
   use mctc_env, only : sp, dp
   implicit none
   private

   public :: wrap_gemv, wrap_symv


   !> Performs one of the matrix-vector operations
   !>
   !>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
   !>
   !> where alpha and beta are scalars, x and y are vectors and A is an
   !> m by n matrix.
   interface wrap_gemv
      module procedure :: wrap_sgemv
      module procedure :: wrap_dgemv
      module procedure :: wrap_sgemv312
      module procedure :: wrap_sgemv321
      module procedure :: wrap_sgemv422
      module procedure :: wrap_dgemv312
      module procedure :: wrap_dgemv321
      module procedure :: wrap_dgemv422
   end interface wrap_gemv

   !> Performs the matrix-vector  operation
   !>
   !>    y := alpha*A*x + beta*y,
   !>
   !> where alpha and beta are scalars, x and y are n element vectors and
   !> A is an n by n symmetric matrix.
   interface wrap_symv
      module procedure :: wrap_ssymv
      module procedure :: wrap_dsymv
   end interface wrap_symv


   !> Performs one of the matrix-vector operations
   !>
   !>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
   !>
   !> where alpha and beta are scalars, x and y are vectors and A is an
   !> m by n matrix.
   interface blas_gemv
      pure subroutine sgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
         import :: sp
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(in) :: x(*)
         real(sp), intent(inout) :: y(*)
         real(sp), intent(in) :: alpha
         real(sp), intent(in) :: beta
         character(len=1), intent(in) :: trans
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine sgemv
      pure subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
         import :: dp
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(in) :: x(*)
         real(dp), intent(inout) :: y(*)
         real(dp), intent(in) :: alpha
         real(dp), intent(in) :: beta
         character(len=1), intent(in) :: trans
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine dgemv
   end interface blas_gemv

   !> Performs the matrix-vector  operation
   !>
   !>    y := alpha*A*x + beta*y,
   !>
   !> where alpha and beta are scalars, x and y are n element vectors and
   !> A is an n by n symmetric matrix.
   interface blas_symv
      pure subroutine ssymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
         import :: sp
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(in) :: x(*)
         real(sp), intent(inout) :: y(*)
         character(len=1), intent(in) :: uplo
         real(sp), intent(in) :: alpha
         real(sp), intent(in) :: beta
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine ssymv
      pure subroutine dsymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
         import :: dp
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(in) :: x(*)
         real(dp), intent(inout) :: y(*)
         character(len=1), intent(in) :: uplo
         real(dp), intent(in) :: alpha
         real(dp), intent(in) :: beta
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine dsymv
   end interface blas_symv


contains


subroutine wrap_sgemv422(amat, xvec, yvec, alpha, beta, trans)
   real(sp), intent(in), contiguous, target :: amat(:, :, :, :)
   real(sp), intent(in), contiguous, target :: xvec(:, :)
   real(sp), intent(inout), contiguous, target :: yvec(:, :)
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(sp), pointer :: aptr(:, :), yptr(:), xptr(:)
   character(len=1) :: tra
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)*size(amat, 4)) => amat
   xptr(1:size(xvec, 1)*size(xvec, 2)) => xvec
   yptr(1:size(yvec, 1)*size(yvec, 2)) => yvec
   call wrap_gemv(aptr, xptr, yptr, alpha, beta, tra)
end subroutine wrap_sgemv422


subroutine wrap_sgemv312(amat, xvec, yvec, alpha, beta, trans)
   real(sp), intent(in), contiguous, target :: amat(:, :, :)
   real(sp), intent(in) :: xvec(:)
   real(sp), intent(inout), contiguous, target :: yvec(:, :)
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(sp), pointer :: aptr(:, :), yptr(:)
   character(len=1) :: tra
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   if (any(tra == ['n', 'N'])) then
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
      yptr(1:size(yvec, 1)*size(yvec, 2)) => yvec
   else
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
      yptr(1:size(yvec, 1) * size(yvec, 2)) => yvec
   end if
   call wrap_gemv(aptr, xvec, yptr, alpha, beta, tra)
end subroutine wrap_sgemv312


subroutine wrap_sgemv321(amat, xvec, yvec, alpha, beta, trans)
   real(sp), intent(in), contiguous, target :: amat(:, :, :)
   real(sp), intent(in), contiguous, target :: xvec(:, :)
   real(sp), intent(inout) :: yvec(:)
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(sp), pointer :: aptr(:, :), xptr(:)
   character(len=1) :: tra
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   if (any(tra == ['n', 'N'])) then
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
      xptr(1:size(xvec, 1)*size(xvec, 2)) => xvec
   else
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
      xptr(1:size(xvec, 1) * size(xvec, 2)) => xvec
   end if
   call wrap_gemv(aptr, xptr, yvec, alpha, beta, tra)
end subroutine wrap_sgemv321


subroutine wrap_dgemv422(amat, xvec, yvec, alpha, beta, trans)
   real(dp), intent(in), contiguous, target :: amat(:, :, :, :)
   real(dp), intent(in), contiguous, target :: xvec(:, :)
   real(dp), intent(inout), contiguous, target :: yvec(:, :)
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(dp), pointer :: aptr(:, :), yptr(:), xptr(:)
   character(len=1) :: tra
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)*size(amat, 4)) => amat
   xptr(1:size(xvec, 1)*size(xvec, 2)) => xvec
   yptr(1:size(yvec, 1)*size(yvec, 2)) => yvec
   call wrap_gemv(aptr, xptr, yptr, alpha, beta, tra)
end subroutine wrap_dgemv422


subroutine wrap_dgemv312(amat, xvec, yvec, alpha, beta, trans)
   real(dp), intent(in), contiguous, target :: amat(:, :, :)
   real(dp), intent(in) :: xvec(:)
   real(dp), intent(inout), contiguous, target :: yvec(:, :)
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(dp), pointer :: aptr(:, :), yptr(:)
   character(len=1) :: tra
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   if (any(tra == ['n', 'N'])) then
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
      yptr(1:size(yvec, 1)*size(yvec, 2)) => yvec
   else
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
      yptr(1:size(yvec, 1) * size(yvec, 2)) => yvec
   end if
   call wrap_gemv(aptr, xvec, yptr, alpha, beta, tra)
end subroutine wrap_dgemv312


subroutine wrap_dgemv321(amat, xvec, yvec, alpha, beta, trans)
   real(dp), intent(in), contiguous, target :: amat(:, :, :)
   real(dp), intent(in), contiguous, target :: xvec(:, :)
   real(dp), intent(inout) :: yvec(:)
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(dp), pointer :: aptr(:, :), xptr(:)
   character(len=1) :: tra
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   if (any(tra == ['n', 'N'])) then
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
      xptr(1:size(xvec, 1)*size(xvec, 2)) => xvec
   else
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
      xptr(1:size(xvec, 1) * size(xvec, 2)) => xvec
   end if
   call wrap_gemv(aptr, xptr, yvec, alpha, beta, tra)
end subroutine wrap_dgemv321


pure subroutine wrap_sgemv(amat, xvec, yvec, alpha, beta, trans)
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(in) :: xvec(:)
   real(sp), intent(inout) :: yvec(:)
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(sp) :: a, b
   character(len=1) :: tra
   integer :: incx, incy, m, n, lda
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_sp
   end if
   if (present(beta)) then
      b = beta
   else
      b = 0
   end if
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   incx = 1
   incy = 1
   lda = max(1, size(amat, 1))
   m = size(amat, 1)
   n = size(amat, 2)
   call blas_gemv(tra, m, n, a, amat, lda, xvec, incx, b, yvec, incy)
end subroutine wrap_sgemv


pure subroutine wrap_dgemv(amat, xvec, yvec, alpha, beta, trans)
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(in) :: xvec(:)
   real(dp), intent(inout) :: yvec(:)
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(dp) :: a, b
   character(len=1) :: tra
   integer :: incx, incy, m, n, lda
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_dp
   end if
   if (present(beta)) then
      b = beta
   else
      b = 0
   end if
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   incx = 1
   incy = 1
   lda = max(1, size(amat, 1))
   m = size(amat, 1)
   n = size(amat, 2)
   call blas_gemv(tra, m, n, a, amat, lda, xvec, incx, b, yvec, incy)
end subroutine wrap_dgemv


pure subroutine wrap_ssymv(amat, xvec, yvec, uplo, alpha, beta)
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(in) :: xvec(:)
   real(sp), intent(inout) :: yvec(:)
   character(len=1), intent(in), optional :: uplo
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   character(len=1) :: ula
   real(sp) :: a, b
   integer :: incx, incy, n, lda
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_sp
   end if
   if (present(beta)) then
      b = beta
   else
      b = 0
   end if
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   incx = 1
   incy = 1
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   call blas_symv(ula, n, a, amat, lda, xvec, incx, b, yvec, incy)
end subroutine wrap_ssymv


pure subroutine wrap_dsymv(amat, xvec, yvec, uplo, alpha, beta)
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(in) :: xvec(:)
   real(dp), intent(inout) :: yvec(:)
   character(len=1), intent(in), optional :: uplo
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   character(len=1) :: ula
   real(dp) :: a, b
   integer :: incx, incy, n, lda
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_dp
   end if
   if (present(beta)) then
      b = beta
   else
      b = 0
   end if
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   incx = 1
   incy = 1
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   call blas_symv(ula, n, a, amat, lda, xvec, incx, b, yvec, incy)
end subroutine wrap_dsymv


end module tblite_blas_level2
