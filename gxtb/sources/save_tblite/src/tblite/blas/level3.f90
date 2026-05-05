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

!> @file tblite/blas/level3.f90
!> Provides interfactes to level 3 BLAS routines

!> High-level interface to level 3 basic linear algebra subprogram operations
module tblite_blas_level3
   use mctc_env, only : sp, dp
   implicit none
   private

   public :: wrap_gemm, wrap_trsm, wrap_symm


   !> Performs one of the matrix-matrix operations
   !>
   !>    C := alpha*op( A )*op( B ) + beta*C,
   !>
   !> where  op( X ) is one of
   !>
   !>    op( X ) = X   or   op( X ) = X**T,
   !>
   !> where alpha and beta are scalars, and A, B and C are matrices
   !> with op( A ) an m by k matrix, op( B ) a k by n matrix and C an m by n matrix.
   interface wrap_gemm
      module procedure :: wrap_sgemm
      module procedure :: wrap_dgemm
      module procedure :: wrap_sgemm323
      module procedure :: wrap_sgemm233
      module procedure :: wrap_sgemm332
      module procedure :: wrap_dgemm323
      module procedure :: wrap_dgemm233
      module procedure :: wrap_dgemm332
   end interface wrap_gemm

   !> Solves one of the matrix equations
   !>
   !>    op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
   !>
   !> where alpha is a scalar, X and B are m by n matrices, A is a unit, or
   !> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
   !>
   !>    op( A ) = A   or   op( A ) = A**T.
   !>
   !> The matrix X is overwritten on B.
   interface wrap_trsm
      module procedure :: wrap_strsm
      module procedure :: wrap_dtrsm
   end interface wrap_trsm

   !> Performs one of the symmetric matrix-matrix operations
   !>
   !>    C := alpha*A*B + beta*C,
   !>
   !> or
   !>
   !>    C := alpha*B*A + beta*C,
   !>
   !> where alpha and beta are scalars, A is a symmetric matrix 
   !> and B and C are m by n matrices.
   interface wrap_symm
      module procedure :: wrap_ssymm
      module procedure :: wrap_dsymm
   end interface wrap_symm

   !> Performs one of the matrix-matrix operations
   !>
   !>    C := alpha*op( A )*op( B ) + beta*C,
   !>
   !> where  op( X ) is one of
   !>
   !>    op( X ) = X   or   op( X ) = X**T,
   !>
   !> alpha and beta are scalars, and A, B and C are matrices, with op( A )
   !> an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
   interface blas_gemm
      pure subroutine sgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, &
            & beta, c, ldc)
         import :: sp
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(in) :: b(ldb, *)
         real(sp), intent(inout) :: c(ldc, *)
         character(len=1), intent(in) :: transa
         character(len=1), intent(in) :: transb
         real(sp), intent(in) :: alpha
         real(sp), intent(in) :: beta
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: k
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
         integer, intent(in) :: ldc
      end subroutine sgemm
      pure subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, &
            & beta, c, ldc)
         import :: dp
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(in) :: b(ldb, *)
         real(dp), intent(inout) :: c(ldc, *)
         character(len=1), intent(in) :: transa
         character(len=1), intent(in) :: transb
         real(dp), intent(in) :: alpha
         real(dp), intent(in) :: beta
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: k
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
         integer, intent(in) :: ldc
      end subroutine dgemm
   end interface blas_gemm

   !> Solves one of the matrix equations
   !>
   !>    op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
   !>
   !> where alpha is a scalar, X and B are m by n matrices, A is a unit, or
   !> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
   !>
   !>    op( A ) = A   or   op( A ) = A**T.
   !>
   !> The matrix X is overwritten on B.
   interface blas_trsm
      pure subroutine strsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
         import :: sp
         integer, intent(in) :: ldb
         integer, intent(in) :: lda
         character(len=1), intent(in) :: side
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: transa
         character(len=1), intent(in) :: diag
         integer, intent(in) :: m
         integer, intent(in) :: n
         real(sp), intent(in) :: alpha
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(inout) :: b(ldb, *)
      end subroutine strsm
      pure subroutine dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
         import :: dp
         integer, intent(in) :: ldb
         integer, intent(in) :: lda
         character(len=1), intent(in) :: side
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: transa
         character(len=1), intent(in) :: diag
         integer, intent(in) :: m
         integer, intent(in) :: n
         real(dp), intent(in) :: alpha
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(inout) :: b(ldb, *)
      end subroutine dtrsm
   end interface blas_trsm

   !> Performs one of the symmetric matrix-matrix operations
   !>
   !>    C := alpha*A*B + beta*C,
   !>
   !> or
   !>
   !>    C := alpha*B*A + beta*C,
   !>
   !> where alpha and beta are scalars, A is a symmetric matrix 
   !> and B and C are m by n matrices.
   interface blas_symm
      pure subroutine ssymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
         import :: sp
         character(len=1), intent(in) :: side, uplo
         integer, intent(in) :: m, n, lda, ldb, ldc
         real(sp), intent(in) :: alpha, beta
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(in) :: b(ldb, *)
         real(sp), intent(inout) :: c(ldc, *)
      end subroutine ssymm
      pure subroutine dsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
         import :: dp
         character(len=1), intent(in) :: side, uplo
         integer, intent(in) :: m, n, lda, ldb, ldc
         real(dp), intent(in) :: alpha, beta
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(in) :: b(ldb, *)
         real(dp), intent(inout) :: c(ldc, *)
      end subroutine dsymm
   end interface blas_symm

contains


pure subroutine wrap_sgemm(amat, bmat, cmat, transa, transb, alpha, beta)
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(in) :: bmat(:, :)
   real(sp), intent(inout) :: cmat(:, :)
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: transb
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   character(len=1) :: tra, trb
   real(sp) :: a, b
   integer :: m, n, k, lda, ldb, ldc
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_sp
   end if
   if (present(beta)) then
      b = beta
   else
      b = 0.0_sp
   end if
   if (present(transa)) then
      tra = transa
   else
      tra = 'n'
   end if
   if (present(transb)) then
      trb = transb
   else
      trb = 'n'
   end if
   if ((tra.eq.'n'.or.tra.eq.'N')) then
      k = size(amat, 2)
   else
      k = size(amat, 1)
   end if
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   ldc = max(1, size(cmat, 1))
   m = size(cmat, 1)
   n = size(cmat, 2)
   call blas_gemm(tra, trb, m, n, k, a, amat, lda, bmat, ldb, b, cmat, ldc)
end subroutine wrap_sgemm


pure subroutine wrap_dgemm(amat, bmat, cmat, transa, transb, alpha, beta)
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(in) :: bmat(:, :)
   real(dp), intent(inout) :: cmat(:, :)
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: transb
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   character(len=1) :: tra, trb
   real(dp) :: a, b
   integer :: m, n, k, lda, ldb, ldc
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_dp
   end if
   if (present(beta)) then
      b = beta
   else
      b = 0.0_dp
   end if
   if (present(transa)) then
      tra = transa
   else
      tra = 'n'
   end if
   if (present(transb)) then
      trb = transb
   else
      trb = 'n'
   end if
   if ((tra.eq.'n'.or.tra.eq.'N')) then
      k = size(amat, 2)
   else
      k = size(amat, 1)
   end if
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   ldc = max(1, size(cmat, 1))
   m = size(cmat, 1)
   n = size(cmat, 2)
   call blas_gemm(tra, trb, m, n, k, a, amat, lda, bmat, ldb, b, cmat, ldc)
end subroutine wrap_dgemm


subroutine wrap_sgemm323(amat, bmat, cmat, transa, transb, alpha, beta)
   real(sp), intent(in), contiguous, target :: amat(:, :, :)
   real(sp), intent(in) :: bmat(:, :)
   real(sp), intent(inout), contiguous, target :: cmat(:, :, :)
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: transb
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   real(sp), pointer :: aptr(:, :), cptr(:, :)
   character(len=1) :: tra
   if (present(transa)) then
      tra = transa
   else
      tra = 'n'
   end if
   if (any(tra == ['n', 'N'])) then
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
   else
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
   end if
   cptr(1:size(cmat, 1)*size(cmat, 2), 1:size(cmat, 3)) => cmat
   call wrap_gemm(aptr, bmat, cptr, tra, transb, alpha, beta)
end subroutine wrap_sgemm323


subroutine wrap_sgemm233(amat, bmat, cmat, transa, transb, alpha, beta)
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(in), contiguous, target :: bmat(:, :, :)
   real(sp), intent(inout), contiguous, target :: cmat(:, :, :)
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: transb
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   real(sp), pointer :: bptr(:, :), cptr(:, :)
   character(len=1) :: trb
   if (present(transb)) then
      trb = transb
   else
      trb = 'n'
   end if
   if (any(trb == ['n', 'N'])) then
      bptr(1:size(bmat, 1), 1:size(bmat, 2)*size(bmat, 3)) => bmat
   else
      bptr(1:size(bmat, 1)*size(bmat, 2), 1:size(bmat, 3)) => bmat
   end if
   cptr(1:size(cmat, 1), 1:size(cmat, 2)*size(cmat, 3)) => cmat
   call wrap_gemm(amat, bptr, cptr, transa, trb, alpha, beta)
end subroutine wrap_sgemm233


subroutine wrap_sgemm332(amat, bmat, cmat, transa, transb, alpha, beta)
   real(sp), intent(in), contiguous, target :: amat(:, :, :)
   real(sp), intent(in), contiguous, target :: bmat(:, :, :)
   real(sp), intent(inout) :: cmat(:, :)
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: transb
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   real(sp), pointer :: aptr(:, :), bptr(:, :)
   character(len=1) :: tra, trb
   if (present(transa)) then
      tra = transa
   else
      tra = 'n'
   end if
   if (present(transb)) then
      trb = transb
   else
      trb = 'n'
   end if
   if (any(tra == ['n', 'N'])) then
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
   else
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
   end if
   if (any(trb == ['n', 'N'])) then
      bptr(1:size(bmat, 1)*size(bmat, 2), 1:size(bmat, 3)) => bmat
   else
      bptr(1:size(bmat, 1), 1:size(bmat, 2)*size(bmat, 3)) => bmat
   end if
   call wrap_gemm(aptr, bptr, cmat, tra, trb, alpha, beta)
end subroutine wrap_sgemm332


subroutine wrap_dgemm323(amat, bmat, cmat, transa, transb, alpha, beta)
   real(dp), intent(in), contiguous, target :: amat(:, :, :)
   real(dp), intent(in) :: bmat(:, :)
   real(dp), intent(inout), contiguous, target :: cmat(:, :, :)
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: transb
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   real(dp), pointer :: aptr(:, :), cptr(:, :)
   character(len=1) :: tra
   if (present(transa)) then
      tra = transa
   else
      tra = 'n'
   end if
   if (any(tra == ['n', 'N'])) then
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
   else
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
   end if
   cptr(1:size(cmat, 1)*size(cmat, 2), 1:size(cmat, 3)) => cmat
   call wrap_gemm(aptr, bmat, cptr, tra, transb, alpha, beta)
end subroutine wrap_dgemm323


subroutine wrap_dgemm233(amat, bmat, cmat, transa, transb, alpha, beta)
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(in), contiguous, target :: bmat(:, :, :)
   real(dp), intent(inout), contiguous, target :: cmat(:, :, :)
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: transb
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   real(dp), pointer :: bptr(:, :), cptr(:, :)
   character(len=1) :: trb
   if (present(transb)) then
      trb = transb
   else
      trb = 'n'
   end if
   if (any(trb == ['n', 'N'])) then
      bptr(1:size(bmat, 1), 1:size(bmat, 2)*size(bmat, 3)) => bmat
   else
      bptr(1:size(bmat, 1)*size(bmat, 2), 1:size(bmat, 3)) => bmat
   end if
   cptr(1:size(cmat, 1), 1:size(cmat, 2)*size(cmat, 3)) => cmat
   call wrap_gemm(amat, bptr, cptr, transa, trb, alpha, beta)
end subroutine wrap_dgemm233


subroutine wrap_dgemm332(amat, bmat, cmat, transa, transb, alpha, beta)
   real(dp), intent(in), contiguous, target :: amat(:, :, :)
   real(dp), intent(in), contiguous, target :: bmat(:, :, :)
   real(dp), intent(inout) :: cmat(:, :)
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: transb
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   real(dp), pointer :: aptr(:, :), bptr(:, :)
   character(len=1) :: tra, trb
   if (present(transa)) then
      tra = transa
   else
      tra = 'n'
   end if
   if (present(transb)) then
      trb = transb
   else
      trb = 'n'
   end if
   if (any(tra == ['n', 'N'])) then
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
   else
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
   end if
   if (any(trb == ['n', 'N'])) then
      bptr(1:size(bmat, 1)*size(bmat, 2), 1:size(bmat, 3)) => bmat
   else
      bptr(1:size(bmat, 1), 1:size(bmat, 2)*size(bmat, 3)) => bmat
   end if
   call wrap_gemm(aptr, bptr, cmat, tra, trb, alpha, beta)
end subroutine wrap_dgemm332


pure subroutine wrap_strsm(amat, bmat, side, uplo, transa, diag, alpha)
   real(sp), contiguous, intent(in) :: amat(:,  :)
   real(sp), contiguous, intent(inout) :: bmat(:, :)
   real(sp), intent(in), optional :: alpha
   character(len=1), intent(in), optional :: side
   character(len=1), intent(in), optional :: uplo
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: diag
   real(sp) :: a
   character(len=1) :: sda, ula, tra, dga
   integer :: m, n, lda, ldb

   a = 1.0_sp
   if (present(alpha)) a = alpha
   dga = 'n'
   if (present(diag)) dga = diag
   sda = 'l'
   if (present(side)) sda = side
   tra = 'n'
   if (present(transa)) tra = transa
   ula = 'u'
   if (present(uplo)) ula = uplo
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   m = size(bmat, 1)
   n = size(bmat, 2)
   call blas_trsm(sda, ula, tra, dga, m, n, a, amat, lda, bmat, ldb)
end subroutine wrap_strsm


pure subroutine wrap_dtrsm(amat, bmat, side, uplo, transa, diag, alpha)
   real(dp), contiguous, intent(in) :: amat(:,  :)
   real(dp), contiguous, intent(inout) :: bmat(:, :)
   real(dp), intent(in), optional :: alpha
   character(len=1), intent(in), optional :: side
   character(len=1), intent(in), optional :: uplo
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: diag
   real(dp) :: a
   character(len=1) :: sda, ula, tra, dga
   integer :: m, n, lda, ldb

   a = 1.0_dp
   if (present(alpha)) a = alpha
   dga = 'n'
   if (present(diag)) dga = diag
   sda = 'l'
   if (present(side)) sda = side
   tra = 'n'
   if (present(transa)) tra = transa
   ula = 'u'
   if (present(uplo)) ula = uplo
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   m = size(bmat, 1)
   n = size(bmat, 2)
   call blas_trsm(sda, ula, tra, dga, m, n, a, amat, lda, bmat, ldb)
end subroutine wrap_dtrsm


pure subroutine wrap_ssymm(amat, bmat, cmat, side, uplo, alpha, beta)
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(in) :: bmat(:, :)
   real(sp), intent(inout) :: cmat(:, :)
   character(len=1), intent(in), optional :: side
   character(len=1), intent(in), optional :: uplo
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta

   character(len=1) :: sda, ula
   real(sp) :: a, b
   integer :: m, n, lda, ldb, ldc

   a = 1.0_sp
   if (present(alpha)) a = alpha
   b = 0.0_sp
   if (present(beta))  b = beta
   sda = 'l'
   if (present(side))  sda = side
   ula = 'u'
   if (present(uplo))  ula = uplo
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   ldc = max(1, size(cmat, 1))
   m   = size(cmat, 1)
   n   = size(cmat, 2)
   call blas_symm(sda, ula, m, n, a, amat, lda, bmat, ldb, b, cmat, ldc)
end subroutine wrap_ssymm


pure subroutine wrap_dsymm(amat, bmat, cmat, side, uplo, alpha, beta)
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(in) :: bmat(:, :)
   real(dp), intent(inout) :: cmat(:, :)
   character(len=1), intent(in), optional :: side
   character(len=1), intent(in), optional :: uplo
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta

   character(len=1) :: sda, ula
   real(dp) :: a, b
   integer :: m, n, lda, ldb, ldc

   a = 1.0_dp
   if (present(alpha)) a = alpha
   b = 0.0_dp
   if (present(beta))  b = beta
   sda = 'l'
   if (present(side))  sda = side
   ula = 'u'
   if (present(uplo))  ula = uplo
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   ldc = max(1, size(cmat, 1))
   m   = size(cmat, 1)
   n   = size(cmat, 2)
   call blas_symm(sda, ula, m, n, a, amat, lda, bmat, ldb, b, cmat, ldc)
end subroutine wrap_dsymm


end module tblite_blas_level3
