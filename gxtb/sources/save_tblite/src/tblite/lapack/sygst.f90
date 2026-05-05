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

!> Reduces a real symmetric-definite generalized eigenproblem to standard form.
module tblite_lapack_sygst
   use mctc_env, only : sp, dp
   implicit none
   private

   public :: wrap_sygst


   !> Reduces a real symmetric-definite generalized eigenproblem to standard form.
   !>
   !> If ITYPE = 1, the problem is A*x = lambda*B*x,
   !> and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)
   !>
   !> If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
   !> B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T*A*L.
   !>
   !> B must have been previously factorized as U**T*U or L*L**T by POTRF.
   interface wrap_sygst
      module procedure :: wrap_ssygst
      module procedure :: wrap_dsygst
   end interface wrap_sygst


   !> Reduces a real symmetric-definite generalized eigenproblem to standard form.
   !>
   !> If ITYPE = 1, the problem is A*x = lambda*B*x,
   !> and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)
   !>
   !> If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
   !> B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T*A*L.
   !>
   !> B must have been previously factorized as U**T*U or L*L**T by POTRF.
   interface lapack_sygst
      pure subroutine ssygst(itype, uplo, n, a, lda, b, ldb, info)
         import :: sp
         real(sp), intent(inout) :: a(lda, *)
         real(sp), intent(in) :: b(ldb, *)
         integer, intent(in) :: itype
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
      end subroutine ssygst
      pure subroutine dsygst(itype, uplo, n, a, lda, b, ldb, info)
         import :: dp
         real(dp), intent(inout) :: a(lda, *)
         real(dp), intent(in) :: b(ldb, *)
         integer, intent(in) :: itype
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
      end subroutine dsygst
   end interface lapack_sygst

contains

pure subroutine wrap_ssygst(amat, bmat, info, itype, uplo)
   real(sp), intent(inout) :: amat(:, :)
   real(sp), intent(in) :: bmat(:, :)
   integer, intent(in), optional :: itype
   character(len=1), intent(in), optional :: uplo
   integer, intent(out) :: info
   character(len=1) :: ula
   integer :: ita, n, lda, ldb

   ita = 1
   if(present(itype)) ita = itype
   ula = 'u'
   if(present(uplo)) ula = uplo
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   n = size(amat, 2)
   call lapack_sygst(ita, ula, n, amat, lda, bmat, ldb, info)
end subroutine wrap_ssygst

pure subroutine wrap_dsygst(amat, bmat, info, itype, uplo)
   real(dp), intent(inout) :: amat(:, :)
   real(dp), intent(in) :: bmat(:, :)
   integer, intent(in), optional :: itype
   character(len=1), intent(in), optional :: uplo
   integer, intent(out) :: info
   character(len=1) :: ula
   integer :: ita, n, lda, ldb

   ita = 1
   if(present(itype)) ita = itype
   ula = 'u'
   if(present(uplo)) ula = uplo
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   n = size(amat, 2)
   call lapack_sygst(ita, ula, n, amat, lda, bmat, ldb, info)
end subroutine wrap_dsygst

end module tblite_lapack_sygst
