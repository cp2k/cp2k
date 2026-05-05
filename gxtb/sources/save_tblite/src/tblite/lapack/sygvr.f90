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

module tblite_lapack_sygvr
   use mctc_env, only : sp, dp, wp, error_type, fatal_error
   use tblite_output_format, only : format_string
   use tblite_scf_diag, only : diag_solver_type
   use tblite_lapack_sygst, only : wrap_sygst
   use tblite_lapack_potrf, only : wrap_potrf
   use tblite_blas_level3, only : wrap_trsm
   implicit none
   private

   public :: new_sygvr


   !> Computes selected eigenvalues and, optionally, eigenvectors
   !> of a real symmetric matrix A.  Eigenvalues and eigenvectors can be
   !> selected by specifying either a range of values or a range of
   !> indices for the desired eigenvalues.
   !>
   !> The desired accuracy of the output can be specified by the input
   !> parameter ABSTOL.
   interface lapack_syevr
      pure subroutine ssyevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, &
            & z, ldz, isuppz, work, lwork, iwork, liwork, info)
         import :: sp
         real(sp), intent(inout) :: a(lda, *)
         real(sp), intent(out) :: w(*)
         character(len=1), intent(in) :: uplo
         real(sp), intent(out) :: z(ldz, *)
         real(sp), intent(in) :: vl
         real(sp), intent(in) :: vu
         integer, intent(in) :: il
         integer, intent(in) :: iu
         integer, intent(out) :: m
         integer, intent(out) :: isuppz(*)
         real(sp), intent(in) :: abstol
         integer, intent(out) :: info
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: range
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldz
         real(sp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         integer, intent(inout) :: iwork(*)
         integer, intent(in) :: liwork
      end subroutine ssyevr
      pure subroutine dsyevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, &
            & z, ldz, isuppz, work, lwork, iwork, liwork, info)
         import :: dp
         real(dp), intent(inout) :: a(lda, *)
         real(dp), intent(out) :: w(*)
         character(len=1), intent(in) :: uplo
         real(dp), intent(out) :: z(ldz, *)
         real(dp), intent(in) :: vl
         real(dp), intent(in) :: vu
         integer, intent(in) :: il
         integer, intent(in) :: iu
         integer, intent(out) :: m
         integer, intent(out) :: isuppz(*)
         real(dp), intent(in) :: abstol
         integer, intent(out) :: info
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: range
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldz
         real(dp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         integer, intent(inout) :: iwork(*)
         integer, intent(in) :: liwork
      end subroutine dsyevr
   end interface lapack_syevr


   type, public, extends(diag_solver_type) :: sygvr_solver
      private
      integer :: n = 0
      integer, allocatable :: iwork(:)
      integer, allocatable :: isuppz(:)
      real(sp), allocatable :: swork(:)
      real(sp), allocatable :: sbmat(:, :)
      real(sp), allocatable :: schole(:)
      real(dp), allocatable :: dwork(:)
      real(dp), allocatable :: dbmat(:, :)
      real(dp), allocatable :: dchole(:)
   contains
      procedure :: solve_sp
      procedure :: solve_dp
   end type sygvr_solver

contains

subroutine new_sygvr(self, overlap, nel, kt)
   type(sygvr_solver), intent(out) :: self
   real(wp), intent(in) :: overlap(:, :)
   real(wp), intent(in) :: nel(:)
   real(wp), intent(in) :: kt
   self%n = size(overlap, 1)
   self%nel = nel
   self%kt = kt
end subroutine new_sygvr

subroutine solve_sp(self, hmat, smat, eval, error)
   class(sygvr_solver), intent(inout) :: self
   real(sp), contiguous, intent(inout) :: hmat(:, :)
   real(sp), contiguous, intent(in) :: smat(:, :)
   real(sp), contiguous, intent(inout) :: eval(:)
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: itype = 1
   character(len=1), parameter :: uplo = 'U'

   logical :: upper
   integer :: info, lswork, liwork, m, ii, jj
   real(sp) :: abstol, vl, vu
   character(len=1) :: trans, uplo_new

   if (self%n == 0) then
      self%n = size(hmat, 1)
   end if
   if (.not.allocated(self%swork)) then
      lswork = query(self%n, uplo, 'S')
      allocate(self%swork(lswork))
   end if
   if (.not.allocated(self%iwork)) then
      liwork = max(1, 10 * self%n)
      allocate(self%iwork(liwork))
   end if
   if (.not.allocated(self%schole)) then
      allocate(self%schole(self%n))
   end if
   if (.not.allocated(self%isuppz)) then
      allocate(self%isuppz(2*self%n))
   end if
   self%sbmat = smat
   lswork = size(self%swork)
   liwork = size(self%iwork)
   upper = (uplo == 'U' .or. uplo == 'u')
   abstol = slamch_s()
   info = 0

   ! Form a Cholesky factorization of B.
   call wrap_potrf(self%sbmat, info, uplo=uplo)
   call handle_info(error, info)
   if (allocated(error)) return

   ! Transform problem to standard eigenvalue problem and solve.
   call wrap_sygst(hmat, self%sbmat, info, itype=itype, uplo=uplo)
   call handle_info(error, info)
   if (allocated(error)) return

   ! Save Cholesky factor in the other triangle of H and chole
   do ii = 1, self%n
      self%schole(ii) = self%sbmat(ii, ii)
   end do
   if (upper) then
      do jj = 1, self%n
         do ii = jj + 1, self%n
            hmat(ii, jj) = self%sbmat(jj, ii)
         end do
      end do
   else
      do jj = 1, self%n
         do ii = 1, jj - 1
            hmat(ii, jj) = self%sbmat(jj, ii)
         end do
      end do
   end if

   call lapack_syevr('V', 'A', uplo, self%n, hmat, self%n, vl, vu, 1, self%n, abstol, &
      & m, eval, self%sbmat, self%n, self%isuppz, self%swork, lswork, self%iwork, liwork, &
      & info)
   call handle_info(error, info)
   if (allocated(error)) return

   ! Backtransform eigenvectors to the original problem.
   do ii = 1, self%n
      hmat(ii, ii) = self%schole(ii)
   end do

   uplo_new = merge('L', 'U', upper)
   upper = .not.upper

   ! For A*x=(lambda)*B*x
   ! backtransform eigenvectors: x = inv(L)'*y or inv(U)*y'
   trans = merge('N', 'T', upper)
   call wrap_trsm(hmat, self%sbmat, side='L', uplo=uplo_new, transa=trans, diag='N')
   do ii = 1, m
      hmat(1:self%n, ii) = self%sbmat(1:self%n, ii)
   end do
   hmat(1:self%n, m+1:self%n) = 0.0_sp

end subroutine solve_sp

subroutine solve_dp(self, hmat, smat, eval, error)
   class(sygvr_solver), intent(inout) :: self
   real(dp), contiguous, intent(inout) :: hmat(:, :)
   real(dp), contiguous, intent(in) :: smat(:, :)
   real(dp), contiguous, intent(inout) :: eval(:)
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: itype = 1
   character(len=1), parameter :: uplo = 'U'

   logical :: upper
   integer :: info, ldwork, liwork, m, ii, jj
   real(dp) :: abstol, vl, vu
   character(len=1) :: trans, uplo_new

   if (self%n == 0) then
      self%n = size(hmat, 1)
   end if
   if (.not.allocated(self%dwork)) then
      ldwork = query(self%n, uplo, 'D')
      allocate(self%dwork(ldwork))
   end if
   if (.not.allocated(self%iwork)) then
      liwork = max(1, 10 * self%n)
      allocate(self%iwork(liwork))
   end if
   if (.not.allocated(self%dchole)) then
      allocate(self%dchole(self%n))
   end if
   if (.not.allocated(self%isuppz)) then
      allocate(self%isuppz(2*self%n))
   end if
   self%dbmat = smat
   ldwork = size(self%dwork)
   liwork = size(self%iwork)
   upper = (uplo == 'U' .or. uplo == 'u')
   abstol = dlamch_s()
   info = 0

   ! Form a Cholesky factorization of B.
   call wrap_potrf(self%dbmat, info, uplo=uplo)
   call handle_info(error, info)
   if (allocated(error)) return

   ! Transform problem to standard eigenvalue problem and solve.
   call wrap_sygst(hmat, self%dbmat, info, itype=itype, uplo=uplo)
   call handle_info(error, info)
   if (allocated(error)) return

   ! Save Cholesky factor in the other triangle of H and chole
   do ii = 1, self%n
      self%dchole(ii) = self%dbmat(ii, ii)
   end do
   if (upper) then
      do jj = 1, self%n
         do ii = jj + 1, self%n
            hmat(ii, jj) = self%dbmat(jj, ii)
         end do
      end do
   else
      do jj = 1, self%n
         do ii = 1, jj - 1
            hmat(ii, jj) = self%dbmat(jj, ii)
         end do
      end do
   end if

   call lapack_syevr('V', 'A', uplo, self%n, hmat, self%n, vl, vu, 1, self%n, abstol, &
      & m, eval, self%dbmat, self%n, self%isuppz, self%dwork, ldwork, self%iwork, liwork, &
      & info)
   call handle_info(error, info)
   if (allocated(error)) return

   ! Backtransform eigenvectors to the original problem.
   do ii = 1, self%n
      hmat(ii, ii) = self%dchole(ii)
   end do

   uplo_new = merge('L', 'U', upper)
   upper = .not.upper

   ! For A*x=(lambda)*B*x
   ! backtransform eigenvectors: x = inv(L)'*y or inv(U)*y'
   trans = merge('N', 'T', upper)
   call wrap_trsm(hmat, self%dbmat, side='L', uplo=uplo_new, transa=trans, diag='N')
   do ii = 1, m
      hmat(1:self%n, ii) = self%dbmat(1:self%n, ii)
   end do
   hmat(1:self%n, m+1:self%n) = 0.0_dp

end subroutine solve_dp


pure function query(n, uplo, prefix) result(lwork)
   interface
      pure integer function ilaenv(ispec, name, opts, n1, n2, n3, n4)
         integer, intent(in) :: ispec
         character(len=1), intent(in) :: name
         character(len=1), intent(in) :: opts
         integer, intent(in) :: n1
         integer, intent(in) :: n2
         integer, intent(in) :: n3
         integer, intent(in) :: n4
      end function ilaenv
   end interface
   integer, intent(in) :: n
   character(len=1), intent(in) :: uplo
   character(len=1), intent(in) :: prefix
   integer :: lwork
   integer :: nb
   nb = ilaenv(1, prefix//'SYTRD', uplo, n, -1, -1, -1)
   nb = max(nb, ilaenv(1, prefix//'ORMTR', uplo, n, -1, -1, -1))
   lwork = max(1, 26 * n, (nb + 1)*n)
end function query

pure function slamch_s() result(rmach)
   real(sp), parameter :: one = 1.0_sp, zero = 0.0_sp, rnd = one
   real(sp), parameter :: eps = merge(epsilon(zero), epsilon(zero) * 0.5_sp, one /= rnd)
   real(sp), parameter :: tinyz = tiny(zero), small = one / huge(zero)
   real(sp), parameter :: sfmin = merge(small*(one + eps), tinyz, small.ge.tinyz)
   real(sp) :: rmach
   rmach = sfmin
end function slamch_s

pure function dlamch_s() result(rmach)
   real(dp), parameter :: one = 1.0_dp, zero = 0.0_dp, rnd = one
   real(dp), parameter :: eps = merge(epsilon(zero), epsilon(zero) * 0.5_dp, one /= rnd)
   real(dp), parameter :: tinyz = tiny(zero), small = one / huge(zero)
   real(dp), parameter :: sfmin = merge(small*(one + eps), tinyz, small.ge.tinyz)
   real(dp) :: rmach
   rmach = sfmin
end function dlamch_s

subroutine handle_info(error, info)
   type(error_type), allocatable, intent(out) :: error
   integer, intent(in) :: info

   if (info /= 0) then
      call fatal_error(error, "(sygvr) failed to solve eigenvalue problem.&
         & info="//format_string(info, '(i0)'))
   end if
end subroutine handle_info

end module tblite_lapack_sygvr
