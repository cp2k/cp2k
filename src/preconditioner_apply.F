!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief computes preconditioners, and implements methods to apply them
!>      currently used in qs_ot
!> \par History
!>      - [UB] 2009-05-13 Adding stable approximate inverse (full and sparse)
!> \author Joost VandeVondele (09.2002)
! **************************************************************************************************
MODULE preconditioner_apply
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_gemm
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_release,&
                                              cp_cfm_to_fm,&
                                              cp_cfm_type,&
                                              cp_fm_to_cfm
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_get_info, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, &
        dbcsr_multiply, dbcsr_release, dbcsr_set, dbcsr_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale,&
                                              cp_fm_scale_and_add
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_restore
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE input_constants,                 ONLY: &
        ot_precond_full_all, ot_precond_full_kinetic, ot_precond_full_single, &
        ot_precond_full_single_inverse, ot_precond_s_inverse, ot_precond_solver_chebyshev, &
        ot_precond_solver_direct, ot_precond_solver_inv_chol, ot_precond_solver_update
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: z_one,&
                                              z_zero
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE preconditioner_types,            ONLY: preconditioner_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'preconditioner_apply'

   PUBLIC :: apply_preconditioner_dbcsr_complex, apply_preconditioner_fm, apply_preconditioner_dbcsr

CONTAINS

! **************************************************************************************************
!> \brief applies a previously created preconditioner to a full matrix
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
! **************************************************************************************************
   SUBROUTINE apply_preconditioner_fm(preconditioner_env, matrix_in, matrix_out)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(cp_fm_type), INTENT(IN)                       :: matrix_in, matrix_out

      CHARACTER(len=*), PARAMETER :: routineN = 'apply_preconditioner_fm'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      SELECT CASE (preconditioner_env%in_use)
      CASE (0)
         CPABORT("No preconditioner in use")
      CASE (ot_precond_full_single)
         CALL apply_full_single(preconditioner_env, matrix_in, matrix_out)
      CASE (ot_precond_full_all)
         CALL apply_full_all(preconditioner_env, matrix_in, matrix_out)
      CASE (ot_precond_full_kinetic, ot_precond_full_single_inverse, ot_precond_s_inverse)
         SELECT CASE (preconditioner_env%solver)
         CASE (ot_precond_solver_inv_chol, ot_precond_solver_update)
            CALL apply_full_single(preconditioner_env, matrix_in, matrix_out)
         CASE (ot_precond_solver_direct)
            CALL apply_full_direct(preconditioner_env, matrix_in, matrix_out)
         CASE (ot_precond_solver_chebyshev)
            CALL apply_chebyshev_fm(preconditioner_env, matrix_in, matrix_out)
         CASE DEFAULT
            CPABORT("Solver not implemented")
         END SELECT
      CASE DEFAULT
         CPABORT("Unknown preconditioner")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE apply_preconditioner_fm

! **************************************************************************************************
!> \brief ...
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
! **************************************************************************************************
   SUBROUTINE apply_preconditioner_dbcsr(preconditioner_env, matrix_in, matrix_out)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(dbcsr_type)                                   :: matrix_in, matrix_out

      CHARACTER(len=*), PARAMETER :: routineN = 'apply_preconditioner_dbcsr'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      SELECT CASE (preconditioner_env%in_use)
      CASE (0)
         CPABORT("No preconditioner in use")
      CASE (ot_precond_full_single)
         CALL apply_single(preconditioner_env, matrix_in, matrix_out)
      CASE (ot_precond_full_all)
         CALL apply_all(preconditioner_env, matrix_in, matrix_out)
      CASE (ot_precond_full_kinetic, ot_precond_full_single_inverse, ot_precond_s_inverse)
         SELECT CASE (preconditioner_env%solver)
         CASE (ot_precond_solver_inv_chol, ot_precond_solver_update)
            CALL apply_single(preconditioner_env, matrix_in, matrix_out)
         CASE (ot_precond_solver_direct)
            CPABORT("Apply_full_direct not supported with ot")
            !CALL apply_full_direct(preconditioner_env, matrix_in, matrix_out)
         CASE (ot_precond_solver_chebyshev)
            CALL apply_chebyshev_dbcsr(preconditioner_env, matrix_in, matrix_out)
         CASE DEFAULT
            CPABORT("Wrong solver")
         END SELECT
      CASE DEFAULT
         CPABORT("Wrong preconditioner")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE apply_preconditioner_dbcsr

! **************************************************************************************************
!> \brief Apply a complex orbital preconditioner to paired real/imaginary DBCSR matrices.
!> \param preconditioner_env complex k-point preconditioner
!> \param matrix_in_re real input channel
!> \param matrix_in_im imaginary input channel
!> \param matrix_out_re real output channel
!> \param matrix_out_im imaginary output channel
! **************************************************************************************************
   SUBROUTINE apply_preconditioner_dbcsr_complex(preconditioner_env, matrix_in_re, matrix_in_im, &
                                                 matrix_out_re, matrix_out_im)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(dbcsr_type)                                   :: matrix_in_re, matrix_in_im, &
                                                            matrix_out_re, matrix_out_im

      CHARACTER(len=*), PARAMETER :: routineN = 'apply_preconditioner_dbcsr_complex'

      INTEGER                                            :: handle, k, n, npre
      TYPE(cp_cfm_type)                                  :: matrix_in, matrix_out
      TYPE(cp_fm_type)                                   :: matrix_in_im_fm, matrix_in_re_fm, &
                                                            matrix_out_im_fm, matrix_out_re_fm

      CALL timeset(routineN, handle)

      SELECT CASE (preconditioner_env%in_use)
      CASE (ot_precond_full_single_inverse, ot_precond_full_kinetic, ot_precond_s_inverse)
      CASE DEFAULT
         CPABORT("Unsupported complex K-point OT preconditioner")
      END SELECT
      CPASSERT(ASSOCIATED(preconditioner_env%complex_fm))
      CALL dbcsr_get_info(matrix_in_re, nfullrows_total=n, nfullcols_total=k)
      CALL cp_cfm_get_info(preconditioner_env%complex_fm, nrow_global=npre)
      CPASSERT(n == npre)

      CALL cp_fm_create(matrix_in_re_fm, preconditioner_env%complex_fm%matrix_struct, &
                        nrow=n, ncol=k, name='complex preconditioner input real')
      CALL cp_fm_create(matrix_in_im_fm, preconditioner_env%complex_fm%matrix_struct, &
                        nrow=n, ncol=k, name='complex preconditioner input imaginary')
      CALL cp_fm_create(matrix_out_re_fm, preconditioner_env%complex_fm%matrix_struct, &
                        nrow=n, ncol=k, name='complex preconditioner output real')
      CALL cp_fm_create(matrix_out_im_fm, preconditioner_env%complex_fm%matrix_struct, &
                        nrow=n, ncol=k, name='complex preconditioner output imaginary')
      CALL copy_dbcsr_to_fm(matrix_in_re, matrix_in_re_fm)
      CALL copy_dbcsr_to_fm(matrix_in_im, matrix_in_im_fm)

      CALL cp_cfm_create(matrix_in, matrix_in_re_fm%matrix_struct, &
                         name='complex preconditioner input')
      CALL cp_cfm_create(matrix_out, matrix_in_re_fm%matrix_struct, &
                         name='complex preconditioner output')
      CALL cp_fm_to_cfm(matrix_in_re_fm, matrix_in_im_fm, matrix_in)
      CALL cp_cfm_gemm('N', 'N', n, k, n, z_one, preconditioner_env%complex_fm, &
                       matrix_in, z_zero, matrix_out)
      CALL cp_cfm_to_fm(matrix_out, matrix_out_re_fm, matrix_out_im_fm)
      CALL dbcsr_set(matrix_out_re, 0.0_dp)
      CALL dbcsr_set(matrix_out_im, 0.0_dp)
      CALL copy_fm_to_dbcsr(matrix_out_re_fm, matrix_out_re)
      CALL copy_fm_to_dbcsr(matrix_out_im_fm, matrix_out_im)

      CALL cp_cfm_release(matrix_out)
      CALL cp_cfm_release(matrix_in)
      CALL cp_fm_release(matrix_out_im_fm)
      CALL cp_fm_release(matrix_out_re_fm)
      CALL cp_fm_release(matrix_in_im_fm)
      CALL cp_fm_release(matrix_in_re_fm)

      CALL timestop(handle)

   END SUBROUTINE apply_preconditioner_dbcsr_complex

! **************************************************************************************************
!> \brief apply to full matrix, complete inversion has already been done
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
! **************************************************************************************************
   SUBROUTINE apply_full_single(preconditioner_env, matrix_in, matrix_out)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(cp_fm_type), INTENT(IN)                       :: matrix_in, matrix_out

      CHARACTER(len=*), PARAMETER                        :: routineN = 'apply_full_single'

      INTEGER                                            :: handle, k, n

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(matrix_in, nrow_global=n, ncol_global=k)
      CALL parallel_gemm('N', 'N', n, k, n, 1.0_dp, preconditioner_env%fm, &
                         matrix_in, 0.0_dp, matrix_out)
      CALL timestop(handle)

   END SUBROUTINE apply_full_single

! **************************************************************************************************
!> \brief apply to dbcsr matrix, complete inversion has already been done
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
! **************************************************************************************************
   SUBROUTINE apply_single(preconditioner_env, matrix_in, matrix_out)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(dbcsr_type)                                   :: matrix_in, matrix_out

      CHARACTER(len=*), PARAMETER                        :: routineN = 'apply_single'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (.NOT. ASSOCIATED(preconditioner_env%dbcsr_matrix)) THEN
         CPABORT("NOT ASSOCIATED preconditioner_env%dbcsr_matrix")
      END IF
      CALL dbcsr_multiply('N', 'N', 1.0_dp, preconditioner_env%dbcsr_matrix, matrix_in, &
                          0.0_dp, matrix_out)

      CALL timestop(handle)

   END SUBROUTINE apply_single

! **************************************************************************************************
!> \brief preconditioner contains the factorization, application done by
!>        solving the linear system
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
! **************************************************************************************************
   SUBROUTINE apply_full_direct(preconditioner_env, matrix_in, matrix_out)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(cp_fm_type), INTENT(IN)                       :: matrix_in, matrix_out

      CHARACTER(len=*), PARAMETER                        :: routineN = 'apply_full_direct'

      INTEGER                                            :: handle, k, n
      TYPE(cp_fm_type)                                   :: work

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(matrix_in, nrow_global=n, ncol_global=k)
      CALL cp_fm_create(work, matrix_in%matrix_struct, name="apply_full_single")
      CALL cp_fm_cholesky_restore(matrix_in, k, preconditioner_env%fm, work,&
           &                      "SOLVE", transa="T")
      CALL cp_fm_cholesky_restore(work, k, preconditioner_env%fm, matrix_out,&
           &                      "SOLVE", transa="N")
      CALL cp_fm_release(work)

      CALL timestop(handle)

   END SUBROUTINE apply_full_direct

! **************************************************************************************************
!> \brief Apply a Chebyshev approximation to the inverse of the stored dense SPD operator.
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
! **************************************************************************************************
   SUBROUTINE apply_chebyshev_fm(preconditioner_env, matrix_in, matrix_out)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(cp_fm_type), INTENT(IN)                       :: matrix_in, matrix_out

      INTEGER                                            :: degree, iteration, k, n
      REAL(KIND=dp)                                      :: delta, rho, rho_previous, sigma, theta
      TYPE(cp_fm_type)                                   :: direction, residual

      CPASSERT(ASSOCIATED(preconditioner_env%fm))
      degree = preconditioner_env%polynomial_degree
      CPASSERT(degree >= 1)
      theta = 0.5_dp*(preconditioner_env%polynomial_max + preconditioner_env%polynomial_min)
      delta = 0.5_dp*(preconditioner_env%polynomial_max - preconditioner_env%polynomial_min)
      sigma = theta/delta
      rho_previous = 1.0_dp/sigma
      CALL cp_fm_get_info(matrix_in, nrow_global=n, ncol_global=k)
      CALL cp_fm_create(direction, matrix_in%matrix_struct, name='Chebyshev direction')
      CALL cp_fm_create(residual, matrix_in%matrix_struct, name='Chebyshev residual')
      CALL cp_fm_to_fm(matrix_in, direction)
      CALL cp_fm_scale(1.0_dp/theta, direction)
      CALL cp_fm_to_fm(direction, matrix_out)
      DO iteration = 2, degree
         CALL cp_fm_to_fm(matrix_in, residual)
         CALL parallel_gemm('N', 'N', n, k, n, -1.0_dp, preconditioner_env%fm, &
                            matrix_out, 1.0_dp, residual)
         rho = 1.0_dp/(2.0_dp*sigma - rho_previous)
         CALL cp_fm_scale(rho*rho_previous, direction)
         CALL cp_fm_scale_and_add(1.0_dp, direction, 2.0_dp*rho/delta, residual)
         CALL cp_fm_scale_and_add(1.0_dp, matrix_out, 1.0_dp, direction)
         rho_previous = rho
      END DO
      CALL cp_fm_release(residual)
      CALL cp_fm_release(direction)

   END SUBROUTINE apply_chebyshev_fm

! **************************************************************************************************
!> \brief Apply a Chebyshev approximation to the inverse of the stored sparse SPD operator.
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
! **************************************************************************************************
   SUBROUTINE apply_chebyshev_dbcsr(preconditioner_env, matrix_in, matrix_out)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(dbcsr_type)                                   :: matrix_in, matrix_out

      INTEGER                                            :: degree, iteration
      REAL(KIND=dp)                                      :: delta, rho, rho_previous, sigma, theta
      TYPE(dbcsr_type)                                   :: direction, residual

      CPASSERT(ASSOCIATED(preconditioner_env%dbcsr_matrix))
      degree = preconditioner_env%polynomial_degree
      CPASSERT(degree >= 1)
      theta = 0.5_dp*(preconditioner_env%polynomial_max + preconditioner_env%polynomial_min)
      delta = 0.5_dp*(preconditioner_env%polynomial_max - preconditioner_env%polynomial_min)
      sigma = theta/delta
      rho_previous = 1.0_dp/sigma
      CALL dbcsr_copy(direction, matrix_in, name='Chebyshev direction')
      CALL dbcsr_copy(residual, matrix_in, name='Chebyshev residual')
      CALL dbcsr_set(matrix_out, 0.0_dp)
      CALL dbcsr_add(matrix_out, direction, 1.0_dp, 1.0_dp/theta)
      CALL dbcsr_set(direction, 0.0_dp)
      CALL dbcsr_add(direction, matrix_in, 1.0_dp, 1.0_dp/theta)
      DO iteration = 2, degree
         CALL dbcsr_set(residual, 0.0_dp)
         CALL dbcsr_add(residual, matrix_in, 1.0_dp, 1.0_dp)
         CALL dbcsr_multiply('N', 'N', -1.0_dp, preconditioner_env%dbcsr_matrix, &
                             matrix_out, 1.0_dp, residual)
         rho = 1.0_dp/(2.0_dp*sigma - rho_previous)
         CALL dbcsr_add(direction, residual, rho*rho_previous, 2.0_dp*rho/delta)
         CALL dbcsr_add(matrix_out, direction, 1.0_dp, 1.0_dp)
         rho_previous = rho
      END DO
      CALL dbcsr_release(residual)
      CALL dbcsr_release(direction)

   END SUBROUTINE apply_chebyshev_dbcsr

! **************************************************************************************************
!> \brief full all to a full matrix
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
! **************************************************************************************************
   SUBROUTINE apply_full_all(preconditioner_env, matrix_in, matrix_out)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(cp_fm_type), INTENT(IN)                       :: matrix_in, matrix_out

      CHARACTER(len=*), PARAMETER                        :: routineN = 'apply_full_all'

      INTEGER                                            :: handle, i, j, k, n, ncol_local, &
                                                            nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: dum
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: local_data
      TYPE(cp_fm_type)                                   :: matrix_tmp

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(matrix_in, nrow_global=n, ncol_global=k)

      CALL cp_fm_create(matrix_tmp, matrix_in%matrix_struct, name="apply_full_all")
      CALL cp_fm_get_info(matrix_tmp, nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices, local_data=local_data)

      !
      CALL parallel_gemm('T', 'N', n, k, n, 1.0_dp, preconditioner_env%fm, &
                         matrix_in, 0.0_dp, matrix_tmp)

      ! do the right scaling
      DO j = 1, ncol_local
      DO i = 1, nrow_local
         dum = 1.0_dp/MAX(preconditioner_env%energy_gap, &
                          preconditioner_env%full_evals(row_indices(i)) - preconditioner_env%occ_evals(col_indices(j)))
         local_data(i, j) = local_data(i, j)*dum
      END DO
      END DO

      ! mult back
      CALL parallel_gemm('N', 'N', n, k, n, 1.0_dp, preconditioner_env%fm, &
                         matrix_tmp, 0.0_dp, matrix_out)

      CALL cp_fm_release(matrix_tmp)

      CALL timestop(handle)

   END SUBROUTINE apply_full_all

! **************************************************************************************************
!> \brief full all to a dbcsr matrix
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
! **************************************************************************************************
   SUBROUTINE apply_all(preconditioner_env, matrix_in, matrix_out)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(dbcsr_type)                                   :: matrix_in, matrix_out

      CHARACTER(len=*), PARAMETER                        :: routineN = 'apply_all'

      INTEGER                                            :: col, col_offset, col_size, handle, i, j, &
                                                            row, row_offset, row_size
      REAL(KIND=dp)                                      :: dum
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: DATA
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_type)                                   :: matrix_tmp

      CALL timeset(routineN, handle)

      CALL dbcsr_copy(matrix_tmp, matrix_in, name="apply_full_all")
      CALL dbcsr_multiply('T', 'N', 1.0_dp, preconditioner_env%dbcsr_matrix, &
                          matrix_in, 0.0_dp, matrix_tmp)
      ! do the right scaling
      CALL dbcsr_iterator_start(iter, matrix_tmp)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, DATA, &
                                        row_size=row_size, col_size=col_size, &
                                        row_offset=row_offset, col_offset=col_offset)
         DO j = 1, col_size
         DO i = 1, row_size
            dum = 1.0_dp/MAX(preconditioner_env%energy_gap, &
                             preconditioner_env%full_evals(row_offset + i - 1) &
                             - preconditioner_env%occ_evals(col_offset + j - 1))
            DATA(i, j) = DATA(i, j)*dum
         END DO
         END DO
      END DO
      CALL dbcsr_iterator_stop(iter)

      ! mult back
      CALL dbcsr_multiply('N', 'N', 1.0_dp, preconditioner_env%dbcsr_matrix, &
                          matrix_tmp, 0.0_dp, matrix_out)
      CALL dbcsr_release(matrix_tmp)
      CALL timestop(handle)

   END SUBROUTINE apply_all

END MODULE preconditioner_apply
