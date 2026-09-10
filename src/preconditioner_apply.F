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
        dbcsr_add, dbcsr_copy, dbcsr_get_info, dbcsr_init_p, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, &
        dbcsr_multiply, dbcsr_release, dbcsr_set, dbcsr_type, dbcsr_type_no_symmetry
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_m_by_n_from_template
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale,&
                                              cp_fm_scale_and_add
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_restore
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE input_constants,                 ONLY: &
        ot_low_rank_base_lattice_fft, ot_low_rank_base_overlap, ot_precond_fermi_low_rank, &
        ot_precond_full_all, ot_precond_full_all_covariant, ot_precond_full_kinetic, &
        ot_precond_full_single, ot_precond_full_single_inverse, ot_precond_s_inverse, &
        ot_precond_solver_chebyshev, ot_precond_solver_direct, ot_precond_solver_inv_chol, &
        ot_precond_solver_update
   USE kinds,                           ONLY: dp
   USE lattice_low_rank_preconditioner, ONLY: apply_lattice_inverse_dense
   USE low_rank_preconditioner_model,   ONLY: low_rank_inverse_weight
   USE mathconstants,                   ONLY: z_one,&
                                              z_zero
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE preconditioner_types,            ONLY: preconditioner_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'preconditioner_apply'

   PUBLIC :: apply_preconditioner_cfm_complex, apply_preconditioner_dbcsr_complex, &
             apply_preconditioner_fm, apply_preconditioner_dbcsr

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
      CASE (ot_precond_fermi_low_rank)
         CALL apply_low_rank_fm(preconditioner_env, matrix_in, matrix_out)
      CASE (ot_precond_full_all_covariant)
         CALL apply_full_all_covariant(preconditioner_env, matrix_in, matrix_out)
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
      CASE (ot_precond_fermi_low_rank)
         CALL apply_low_rank_dbcsr(preconditioner_env, matrix_in, matrix_out)
      CASE (ot_precond_full_all_covariant)
         CALL apply_all_covariant(preconditioner_env, matrix_in, matrix_out)
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
!> \brief Apply a complex k-point orbital preconditioner.
!> \param preconditioner_env complex preconditioner storage
!> \param matrix_in complex input channel
!> \param matrix_out complex output channel
! **************************************************************************************************
   SUBROUTINE apply_preconditioner_cfm_complex(preconditioner_env, matrix_in, matrix_out)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(cp_cfm_type), INTENT(IN)                      :: matrix_in, matrix_out

      CHARACTER(len=*), PARAMETER :: routineN = 'apply_preconditioner_cfm_complex'

      COMPLEX(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: local_data
      INTEGER                                            :: handle, i, j, k, n, ncol_local, npre, &
                                                            nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: scale
      TYPE(cp_cfm_type)                                  :: matrix_spectral

      CALL timeset(routineN, handle)

      SELECT CASE (preconditioner_env%in_use)
      CASE (ot_precond_full_all, ot_precond_full_single, ot_precond_full_single_inverse, &
            ot_precond_full_kinetic, ot_precond_s_inverse)
      CASE DEFAULT
         CPABORT("Unsupported complex K-point OT preconditioner")
      END SELECT
      CPASSERT(ASSOCIATED(preconditioner_env%complex_fm))
      CALL cp_cfm_get_info(matrix_in, nrow_global=n, ncol_global=k)
      CALL cp_cfm_get_info(preconditioner_env%complex_fm, nrow_global=npre)
      CPASSERT(n == npre)

      IF (preconditioner_env%in_use == ot_precond_full_all) THEN
         CPASSERT(ASSOCIATED(preconditioner_env%full_evals))
         CPASSERT(ASSOCIATED(preconditioner_env%occ_evals))
         CPASSERT(n == SIZE(preconditioner_env%full_evals))
         CPASSERT(k == SIZE(preconditioner_env%occ_evals))
         CALL cp_cfm_create(matrix_spectral, matrix_in%matrix_struct, &
                            name='complex FULL_ALL spectral input')
         CALL cp_cfm_gemm('C', 'N', n, k, n, z_one, preconditioner_env%complex_fm, &
                          matrix_in, z_zero, matrix_spectral)
         CALL cp_cfm_get_info(matrix_spectral, nrow_local=nrow_local, ncol_local=ncol_local, &
                              row_indices=row_indices, col_indices=col_indices, local_data=local_data)
         DO j = 1, ncol_local
            DO i = 1, nrow_local
               scale = 1.0_dp/MAX(preconditioner_env%energy_gap, &
                                  preconditioner_env%full_evals(row_indices(i)) - &
                                  preconditioner_env%occ_evals(col_indices(j)))
               local_data(i, j) = scale*local_data(i, j)
            END DO
         END DO
         CALL cp_cfm_gemm('N', 'N', n, k, n, z_one, preconditioner_env%complex_fm, &
                          matrix_spectral, z_zero, matrix_out)
         CALL cp_cfm_release(matrix_spectral)
      ELSE
         CALL cp_cfm_gemm('N', 'N', n, k, n, z_one, preconditioner_env%complex_fm, &
                          matrix_in, z_zero, matrix_out)
      END IF

      CALL timestop(handle)

   END SUBROUTINE apply_preconditioner_cfm_complex

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

      INTEGER                                            :: handle, k, n
      TYPE(cp_cfm_type)                                  :: matrix_in, matrix_out
      TYPE(cp_fm_type)                                   :: matrix_in_im_fm, matrix_in_re_fm, &
                                                            matrix_out_im_fm, matrix_out_re_fm

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(matrix_in_re, nfullrows_total=n, nfullcols_total=k)
      CPASSERT(ASSOCIATED(preconditioner_env%complex_fm))

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
      CALL apply_preconditioner_cfm_complex(preconditioner_env, matrix_in, matrix_out)
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
!> \brief Applies the overlap-inverse base and a bounded low-rank spectral correction.
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
! **************************************************************************************************
   SUBROUTINE apply_low_rank_fm(preconditioner_env, matrix_in, matrix_out)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(cp_fm_type), INTENT(IN)                       :: matrix_in, matrix_out

      CHARACTER(len=*), PARAMETER                        :: routineN = 'apply_low_rank_fm'

      INTEGER                                            :: handle, i, j, k, n, ncol_local, nocc, &
                                                            nrow_local, rank_used
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: correction
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: local_data
      TYPE(cp_fm_type)                                   :: matrix_tmp

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(preconditioner_env%fm))
      CPASSERT(ASSOCIATED(preconditioner_env%full_evals))
      CPASSERT(ASSOCIATED(preconditioner_env%occ_evals))

      CALL cp_fm_get_info(matrix_in, nrow_global=n, ncol_global=k)
      nocc = SIZE(preconditioner_env%occ_evals)
      CPASSERT(n == SIZE(preconditioner_env%full_evals))
      rank_used = MIN(preconditioner_env%spectral_rank, n - nocc)

      ! With the exact base V V^T = S^-1. The optional lattice base replaces
      ! S^-1 by the inverse of its block-circulant translational projection:
      ! P = B^-1/window + V_r diag(1/max(gap, epsilon_a-mu) - 1/window) V_r^T.
      ! The common reference mu makes the operator right-covariant for any number
      ! of orbital columns. The lattice path retains only non-negative corrections.
      SELECT CASE (preconditioner_env%low_rank_base)
      CASE (ot_low_rank_base_overlap)
         CPASSERT(ASSOCIATED(preconditioner_env%base_fm))
         CALL parallel_gemm('N', 'N', n, k, n, preconditioner_env%spectral_base_scale, &
                            preconditioner_env%base_fm, matrix_in, 0.0_dp, matrix_out)
      CASE (ot_low_rank_base_lattice_fft)
         CALL apply_lattice_base_fm(preconditioner_env, matrix_in, matrix_out)
      CASE DEFAULT
         CPABORT("Unknown FERMI_LOW_RANK base operator during application")
      END SELECT

      IF (rank_used > 0) THEN
         CALL cp_fm_create(matrix_tmp, matrix_in%matrix_struct, name=routineN, &
                           nrow=rank_used, ncol=k)
         CALL cp_fm_get_info(matrix_tmp, nrow_local=nrow_local, ncol_local=ncol_local, &
                             row_indices=row_indices, col_indices=col_indices, local_data=local_data)

         CALL parallel_gemm('T', 'N', rank_used, k, n, 1.0_dp, preconditioner_env%fm, &
                            matrix_in, 0.0_dp, matrix_tmp, a_first_col=nocc + 1)

         DO j = 1, ncol_local
            DO i = 1, nrow_local
               correction = low_rank_inverse_weight( &
                            preconditioner_env%full_evals(nocc + row_indices(i)), &
                            preconditioner_env%spectral_reference, &
                            preconditioner_env%energy_gap) &
                            - preconditioner_env%spectral_base_scale
               local_data(i, j) = local_data(i, j)*correction
            END DO
         END DO

         CALL parallel_gemm('N', 'N', n, k, rank_used, 1.0_dp, preconditioner_env%fm, &
                            matrix_tmp, 1.0_dp, matrix_out, a_first_col=nocc + 1)
         CALL cp_fm_release(matrix_tmp)
      END IF

      CALL timestop(handle)

   END SUBROUTINE apply_low_rank_fm

! **************************************************************************************************
!> \brief Apply the replicated-cell FFT base to a distributed full matrix.
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
! **************************************************************************************************
   SUBROUTINE apply_lattice_base_fm(preconditioner_env, matrix_in, matrix_out)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(cp_fm_type), INTENT(IN)                       :: matrix_in, matrix_out

      CHARACTER(LEN=*), PARAMETER :: routineN = 'apply_lattice_base_fm'

      INTEGER                                            :: handle, i, j, k, n, ncol_local, &
                                                            nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: used_fft
      REAL(KIND=dp)                                      :: imaginary_residual, scale
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: global_input, global_output
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: input_data, output_data

      CALL timeset(routineN, handle)
      CPASSERT(ALLOCATED(preconditioner_env%lattice_inverse_k))
      CALL cp_fm_get_info(matrix_in, nrow_global=n, ncol_global=k, nrow_local=nrow_local, &
                          ncol_local=ncol_local, row_indices=row_indices, &
                          col_indices=col_indices, local_data=input_data)
      CALL cp_fm_get_info(matrix_out, local_data=output_data)
      ALLOCATE (global_input(n, k), source=0.0_dp)
      ALLOCATE (global_output(n, k))
      DO j = 1, ncol_local
         DO i = 1, nrow_local
            global_input(row_indices(i), col_indices(j)) = input_data(i, j)
         END DO
      END DO
      CALL preconditioner_env%para_env%sum(global_input)
      CALL apply_lattice_inverse_dense(preconditioner_env%lattice_inverse_k, &
                                       preconditioner_env%lattice_dims, global_input, global_output, &
                                       used_fft, imaginary_residual)
      IF (.NOT. used_fft) CPABORT("FERMI_LOW_RANK lattice inverse failed to use FFTs")
      scale = MAX(1.0_dp, MAXVAL(ABS(global_output)))
      IF (imaginary_residual > 1.0E-10_dp*scale) THEN
         CPABORT("FERMI_LOW_RANK lattice inverse produced a non-real Gamma-point result")
      END IF
      DO j = 1, ncol_local
         DO i = 1, nrow_local
            output_data(i, j) = preconditioner_env%spectral_base_scale* &
                                global_output(row_indices(i), col_indices(j))
         END DO
      END DO
      DEALLOCATE (global_input, global_output)
      CALL timestop(handle)

   END SUBROUTINE apply_lattice_base_fm

! **************************************************************************************************
!> \brief DBCSR variant of the bounded low-rank spectral correction.
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
! **************************************************************************************************
   SUBROUTINE apply_low_rank_dbcsr(preconditioner_env, matrix_in, matrix_out)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(dbcsr_type)                                   :: matrix_in, matrix_out

      CHARACTER(len=*), PARAMETER :: routineN = 'apply_low_rank_dbcsr'

      INTEGER                                            :: col, col_offset, col_size, handle, i, j, &
                                                            k, n, nocc, rank_used, row, &
                                                            row_offset, row_size
      REAL(KIND=dp)                                      :: correction
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: DATA
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_type)                                   :: matrix_tmp

      CALL timeset(routineN, handle)

      IF (preconditioner_env%low_rank_base == ot_low_rank_base_lattice_fft) THEN
         CALL apply_low_rank_dbcsr_via_fm(preconditioner_env, matrix_in, matrix_out)
         CALL timestop(handle)
         RETURN
      END IF

      CPASSERT(ASSOCIATED(preconditioner_env%base_dbcsr_matrix))
      CPASSERT(ASSOCIATED(preconditioner_env%dbcsr_matrix))
      CPASSERT(ASSOCIATED(preconditioner_env%full_evals))
      CPASSERT(ASSOCIATED(preconditioner_env%occ_evals))

      CALL dbcsr_get_info(matrix_in, nfullrows_total=n, nfullcols_total=k)
      nocc = SIZE(preconditioner_env%occ_evals)
      CPASSERT(n == SIZE(preconditioner_env%full_evals))
      rank_used = MIN(preconditioner_env%spectral_rank, n - nocc)

      ! DBCSR form of the same overlap-inverse base plus positive low-rank update.
      CALL dbcsr_multiply('N', 'N', preconditioner_env%spectral_base_scale, &
                          preconditioner_env%base_dbcsr_matrix, matrix_in, 0.0_dp, matrix_out)

      IF (rank_used > 0) THEN
         CALL dbcsr_copy(matrix_tmp, matrix_in, name=routineN)
         CALL dbcsr_set(matrix_tmp, 0.0_dp)
         CALL dbcsr_multiply('T', 'N', 1.0_dp, preconditioner_env%dbcsr_matrix, &
                             matrix_in, 0.0_dp, matrix_tmp, &
                             first_row=nocc + 1, last_row=nocc + rank_used)

         CALL dbcsr_iterator_start(iter, matrix_tmp)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, row, col, DATA, &
                                           row_size=row_size, col_size=col_size, &
                                           row_offset=row_offset, col_offset=col_offset)
            DO j = 1, col_size
               DO i = 1, row_size
                  IF (row_offset + i - 1 >= nocc + 1 .AND. &
                      row_offset + i - 1 <= nocc + rank_used) THEN
                     correction = low_rank_inverse_weight( &
                                  preconditioner_env%full_evals(row_offset + i - 1), &
                                  preconditioner_env%spectral_reference, &
                                  preconditioner_env%energy_gap) &
                                  - preconditioner_env%spectral_base_scale
                     DATA(i, j) = DATA(i, j)*correction
                  END IF
               END DO
            END DO
         END DO
         CALL dbcsr_iterator_stop(iter)

         CALL dbcsr_multiply('N', 'N', 1.0_dp, preconditioner_env%dbcsr_matrix, &
                             matrix_tmp, 1.0_dp, matrix_out, &
                             first_k=nocc + 1, last_k=nocc + rank_used)
         CALL dbcsr_release(matrix_tmp)
      END IF

      CALL timestop(handle)

   END SUBROUTINE apply_low_rank_dbcsr

! **************************************************************************************************
!> \brief Bridge DBCSR OT gradients to the full-matrix lattice-FFT implementation.
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
! **************************************************************************************************
   SUBROUTINE apply_low_rank_dbcsr_via_fm(preconditioner_env, matrix_in, matrix_out)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(dbcsr_type)                                   :: matrix_in, matrix_out

      INTEGER                                            :: k, n
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: input_fm, output_fm

      NULLIFY (fm_struct)
      CALL dbcsr_get_info(matrix_in, nfullrows_total=n, nfullcols_total=k)
      CALL cp_fm_struct_create(fm_struct, nrow_global=n, ncol_global=k, &
                               context=preconditioner_env%ctxt, &
                               para_env=preconditioner_env%para_env)
      CALL cp_fm_create(input_fm, fm_struct, name="lattice low-rank input")
      CALL cp_fm_create(output_fm, fm_struct, name="lattice low-rank output")
      CALL cp_fm_struct_release(fm_struct)
      CALL copy_dbcsr_to_fm(matrix_in, input_fm)
      CALL apply_low_rank_fm(preconditioner_env, input_fm, output_fm)
      CALL copy_fm_to_dbcsr(output_fm, matrix_out)
      CALL cp_fm_release(input_fm)
      CALL cp_fm_release(output_fm)

   END SUBROUTINE apply_low_rank_dbcsr_via_fm

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

! **************************************************************************************************
!> \brief Apply the state-selective inverse in a private canonical occupied gauge.
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
! **************************************************************************************************
   SUBROUTINE apply_full_all_covariant(preconditioner_env, matrix_in, matrix_out)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(cp_fm_type), INTENT(IN)                       :: matrix_in, matrix_out

      CHARACTER(len=*), PARAMETER :: routineN = 'apply_full_all_covariant'

      INTEGER                                            :: handle, k, n
      TYPE(cp_fm_type)                                   :: canonical_in, canonical_out

      CALL timeset(routineN, handle)
      CPASSERT(ASSOCIATED(preconditioner_env%occ_rotation))
      CALL cp_fm_get_info(matrix_in, nrow_global=n, ncol_global=k)

      CALL cp_fm_create(canonical_in, matrix_in%matrix_struct, name=routineN)
      CALL cp_fm_create(canonical_out, matrix_out%matrix_struct, name=routineN)
      CALL parallel_gemm('N', 'N', n, k, k, 1.0_dp, matrix_in, &
                         preconditioner_env%occ_rotation, 0.0_dp, canonical_in)
      CALL apply_full_all(preconditioner_env, canonical_in, canonical_out)
      CALL parallel_gemm('N', 'T', n, k, k, 1.0_dp, canonical_out, &
                         preconditioner_env%occ_rotation, 0.0_dp, matrix_out)
      CALL cp_fm_release(canonical_in)
      CALL cp_fm_release(canonical_out)
      CALL timestop(handle)

   END SUBROUTINE apply_full_all_covariant

! **************************************************************************************************
!> \brief DBCSR form of the rotation-covariant state-selective inverse.
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
! **************************************************************************************************
   SUBROUTINE apply_all_covariant(preconditioner_env, matrix_in, matrix_out)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(dbcsr_type)                                   :: matrix_in, matrix_out

      CHARACTER(len=*), PARAMETER :: routineN = 'apply_all_covariant'

      INTEGER                                            :: handle, k
      TYPE(dbcsr_type)                                   :: canonical_in, canonical_out

      CALL timeset(routineN, handle)
      CPASSERT(ASSOCIATED(preconditioner_env%occ_rotation))
      CALL dbcsr_get_info(matrix_in, nfullcols_total=k)

      IF (.NOT. ASSOCIATED(preconditioner_env%occ_rotation_dbcsr)) THEN
         CALL dbcsr_init_p(preconditioner_env%occ_rotation_dbcsr)
         CALL cp_dbcsr_m_by_n_from_template(preconditioner_env%occ_rotation_dbcsr, &
                                            template=matrix_in, m=k, n=k, &
                                            sym=dbcsr_type_no_symmetry)
         CALL copy_fm_to_dbcsr(preconditioner_env%occ_rotation, &
                               preconditioner_env%occ_rotation_dbcsr)
      END IF

      CALL dbcsr_copy(canonical_in, matrix_in, name=routineN)
      CALL dbcsr_set(canonical_in, 0.0_dp)
      CALL dbcsr_copy(canonical_out, matrix_out, name=routineN)
      CALL dbcsr_set(canonical_out, 0.0_dp)
      CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_in, &
                          preconditioner_env%occ_rotation_dbcsr, 0.0_dp, canonical_in)
      CALL apply_all(preconditioner_env, canonical_in, canonical_out)
      CALL dbcsr_multiply('N', 'T', 1.0_dp, canonical_out, &
                          preconditioner_env%occ_rotation_dbcsr, 0.0_dp, matrix_out)
      CALL dbcsr_release(canonical_in)
      CALL dbcsr_release(canonical_out)
      CALL timestop(handle)

   END SUBROUTINE apply_all_covariant

END MODULE preconditioner_apply
