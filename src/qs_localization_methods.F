!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Localization methods such as 2x2 Jacobi rotations
!>                                   Steepest Decents
!>                                   Conjugate Gradient
!> \par History
!>      Initial parallellization of jacobi (JVDV 07.2003)
!>      direct minimization using exponential parametrization (JVDV 09.2003)
!>      crazy rotations go fast (JVDV 10.2003)
!> \author CJM (04.2003)
! **************************************************************************************************
MODULE qs_localization_methods
   USE cell_types,                      ONLY: cell_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_column_scale,&
                                              cp_cfm_rot_cols,&
                                              cp_cfm_rot_rows,&
                                              cp_cfm_scale,&
                                              cp_cfm_scale_and_add,&
                                              cp_cfm_schur_product,&
                                              cp_cfm_trace
   USE cp_cfm_diag,                     ONLY: cp_cfm_heevd
   USE cp_cfm_types,                    ONLY: &
        cp_cfm_create, cp_cfm_get_element, cp_cfm_get_info, cp_cfm_get_submatrix, cp_cfm_release, &
        cp_cfm_set_all, cp_cfm_set_submatrix, cp_cfm_to_cfm, cp_cfm_to_fm, cp_cfm_type, &
        cp_fm_to_cfm
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply
   USE cp_external_control,             ONLY: external_control
   USE cp_fm_basic_linalg,              ONLY: cp_fm_frobenius_norm,&
                                              cp_fm_pdgeqpf,&
                                              cp_fm_pdorgqr,&
                                              cp_fm_scale,&
                                              cp_fm_scale_and_add,&
                                              cp_fm_trace,&
                                              cp_fm_transpose,&
                                              cp_fm_triangular_multiply
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose
   USE cp_fm_diag,                      ONLY: cp_fm_syevd
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_get,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: &
        cp_fm_create, cp_fm_get_element, cp_fm_get_info, cp_fm_get_submatrix, cp_fm_init_random, &
        cp_fm_maxabsrownorm, cp_fm_maxabsval, cp_fm_release, cp_fm_set_all, cp_fm_set_submatrix, &
        cp_fm_to_fm, cp_fm_to_fm_submat, cp_fm_type
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit,&
                                              cp_logger_get_default_unit_nr
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE kahan_sum,                       ONLY: accurate_sum
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE mathconstants,                   ONLY: pi,&
                                              twopi
   USE matrix_exp,                      ONLY: exp_pade_real,&
                                              get_nsquare_norder
   USE message_passing,                 ONLY: mp_para_env_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PUBLIC :: initialize_weights, crazy_rotations, &
             direct_mini, rotate_orbitals, approx_l1_norm_sd, jacobi_rotations, scdm_qrfact, zij_matrix, &
             jacobi_cg_edf_ls

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_localization_methods'

   PRIVATE

   TYPE set_c_1d_type
      COMPLEX(KIND=dp), POINTER, DIMENSION(:) :: c_array
   END TYPE
   TYPE set_c_2d_type
      COMPLEX(KIND=dp), POINTER, DIMENSION(:, :) :: c_array
   END TYPE

CONTAINS
! **************************************************************************************************
!> \brief ...
!> \param C ...
!> \param iterations ...
!> \param eps ...
!> \param converged ...
!> \param sweeps ...
! **************************************************************************************************
   SUBROUTINE approx_l1_norm_sd(C, iterations, eps, converged, sweeps)
      TYPE(cp_fm_type), INTENT(IN)                       :: C
      INTEGER, INTENT(IN)                                :: iterations
      REAL(KIND=dp), INTENT(IN)                          :: eps
      LOGICAL, INTENT(INOUT)                             :: converged
      INTEGER, INTENT(INOUT)                             :: sweeps

      CHARACTER(len=*), PARAMETER                        :: routineN = 'approx_l1_norm_sd'
      INTEGER, PARAMETER                                 :: taylor_order = 100
      REAL(KIND=dp), PARAMETER                           :: alpha = 0.1_dp, f2_eps = 0.01_dp

      INTEGER                                            :: handle, i, istep, k, n, ncol_local, &
                                                            nrow_local, output_unit, p
      REAL(KIND=dp)                                      :: expfactor, f2, f2old, gnorm, tnorm
      TYPE(cp_blacs_env_type), POINTER                   :: context
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_k_k
      TYPE(cp_fm_type)                                   :: CTmp, G, Gp1, Gp2, U
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      NULLIFY (context, para_env, fm_struct_k_k)

      output_unit = cp_logger_get_default_io_unit()

      CALL cp_fm_struct_get(C%matrix_struct, nrow_global=n, ncol_global=k, &
                            nrow_local=nrow_local, ncol_local=ncol_local, &
                            para_env=para_env, context=context)
      CALL cp_fm_struct_create(fm_struct_k_k, para_env=para_env, context=context, &
                               nrow_global=k, ncol_global=k)
      CALL cp_fm_create(CTmp, C%matrix_struct)
      CALL cp_fm_create(U, fm_struct_k_k)
      CALL cp_fm_create(G, fm_struct_k_k)
      CALL cp_fm_create(Gp1, fm_struct_k_k)
      CALL cp_fm_create(Gp2, fm_struct_k_k)
      !
      ! printing
      IF (output_unit > 0) THEN
         WRITE (output_unit, '(1X)')
         WRITE (output_unit, '(2X,A)') '-----------------------------------------------------------------------------'
         WRITE (output_unit, '(A,I5)') '      Nbr iterations =', iterations
         WRITE (output_unit, '(A,E10.2)') '     eps convergence =', eps
         WRITE (output_unit, '(A,I5)') '    Max Taylor order =', taylor_order
         WRITE (output_unit, '(A,E10.2)') '              f2 eps =', f2_eps
         WRITE (output_unit, '(A,E10.2)') '               alpha =', alpha
         WRITE (output_unit, '(A)') '     iteration    approx_l1_norm    g_norm   rel_err'
      END IF
      !
      f2old = 0.0_dp
      converged = .FALSE.
      !
      ! Start the steepest descent
      DO istep = 1, iterations
         !
         !-------------------------------------------------------------------
         ! compute f_2
         ! f_2(x)=(x^2+eps)^1/2
         f2 = 0.0_dp
         DO p = 1, ncol_local ! p
            DO i = 1, nrow_local ! i
               f2 = f2 + SQRT(C%local_data(i, p)**2 + f2_eps)
            END DO
         END DO
         CALL C%matrix_struct%para_env%sum(f2)
         !write(*,*) 'qs_localize: f_2=',f2
         !-------------------------------------------------------------------
         ! compute the derivative of f_2
         ! f_2(x)=(x^2+eps)^1/2
         DO p = 1, ncol_local ! p
            DO i = 1, nrow_local ! i
               CTmp%local_data(i, p) = C%local_data(i, p)/SQRT(C%local_data(i, p)**2 + f2_eps)
            END DO
         END DO
         CALL parallel_gemm('T', 'N', k, k, n, 1.0_dp, CTmp, C, 0.0_dp, G)
         ! antisymmetrize
         CALL cp_fm_transpose(G, U)
         CALL cp_fm_scale_and_add(-0.5_dp, G, 0.5_dp, U)
         !
         !-------------------------------------------------------------------
         !
         gnorm = cp_fm_frobenius_norm(G)
         !write(*,*) 'qs_localize: norm(G)=',gnorm
         !
         ! rescale for steepest descent
         CALL cp_fm_scale(-alpha, G)
         !
         ! compute unitary transform
         ! zeroth order
         CALL cp_fm_set_all(U, 0.0_dp, 1.0_dp)
         ! first order
         expfactor = 1.0_dp
         CALL cp_fm_scale_and_add(1.0_dp, U, expfactor, G)
         tnorm = cp_fm_frobenius_norm(G)
         !write(*,*) 'Taylor expansion i=',1,' norm(X^i)/i!=',tnorm
         IF (tnorm .GT. 1.0E-10_dp) THEN
            ! other orders
            CALL cp_fm_to_fm(G, Gp1)
            DO i = 2, taylor_order
               ! new power of G
               CALL parallel_gemm('N', 'N', k, k, k, 1.0_dp, G, Gp1, 0.0_dp, Gp2)
               CALL cp_fm_to_fm(Gp2, Gp1)
               ! add to the taylor expansion so far
               expfactor = expfactor/REAL(i, KIND=dp)
               CALL cp_fm_scale_and_add(1.0_dp, U, expfactor, Gp1)
               tnorm = cp_fm_frobenius_norm(Gp1)
               !write(*,*) 'Taylor expansion i=',i,' norm(X^i)/i!=',tnorm*expfactor
               IF (tnorm*expfactor .LT. 1.0E-10_dp) EXIT
            END DO
         END IF
         !
         ! incrementaly rotate the MOs
         CALL parallel_gemm('N', 'N', n, k, k, 1.0_dp, C, U, 0.0_dp, CTmp)
         CALL cp_fm_to_fm(CTmp, C)
         !
         ! printing
         IF (output_unit .GT. 0) THEN
            WRITE (output_unit, '(10X,I4,E18.10,2E10.2)') istep, f2, gnorm, ABS((f2 - f2old)/f2)
         END IF
         !
         ! Are we done?
         sweeps = istep
         !IF(gnorm.LE.grad_thresh.AND.ABS((f2-f2old)/f2).LE.f2_thresh.AND.istep.GT.1) THEN
         IF (ABS((f2 - f2old)/f2) .LE. eps .AND. istep .GT. 1) THEN
            converged = .TRUE.
            EXIT
         END IF
         f2old = f2
      END DO
      !
      ! here we should do one refine step to enforce C'*S*C=1 for any case
      !
      ! Print the final result
      IF (output_unit .GT. 0) WRITE (output_unit, '(A,E16.10)') ' sparseness function f2 = ', f2
      !
      ! sparsity
      !DO i=1,size(thresh,1)
      !   gnorm = 0.0_dp
      !   DO o=1,ncol_local
      !      DO p=1,nrow_local
      !         IF(ABS(C%local_data(p,o)).GT.thresh(i)) THEN
      !            gnorm = gnorm + 1.0_dp
      !         ENDIF
      !      ENDDO
      !   ENDDO
      !   CALL C%matrix_struct%para_env%sum(gnorm)
      !   IF(output_unit.GT.0) THEN
      !      WRITE(output_unit,*) 'qs_localize: ratio2=',gnorm / ( REAL(k,KIND=dp)*REAL(n,KIND=dp) ),thresh(i)
      !   ENDIF
      !ENDDO
      !
      ! deallocate
      CALL cp_fm_struct_release(fm_struct_k_k)
      CALL cp_fm_release(CTmp)
      CALL cp_fm_release(U)
      CALL cp_fm_release(G)
      CALL cp_fm_release(Gp1)
      CALL cp_fm_release(Gp2)

      CALL timestop(handle)

   END SUBROUTINE approx_l1_norm_sd
! **************************************************************************************************
!> \brief ...
!> \param cell ...
!> \param weights ...
! **************************************************************************************************
   SUBROUTINE initialize_weights(cell, weights)

      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), DIMENSION(:)                        :: weights

      REAL(KIND=dp), DIMENSION(3, 3)                     :: metric

      CPASSERT(ASSOCIATED(cell))

      metric = 0.0_dp
      CALL dgemm('T', 'N', 3, 3, 3, 1._dp, cell%hmat(:, :), 3, cell%hmat(:, :), 3, 0.0_dp, metric(:, :), 3)

      weights(1) = METRIC(1, 1) - METRIC(1, 2) - METRIC(1, 3)
      weights(2) = METRIC(2, 2) - METRIC(1, 2) - METRIC(2, 3)
      weights(3) = METRIC(3, 3) - METRIC(1, 3) - METRIC(2, 3)
      weights(4) = METRIC(1, 2)
      weights(5) = METRIC(1, 3)
      weights(6) = METRIC(2, 3)

   END SUBROUTINE initialize_weights

! **************************************************************************************************
!> \brief wrapper for the jacobi routines, should be removed if jacobi_rot_para
!>        can deal with serial para_envs.
!> \param weights ...
!> \param zij ...
!> \param vectors ...
!> \param para_env ...
!> \param max_iter ...
!> \param eps_localization ...
!> \param sweeps ...
!> \param out_each ...
!> \param target_time ...
!> \param start_time ...
!> \par History
!> \author Joost VandeVondele (02.2010)
! **************************************************************************************************
   SUBROUTINE jacobi_rotations(weights, zij, vectors, para_env, max_iter, &
                               eps_localization, sweeps, out_each, target_time, start_time)

      REAL(KIND=dp), INTENT(IN)                          :: weights(:)
      TYPE(cp_fm_type), INTENT(IN)                       :: zij(:, :), vectors
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: max_iter
      REAL(KIND=dp), INTENT(IN)                          :: eps_localization
      INTEGER                                            :: sweeps
      INTEGER, INTENT(IN)                                :: out_each
      REAL(dp)                                           :: target_time, start_time

      IF (para_env%num_pe == 1) THEN
         CALL jacobi_rotations_serial(weights, zij, vectors, max_iter, eps_localization, sweeps, out_each)
      ELSE
         CALL jacobi_rot_para(weights, zij, vectors, para_env, max_iter, eps_localization, &
                              sweeps, out_each, target_time, start_time)
      END IF

   END SUBROUTINE jacobi_rotations

! **************************************************************************************************
!> \brief this routine, private to the module is a serial backup, till we have jacobi_rot_para to work in serial
!>        while the routine below works in parallel, it is too slow to be useful
!> \param weights ...
!> \param zij ...
!> \param vectors ...
!> \param max_iter ...
!> \param eps_localization ...
!> \param sweeps ...
!> \param out_each ...
! **************************************************************************************************
   SUBROUTINE jacobi_rotations_serial(weights, zij, vectors, max_iter, eps_localization, sweeps, &
                                      out_each)
      REAL(KIND=dp), INTENT(IN)                          :: weights(:)
      TYPE(cp_fm_type), INTENT(IN)                       :: zij(:, :), vectors
      INTEGER, INTENT(IN)                                :: max_iter
      REAL(KIND=dp), INTENT(IN)                          :: eps_localization
      INTEGER                                            :: sweeps
      INTEGER, INTENT(IN)                                :: out_each

      CHARACTER(len=*), PARAMETER :: routineN = 'jacobi_rotations_serial'

      COMPLEX(KIND=dp), POINTER                          :: mii(:), mij(:), mjj(:)
      INTEGER                                            :: dim2, handle, idim, istate, jstate, &
                                                            nstate, unit_nr
      REAL(KIND=dp)                                      :: ct, st, t1, t2, theta, tolerance
      TYPE(cp_cfm_type)                                  :: c_rmat
      TYPE(cp_cfm_type), ALLOCATABLE, DIMENSION(:)       :: c_zij
      TYPE(cp_fm_type)                                   :: rmat

      CALL timeset(routineN, handle)

      dim2 = SIZE(zij, 2)
      ALLOCATE (c_zij(dim2))
      NULLIFY (mii, mij, mjj)
      ALLOCATE (mii(dim2), mij(dim2), mjj(dim2))

      CALL cp_fm_create(rmat, zij(1, 1)%matrix_struct)
      CALL cp_fm_set_all(rmat, 0._dp, 1._dp)

      CALL cp_cfm_create(c_rmat, zij(1, 1)%matrix_struct)
      CALL cp_cfm_set_all(c_rmat, (0._dp, 0._dp), (1._dp, 0._dp))
      DO idim = 1, dim2
         CALL cp_cfm_create(c_zij(idim), zij(1, 1)%matrix_struct)
         c_zij(idim)%local_data = CMPLX(zij(1, idim)%local_data, &
                                        zij(2, idim)%local_data, dp)
      END DO

      CALL cp_fm_get_info(rmat, nrow_global=nstate)
      tolerance = 1.0e10_dp

      sweeps = 0
      unit_nr = -1
      IF (rmat%matrix_struct%para_env%is_source()) THEN
         unit_nr = cp_logger_get_default_unit_nr()
         WRITE (unit_nr, '(T4,A )') " Localization by iterative Jacobi rotation"
      END IF

      ! do jacobi sweeps until converged
      DO WHILE (tolerance >= eps_localization .AND. sweeps < max_iter)
         sweeps = sweeps + 1
         t1 = m_walltime()

         DO istate = 1, nstate
            DO jstate = istate + 1, nstate
               DO idim = 1, dim2
                  CALL cp_cfm_get_element(c_zij(idim), istate, istate, mii(idim))
                  CALL cp_cfm_get_element(c_zij(idim), istate, jstate, mij(idim))
                  CALL cp_cfm_get_element(c_zij(idim), jstate, jstate, mjj(idim))
               END DO
               CALL get_angle(mii, mjj, mij, weights, theta)
               st = SIN(theta)
               ct = COS(theta)
               CALL rotate_zij(istate, jstate, st, ct, c_zij)

               CALL rotate_rmat(istate, jstate, st, ct, c_rmat)
            END DO
         END DO

         CALL check_tolerance(c_zij, weights, tolerance)

         t2 = m_walltime()
         IF (unit_nr > 0 .AND. MODULO(sweeps, out_each) == 0) THEN
            WRITE (unit_nr, '(T4,A,I7,T30,A,E12.4,T60,A,F8.3)') &
               "Iteration:", sweeps, "Tolerance:", tolerance, "Time:", t2 - t1
            CALL m_flush(unit_nr)
         END IF

      END DO

      DO idim = 1, dim2
         zij(1, idim)%local_data = REAL(c_zij(idim)%local_data, dp)
         zij(2, idim)%local_data = AIMAG(c_zij(idim)%local_data)
         CALL cp_cfm_release(c_zij(idim))
      END DO
      DEALLOCATE (c_zij)
      DEALLOCATE (mii, mij, mjj)
      rmat%local_data = REAL(c_rmat%local_data, dp)
      CALL cp_cfm_release(c_rmat)
      CALL rotate_orbitals(rmat, vectors)
      CALL cp_fm_release(rmat)

      CALL timestop(handle)

   END SUBROUTINE jacobi_rotations_serial
! **************************************************************************************************
!> \brief very similar to jacobi_rotations_serial with some extra output options
!> \param weights ...
!> \param c_zij ...
!> \param max_iter ...
!> \param c_rmat ...
!> \param eps_localization ...
!> \param tol_out ...
!> \param jsweeps ...
!> \param out_each ...
!> \param c_zij_out ...
!> \param grad_final ...
! **************************************************************************************************
   SUBROUTINE jacobi_rotations_serial_1(weights, c_zij, max_iter, c_rmat, eps_localization, &
                                        tol_out, jsweeps, out_each, c_zij_out, grad_final)
      REAL(KIND=dp), INTENT(IN)                          :: weights(:)
      TYPE(cp_cfm_type), INTENT(IN)                      :: c_zij(:)
      INTEGER, INTENT(IN)                                :: max_iter
      TYPE(cp_cfm_type), INTENT(IN)                      :: c_rmat
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: eps_localization
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: tol_out
      INTEGER, INTENT(OUT), OPTIONAL                     :: jsweeps
      INTEGER, INTENT(IN), OPTIONAL                      :: out_each
      TYPE(cp_cfm_type), INTENT(IN), OPTIONAL            :: c_zij_out(:)
      TYPE(cp_fm_type), INTENT(OUT), OPTIONAL, POINTER   :: grad_final

      CHARACTER(len=*), PARAMETER :: routineN = 'jacobi_rotations_serial_1'

      COMPLEX(KIND=dp)                                   :: mzii
      COMPLEX(KIND=dp), POINTER                          :: mii(:), mij(:), mjj(:)
      INTEGER                                            :: dim2, handle, idim, istate, jstate, &
                                                            nstate, sweeps, unit_nr
      REAL(KIND=dp)                                      :: alpha, avg_spread_ii, ct, spread_ii, st, &
                                                            sum_spread_ii, t1, t2, theta, tolerance
      TYPE(cp_cfm_type)                                  :: c_rmat_local
      TYPE(cp_cfm_type), ALLOCATABLE                     :: c_zij_local(:)

      CALL timeset(routineN, handle)

      dim2 = SIZE(c_zij)
      NULLIFY (mii, mij, mjj)
      ALLOCATE (mii(dim2), mij(dim2), mjj(dim2))

      ALLOCATE (c_zij_local(dim2))
      CALL cp_cfm_create(c_rmat_local, c_rmat%matrix_struct)
      CALL cp_cfm_set_all(c_rmat_local, (0.0_dp, 0.0_dp), (1.0_dp, 0.0_dp))
      DO idim = 1, dim2
         CALL cp_cfm_create(c_zij_local(idim), c_zij(idim)%matrix_struct)
         c_zij_local(idim)%local_data = c_zij(idim)%local_data
      END DO

      CALL cp_cfm_get_info(c_rmat_local, nrow_global=nstate)
      tolerance = 1.0e10_dp

      IF (PRESENT(grad_final)) CALL cp_fm_set_all(grad_final, 0.0_dp)

      sweeps = 0
      IF (PRESENT(out_each)) THEN
         unit_nr = -1
         IF (c_rmat_local%matrix_struct%para_env%is_source()) THEN
            unit_nr = cp_logger_get_default_unit_nr()
         END IF
         alpha = 0.0_dp
         DO idim = 1, dim2
            alpha = alpha + weights(idim)
         END DO
      END IF

      ! do jacobi sweeps until converged
      DO WHILE (sweeps < max_iter)
         sweeps = sweeps + 1
         IF (PRESENT(eps_localization)) THEN
            IF (tolerance < eps_localization) EXIT
         END IF
         IF (PRESENT(out_each)) t1 = m_walltime()

         DO istate = 1, nstate
            DO jstate = istate + 1, nstate
               DO idim = 1, dim2
                  CALL cp_cfm_get_element(c_zij_local(idim), istate, istate, mii(idim))
                  CALL cp_cfm_get_element(c_zij_local(idim), istate, jstate, mij(idim))
                  CALL cp_cfm_get_element(c_zij_local(idim), jstate, jstate, mjj(idim))
               END DO
               CALL get_angle(mii, mjj, mij, weights, theta)
               st = SIN(theta)
               ct = COS(theta)
               CALL rotate_zij(istate, jstate, st, ct, c_zij_local)

               CALL rotate_rmat(istate, jstate, st, ct, c_rmat_local)
            END DO
         END DO

         IF (PRESENT(grad_final)) THEN
            CALL check_tolerance(c_zij_local, weights, tolerance, grad=grad_final)
         ELSE
            CALL check_tolerance(c_zij_local, weights, tolerance)
         END IF
         IF (PRESENT(tol_out)) tol_out = tolerance

         IF (PRESENT(out_each)) THEN
            t2 = m_walltime()
            IF (unit_nr > 0 .AND. MODULO(sweeps, out_each) == 0) THEN
               sum_spread_ii = 0.0_dp
               DO istate = 1, nstate
                  spread_ii = 0.0_dp
                  DO idim = 1, dim2
                     CALL cp_cfm_get_element(c_zij_local(idim), istate, istate, mzii)
                     spread_ii = spread_ii + weights(idim)* &
                                 ABS(mzii)**2/twopi/twopi
                  END DO
                  sum_spread_ii = sum_spread_ii + spread_ii
               END DO
               sum_spread_ii = alpha*nstate/twopi/twopi - sum_spread_ii
               avg_spread_ii = sum_spread_ii/nstate
               WRITE (unit_nr, '(T4,A,T26,A,T48,A,T64,A)') &
                  "Iteration", "Avg. Spread_ii", "Tolerance", "Time"
               WRITE (unit_nr, '(T4,I7,T20,F20.10,T45,E12.4,T60,F8.3)') &
                  sweeps, avg_spread_ii, tolerance, t2 - t1
               CALL m_flush(unit_nr)
            END IF
            IF (PRESENT(jsweeps)) jsweeps = sweeps
         END IF

      END DO

      IF (PRESENT(c_zij_out)) THEN
         DO idim = 1, dim2
            CALL cp_cfm_to_cfm(c_zij_local(idim), c_zij_out(idim))
         END DO
      END IF
      CALL cp_cfm_to_cfm(c_rmat_local, c_rmat)

      DEALLOCATE (mii, mij, mjj)
      DO idim = 1, dim2
         CALL cp_cfm_release(c_zij_local(idim))
      END DO
      DEALLOCATE (c_zij_local)
      CALL cp_cfm_release(c_rmat_local)

      CALL timestop(handle)

   END SUBROUTINE jacobi_rotations_serial_1
! **************************************************************************************************
!> \brief combine jacobi rotations (serial) and conjugate gradient with golden section line search
!>        for partially occupied wannier functions
!> \param para_env ...
!> \param weights ...
!> \param zij ...
!> \param vectors ...
!> \param max_iter ...
!> \param eps_localization ...
!> \param iter ...
!> \param out_each ...
!> \param nextra ...
!> \param do_cg ...
!> \param nmo ...
!> \param vectors_2 ...
!> \param mos_guess ...
! **************************************************************************************************
   SUBROUTINE jacobi_cg_edf_ls(para_env, weights, zij, vectors, max_iter, eps_localization, &
                               iter, out_each, nextra, do_cg, nmo, vectors_2, mos_guess)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      REAL(KIND=dp), INTENT(IN)                          :: weights(:)
      TYPE(cp_fm_type), INTENT(IN)                       :: zij(:, :), vectors
      INTEGER, INTENT(IN)                                :: max_iter
      REAL(KIND=dp), INTENT(IN)                          :: eps_localization
      INTEGER                                            :: iter
      INTEGER, INTENT(IN)                                :: out_each, nextra
      LOGICAL, INTENT(IN)                                :: do_cg
      INTEGER, INTENT(IN), OPTIONAL                      :: nmo
      TYPE(cp_fm_type), INTENT(IN), OPTIONAL             :: vectors_2, mos_guess

      CHARACTER(len=*), PARAMETER                        :: routineN = 'jacobi_cg_edf_ls'
      COMPLEX(KIND=dp), PARAMETER                        :: cone = (1.0_dp, 0.0_dp), &
                                                            czero = (0.0_dp, 0.0_dp)
      REAL(KIND=dp), PARAMETER                           :: gold_sec = 0.3819_dp

      COMPLEX(KIND=dp)                                   :: cnorm2_Gct, cnorm2_Gct_cross, mzii
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: tmp_cmat
      COMPLEX(KIND=dp), DIMENSION(:), POINTER            :: arr_zii
      COMPLEX(KIND=dp), DIMENSION(:, :), POINTER         :: matrix_zii
      INTEGER :: dim2, handle, icinit, idim, istate, line_search_count, line_searches, lsl, lsm, &
         lsr, miniter, nao, ndummy, nocc, norextra, northo, nstate, unit_nr
      INTEGER, DIMENSION(1)                              :: iloc
      LOGICAL                                            :: do_cinit_mo, do_cinit_random, &
                                                            do_U_guess_mo, new_direction
      REAL(KIND=dp) :: alpha, avg_spread_ii, beta, beta_pr, ds, ds_min, mintol, norm, norm2_Gct, &
         norm2_Gct_cross, norm2_old, spread_ii, spread_sum, sum_spread_ii, t1, tol, tolc, weight
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: sum_spread
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: tmp_mat, tmp_mat_1
      REAL(KIND=dp), DIMENSION(50)                       :: energy, pos
      REAL(KIND=dp), DIMENSION(:), POINTER               :: tmp_arr
      TYPE(cp_blacs_env_type), POINTER                   :: context
      TYPE(cp_cfm_type)                                  :: c_tilde, ctrans_lambda, Gct_old, &
                                                            grad_ctilde, skc, tmp_cfm, tmp_cfm_1, &
                                                            tmp_cfm_2, U, UL, V, VL, zdiag
      TYPE(cp_cfm_type), ALLOCATABLE, DIMENSION(:)       :: c_zij, zij_0
      TYPE(cp_fm_struct_type), POINTER                   :: tmp_fm_struct
      TYPE(cp_fm_type)                                   :: id_nextra, matrix_U, matrix_V, &
                                                            matrix_V_all, rmat, tmp_fm, vectors_all

      CALL timeset(routineN, handle)

      dim2 = SIZE(zij, 2)
      NULLIFY (context)
      NULLIFY (matrix_zii, arr_zii)
      NULLIFY (tmp_fm_struct)
      NULLIFY (tmp_arr)

      ALLOCATE (c_zij(dim2))

      CALL cp_fm_get_info(zij(1, 1), nrow_global=nstate)

      ALLOCATE (sum_spread(nstate))
      ALLOCATE (matrix_zii(nstate, dim2))
      matrix_zii = czero
      sum_spread = 0.0_dp

      alpha = 0.0_dp
      DO idim = 1, dim2
         alpha = alpha + weights(idim)
         CALL cp_cfm_create(c_zij(idim), zij(1, 1)%matrix_struct)
         c_zij(idim)%local_data = CMPLX(zij(1, idim)%local_data, &
                                        zij(2, idim)%local_data, dp)
      END DO

      ALLOCATE (zij_0(dim2))

      CALL cp_cfm_create(U, zij(1, 1)%matrix_struct)
      CALL cp_fm_create(matrix_U, zij(1, 1)%matrix_struct)

      CALL cp_cfm_set_all(U, czero, cone)
      CALL cp_fm_set_all(matrix_U, 0.0_dp, 1.0_dp)

      CALL cp_fm_get_info(vectors, nrow_global=nao)
      IF (nextra > 0) THEN
         IF (PRESENT(mos_guess)) THEN
            do_cinit_random = .FALSE.
            do_cinit_mo = .TRUE.
            CALL cp_fm_get_info(mos_guess, ncol_global=ndummy)
         ELSE
            do_cinit_random = .TRUE.
            do_cinit_mo = .FALSE.
            ndummy = nstate
         END IF

         IF (do_cinit_random) THEN
            icinit = 1
            do_U_guess_mo = .FALSE.
         ELSEIF (do_cinit_mo) THEN
            icinit = 2
            do_U_guess_mo = .TRUE.
         END IF

         nocc = nstate - nextra
         northo = nmo - nocc
         norextra = nmo - nstate
         CALL cp_fm_struct_get(zij(1, 1)%matrix_struct, context=context)

         ALLOCATE (tmp_cmat(nstate, nstate))
         CALL cp_fm_struct_create(tmp_fm_struct, nrow_global=nmo, ncol_global=nmo, &
                                  para_env=para_env, context=context)
         DO idim = 1, dim2
            CALL cp_cfm_create(zij_0(idim), tmp_fm_struct)
            CALL cp_cfm_set_all(zij_0(idim), czero, cone)
            CALL cp_cfm_get_submatrix(c_zij(idim), tmp_cmat)
            CALL cp_cfm_set_submatrix(zij_0(idim), tmp_cmat)
         END DO
         CALL cp_fm_struct_release(tmp_fm_struct)
         DEALLOCATE (tmp_cmat)

         CALL cp_fm_struct_create(tmp_fm_struct, nrow_global=nmo, ncol_global=nstate, &
                                  para_env=para_env, context=context)
         CALL cp_cfm_create(V, tmp_fm_struct)
         CALL cp_fm_create(matrix_V, tmp_fm_struct)
         CALL cp_cfm_create(zdiag, tmp_fm_struct)
         CALL cp_fm_create(rmat, tmp_fm_struct)
         CALL cp_fm_struct_release(tmp_fm_struct)
         CALL cp_cfm_set_all(V, czero, cone)
         CALL cp_fm_set_all(matrix_V, 0.0_dp, 1.0_dp)

         CALL cp_fm_struct_create(tmp_fm_struct, nrow_global=nmo, ncol_global=ndummy, &
                                  para_env=para_env, context=context)
         CALL cp_fm_create(matrix_V_all, tmp_fm_struct)
         CALL cp_fm_struct_release(tmp_fm_struct)
         CALL cp_fm_set_all(matrix_V_all, 0._dp, 1._dp)

         ALLOCATE (arr_zii(nstate))

         CALL cp_fm_struct_create(tmp_fm_struct, nrow_global=northo, ncol_global=nextra, &
                                  para_env=para_env, context=context)
         CALL cp_cfm_create(c_tilde, tmp_fm_struct)
         CALL cp_cfm_create(grad_ctilde, tmp_fm_struct)
         CALL cp_cfm_create(Gct_old, tmp_fm_struct)
         CALL cp_cfm_create(skc, tmp_fm_struct)
         CALL cp_fm_struct_release(tmp_fm_struct)
         CALL cp_cfm_set_all(c_tilde, czero)
         CALL cp_cfm_set_all(Gct_old, czero)
         CALL cp_cfm_set_all(skc, czero)

         CALL cp_fm_struct_create(tmp_fm_struct, nrow_global=northo, ncol_global=nstate, &
                                  para_env=para_env, context=context)
         CALL cp_cfm_create(VL, tmp_fm_struct)
         CALL cp_cfm_set_all(VL, czero)
         CALL cp_fm_struct_release(tmp_fm_struct)

         CALL cp_fm_struct_create(tmp_fm_struct, nrow_global=nextra, ncol_global=nextra, &
                                  para_env=para_env, context=context)
         CALL cp_fm_create(id_nextra, tmp_fm_struct)
         CALL cp_cfm_create(ctrans_lambda, tmp_fm_struct)
         CALL cp_fm_struct_release(tmp_fm_struct)
         CALL cp_cfm_set_all(ctrans_lambda, czero)
         CALL cp_fm_set_all(id_nextra, 0.0_dp, 1.0_dp)

         CALL cp_fm_struct_create(tmp_fm_struct, nrow_global=nextra, ncol_global=nstate, &
                                  para_env=para_env, context=context)
         CALL cp_cfm_create(UL, tmp_fm_struct)
         CALL cp_fm_struct_release(tmp_fm_struct)
         CALL cp_cfm_set_all(UL, czero)

         CALL cp_fm_struct_create(tmp_fm_struct, nrow_global=nao, ncol_global=nmo, &
                                  para_env=para_env, context=context)
         CALL cp_fm_create(vectors_all, tmp_fm_struct)
         CALL cp_fm_struct_release(tmp_fm_struct)
         ALLOCATE (tmp_mat(nao, nstate))
         CALL cp_fm_get_submatrix(vectors, tmp_mat)
         CALL cp_fm_set_submatrix(vectors_all, tmp_mat, 1, 1, nao, nstate)
         DEALLOCATE (tmp_mat)
         ALLOCATE (tmp_mat(nao, norextra))
         CALL cp_fm_get_submatrix(vectors_2, tmp_mat)
         CALL cp_fm_set_submatrix(vectors_all, tmp_mat, 1, nstate + 1, nao, norextra)
         DEALLOCATE (tmp_mat)

         ! initialize c_tilde
         SELECT CASE (icinit)
         CASE (1) ! random coefficients
            !WRITE (*, *) "RANDOM INITIAL GUESS FOR C"
            CALL cp_fm_create(tmp_fm, c_tilde%matrix_struct)
            CALL cp_fm_init_random(tmp_fm, nextra)
            CALL ortho_vectors(tmp_fm)
            c_tilde%local_data = tmp_fm%local_data
            CALL cp_fm_release(tmp_fm)
            ALLOCATE (tmp_cmat(northo, nextra))
            CALL cp_cfm_get_submatrix(c_tilde, tmp_cmat)
            CALL cp_cfm_set_submatrix(V, tmp_cmat, nocc + 1, nocc + 1, northo, nextra)
            DEALLOCATE (tmp_cmat)
         CASE (2) ! MO based coeffs
            CALL parallel_gemm("T", "N", nmo, ndummy, nao, 1.0_dp, vectors_all, mos_guess, 0.0_dp, matrix_V_all)
            ALLOCATE (tmp_arr(nmo))
            ALLOCATE (tmp_mat(nmo, ndummy))
            ALLOCATE (tmp_mat_1(nmo, nstate))
            ! normalize matrix_V_all
            CALL cp_fm_get_submatrix(matrix_V_all, tmp_mat)
            DO istate = 1, ndummy
               tmp_arr(:) = tmp_mat(:, istate)
               norm = norm2(tmp_arr)
               tmp_arr(:) = tmp_arr(:)/norm
               tmp_mat(:, istate) = tmp_arr(:)
            END DO
            CALL cp_fm_set_submatrix(matrix_V_all, tmp_mat)
            CALL cp_fm_get_submatrix(matrix_V_all, tmp_mat_1, 1, 1, nmo, nstate)
            CALL cp_fm_set_submatrix(matrix_V, tmp_mat_1)
            DEALLOCATE (tmp_arr, tmp_mat, tmp_mat_1)
            CALL cp_fm_to_cfm(msourcer=matrix_V, mtarget=V)
            ALLOCATE (tmp_mat(northo, ndummy))
            ALLOCATE (tmp_mat_1(northo, nextra))
            CALL cp_fm_get_submatrix(matrix_V_all, tmp_mat, nocc + 1, 1, northo, ndummy)
            ALLOCATE (tmp_arr(ndummy))
            tmp_arr = 0.0_dp
            DO istate = 1, ndummy
               tmp_arr(istate) = norm2(tmp_mat(:, istate))
            END DO
            ! find edfs
            DO istate = 1, nextra
               iloc = MAXLOC(tmp_arr)
               tmp_mat_1(:, istate) = tmp_mat(:, iloc(1))
               tmp_arr(iloc(1)) = 0.0_dp
            END DO

            DEALLOCATE (tmp_arr, tmp_mat)

            CALL cp_fm_struct_create(tmp_fm_struct, nrow_global=northo, ncol_global=nextra, &
                                     para_env=para_env, context=context)
            CALL cp_fm_create(tmp_fm, tmp_fm_struct)
            CALL cp_fm_struct_release(tmp_fm_struct)
            CALL cp_fm_set_submatrix(tmp_fm, tmp_mat_1)
            DEALLOCATE (tmp_mat_1)
            CALL ortho_vectors(tmp_fm)
            CALL cp_fm_to_cfm(msourcer=tmp_fm, mtarget=c_tilde)
            CALL cp_fm_release(tmp_fm)
            ! initialize U
            IF (do_U_guess_mo) THEN
               ALLOCATE (tmp_cmat(nocc, nstate))
               CALL cp_cfm_get_submatrix(V, tmp_cmat, 1, 1, nocc, nstate)
               CALL cp_cfm_set_submatrix(U, tmp_cmat, 1, 1, nocc, nstate)
               DEALLOCATE (tmp_cmat)
               ALLOCATE (tmp_cmat(northo, nstate))
               CALL cp_cfm_get_submatrix(V, tmp_cmat, nocc + 1, 1, northo, nstate)
               CALL cp_cfm_set_submatrix(VL, tmp_cmat, 1, 1, northo, nstate)
               DEALLOCATE (tmp_cmat)
               CALL parallel_gemm("C", "N", nextra, nstate, northo, cone, c_tilde, VL, czero, UL)
               ALLOCATE (tmp_cmat(nextra, nstate))
               CALL cp_cfm_get_submatrix(UL, tmp_cmat, 1, 1, nextra, nstate)
               CALL cp_cfm_set_submatrix(U, tmp_cmat, nocc + 1, 1, nextra, nstate)
               DEALLOCATE (tmp_cmat)
               CALL cp_fm_create(tmp_fm, U%matrix_struct)
               tmp_fm%local_data = REAL(U%local_data, KIND=dp)
               CALL ortho_vectors(tmp_fm)
               CALL cp_fm_to_cfm(msourcer=tmp_fm, mtarget=U)
               CALL cp_fm_release(tmp_fm)
               CALL cp_cfm_to_fm(U, matrix_U)
            END IF
            ! reevaluate V
            ALLOCATE (tmp_cmat(nocc, nstate))
            CALL cp_cfm_get_submatrix(U, tmp_cmat, 1, 1, nocc, nstate)
            CALL cp_cfm_set_submatrix(V, tmp_cmat, 1, 1, nocc, nstate)
            DEALLOCATE (tmp_cmat)
            ALLOCATE (tmp_cmat(nextra, nstate))
            CALL cp_cfm_get_submatrix(U, tmp_cmat, nocc + 1, 1, nextra, nstate)
            CALL cp_cfm_set_submatrix(UL, tmp_cmat, 1, 1, nextra, nstate)
            DEALLOCATE (tmp_cmat)
            CALL parallel_gemm("N", "N", northo, nstate, nextra, cone, c_tilde, UL, czero, VL)
            ALLOCATE (tmp_cmat(northo, nstate))
            CALL cp_cfm_get_submatrix(VL, tmp_cmat)
            CALL cp_cfm_set_submatrix(V, tmp_cmat, nocc + 1, 1, northo, nstate)
            DEALLOCATE (tmp_cmat)
         END SELECT
      ELSE
         DO idim = 1, dim2
            CALL cp_cfm_create(zij_0(idim), zij(1, 1)%matrix_struct)
            CALL cp_cfm_to_cfm(c_zij(idim), zij_0(idim))
         END DO
         CALL cp_fm_create(rmat, zij(1, 1)%matrix_struct)
         CALL cp_fm_set_all(rmat, 0._dp, 1._dp)
      END IF

      unit_nr = -1
      IF (rmat%matrix_struct%para_env%is_source()) THEN
         unit_nr = cp_logger_get_default_unit_nr()
         WRITE (unit_nr, '(T4,A )') " Localization by combined Jacobi rotations and Non-Linear Conjugate Gradient"
      END IF

      norm2_old = 1.0E30_dp
      ds_min = 1.0_dp
      new_direction = .TRUE.
      iter = 0
      line_searches = 0
      line_search_count = 0
      tol = 1.0E+20_dp
      mintol = 1.0E+10_dp
      miniter = 0

      !IF (nextra > 0) WRITE(*,*) 'random_guess, MO_guess, U_guess, conjugate_gradient: ', &
      !                            do_cinit_random, do_cinit_mo, do_U_guess_mo, do_cg

      ! do conjugate gradient until converged
      DO WHILE (iter < max_iter)
         iter = iter + 1
         !WRITE(*,*) 'iter = ', iter
         t1 = m_walltime()

         IF (iter > 1) THEN
            ! comput U
            CALL cp_cfm_create(tmp_cfm, zij(1, 1)%matrix_struct)
            CALL cp_cfm_create(tmp_cfm_2, zij(1, 1)%matrix_struct)
            IF (para_env%num_pe == 1) THEN
               CALL jacobi_rotations_serial_1(weights, c_zij, 1, tmp_cfm_2, tol_out=tol)
            ELSE
               CALL jacobi_rot_para_1(weights, c_zij, para_env, 1, tmp_cfm_2, tol_out=tol)
            END IF
            CALL parallel_gemm('N', 'N', nstate, nstate, nstate, cone, U, tmp_cfm_2, czero, tmp_cfm)
            CALL cp_cfm_to_cfm(tmp_cfm, U)
            CALL cp_cfm_release(tmp_cfm)
            CALL cp_cfm_release(tmp_cfm_2)
         END IF

         IF (nextra > 0) THEN
            ALLOCATE (tmp_cmat(nextra, nstate))
            CALL cp_cfm_get_submatrix(U, tmp_cmat, nocc + 1, 1, nextra, nstate)
            CALL cp_cfm_set_submatrix(UL, tmp_cmat)
            DEALLOCATE (tmp_cmat)
            IF (iter > 1) THEN
               ! orthonormalize c_tilde
               CALL cp_fm_create(tmp_fm, c_tilde%matrix_struct)
               tmp_fm%local_data = REAL(c_tilde%local_data, KIND=dp)
               CALL ortho_vectors(tmp_fm)
               CALL cp_fm_to_cfm(msourcer=tmp_fm, mtarget=c_tilde)
               CALL cp_fm_release(tmp_fm)

               ALLOCATE (tmp_cmat(nocc, nstate))
               CALL cp_cfm_get_submatrix(U, tmp_cmat, 1, 1, nocc, nstate)
               CALL cp_cfm_set_submatrix(V, tmp_cmat, 1, 1, nocc, nstate)
               DEALLOCATE (tmp_cmat)
               CALL parallel_gemm("N", "N", northo, nstate, nextra, cone, c_tilde, UL, czero, VL)
               ALLOCATE (tmp_cmat(northo, nstate))
               CALL cp_cfm_get_submatrix(VL, tmp_cmat)
               CALL cp_cfm_set_submatrix(V, tmp_cmat, nocc + 1, 1, northo, nstate)
               DEALLOCATE (tmp_cmat)
            END IF

            ! reset if new_direction
            IF (new_direction .AND. MOD(line_searches, 20) .EQ. 5) THEN
               CALL cp_cfm_set_all(skc, czero)
               CALL cp_cfm_set_all(Gct_old, czero)
               norm2_old = 1.0E30_dp
            END IF

            CALL cp_cfm_create(tmp_cfm, V%matrix_struct)
            CALL cp_cfm_to_cfm(V, tmp_cfm)
            CALL cp_cfm_create(tmp_cfm_1, V%matrix_struct)
            ndummy = nmo
         ELSE
            CALL cp_cfm_create(tmp_cfm, zij(1, 1)%matrix_struct)
            CALL cp_cfm_to_cfm(U, tmp_cfm)
            CALL cp_cfm_create(tmp_cfm_1, zij(1, 1)%matrix_struct)
            ndummy = nstate
         END IF
         ! update z_ij
         DO idim = 1, dim2
            ! 'tmp_cfm_1 = zij_0*tmp_cfm'
            CALL parallel_gemm("N", "N", ndummy, nstate, ndummy, cone, zij_0(idim), &
                               tmp_cfm, czero, tmp_cfm_1)
            ! 'c_zij = tmp_cfm_dagg*tmp_cfm_1'
            CALL parallel_gemm("C", "N", nstate, nstate, ndummy, cone, tmp_cfm, tmp_cfm_1, &
                               czero, c_zij(idim))
         END DO
         CALL cp_cfm_release(tmp_cfm)
         CALL cp_cfm_release(tmp_cfm_1)
         ! compute spread
         DO istate = 1, nstate
            spread_ii = 0.0_dp
            DO idim = 1, dim2
               CALL cp_cfm_get_element(c_zij(idim), istate, istate, mzii)
               spread_ii = spread_ii + weights(idim)* &
                           ABS(mzii)**2/twopi/twopi
               matrix_zii(istate, idim) = mzii
            END DO
            !WRITE(*,*) 'spread_ii', spread_ii
            sum_spread(istate) = spread_ii
         END DO
         CALL c_zij(1)%matrix_struct%para_env%sum(spread_ii)
         spread_sum = accurate_sum(sum_spread)

         IF (nextra > 0) THEN
            ! update c_tilde
            CALL cp_cfm_set_all(zdiag, czero)
            CALL cp_cfm_set_all(grad_ctilde, czero)
            CALL cp_cfm_create(tmp_cfm, V%matrix_struct)
            CALL cp_cfm_set_all(tmp_cfm, czero)
            CALL cp_cfm_create(tmp_cfm_1, V%matrix_struct)
            CALL cp_cfm_set_all(tmp_cfm_1, czero)
            ALLOCATE (tmp_cmat(northo, nstate))
            DO idim = 1, dim2
               weight = weights(idim)
               arr_zii = matrix_zii(:, idim)
               ! tmp_cfm = zij_0*V
               CALL parallel_gemm("N", "N", nmo, nstate, nmo, cone, &
                                  zij_0(idim), V, czero, tmp_cfm)
               ! tmp_cfm = tmp_cfm*diag_zij_dagg
               CALL cp_cfm_column_scale(tmp_cfm, CONJG(arr_zii))
               ! tmp_cfm_1 = tmp_cfm*U_dagg
               CALL parallel_gemm("N", "C", nmo, nstate, nstate, cone, tmp_cfm, &
                                  U, czero, tmp_cfm_1)
               CALL cp_cfm_scale(weight, tmp_cfm_1)
               ! zdiag = zdiag + tmp_cfm_1'
               CALL cp_cfm_scale_and_add(cone, zdiag, cone, tmp_cfm_1)

               ! tmp_cfm = zij_0_dagg*V
               CALL parallel_gemm("C", "N", nmo, nstate, nmo, cone, &
                                  zij_0(idim), V, czero, tmp_cfm)

               ! tmp_cfm = tmp_cfm*diag_zij
               CALL cp_cfm_column_scale(tmp_cfm, arr_zii)
               ! tmp_cfm_1 = tmp_cfm*U_dagg
               CALL parallel_gemm("N", "C", nmo, nstate, nstate, cone, tmp_cfm, &
                                  U, czero, tmp_cfm_1)
               CALL cp_cfm_scale(weight, tmp_cfm_1)
               ! zdiag = zdiag + tmp_cfm_1'
               CALL cp_cfm_scale_and_add(cone, zdiag, cone, tmp_cfm_1)
            END DO ! idim
            CALL cp_cfm_release(tmp_cfm)
            CALL cp_cfm_release(tmp_cfm_1)
            DEALLOCATE (tmp_cmat)
            ALLOCATE (tmp_cmat(northo, nextra))
            CALL cp_cfm_get_submatrix(zdiag, tmp_cmat, nocc + 1, nocc + 1, &
                                      northo, nextra, .FALSE.)
            ! 'grad_ctilde'
            CALL cp_cfm_set_submatrix(grad_ctilde, tmp_cmat)
            DEALLOCATE (tmp_cmat)
            ! ctrans_lambda = c_tilde_dagg*grad_ctilde
            CALL parallel_gemm("C", "N", nextra, nextra, northo, cone, c_tilde, grad_ctilde, czero, ctrans_lambda)
            !WRITE(*,*) "norm(ctrans_lambda) = ", cp_cfm_norm(ctrans_lambda, "F")
            ! 'grad_ctilde = - c_tilde*ctrans_lambda + grad_ctilde'
            CALL parallel_gemm("N", "N", northo, nextra, nextra, -cone, c_tilde, ctrans_lambda, cone, grad_ctilde)
         END IF ! nextra > 0

         ! tolerance
         IF (nextra > 0) THEN
            tolc = 0.0_dp
            CALL cp_fm_create(tmp_fm, grad_ctilde%matrix_struct)
            CALL cp_cfm_to_fm(grad_ctilde, tmp_fm)
            CALL cp_fm_maxabsval(tmp_fm, tolc)
            CALL cp_fm_release(tmp_fm)
            !WRITE(*,*) 'tolc = ', tolc
            tol = tol + tolc
         END IF
         !WRITE(*,*) 'tol = ', tol

         IF (nextra > 0) THEN
            !WRITE(*,*) 'new_direction: ', new_direction
            IF (new_direction) THEN
               line_searches = line_searches + 1
               IF (mintol > tol) THEN
                  mintol = tol
                  miniter = iter
               END IF

               IF (unit_nr > 0 .AND. MODULO(iter, out_each) == 0) THEN
                  sum_spread_ii = alpha*nstate/twopi/twopi - spread_sum
                  avg_spread_ii = sum_spread_ii/nstate
                  WRITE (unit_nr, '(T4,A,T26,A,T48,A)') &
                     "Iteration", "Avg. Spread_ii", "Tolerance"
                  WRITE (unit_nr, '(T4,I7,T20,F20.10,T45,E12.4)') &
                     iter, avg_spread_ii, tol
                  CALL m_flush(unit_nr)
               END IF
               IF (tol < eps_localization) EXIT

               IF (do_cg) THEN
                  cnorm2_Gct = czero
                  cnorm2_Gct_cross = czero
                  CALL cp_cfm_trace(grad_ctilde, Gct_old, cnorm2_Gct_cross)
                  norm2_Gct_cross = REAL(cnorm2_Gct_cross, KIND=dp)
                  Gct_old%local_data = grad_ctilde%local_data
                  CALL cp_cfm_trace(grad_ctilde, Gct_old, cnorm2_Gct)
                  norm2_Gct = REAL(cnorm2_Gct, KIND=dp)
                  ! compute beta_pr
                  beta_pr = (norm2_Gct - norm2_Gct_cross)/norm2_old
                  norm2_old = norm2_Gct
                  beta = MAX(0.0_dp, beta_pr)
                  !WRITE(*,*) 'beta = ', beta
                  ! compute skc / ska = beta * skc / ska + grad_ctilde / G
                  CALL cp_cfm_scale(beta, skc)
                  CALL cp_cfm_scale_and_add(cone, skc, cone, Gct_old)
                  CALL cp_cfm_trace(skc, Gct_old, cnorm2_Gct_cross)
                  norm2_Gct_cross = REAL(cnorm2_Gct_cross, KIND=dp)
                  IF (norm2_Gct_cross .LE. 0.0_dp) THEN ! back to steepest ascent
                     CALL cp_cfm_scale_and_add(czero, skc, cone, Gct_old)
                  END IF
               ELSE
                  CALL cp_cfm_scale_and_add(czero, skc, cone, grad_ctilde)
               END IF
               line_search_count = 0
            END IF

            line_search_count = line_search_count + 1
            !WRITE(*,*) 'line_search_count = ', line_search_count
            energy(line_search_count) = spread_sum

            ! gold line search
            new_direction = .FALSE.
            IF (line_search_count .EQ. 1) THEN
               lsl = 1
               lsr = 0
               lsm = 1
               pos(1) = 0.0_dp
               pos(2) = ds_min/gold_sec
               ds = pos(2)
            ELSE
               IF (line_search_count .EQ. 50) THEN
                  IF (ABS(energy(line_search_count) - energy(line_search_count - 1)) < 1.0E-4_dp) THEN
                     CPWARN("Line search failed to converge properly")
                     ds_min = 0.1_dp
                     new_direction = .TRUE.
                     ds = pos(line_search_count)
                     line_search_count = 0
                  ELSE
                     CPABORT("No. of line searches exceeds 50")
                  END IF
               ELSE
                  IF (lsr .EQ. 0) THEN
                     IF (energy(line_search_count - 1) .GT. energy(line_search_count)) THEN
                        lsr = line_search_count
                        pos(line_search_count + 1) = pos(lsm) + (pos(lsr) - pos(lsm))*gold_sec
                     ELSE
                        lsl = lsm
                        lsm = line_search_count
                        pos(line_search_count + 1) = pos(line_search_count)/gold_sec
                     END IF
                  ELSE
                     IF (pos(line_search_count) .LT. pos(lsm)) THEN
                        IF (energy(line_search_count) .GT. energy(lsm)) THEN
                           lsr = lsm
                           lsm = line_search_count
                        ELSE
                           lsl = line_search_count
                        END IF
                     ELSE
                        IF (energy(line_search_count) .GT. energy(lsm)) THEN
                           lsl = lsm
                           lsm = line_search_count
                        ELSE
                           lsr = line_search_count
                        END IF
                     END IF
                     IF (pos(lsr) - pos(lsm) .GT. pos(lsm) - pos(lsl)) THEN
                        pos(line_search_count + 1) = pos(lsm) + gold_sec*(pos(lsr) - pos(lsm))
                     ELSE
                        pos(line_search_count + 1) = pos(lsl) + gold_sec*(pos(lsm) - pos(lsl))
                     END IF
                     IF ((pos(lsr) - pos(lsl)) .LT. 1.0E-3_dp*pos(lsr)) THEN
                        new_direction = .TRUE.
                     END IF
                  END IF ! lsr .eq. 0
               END IF ! line_search_count .eq. 50
               ! now go to the suggested point
               ds = pos(line_search_count + 1) - pos(line_search_count)
               !WRITE(*,*) 'lsl, lsr, lsm, ds = ', lsl, lsr, lsm, ds
               IF ((ABS(ds) < 1.0E-10_dp) .AND. (lsl == 1)) THEN
                  new_direction = .TRUE.
                  ds_min = 0.5_dp/alpha
               ELSEIF (ABS(ds) > 10.0_dp) THEN
                  new_direction = .TRUE.
                  ds_min = 0.5_dp/alpha
               ELSE
                  ds_min = pos(line_search_count + 1)
               END IF
            END IF ! first step
            ! 'c_tilde = c_tilde + d*skc'
            CALL cp_cfm_scale(ds, skc)
            CALL cp_cfm_scale_and_add(cone, c_tilde, cone, skc)
         ELSE
            IF (mintol > tol) THEN
               mintol = tol
               miniter = iter
            END IF
            IF (unit_nr > 0 .AND. MODULO(iter, out_each) == 0) THEN
               sum_spread_ii = alpha*nstate/twopi/twopi - spread_sum
               avg_spread_ii = sum_spread_ii/nstate
               WRITE (unit_nr, '(T4,A,T26,A,T48,A)') &
                  "Iteration", "Avg. Spread_ii", "Tolerance"
               WRITE (unit_nr, '(T4,I7,T20,F20.10,T45,E12.4)') &
                  iter, avg_spread_ii, tol
               CALL m_flush(unit_nr)
            END IF
            IF (tol < eps_localization) EXIT
         END IF ! nextra > 0

      END DO ! iteration

      IF ((unit_nr > 0) .AND. (iter == max_iter)) THEN
         WRITE (unit_nr, '(T4,A,T4,A)') "Min. Itr.", "Min. Tol."
         WRITE (unit_nr, '(T4,I7,T4,E12.4)') miniter, mintol
         CALL m_flush(unit_nr)
      END IF

      CALL cp_cfm_to_fm(U, matrix_U)

      IF (nextra > 0) THEN
         rmat%local_data = REAL(V%local_data, KIND=dp)
         CALL rotate_orbitals_edf(rmat, vectors_all, vectors)

         CALL cp_cfm_release(c_tilde)
         CALL cp_cfm_release(grad_ctilde)
         CALL cp_cfm_release(Gct_old)
         CALL cp_cfm_release(skc)
         CALL cp_cfm_release(UL)
         CALL cp_cfm_release(zdiag)
         CALL cp_cfm_release(ctrans_lambda)
         CALL cp_fm_release(id_nextra)
         CALL cp_fm_release(vectors_all)
         CALL cp_cfm_release(V)
         CALL cp_fm_release(matrix_V)
         CALL cp_fm_release(matrix_V_all)
         CALL cp_cfm_release(VL)
         DEALLOCATE (arr_zii)
      ELSE
         rmat%local_data = matrix_U%local_data
         CALL rotate_orbitals(rmat, vectors)
      END IF
      DO idim = 1, dim2
         CALL cp_cfm_release(zij_0(idim))
      END DO
      DEALLOCATE (zij_0)

      DO idim = 1, dim2
         zij(1, idim)%local_data = REAL(c_zij(idim)%local_data, dp)
         zij(2, idim)%local_data = AIMAG(c_zij(idim)%local_data)
         CALL cp_cfm_release(c_zij(idim))
      END DO
      DEALLOCATE (c_zij)
      CALL cp_fm_release(rmat)
      CALL cp_cfm_release(U)
      CALL cp_fm_release(matrix_U)
      DEALLOCATE (matrix_zii, sum_spread)

      CALL timestop(handle)

   END SUBROUTINE jacobi_cg_edf_ls

! **************************************************************************************************
!> \brief ...
!> \param vmatrix ...
! **************************************************************************************************
   SUBROUTINE ortho_vectors(vmatrix)

      TYPE(cp_fm_type), INTENT(IN)                       :: vmatrix

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'ortho_vectors'

      INTEGER                                            :: handle, n, ncol
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(cp_fm_type)                                   :: overlap_vv

      CALL timeset(routineN, handle)

      NULLIFY (fm_struct_tmp)

      CALL cp_fm_get_info(matrix=vmatrix, nrow_global=n, ncol_global=ncol)

      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ncol, ncol_global=ncol, &
                               para_env=vmatrix%matrix_struct%para_env, &
                               context=vmatrix%matrix_struct%context)
      CALL cp_fm_create(overlap_vv, fm_struct_tmp, "overlap_vv")
      CALL cp_fm_struct_release(fm_struct_tmp)

      CALL parallel_gemm('T', 'N', ncol, ncol, n, 1.0_dp, vmatrix, vmatrix, 0.0_dp, overlap_vv)
      CALL cp_fm_cholesky_decompose(overlap_vv)
      CALL cp_fm_triangular_multiply(overlap_vv, vmatrix, n_cols=ncol, side='R', invert_tr=.TRUE.)

      CALL cp_fm_release(overlap_vv)

      CALL timestop(handle)

   END SUBROUTINE ortho_vectors

! **************************************************************************************************
!> \brief ...
!> \param istate ...
!> \param jstate ...
!> \param st ...
!> \param ct ...
!> \param zij ...
! **************************************************************************************************
   SUBROUTINE rotate_zij(istate, jstate, st, ct, zij)
      INTEGER, INTENT(IN)                                :: istate, jstate
      REAL(KIND=dp), INTENT(IN)                          :: st, ct
      TYPE(cp_cfm_type)                                  :: zij(:)

      INTEGER                                            :: id

! Locals

      DO id = 1, SIZE(zij, 1)
         CALL cp_cfm_rot_cols(zij(id), istate, jstate, ct, st)
         CALL cp_cfm_rot_rows(zij(id), istate, jstate, ct, st)
      END DO

   END SUBROUTINE rotate_zij
! **************************************************************************************************
!> \brief ...
!> \param istate ...
!> \param jstate ...
!> \param st ...
!> \param ct ...
!> \param rmat ...
! **************************************************************************************************
   SUBROUTINE rotate_rmat(istate, jstate, st, ct, rmat)
      INTEGER, INTENT(IN)                                :: istate, jstate
      REAL(KIND=dp), INTENT(IN)                          :: st, ct
      TYPE(cp_cfm_type), INTENT(IN)                      :: rmat

      CALL cp_cfm_rot_cols(rmat, istate, jstate, ct, st)

   END SUBROUTINE rotate_rmat
! **************************************************************************************************
!> \brief ...
!> \param mii ...
!> \param mjj ...
!> \param mij ...
!> \param weights ...
!> \param theta ...
!> \param grad_ij ...
!> \param step ...
! **************************************************************************************************
   SUBROUTINE get_angle(mii, mjj, mij, weights, theta, grad_ij, step)
      COMPLEX(KIND=dp), POINTER                          :: mii(:), mjj(:), mij(:)
      REAL(KIND=dp), INTENT(IN)                          :: weights(:)
      REAL(KIND=dp), INTENT(OUT)                         :: theta
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: grad_ij, step

      COMPLEX(KIND=dp)                                   :: z11, z12, z22
      INTEGER                                            :: dim_m, idim
      REAL(KIND=dp)                                      :: a12, b12, d2, ratio

      a12 = 0.0_dp
      b12 = 0.0_dp
      dim_m = SIZE(mii)
      DO idim = 1, dim_m
         z11 = mii(idim)
         z22 = mjj(idim)
         z12 = mij(idim)
         a12 = a12 + weights(idim)*REAL(CONJG(z12)*(z11 - z22), KIND=dp)
         b12 = b12 + weights(idim)*REAL((z12*CONJG(z12) - &
                                         0.25_dp*(z11 - z22)*(CONJG(z11) - CONJG(z22))), KIND=dp)
      END DO
      IF (ABS(b12) > 1.e-10_dp) THEN
         ratio = -a12/b12
         theta = 0.25_dp*ATAN(ratio)
      ELSEIF (ABS(b12) < 1.e-10_dp) THEN
         b12 = 0.0_dp
         theta = 0.0_dp
      ELSE
         theta = 0.25_dp*pi
      END IF
      IF (PRESENT(grad_ij)) theta = theta + step*grad_ij
! Check second derivative info
      d2 = a12*SIN(4._dp*theta) - b12*COS(4._dp*theta)
      IF (d2 <= 0._dp) THEN ! go to the maximum, not the minimum
         IF (theta > 0.0_dp) THEN ! make theta as small as possible
            theta = theta - 0.25_dp*pi
         ELSE
            theta = theta + 0.25_dp*pi
         END IF
      END IF
   END SUBROUTINE get_angle
! **************************************************************************************************
!> \brief ...
!> \param zij ...
!> \param weights ...
!> \param tolerance ...
!> \param grad ...
! **************************************************************************************************
   SUBROUTINE check_tolerance(zij, weights, tolerance, grad)
      TYPE(cp_cfm_type)                                  :: zij(:)
      REAL(KIND=dp), INTENT(IN)                          :: weights(:)
      REAL(KIND=dp), INTENT(OUT)                         :: tolerance
      TYPE(cp_fm_type), INTENT(OUT), OPTIONAL            :: grad

      CHARACTER(len=*), PARAMETER                        :: routineN = 'check_tolerance'

      INTEGER                                            :: handle
      TYPE(cp_fm_type)                                   :: force

      CALL timeset(routineN, handle)

! compute gradient at t=0

      CALL cp_fm_create(force, zij(1)%matrix_struct)
      CALL cp_fm_set_all(force, 0._dp)
      CALL grad_at_0(zij, weights, force)
      CALL cp_fm_maxabsval(force, tolerance)
      IF (PRESENT(grad)) CALL cp_fm_to_fm(force, grad)
      CALL cp_fm_release(force)

      CALL timestop(handle)

   END SUBROUTINE check_tolerance
! **************************************************************************************************
!> \brief ...
!> \param rmat ...
!> \param vectors ...
! **************************************************************************************************
   SUBROUTINE rotate_orbitals(rmat, vectors)
      TYPE(cp_fm_type), INTENT(IN)                       :: rmat, vectors

      INTEGER                                            :: k, n
      TYPE(cp_fm_type)                                   :: wf

      CALL cp_fm_create(wf, vectors%matrix_struct)
      CALL cp_fm_get_info(vectors, nrow_global=n, ncol_global=k)
      CALL parallel_gemm("N", "N", n, k, k, 1.0_dp, vectors, rmat, 0.0_dp, wf)
      CALL cp_fm_to_fm(wf, vectors)
      CALL cp_fm_release(wf)
   END SUBROUTINE rotate_orbitals
! **************************************************************************************************
!> \brief ...
!> \param rmat ...
!> \param vec_all ...
!> \param vectors ...
! **************************************************************************************************
   SUBROUTINE rotate_orbitals_edf(rmat, vec_all, vectors)
      TYPE(cp_fm_type), INTENT(IN)                       :: rmat, vec_all, vectors

      INTEGER                                            :: k, l, n
      TYPE(cp_fm_type)                                   :: wf

      CALL cp_fm_create(wf, vectors%matrix_struct)
      CALL cp_fm_get_info(vec_all, nrow_global=n, ncol_global=k)
      CALL cp_fm_get_info(rmat, ncol_global=l)

      CALL parallel_gemm("N", "N", n, l, k, 1.0_dp, vec_all, rmat, 0.0_dp, wf)
      CALL cp_fm_to_fm(wf, vectors)
      CALL cp_fm_release(wf)
   END SUBROUTINE rotate_orbitals_edf
! **************************************************************************************************
!> \brief ...
!> \param diag ...
!> \param weights ...
!> \param matrix ...
!> \param ndim ...
! **************************************************************************************************
   SUBROUTINE gradsq_at_0(diag, weights, matrix, ndim)
      COMPLEX(KIND=dp), DIMENSION(:, :), POINTER         :: diag
      REAL(KIND=dp), INTENT(IN)                          :: weights(:)
      TYPE(cp_fm_type), INTENT(IN)                       :: matrix
      INTEGER, INTENT(IN)                                :: ndim

      COMPLEX(KIND=dp)                                   :: zii, zjj
      INTEGER                                            :: idim, istate, jstate, ncol_local, &
                                                            nrow_global, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: gradsq_ij

      CALL cp_fm_get_info(matrix, nrow_local=nrow_local, &
                          ncol_local=ncol_local, nrow_global=nrow_global, &
                          row_indices=row_indices, col_indices=col_indices)

      DO istate = 1, nrow_local
         DO jstate = 1, ncol_local
! get real and imaginary parts
            gradsq_ij = 0.0_dp
            DO idim = 1, ndim
               zii = diag(row_indices(istate), idim)
               zjj = diag(col_indices(jstate), idim)
               gradsq_ij = gradsq_ij + weights(idim)* &
                           4.0_dp*REAL((CONJG(zii)*zii + CONJG(zjj)*zjj), KIND=dp)
            END DO
            matrix%local_data(istate, jstate) = gradsq_ij
         END DO
      END DO
   END SUBROUTINE gradsq_at_0
! **************************************************************************************************
!> \brief ...
!> \param matrix_p ...
!> \param weights ...
!> \param matrix ...
! **************************************************************************************************
   SUBROUTINE grad_at_0(matrix_p, weights, matrix)
      TYPE(cp_cfm_type)                                  :: matrix_p(:)
      REAL(KIND=dp), INTENT(IN)                          :: weights(:)
      TYPE(cp_fm_type), INTENT(IN)                       :: matrix

      COMPLEX(KIND=dp)                                   :: zii, zij, zjj
      COMPLEX(KIND=dp), DIMENSION(:, :), POINTER         :: diag
      INTEGER                                            :: dim_m, idim, istate, jstate, ncol_local, &
                                                            nrow_global, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: grad_ij

      NULLIFY (diag)
      CALL cp_fm_get_info(matrix, nrow_local=nrow_local, &
                          ncol_local=ncol_local, nrow_global=nrow_global, &
                          row_indices=row_indices, col_indices=col_indices)
      dim_m = SIZE(matrix_p, 1)
      ALLOCATE (diag(nrow_global, dim_m))

      DO idim = 1, dim_m
         DO istate = 1, nrow_global
            CALL cp_cfm_get_element(matrix_p(idim), istate, istate, diag(istate, idim))
         END DO
      END DO

      DO istate = 1, nrow_local
         DO jstate = 1, ncol_local
! get real and imaginary parts
            grad_ij = 0.0_dp
            DO idim = 1, dim_m
               zii = diag(row_indices(istate), idim)
               zjj = diag(col_indices(jstate), idim)
               zij = matrix_p(idim)%local_data(istate, jstate)
               grad_ij = grad_ij + weights(idim)* &
                         REAL(4.0_dp*CONJG(zij)*(zjj - zii), dp)
            END DO
            matrix%local_data(istate, jstate) = grad_ij
         END DO
      END DO
      DEALLOCATE (diag)
   END SUBROUTINE grad_at_0

! return energy and maximum gradient in the current point
! **************************************************************************************************
!> \brief ...
!> \param weights ...
!> \param zij ...
!> \param tolerance ...
!> \param value ...
! **************************************************************************************************
   SUBROUTINE check_tolerance_new(weights, zij, tolerance, value)
      REAL(KIND=dp), INTENT(IN)                          :: weights(:)
      TYPE(cp_fm_type), INTENT(IN)                       :: zij(:, :)
      REAL(KIND=dp)                                      :: tolerance, value

      COMPLEX(KIND=dp)                                   :: kii, kij, kjj
      COMPLEX(KIND=dp), DIMENSION(:, :), POINTER         :: diag
      INTEGER                                            :: idim, istate, jstate, ncol_local, &
                                                            nrow_global, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: grad_ij, ra, rb

      NULLIFY (diag)
      CALL cp_fm_get_info(zij(1, 1), nrow_local=nrow_local, &
                          ncol_local=ncol_local, nrow_global=nrow_global, &
                          row_indices=row_indices, col_indices=col_indices)
      ALLOCATE (diag(nrow_global, SIZE(zij, 2)))
      value = 0.0_dp
      DO idim = 1, SIZE(zij, 2)
         DO istate = 1, nrow_global
            CALL cp_fm_get_element(zij(1, idim), istate, istate, ra)
            CALL cp_fm_get_element(zij(2, idim), istate, istate, rb)
            diag(istate, idim) = CMPLX(ra, rb, dp)
            value = value + weights(idim) - weights(idim)*ABS(diag(istate, idim))**2
         END DO
      END DO
      tolerance = 0.0_dp
      DO istate = 1, nrow_local
         DO jstate = 1, ncol_local
            grad_ij = 0.0_dp
            DO idim = 1, SIZE(zij, 2)
               kii = diag(row_indices(istate), idim)
               kjj = diag(col_indices(jstate), idim)
               ra = zij(1, idim)%local_data(istate, jstate)
               rb = zij(2, idim)%local_data(istate, jstate)
               kij = CMPLX(ra, rb, dp)
               grad_ij = grad_ij + weights(idim)* &
                         REAL(4.0_dp*CONJG(kij)*(kjj - kii), dp)
            END DO
            tolerance = MAX(ABS(grad_ij), tolerance)
         END DO
      END DO
      CALL zij(1, 1)%matrix_struct%para_env%max(tolerance)

      DEALLOCATE (diag)

   END SUBROUTINE check_tolerance_new

! **************************************************************************************************
!> \brief yet another crazy try, computes the angles needed to rotate the orbitals first
!>        and rotates them all at the same time (hoping for the best of course)
!> \param weights ...
!> \param zij ...
!> \param vectors ...
!> \param max_iter ...
!> \param max_crazy_angle ...
!> \param crazy_scale ...
!> \param crazy_use_diag ...
!> \param eps_localization ...
!> \param iterations ...
!> \param converged ...
! **************************************************************************************************
   SUBROUTINE crazy_rotations(weights, zij, vectors, max_iter, max_crazy_angle, crazy_scale, crazy_use_diag, &
                              eps_localization, iterations, converged)
      REAL(KIND=dp), INTENT(IN)                          :: weights(:)
      TYPE(cp_fm_type), INTENT(IN)                       :: zij(:, :), vectors
      INTEGER, INTENT(IN)                                :: max_iter
      REAL(KIND=dp), INTENT(IN)                          :: max_crazy_angle
      REAL(KIND=dp)                                      :: crazy_scale
      LOGICAL                                            :: crazy_use_diag
      REAL(KIND=dp), INTENT(IN)                          :: eps_localization
      INTEGER                                            :: iterations
      LOGICAL, INTENT(out), OPTIONAL                     :: converged

      CHARACTER(len=*), PARAMETER                        :: routineN = 'crazy_rotations'
      COMPLEX(KIND=dp), PARAMETER                        :: cone = (1.0_dp, 0.0_dp), &
                                                            czero = (0.0_dp, 0.0_dp)

      COMPLEX(KIND=dp), DIMENSION(:), POINTER            :: evals_exp
      COMPLEX(KIND=dp), DIMENSION(:, :), POINTER         :: diag_z
      COMPLEX(KIND=dp), POINTER                          :: mii(:), mij(:), mjj(:)
      INTEGER                                            :: dim2, handle, i, icol, idim, irow, &
                                                            method, ncol_global, ncol_local, &
                                                            norder, nrow_global, nrow_local, &
                                                            nsquare, unit_nr
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: do_emd
      REAL(KIND=dp)                                      :: eps_exp, limit_crazy_angle, maxeval, &
                                                            norm, ra, rb, theta, tolerance, value
      REAL(KIND=dp), DIMENSION(:), POINTER               :: evals
      TYPE(cp_cfm_type)                                  :: cmat_A, cmat_R, cmat_t1
      TYPE(cp_fm_type)                                   :: mat_R, mat_t, mat_theta, mat_U

      CALL timeset(routineN, handle)
      NULLIFY (row_indices, col_indices)
      CALL cp_fm_get_info(zij(1, 1), nrow_global=nrow_global, &
                          ncol_global=ncol_global, &
                          row_indices=row_indices, col_indices=col_indices, &
                          nrow_local=nrow_local, ncol_local=ncol_local)

      limit_crazy_angle = max_crazy_angle

      NULLIFY (diag_z, evals, evals_exp, mii, mij, mjj)
      dim2 = SIZE(zij, 2)
      ALLOCATE (diag_z(nrow_global, dim2))
      ALLOCATE (evals(nrow_global))
      ALLOCATE (evals_exp(nrow_global))

      CALL cp_cfm_create(cmat_A, zij(1, 1)%matrix_struct)
      CALL cp_cfm_create(cmat_R, zij(1, 1)%matrix_struct)
      CALL cp_cfm_create(cmat_t1, zij(1, 1)%matrix_struct)

      CALL cp_fm_create(mat_U, zij(1, 1)%matrix_struct)
      CALL cp_fm_create(mat_t, zij(1, 1)%matrix_struct)
      CALL cp_fm_create(mat_R, zij(1, 1)%matrix_struct)

      CALL cp_fm_create(mat_theta, zij(1, 1)%matrix_struct)

      CALL cp_fm_set_all(mat_R, 0.0_dp, 1.0_dp)
      CALL cp_fm_set_all(mat_t, 0.0_dp)
      ALLOCATE (mii(dim2), mij(dim2), mjj(dim2))
      DO idim = 1, dim2
         CALL cp_fm_scale_and_add(1.0_dp, mat_t, weights(idim), zij(1, idim))
         CALL cp_fm_scale_and_add(1.0_dp, mat_t, weights(idim), zij(2, idim))
      END DO
      CALL cp_fm_syevd(mat_t, mat_U, evals)
      DO idim = 1, dim2
         ! rotate z's
         CALL parallel_gemm('N', 'N', nrow_global, nrow_global, nrow_global, 1.0_dp, zij(1, idim), mat_U, 0.0_dp, mat_t)
         CALL parallel_gemm('T', 'N', nrow_global, nrow_global, nrow_global, 1.0_dp, mat_U, mat_t, 0.0_dp, zij(1, idim))
         CALL parallel_gemm('N', 'N', nrow_global, nrow_global, nrow_global, 1.0_dp, zij(2, idim), mat_U, 0.0_dp, mat_t)
         CALL parallel_gemm('T', 'N', nrow_global, nrow_global, nrow_global, 1.0_dp, mat_U, mat_t, 0.0_dp, zij(2, idim))
      END DO
      ! collect rotation matrix
      CALL parallel_gemm('N', 'N', nrow_global, nrow_global, nrow_global, 1.0_dp, mat_R, mat_U, 0.0_dp, mat_t)
      CALL cp_fm_to_fm(mat_t, mat_R)

      unit_nr = -1
      IF (cmat_A%matrix_struct%para_env%is_source()) THEN
         unit_nr = cp_logger_get_default_unit_nr()
         WRITE (unit_nr, '(T2,A7,A6,1X,A20,A12,A12,A12)') &
            "CRAZY| ", "Iter", "value    ", "gradient", "Max. eval", "limit"
      END IF

      iterations = 0
      tolerance = 1.0_dp

      DO
         iterations = iterations + 1
         DO idim = 1, dim2
            DO i = 1, nrow_global
               CALL cp_fm_get_element(zij(1, idim), i, i, ra)
               CALL cp_fm_get_element(zij(2, idim), i, i, rb)
               diag_z(i, idim) = CMPLX(ra, rb, dp)
            END DO
         END DO
         DO irow = 1, nrow_local
            DO icol = 1, ncol_local
               DO idim = 1, dim2
                  ra = zij(1, idim)%local_data(irow, icol)
                  rb = zij(2, idim)%local_data(irow, icol)
                  mij(idim) = CMPLX(ra, rb, dp)
                  mii(idim) = diag_z(row_indices(irow), idim)
                  mjj(idim) = diag_z(col_indices(icol), idim)
               END DO
               IF (row_indices(irow) .NE. col_indices(icol)) THEN
                  CALL get_angle(mii, mjj, mij, weights, theta)
                  theta = crazy_scale*theta
                  IF (theta .GT. limit_crazy_angle) theta = limit_crazy_angle
                  IF (theta .LT. -limit_crazy_angle) theta = -limit_crazy_angle
                  IF (crazy_use_diag) THEN
                     cmat_A%local_data(irow, icol) = -CMPLX(0.0_dp, theta, dp)
                  ELSE
                     mat_theta%local_data(irow, icol) = -theta
                  END IF
               ELSE
                  IF (crazy_use_diag) THEN
                     cmat_A%local_data(irow, icol) = czero
                  ELSE
                     mat_theta%local_data(irow, icol) = 0.0_dp
                  END IF
               END IF
            END DO
         END DO

         ! construct rotation matrix U based on A using diagonalization
         ! alternatively, exp based on repeated squaring could be faster
         IF (crazy_use_diag) THEN
            CALL cp_cfm_heevd(cmat_A, cmat_R, evals)
            maxeval = MAXVAL(ABS(evals))
            evals_exp(:) = EXP((0.0_dp, -1.0_dp)*evals(:))
            CALL cp_cfm_to_cfm(cmat_R, cmat_t1)
            CALL cp_cfm_column_scale(cmat_t1, evals_exp)
            CALL parallel_gemm('N', 'C', nrow_global, nrow_global, nrow_global, cone, &
                               cmat_t1, cmat_R, czero, cmat_A)
            mat_U%local_data = REAL(cmat_A%local_data, KIND=dp) ! U is a real matrix
         ELSE
            do_emd = .FALSE.
            method = 2
            eps_exp = 1.0_dp*EPSILON(eps_exp)
            CALL cp_fm_maxabsrownorm(mat_theta, norm)
            maxeval = norm ! an upper bound
            CALL get_nsquare_norder(norm, nsquare, norder, eps_exp, method, do_emd)
            CALL exp_pade_real(mat_U, mat_theta, nsquare, norder)
         END IF

         DO idim = 1, dim2
            ! rotate z's
            CALL parallel_gemm('N', 'N', nrow_global, nrow_global, nrow_global, 1.0_dp, zij(1, idim), mat_U, 0.0_dp, mat_t)
            CALL parallel_gemm('T', 'N', nrow_global, nrow_global, nrow_global, 1.0_dp, mat_U, mat_t, 0.0_dp, zij(1, idim))
            CALL parallel_gemm('N', 'N', nrow_global, nrow_global, nrow_global, 1.0_dp, zij(2, idim), mat_U, 0.0_dp, mat_t)
            CALL parallel_gemm('T', 'N', nrow_global, nrow_global, nrow_global, 1.0_dp, mat_U, mat_t, 0.0_dp, zij(2, idim))
         END DO
         ! collect rotation matrix
         CALL parallel_gemm('N', 'N', nrow_global, nrow_global, nrow_global, 1.0_dp, mat_R, mat_U, 0.0_dp, mat_t)
         CALL cp_fm_to_fm(mat_t, mat_R)

         CALL check_tolerance_new(weights, zij, tolerance, value)

         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A7,I6,1X,G20.15,E12.4,E12.4,E12.4)') &
               "CRAZY| ", iterations, value, tolerance, maxeval, limit_crazy_angle
            CALL m_flush(unit_nr)
         END IF
         IF (tolerance .LT. eps_localization .OR. iterations .GE. max_iter) EXIT
      END DO

      IF (PRESENT(converged)) converged = (tolerance .LT. eps_localization)

      CALL cp_cfm_release(cmat_A)
      CALL cp_cfm_release(cmat_R)
      CALL cp_cfm_release(cmat_T1)

      CALL cp_fm_release(mat_U)
      CALL cp_fm_release(mat_T)
      CALL cp_fm_release(mat_theta)

      CALL rotate_orbitals(mat_R, vectors)

      CALL cp_fm_release(mat_R)
      DEALLOCATE (evals_exp, evals, diag_z)
      DEALLOCATE (mii, mij, mjj)

      CALL timestop(handle)

   END SUBROUTINE crazy_rotations

! **************************************************************************************************
!> \brief use the exponential parametrization as described in to perform a direct mini
!>        Gerd Berghold et al. PRB 61 (15), pag. 10040 (2000)
!> none of the input is modified for the time being, just finds the rotations
!> that minimizes, and throws it away afterwards :-)
!> apart from being expensive and not cleaned, this works fine
!> useful to try different spread functionals
!> \param weights ...
!> \param zij ...
!> \param vectors ...
!> \param max_iter ...
!> \param eps_localization ...
!> \param iterations ...
! **************************************************************************************************
   SUBROUTINE direct_mini(weights, zij, vectors, max_iter, eps_localization, iterations)
      REAL(KIND=dp), INTENT(IN)                          :: weights(:)
      TYPE(cp_fm_type), INTENT(IN)                       :: zij(:, :), vectors
      INTEGER, INTENT(IN)                                :: max_iter
      REAL(KIND=dp), INTENT(IN)                          :: eps_localization
      INTEGER                                            :: iterations

      CHARACTER(len=*), PARAMETER                        :: routineN = 'direct_mini'
      COMPLEX(KIND=dp), PARAMETER                        :: cone = (1.0_dp, 0.0_dp), &
                                                            czero = (0.0_dp, 0.0_dp)
      REAL(KIND=dp), PARAMETER                           :: gold_sec = 0.3819_dp

      COMPLEX(KIND=dp)                                   :: lk, ll, tmp
      COMPLEX(KIND=dp), DIMENSION(:), POINTER            :: evals_exp
      COMPLEX(KIND=dp), DIMENSION(:, :), POINTER         :: diag_z
      INTEGER                                            :: handle, i, icol, idim, irow, &
                                                            line_search_count, line_searches, lsl, &
                                                            lsm, lsr, n, ncol_local, ndim, &
                                                            nrow_local, output_unit
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: new_direction
      REAL(KIND=dp)                                      :: a, b, beta_pr, c, denom, ds, ds_min, fa, &
                                                            fb, fc, nom, normg, normg_cross, &
                                                            normg_old, npos, omega, tol, val, x0, &
                                                            x1, xa, xb, xc
      REAL(KIND=dp), DIMENSION(150)                      :: energy, grad, pos
      REAL(KIND=dp), DIMENSION(:), POINTER               :: evals, fval, fvald
      TYPE(cp_cfm_type)                                  :: cmat_A, cmat_B, cmat_M, cmat_R, cmat_t1, &
                                                            cmat_t2, cmat_U
      TYPE(cp_cfm_type), ALLOCATABLE, DIMENSION(:)       :: c_zij
      TYPE(cp_fm_type)                                   :: matrix_A, matrix_G, matrix_G_old, &
                                                            matrix_G_search, matrix_H, matrix_R, &
                                                            matrix_T

      NULLIFY (evals, evals_exp, diag_z, fval, fvald)

      CALL timeset(routineN, handle)
      output_unit = cp_logger_get_default_io_unit()

      n = zij(1, 1)%matrix_struct%nrow_global
      ndim = (SIZE(zij, 2))

      IF (output_unit > 0) THEN
         WRITE (output_unit, '(T4,A )') "Localization by direct minimization of the functional; "
         WRITE (output_unit, '(T5,2A13,A20,A20,A10 )') " Line search ", " Iteration ", " Functional ", " Tolerance ", " ds Min "
      END IF

      ALLOCATE (evals(n), evals_exp(n), diag_z(n, ndim), fval(n), fvald(n))
      ALLOCATE (c_zij(ndim))

      ! create the three complex matrices Z
      DO idim = 1, ndim
         CALL cp_cfm_create(c_zij(idim), zij(1, 1)%matrix_struct)
         c_zij(idim)%local_data = CMPLX(zij(1, idim)%local_data, &
                                        zij(2, idim)%local_data, dp)
      END DO

      CALL cp_fm_create(matrix_A, zij(1, 1)%matrix_struct)
      CALL cp_fm_create(matrix_G, zij(1, 1)%matrix_struct)
      CALL cp_fm_create(matrix_T, zij(1, 1)%matrix_struct)
      CALL cp_fm_create(matrix_H, zij(1, 1)%matrix_struct)
      CALL cp_fm_create(matrix_G_search, zij(1, 1)%matrix_struct)
      CALL cp_fm_create(matrix_G_old, zij(1, 1)%matrix_struct)
      CALL cp_fm_create(matrix_R, zij(1, 1)%matrix_struct)
      CALL cp_fm_set_all(matrix_R, 0.0_dp, 1.0_dp)

      CALL cp_fm_set_all(matrix_A, 0.0_dp)
!    CALL cp_fm_init_random ( matrix_A )

      CALL cp_cfm_create(cmat_A, zij(1, 1)%matrix_struct)
      CALL cp_cfm_create(cmat_U, zij(1, 1)%matrix_struct)
      CALL cp_cfm_create(cmat_R, zij(1, 1)%matrix_struct)
      CALL cp_cfm_create(cmat_t1, zij(1, 1)%matrix_struct)
      CALL cp_cfm_create(cmat_t2, zij(1, 1)%matrix_struct)
      CALL cp_cfm_create(cmat_B, zij(1, 1)%matrix_struct)
      CALL cp_cfm_create(cmat_M, zij(1, 1)%matrix_struct)

      CALL cp_cfm_get_info(cmat_B, nrow_local=nrow_local, ncol_local=ncol_local, &
                           row_indices=row_indices, col_indices=col_indices)

      CALL cp_fm_set_all(matrix_G_old, 0.0_dp)
      CALL cp_fm_set_all(matrix_G_search, 0.0_dp)
      normg_old = 1.0E30_dp
      ds_min = 1.0_dp
      new_direction = .TRUE.
      Iterations = 0
      line_searches = 0
      line_search_count = 0
      DO
         iterations = iterations + 1
         ! compute U,R,evals given A
         cmat_A%local_data = CMPLX(0.0_dp, matrix_A%local_data, dp) ! cmat_A is hermitian, evals are reals
         CALL cp_cfm_heevd(cmat_A, cmat_R, evals)
         evals_exp(:) = EXP((0.0_dp, -1.0_dp)*evals(:))
         CALL cp_cfm_to_cfm(cmat_R, cmat_t1)
         CALL cp_cfm_column_scale(cmat_t1, evals_exp)
         CALL parallel_gemm('N', 'C', n, n, n, cone, cmat_t1, cmat_R, czero, cmat_U)
         cmat_U%local_data = REAL(cmat_U%local_data, KIND=dp) ! enforce numerics, U is a real matrix

         IF (new_direction .AND. MOD(line_searches, 20) .EQ. 5) THEN ! reset with A .eq. 0
            DO idim = 1, ndim
               CALL parallel_gemm('N', 'N', n, n, n, cone, c_zij(idim), cmat_U, czero, cmat_t1)
               CALL parallel_gemm('C', 'N', n, n, n, cone, cmat_U, cmat_t1, czero, c_zij(idim))
            END DO
            ! collect rotation matrix
            matrix_H%local_data = REAL(cmat_U%local_data, KIND=dp)
            CALL parallel_gemm('N', 'N', n, n, n, 1.0_dp, matrix_R, matrix_H, 0.0_dp, matrix_T)
            CALL cp_fm_to_fm(matrix_T, matrix_R)

            CALL cp_cfm_set_all(cmat_U, czero, cone)
            CALL cp_cfm_set_all(cmat_R, czero, cone)
            CALL cp_cfm_set_all(cmat_A, czero)
            CALL cp_fm_set_all(matrix_A, 0.0_dp)
            evals(:) = 0.0_dp
            evals_exp(:) = EXP((0.0_dp, -1.0_dp)*evals(:))
            CALL cp_fm_set_all(matrix_G_old, 0.0_dp)
            CALL cp_fm_set_all(matrix_G_search, 0.0_dp)
            normg_old = 1.0E30_dp
         END IF

         ! compute Omega and M
         CALL cp_cfm_set_all(cmat_M, czero)
         omega = 0.0_dp
         DO idim = 1, ndim
            CALL parallel_gemm('N', 'N', n, n, n, cone, c_zij(idim), cmat_U, czero, cmat_t1) ! t1=ZU
            CALL parallel_gemm('C', 'N', n, n, n, cone, cmat_U, cmat_t1, czero, cmat_t2) ! t2=(U^T)ZU
            DO i = 1, n
               CALL cp_cfm_get_element(cmat_t2, i, i, diag_z(i, idim))
               SELECT CASE (2) ! allows for selection of different spread functionals
               CASE (1)
                  fval(i) = -weights(idim)*LOG(ABS(diag_z(i, idim))**2)
                  fvald(i) = -weights(idim)/(ABS(diag_z(i, idim))**2)
               CASE (2) ! corresponds to the jacobi setup
                  fval(i) = weights(idim) - weights(idim)*ABS(diag_z(i, idim))**2
                  fvald(i) = -weights(idim)
               END SELECT
               omega = omega + fval(i)
            END DO
            DO icol = 1, ncol_local
               DO irow = 1, nrow_local
                  tmp = cmat_t1%local_data(irow, icol)*CONJG(diag_z(col_indices(icol), idim))
                  cmat_M%local_data(irow, icol) = cmat_M%local_data(irow, icol) &
                                                  + 4.0_dp*fvald(col_indices(icol))*REAL(tmp, KIND=dp)
               END DO
            END DO
         END DO

         ! compute Hessian diagonal approximation for the preconditioner
         IF (.TRUE.) THEN
            CALL gradsq_at_0(diag_z, weights, matrix_H, ndim)
         ELSE
            CALL cp_fm_set_all(matrix_H, 1.0_dp)
         END IF

         ! compute B
         DO icol = 1, ncol_local
            DO irow = 1, nrow_local
               ll = (0.0_dp, -1.0_dp)*evals(row_indices(irow))
               lk = (0.0_dp, -1.0_dp)*evals(col_indices(icol))
               IF (ABS(ll - lk) .LT. 0.5_dp) THEN ! use a series expansion to avoid loss of precision
                  tmp = 1.0_dp
                  cmat_B%local_data(irow, icol) = 0.0_dp
                  DO i = 1, 16
                     cmat_B%local_data(irow, icol) = cmat_B%local_data(irow, icol) + tmp
                     tmp = tmp*(ll - lk)/(i + 1)
                  END DO
                  cmat_B%local_data(irow, icol) = cmat_B%local_data(irow, icol)*EXP(lk)
               ELSE
                  cmat_B%local_data(irow, icol) = (EXP(lk) - EXP(ll))/(lk - ll)
               END IF
            END DO
         END DO
         ! compute gradient matrix_G

         CALL parallel_gemm('C', 'N', n, n, n, cone, cmat_M, cmat_R, czero, cmat_t1) ! t1=(M^T)(R^T)
         CALL parallel_gemm('C', 'N', n, n, n, cone, cmat_R, cmat_t1, czero, cmat_t2) ! t2=(R)t1
         CALL cp_cfm_schur_product(cmat_t2, cmat_B, cmat_t1)
         CALL parallel_gemm('N', 'C', n, n, n, cone, cmat_t1, cmat_R, czero, cmat_t2)
         CALL parallel_gemm('N', 'N', n, n, n, cone, cmat_R, cmat_t2, czero, cmat_t1)
         matrix_G%local_data = REAL(cmat_t1%local_data, KIND=dp)
         CALL cp_fm_transpose(matrix_G, matrix_T)
         CALL cp_fm_scale_and_add(-1.0_dp, matrix_G, 1.0_dp, matrix_T)
         CALL cp_fm_maxabsval(matrix_G, tol)

         ! from here on, minimizing technology
         IF (new_direction) THEN
            ! energy converged up to machine precision ?
            line_searches = line_searches + 1
            ! DO i=1,line_search_count
            !   write(15,*) pos(i),energy(i)
            ! ENDDO
            ! write(15,*) ""
            ! CALL m_flush(15)
            !write(16,*) evals(:)
            !write(17,*) matrix_A%local_data(:,:)
            !write(18,*) matrix_G%local_data(:,:)
            IF (output_unit > 0) THEN
               WRITE (output_unit, '(T5,I10,T18,I10,T31,2F20.6,F10.3)') line_searches, Iterations, Omega, tol, ds_min
               CALL m_flush(output_unit)
            END IF
            IF (tol < eps_localization .OR. iterations > max_iter) EXIT

            IF (.TRUE.) THEN ! do conjugate gradient CG
               CALL cp_fm_trace(matrix_G, matrix_G_old, normg_cross)
               normg_cross = normg_cross*0.5_dp ! takes into account the fact that A is antisymmetric
               ! apply the preconditioner
               DO icol = 1, ncol_local
                  DO irow = 1, nrow_local
                     matrix_G_old%local_data(irow, icol) = matrix_G%local_data(irow, icol)/matrix_H%local_data(irow, icol)
                  END DO
               END DO
               CALL cp_fm_trace(matrix_G, matrix_G_old, normg)
               normg = normg*0.5_dp
               beta_pr = (normg - normg_cross)/normg_old
               normg_old = normg
               beta_pr = MAX(beta_pr, 0.0_dp)
               CALL cp_fm_scale_and_add(beta_pr, matrix_G_search, -1.0_dp, matrix_G_old)
               CALL cp_fm_trace(matrix_G_search, matrix_G_old, normg_cross)
               IF (normg_cross .GE. 0) THEN ! back to SD
                  IF (matrix_A%matrix_struct%para_env%is_source()) THEN
                     WRITE (cp_logger_get_default_unit_nr(), *) "!"
                  END IF
                  beta_pr = 0.0_dp
                  CALL cp_fm_scale_and_add(beta_pr, matrix_G_search, -1.0_dp, matrix_G_old)
               END IF
            ELSE ! SD
               CALL cp_fm_scale_and_add(0.0_dp, matrix_G_search, -1.0_dp, matrix_G)
            END IF
            ! ds_min=1.0E-4_dp
            line_search_count = 0
         END IF
         line_search_count = line_search_count + 1
         energy(line_search_count) = Omega

         ! line search section
         SELECT CASE (3)
         CASE (1) ! two point line search
            SELECT CASE (line_search_count)
            CASE (1)
               pos(1) = 0.0_dp
               pos(2) = ds_min
               CALL cp_fm_trace(matrix_G, matrix_G_search, grad(1))
               grad(1) = grad(1)/2.0_dp
               new_direction = .FALSE.
            CASE (2)
               new_direction = .TRUE.
               x0 = pos(1) ! 0.0_dp
               c = energy(1)
               b = grad(1)
               x1 = pos(2)
               a = (energy(2) - b*x1 - c)/(x1**2)
               IF (a .LE. 0.0_dp) a = 1.0E-15_dp
               npos = -b/(2.0_dp*a)
               val = a*npos**2 + b*npos + c
               IF (val .LT. energy(1) .AND. val .LE. energy(2)) THEN
                  ! we go to a minimum, but ...
                  ! we take a guard against too large steps
                  pos(3) = MIN(npos, MAXVAL(pos(1:2))*4.0_dp)
               ELSE ! just take an extended step
                  pos(3) = MAXVAL(pos(1:2))*2.0_dp
               END IF
            END SELECT
         CASE (2) ! 3 point line search
            SELECT CASE (line_search_count)
            CASE (1)
               new_direction = .FALSE.
               pos(1) = 0.0_dp
               pos(2) = ds_min*0.8_dp
            CASE (2)
               new_direction = .FALSE.
               IF (energy(2) .GT. energy(1)) THEN
                  pos(3) = ds_min*0.7_dp
               ELSE
                  pos(3) = ds_min*1.4_dp
               END IF
            CASE (3)
               new_direction = .TRUE.
               xa = pos(1)
               xb = pos(2)
               xc = pos(3)
               fa = energy(1)
               fb = energy(2)
               fc = energy(3)
               nom = (xb - xa)**2*(fb - fc) - (xb - xc)**2*(fb - fa)
               denom = (xb - xa)*(fb - fc) - (xb - xc)*(fb - fa)
               IF (ABS(denom) .LE. 1.0E-18_dp*MAX(ABS(fb - fc), ABS(fb - fa))) THEN
                  npos = xb
               ELSE
                  npos = xb - 0.5_dp*nom/denom ! position of the stationary point
               END IF
               val = (npos - xa)*(npos - xb)*fc/((xc - xa)*(xc - xb)) + &
                     (npos - xb)*(npos - xc)*fa/((xa - xb)*(xa - xc)) + &
                     (npos - xc)*(npos - xa)*fb/((xb - xc)*(xb - xa))
               IF (val .LT. fa .AND. val .LE. fb .AND. val .LE. fc) THEN ! OK, we go to a minimum
                  ! we take a guard against too large steps
                  pos(4) = MAX(MAXVAL(pos(1:3))*0.01_dp, &
                               MIN(npos, MAXVAL(pos(1:3))*4.0_dp))
               ELSE ! just take an extended step
                  pos(4) = MAXVAL(pos(1:3))*2.0_dp
               END IF
            END SELECT
         CASE (3) ! golden section hunt
            new_direction = .FALSE.
            IF (line_search_count .EQ. 1) THEN
               lsl = 1
               lsr = 0
               lsm = 1
               pos(1) = 0.0_dp
               pos(2) = ds_min/gold_sec
            ELSE
               IF (line_search_count .EQ. 150) CPABORT("Too many")
               IF (lsr .EQ. 0) THEN
                  IF (energy(line_search_count - 1) .LT. energy(line_search_count)) THEN
                     lsr = line_search_count
                     pos(line_search_count + 1) = pos(lsm) + (pos(lsr) - pos(lsm))*gold_sec
                  ELSE
                     lsl = lsm
                     lsm = line_search_count
                     pos(line_search_count + 1) = pos(line_search_count)/gold_sec
                  END IF
               ELSE
                  IF (pos(line_search_count) .LT. pos(lsm)) THEN
                     IF (energy(line_search_count) .LT. energy(lsm)) THEN
                        lsr = lsm
                        lsm = line_search_count
                     ELSE
                        lsl = line_search_count
                     END IF
                  ELSE
                     IF (energy(line_search_count) .LT. energy(lsm)) THEN
                        lsl = lsm
                        lsm = line_search_count
                     ELSE
                        lsr = line_search_count
                     END IF
                  END IF
                  IF (pos(lsr) - pos(lsm) .GT. pos(lsm) - pos(lsl)) THEN
                     pos(line_search_count + 1) = pos(lsm) + gold_sec*(pos(lsr) - pos(lsm))
                  ELSE
                     pos(line_search_count + 1) = pos(lsl) + gold_sec*(pos(lsm) - pos(lsl))
                  END IF
                  IF ((pos(lsr) - pos(lsl)) .LT. 1.0E-3_dp*pos(lsr)) THEN
                     new_direction = .TRUE.
                  END IF
               END IF ! lsr .eq. 0
            END IF ! first step
         END SELECT
         ! now go to the suggested point
         ds_min = pos(line_search_count + 1)
         ds = pos(line_search_count + 1) - pos(line_search_count)
         CALL cp_fm_scale_and_add(1.0_dp, matrix_A, ds, matrix_G_search)
      END DO

      ! collect rotation matrix
      matrix_H%local_data = REAL(cmat_U%local_data, KIND=dp)
      CALL parallel_gemm('N', 'N', n, n, n, 1.0_dp, matrix_R, matrix_H, 0.0_dp, matrix_T)
      CALL cp_fm_to_fm(matrix_T, matrix_R)
      CALL rotate_orbitals(matrix_R, vectors)
      CALL cp_fm_release(matrix_R)

      CALL cp_fm_release(matrix_A)
      CALL cp_fm_release(matrix_G)
      CALL cp_fm_release(matrix_H)
      CALL cp_fm_release(matrix_T)
      CALL cp_fm_release(matrix_G_search)
      CALL cp_fm_release(matrix_G_old)
      CALL cp_cfm_release(cmat_A)
      CALL cp_cfm_release(cmat_U)
      CALL cp_cfm_release(cmat_R)
      CALL cp_cfm_release(cmat_t1)
      CALL cp_cfm_release(cmat_t2)
      CALL cp_cfm_release(cmat_B)
      CALL cp_cfm_release(cmat_M)

      DEALLOCATE (evals, evals_exp, fval, fvald)

      DO idim = 1, SIZE(c_zij)
         zij(1, idim)%local_data = REAL(c_zij(idim)%local_data, dp)
         zij(2, idim)%local_data = AIMAG(c_zij(idim)%local_data)
         CALL cp_cfm_release(c_zij(idim))
      END DO
      DEALLOCATE (c_zij)
      DEALLOCATE (diag_z)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief Parallel algorithm for jacobi rotations
!> \param weights ...
!> \param zij ...
!> \param vectors ...
!> \param para_env ...
!> \param max_iter ...
!> \param eps_localization ...
!> \param sweeps ...
!> \param out_each ...
!> \param target_time ...
!> \param start_time ...
!> \par History
!>      use allgather for improved performance
!> \author MI (11.2009)
! **************************************************************************************************
   SUBROUTINE jacobi_rot_para(weights, zij, vectors, para_env, max_iter, eps_localization, &
                              sweeps, out_each, target_time, start_time)

      REAL(KIND=dp), INTENT(IN)                          :: weights(:)
      TYPE(cp_fm_type), INTENT(IN)                       :: zij(:, :), vectors
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: max_iter
      REAL(KIND=dp), INTENT(IN)                          :: eps_localization
      INTEGER                                            :: sweeps
      INTEGER, INTENT(IN)                                :: out_each
      REAL(dp)                                           :: target_time, start_time

      CHARACTER(len=*), PARAMETER                        :: routineN = 'jacobi_rot_para'

      INTEGER                                            :: dim2, handle, i, idim, ii, ilow1, ip, j, &
                                                            nblock, nblock_max, ns_me, nstate, &
                                                            output_unit
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: ns_bound
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: rotmat, z_ij_loc_im, z_ij_loc_re
      REAL(KIND=dp)                                      :: xstate
      TYPE(cp_fm_type)                                   :: rmat
      TYPE(set_c_2d_type), DIMENSION(:), POINTER         :: cz_ij_loc

      CALL timeset(routineN, handle)

      output_unit = cp_logger_get_default_io_unit()

      NULLIFY (cz_ij_loc)

      dim2 = SIZE(zij, 2)

      CALL cp_fm_create(rmat, zij(1, 1)%matrix_struct)
      CALL cp_fm_set_all(rmat, 0._dp, 1._dp)

      CALL cp_fm_get_info(rmat, nrow_global=nstate)

      ! Distribution of the states (XXXXX safe against more pe than states ??? XXXXX)
      xstate = REAL(nstate, dp)/REAL(para_env%num_pe, dp)
      ALLOCATE (ns_bound(0:para_env%num_pe - 1, 2))
      DO ip = 1, para_env%num_pe
         ns_bound(ip - 1, 1) = MIN(nstate, NINT(xstate*(ip - 1))) + 1
         ns_bound(ip - 1, 2) = MIN(nstate, NINT(xstate*ip))
      END DO
      nblock_max = 0
      DO ip = 0, para_env%num_pe - 1
         nblock = ns_bound(ip, 2) - ns_bound(ip, 1) + 1
         nblock_max = MAX(nblock_max, nblock)
      END DO

      ! otbtain local part of the matrix (could be made faster, but is likely irrelevant).
      ALLOCATE (z_ij_loc_re(nstate, nblock_max))
      ALLOCATE (z_ij_loc_im(nstate, nblock_max))
      ALLOCATE (cz_ij_loc(dim2))
      DO idim = 1, dim2
         DO ip = 0, para_env%num_pe - 1
            nblock = ns_bound(ip, 2) - ns_bound(ip, 1) + 1
            CALL cp_fm_get_submatrix(zij(1, idim), z_ij_loc_re, 1, ns_bound(ip, 1), nstate, nblock)
            CALL cp_fm_get_submatrix(zij(2, idim), z_ij_loc_im, 1, ns_bound(ip, 1), nstate, nblock)
            IF (para_env%mepos == ip) THEN
               ALLOCATE (cz_ij_loc(idim)%c_array(nstate, nblock))
               DO i = 1, nblock
                  DO j = 1, nstate
                     cz_ij_loc(idim)%c_array(j, i) = CMPLX(z_ij_loc_re(j, i), z_ij_loc_im(j, i), dp)
                  END DO
               END DO
            END IF
         END DO ! ip
      END DO
      DEALLOCATE (z_ij_loc_re)
      DEALLOCATE (z_ij_loc_im)

      ALLOCATE (rotmat(nstate, 2*nblock_max))

      CALL jacobi_rot_para_core(weights, para_env, max_iter, sweeps, out_each, dim2, nstate, nblock_max, ns_bound, &
                                cz_ij_loc, rotmat, output_unit, eps_localization=eps_localization, &
                                target_time=target_time, start_time=start_time)

      ilow1 = ns_bound(para_env%mepos, 1)
      ns_me = ns_bound(para_env%mepos, 2) - ns_bound(para_env%mepos, 1) + 1
      ALLOCATE (z_ij_loc_re(nstate, nblock_max))
      ALLOCATE (z_ij_loc_im(nstate, nblock_max))
      DO idim = 1, dim2
         DO ip = 0, para_env%num_pe - 1
            z_ij_loc_re = 0.0_dp
            z_ij_loc_im = 0.0_dp
            nblock = ns_bound(ip, 2) - ns_bound(ip, 1) + 1
            IF (ip == para_env%mepos) THEN
               ns_me = nblock
               DO i = 1, ns_me
                  ii = ilow1 + i - 1
                  DO j = 1, nstate
                     z_ij_loc_re(j, i) = REAL(cz_ij_loc(idim)%c_array(j, i), dp)
                     z_ij_loc_im(j, i) = AIMAG(cz_ij_loc(idim)%c_array(j, i))
                  END DO
               END DO
            END IF
            CALL para_env%bcast(z_ij_loc_re, ip)
            CALL para_env%bcast(z_ij_loc_im, ip)
            CALL cp_fm_set_submatrix(zij(1, idim), z_ij_loc_re, 1, ns_bound(ip, 1), nstate, nblock)
            CALL cp_fm_set_submatrix(zij(2, idim), z_ij_loc_im, 1, ns_bound(ip, 1), nstate, nblock)
         END DO ! ip
      END DO

      DO ip = 0, para_env%num_pe - 1
         z_ij_loc_re = 0.0_dp
         nblock = ns_bound(ip, 2) - ns_bound(ip, 1) + 1
         IF (ip == para_env%mepos) THEN
            ns_me = nblock
            DO i = 1, ns_me
               ii = ilow1 + i - 1
               DO j = 1, nstate
                  z_ij_loc_re(j, i) = rotmat(j, i)
               END DO
            END DO
         END IF
         CALL para_env%bcast(z_ij_loc_re, ip)
         CALL cp_fm_set_submatrix(rmat, z_ij_loc_re, 1, ns_bound(ip, 1), nstate, nblock)
      END DO

      DEALLOCATE (z_ij_loc_re)
      DEALLOCATE (z_ij_loc_im)
      DO idim = 1, dim2
         DEALLOCATE (cz_ij_loc(idim)%c_array)
      END DO
      DEALLOCATE (cz_ij_loc)

      CALL para_env%sync()
      CALL rotate_orbitals(rmat, vectors)
      CALL cp_fm_release(rmat)

      DEALLOCATE (rotmat)
      DEALLOCATE (ns_bound)

      CALL timestop(handle)

   END SUBROUTINE jacobi_rot_para

! **************************************************************************************************
!> \brief almost identical to 'jacobi_rot_para' but with different inout variables
!> \param weights ...
!> \param czij ...
!> \param para_env ...
!> \param max_iter ...
!> \param rmat ...
!> \param tol_out ...
!> \author Soumya Ghosh (08/21)
! **************************************************************************************************
   SUBROUTINE jacobi_rot_para_1(weights, czij, para_env, max_iter, rmat, tol_out)

      REAL(KIND=dp), INTENT(IN)                          :: weights(:)
      TYPE(cp_cfm_type), INTENT(IN)                      :: czij(:)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: max_iter
      TYPE(cp_cfm_type), INTENT(IN)                      :: rmat
      REAL(dp), INTENT(OUT), OPTIONAL                    :: tol_out

      CHARACTER(len=*), PARAMETER                        :: routineN = 'jacobi_rot_para_1'

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: czij_array
      INTEGER                                            :: dim2, handle, i, idim, ii, ilow1, ip, j, &
                                                            nblock, nblock_max, ns_me, nstate, &
                                                            sweeps
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: ns_bound
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: rotmat, z_ij_loc_re
      REAL(KIND=dp)                                      :: xstate
      TYPE(set_c_2d_type), DIMENSION(:), POINTER         :: cz_ij_loc

      CALL timeset(routineN, handle)

      dim2 = SIZE(czij)

      CALL cp_cfm_set_all(rmat, CMPLX(0._dp, 0._dp, dp), CMPLX(1._dp, 0._dp, dp))

      CALL cp_cfm_get_info(rmat, nrow_global=nstate)

      ! Distribution of the states (XXXXX safe against more pe than states ??? XXXXX)
      xstate = REAL(nstate, dp)/REAL(para_env%num_pe, dp)
      ALLOCATE (ns_bound(0:para_env%num_pe - 1, 2))
      DO ip = 1, para_env%num_pe
         ns_bound(ip - 1, 1) = MIN(nstate, NINT(xstate*(ip - 1))) + 1
         ns_bound(ip - 1, 2) = MIN(nstate, NINT(xstate*ip))
      END DO
      nblock_max = 0
      DO ip = 0, para_env%num_pe - 1
         nblock = ns_bound(ip, 2) - ns_bound(ip, 1) + 1
         nblock_max = MAX(nblock_max, nblock)
      END DO

      ! otbtain local part of the matrix (could be made faster, but is likely irrelevant).
      ALLOCATE (czij_array(nstate, nblock_max))
      ALLOCATE (cz_ij_loc(dim2))
      DO idim = 1, dim2
         DO ip = 0, para_env%num_pe - 1
            nblock = ns_bound(ip, 2) - ns_bound(ip, 1) + 1
            ! cfm --> allocatable
            CALL cp_cfm_get_submatrix(czij(idim), czij_array, 1, ns_bound(ip, 1), nstate, nblock)
            IF (para_env%mepos == ip) THEN
               ns_me = nblock
               ALLOCATE (cz_ij_loc(idim)%c_array(nstate, ns_me))
               DO i = 1, ns_me
                  DO j = 1, nstate
                     cz_ij_loc(idim)%c_array(j, i) = czij_array(j, i)
                  END DO
               END DO
            END IF
         END DO ! ip
      END DO
      DEALLOCATE (czij_array)

      ALLOCATE (rotmat(nstate, 2*nblock_max))

      CALL jacobi_rot_para_core(weights, para_env, max_iter, sweeps, 1, dim2, nstate, nblock_max, ns_bound, &
                                cz_ij_loc, rotmat, 0, tol_out=tol_out)

      ilow1 = ns_bound(para_env%mepos, 1)
      ns_me = ns_bound(para_env%mepos, 2) - ns_bound(para_env%mepos, 1) + 1
      ALLOCATE (z_ij_loc_re(nstate, nblock_max))

      DO ip = 0, para_env%num_pe - 1
         z_ij_loc_re = 0.0_dp
         nblock = ns_bound(ip, 2) - ns_bound(ip, 1) + 1
         IF (ip == para_env%mepos) THEN
            ns_me = nblock
            DO i = 1, ns_me
               ii = ilow1 + i - 1
               DO j = 1, nstate
                  z_ij_loc_re(j, i) = rotmat(j, i)
               END DO
            END DO
         END IF
         CALL para_env%bcast(z_ij_loc_re, ip)
         CALL cp_cfm_set_submatrix(rmat, CMPLX(z_ij_loc_re, 0.0_dp, dp), 1, ns_bound(ip, 1), nstate, nblock)
      END DO

      DEALLOCATE (z_ij_loc_re)
      DO idim = 1, dim2
         DEALLOCATE (cz_ij_loc(idim)%c_array)
      END DO
      DEALLOCATE (cz_ij_loc)

      CALL para_env%sync()

      DEALLOCATE (rotmat)
      DEALLOCATE (ns_bound)

      CALL timestop(handle)

   END SUBROUTINE jacobi_rot_para_1

! **************************************************************************************************
!> \brief Parallel algorithm for jacobi rotations
!> \param weights ...
!> \param para_env ...
!> \param max_iter ...
!> \param sweeps ...
!> \param out_each ...
!> \param dim2 ...
!> \param nstate ...
!> \param nblock_max ...
!> \param ns_bound ...
!> \param cz_ij_loc ...
!> \param rotmat ...
!> \param output_unit ...
!> \param tol_out ...
!> \param eps_localization ...
!> \param target_time ...
!> \param start_time ...
!> \par History
!>      split out to reuse with different input types
!> \author HF (05.2022)
! **************************************************************************************************
   SUBROUTINE jacobi_rot_para_core(weights, para_env, max_iter, sweeps, out_each, dim2, nstate, nblock_max, &
                                   ns_bound, cz_ij_loc, rotmat, output_unit, tol_out, eps_localization, target_time, start_time)

      REAL(KIND=dp), INTENT(IN)                          :: weights(:)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: max_iter
      INTEGER, INTENT(OUT)                               :: sweeps
      INTEGER, INTENT(IN)                                :: out_each, dim2, nstate, nblock_max
      INTEGER, DIMENSION(0:, :), INTENT(IN)              :: ns_bound
      TYPE(set_c_2d_type), DIMENSION(:), POINTER         :: cz_ij_loc
      REAL(dp), DIMENSION(:, :), INTENT(OUT)             :: rotmat
      INTEGER, INTENT(IN)                                :: output_unit
      REAL(dp), INTENT(OUT), OPTIONAL                    :: tol_out
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: eps_localization
      REAL(dp), OPTIONAL                                 :: target_time, start_time

      COMPLEX(KIND=dp)                                   :: zi, zj
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: c_array_me, c_array_partner
      COMPLEX(KIND=dp), POINTER                          :: mii(:), mij(:), mjj(:)
      INTEGER :: i, idim, ii, ik, il1, il2, il_recv, il_recv_partner, ilow1, ilow2, ip, ip_has_i, &
         ip_partner, ip_recv_from, ip_recv_partner, ipair, iperm, istat, istate, iu1, iu2, iup1, &
         iup2, j, jj, jstate, k, kk, lsweep, n1, n2, npair, nperm, ns_me, ns_partner, &
         ns_recv_from, ns_recv_partner
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: rcount, rdispl
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: list_pair
      LOGICAL                                            :: should_stop
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: gmat, rmat_loc, rmat_recv, rmat_send
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: rmat_recv_all
      REAL(KIND=dp)                                      :: ct, func, gmax, grad, ri, rj, st, t1, &
                                                            t2, theta, tolerance, zc, zr
      TYPE(set_c_1d_type), DIMENSION(:), POINTER         :: zdiag_all, zdiag_me
      TYPE(set_c_2d_type), DIMENSION(:), POINTER         :: xyz_mix, xyz_mix_ns

      NULLIFY (zdiag_all, zdiag_me)
      NULLIFY (xyz_mix, xyz_mix_ns)
      NULLIFY (mii, mij, mjj)

      ALLOCATE (mii(dim2), mij(dim2), mjj(dim2))

      ALLOCATE (rcount(para_env%num_pe), STAT=istat)
      ALLOCATE (rdispl(para_env%num_pe), STAT=istat)

      tolerance = 1.0e10_dp
      sweeps = 0

      ! number of processor pairs and number of permutations
      npair = (para_env%num_pe + 1)/2
      nperm = para_env%num_pe - MOD(para_env%num_pe + 1, 2)
      ALLOCATE (list_pair(2, npair))

      ! initialize rotation matrix
      rotmat = 0.0_dp
      DO i = ns_bound(para_env%mepos, 1), ns_bound(para_env%mepos, 2)
         ii = i - ns_bound(para_env%mepos, 1) + 1
         rotmat(i, ii) = 1.0_dp
      END DO

      ALLOCATE (xyz_mix(dim2))
      ALLOCATE (xyz_mix_ns(dim2))
      ALLOCATE (zdiag_me(dim2))
      ALLOCATE (zdiag_all(dim2))

      ns_me = ns_bound(para_env%mepos, 2) - ns_bound(para_env%mepos, 1) + 1
      IF (ns_me /= 0) THEN
         ALLOCATE (c_array_me(nstate, ns_me, dim2))
         DO idim = 1, dim2
            ALLOCATE (xyz_mix_ns(idim)%c_array(nstate, ns_me))
         END DO
         ALLOCATE (gmat(nstate, ns_me))
      END IF

      DO idim = 1, dim2
         ALLOCATE (zdiag_me(idim)%c_array(nblock_max))
         zdiag_me(idim)%c_array = CMPLX(0.0_dp, 0.0_dp, dp)
         ALLOCATE (zdiag_all(idim)%c_array(para_env%num_pe*nblock_max))
         zdiag_all(idim)%c_array = CMPLX(0.0_dp, 0.0_dp, dp)
      END DO
      ALLOCATE (rmat_recv(nblock_max*2, nblock_max))
      ALLOCATE (rmat_send(nblock_max*2, nblock_max))

      ! buffer for message passing
      ALLOCATE (rmat_recv_all(nblock_max*2, nblock_max, 0:para_env%num_pe - 1))

      IF (output_unit > 0) THEN
         WRITE (output_unit, '(T4,A )') " Localization by iterative distributed Jacobi rotation"
         WRITE (output_unit, '(T20,A12,T32, A22,T60, A12,A8 )') "Iteration", "Functional", "Tolerance", " Time "
      END IF

      DO lsweep = 1, max_iter + 1
         sweeps = lsweep
         IF (sweeps == max_iter + 1) THEN
            IF (output_unit > 0) THEN
               WRITE (output_unit, *) ' LOCALIZATION! loop did not converge within the maximum number of iterations.'
               WRITE (output_unit, *) '               Present  Max. gradient = ', tolerance
            END IF
            EXIT
         END IF
         t1 = m_walltime()

         DO iperm = 1, nperm

            ! fix partners for this permutation, and get the number of states
            CALL eberlein(iperm, para_env, list_pair)
            ip_partner = -1
            ns_partner = 0
            DO ipair = 1, npair
               IF (list_pair(1, ipair) == para_env%mepos) THEN
                  ip_partner = list_pair(2, ipair)
                  EXIT
               ELSE IF (list_pair(2, ipair) == para_env%mepos) THEN
                  ip_partner = list_pair(1, ipair)
                  EXIT
               END IF
            END DO
            IF (ip_partner >= 0) THEN
               ns_partner = ns_bound(ip_partner, 2) - ns_bound(ip_partner, 1) + 1
            ELSE
               ns_partner = 0
            END IF

            ! if there is a non-zero block connecting two partners, jacobi-sweep it.
            IF (ns_partner*ns_me /= 0) THEN

               ALLOCATE (rmat_loc(ns_me + ns_partner, ns_me + ns_partner))
               rmat_loc = 0.0_dp
               DO i = 1, ns_me + ns_partner
                  rmat_loc(i, i) = 1.0_dp
               END DO

               ALLOCATE (c_array_partner(nstate, ns_partner, dim2))

               DO idim = 1, dim2
                  ALLOCATE (xyz_mix(idim)%c_array(ns_me + ns_partner, ns_me + ns_partner))
                  DO i = 1, ns_me
                     c_array_me(1:nstate, i, idim) = cz_ij_loc(idim)%c_array(1:nstate, i)
                  END DO
               END DO

               CALL para_env%sendrecv(msgin=c_array_me, dest=ip_partner, &
                                      msgout=c_array_partner, source=ip_partner)

               n1 = ns_me
               n2 = ns_partner
               ilow1 = ns_bound(para_env%mepos, 1)
               iup1 = ns_bound(para_env%mepos, 1) + n1 - 1
               ilow2 = ns_bound(ip_partner, 1)
               iup2 = ns_bound(ip_partner, 1) + n2 - 1
               IF (ns_bound(para_env%mepos, 1) < ns_bound(ip_partner, 1)) THEN
                  il1 = 1
                  iu1 = n1
                  iu1 = n1
                  il2 = 1 + n1
                  iu2 = n1 + n2
               ELSE
                  il1 = 1 + n2
                  iu1 = n1 + n2
                  iu1 = n1 + n2
                  il2 = 1
                  iu2 = n2
               END IF

               DO idim = 1, dim2
                  DO i = 1, n1
                     xyz_mix(idim)%c_array(il1:iu1, il1 + i - 1) = c_array_me(ilow1:iup1, i, idim)
                     xyz_mix(idim)%c_array(il2:iu2, il1 + i - 1) = c_array_me(ilow2:iup2, i, idim)
                  END DO
                  DO i = 1, n2
                     xyz_mix(idim)%c_array(il2:iu2, il2 + i - 1) = c_array_partner(ilow2:iup2, i, idim)
                     xyz_mix(idim)%c_array(il1:iu1, il2 + i - 1) = c_array_partner(ilow1:iup1, i, idim)
                  END DO
               END DO

               DO istate = 1, n1 + n2
                  DO jstate = istate + 1, n1 + n2
                     DO idim = 1, dim2
                        mii(idim) = xyz_mix(idim)%c_array(istate, istate)
                        mij(idim) = xyz_mix(idim)%c_array(istate, jstate)
                        mjj(idim) = xyz_mix(idim)%c_array(jstate, jstate)
                     END DO
                     CALL get_angle(mii, mjj, mij, weights, theta)
                     st = SIN(theta)
                     ct = COS(theta)
                     DO idim = 1, dim2
                        DO i = 1, n1 + n2
                           zi = ct*xyz_mix(idim)%c_array(i, istate) + st*xyz_mix(idim)%c_array(i, jstate)
                           zj = -st*xyz_mix(idim)%c_array(i, istate) + ct*xyz_mix(idim)%c_array(i, jstate)
                           xyz_mix(idim)%c_array(i, istate) = zi
                           xyz_mix(idim)%c_array(i, jstate) = zj
                        END DO
                        DO i = 1, n1 + n2
                           zi = ct*xyz_mix(idim)%c_array(istate, i) + st*xyz_mix(idim)%c_array(jstate, i)
                           zj = -st*xyz_mix(idim)%c_array(istate, i) + ct*xyz_mix(idim)%c_array(jstate, i)
                           xyz_mix(idim)%c_array(istate, i) = zi
                           xyz_mix(idim)%c_array(jstate, i) = zj
                        END DO
                     END DO

                     DO i = 1, n1 + n2
                        ri = ct*rmat_loc(i, istate) + st*rmat_loc(i, jstate)
                        rj = ct*rmat_loc(i, jstate) - st*rmat_loc(i, istate)
                        rmat_loc(i, istate) = ri
                        rmat_loc(i, jstate) = rj
                     END DO
                  END DO
               END DO

               k = nblock_max + 1
               CALL para_env%sendrecv(rotmat(1:nstate, 1:ns_me), ip_partner, &
                                      rotmat(1:nstate, k:k + n2 - 1), ip_partner)

               IF (ilow1 < ilow2) THEN
                  ! no longer compiles in official sdgb:
                  !CALL dgemm("N", "N", nstate, n1, n2, 1.0_dp, rotmat(1, k), nstate, rmat_loc(1 + n1, 1), n1 + n2, 0.0_dp, gmat, nstate)
                  ! probably inefficient:
                  CALL dgemm("N", "N", nstate, n1, n2, 1.0_dp, rotmat(1:, k:), nstate, rmat_loc(1 + n1:, 1:n1), &
                             n2, 0.0_dp, gmat(:, :), nstate)
                  CALL dgemm("N", "N", nstate, n1, n1, 1.0_dp, rotmat(1:, 1:), nstate, rmat_loc(1:, 1:), &
                             n1 + n2, 1.0_dp, gmat(:, :), nstate)
               ELSE
                  CALL dgemm("N", "N", nstate, n1, n2, 1.0_dp, rotmat(1:, k:), nstate, &
                             rmat_loc(1:, n2 + 1:), n1 + n2, 0.0_dp, gmat(:, :), nstate)
                  ! no longer compiles in official sdgb:
                  !CALL dgemm("N", "N", nstate, n1, n1, 1.0_dp, rotmat(1, 1), nstate, rmat_loc(n2 + 1, n2 + 1), n1 + n2, 1.0_dp, gmat, nstate)
                  ! probably inefficient:
                  CALL dgemm("N", "N", nstate, n1, n1, 1.0_dp, rotmat(1:, 1:), nstate, rmat_loc(n2 + 1:, n2 + 1:), &
                             n1, 1.0_dp, gmat(:, :), nstate)
               END IF

               CALL dcopy(nstate*n1, gmat(1, 1), 1, rotmat(1, 1), 1)

               DO idim = 1, dim2
                  DO i = 1, n1
                     xyz_mix_ns(idim)%c_array(1:nstate, i) = CMPLX(0.0_dp, 0.0_dp, dp)
                  END DO

                  DO istate = 1, n1
                     DO jstate = 1, nstate
                        DO i = 1, n2
                           xyz_mix_ns(idim)%c_array(jstate, istate) = &
                              xyz_mix_ns(idim)%c_array(jstate, istate) + &
                              c_array_partner(jstate, i, idim)*rmat_loc(il2 + i - 1, il1 + istate - 1)
                        END DO
                     END DO
                  END DO
                  DO istate = 1, n1
                     DO jstate = 1, nstate
                        DO i = 1, n1
                           xyz_mix_ns(idim)%c_array(jstate, istate) = xyz_mix_ns(idim)%c_array(jstate, istate) + &
                                                                 c_array_me(jstate, i, idim)*rmat_loc(il1 + i - 1, il1 + istate - 1)
                        END DO
                     END DO
                  END DO
               END DO ! idim

               DEALLOCATE (c_array_partner)

            ELSE ! save my data
               DO idim = 1, dim2
                  DO i = 1, ns_me
                     xyz_mix_ns(idim)%c_array(1:nstate, i) = cz_ij_loc(idim)%c_array(1:nstate, i)
                  END DO
               END DO
            END IF

            DO idim = 1, dim2
               DO i = 1, ns_me
                  cz_ij_loc(idim)%c_array(1:nstate, i) = CMPLX(0.0_dp, 0.0_dp, dp)
               END DO
            END DO

            IF (ns_partner*ns_me /= 0) THEN
               ! transpose rotation matrix rmat_loc
               DO i = 1, ns_me + ns_partner
                  DO j = i + 1, ns_me + ns_partner
                     ri = rmat_loc(i, j)
                     rmat_loc(i, j) = rmat_loc(j, i)
                     rmat_loc(j, i) = ri
                  END DO
               END DO

               ! prepare for distribution
               DO i = 1, n1
                  rmat_send(1:n1, i) = rmat_loc(il1:iu1, il1 + i - 1)
               END DO
               ik = nblock_max
               DO i = 1, n2
                  rmat_send(ik + 1:ik + n1, i) = rmat_loc(il1:iu1, il2 + i - 1)
               END DO
            ELSE
               rmat_send = 0.0_dp
            END IF

            ! collect data from all tasks (this takes some significant time)
            CALL para_env%allgather(rmat_send, rmat_recv_all)

            ! update blocks everywhere
            DO ip = 0, para_env%num_pe - 1

               ip_recv_from = MOD(para_env%mepos - IP + para_env%num_pe, para_env%num_pe)
               rmat_recv(:, :) = rmat_recv_all(:, :, ip_recv_from)

               ns_recv_from = ns_bound(ip_recv_from, 2) - ns_bound(ip_recv_from, 1) + 1

               IF (ns_me /= 0) THEN
                  IF (ns_recv_from /= 0) THEN
                     !look for the partner of ip_recv_from
                     ip_recv_partner = -1
                     ns_recv_partner = 0
                     DO ipair = 1, npair
                        IF (list_pair(1, ipair) == ip_recv_from) THEN
                           ip_recv_partner = list_pair(2, ipair)
                           EXIT
                        ELSE IF (list_pair(2, ipair) == ip_recv_from) THEN
                           ip_recv_partner = list_pair(1, ipair)
                           EXIT
                        END IF
                     END DO

                     IF (ip_recv_partner >= 0) THEN
                        ns_recv_partner = ns_bound(ip_recv_partner, 2) - ns_bound(ip_recv_partner, 1) + 1
                     END IF
                     IF (ns_recv_partner > 0) THEN
                        il1 = ns_bound(para_env%mepos, 1)
                        il_recv = ns_bound(ip_recv_from, 1)
                        il_recv_partner = ns_bound(ip_recv_partner, 1)
                        ik = nblock_max

                        DO idim = 1, dim2
                           DO i = 1, ns_recv_from
                              ii = il_recv + i - 1
                              DO j = 1, ns_me
                                 jj = j
                                 DO k = 1, ns_recv_from
                                    kk = il_recv + k - 1
                                    cz_ij_loc(idim)%c_array(ii, jj) = cz_ij_loc(idim)%c_array(ii, jj) + &
                                                                      rmat_recv(i, k)*xyz_mix_ns(idim)%c_array(kk, j)
                                 END DO
                              END DO
                           END DO
                           DO i = 1, ns_recv_from
                              ii = il_recv + i - 1
                              DO j = 1, ns_me
                                 jj = j
                                 DO k = 1, ns_recv_partner
                                    kk = il_recv_partner + k - 1
                                    cz_ij_loc(idim)%c_array(ii, jj) = cz_ij_loc(idim)%c_array(ii, jj) + &
                                                                      rmat_recv(ik + i, k)*xyz_mix_ns(idim)%c_array(kk, j)
                                 END DO
                              END DO
                           END DO
                        END DO ! idim
                     ELSE
                        il1 = ns_bound(para_env%mepos, 1)
                        il_recv = ns_bound(ip_recv_from, 1)
                        DO idim = 1, dim2
                           DO j = 1, ns_me
                              jj = j
                              DO i = 1, ns_recv_from
                                 ii = il_recv + i - 1
                                 cz_ij_loc(idim)%c_array(ii, jj) = xyz_mix_ns(idim)%c_array(ii, j)
                              END DO
                           END DO
                        END DO ! idim
                     END IF
                  END IF
               END IF ! ns_me
            END DO ! ip

            IF (ns_partner*ns_me /= 0) THEN
               DEALLOCATE (rmat_loc)
               DO idim = 1, dim2
                  DEALLOCATE (xyz_mix(idim)%c_array)
               END DO
            END IF

         END DO ! iperm

         ! calculate the max gradient
         DO idim = 1, dim2
            DO i = ns_bound(para_env%mepos, 1), ns_bound(para_env%mepos, 2)
               ii = i - ns_bound(para_env%mepos, 1) + 1
               zdiag_me(idim)%c_array(ii) = cz_ij_loc(idim)%c_array(i, ii)
               zdiag_me(idim)%c_array(ii) = cz_ij_loc(idim)%c_array(i, ii)
            END DO
            rcount(:) = SIZE(zdiag_me(idim)%c_array)
            rdispl(1) = 0
            DO ip = 2, para_env%num_pe
               rdispl(ip) = rdispl(ip - 1) + rcount(ip - 1)
            END DO
            ! collect all the diagonal elements in a replicated 1d array
            CALL para_env%allgatherv(zdiag_me(idim)%c_array, zdiag_all(idim)%c_array, rcount, rdispl)
         END DO

         gmax = 0.0_dp
         DO j = ns_bound(para_env%mepos, 1), ns_bound(para_env%mepos, 2)
            k = j - ns_bound(para_env%mepos, 1) + 1
            DO i = 1, j - 1
               ! find the location of the diagonal element (i,i)
               DO ip = 0, para_env%num_pe - 1
                  IF (i >= ns_bound(ip, 1) .AND. i <= ns_bound(ip, 2)) THEN
                     ip_has_i = ip
                     EXIT
                  END IF
               END DO
               ii = nblock_max*ip_has_i + i - ns_bound(ip_has_i, 1) + 1
               ! mepos has the diagonal element (j,j), as well as the off diagonal (i,j)
               jj = nblock_max*para_env%mepos + j - ns_bound(para_env%mepos, 1) + 1
               grad = 0.0_dp
               DO idim = 1, dim2
                  zi = zdiag_all(idim)%c_array(ii)
                  zj = zdiag_all(idim)%c_array(jj)
                  grad = grad + weights(idim)*REAL(4.0_dp*CONJG(cz_ij_loc(idim)%c_array(i, k))*(zj - zi), dp)
               END DO
               gmax = MAX(gmax, ABS(grad))
            END DO
         END DO

         CALL para_env%max(gmax)
         tolerance = gmax
         IF (PRESENT(tol_out)) tol_out = tolerance

         func = 0.0_dp
         DO i = ns_bound(para_env%mepos, 1), ns_bound(para_env%mepos, 2)
            k = i - ns_bound(para_env%mepos, 1) + 1
            DO idim = 1, dim2
               zr = REAL(cz_ij_loc(idim)%c_array(i, k), dp)
               zc = AIMAG(cz_ij_loc(idim)%c_array(i, k))
               func = func + weights(idim)*(1.0_dp - (zr*zr + zc*zc))/twopi/twopi
            END DO
         END DO
         CALL para_env%sum(func)
         t2 = m_walltime()

         IF (output_unit > 0 .AND. MODULO(sweeps, out_each) == 0) THEN
            WRITE (output_unit, '(T20,I12,T35,F20.10,T60,E12.4,F8.3)') sweeps, func, tolerance, t2 - t1
            CALL m_flush(output_unit)
         END IF
         IF (PRESENT(eps_localization)) THEN
            IF (tolerance < eps_localization) EXIT
         END IF
         IF (PRESENT(target_time) .AND. PRESENT(start_time)) THEN
            CALL external_control(should_stop, "LOC", target_time=target_time, start_time=start_time)
            IF (should_stop) EXIT
         END IF

      END DO ! lsweep

      ! buffer for message passing
      DEALLOCATE (rmat_recv_all)

      DEALLOCATE (rmat_recv)
      DEALLOCATE (rmat_send)
      IF (ns_me > 0) THEN
         DEALLOCATE (c_array_me)
      END IF
      DO idim = 1, dim2
         DEALLOCATE (zdiag_me(idim)%c_array)
         DEALLOCATE (zdiag_all(idim)%c_array)
      END DO
      DEALLOCATE (zdiag_me)
      DEALLOCATE (zdiag_all)
      DEALLOCATE (xyz_mix)
      DO idim = 1, dim2
         IF (ns_me /= 0) THEN
            DEALLOCATE (xyz_mix_ns(idim)%c_array)
         END IF
      END DO
      DEALLOCATE (xyz_mix_ns)
      IF (ns_me /= 0) THEN
         DEALLOCATE (gmat)
      END IF
      DEALLOCATE (mii)
      DEALLOCATE (mij)
      DEALLOCATE (mjj)
      DEALLOCATE (list_pair)

   END SUBROUTINE jacobi_rot_para_core

! **************************************************************************************************
!> \brief ...
!> \param iperm ...
!> \param para_env ...
!> \param list_pair ...
! **************************************************************************************************
   SUBROUTINE eberlein(iperm, para_env, list_pair)
      INTEGER, INTENT(IN)                                :: iperm
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, DIMENSION(:, :)                           :: list_pair

      INTEGER                                            :: i, ii, jj, npair

      npair = (para_env%num_pe + 1)/2
      IF (iperm == 1) THEN
!..set up initial ordering
         DO I = 0, para_env%num_pe - 1
            II = ((i + 1) + 1)/2
            JJ = MOD((i + 1) + 1, 2) + 1
            list_pair(JJ, II) = i
         END DO
         IF (MOD(para_env%num_pe, 2) == 1) list_pair(2, npair) = -1
      ELSEIF (MOD(iperm, 2) == 0) THEN
!..a type shift
         jj = list_pair(1, npair)
         DO I = npair, 3, -1
            list_pair(1, I) = list_pair(1, I - 1)
         END DO
         list_pair(1, 2) = list_pair(2, 1)
         list_pair(2, 1) = jj
      ELSE
!..b type shift
         jj = list_pair(2, 1)
         DO I = 1, npair - 1
            list_pair(2, I) = list_pair(2, I + 1)
         END DO
         list_pair(2, npair) = jj
      END IF

   END SUBROUTINE eberlein

! **************************************************************************************************
!> \brief ...
!> \param vectors ...
!> \param op_sm_set ...
!> \param zij_fm_set ...
! **************************************************************************************************
   SUBROUTINE zij_matrix(vectors, op_sm_set, zij_fm_set)

      TYPE(cp_fm_type), INTENT(IN)                       :: vectors
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: op_sm_set
      TYPE(cp_fm_type), DIMENSION(:, :), INTENT(IN)      :: zij_fm_set

      CHARACTER(len=*), PARAMETER                        :: routineN = 'zij_matrix'

      INTEGER                                            :: handle, i, j, nao, nmoloc
      TYPE(cp_fm_type)                                   :: opvec

      CALL timeset(routineN, handle)

      ! get rows and cols of the input
      CALL cp_fm_get_info(vectors, nrow_global=nao, ncol_global=nmoloc)
      ! replicate the input kind of matrix
      CALL cp_fm_create(opvec, vectors%matrix_struct)

      ! Compute zij here
      DO i = 1, SIZE(zij_fm_set, 2)
         DO j = 1, SIZE(zij_fm_set, 1)
            CALL cp_fm_set_all(zij_fm_set(j, i), 0.0_dp)
            CALL cp_dbcsr_sm_fm_multiply(op_sm_set(j, i)%matrix, vectors, opvec, ncol=nmoloc)
            CALL parallel_gemm("T", "N", nmoloc, nmoloc, nao, 1.0_dp, vectors, opvec, 0.0_dp, &
                               zij_fm_set(j, i))
         END DO
      END DO

      CALL cp_fm_release(opvec)
      CALL timestop(handle)

   END SUBROUTINE zij_matrix

! **************************************************************************************************
!> \brief ...
!> \param vectors ...
! **************************************************************************************************
   SUBROUTINE scdm_qrfact(vectors)

      TYPE(cp_fm_type), INTENT(IN)                       :: vectors

      CHARACTER(len=*), PARAMETER                        :: routineN = 'scdm_qrfact'

      INTEGER                                            :: handle, ncolT, nrowT
      REAL(KIND=dp), DIMENSION(:), POINTER               :: tau
      TYPE(cp_fm_struct_type), POINTER                   :: cstruct
      TYPE(cp_fm_type)                                   :: CTp, Qf, tmp

      CALL timeset(routineN, handle)

      ! Create Transpose of Coefficient Matrix vectors
      nrowT = vectors%matrix_struct%ncol_global
      ncolT = vectors%matrix_struct%nrow_global

      CALL cp_fm_struct_create(cstruct, vectors%matrix_struct%para_env, &
                               vectors%matrix_struct%context, nrowT, ncolT, vectors%matrix_struct%ncol_block, &
                               vectors%matrix_struct%nrow_block, first_p_pos=vectors%matrix_struct%first_p_pos)
      CALL cp_fm_create(CTp, cstruct)
      CALL cp_fm_struct_release(cstruct)

      ALLOCATE (tau(nrowT))

      CALL cp_fm_transpose(vectors, CTp)

      ! Get QR decomposition of CTs
      CALL cp_fm_pdgeqpf(CTp, tau, nrowT, ncolT, 1, 1)

      ! Construction of Q from the scalapack output
      CALL cp_fm_struct_create(cstruct, para_env=CTp%matrix_struct%para_env, &
                               context=CTp%matrix_struct%context, nrow_global=CTp%matrix_struct%nrow_global, &
                               ncol_global=CTp%matrix_struct%nrow_global)
      CALL cp_fm_create(Qf, cstruct)
      CALL cp_fm_struct_release(cstruct)
      CALL cp_fm_to_fm_submat(CTp, Qf, nrowT, nrowT, 1, 1, 1, 1)

      ! Get Q
      CALL cp_fm_pdorgqr(Qf, tau, nrowT, 1, 1)

      ! Transform original coefficient matrix vectors
      CALL cp_fm_create(tmp, vectors%matrix_struct)
      CALL cp_fm_set_all(tmp, 0.0_dp, 1.0_dp)
      CALL cp_fm_to_fm(vectors, tmp)
      CALL parallel_gemm('N', 'N', ncolT, nrowT, nrowT, 1.0_dp, tmp, Qf, 0.0_dp, vectors)

      ! Cleanup
      CALL cp_fm_release(CTp)
      CALL cp_fm_release(tmp)
      CALL cp_fm_release(Qf)
      DEALLOCATE (tau)

      CALL timestop(handle)

   END SUBROUTINE scdm_qrfact

END MODULE qs_localization_methods
