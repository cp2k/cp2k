!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Constrained fitting of a real-space density by an AO density matrix.
!>        This is a dense reference implementation with cubic spectral projections.
!> \author CP2K developers group
! **************************************************************************************************
MODULE qs_density_fit
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_copy,&
                                              dbcsr_create,&
                                              dbcsr_get_info,&
                                              dbcsr_p_type,&
                                              dbcsr_release,&
                                              dbcsr_set,&
                                              dbcsr_type_symmetric
   USE cp_dbcsr_contrib,                ONLY: dbcsr_dot
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                              cp_fm_scale_and_add,&
                                              cp_fm_trace,&
                                              cp_fm_uplo_to_full
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                              cp_fm_cholesky_restore
   USE cp_fm_diag,                      ONLY: cp_fm_syevd
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_copy,&
                                              pw_integral_ab,&
                                              pw_integrate_function
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_collocate_density,            ONLY: calculate_rho_elec
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integrate_potential,          ONLY: integrate_v_rspace
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_density_fit'

   PUBLIC :: fit_constrained_density

CONTAINS

! **************************************************************************************************
!> \brief Fits the density currently stored in rho_struct with a spin-restricted AO density matrix.
!>        The density matrix is projected in an orthonormal AO representation onto
!>        0 <= P <= 2 and Tr(P) = Ne after every gradient step.
!> \param qs_env QS environment
!> \param rho_struct density structure containing the target density on entry and the fit on exit
!> \param max_iter maximum number of projected-gradient iterations
!> \param eps_rms requested RMS density residual
!> \param step_size initial spectral projected-gradient step
!> \param max_backtrack maximum number of step halvings per iteration
! **************************************************************************************************
   SUBROUTINE fit_constrained_density(qs_env, rho_struct, max_iter, eps_rms, step_size, max_backtrack)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_rho_type), INTENT(INOUT)                   :: rho_struct
      INTEGER, INTENT(IN)                                :: max_iter
      REAL(KIND=dp), INTENT(IN)                          :: eps_rms, step_size
      INTEGER, INTENT(IN)                                :: max_backtrack

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fit_constrained_density'

      INTEGER                                            :: backtrack, handle, info, iteration, nao, &
                                                            nelectron, unit_nr
      LOGICAL                                            :: accepted, converged, do_kpoints, &
                                                            have_previous
      REAL(KIND=dp) :: alpha, direction_derivative, direction_norm2, fitted_grid_trace, &
         matrix_trace, objective, objective_trial, rms, rms_trial, ss, step_length, sy
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues
      REAL(KIND=dp), DIMENSION(:), POINTER               :: tot_rho_r
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type) :: density_orth, direction_orth, eigenvectors, gradient_orth, overlap_chol, &
         previous_density, previous_gradient, trial_orth, work1, work2
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type)                                 :: gradient_ao
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, rho_ao
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_g
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_r3d_rs_type)                               :: residual_rspace, target_rspace
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      CALL timeset(routineN, handle)

      NULLIFY (auxbas_pw_pool, blacs_env, dft_control, fm_struct, gradient_ao%matrix, &
               ks_env, matrix_s, para_env, pw_env, rho_ao, rho_g, rho_r, tot_rho_r)

      CALL get_qs_env(qs_env, blacs_env=blacs_env, dft_control=dft_control, &
                      do_kpoints=do_kpoints, ks_env=ks_env, matrix_s=matrix_s, &
                      nelectron_total=nelectron, para_env=para_env, pw_env=pw_env)

      IF (do_kpoints) CPABORT("Harris CUBE_FIT is currently available only at the Gamma point")
      IF (dft_control%nspins /= 1) THEN
         CPABORT("Harris CUBE_FIT currently requires a spin-restricted calculation")
      END IF
      IF (dft_control%qs_control%gapw) CPABORT("Harris CUBE_FIT currently supports GPW only")

      CALL qs_rho_get(rho_struct, rho_ao=rho_ao, rho_r=rho_r, rho_g=rho_g, tot_rho_r=tot_rho_r)
      CPASSERT(ASSOCIATED(rho_ao) .AND. SIZE(rho_ao) == 1)
      CPASSERT(ASSOCIATED(rho_r) .AND. SIZE(rho_r) == 1)
      CPASSERT(ASSOCIATED(rho_g) .AND. SIZE(rho_g) == 1)

      CALL dbcsr_get_info(matrix_s(1)%matrix, nfullrows_total=nao)
      IF (nelectron < 0 .OR. nelectron > 2*nao) THEN
         CPABORT("Electron count is incompatible with the constrained AO density-matrix bounds")
      END IF

      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
      CALL auxbas_pw_pool%create_pw(target_rspace)
      CALL auxbas_pw_pool%create_pw(residual_rspace)
      CALL pw_copy(rho_r(1), target_rspace)

      ALLOCATE (gradient_ao%matrix)
      CALL dbcsr_create(gradient_ao%matrix, template=matrix_s(1)%matrix, &
                        name="Harris density-fit gradient", matrix_type=dbcsr_type_symmetric)
      CALL dbcsr_copy(gradient_ao%matrix, matrix_s(1)%matrix)
      CALL dbcsr_set(gradient_ao%matrix, 0.0_dp)

      CALL cp_fm_struct_create(fm_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=nao, ncol_global=nao)
      CALL cp_fm_create(overlap_chol, fm_struct, name="density-fit overlap")
      CALL cp_fm_create(density_orth, fm_struct, name="density fit in orthogonal basis")
      CALL cp_fm_create(direction_orth, fm_struct, name="density-fit search direction")
      CALL cp_fm_create(gradient_orth, fm_struct, name="density-fit gradient in orthogonal basis")
      CALL cp_fm_create(previous_density, fm_struct, name="previous density-fit matrix")
      CALL cp_fm_create(previous_gradient, fm_struct, name="previous density-fit gradient")
      CALL cp_fm_create(trial_orth, fm_struct, name="density-fit trial matrix")
      CALL cp_fm_create(eigenvectors, fm_struct, name="density-fit eigenvectors")
      CALL cp_fm_create(work1, fm_struct, name="density-fit work matrix 1")
      CALL cp_fm_create(work2, fm_struct, name="density-fit work matrix 2")
      CALL cp_fm_struct_release(fm_struct)
      ALLOCATE (eigenvalues(nao))

      CALL copy_dbcsr_to_fm(matrix_s(1)%matrix, overlap_chol)
      CALL cp_fm_uplo_to_full(overlap_chol, work1)
      CALL cp_fm_cholesky_decompose(overlap_chol, info_out=info)
      IF (info /= 0) CPABORT("Overlap Cholesky decomposition failed in Harris CUBE_FIT")

      ! X = U P U^T for S = U^T U.  Tr(X) is the electron count.
      CALL copy_dbcsr_to_fm(rho_ao(1)%matrix, work1)
      CALL cp_fm_uplo_to_full(work1, work2)
      CALL cp_fm_cholesky_restore(work1, nao, overlap_chol, work2, "MULTIPLY", pos="LEFT")
      CALL cp_fm_cholesky_restore(work2, nao, overlap_chol, density_orth, &
                                  "MULTIPLY", pos="RIGHT", transa="T")
      CALL project_density_matrix(density_orth, eigenvectors, work1, eigenvalues, &
                                  REAL(nelectron, dp), 2.0_dp)
      CALL orthogonal_to_ao(density_orth, overlap_chol, work1, work2, nao)
      CALL copy_fm_to_dbcsr(work2, rho_ao(1)%matrix)
      CALL evaluate_density_fit(qs_env, ks_env, rho_ao, rho_r, rho_g, tot_rho_r, &
                                target_rspace, residual_rspace, objective, rms)

      logger => cp_get_default_logger()
      IF (logger%para_env%is_source()) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      END IF
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, "(/,T3,A,I0)") "HARRIS| Constrained AO density fit; basis functions: ", nao
         WRITE (unit_nr, "(T3,A,I0)") "HARRIS| Constrained electron count: ", nelectron
         WRITE (unit_nr, "(T3,A,ES16.8)") "HARRIS| Initial density RMS error: ", rms
         WRITE (unit_nr, "(T3,A)") "HARRIS|  Iteration       RMS error          Objective       Step"
      END IF

      converged = rms <= eps_rms
      alpha = step_size
      have_previous = .FALSE.
      iteration = 0
      DO WHILE (iteration < max_iter .AND. .NOT. converged)
         iteration = iteration + 1

         CALL dbcsr_set(gradient_ao%matrix, 0.0_dp)
         CALL integrate_v_rspace(v_rspace=residual_rspace, hmat=gradient_ao, &
                                 qs_env=qs_env, calculate_forces=.FALSE.)
         CALL copy_dbcsr_to_fm(gradient_ao%matrix, work1)
         CALL cp_fm_uplo_to_full(work1, work2)

         ! dF/dX = U^-T (dF/dP) U^-1.
         CALL cp_fm_cholesky_restore(work1, nao, overlap_chol, work2, &
                                     "SOLVE", pos="LEFT", transa="T")
         CALL cp_fm_cholesky_restore(work2, nao, overlap_chol, gradient_orth, &
                                     "SOLVE", pos="RIGHT")

         ! Barzilai-Borwein spectral step for the convex projected-gradient iteration.
         IF (have_previous) THEN
            CALL cp_fm_to_fm(density_orth, work1)
            CALL cp_fm_scale_and_add(1.0_dp, work1, beta=-1.0_dp, matrix_b=previous_density)
            CALL cp_fm_to_fm(gradient_orth, work2)
            CALL cp_fm_scale_and_add(1.0_dp, work2, beta=-1.0_dp, matrix_b=previous_gradient)
            CALL cp_fm_trace(work1, work1, ss)
            CALL cp_fm_trace(work1, work2, sy)
            IF (sy > 100.0_dp*EPSILON(sy)*MAX(1.0_dp, ss)) THEN
               alpha = MIN(100.0_dp*step_size, MAX(1.0E-6_dp*step_size, ss/sy))
            END IF
         END IF
         CALL cp_fm_to_fm(density_orth, previous_density)
         CALL cp_fm_to_fm(gradient_orth, previous_gradient)
         have_previous = .TRUE.

         ! Project once.  Backtracking then stays on the feasible line segment between
         ! the accepted matrix and its projected spectral-gradient point.
         CALL cp_fm_to_fm(density_orth, trial_orth)
         CALL cp_fm_scale_and_add(1.0_dp, trial_orth, beta=-alpha, matrix_b=gradient_orth)
         CALL project_density_matrix(trial_orth, eigenvectors, work1, eigenvalues, &
                                     REAL(nelectron, dp), 2.0_dp)
         CALL cp_fm_to_fm(trial_orth, direction_orth)
         CALL cp_fm_scale_and_add(1.0_dp, direction_orth, beta=-1.0_dp, matrix_b=density_orth)
         CALL cp_fm_trace(gradient_orth, direction_orth, direction_derivative)
         CALL cp_fm_trace(direction_orth, direction_orth, direction_norm2)
         IF (direction_norm2 <= 100.0_dp*EPSILON(direction_norm2)) THEN
            IF (unit_nr > 0) WRITE (unit_nr, "(T3,A)") &
               "HARRIS| Density fit stopped: projected gradient reached a stationary point"
            converged = .TRUE.
            EXIT
         END IF

         accepted = .FALSE.
         step_length = 1.0_dp
         DO backtrack = 0, max_backtrack
            CALL cp_fm_to_fm(density_orth, trial_orth)
            CALL cp_fm_scale_and_add(1.0_dp, trial_orth, beta=step_length, matrix_b=direction_orth)
            CALL orthogonal_to_ao(trial_orth, overlap_chol, work1, work2, nao)
            CALL copy_fm_to_dbcsr(work2, rho_ao(1)%matrix)
            CALL evaluate_density_fit(qs_env, ks_env, rho_ao, rho_r, rho_g, tot_rho_r, &
                                      target_rspace, residual_rspace, objective_trial, rms_trial)
            IF (objective_trial <= objective + 1.0E-4_dp*step_length*direction_derivative) THEN
               accepted = .TRUE.
               EXIT
            END IF
            step_length = 0.5_dp*step_length
         END DO

         IF (.NOT. accepted) THEN
            ! Leave both the AO matrix and the grids at the last accepted point.
            CALL orthogonal_to_ao(density_orth, overlap_chol, work1, work2, nao)
            CALL copy_fm_to_dbcsr(work2, rho_ao(1)%matrix)
            CALL evaluate_density_fit(qs_env, ks_env, rho_ao, rho_r, rho_g, tot_rho_r, &
                                      target_rspace, residual_rspace, objective, rms)
            IF (unit_nr > 0) WRITE (unit_nr, "(T3,A)") &
               "HARRIS| Density fit stopped: projected line search reached a stationary point"
            converged = .TRUE.
            EXIT
         END IF

         CALL cp_fm_to_fm(trial_orth, density_orth)
         objective = objective_trial
         rms = rms_trial
         IF (unit_nr > 0) WRITE (unit_nr, "(T3,A,I7,3ES19.9)") &
            "HARRIS| ", iteration, rms, objective, alpha*step_length
         converged = rms <= eps_rms
      END DO

      CALL dbcsr_dot(rho_ao(1)%matrix, matrix_s(1)%matrix, matrix_trace)
      ! pw_integrate_function performs an MPI reduction and therefore has to be
      ! called collectively, not only by the rank that owns the output unit.
      fitted_grid_trace = pw_integrate_function(rho_r(1), isign=1)
      IF (unit_nr > 0) THEN
         IF (rms <= eps_rms) THEN
            WRITE (unit_nr, "(T3,A,I0,A)") "HARRIS| Density fit reached its RMS target after ", iteration, " iterations"
         ELSE IF (iteration >= max_iter) THEN
            WRITE (unit_nr, "(T3,A,I0,A)") "HARRIS| Density fit stopped after ", iteration, " iterations"
         END IF
         WRITE (unit_nr, "(T3,A,ES16.8)") "HARRIS| Final density RMS error: ", rms
         WRITE (unit_nr, "(T3,A,F20.10)") "HARRIS| AO density-matrix electron count: ", matrix_trace
         WRITE (unit_nr, "(T3,A,F20.10)") "HARRIS| Fitted grid electron count: ", fitted_grid_trace
      END IF

      DEALLOCATE (eigenvalues)
      CALL cp_fm_release(overlap_chol)
      CALL cp_fm_release(density_orth)
      CALL cp_fm_release(direction_orth)
      CALL cp_fm_release(gradient_orth)
      CALL cp_fm_release(previous_density)
      CALL cp_fm_release(previous_gradient)
      CALL cp_fm_release(trial_orth)
      CALL cp_fm_release(eigenvectors)
      CALL cp_fm_release(work1)
      CALL cp_fm_release(work2)
      CALL dbcsr_release(gradient_ao%matrix)
      DEALLOCATE (gradient_ao%matrix)
      CALL auxbas_pw_pool%give_back_pw(residual_rspace)
      CALL auxbas_pw_pool%give_back_pw(target_rspace)

      CALL timestop(handle)

   END SUBROUTINE fit_constrained_density

! **************************************************************************************************
!> \brief Projects a symmetric matrix onto eigenvalue bounds and a prescribed trace.
!> \param matrix ...
!> \param eigenvectors ...
!> \param work ...
!> \param eigenvalues ...
!> \param trace_target ...
!> \param max_occupation ...
! **************************************************************************************************
   SUBROUTINE project_density_matrix(matrix, eigenvectors, work, eigenvalues, trace_target, max_occupation)
      TYPE(cp_fm_type), INTENT(INOUT)                    :: matrix, eigenvectors, work
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: eigenvalues
      REAL(KIND=dp), INTENT(IN)                          :: trace_target, max_occupation

      INTEGER                                            :: i, iteration, n
      REAL(KIND=dp)                                      :: tau, tau_high, tau_low, trace_value

      n = SIZE(eigenvalues)
      CPASSERT(trace_target >= 0.0_dp .AND. trace_target <= max_occupation*REAL(n, dp))

      CALL cp_fm_syevd(matrix, eigenvectors, eigenvalues)

      IF (trace_target <= 0.0_dp) THEN
         eigenvalues = 0.0_dp
      ELSE IF (trace_target >= max_occupation*REAL(n, dp)) THEN
         eigenvalues = max_occupation
      ELSE
         tau_low = MINVAL(eigenvalues) - max_occupation
         tau_high = MAXVAL(eigenvalues)
         DO iteration = 1, 200
            tau = 0.5_dp*(tau_low + tau_high)
            trace_value = SUM(MIN(max_occupation, MAX(0.0_dp, eigenvalues - tau)))
            IF (trace_value > trace_target) THEN
               tau_low = tau
            ELSE
               tau_high = tau
            END IF
            IF (tau_high - tau_low <= 10.0_dp*EPSILON(tau)*MAX(1.0_dp, ABS(tau))) EXIT
         END DO
         tau = 0.5_dp*(tau_low + tau_high)
         DO i = 1, n
            eigenvalues(i) = MIN(max_occupation, MAX(0.0_dp, eigenvalues(i) - tau))
         END DO
      END IF

      CALL cp_fm_to_fm(eigenvectors, work)
      CALL cp_fm_column_scale(work, eigenvalues)
      CALL parallel_gemm('N', 'T', n, n, n, 1.0_dp, work, eigenvectors, 0.0_dp, matrix)

   END SUBROUTINE project_density_matrix

! **************************************************************************************************
!> \brief Transforms X in the orthonormal representation back to P = U^-1 X U^-T.
!> \param density_orth ...
!> \param overlap_chol ...
!> \param work1 ...
!> \param density_ao ...
!> \param nao ...
! **************************************************************************************************
   SUBROUTINE orthogonal_to_ao(density_orth, overlap_chol, work1, density_ao, nao)
      TYPE(cp_fm_type), INTENT(IN)                       :: density_orth, overlap_chol
      TYPE(cp_fm_type), INTENT(INOUT)                    :: work1, density_ao
      INTEGER, INTENT(IN)                                :: nao

      CALL cp_fm_cholesky_restore(density_orth, nao, overlap_chol, work1, "SOLVE", pos="LEFT")
      CALL cp_fm_cholesky_restore(work1, nao, overlap_chol, density_ao, &
                                  "SOLVE", pos="RIGHT", transa="T")

   END SUBROUTINE orthogonal_to_ao

! **************************************************************************************************
!> \brief Collocates the current AO density and evaluates its least-squares residual.
!> \param qs_env ...
!> \param ks_env ...
!> \param rho_ao ...
!> \param rho_r ...
!> \param rho_g ...
!> \param tot_rho_r ...
!> \param target_rspace ...
!> \param residual_rspace ...
!> \param objective ...
!> \param rms ...
! **************************************************************************************************
   SUBROUTINE evaluate_density_fit(qs_env, ks_env, rho_ao, rho_r, rho_g, tot_rho_r, &
                                   target_rspace, residual_rspace, objective, rms)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_g
      REAL(KIND=dp), DIMENSION(:), POINTER               :: tot_rho_r
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: target_rspace
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: residual_rspace
      REAL(KIND=dp), INTENT(OUT)                         :: objective, rms

      MARK_USED(qs_env)

      CALL calculate_rho_elec(matrix_p=rho_ao(1)%matrix, rho=rho_r(1), &
                              rho_gspace=rho_g(1), total_rho=tot_rho_r(1), ks_env=ks_env)
      CALL pw_copy(rho_r(1), residual_rspace)
      CALL pw_axpy(target_rspace, residual_rspace, alpha=-1.0_dp, beta=1.0_dp)
      objective = 0.5_dp*pw_integral_ab(residual_rspace, residual_rspace)
      rms = SQRT(2.0_dp*objective/rho_r(1)%pw_grid%vol)

   END SUBROUTINE evaluate_density_fit

END MODULE qs_density_fit
