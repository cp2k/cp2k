!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief localize wavefunctions
!>      linear response scf
!> \par History
!>      created 07-2005 [MI]
!> \author MI
! **************************************************************************************************
MODULE qs_linres_methods
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply
   USE cp_external_control,             ONLY: external_control
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add,&
                                              cp_fm_trace
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_set_submatrix,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_fm_vect,                      ONLY: cp_fm_vect_dealloc
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_generate_filename,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_checksum,&
                                              dbcsr_copy,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_type
   USE input_constants,                 ONLY: do_loc_none,&
                                              op_loc_berry,&
                                              ot_precond_none,&
                                              ot_precond_solver_default,&
                                              state_loc_all
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE message_passing,                 ONLY: mp_bcast
   USE preconditioner,                  ONLY: apply_preconditioner,&
                                              make_preconditioner
   USE qs_2nd_kernel_ao,                ONLY: build_dm_response
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_linres_kernel,                ONLY: apply_op_2
   USE qs_linres_types,                 ONLY: linres_control_type
   USE qs_loc_methods,                  ONLY: qs_loc_driver
   USE qs_loc_types,                    ONLY: get_qs_loc_env,&
                                              localized_wfn_control_type,&
                                              qs_loc_env_create,&
                                              qs_loc_env_new_type,&
                                              qs_loc_env_release,&
                                              qs_loc_env_retain
   USE qs_loc_utils,                    ONLY: loc_write_restart,&
                                              qs_loc_control_init,&
                                              qs_loc_init
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_p_env_methods,                ONLY: p_env_check_i_alloc,&
                                              p_env_update_rho
   USE qs_p_env_types,                  ONLY: qs_p_env_type
   USE qs_rho_methods,                  ONLY: qs_rho_update_rho
   USE qs_rho_types,                    ONLY: qs_rho_type
   USE string_utilities,                ONLY: xstring
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! *** Public subroutines ***
   PUBLIC :: linres_localize, linres_solver
   PUBLIC :: linres_write_restart, linres_read_restart
   PUBLIC :: build_dm_response
   PUBLIC :: p_env_check_i_alloc

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_methods'

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Find the centers and spreads of the wfn,
!>      if required apply a localization algorithm
!> \param qs_env ...
!> \param linres_control ...
!> \param nspins ...
!> \param centers_only ...
!> \par History
!>      07.2005 created [MI]
!> \author MI
! **************************************************************************************************
   SUBROUTINE linres_localize(qs_env, linres_control, nspins, centers_only)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(linres_control_type), POINTER                 :: linres_control
      INTEGER, INTENT(IN)                                :: nspins
      LOGICAL, INTENT(IN), OPTIONAL                      :: centers_only

      INTEGER                                            :: iounit, ispin, istate, nmoloc(2)
      LOGICAL                                            :: my_centers_only
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: mos_localized
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(localized_wfn_control_type), POINTER          :: localized_wfn_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(qs_loc_env_new_type), POINTER                 :: qs_loc_env
      TYPE(section_vals_type), POINTER                   :: loc_print_section, loc_section, &
                                                            lr_section

      NULLIFY (logger, lr_section, loc_section, loc_print_section, localized_wfn_control)
      logger => cp_get_default_logger()
      lr_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%LINRES")
      loc_section => section_vals_get_subs_vals(lr_section, "LOCALIZE")
      loc_print_section => section_vals_get_subs_vals(lr_section, "LOCALIZE%PRINT")
      iounit = cp_print_key_unit_nr(logger, lr_section, "PRINT%PROGRAM_RUN_INFO", &
                                    extension=".linresLog")
      my_centers_only = .FALSE.
      IF (PRESENT(centers_only)) my_centers_only = centers_only

      NULLIFY (mos, mo_coeff, qs_loc_env, mos_localized)
      CALL get_qs_env(qs_env=qs_env, mos=mos)
      CALL qs_loc_env_create(qs_loc_env)
      CALL qs_loc_env_retain(qs_loc_env)
      linres_control%qs_loc_env => qs_loc_env
      CALL qs_loc_env_release(qs_loc_env)
      qs_loc_env => linres_control%qs_loc_env
      CALL qs_loc_control_init(qs_loc_env, loc_section, do_homo=.TRUE.)
      CALL get_qs_loc_env(qs_loc_env, localized_wfn_control=localized_wfn_control)

      ALLOCATE (mos_localized(nspins))
      DO ispin = 1, nspins
         CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff)
         CALL cp_fm_create(mos_localized(ispin)%matrix, mo_coeff%matrix_struct)
         CALL cp_fm_to_fm(mo_coeff, mos_localized(ispin)%matrix)
      END DO

      nmoloc(1:2) = 0
      IF (my_centers_only) THEN
         localized_wfn_control%set_of_states = state_loc_all
         localized_wfn_control%localization_method = do_loc_none
         localized_wfn_control%operator_type = op_loc_berry
      END IF

      CALL qs_loc_init(qs_env, qs_loc_env, loc_section, mos_localized=mos_localized, &
                       do_homo=.TRUE.)

      ! The orbital centers are stored in linres_control%localized_wfn_control
      DO ispin = 1, nspins
         CALL qs_loc_driver(qs_env, qs_loc_env, loc_print_section, myspin=ispin, &
                            ext_mo_coeff=mos_localized(ispin)%matrix)
         CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff)
         CALL cp_fm_to_fm(mos_localized(ispin)%matrix, mo_coeff)
      END DO

      CALL loc_write_restart(qs_loc_env, loc_print_section, mos, &
                             mos_localized, do_homo=.TRUE.)
      CALL cp_fm_vect_dealloc(mos_localized)

      ! Write Centers and Spreads on std out
      IF (iounit > 0) THEN
         DO ispin = 1, nspins
            WRITE (iounit, "(/,T2,A,I2)") &
               "WANNIER CENTERS for spin ", ispin
            WRITE (iounit, "(/,T18,A,3X,A)") &
               "--------------- Centers --------------- ", &
               "--- Spreads ---"
            DO istate = 1, SIZE(localized_wfn_control%centers_set(ispin)%array, 2)
               WRITE (iounit, "(T5,A6,I6,2X,3f12.6,5X,f12.6)") &
                  'state ', istate, localized_wfn_control%centers_set(ispin)%array(1:3, istate), &
                  localized_wfn_control%centers_set(ispin)%array(4, istate)
            END DO
         END DO
         CALL m_flush(iounit)
      END IF

   END SUBROUTINE linres_localize

! **************************************************************************************************
!> \brief scf loop to optimize the first order wavefunctions (psi1)
!>      given a perturbation as an operator applied to the ground
!>      state orbitals (h1_psi0)
!> \param p_env ...
!> \param qs_env ...
!> \param psi1 ...
!> \param h1_psi0 ...
!> \param psi0_order ...
!> \param iounit ...
!> \param should_stop ...
!> \par History
!>      07.2005 created [MI]
!> \author MI
! **************************************************************************************************
   SUBROUTINE linres_solver(p_env, qs_env, psi1, h1_psi0, psi0_order, iounit, should_stop)
      TYPE(qs_p_env_type)                                :: p_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_p_type), DIMENSION(:)                   :: psi1, h1_psi0, psi0_order
      INTEGER, INTENT(IN)                                :: iounit
      LOGICAL, INTENT(OUT)                               :: should_stop

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'linres_solver'

      INTEGER                                            :: handle, ispin, iter, maxnmo, maxnmo_o, &
                                                            nao, ncol, nmo, nspins
      LOGICAL                                            :: restart
      REAL(dp)                                           :: alpha, beta, norm_res, t1, t2
      REAL(dp), DIMENSION(:), POINTER                    :: tr_pAp, tr_rz0, tr_rz00, tr_rz1
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: Ap, chc, mo_coeff_array, p, r, Sc, z
      TYPE(cp_fm_struct_type), POINTER                   :: tmp_fm_struct
      TYPE(cp_fm_type), POINTER                          :: buf, mo_coeff
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s, matrix_t
      TYPE(dbcsr_type), POINTER                          :: matrix_x
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(linres_control_type), POINTER                 :: linres_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos

      CALL timeset(routineN, handle)

      NULLIFY (dft_control, linres_control, matrix_s, matrix_t, matrix_ks, para_env)
      NULLIFY (buf, mos, tmp_fm_struct, mo_coeff)

      t1 = m_walltime()

      CALL get_qs_env(qs_env=qs_env, &
                      matrix_ks=matrix_ks, &
                      matrix_s=matrix_s, &
                      kinetic=matrix_t, &
                      dft_control=dft_control, &
                      linres_control=linres_control, &
                      para_env=para_env, &
                      mos=mos)

      nspins = dft_control%nspins
      CALL get_mo_set(mos(1)%mo_set, nao=nao)
      maxnmo = 0
      maxnmo_o = 0
      DO ispin = 1, nspins
         CALL get_mo_set(mos(ispin)%mo_set, nmo=ncol)
         maxnmo = MAX(maxnmo, ncol)
         CALL cp_fm_get_info(psi0_order(ispin)%matrix, ncol_global=ncol)
         maxnmo_o = MAX(maxnmo_o, ncol)
      END DO
      !
      CALL check_p_env_init(p_env, linres_control, nspins)
      !
      ! allocate the vectors
      ALLOCATE (tr_pAp(nspins), tr_rz0(nspins), tr_rz00(nspins), tr_rz1(nspins), &
                r(nspins), p(nspins), z(nspins), Ap(nspins), mo_coeff_array(nspins))
      DO ispin = 1, nspins
         CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff)
         mo_coeff_array(ispin)%matrix => mo_coeff
      END DO
      !
      DO ispin = 1, nspins
         NULLIFY (r(ispin)%matrix, p(ispin)%matrix, z(ispin)%matrix, Ap(ispin)%matrix)
         CALL cp_fm_create(r(ispin)%matrix, psi1(ispin)%matrix%matrix_struct)
         CALL cp_fm_create(p(ispin)%matrix, psi1(ispin)%matrix%matrix_struct)
         CALL cp_fm_create(z(ispin)%matrix, psi1(ispin)%matrix%matrix_struct)
         CALL cp_fm_create(Ap(ispin)%matrix, psi1(ispin)%matrix%matrix_struct)
      END DO
      !
      ! compute S*C0, C0_order'*H*C0_order (this should be done once for all)
      ALLOCATE (chc(nspins), Sc(nspins))
      DO ispin = 1, nspins
         CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff, nmo=nmo)
         CALL cp_fm_create(Sc(ispin)%matrix, mo_coeff%matrix_struct)
         NULLIFY (tmp_fm_struct)
         CALL cp_fm_struct_create(tmp_fm_struct, nrow_global=nmo, &
                                  ncol_global=nmo, para_env=para_env, &
                                  context=mo_coeff%matrix_struct%context)
         CALL cp_fm_create(chc(ispin)%matrix, tmp_fm_struct)
         CALL cp_fm_struct_release(tmp_fm_struct)
      END DO
      !
      DO ispin = 1, nspins
         !
         ! C0_order' * H * C0_order
         mo_coeff => psi0_order(ispin)%matrix
         CALL cp_fm_create(buf, mo_coeff%matrix_struct)
         CALL cp_fm_get_info(mo_coeff, ncol_global=ncol)
         CALL cp_dbcsr_sm_fm_multiply(matrix_ks(ispin)%matrix, mo_coeff, buf, ncol)
         CALL cp_gemm('T', 'N', ncol, ncol, nao, -1.0_dp, mo_coeff, buf, 0.0_dp, chc(ispin)%matrix)
         CALL cp_fm_release(buf)
         !
         ! S * C0
         CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff)
         CALL cp_fm_get_info(mo_coeff, ncol_global=ncol)
         CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, mo_coeff, Sc(ispin)%matrix, ncol)
      END DO
      !
      ! header
      IF (iounit > 0) THEN
         WRITE (iounit, "(/,T3,A,T16,A,T25,A,T38,A,T52,A,T72,A,/,T3,A)") &
            "Iteration", "Method", "Restart", "Stepsize", "Convergence", "Time", &
            REPEAT("-", 78)
      END IF
      !
      ! orthogonalize x with respect to the psi0
      CALL preortho(psi1, mo_coeff_array, Sc)
      !
      ! build the preconditioner
      IF (linres_control%preconditioner_type /= ot_precond_none) THEN
         IF (p_env%new_preconditioner) THEN
            DO ispin = 1, nspins
               IF (ASSOCIATED(matrix_t)) THEN
                  CALL make_preconditioner(p_env%preconditioner(ispin), &
                                           linres_control%preconditioner_type, ot_precond_solver_default, &
                                           matrix_ks(ispin)%matrix, matrix_s(1)%matrix, matrix_t(1)%matrix, &
                                           mos(ispin)%mo_set, linres_control%energy_gap)
               ELSE
                  NULLIFY (matrix_x)
                  CALL make_preconditioner(p_env%preconditioner(ispin), &
                                           linres_control%preconditioner_type, ot_precond_solver_default, &
                                           matrix_ks(ispin)%matrix, matrix_s(1)%matrix, matrix_x, &
                                           mos(ispin)%mo_set, linres_control%energy_gap)
               END IF
            END DO
            p_env%new_preconditioner = .FALSE.
         END IF
      END IF
      !
      ! initialization of the linear solver
      !
      ! A * x0
      CALL apply_op(qs_env, p_env, psi0_order, psi1, Ap, chc)
      !
      !
      ! r_0 = b - Ax0
      DO ispin = 1, nspins
         CALL cp_fm_to_fm(h1_psi0(ispin)%matrix, r(ispin)%matrix)
         CALL cp_fm_scale_and_add(-1.0_dp, r(ispin)%matrix, -1.0_dp, Ap(ispin)%matrix)
      END DO
      !
      ! proj r
      CALL postortho(r, mo_coeff_array, Sc)
      !
      ! preconditioner
      linres_control%flag = ""
      IF (linres_control%preconditioner_type .EQ. ot_precond_none) THEN
         !
         ! z_0 = r_0
         DO ispin = 1, nspins
            CALL cp_fm_to_fm(r(ispin)%matrix, z(ispin)%matrix)
         END DO
         linres_control%flag = "CG"
      ELSE
         !
         ! z_0 = M * r_0
         DO ispin = 1, nspins
            CALL apply_preconditioner(p_env%preconditioner(ispin), r(ispin)%matrix, &
                                      z(ispin)%matrix)
         END DO
         linres_control%flag = "PCG"
      END IF
      !
      DO ispin = 1, nspins
         !
         ! p_0 = z_0
         CALL cp_fm_to_fm(z(ispin)%matrix, p(ispin)%matrix)
         !
         ! trace(r_0 * z_0)
         CALL cp_fm_trace(r(ispin)%matrix, z(ispin)%matrix, tr_rz0(ispin))
      END DO
      IF (SUM(tr_rz0) < 0.0_dp) CPABORT("tr(r_j*z_j) < 0")
      norm_res = ABS(SUM(tr_rz0))/SQRT(REAL(nspins*nao*maxnmo_o, dp))
      !
      alpha = 0.0_dp
      restart = .FALSE.
      should_stop = .FALSE.
      iteration: DO iter = 1, linres_control%max_iter
         !
         ! check convergence
         linres_control%converged = .FALSE.
         IF (norm_res .LT. linres_control%eps) THEN
            linres_control%converged = .TRUE.
         END IF
         !
         t2 = m_walltime()
         IF (iter .EQ. 1 .OR. MOD(iter, 1) .EQ. 0 .OR. linres_control%converged &
             .OR. restart .OR. should_stop) THEN
            IF (iounit > 0) THEN
               WRITE (iounit, "(T5,I5,T18,A3,T28,L1,T38,1E8.2,T48,F16.10,T68,F8.2)") &
                  iter, linres_control%flag, restart, alpha, norm_res, t2 - t1
               CALL m_flush(iounit)
            END IF
         END IF
         !
         IF (linres_control%converged) THEN
            IF (iounit > 0) THEN
               WRITE (iounit, "(T3,A,I4,A)") "The linear solver converged in ", iter, " iterations."
               CALL m_flush(iounit)
            END IF
            EXIT iteration
         ELSE IF (should_stop) THEN
            IF (iounit > 0) THEN
               WRITE (iounit, "(T3,A,I4,A)") "The linear solver did NOT converge! External stop"
               CALL m_flush(iounit)
            END IF
            EXIT iteration
         END IF
         !
         ! Max number of iteration reached
         IF (iter == linres_control%max_iter) THEN
            IF (iounit > 0) THEN
               WRITE (iounit, "(T3,A)") &
                  "The linear solver didn't converge! Maximum number of iterations reached."
               CALL m_flush(iounit)
            END IF
            linres_control%converged = .FALSE.
         END IF
         !
         ! Apply the operators that do not depend on the perturbation
         CALL apply_op(qs_env, p_env, psi0_order, p, Ap, chc)
         !
         ! proj Ap onto the virtual subspace
         CALL postortho(Ap, mo_coeff_array, Sc)
         !
         DO ispin = 1, nspins
            !
            ! tr(Ap_j*p_j)
            CALL cp_fm_trace(Ap(ispin)%matrix, p(ispin)%matrix, tr_pAp(ispin))
         END DO
         !
         ! alpha = tr(r_j*z_j) / tr(Ap_j*p_j)
         IF (SUM(tr_pAp) < 1.0e-10_dp) THEN
            alpha = 1.0_dp
         ELSE
            alpha = SUM(tr_rz0)/SUM(tr_pAp)
         END IF
         DO ispin = 1, nspins
            !
            ! x_j+1 = x_j + alpha * p_j
            CALL cp_fm_scale_and_add(1.0_dp, psi1(ispin)%matrix, alpha, p(ispin)%matrix)
         END DO
         !
         ! need to recompute the residue
         restart = .FALSE.
         IF (MOD(iter, linres_control%restart_every) .EQ. 0) THEN
            !
            ! r_j+1 = b - A * x_j+1
            CALL apply_op(qs_env, p_env, psi0_order, psi1, Ap, chc)
            !
            DO ispin = 1, nspins
               CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff)
               CALL cp_fm_to_fm(h1_psi0(ispin)%matrix, r(ispin)%matrix)
               CALL cp_fm_scale_and_add(-1.0_dp, r(ispin)%matrix, -1.0_dp, Ap(ispin)%matrix)
            END DO
            CALL postortho(r, mo_coeff_array, Sc)
            !
            restart = .TRUE.
         ELSE
            ! proj Ap onto the virtual subspace
            CALL postortho(Ap, mo_coeff_array, Sc)
            !
            ! r_j+1 = r_j - alpha * Ap_j
            DO ispin = 1, nspins
               CALL cp_fm_scale_and_add(1.0_dp, r(ispin)%matrix, -alpha, Ap(ispin)%matrix)
            END DO
            restart = .FALSE.
         END IF
         !
         ! preconditioner
         linres_control%flag = ""
         IF (linres_control%preconditioner_type .EQ. ot_precond_none) THEN
            !
            ! z_j+1 = r_j+1
            DO ispin = 1, nspins
               CALL cp_fm_to_fm(r(ispin)%matrix, z(ispin)%matrix)
            END DO
            linres_control%flag = "CG"
         ELSE
            !
            ! z_j+1 = M * r_j+1
            DO ispin = 1, nspins
               CALL apply_preconditioner(p_env%preconditioner(ispin), r(ispin)%matrix, &
                                         z(ispin)%matrix)
            END DO
            linres_control%flag = "PCG"
         END IF
         !
         DO ispin = 1, nspins
            !
            ! tr(r_j+1*z_j+1)
            CALL cp_fm_trace(r(ispin)%matrix, z(ispin)%matrix, tr_rz1(ispin))
         END DO
         IF (SUM(tr_rz1) < 0.0_dp) CPABORT("tr(r_j+1*z_j+1) < 0")
         norm_res = SUM(tr_rz1)/SQRT(REAL(nspins*nao*maxnmo_o, dp))
         !
         ! beta = tr(r_j+1*z_j+1) / tr(r_j*z_j)
         IF (SUM(tr_rz0) < 1.0e-10_dp) THEN
            beta = 0.0_dp
         ELSE
            beta = SUM(tr_rz1)/SUM(tr_rz0)
         END IF
         DO ispin = 1, nspins
            !
            ! p_j+1 = z_j+1 + beta * p_j
            CALL cp_fm_scale_and_add(beta, p(ispin)%matrix, 1.0_dp, z(ispin)%matrix)
            tr_rz00(ispin) = tr_rz0(ispin)
            tr_rz0(ispin) = tr_rz1(ispin)
         END DO
         !
         ! Can we exit the SCF loop?
         CALL external_control(should_stop, "LINRES", target_time=qs_env%target_time, &
                               start_time=qs_env%start_time)

      END DO iteration
      !
      ! proj psi1
      CALL preortho(psi1, mo_coeff_array, Sc)
      !
      !
      ! clean up
      DO ispin = 1, nspins
         CALL cp_fm_release(r(ispin)%matrix)
         CALL cp_fm_release(p(ispin)%matrix)
         CALL cp_fm_release(z(ispin)%matrix)
         CALL cp_fm_release(Ap(ispin)%matrix)
         !
         CALL cp_fm_release(Sc(ispin)%matrix)
         CALL cp_fm_release(chc(ispin)%matrix)
      END DO
      DEALLOCATE (tr_pAp, tr_rz0, tr_rz00, tr_rz1, r, p, z, Ap, Sc, chc, mo_coeff_array)
      !
      CALL timestop(handle)
      !
   END SUBROUTINE linres_solver

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param p_env ...
!> \param c0 ...
!> \param v ...
!> \param Av ...
!> \param chc ...
! **************************************************************************************************
   SUBROUTINE apply_op(qs_env, p_env, c0, v, Av, chc)
      !
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      TYPE(qs_p_env_type)                                :: p_env
      TYPE(cp_fm_p_type), DIMENSION(:)                   :: c0, v, Av, chc

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'apply_op'

      INTEGER                                            :: handle, ispin, nspins
      REAL(dp)                                           :: chksum
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(linres_control_type), POINTER                 :: linres_control
      TYPE(qs_rho_type), POINTER                         :: rho

      CALL timeset(routineN, handle)

      NULLIFY (dft_control, matrix_ks, matrix_s, linres_control)
      CALL get_qs_env(qs_env=qs_env, &
                      matrix_ks=matrix_ks, &
                      matrix_s=matrix_s, &
                      dft_control=dft_control, &
                      linres_control=linres_control)

      nspins = dft_control%nspins

      ! apply the uncoupled operator
      DO ispin = 1, nspins
         CALL apply_op_1(v(ispin)%matrix, Av(ispin)%matrix, matrix_ks(ispin)%matrix, &
                         matrix_s(1)%matrix, chc(ispin)%matrix)
      END DO

      IF (linres_control%do_kernel) THEN

         ! build DM, refill p1, build_dm_response keeps sparse structure
         DO ispin = 1, nspins
            CALL dbcsr_copy(p_env%p1(ispin)%matrix, matrix_s(1)%matrix)
         END DO
         CALL build_dm_response(c0, v, p_env%p1)

         chksum = 0.0_dp
         DO ispin = 1, nspins
            chksum = chksum + dbcsr_checksum(p_env%p1(ispin)%matrix)
         END DO

         ! skip the kernel if the DM is very small
         IF (chksum .GT. 1.0E-14_dp) THEN

            CALL p_env_check_i_alloc(p_env, qs_env)

            CALL p_env_update_rho(p_env, qs_env)

            CALL get_qs_env(qs_env, rho=rho) ! that could be called before
            CALL qs_rho_update_rho(rho, qs_env=qs_env) ! that could be called before

            DO ispin = 1, nspins
               CALL dbcsr_set(p_env%kpp1(ispin)%matrix, 0.0_dp)
               IF (ASSOCIATED(p_env%kpp1_admm)) CALL dbcsr_set(p_env%kpp1_admm(ispin)%matrix, 0.0_dp)
            END DO

            CALL apply_op_2(qs_env, p_env, c0, Av)

         END IF

      END IF

      CALL timestop(handle)

   END SUBROUTINE apply_op

! **************************************************************************************************
!> \brief ...
!> \param v ...
!> \param Av ...
!> \param matrix_ks ...
!> \param matrix_s ...
!> \param chc ...
! **************************************************************************************************
   SUBROUTINE apply_op_1(v, Av, matrix_ks, matrix_s, chc)
      !
      TYPE(cp_fm_type), INTENT(IN), POINTER              :: v, Av
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_ks, matrix_s
      TYPE(cp_fm_type), INTENT(IN), POINTER              :: chc

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'apply_op_1'

      INTEGER                                            :: handle, ncol, nrow
      TYPE(cp_fm_type), POINTER                          :: buf

      CALL timeset(routineN, handle)
      !
      NULLIFY (buf)
      !
      CALL cp_fm_create(buf, v%matrix_struct)
      !
      CALL cp_fm_get_info(v, ncol_global=ncol, nrow_global=nrow)
      ! H * v
      CALL cp_dbcsr_sm_fm_multiply(matrix_ks, v, Av, ncol)
      ! v * e  (chc already multiplied by -1)
      CALL cp_gemm('N', 'N', nrow, ncol, ncol, 1.0_dp, v, chc, 0.0_dp, buf)
      ! S * ve
      CALL cp_dbcsr_sm_fm_multiply(matrix_s, buf, Av, ncol, alpha=1.0_dp, beta=1.0_dp)
      !Results is H*C1 - S*<iHj>*C1
      !
      CALL cp_fm_release(buf)
      !
      CALL timestop(handle)
      !
   END SUBROUTINE apply_op_1

!MERGE
! **************************************************************************************************
!> \brief ...
!> \param v ...
!> \param psi0 ...
!> \param S_psi0 ...
! **************************************************************************************************
   SUBROUTINE preortho(v, psi0, S_psi0)
      !v = (I-PS)v
      !
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(INOUT)    :: v, psi0, S_psi0

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'preortho'

      INTEGER                                            :: handle, ispin, mp, mt, mv, np, nspins, &
                                                            nt, nv
      TYPE(cp_fm_struct_type), POINTER                   :: tmp_fm_struct
      TYPE(cp_fm_type), POINTER                          :: buf

      CALL timeset(routineN, handle)
      !
      NULLIFY (buf, tmp_fm_struct)
      !
      nspins = SIZE(v, 1)
      !
      DO ispin = 1, nspins
         CALL cp_fm_get_info(v(ispin)%matrix, ncol_global=mv, nrow_global=nv)
         CALL cp_fm_get_info(psi0(ispin)%matrix, ncol_global=mp, nrow_global=np)
         !
         CALL cp_fm_struct_create(tmp_fm_struct, nrow_global=nv, ncol_global=mp, &
                                  para_env=v(ispin)%matrix%matrix_struct%para_env, &
                                  context=v(ispin)%matrix%matrix_struct%context)
         CALL cp_fm_create(buf, tmp_fm_struct)
         CALL cp_fm_struct_release(tmp_fm_struct)
         !
         CALL cp_fm_get_info(buf, ncol_global=mt, nrow_global=nt)
         CPASSERT(nv == np)
         CPASSERT(mt >= mv)
         CPASSERT(mt >= mp)
         CPASSERT(nt == nv)
         !
         ! buf = v' * S_psi0
         CALL cp_gemm('T', 'N', mv, mp, nv, 1.0_dp, v(ispin)%matrix, S_psi0(ispin)%matrix, 0.0_dp, buf)
         ! v = v - psi0 * buf'
         CALL cp_gemm('N', 'T', nv, mv, mp, -1.0_dp, psi0(ispin)%matrix, buf, 1.0_dp, v(ispin)%matrix)
         !
         CALL cp_fm_release(buf)
      END DO
      !
      CALL timestop(handle)
      !
   END SUBROUTINE preortho

! **************************************************************************************************
!> \brief projects first index of v onto the virtual subspace
!> \param v matrix to be projected
!> \param psi0 matrix with occupied orbitals
!> \param S_psi0 matrix containing product of metric and occupied orbitals
! **************************************************************************************************
   SUBROUTINE postortho(v, psi0, S_psi0)
      !v = (I-SP)v
      !
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(INOUT)    :: v, psi0, S_psi0

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'postortho'

      INTEGER                                            :: handle, ispin, mp, mt, mv, np, nspins, &
                                                            nt, nv
      TYPE(cp_fm_struct_type), POINTER                   :: tmp_fm_struct
      TYPE(cp_fm_type), POINTER                          :: buf

      CALL timeset(routineN, handle)
      !
      NULLIFY (buf, tmp_fm_struct)
      !
      nspins = SIZE(v, 1)
      !
      DO ispin = 1, nspins
         CALL cp_fm_get_info(v(ispin)%matrix, ncol_global=mv, nrow_global=nv)
         CALL cp_fm_get_info(psi0(ispin)%matrix, ncol_global=mp, nrow_global=np)
         !
         CALL cp_fm_struct_create(tmp_fm_struct, nrow_global=nv, ncol_global=mp, &
                                  para_env=v(ispin)%matrix%matrix_struct%para_env, &
                                  context=v(ispin)%matrix%matrix_struct%context)
         CALL cp_fm_create(buf, tmp_fm_struct)
         CALL cp_fm_struct_release(tmp_fm_struct)
         !
         CALL cp_fm_get_info(buf, ncol_global=mt, nrow_global=nt)
         CPASSERT(nv == np)
         CPASSERT(mt >= mv)
         CPASSERT(mt >= mp)
         CPASSERT(nt == nv)
         !
         ! buf = v' * psi0
         CALL cp_gemm('T', 'N', mv, mp, nv, 1.0_dp, v(ispin)%matrix, psi0(ispin)%matrix, 0.0_dp, buf)
         ! v = v - S_psi0 * buf'
         CALL cp_gemm('N', 'T', nv, mv, mp, -1.0_dp, S_psi0(ispin)%matrix, buf, 1.0_dp, v(ispin)%matrix)
         !
         CALL cp_fm_release(buf)
      END DO
      !
      CALL timestop(handle)
      !
   END SUBROUTINE postortho

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param linres_section ...
!> \param vec ...
!> \param ivec ...
!> \param tag ...
!> \param ind ...
! **************************************************************************************************
   SUBROUTINE linres_write_restart(qs_env, linres_section, vec, ivec, tag, ind)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: linres_section
      TYPE(cp_fm_p_type), DIMENSION(:)                   :: vec
      INTEGER, INTENT(IN)                                :: ivec
      CHARACTER(LEN=*)                                   :: tag
      INTEGER, INTENT(IN), OPTIONAL                      :: ind

      CHARACTER(LEN=*), PARAMETER :: routineN = 'linres_write_restart'

      CHARACTER(LEN=default_path_length)                 :: filename
      CHARACTER(LEN=default_string_length)               :: my_middle, my_pos, my_status
      INTEGER                                            :: handle, i, i_block, ia, ie, iounit, &
                                                            ispin, j, max_block, nao, nmo, nspins, &
                                                            rst_unit
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: vecbuffer
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(section_vals_type), POINTER                   :: print_key

      NULLIFY (logger, mo_coeff, mos, para_env, print_key, vecbuffer)

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()

      IF (BTEST(cp_print_key_should_output(logger%iter_info, linres_section, "PRINT%RESTART", &
                                           used_print_key=print_key), &
                cp_p_file)) THEN

         iounit = cp_print_key_unit_nr(logger, linres_section, &
                                       "PRINT%PROGRAM_RUN_INFO", extension=".Log")

         CALL get_qs_env(qs_env=qs_env, &
                         mos=mos, &
                         para_env=para_env)

         nspins = SIZE(mos)

         my_status = "REPLACE"
         my_pos = "REWIND"
         CALL XSTRING(tag, ia, ie)
         IF (PRESENT(ind)) THEN
            my_middle = "RESTART-"//tag(ia:ie)//TRIM(ADJUSTL(cp_to_string(ivec)))
         ELSE
            my_middle = "RESTART-"//tag(ia:ie)
            IF (ivec > 1) THEN
               my_status = "OLD"
               my_pos = "APPEND"
            END IF
         END IF
         rst_unit = cp_print_key_unit_nr(logger, linres_section, "PRINT%RESTART", &
                                         extension=".lr", middle_name=TRIM(my_middle), file_status=TRIM(my_status), &
                                         file_position=TRIM(my_pos), file_action="WRITE", file_form="UNFORMATTED")

         filename = cp_print_key_generate_filename(logger, print_key, &
                                                   extension=".lr", middle_name=TRIM(my_middle), my_local=.FALSE.)

         IF (iounit > 0) THEN
            WRITE (UNIT=iounit, FMT="(T2,A)") &
               "LINRES| Writing response functions to the restart file <"//TRIM(ADJUSTL(filename))//">"
         END IF

         !
         ! write data to file
         ! use the scalapack block size as a default for buffering columns
         CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff)
         CALL cp_fm_get_info(mo_coeff, nrow_global=nao, ncol_block=max_block)
         ALLOCATE (vecbuffer(nao, max_block))

         IF (PRESENT(ind)) THEN
            IF (rst_unit > 0) WRITE (rst_unit) ind, ivec, nspins, nao
         ELSE
            IF (rst_unit > 0) WRITE (rst_unit) ivec, nspins, nao
         END IF

         DO ispin = 1, nspins
            CALL cp_fm_get_info(vec(ispin)%matrix, ncol_global=nmo)

            IF (rst_unit > 0) WRITE (rst_unit) nmo

            DO i = 1, nmo, MAX(max_block, 1)
               i_block = MIN(max_block, nmo - i + 1)
               CALL cp_fm_get_submatrix(vec(ispin)%matrix, vecbuffer, 1, i, nao, i_block)
               ! doing this in one write would increase efficiency, but breaks RESTART compatibility.
               ! to old ones, and in cases where max_block is different between runs, as might happen during
               ! restarts with a different number of CPUs
               DO j = 1, i_block
                  IF (rst_unit > 0) WRITE (rst_unit) vecbuffer(1:nao, j)
               END DO
            END DO
         END DO

         DEALLOCATE (vecbuffer)

         CALL cp_print_key_finished_output(rst_unit, logger, linres_section, &
                                           "PRINT%RESTART")
      END IF

      CALL timestop(handle)

   END SUBROUTINE linres_write_restart

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param linres_section ...
!> \param vec ...
!> \param ivec ...
!> \param tag ...
!> \param ind ...
! **************************************************************************************************
   SUBROUTINE linres_read_restart(qs_env, linres_section, vec, ivec, tag, ind)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: linres_section
      TYPE(cp_fm_p_type), DIMENSION(:)                   :: vec
      INTEGER, INTENT(IN)                                :: ivec
      CHARACTER(LEN=*)                                   :: tag
      INTEGER, INTENT(INOUT), OPTIONAL                   :: ind

      CHARACTER(LEN=*), PARAMETER :: routineN = 'linres_read_restart'

      CHARACTER(LEN=default_path_length)                 :: filename
      CHARACTER(LEN=default_string_length)               :: my_middle
      INTEGER :: group, handle, i, i_block, ia, ie, iostat, iounit, ispin, iv, iv1, ivec_tmp, j, &
         max_block, n_rep_val, nao, nao_tmp, nmo, nmo_tmp, nspins, nspins_tmp, rst_unit, source
      LOGICAL                                            :: file_exists
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: vecbuffer
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(section_vals_type), POINTER                   :: print_key

      file_exists = .FALSE.

      CALL timeset(routineN, handle)

      NULLIFY (mos, para_env, logger, print_key, vecbuffer)
      logger => cp_get_default_logger()

      iounit = cp_print_key_unit_nr(logger, linres_section, &
                                    "PRINT%PROGRAM_RUN_INFO", extension=".Log")

      CALL get_qs_env(qs_env=qs_env, &
                      para_env=para_env, &
                      mos=mos)

      nspins = SIZE(mos)
      group = para_env%group
      source = para_env%source !ionode???

      rst_unit = -1
      IF (para_env%ionode) THEN
         CALL section_vals_val_get(linres_section, "WFN_RESTART_FILE_NAME", &
                                   n_rep_val=n_rep_val)

         CALL XSTRING(tag, ia, ie)
         IF (PRESENT(ind)) THEN
            my_middle = "RESTART-"//tag(ia:ie)//TRIM(ADJUSTL(cp_to_string(ivec)))
         ELSE
            my_middle = "RESTART-"//tag(ia:ie)
         END IF

         IF (n_rep_val > 0) THEN
            CALL section_vals_val_get(linres_section, "WFN_RESTART_FILE_NAME", c_val=filename)
            CALL xstring(filename, ia, ie)
            filename = filename(ia:ie)//TRIM(my_middle)//".lr"
         ELSE
            ! try to read from the filename that is generated automatically from the printkey
            print_key => section_vals_get_subs_vals(linres_section, "PRINT%RESTART")
            filename = cp_print_key_generate_filename(logger, print_key, &
                                                      extension=".lr", middle_name=TRIM(my_middle), my_local=.FALSE.)
         END IF
         INQUIRE (FILE=filename, exist=file_exists)
         !
         ! open file
         IF (file_exists) THEN
            CALL open_file(file_name=TRIM(filename), &
                           file_action="READ", &
                           file_form="UNFORMATTED", &
                           file_position="REWIND", &
                           file_status="OLD", &
                           unit_number=rst_unit)

            IF (iounit > 0) WRITE (iounit, "(T2,A)") &
               "LINRES| Reading response wavefunctions from the restart file <"//TRIM(ADJUSTL(filename))//">"
         ELSE
            IF (iounit > 0) WRITE (iounit, "(T2,A)") &
               "LINRES| Restart file  <"//TRIM(ADJUSTL(filename))//"> not found"
         END IF
      END IF

      CALL mp_bcast(file_exists, source, group)

      IF (file_exists) THEN

         CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff)
         CALL cp_fm_get_info(mo_coeff, nrow_global=nao, ncol_block=max_block)

         ALLOCATE (vecbuffer(nao, max_block))
         !
         ! read headers
         IF (PRESENT(ind)) THEN
            iv1 = ivec
         ELSE
            iv1 = 1
         END IF
         DO iv = iv1, ivec

            IF (PRESENT(ind)) THEN
               IF (rst_unit > 0) READ (rst_unit, IOSTAT=iostat) ind, ivec_tmp, nspins_tmp, nao_tmp
               CALL mp_bcast(iostat, source, group)
               CALL mp_bcast(ind, source, group)
            ELSE
               IF (rst_unit > 0) READ (rst_unit, IOSTAT=iostat) ivec_tmp, nspins_tmp, nao_tmp
               CALL mp_bcast(iostat, source, group)
            END IF

            IF (iostat .NE. 0) EXIT
            CALL mp_bcast(ivec_tmp, source, group)
            CALL mp_bcast(nspins_tmp, source, group)
            CALL mp_bcast(nao_tmp, source, group)

            ! check that the number nao, nmo and nspins are
            ! the same as in the current mos
            IF (nspins_tmp .NE. nspins) CPABORT("nspins not consistent")
            IF (nao_tmp .NE. nao) CPABORT("nao not consistent")
            !
            DO ispin = 1, nspins
               CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff)
               CALL cp_fm_get_info(mo_coeff, ncol_global=nmo)
               !
               IF (rst_unit > 0) READ (rst_unit) nmo_tmp
               CALL mp_bcast(nmo_tmp, source, group)
               IF (nmo_tmp .NE. nmo) CPABORT("nmo not consistent")
               !
               ! read the response
               DO i = 1, nmo, MAX(max_block, 1)
                  i_block = MIN(max_block, nmo - i + 1)
                  DO j = 1, i_block
                     IF (rst_unit > 0) READ (rst_unit) vecbuffer(1:nao, j)
                  END DO
                  IF (iv .EQ. ivec_tmp) THEN
                     CALL mp_bcast(vecbuffer, source, group)
                     CALL cp_fm_set_submatrix(vec(ispin)%matrix, vecbuffer, 1, i, nao, i_block)
                  END IF
               END DO
            END DO
            IF (ivec .EQ. ivec_tmp) EXIT
         END DO

         IF (iostat /= 0) THEN
            IF (iounit > 0) WRITE (iounit, "(T2,A)") &
               "LINRES| Restart file <"//TRIM(ADJUSTL(filename))//"> not found"
         END IF

         DEALLOCATE (vecbuffer)

      END IF

      IF (para_env%ionode) THEN
         IF (file_exists) CALL close_file(unit_number=rst_unit)
      END IF

      CALL timestop(handle)

   END SUBROUTINE linres_read_restart

! **************************************************************************************************
!> \brief ...
!> \param p_env ...
!> \param linres_control ...
!> \param nspins ...
! **************************************************************************************************
   SUBROUTINE check_p_env_init(p_env, linres_control, nspins)
      !
      TYPE(qs_p_env_type)                                :: p_env
      TYPE(linres_control_type), INTENT(IN)              :: linres_control
      INTEGER, INTENT(IN)                                :: nspins

      INTEGER                                            :: ispin, ncol, nrow

      IF (linres_control%preconditioner_type /= ot_precond_none) THEN
         CPASSERT(ASSOCIATED(p_env%preconditioner))
         DO ispin = 1, nspins
            CALL cp_fm_get_info(p_env%PS_psi0(ispin)%matrix, nrow_global=nrow, ncol_global=ncol)
            CPASSERT(nrow == p_env%n_ao(ispin))
            CPASSERT(ncol == p_env%n_mo(ispin))
         END DO
      END IF

   END SUBROUTINE check_p_env_init

END MODULE qs_linres_methods
