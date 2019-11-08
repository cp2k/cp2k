!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief localize wavefunctions
!>      linear response scf
!> \par History
!>      created 07-2005 [MI]
!> \author MI
! **************************************************************************************************
MODULE qs_linres_methods
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_plus_fm_fm_t,&
                                              cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set
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
                                              dbcsr_filter,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_type
   USE hartree_local_methods,           ONLY: Vh_1c_gg_integrals
   USE input_constants,                 ONLY: do_loc_none,&
                                              op_loc_berry,&
                                              ot_precond_none,&
                                              ot_precond_solver_default,&
                                              state_loc_all
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_get_lval,&
                                              section_get_rval,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE lri_environment_types,           ONLY: lri_density_type,&
                                              lri_environment_type,&
                                              lri_kind_type
   USE lri_ks_methods,                  ONLY: calculate_lri_ks_matrix
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE message_passing,                 ONLY: mp_bcast,&
                                              mp_sum
   USE mulliken,                        ONLY: ao_charges
   USE particle_types,                  ONLY: particle_type
   USE preconditioner,                  ONLY: apply_preconditioner,&
                                              make_preconditioner
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_copy,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_poisson_methods,              ONLY: pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_p_type,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_create,&
                                              pw_p_type,&
                                              pw_release,&
                                              pw_retain
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_gapw_densities,               ONLY: prepare_gapw_den
   USE qs_integrate_potential,          ONLY: integrate_v_rspace,&
                                              integrate_v_rspace_diagonal,&
                                              integrate_v_rspace_one_center
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_kpp1_env_types,               ONLY: qs_kpp1_env_type
   USE qs_ks_atom,                      ONLY: update_ks_atom
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
   USE qs_p_env_types,                  ONLY: qs_p_env_type
   USE qs_rho0_ggrid,                   ONLY: integrate_vhg0_rspace
   USE qs_rho_methods,                  ONLY: qs_rho_rebuild,&
                                              qs_rho_update_rho
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_vxc_atom,                     ONLY: calculate_xc_2nd_deriv_atom
   USE string_utilities,                ONLY: xstring
   USE xc,                              ONLY: xc_calc_2nd_deriv,&
                                              xc_prep_2nd_deriv,&
                                              xc_vxc_pw_create
   USE xc_derivatives,                  ONLY: xc_functionals_get_needs
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
   USE xc_rho_set_types,                ONLY: xc_rho_set_create,&
                                              xc_rho_set_get,&
                                              xc_rho_set_release,&
                                              xc_rho_set_type,&
                                              xc_rho_set_update
   USE xtb_ehess,                       ONLY: xtb_coulomb_hessian
   USE xtb_types,                       ONLY: get_xtb_atom_param,&
                                              xtb_atom_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! *** Public subroutines ***
   PUBLIC :: linres_localize, linres_solver
   PUBLIC :: linres_write_restart, linres_read_restart

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

      CHARACTER(LEN=*), PARAMETER :: routineN = 'linres_localize', &
         routineP = moduleN//':'//routineN

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
      ENDIF

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
      TYPE(qs_p_env_type), POINTER                       :: p_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: psi1, h1_psi0, psi0_order
      INTEGER, INTENT(IN)                                :: iounit
      LOGICAL, INTENT(OUT)                               :: should_stop

      CHARACTER(LEN=*), PARAMETER :: routineN = 'linres_solver', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, iter, maxnmo, maxnmo_o, &
                                                            nao, ncol, nmo, nspins
      LOGICAL                                            :: restart
      REAL(dp)                                           :: norm_res, t1, t2
      REAL(dp), DIMENSION(:), POINTER                    :: alpha, beta, tr_pAp, tr_rz0, tr_rz00, &
                                                            tr_rz1
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
      NULLIFY (Ap, r, p, z, buf, mos, tmp_fm_struct, mo_coeff)
      NULLIFY (Sc, chc)

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
      ENDDO
      !
      CALL check_p_env_init(p_env, linres_control, nspins)
      !
      ! allocate the vectors
      ALLOCATE (alpha(nspins), beta(nspins), tr_pAp(nspins), tr_rz0(nspins), tr_rz00(nspins), tr_rz1(nspins), &
                r(nspins), p(nspins), z(nspins), Ap(nspins), mo_coeff_array(nspins))
      DO ispin = 1, nspins
         CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff)
         mo_coeff_array(ispin)%matrix => mo_coeff
      ENDDO
      !
      DO ispin = 1, nspins
         NULLIFY (r(ispin)%matrix, p(ispin)%matrix, z(ispin)%matrix, Ap(ispin)%matrix)
         CALL cp_fm_create(r(ispin)%matrix, psi1(ispin)%matrix%matrix_struct)
         CALL cp_fm_create(p(ispin)%matrix, psi1(ispin)%matrix%matrix_struct)
         CALL cp_fm_create(z(ispin)%matrix, psi1(ispin)%matrix%matrix_struct)
         CALL cp_fm_create(Ap(ispin)%matrix, psi1(ispin)%matrix%matrix_struct)
      ENDDO
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
      ENDDO
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
      ENDDO
      !
      !
      !
      ! header
      IF (iounit > 0) THEN
         WRITE (iounit, "(/,T3,A,T16,A,T25,A,T38,A,T52,A,T72,A,/,T3,A)") &
            "Iteration", "Method", "Restart", "Stepsize", "Convergence", "Time", &
            REPEAT("-", 80)
      ENDIF
      !
      ! orthogonalize x with respect to the psi0
      CALL preortho(psi1, mo_coeff_array, Sc)
      !
      ! build the preconditioner
      IF (linres_control%preconditioner_type /= ot_precond_none) THEN
         IF (p_env%new_preconditioner) THEN
            p_env%os_valid = .FALSE.
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
            ENDDO
            p_env%new_preconditioner = .FALSE.
         ENDIF
      ENDIF
      !
      ! initalization of the linear solver
      !
      ! A * x0
      CALL apply_op(qs_env, p_env, psi0_order, psi1, Ap, chc)
      !
      !
      ! r_0 = b - Ax0
      DO ispin = 1, nspins
         CALL cp_fm_to_fm(h1_psi0(ispin)%matrix, r(ispin)%matrix)
         CALL cp_fm_scale_and_add(-1.0_dp, r(ispin)%matrix, -1.0_dp, Ap(ispin)%matrix)
      ENDDO
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
         ENDDO
         linres_control%flag = "CG"
      ELSE
         !
         ! z_0 = M * r_0
         DO ispin = 1, nspins
            CALL apply_preconditioner(p_env%preconditioner(ispin), r(ispin)%matrix, &
                                      z(ispin)%matrix)
         ENDDO
         linres_control%flag = "PCG"
      ENDIF
      !
      norm_res = 0.0_dp
      DO ispin = 1, nspins
         !
         ! p_0 = z_0
         CALL cp_fm_to_fm(z(ispin)%matrix, p(ispin)%matrix)
         !
         ! trace(r_0 * z_0)
         CALL cp_fm_trace(r(ispin)%matrix, z(ispin)%matrix, tr_rz0(ispin))
         IF (tr_rz0(ispin) .LT. 0.0_dp) CPABORT("tr(r_j*z_j) < 0")
         norm_res = MAX(norm_res, ABS(tr_rz0(ispin))/SQRT(REAL(nao*maxnmo_o, dp)))
      ENDDO
      !
      alpha(:) = 0.0_dp
      restart = .FALSE.
      should_stop = .FALSE.
      iteration: DO iter = 1, linres_control%max_iter
         !
         ! check convergence
         linres_control%converged = .FALSE.
         IF (norm_res .LT. linres_control%eps) THEN
            linres_control%converged = .TRUE.
         ENDIF
         !
         t2 = m_walltime()
         IF (iter .EQ. 1 .OR. MOD(iter, 1) .EQ. 0 .OR. linres_control%converged .OR. restart .OR. should_stop) THEN
            IF (iounit > 0) THEN
               WRITE (iounit, "(T5,I5,T18,A3,T28,L1,T38,1E8.2,T48,F16.10,T68,F8.2)") &
                  iter, linres_control%flag, restart, MAXVAL(alpha), norm_res, t2-t1
               CALL m_flush(iounit)
            ENDIF
         ENDIF
         !
         IF (linres_control%converged) THEN
            IF (iounit > 0) THEN
               WRITE (iounit, "(/,T2,A,I4,A,/)") "The linear solver converged in ", iter, " iterations."
               CALL m_flush(iounit)
            ENDIF
            EXIT iteration
         ELSE IF (should_stop) THEN
            IF (iounit > 0) THEN
               WRITE (iounit, "(/,T2,A,I4,A,/)") "The linear solver did NOT converge! External stop"
               CALL m_flush(iounit)
            END IF
            EXIT iteration
         ENDIF
         !
         ! Max number of iteration reached
         IF (iter == linres_control%max_iter) THEN
            IF (iounit > 0) THEN
               WRITE (iounit, "(/,T2,A/)") &
                  "The linear solver didnt converge! Maximum number of iterations reached."
               CALL m_flush(iounit)
            ENDIF
            linres_control%converged = .FALSE.
         ENDIF
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
            IF (tr_pAp(ispin) .LT. 0.0_dp) THEN
               ! try to fix it by getting rid of the preconditioner
               IF (iter > 1) THEN
                  CALL cp_fm_scale_and_add(beta(ispin), p(ispin)%matrix, -1.0_dp, z(ispin)%matrix)
                  CALL cp_fm_trace(r(ispin)%matrix, r(ispin)%matrix, tr_rz1(ispin))
                  beta(ispin) = tr_rz1(ispin)/tr_rz00(ispin)
                  CALL cp_fm_scale_and_add(beta(ispin), p(ispin)%matrix, 1.0_dp, r(ispin)%matrix)
                  tr_rz0(ispin) = tr_rz1(ispin)
               ELSE
                  CALL cp_fm_to_fm(r(ispin)%matrix, p(ispin)%matrix)
                  CALL cp_fm_trace(r(ispin)%matrix, r(ispin)%matrix, tr_rz0(ispin))
               END IF
               linres_control%flag = "CG"

               CALL apply_op(qs_env, p_env, psi0_order, p, Ap, chc)
               CALL postortho(Ap, mo_coeff_array, Sc)
               CALL cp_fm_trace(Ap(ispin)%matrix, p(ispin)%matrix, tr_pAp(ispin))
               CPABORT("tr(Ap_j*p_j) < 0")
            END IF
            !
            ! alpha = tr(r_j*z_j) / tr(Ap_j*p_j)
            IF (tr_pAp(ispin) .LT. 1.0e-10_dp) THEN
               alpha(ispin) = 1.0_dp
            ELSE
               alpha(ispin) = tr_rz0(ispin)/tr_pAp(ispin)
            ENDIF
            !
            ! x_j+1 = x_j + alpha * p_j
            CALL cp_fm_scale_and_add(1.0_dp, psi1(ispin)%matrix, alpha(ispin), p(ispin)%matrix)
         ENDDO
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
            ENDDO
            CALL postortho(r, mo_coeff_array, Sc)
            !
            restart = .TRUE.
         ELSE
            ! proj Ap onto the virtual subspace
            CALL postortho(Ap, mo_coeff_array, Sc)
            !
            ! r_j+1 = r_j - alpha * Ap_j
            DO ispin = 1, nspins
               CALL cp_fm_scale_and_add(1.0_dp, r(ispin)%matrix, -alpha(ispin), Ap(ispin)%matrix)
            ENDDO
            restart = .FALSE.
         ENDIF
         !
         ! preconditioner
         linres_control%flag = ""
         IF (linres_control%preconditioner_type .EQ. ot_precond_none) THEN
            !
            ! z_j+1 = r_j+1
            DO ispin = 1, nspins
               CALL cp_fm_to_fm(r(ispin)%matrix, z(ispin)%matrix)
            ENDDO
            linres_control%flag = "CG"
         ELSE
            !
            ! z_j+1 = M * r_j+1
            DO ispin = 1, nspins
               CALL apply_preconditioner(p_env%preconditioner(ispin), r(ispin)%matrix, &
                                         z(ispin)%matrix)
            ENDDO
            linres_control%flag = "PCG"
         ENDIF
         !
         norm_res = 0.0_dp
         DO ispin = 1, nspins
            !
            ! tr(r_j+1*z_j+1)
            CALL cp_fm_trace(r(ispin)%matrix, z(ispin)%matrix, tr_rz1(ispin))
            IF (tr_rz1(ispin) .LT. 0.0_dp) CPABORT("tr(r_j+1*z_j+1) < 0")
            norm_res = MAX(norm_res, tr_rz1(ispin)/SQRT(REAL(nao*maxnmo_o, dp)))
            !
            ! beta = tr(r_j+1*z_j+1) / tr(r_j*z_j)
            IF (tr_rz0(ispin) .LT. 1.0e-10_dp) THEN
               beta(ispin) = 0.0_dp
            ELSE
               beta(ispin) = tr_rz1(ispin)/tr_rz0(ispin)
            ENDIF
            !
            ! p_j+1 = z_j+1 + beta * p_j
            CALL cp_fm_scale_and_add(beta(ispin), p(ispin)%matrix, 1.0_dp, z(ispin)%matrix)
            tr_rz00(ispin) = tr_rz0(ispin)
            tr_rz0(ispin) = tr_rz1(ispin)
         ENDDO

         ! Can we exit the SCF loop?
         CALL external_control(should_stop, "LINRES", target_time=qs_env%target_time, &
                               start_time=qs_env%start_time)

      ENDDO iteration
      !
      ! proj psi1
      CALL preortho(psi1, mo_coeff_array, Sc)
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
      ENDDO
      DEALLOCATE (alpha, beta, tr_pAp, tr_rz0, tr_rz00, tr_rz1, r, p, z, Ap, Sc, chc, mo_coeff_array)
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
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_p_env_type), POINTER                       :: p_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: c0, v, Av, chc

      CHARACTER(LEN=*), PARAMETER :: routineN = 'apply_op', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nspins
      REAL(dp)                                           :: chksum
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s, rho1_ao
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(linres_control_type), POINTER                 :: linres_control
      TYPE(qs_rho_type), POINTER                         :: rho

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(v))
      CPASSERT(ASSOCIATED(Av))
      CPASSERT(ASSOCIATED(chc))

      NULLIFY (dft_control, matrix_ks, matrix_s, linres_control, rho1_ao)
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
      ENDDO

      IF (linres_control%do_kernel) THEN

         ! build DM, refill p1, build_dm_response keeps sparse structure
         DO ispin = 1, nspins
            CALL dbcsr_copy(p_env%p1(ispin)%matrix, matrix_s(1)%matrix)
         ENDDO
         CALL build_dm_response(c0, v, p_env%p1)
         DO ispin = 1, nspins
            CALL dbcsr_filter(p_env%p1(ispin)%matrix, linres_control%eps_filter)
         ENDDO

         chksum = 0.0_dp
         DO ispin = 1, nspins
            chksum = chksum+dbcsr_checksum(p_env%p1(ispin)%matrix)
         ENDDO

         ! skip the kernel if the DM is very small
         IF (chksum .GT. 1.0E-14_dp) THEN

            CALL p_env_check_i_alloc(p_env, qs_env)

            CALL qs_rho_get(p_env%rho1, rho_ao=rho1_ao)
            DO ispin = 1, nspins
               CALL dbcsr_copy(rho1_ao(ispin)%matrix, p_env%p1(ispin)%matrix)
            ENDDO

            CALL qs_rho_update_rho(rho_struct=p_env%rho1, &
                                   local_rho_set=p_env%local_rho_set, &
                                   qs_env=qs_env)

            CALL get_qs_env(qs_env, rho=rho) ! that could be called before
            CALL qs_rho_update_rho(rho, qs_env=qs_env) ! that could be called before

            CALL apply_op_2(qs_env, p_env, c0, v, Av, chc)

         ENDIF

      ENDIF

      CALL timestop(handle)

   END SUBROUTINE apply_op

! **************************************************************************************************
!> \brief ...
!> \param c0 ...
!> \param c1 ...
!> \param dm ...
! **************************************************************************************************
   SUBROUTINE build_dm_response(c0, c1, dm)
      !
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: c0, c1
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: dm

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_dm_response', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ispin, ncol, nspins

      CPASSERT(ASSOCIATED(c0))
      CPASSERT(ASSOCIATED(c1))
      CPASSERT(ASSOCIATED(dm))

      nspins = SIZE(dm, 1)

      DO ispin = 1, nspins
         CALL dbcsr_set(dm(ispin)%matrix, 0.0_dp)
         CALL cp_fm_get_info(c0(ispin)%matrix, ncol_global=ncol)
         CALL cp_dbcsr_plus_fm_fm_t(dm(ispin)%matrix, &
                                    matrix_v=c0(ispin)%matrix, &
                                    matrix_g=c1(ispin)%matrix, &
                                    ncol=ncol, alpha=1.0_dp)
         CALL cp_dbcsr_plus_fm_fm_t(dm(ispin)%matrix, &
                                    matrix_v=c1(ispin)%matrix, &
                                    matrix_g=c0(ispin)%matrix, &
                                    ncol=ncol, alpha=1.0_dp)
      ENDDO

   END SUBROUTINE build_dm_response

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
      TYPE(cp_fm_type), POINTER                          :: v, Av
      TYPE(dbcsr_type), POINTER                          :: matrix_ks, matrix_s
      TYPE(cp_fm_type), POINTER                          :: chc

      CHARACTER(LEN=*), PARAMETER :: routineN = 'apply_op_1', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ncol, nrow
      TYPE(cp_fm_type), POINTER                          :: buf

      CALL timeset(routineN, handle)
      !
      CPASSERT(ASSOCIATED(v))
      CPASSERT(ASSOCIATED(Av))
      CPASSERT(ASSOCIATED(matrix_ks))
      CPASSERT(ASSOCIATED(matrix_s))
      CPASSERT(ASSOCIATED(chc))
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

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param p_env ...
!> \param c0 ...
!> \param v ...
!> \param Av ...
!> \param chc ...
! **************************************************************************************************
   SUBROUTINE apply_op_2(qs_env, p_env, c0, v, Av, chc)
      !
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_p_env_type), POINTER                       :: p_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: c0, v, Av, chc

      CHARACTER(len=*), PARAMETER :: routineN = 'apply_op_2', routineP = moduleN//':'//routineN

      TYPE(dft_control_type), POINTER                    :: dft_control

      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control)
      IF (dft_control%qs_control%semi_empirical) THEN
         CPABORT("Linear response not available with SE methods")
      ELSEIF (dft_control%qs_control%dftb) THEN
         CPABORT("Linear response not available with DFTB")
      ELSEIF (dft_control%qs_control%xtb) THEN
         CALL apply_op_2_xtb(qs_env, p_env, c0, v, Av, chc)
      ELSE
         CALL apply_op_2_dft(qs_env, p_env, c0, v, Av, chc)
      END IF

   END SUBROUTINE apply_op_2

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param p_env ...
!> \param c0 ...
!> \param v ...
!> \param Av ...
!> \param chc ...
! **************************************************************************************************
   SUBROUTINE apply_op_2_dft(qs_env, p_env, c0, v, Av, chc)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_p_env_type), POINTER                       :: p_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: c0, v, Av, chc

      CHARACTER(len=*), PARAMETER :: routineN = 'apply_op_2_dft', routineP = moduleN//':'//routineN
      REAL(KIND=dp), PARAMETER                           :: h = 0.001_dp

      INTEGER                                            :: handle, ikind, ispin, ncol, nkind, ns, &
                                                            nspins
      INTEGER, DIMENSION(2, 3)                           :: bo
      LOGICAL                                            :: deriv2_analytic, gapw, gapw_xc, &
                                                            lr_triplet, lrigpw, lsd
      REAL(KIND=dp)                                      :: energy_hartree, energy_hartree_1c, exc, &
                                                            fac
      REAL(KIND=dp), DIMENSION(3, 3)                     :: virial_xc
      REAL(kind=dp), DIMENSION(:, :, :), POINTER         :: rho3, rhoa, rhob
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: k1mat, rho1_ao, rho_ao
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ksmat, psmat
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(linres_control_type), POINTER                 :: linres_control
      TYPE(lri_density_type), POINTER                    :: lri_density
      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(lri_kind_type), DIMENSION(:), POINTER         :: lri_v_int
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: rho1_tot_gspace, v_hartree_gspace, &
                                                            v_hartree_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER :: rho1_g, rho1_g_pw, rho1_r, rho1_r_pw, rho_g, &
         rho_r, tau, tau_pw, v_rspace_new, v_xc, vxc_rho_1, vxc_rho_2, vxc_rho_3, vxc_rho_4
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_rho_type), POINTER                         :: rho, rho1, rho1_xc
      TYPE(section_vals_type), POINTER                   :: input, scf_section, xc_fun_section, &
                                                            xc_section
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type), POINTER                     :: rho1_set

      CALL timeset(routineN, handle)

      NULLIFY (auxbas_pw_pool, pw_pools, pw_env, v_rspace_new, &
               rho1_r, rho1_g_pw, tau_pw, v_xc, rho1_set, rho1_ao, rho_ao, &
               poisson_env, input, scf_section, rho, dft_control, logger, rho1_g)
      logger => cp_get_default_logger()

      energy_hartree = 0.0_dp
      energy_hartree_1c = 0.0_dp

      CPASSERT(ASSOCIATED(c0))
      CPASSERT(ASSOCIATED(v))
      CPASSERT(ASSOCIATED(Av))
      CPASSERT(ASSOCIATED(chc))

      CPASSERT(ASSOCIATED(p_env%kpp1_env))
      CPASSERT(ASSOCIATED(p_env%kpp1))
      rho1 => p_env%rho1
      rho1_xc => p_env%rho1_xc
      CPASSERT(ASSOCIATED(rho1))

      CPASSERT(p_env%kpp1_env%ref_count > 0)

      CALL get_qs_env(qs_env=qs_env, &
                      pw_env=pw_env, &
                      input=input, &
                      rho=rho, &
                      linres_control=linres_control, &
                      dft_control=dft_control)

      lrigpw = dft_control%qs_control%lrigpw
      IF (lrigpw) THEN
         CALL get_qs_env(qs_env, &
                         lri_env=lri_env, &
                         lri_density=lri_density, &
                         atomic_kind_set=atomic_kind_set)
      ENDIF

      CALL qs_rho_get(rho, rho_ao=rho_ao)

      lr_triplet = linres_control%lr_triplet
      CALL kpp1_check_i_alloc(p_env%kpp1_env, qs_env, lr_triplet)
      gapw = dft_control%qs_control%gapw
      gapw_xc = dft_control%qs_control%gapw_xc
      IF (gapw_xc) THEN
         CPASSERT(ASSOCIATED(rho1_xc))
      END IF

      nspins = SIZE(p_env%kpp1)
      lsd = (nspins == 2)

      xc_section => section_vals_get_subs_vals(input, "DFT%XC")
      scf_section => section_vals_get_subs_vals(input, "DFT%SCF")

      p_env%kpp1_env%iter = p_env%kpp1_env%iter+1

      ! gets the tmp grids
      CPASSERT(ASSOCIATED(pw_env))
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                      pw_pools=pw_pools, poisson_env=poisson_env)
      ALLOCATE (v_rspace_new(nspins))
      CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_gspace%pw, &
                             use_data=COMPLEXDATA1D, &
                             in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_rspace%pw, &
                             use_data=REALDATA3D, &
                             in_space=REALSPACE)

      IF (gapw .OR. gapw_xc) &
         CALL prepare_gapw_den(qs_env, p_env%local_rho_set, do_rho0=(.NOT. gapw_xc))

      ! *** calculate the hartree potential on the total density ***
      CALL pw_pool_create_pw(auxbas_pw_pool, rho1_tot_gspace%pw, &
                             use_data=COMPLEXDATA1D, &
                             in_space=RECIPROCALSPACE)

      CALL qs_rho_get(rho1, rho_g=rho1_g)
      CALL pw_copy(rho1_g(1)%pw, rho1_tot_gspace%pw)
      DO ispin = 2, nspins
         CALL pw_axpy(rho1_g(ispin)%pw, rho1_tot_gspace%pw)
      END DO
      IF (gapw) &
         CALL pw_axpy(p_env%local_rho_set%rho0_mpole%rho0_s_gs%pw, rho1_tot_gspace%pw)

      IF (.NOT. (nspins == 1 .AND. lr_triplet)) THEN
         CALL pw_poisson_solve(poisson_env, rho1_tot_gspace%pw, &
                               energy_hartree, &
                               v_hartree_gspace%pw)
         CALL pw_transfer(v_hartree_gspace%pw, v_hartree_rspace%pw)
      ENDIF

      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho1_tot_gspace%pw)

      ! *** calculate the xc potential ***
      IF (gapw_xc) THEN
         CALL qs_rho_get(rho1_xc, rho_r=rho1_r)
      ELSE
         CALL qs_rho_get(rho1, rho_r=rho1_r)
      END IF

      IF (nspins == 1 .AND. (lr_triplet)) THEN

         lsd = .TRUE.
         ALLOCATE (rho1_r_pw(2))
         DO ispin = 1, 2
            NULLIFY (rho1_r_pw(ispin)%pw)
            CALL pw_create(rho1_r_pw(ispin)%pw, rho1_r(1)%pw%pw_grid, &
                           rho1_r(1)%pw%in_use, rho1_r(1)%pw%in_space)
            CALL pw_transfer(rho1_r(1)%pw, rho1_r_pw(ispin)%pw)
         END DO

      ELSE

         ALLOCATE (rho1_r_pw(nspins))
         DO ispin = 1, nspins
            rho1_r_pw(ispin)%pw => rho1_r(ispin)%pw
            CALL pw_retain(rho1_r_pw(ispin)%pw)
         END DO

      END IF

      NULLIFY (tau_pw)

      deriv2_analytic = section_get_lval(xc_section, "2ND_DERIV_ANALYTICAL")

      IF (deriv2_analytic) THEN
         !------!
         ! rho1 !
         !------!
         bo = rho1_r(1)%pw%pw_grid%bounds_local
         ! create the place where to store the argument for the functionals
         CALL xc_rho_set_create(rho1_set, bo, &
                                rho_cutoff=section_get_rval(xc_section, "DENSITY_CUTOFF"), &
                                drho_cutoff=section_get_rval(xc_section, "GRADIENT_CUTOFF"), &
                                tau_cutoff=section_get_rval(xc_section, "TAU_CUTOFF"))

         xc_fun_section => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
         needs = xc_functionals_get_needs(xc_fun_section, lsd, .TRUE.)

         ! calculate the arguments needed by the functionals
         CALL xc_rho_set_update(rho1_set, rho1_r_pw, rho1_g_pw, tau_pw, needs, &
                                section_get_ival(xc_section, "XC_GRID%XC_DERIV"), &
                                section_get_ival(xc_section, "XC_GRID%XC_SMOOTH_RHO"), &
                                auxbas_pw_pool)

         ALLOCATE (v_xc(nspins))
         DO ispin = 1, nspins
            NULLIFY (v_xc(ispin)%pw)
            CALL pw_pool_create_pw(auxbas_pw_pool, v_xc(ispin)%pw, &
                                   use_data=REALDATA3D, &
                                   in_space=REALSPACE)
            CALL pw_zero(v_xc(ispin)%pw)
         END DO

         fac = 0._dp
         IF (nspins == 1) THEN
            IF (lr_triplet) fac = -1.0_dp
         END IF

         CALL xc_calc_2nd_deriv(v_xc, p_env%kpp1_env%deriv_set, p_env%kpp1_env%rho_set, &
                                rho1_set, auxbas_pw_pool, xc_section=xc_section, &
                                tddfpt_fac=fac)

         DO ispin = 1, nspins
            v_rspace_new(ispin)%pw => v_xc(ispin)%pw
         END DO
         DEALLOCATE (v_xc)

         IF (gapw) CALL calculate_xc_2nd_deriv_atom(p_env, qs_env, xc_section, do_tddft=.FALSE., &
                                                    do_triplet=lr_triplet)

         CALL xc_rho_set_release(rho1_set)

      ELSE
         NULLIFY (tau)
         ALLOCATE (v_xc(nspins))
         DO ispin = 1, nspins
            NULLIFY (v_xc(ispin)%pw)
            CALL pw_pool_create_pw(auxbas_pw_pool, v_xc(ispin)%pw, &
                                   use_data=REALDATA3D, &
                                   in_space=REALSPACE)
            CALL pw_zero(v_xc(ispin)%pw)
         END DO
         ! finite differences
         CPASSERT(.NOT. gapw)
         IF (nspins == 2) THEN
            NULLIFY (vxc_rho_1, vxc_rho_2, rho_g)
            ALLOCATE (rho_r(2))
            DO ispin = 1, nspins
               NULLIFY (rho_r(ispin)%pw)
               CALL pw_pool_create_pw(auxbas_pw_pool, rho_r(ispin)%pw, in_space=REALSPACE, use_data=REALDATA3D)
            END DO
            CALL xc_rho_set_get(p_env%kpp1_env%rho_set, rhoa=rhoa, rhob=rhob)
            rho_r(1)%pw%cr3d(:, :, :) = rhoa(:, :, :)+0.5_dp*h*rho1_r(1)%pw%cr3d(:, :, :)
            rho_r(2)%pw%cr3d(:, :, :) = rhob(:, :, :)+0.5_dp*h*rho1_r(2)%pw%cr3d(:, :, :)
            CALL xc_vxc_pw_create(vxc_rho_1, tau_pw, exc, rho_r, rho_g, tau, xc_section, &
                                  auxbas_pw_pool, .FALSE., virial_xc)
            rho_r(1)%pw%cr3d(:, :, :) = rhoa(:, :, :)-0.5_dp*h*rho1_r(1)%pw%cr3d(:, :, :)
            rho_r(2)%pw%cr3d(:, :, :) = rhob(:, :, :)-0.5_dp*h*rho1_r(2)%pw%cr3d(:, :, :)
            CALL xc_vxc_pw_create(vxc_rho_2, tau_pw, exc, rho_r, rho_g, tau, xc_section, &
                                  auxbas_pw_pool, .FALSE., virial_xc)
            v_xc(1)%pw%cr3d(:, :, :) = (vxc_rho_1(1)%pw%cr3d(:, :, :)-vxc_rho_2(1)%pw%cr3d(:, :, :))/h
            v_xc(2)%pw%cr3d(:, :, :) = (vxc_rho_1(2)%pw%cr3d(:, :, :)-vxc_rho_2(2)%pw%cr3d(:, :, :))/h
            DO ispin = 1, nspins
               CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_r(ispin)%pw)
            END DO
            DEALLOCATE (rho_r)
            DO ispin = 1, nspins
               CALL pw_release(vxc_rho_1(ispin)%pw)
               CALL pw_release(vxc_rho_2(ispin)%pw)
            END DO
            DEALLOCATE (vxc_rho_1, vxc_rho_2)
         ELSE IF (nspins == 1 .AND. (lr_triplet)) THEN
            NULLIFY (vxc_rho_1, vxc_rho_2, vxc_rho_3, vxc_rho_4, rho_g)
            ALLOCATE (rho_r(2))
            DO ispin = 1, 2
               NULLIFY (rho_r(ispin)%pw)
               CALL pw_pool_create_pw(auxbas_pw_pool, rho_r(ispin)%pw, in_space=REALSPACE, use_data=REALDATA3D)
            END DO
            CALL xc_rho_set_get(p_env%kpp1_env%rho_set, rhoa=rhoa, rhob=rhob)
            ! K(alpha,alpha)
            rho_r(1)%pw%cr3d(:, :, :) = rhoa(:, :, :)+0.5_dp*h*rho1_r(1)%pw%cr3d(:, :, :)
            rho_r(2)%pw%cr3d(:, :, :) = rhob(:, :, :)
            CALL xc_vxc_pw_create(vxc_rho_1, tau_pw, exc, rho_r, rho_g, tau, xc_section, &
                                  auxbas_pw_pool, .FALSE., virial_xc)
            rho_r(1)%pw%cr3d(:, :, :) = rhoa(:, :, :)-0.5_dp*h*rho1_r(1)%pw%cr3d(:, :, :)
            rho_r(2)%pw%cr3d(:, :, :) = rhob(:, :, :)
            CALL xc_vxc_pw_create(vxc_rho_2, tau_pw, exc, rho_r, rho_g, tau, xc_section, &
                                  auxbas_pw_pool, .FALSE., virial_xc)
            v_xc(1)%pw%cr3d(:, :, :) = (vxc_rho_1(1)%pw%cr3d(:, :, :)-vxc_rho_2(1)%pw%cr3d(:, :, :))/h
            ! K(alpha,beta)
            rho_r(1)%pw%cr3d(:, :, :) = rhoa(:, :, :)
            rho_r(2)%pw%cr3d(:, :, :) = rhob(:, :, :)+0.5_dp*h*rho1_r(1)%pw%cr3d(:, :, :)
            CALL xc_vxc_pw_create(vxc_rho_3, tau_pw, exc, rho_r, rho_g, tau, xc_section, &
                                  auxbas_pw_pool, .FALSE., virial_xc)
            rho_r(1)%pw%cr3d(:, :, :) = rhoa(:, :, :)
            rho_r(2)%pw%cr3d(:, :, :) = rhob(:, :, :)-0.5_dp*h*rho1_r(1)%pw%cr3d(:, :, :)
            CALL xc_vxc_pw_create(vxc_rho_4, tau_pw, exc, rho_r, rho_g, tau, xc_section, &
                                  auxbas_pw_pool, .FALSE., virial_xc)
            v_xc(1)%pw%cr3d(:, :, :) = v_xc(1)%pw%cr3d(:, :, :)- &
                                       (vxc_rho_3(1)%pw%cr3d(:, :, :)-vxc_rho_4(1)%pw%cr3d(:, :, :))/h
            DO ispin = 1, 2
               CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_r(ispin)%pw)
            END DO
            DEALLOCATE (rho_r)
            DO ispin = 1, 2
               CALL pw_release(vxc_rho_1(ispin)%pw)
               CALL pw_release(vxc_rho_2(ispin)%pw)
               CALL pw_release(vxc_rho_3(ispin)%pw)
               CALL pw_release(vxc_rho_4(ispin)%pw)
            END DO
            DEALLOCATE (vxc_rho_1, vxc_rho_2, vxc_rho_3, vxc_rho_4)
         ELSE
            NULLIFY (vxc_rho_1, vxc_rho_2, rho_r, rho_g)
            ALLOCATE (rho_r(1))
            NULLIFY (rho_r(1)%pw)
            CALL pw_pool_create_pw(auxbas_pw_pool, rho_r(1)%pw, in_space=REALSPACE, use_data=REALDATA3D)
            CALL xc_rho_set_get(p_env%kpp1_env%rho_set, rho=rho3)
            rho_r(1)%pw%cr3d(:, :, :) = rho3(:, :, :)+0.5_dp*h*rho1_r(1)%pw%cr3d(:, :, :)
            CALL xc_vxc_pw_create(vxc_rho_1, tau_pw, exc, rho_r, rho_g, tau, xc_section, &
                                  auxbas_pw_pool, .FALSE., virial_xc)
            rho_r(1)%pw%cr3d(:, :, :) = rho3(:, :, :)-0.5_dp*h*rho1_r(1)%pw%cr3d(:, :, :)
            CALL xc_vxc_pw_create(vxc_rho_2, tau_pw, exc, rho_r, rho_g, tau, xc_section, &
                                  auxbas_pw_pool, .FALSE., virial_xc)
            v_xc(1)%pw%cr3d(:, :, :) = (vxc_rho_1(1)%pw%cr3d(:, :, :)-vxc_rho_2(1)%pw%cr3d(:, :, :))/h
            CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_r(1)%pw)
            DEALLOCATE (rho_r)
            CALL pw_release(vxc_rho_1(1)%pw)
            CALL pw_release(vxc_rho_2(1)%pw)
            DEALLOCATE (vxc_rho_1, vxc_rho_2)
         END IF
         DO ispin = 1, nspins
            v_rspace_new(ispin)%pw => v_xc(ispin)%pw
         END DO
         DEALLOCATE (v_xc)
      END IF

      DO ispin = 1, SIZE(rho1_r_pw)
         CALL pw_release(rho1_r_pw(ispin)%pw)
      END DO
      DEALLOCATE (rho1_r_pw)

      !-------------------------------!
      ! Add both hartree and xc terms !
      !-------------------------------!
      DO ispin = 1, nspins

         IF (gapw_xc) THEN
            ! XC and Hartree are integrated separatedly
            ! XC uses the sofft basis set only
            v_rspace_new(ispin)%pw%cr3d = v_rspace_new(ispin)%pw%cr3d* &
                                          v_rspace_new(ispin)%pw%pw_grid%dvol

            IF (nspins == 1) THEN

               IF (.NOT. (lr_triplet)) THEN

                  v_rspace_new(1)%pw%cr3d = 2.0_dp*v_rspace_new(1)%pw%cr3d

               END IF
               CALL qs_rho_get(rho1, rho_ao=rho1_ao)
               ! remove kpp1_env%v_ao and work directly on k_p_p1 ?
               CALL dbcsr_set(p_env%kpp1_env%v_ao(ispin)%matrix, 0.0_dp)
               CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin), &
                                       pmat=rho1_ao(ispin), &
                                       hmat=p_env%kpp1_env%v_ao(ispin), &
                                       qs_env=qs_env, &
                                       calculate_forces=.FALSE., gapw=gapw_xc)

               ! add hartree only for SINGLETS
               IF (.NOT. lr_triplet) THEN
                  v_hartree_rspace%pw%cr3d = v_hartree_rspace%pw%cr3d* &
                                             v_hartree_rspace%pw%pw_grid%dvol
                  v_rspace_new(1)%pw%cr3d = 2.0_dp*v_hartree_rspace%pw%cr3d

                  CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin), &
                                          pmat=rho_ao(ispin), &
                                          hmat=p_env%kpp1_env%v_ao(ispin), &
                                          qs_env=qs_env, &
                                          calculate_forces=.FALSE., gapw=gapw)
               END IF
            ELSE
               ! remove kpp1_env%v_ao and work directly on k_p_p1 ?
               CALL dbcsr_set(p_env%kpp1_env%v_ao(ispin)%matrix, 0.0_dp)
               CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin), &
                                       pmat=rho_ao(ispin), &
                                       hmat=p_env%kpp1_env%v_ao(ispin), &
                                       qs_env=qs_env, &
                                       calculate_forces=.FALSE., gapw=gapw_xc)

               IF (ispin == 1) THEN
                  v_hartree_rspace%pw%cr3d = v_hartree_rspace%pw%cr3d* &
                                             v_hartree_rspace%pw%pw_grid%dvol
               END IF
               v_rspace_new(ispin)%pw%cr3d = v_hartree_rspace%pw%cr3d
               CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin), &
                                       pmat=rho_ao(ispin), &
                                       hmat=p_env%kpp1_env%v_ao(ispin), &
                                       qs_env=qs_env, &
                                       calculate_forces=.FALSE., gapw=gapw)
            END IF

         ELSE

            v_rspace_new(ispin)%pw%cr3d = v_rspace_new(ispin)%pw%cr3d* &
                                          v_rspace_new(ispin)%pw%pw_grid%dvol

            IF (nspins == 1) THEN

               IF (.NOT. (lr_triplet)) THEN
                  v_rspace_new(1)%pw%cr3d = 2.0_dp*v_rspace_new(1)%pw%cr3d

               ENDIF

               ! add hartree only for SINGLETS
               !IF (res_etype == tddfpt_singlet) THEN
               IF (.NOT. lr_triplet) THEN
                  v_hartree_rspace%pw%cr3d = v_hartree_rspace%pw%cr3d* &
                                             v_hartree_rspace%pw%pw_grid%dvol
                  v_rspace_new(1)%pw%cr3d = v_rspace_new(1)%pw%cr3d+ &
                                            2.0_dp*v_hartree_rspace%pw%cr3d
               END IF
            ELSE
               IF (ispin == 1) THEN
                  v_hartree_rspace%pw%cr3d = v_hartree_rspace%pw%cr3d* &
                                             v_hartree_rspace%pw%pw_grid%dvol
               END IF
               v_rspace_new(ispin)%pw%cr3d = v_rspace_new(ispin)%pw%cr3d+ &
                                             v_hartree_rspace%pw%cr3d
            END IF

            ! remove kpp1_env%v_ao and work directly on k_p_p1 ?
            CALL dbcsr_set(p_env%kpp1_env%v_ao(ispin)%matrix, 0.0_dp)

            IF (lrigpw) THEN
               lri_v_int => lri_density%lri_coefs(ispin)%lri_kinds
               CALL get_qs_env(qs_env, para_env=para_env, nkind=nkind)
               DO ikind = 1, nkind
                  lri_v_int(ikind)%v_int = 0.0_dp
               END DO
               CALL integrate_v_rspace_one_center(v_rspace_new(ispin), qs_env, &
                                                  lri_v_int, .FALSE., "LRI_AUX")
               DO ikind = 1, nkind
                  CALL mp_sum(lri_v_int(ikind)%v_int, para_env%group)
               END DO
               ALLOCATE (k1mat(1))
               k1mat(1)%matrix => p_env%kpp1_env%v_ao(ispin)%matrix
               IF (lri_env%exact_1c_terms) THEN
                  CALL integrate_v_rspace_diagonal(v_rspace_new(ispin), k1mat(1)%matrix, &
                                                   rho_ao(ispin)%matrix, qs_env, .FALSE., "ORB")
               END IF
               CALL calculate_lri_ks_matrix(lri_env, lri_v_int, k1mat, atomic_kind_set)
               DEALLOCATE (k1mat)
            ELSE
               CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin), &
                                       pmat=rho_ao(ispin), &
                                       hmat=p_env%kpp1_env%v_ao(ispin), &
                                       qs_env=qs_env, &
                                       calculate_forces=.FALSE., gapw=gapw)
            END IF

         END IF

         CALL dbcsr_copy(p_env%kpp1(ispin)%matrix, p_env%kpp1_env%v_ao(ispin)%matrix)
      END DO

      IF (gapw) THEN

         IF (.NOT. ((nspins == 1 .AND. lr_triplet))) THEN
            CALL Vh_1c_gg_integrals(qs_env, energy_hartree_1c, tddft=.TRUE., do_triplet=lr_triplet, &
                                    p_env=p_env)

            CALL integrate_vhg0_rspace(qs_env, v_hartree_rspace, &
                                       .FALSE., tddft=.TRUE., do_triplet=lr_triplet, &
                                       p_env=p_env)
         END IF

         ! ***  Add single atom contributions to the KS matrix ***
         ! remap pointer
         ns = SIZE(p_env%kpp1)
         ksmat(1:ns, 1:1) => p_env%kpp1(1:ns)
         ns = SIZE(rho_ao)
         psmat(1:ns, 1:1) => rho_ao(1:ns)

         CALL update_ks_atom(qs_env, ksmat, psmat, .FALSE., .TRUE., p_env)

      END IF

      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_gspace%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_rspace%pw)
      DO ispin = 1, nspins
         CALL pw_pool_give_back_pw(auxbas_pw_pool, v_rspace_new(ispin)%pw)
      END DO
      DEALLOCATE (v_rspace_new)

      DO ispin = 1, nspins
         CALL cp_fm_get_info(c0(ispin)%matrix, ncol_global=ncol)
         CALL cp_dbcsr_sm_fm_multiply(p_env%kpp1(ispin)%matrix, &
                                      c0(ispin)%matrix, &
                                      Av(ispin)%matrix, &
                                      ncol=ncol, alpha=1.0_dp, beta=1.0_dp)
      ENDDO
      !
      CALL timestop(handle)
      !
   END SUBROUTINE apply_op_2_dft

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param p_env ...
!> \param c0 ...
!> \param v ...
!> \param Av ...
!> \param chc ...
! **************************************************************************************************
   SUBROUTINE apply_op_2_xtb(qs_env, p_env, c0, v, Av, chc)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_p_env_type), POINTER                       :: p_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: c0, v, Av, chc

      CHARACTER(len=*), PARAMETER :: routineN = 'apply_op_2_xtb', routineP = moduleN//':'//routineN

      INTEGER                                            :: atom_a, handle, iatom, ikind, is, ispin, &
                                                            na, natom, natorb, ncol, nkind, ns, &
                                                            nsgf, nspins
      INTEGER, DIMENSION(25)                             :: lao
      INTEGER, DIMENSION(5)                              :: occ
      LOGICAL                                            :: lr_triplet
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: mcharge, mcharge1
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: aocg, aocg1, charges, charges1
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: p_matrix, rho_ao
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_p, matrix_p1, matrix_s
      TYPE(dbcsr_type), POINTER                          :: s_matrix
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(linres_control_type), POINTER                 :: linres_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho, rho1
      TYPE(xtb_atom_type), POINTER                       :: xtb_kind

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(c0))
      CPASSERT(ASSOCIATED(v))
      CPASSERT(ASSOCIATED(Av))
      CPASSERT(ASSOCIATED(chc))

      CPASSERT(ASSOCIATED(p_env%kpp1_env))
      CPASSERT(ASSOCIATED(p_env%kpp1))

      rho1 => p_env%rho1
      CPASSERT(ASSOCIATED(rho1))

      CPASSERT(p_env%kpp1_env%ref_count > 0)

      CALL get_qs_env(qs_env=qs_env, &
                      pw_env=pw_env, &
                      para_env=para_env, &
                      rho=rho, &
                      linres_control=linres_control, &
                      dft_control=dft_control)

      CALL qs_rho_get(rho, rho_ao=rho_ao)

      lr_triplet = linres_control%lr_triplet
      CPASSERT(.NOT. lr_triplet)

      nspins = SIZE(p_env%kpp1)
      p_env%kpp1_env%iter = p_env%kpp1_env%iter+1

      DO ispin = 1, nspins
         CALL dbcsr_set(p_env%kpp1(ispin)%matrix, 0.0_dp)
      ENDDO

      IF (dft_control%qs_control%xtb_control%coulomb_interaction) THEN
         ! Mulliken charges
         CALL get_qs_env(qs_env, particle_set=particle_set, matrix_s_kp=matrix_s)
         natom = SIZE(particle_set)
         CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
         CALL qs_rho_get(rho1, rho_ao_kp=matrix_p1)
         ALLOCATE (mcharge(natom), charges(natom, 5))
         ALLOCATE (mcharge1(natom), charges1(natom, 5))
         charges = 0.0_dp
         charges1 = 0.0_dp
         CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set)
         nkind = SIZE(atomic_kind_set)
         CALL get_qs_kind_set(qs_kind_set, maxsgf=nsgf)
         ALLOCATE (aocg(nsgf, natom))
         aocg = 0.0_dp
         ALLOCATE (aocg1(nsgf, natom))
         aocg1 = 0.0_dp
         p_matrix => matrix_p(:, 1)
         s_matrix => matrix_s(1, 1)%matrix
         CALL ao_charges(matrix_p, matrix_s, aocg, para_env)
         CALL ao_charges(matrix_p1, matrix_s, aocg1, para_env)
         DO ikind = 1, nkind
            CALL get_atomic_kind(atomic_kind_set(ikind), natom=na)
            CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_kind)
            CALL get_xtb_atom_param(xtb_kind, natorb=natorb, lao=lao, occupation=occ)
            DO iatom = 1, na
               atom_a = atomic_kind_set(ikind)%atom_list(iatom)
               charges(atom_a, :) = REAL(occ(:), KIND=dp)
               DO is = 1, natorb
                  ns = lao(is)+1
                  charges(atom_a, ns) = charges(atom_a, ns)-aocg(is, atom_a)
                  charges1(atom_a, ns) = charges1(atom_a, ns)-aocg1(is, atom_a)
               END DO
            END DO
         END DO
         DEALLOCATE (aocg, aocg1)
         DO iatom = 1, natom
            mcharge(iatom) = SUM(charges(iatom, :))
            mcharge1(iatom) = SUM(charges1(iatom, :))
         END DO
         ! Coulomb Kernel
         CALL xtb_coulomb_hessian(qs_env, p_env%kpp1, charges1, mcharge1, mcharge)
         !
         DEALLOCATE (charges, mcharge, charges1, mcharge1)
      END IF

      DO ispin = 1, nspins
         CALL cp_fm_get_info(c0(ispin)%matrix, ncol_global=ncol)
         CALL cp_dbcsr_sm_fm_multiply(p_env%kpp1(ispin)%matrix, &
                                      c0(ispin)%matrix, &
                                      Av(ispin)%matrix, &
                                      ncol=ncol, alpha=1.0_dp, beta=1.0_dp)
      ENDDO

      CALL timestop(handle)

   END SUBROUTINE apply_op_2_xtb

! **************************************************************************************************
!> \brief ...
!> \param p_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE p_env_check_i_alloc(p_env, qs_env)
      TYPE(qs_p_env_type), POINTER                       :: p_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'p_env_check_i_alloc', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=25)                                  :: name
      INTEGER                                            :: handle, ispin, nspins
      LOGICAL                                            :: gapw_xc
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control

      CALL timeset(routineN, handle)

      NULLIFY (dft_control, matrix_s)

      CPASSERT(ASSOCIATED(p_env))
      CPASSERT(p_env%ref_count > 0)
      CALL get_qs_env(qs_env, dft_control=dft_control)
      gapw_xc = dft_control%qs_control%gapw_xc
      IF (.NOT. ASSOCIATED(p_env%kpp1)) THEN
         CALL get_qs_env(qs_env, matrix_s=matrix_s)
         nspins = dft_control%nspins

         CALL dbcsr_allocate_matrix_set(p_env%kpp1, nspins)
         name = "p_env"//cp_to_string(p_env%id_nr)//"%kpp1-"
         !CALL compress(name,full=.TRUE.)
         DO ispin = 1, nspins
            ALLOCATE (p_env%kpp1(ispin)%matrix)
            CALL dbcsr_copy(p_env%kpp1(ispin)%matrix, matrix_s(1)%matrix, &
                            name=TRIM(name)//ADJUSTL(cp_to_string(ispin)))
            CALL dbcsr_set(p_env%kpp1(ispin)%matrix, 0.0_dp)
         END DO

         CALL qs_rho_rebuild(p_env%rho1, qs_env=qs_env)
         IF (gapw_xc) THEN
            CALL qs_rho_rebuild(p_env%rho1_xc, qs_env=qs_env)
         END IF
      END IF

      IF (.NOT. ASSOCIATED(p_env%rho1)) THEN
         CALL qs_rho_rebuild(p_env%rho1, qs_env=qs_env)
         IF (gapw_xc) THEN
            CALL qs_rho_rebuild(p_env%rho1_xc, qs_env=qs_env)
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE p_env_check_i_alloc

! **************************************************************************************************
!> \brief ...
!> \param kpp1_env ...
!> \param qs_env ...
!> \param lr_triplet ...
! **************************************************************************************************
   SUBROUTINE kpp1_check_i_alloc(kpp1_env, qs_env, lr_triplet)

      TYPE(qs_kpp1_env_type), POINTER                    :: kpp1_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: lr_triplet

      CHARACTER(len=*), PARAMETER :: routineN = 'kpp1_check_i_alloc', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ispin, nspins
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: my_rho_r, rho_r
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: input, xc_section

      NULLIFY (pw_env, auxbas_pw_pool, matrix_s, rho, rho_r, input)

      CPASSERT(ASSOCIATED(kpp1_env))
      CPASSERT(kpp1_env%ref_count > 0)

      CALL get_qs_env(qs_env, pw_env=pw_env, &
                      matrix_s=matrix_s, input=input, rho=rho)

      CALL qs_rho_get(rho, rho_r=rho_r)
      nspins = SIZE(rho_r)

      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      IF (.NOT. ASSOCIATED(kpp1_env%v_rspace)) THEN
         ALLOCATE (kpp1_env%v_rspace(nspins))
         DO ispin = 1, nspins
            CALL pw_pool_create_pw(auxbas_pw_pool, &
                                   kpp1_env%v_rspace(ispin)%pw, &
                                   use_data=REALDATA3D, in_space=REALSPACE)
         END DO
      END IF

      IF (.NOT. ASSOCIATED(kpp1_env%v_ao)) THEN
         CALL dbcsr_allocate_matrix_set(kpp1_env%v_ao, nspins)
         DO ispin = 1, nspins
            ALLOCATE (kpp1_env%v_ao(ispin)%matrix)
            CALL dbcsr_copy(kpp1_env%v_ao(ispin)%matrix, matrix_s(1)%matrix, &
                            name="kpp1%v_ao-"//ADJUSTL(cp_to_string(ispin)))
         END DO
      END IF

      IF (.NOT. ASSOCIATED(kpp1_env%deriv_set)) THEN

         IF (nspins == 1 .AND. lr_triplet) THEN
            ALLOCATE (my_rho_r(2))
            DO ispin = 1, 2
               CALL pw_pool_create_pw(auxbas_pw_pool, my_rho_r(ispin)%pw, &
                                      use_data=rho_r(1)%pw%in_use, in_space=rho_r(1)%pw%in_space)
               my_rho_r(ispin)%pw%cr3d = 0.5_dp*rho_r(1)%pw%cr3d
            END DO
         ELSE
            ALLOCATE (my_rho_r(SIZE(rho_r)))
            DO ispin = 1, SIZE(rho_r)
               my_rho_r(ispin)%pw => rho_r(ispin)%pw
               CALL pw_retain(my_rho_r(ispin)%pw)
            END DO
         END IF

         xc_section => section_vals_get_subs_vals(input, "DFT%XC")

         CALL xc_prep_2nd_deriv(kpp1_env%deriv_set, kpp1_env%rho_set, &
                                my_rho_r, auxbas_pw_pool, &
                                xc_section=xc_section)

         DO ispin = 1, SIZE(my_rho_r)
            CALL pw_release(my_rho_r(ispin)%pw)
         ENDDO
         DEALLOCATE (my_rho_r)
      ENDIF

   END SUBROUTINE kpp1_check_i_alloc

! **************************************************************************************************
!> \brief ...
!> \param v ...
!> \param psi0 ...
!> \param S_psi0 ...
! **************************************************************************************************
   SUBROUTINE preortho(v, psi0, S_psi0)
      !v = (I-PS)v
      !
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: v, psi0, S_psi0

      CHARACTER(LEN=*), PARAMETER :: routineN = 'preortho', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, mp, mt, mv, np, nspins, &
                                                            nt, nv
      TYPE(cp_fm_struct_type), POINTER                   :: tmp_fm_struct
      TYPE(cp_fm_type), POINTER                          :: buf

      CALL timeset(routineN, handle)
      !
      CPASSERT(ASSOCIATED(v))
      CPASSERT(ASSOCIATED(S_psi0))
      CPASSERT(ASSOCIATED(psi0))
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
      ENDDO
      !
      CALL timestop(handle)
      !
   END SUBROUTINE preortho

! **************************************************************************************************
!> \brief ...
!> \param v ...
!> \param psi0 ...
!> \param S_psi0 ...
! **************************************************************************************************
   SUBROUTINE postortho(v, psi0, S_psi0)
      !v = (I-SP)v
      !
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: v, psi0, S_psi0

      CHARACTER(LEN=*), PARAMETER :: routineN = 'postortho', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, mp, mt, mv, np, nspins, &
                                                            nt, nv
      TYPE(cp_fm_struct_type), POINTER                   :: tmp_fm_struct
      TYPE(cp_fm_type), POINTER                          :: buf

      CALL timeset(routineN, handle)
      !
      CPASSERT(ASSOCIATED(v))
      CPASSERT(ASSOCIATED(S_psi0))
      CPASSERT(ASSOCIATED(psi0))
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
      ENDDO
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
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: vec
      INTEGER, INTENT(IN)                                :: ivec
      CHARACTER(LEN=*)                                   :: tag
      INTEGER, INTENT(IN), OPTIONAL                      :: ind

      CHARACTER(LEN=*), PARAMETER :: routineN = 'linres_write_restart', &
         routineP = moduleN//':'//routineN

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
               i_block = MIN(max_block, nmo-i+1)
               CALL cp_fm_get_submatrix(vec(ispin)%matrix, vecbuffer, 1, i, nao, i_block)
               ! doing this in one write would increase efficiency, but breaks RESTART compatibility.
               ! to old ones, and in cases where max_block is different between runs, as might happen during
               ! restarts with a different number of CPUs
               DO j = 1, i_block
                  IF (rst_unit > 0) WRITE (rst_unit) vecbuffer(1:nao, j)
               ENDDO
            ENDDO
         ENDDO

         DEALLOCATE (vecbuffer)

         CALL cp_print_key_finished_output(rst_unit, logger, linres_section, &
                                           "PRINT%RESTART")
      ENDIF

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
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: vec
      INTEGER, INTENT(IN)                                :: ivec
      CHARACTER(LEN=*)                                   :: tag
      INTEGER, INTENT(INOUT), OPTIONAL                   :: ind

      CHARACTER(LEN=*), PARAMETER :: routineN = 'linres_read_restart', &
         routineP = moduleN//':'//routineN

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
         ENDIF
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
         ENDIF
      ENDIF

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
                  i_block = MIN(max_block, nmo-i+1)
                  DO j = 1, i_block
                     IF (rst_unit > 0) READ (rst_unit) vecbuffer(1:nao, j)
                  ENDDO
                  IF (iv .EQ. ivec_tmp) THEN
                     CALL mp_bcast(vecbuffer, source, group)
                     CALL cp_fm_set_submatrix(vec(ispin)%matrix, vecbuffer, 1, i, nao, i_block)
                  ENDIF
               ENDDO
            ENDDO
            IF (ivec .EQ. ivec_tmp) EXIT
         ENDDO

         IF (iostat /= 0) THEN
            IF (iounit > 0) WRITE (iounit, "(T2,A)") &
               "LINRES| Restart file <"//TRIM(ADJUSTL(filename))//"> not found"
         ENDIF

         DEALLOCATE (vecbuffer)

      ENDIF

      IF (para_env%ionode) THEN
         IF (file_exists) CALL close_file(unit_number=rst_unit)
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE linres_read_restart

! **************************************************************************************************

! **************************************************************************************************
!> \brief ...
!> \param p_env ...
!> \param linres_control ...
!> \param nspins ...
! **************************************************************************************************
   SUBROUTINE check_p_env_init(p_env, linres_control, nspins)
      !
      TYPE(qs_p_env_type), POINTER                       :: p_env
      TYPE(linres_control_type), POINTER                 :: linres_control
      INTEGER, INTENT(IN)                                :: nspins

      CHARACTER(LEN=*), PARAMETER :: routineN = 'check_p_env_init', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ispin, ncol, nrow

      p_env%iter = 0
      p_env%only_energy = .FALSE.
      p_env%ls_count = 0

      p_env%ls_pos = 0.0_dp
      p_env%ls_energy = 0.0_dp
      p_env%ls_grad = 0.0_dp
      p_env%gnorm_old = 1.0_dp

      IF (linres_control%preconditioner_type /= ot_precond_none) THEN
         CPASSERT(ASSOCIATED(p_env%preconditioner))
         DO ispin = 1, nspins
            CALL cp_fm_get_info(p_env%PS_psi0(ispin)%matrix, nrow_global=nrow, ncol_global=ncol)
            CPASSERT(nrow == p_env%n_ao(ispin))
            CPASSERT(ncol == p_env%n_mo(ispin))
         ENDDO
      ENDIF

   END SUBROUTINE check_p_env_init

END MODULE qs_linres_methods
