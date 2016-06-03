!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!
MODULE qs_tddfpt2_methods
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE bibliography,                    ONLY: Iannuzzi2005,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type,&
                                              tddfpt2_control_type
   USE cp_dbcsr_interface,              ONLY: cp_dbcsr_allocate_matrix_set,&
                                              cp_dbcsr_copy,&
                                              cp_dbcsr_deallocate_matrix,&
                                              cp_dbcsr_init_p,&
                                              cp_dbcsr_p_type,&
                                              cp_dbcsr_set,&
                                              cp_dbcsr_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_sm_fm_multiply
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                              cp_fm_scale,&
                                              cp_fm_scale_and_add,&
                                              cp_fm_schur_product,&
                                              cp_fm_symm,&
                                              cp_fm_trace,&
                                              cp_fm_triangular_invert
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose
   USE cp_fm_diag,                      ONLY: choose_eigv_solver
   USE cp_fm_pool_types,                ONLY: cp_fm_pool_type,&
                                              fm_pool_create,&
                                              fm_pool_release
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: &
        cp_fm_create, cp_fm_get_info, cp_fm_get_submatrix, cp_fm_p_type, cp_fm_release, &
        cp_fm_set_all, cp_fm_set_submatrix, cp_fm_to_fm, cp_fm_to_fm_submat, cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                              cp_iterate,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_unit_nr,&
                                              cp_rm_iter_level
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE input_constants,                 ONLY: cholesky_dbcsr,&
                                              cholesky_inverse,&
                                              cholesky_off,&
                                              cholesky_restore
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_get_rval,&
                                              section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: evolt
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_scale,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_poisson_methods,              ONLY: pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_p_type
   USE qs_collocate_density,            ONLY: calculate_wavefunction
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integrate_potential,          ONLY: integrate_v_rspace
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                              deallocate_mo_set,&
                                              get_mo_set,&
                                              init_mo_set,&
                                              mo_set_p_type,&
                                              mo_set_type
   USE qs_rho_methods,                  ONLY: qs_rho_rebuild,&
                                              qs_rho_update_rho
   USE qs_rho_types,                    ONLY: qs_rho_create,&
                                              qs_rho_get,&
                                              qs_rho_release,&
                                              qs_rho_set,&
                                              qs_rho_type
   USE qs_scf_methods,                  ONLY: eigensolver
   USE string_utilities,                ONLY: integer_to_string
   USE util,                            ONLY: sort
   USE xc,                              ONLY: xc_calc_2nd_deriv,&
                                              xc_prep_2nd_deriv
   USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                              xc_dset_release
   USE xc_derivatives,                  ONLY: xc_functionals_get_needs
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
   USE xc_rho_set_types,                ONLY: xc_rho_set_create,&
                                              xc_rho_set_release,&
                                              xc_rho_set_type,&
                                              xc_rho_set_update
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: tddfpt
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt2_methods'

CONTAINS

! **************************************************************************************************
!> \brief Perform TDDFPT calculation.
!> \param qs_env  Quickstep environment
!> \par History
!>    * 05.2016 created [Sergey Chulkov]
!> \note Based on the subroutines apply_op(), iterative_solver(), tddfpt_env_init(), and
!>       tddfpt_env_deallocate(). Looks a bit complex at the moment, but the idea was to eliminate
!>       t_env / p_env wherever possible (thus reducing the amount of spaghetti code).
! **************************************************************************************************
   SUBROUTINE tddfpt(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt', routineP = moduleN//':'//routineN
      INTEGER, PARAMETER                                 :: ispin = 1

      CHARACTER(len=20)                                  :: navail_str, nstates_str
      INTEGER :: energy_unit, handle, istate, iter, jstate, log_unit, nao, niters, nmo_occ_active, &
         nmo_occ_total, nspins_save, nstates, nstates_total
      LOGICAL                                            :: do_hfx
      REAL(kind=dp)                                      :: C_hf, conv, t1, t2
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: evals_last, evals_occ, evals_prev, &
                                                            evals_virt
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: action_matrix
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_ks, matrix_s, rho_ij_ao
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: evecs, evecs_new, K_ij_munu
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: C0_C0T_S, evecs_a_fm, matrix_a_fm, &
                                                            mos_occ, mos_virt, rho_ij_fm
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_ij_r
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_rho_type), POINTER                         :: rho_ij_struct
      TYPE(section_vals_type), POINTER                   :: hfx_section, input, tddfpt_section, &
                                                            xc_section
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control
      TYPE(xc_derivative_set_type), POINTER              :: xc_deriv_set
      TYPE(xc_rho_set_type), POINTER                     :: xc_rho_set

      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()

      CALL cite_reference(Iannuzzi2005)

      NULLIFY (blacs_env, dft_control, input, matrix_ks, matrix_s, mos, para_env, pw_env)
      CALL get_qs_env(qs_env, blacs_env=blacs_env, dft_control=dft_control, &
                      input=input, matrix_ks=matrix_ks, matrix_s=matrix_s, mos=mos, para_env=para_env, pw_env=pw_env)
      tddfpt_control => dft_control%tddfpt2_control

      IF (tddfpt_control%nstates <= 0) THEN
         CALL integer_to_string(tddfpt_control%nstates, nstates_str)
         CALL cp_warn(__LOCATION__, "TDDFPT calculation was requestes for "// &
                      TRIM(nstates_str)//" excited states: nothing to do.")
         CALL timestop(handle)
         RETURN
      END IF

      tddfpt_section => section_vals_get_subs_vals(input, "DFT%TDDFPT2")
      xc_section => section_vals_get_subs_vals(input, "DFT%XC")
      hfx_section => section_vals_get_subs_vals(xc_section, "HF")

      CALL section_vals_get(hfx_section, explicit=do_hfx)
      IF (do_hfx) THEN
         CALL section_vals_val_get(hfx_section, "FRACTION", r_val=C_hf)
         do_hfx = (C_hf /= 0.0_dp)
      END IF

      IF (dft_control%nspins > 1) &
         CPABORT("TDDFPT is currently implemented for restricted molecular orbitals only")

      nspins_save = dft_control%nspins
      dft_control%nspins = 1

      CALL get_mo_set(mos(ispin)%mo_set, nao=nao)

      ! generate all virtual molecular orbitals by diagonalising the Kohn-Sham matrix once again
      !
      ! at this point the scf_env has been partially deallocated;
      ! so do not rely on it and assume scf_env%cholesky_method == cholesky_restore
      NULLIFY (mos_occ, mos_virt)
      CALL tddfpt_lumos(mos_occ=mos_occ, evals_occ=evals_occ, &
                        mos_virt=mos_virt, evals_virt=evals_virt, &
                        mo_set=mos(ispin)%mo_set, matrix_ks=matrix_ks(ispin)%matrix, &
                        matrix_s=matrix_s(1)%matrix, cholesky_method=cholesky_restore, &
                        blacs_env=blacs_env, para_env=para_env)
      ! keep the total number of occupied molecular orbitals for future use
      nmo_occ_total = SIZE(evals_occ)

      ! check we have enough molecular orbitals to contruct the requested number of excited states
      IF (INT(tddfpt_control%nstates, kind=int_8) > SIZE(evals_occ, kind=int_8)*SIZE(evals_virt, kind=int_8)) THEN
         CALL integer_to_string(tddfpt_control%nstates, nstates_str)
         CALL integer_to_string(SIZE(evals_occ)*SIZE(evals_virt), navail_str)
         CALL cp_warn(__LOCATION__, TRIM(nstates_str)//" excited states were requested, however only "// &
                      TRIM(navail_str)//" excited states can be constructed using the atomic basis set given.")
         tddfpt_control%nstates = SIZE(evals_occ)*SIZE(evals_virt)
      END IF

      ! adjust the number of trial vectors
      IF (tddfpt_control%added_states < 0 .OR. &
          INT(tddfpt_control%added_states, kind=int_8)+INT(tddfpt_control%nstates, kind=int_8) > &
          SIZE(evals_occ, kind=int_8)*SIZE(evals_virt, kind=int_8)) THEN

         tddfpt_control%added_states = INT(SIZE(evals_occ, kind=int_8)*SIZE(evals_virt, kind=int_8)- &
                                           INT(tddfpt_control%nstates, kind=int_8))
      END IF

      ! ++ matrix to hold the product : C_0 * C_0^T * S == P_0 * S (needed for orthogonalisation)
      NULLIFY (fm_struct, C0_C0T_S)
      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, para_env=para_env, context=blacs_env)
      CALL cp_fm_create(C0_C0T_S, fm_struct)

      ! ++ initialise xc_rho_set and xc_deriv_set (needed for kernel in adiabatic approximation)
      NULLIFY (pw_env)
      CALL get_qs_env(qs_env, pw_env=pw_env)

      NULLIFY (auxbas_pw_pool)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      NULLIFY (rho_ij_fm)
      CALL cp_fm_create(rho_ij_fm, fm_struct)

      NULLIFY (rho_ij_ao)
      CALL cp_dbcsr_allocate_matrix_set(rho_ij_ao, 1)
      CALL cp_dbcsr_init_p(rho_ij_ao(1)%matrix)
      CALL cp_dbcsr_copy(rho_ij_ao(1)%matrix, matrix_s(1)%matrix)

      NULLIFY (rho_ij_struct, rho_ij_r)
      CALL qs_rho_create(rho_ij_struct)
      CALL qs_rho_set(rho_ij_struct, rho_ao=rho_ij_ao)
      CALL qs_rho_rebuild(rho_ij_struct, qs_env, rebuild_ao=.FALSE., rebuild_grids=.TRUE.)
      CALL qs_rho_get(rho_ij_struct, rho_r=rho_ij_r)

      ! alpha == 2.0 due to spin-restricted density
      CALL cp_gemm('N', 'T', nao, nao, nmo_occ_total, 2.0_dp, mos_occ, mos_occ, 0.0_dp, rho_ij_fm)
      CALL copy_fm_to_dbcsr(rho_ij_fm, rho_ij_ao(1)%matrix, keep_sparsity=.TRUE.)
      CALL qs_rho_update_rho(rho_ij_struct, qs_env)

      NULLIFY (xc_rho_set, xc_deriv_set)

      CALL xc_prep_2nd_deriv(xc_deriv_set, xc_rho_set, rho_ij_r, auxbas_pw_pool, xc_section=xc_section)

      CALL qs_rho_release(rho_ij_struct)
      CALL cp_fm_release(rho_ij_fm)

      CALL cp_fm_struct_release(fm_struct)

      log_unit = cp_print_key_unit_nr(logger, tddfpt_section, "PRINT%GUESS_VECTORS", extension=".tddfptLog")

      nstates = tddfpt_control%nstates
      ALLOCATE (evals_prev(nstates))

      ! allocate storage for trial vectors and guess them using virtual orbitals from the ground-state DFT calculation
      CALL tddfpt_guess(nstates_active=tddfpt_control%nstates+tddfpt_control%added_states, &
                        evecs=evecs, evals=evals_prev, mos_occ=mos_occ, evals_occ=evals_occ, &
                        added_occ_mos=tddfpt_control%added_occ_mos, mos_virt=mos_virt, &
                        evals_virt=evals_virt, threshold=tddfpt_control%degenerate_eps, &
                        add_degenerate=tddfpt_control%add_degenerate, blacs_env=blacs_env, &
                        para_env=para_env, log_unit=log_unit)
      ! NOTE: at this point all occupied MOs outside the active space have been dropped, so the variables
      !       mos_occ / evals_occ hold ONLY ACTIVE occupied MOs and their orbital energies.
      !       It helps to save huge amount of computational time by limiting the number of
      !       four-centre exact exchange integrals that are actually contribute to the action matrix.

      nstates_total = SIZE(evecs)
      ALLOCATE (evals_last(nstates_total))

      CALL cp_print_key_finished_output(log_unit, logger, tddfpt_section, "PRINT%GUESS_VECTORS")

      ! the number of occupied molecular orbitals in the active space
      nmo_occ_active = SIZE(evals_occ)

      ! ++ temporary matrices needed to update trial wectors
      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nmo_occ_active, para_env=para_env, context=blacs_env)
      ALLOCATE (evecs_new(nstates_total))
      DO istate = 1, nstates_total
         NULLIFY (evecs_new(istate)%matrix)
         CALL cp_fm_create(evecs_new(istate)%matrix, fm_struct)
      END DO
      CALL cp_fm_struct_release(fm_struct)

      ! ++ action matrix and its eigenvectors
      ALLOCATE (action_matrix(nstates_total, nstates_total))
      NULLIFY (matrix_a_fm, evecs_a_fm)
      CALL cp_fm_struct_create(fm_struct, nrow_global=nstates_total, ncol_global=nstates_total, &
                               para_env=para_env, context=blacs_env)
      CALL cp_fm_create(matrix_a_fm, fm_struct)
      CALL cp_fm_create(evecs_a_fm, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      ! ++ storage for half-transformed four-centre integrals (i j || mu nu),
      ! where i <= j are occupied MOs, and mu, nu are atomic basis functions
      IF (do_hfx) THEN
         CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, para_env=para_env, context=blacs_env)

         jstate = (nmo_occ_active*(nmo_occ_active+1))/2

         ALLOCATE (K_ij_munu(jstate))
         DO istate = 1, jstate
            NULLIFY (K_ij_munu(istate)%matrix)
            CALL cp_fm_create(K_ij_munu(istate)%matrix, fm_struct)
         END DO

         CALL cp_fm_struct_release(fm_struct)
      END IF

      ! compute P_0 * S
      CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, mos_occ, evecs_new(1)%matrix, &
                                   ncol=nmo_occ_active, alpha=1.0_dp, beta=0.0_dp)
      CALL cp_gemm('N', 'T', nao, nao, nmo_occ_active, 1.0_dp, mos_occ, evecs_new(1)%matrix, 0.0_dp, C0_C0T_S)

      ! start the wave-function optimization
      log_unit = cp_print_key_unit_nr(logger, tddfpt_section, "PRINT%ITERATION_INFO", extension=".tddfptLog")
      energy_unit = cp_print_key_unit_nr(logger, tddfpt_section, "PRINT%DETAILED_ENERGY", extension=".tddfptLog")

      IF (log_unit > 0) THEN
         WRITE (log_unit, '(/,1X,A,/)') 'TDDFT WAVEFUNCTION OPTIMIZATION'
         WRITE (log_unit, '(5X,A,T15,A,T24,A,T40,A)') "Step", "Time", "Convergence", "Conv. states"
         WRITE (log_unit, '(1X,50("-"))')
      END IF

      CALL cp_add_iter_level(logger%iter_info, "TDDFT_SCF")

      niters = tddfpt_control%niters
      t1 = m_walltime()

      ! compute and store half-transformed four-centre integrals needed for exact exchange
      IF (do_hfx) &
         CALL tddfpt_compute_hfx(K_ij_munu=K_ij_munu, mos_occ=mos_occ, qs_env=qs_env, matrix_s=matrix_s(1)%matrix)

      DO iter = 1, niters
         CALL cp_iterate(logger%iter_info, last=.FALSE., iter_nr=iter)

         ! initialise the upper-triangular 'action_matrix' by computing (for i <= j)
         ! trace { C_{1,i}^T * ( KS * C_{1,j} - S * C_{1,j} * diag(evals_occ) ) }
         CALL tddfpt_init_by_energy_diff(a_kv=action_matrix, evecs=evecs, evals_occ=evals_occ, &
                                         matrix_ks=matrix_ks(ispin)%matrix, matrix_s=matrix_s(1)%matrix)

         ! add contributions from coulomb terms
         ! alpha == 2 due to singlet state
         CALL tddfpt_apply_coulomb(a_kv=action_matrix, alpha=2.0_dp, evecs=evecs, mos_occ=mos_occ, &
                                   qs_env=qs_env, matrix_s=matrix_s(1)%matrix)

         ! add contributions from the adiabatic TDDFT kernel
         ! alpha == 2 due to relation K_alpha,alpha == K_alpha,beta for singlet states
         CALL tddfpt_apply_kernel(a_kv=action_matrix, alpha=2.0_dp, evecs=evecs, mos_occ=mos_occ, &
                                  xc_rho_set=xc_rho_set, xc_deriv_set=xc_deriv_set, &
                                  qs_env=qs_env, matrix_s=matrix_s(1)%matrix, xc_section=xc_section)

         ! add contribution from exact exchange terms
         IF (do_hfx) &
            CALL tddfpt_apply_hfx(a_kv=action_matrix, K_ij_munu=K_ij_munu, alpha=C_hf, evecs=evecs, qs_env=qs_env)

         ! solve the eigenproblem for action matrix
         CALL cp_fm_set_submatrix(matrix_a_fm, action_matrix)
         CALL choose_eigv_solver(matrix_a_fm, evecs_a_fm, evals_last)

         t2 = m_walltime()

         IF (energy_unit > 0) THEN
            WRITE (energy_unit, '(/,4X,A,T14,A,T36,A)') "State", "Exc. energy (eV)", "Convergence (eV)"
            DO istate = 1, nstates
               WRITE (energy_unit, '(1X,I8,T12,F14.7,T38,ES11.4)') istate, &
                  evals_last(istate)*evolt, (evals_last(istate)-evals_prev(istate))*evolt
            END DO
            WRITE (energy_unit, *)
            CALL m_flush(energy_unit)
         END IF

         conv = MAXVAL(ABS(evals_last(1:nstates)-evals_prev(1:nstates)))

         IF (log_unit > 0) THEN
            jstate = 0
            DO istate = 1, nstates
               IF (ABS(evals_last(istate)-evals_prev(istate)) <= tddfpt_control%conv) &
                  jstate = jstate+1
            END DO

            WRITE (log_unit, '(1X,I8,T12,F7.1,T24,ES11.4,T42,I8)') iter, t2-t1, conv, jstate
            CALL m_flush(log_unit)
         END IF

         t1 = t2
         evals_prev(1:nstates) = evals_last(1:nstates)

         IF (conv <= tddfpt_control%conv) THEN
            CALL cp_iterate(logger%iter_info, last=.TRUE., iter_nr=iter)
            EXIT
         END IF

         CALL cp_fm_get_submatrix(evecs_a_fm, action_matrix)
         DO istate = 1, nstates_total
            CALL cp_fm_set_all(evecs_new(istate)%matrix, 0.0_dp)

            DO jstate = 1, nstates_total
               CALL cp_fm_scale_and_add(1.0_dp, evecs_new(istate)%matrix, &
                                        action_matrix(jstate, istate), evecs(jstate)%matrix)
            END DO
         END DO

         DO istate = 1, nstates_total
            CALL tddfpt_orthogonalize_psi1_psi0(evecs_new(istate)%matrix, C0_C0T_S)
            CALL tddfpt_orthonormalize_psi1_psi1(evecs_new, matrix_s(1)%matrix, istate)

            CALL cp_fm_to_fm(evecs_new(istate)%matrix, evecs(istate)%matrix)
         END DO

      END DO

      CALL cp_rm_iter_level(logger%iter_info, "TDDFT_SCF")

      IF (log_unit > 0) THEN
         IF (iter <= niters) THEN
            CALL integer_to_string(iter, nstates_str)
            WRITE (log_unit, '(/,1X,A)') "*** TDDFT run converged in "//TRIM(nstates_str)//" iteration(s) ***"
         ELSE
            CALL integer_to_string(iter-1, nstates_str)
            WRITE (log_unit, '(/,1X,A)') "*** TDDFT run did NOT converge after "//TRIM(nstates_str)//" iteration(s) ***"
         END IF
      END IF

      CALL cp_print_key_finished_output(energy_unit, logger, tddfpt_section, "PRINT%DETAILED_ENERGY")
      CALL cp_print_key_finished_output(log_unit, logger, tddfpt_section, "PRINT%ITERATION_INFO")

      ! -- clean up all useless stuff
      IF (ALLOCATED(K_ij_munu)) THEN
         DO istate = SIZE(K_ij_munu), 1, -1
            CALL cp_fm_release(K_ij_munu(istate)%matrix)
         END DO
         DEALLOCATE (K_ij_munu)
      END IF

      CALL xc_dset_release(xc_deriv_set)
      CALL xc_rho_set_release(xc_rho_set)

      CALL cp_fm_release(C0_C0T_S)

      CALL cp_fm_release(evecs_a_fm)
      CALL cp_fm_release(matrix_a_fm)
      DEALLOCATE (action_matrix)
      DEALLOCATE (evals_prev)

      DO istate = SIZE(evecs_new), 1, -1
         CALL cp_fm_release(evecs_new(istate)%matrix)
      END DO
      DEALLOCATE (evecs_new)

      ! print summary information
      log_unit = cp_logger_get_default_io_unit()
      CALL tddfpt_print_summary(log_unit=log_unit, nstates=nstates, evecs=evecs, evals=evals_last, &
                                nmo_occ_total=nmo_occ_total, mos_virt=mos_virt, matrix_s=matrix_s, &
                                min_amplitude=tddfpt_control%min_excitation_amplitude, mos_occ=mos_occ, &
                                evals_virt=evals_virt, evals_occ=evals_occ)

      DEALLOCATE (evals_virt)
      DEALLOCATE (evals_occ)
      CALL cp_fm_release(mos_occ)

      DEALLOCATE (evals_last)

      DO istate = SIZE(evecs), 1, -1
         CALL cp_fm_release(evecs(istate)%matrix)
      END DO
      DEALLOCATE (evecs)

      CALL cp_fm_release(mos_virt)

      dft_control%nspins = nspins_save
      CALL timestop(handle)
   END SUBROUTINE tddfpt

! **************************************************************************************************
!> \brief Generate all virtual molecular orbitals for a given spin by diagonalising
!>        the corresponding Kohn-Sham matrix.
!> \param mos_occ         occupied molecular orbitals (allocated on exit)
!> \param evals_occ       occupied orbital energies (allocated on exit)
!> \param mos_virt        virtual molecular orbitals (allocated on exit)
!> \param evals_virt      virtual orbital energies (allocated on exit)
!> \param mo_set          ground state molecular orbitals for a given spin
!> \param matrix_ks       Kohn-Sham matrix for a given spin
!> \param matrix_s        overlap matrix
!> \param cholesky_method Cholesky method to compute the inverse overlap matrix
!> \param blacs_env       BLACS environment for newly allocated distributed matrices
!> \param para_env        parallel environment for newly allocated distributed matrices
!> \par History
!>    * 05.2016 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_lumos(mos_occ, evals_occ, mos_virt, evals_virt, &
                           mo_set, matrix_ks, matrix_s, cholesky_method, blacs_env, para_env)
      TYPE(cp_fm_type), POINTER                          :: mos_occ
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(inout)                                   :: evals_occ
      TYPE(cp_fm_type), POINTER                          :: mos_virt
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(inout)                                   :: evals_virt
      TYPE(mo_set_type), POINTER                         :: mo_set
      TYPE(cp_dbcsr_type), POINTER                       :: matrix_ks, matrix_s
      INTEGER, INTENT(in)                                :: cholesky_method
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_lumos', routineP = moduleN//':'//routineN

      INTEGER                                            :: cholesky_method_inout, handle, homo, &
                                                            iocc, nao, nelectrons, nmo_occ_actual, &
                                                            nmo_scf, nmo_virt_actual
      REAL(kind=dp)                                      :: maxocc
      REAL(kind=dp), DIMENSION(:), POINTER               :: mo_evals_extended, mo_occ_extended, &
                                                            mo_occ_scf
      TYPE(cp_fm_pool_type), POINTER                     :: ao_ao_fm_pool
      TYPE(cp_fm_struct_type), POINTER                   :: ao_ao_fm_struct, ao_mo_occ_fm_struct, &
                                                            ao_mo_virt_fm_struct
      TYPE(cp_fm_type), POINTER                          :: matrix_ks_fm, mo_coeff_extended, &
                                                            ortho_fm, work_fm
      TYPE(mo_set_type), POINTER                         :: mos_extended

      CALL timeset(routineN, handle)

      CPASSERT((.NOT. ASSOCIATED(mos_occ)) .AND. (.NOT. ASSOCIATED(mos_virt)))
      CPASSERT((.NOT. ALLOCATED(evals_occ)) .AND. (.NOT. ALLOCATED(evals_virt)))

      CALL get_mo_set(mo_set, nao=nao, nmo=nmo_scf, homo=homo, maxocc=maxocc, nelectron=nelectrons, occupation_numbers=mo_occ_scf)

      nmo_occ_actual = homo
      nmo_virt_actual = nao-homo

      IF (nmo_virt_actual <= 0) &
         CALL cp_abort(__LOCATION__, &
                       'Unable to generate virtual molecular orbitals. Please use a larger atomic basis set.')

      ! structures of distributed dense matrices
      NULLIFY (ao_ao_fm_struct, ao_mo_occ_fm_struct, ao_mo_virt_fm_struct)
      CALL cp_fm_struct_create(ao_ao_fm_struct, nrow_global=nao, ncol_global=nao, &
                               para_env=para_env, context=blacs_env)
      CALL cp_fm_struct_create(ao_mo_occ_fm_struct, nrow_global=nao, ncol_global=nmo_occ_actual, &
                               para_env=para_env, context=blacs_env)
      CALL cp_fm_struct_create(ao_mo_virt_fm_struct, nrow_global=nao, ncol_global=nmo_virt_actual, &
                               para_env=para_env, context=blacs_env)

      ! storage space
      CALL cp_fm_create(mos_occ, ao_mo_occ_fm_struct)
      CALL cp_fm_create(mos_virt, ao_mo_virt_fm_struct)
      ALLOCATE (evals_occ(nmo_occ_actual))
      ALLOCATE (evals_virt(nmo_virt_actual))

      ! set of molecular orbitals
      NULLIFY (ao_ao_fm_pool, mos_extended)
      CALL fm_pool_create(ao_ao_fm_pool, ao_ao_fm_struct)
      CALL allocate_mo_set(mos_extended, nao, nao, nelectrons, &
                           REAL(nelectrons, dp), maxocc, flexible_electron_count=0.0_dp)
      CALL init_mo_set(mos_extended, fm_pool=ao_ao_fm_pool, name="mos-extended")
      CALL fm_pool_release(ao_ao_fm_pool)

      CALL get_mo_set(mos_extended, mo_coeff=mo_coeff_extended, &
                      eigenvalues=mo_evals_extended, occupation_numbers=mo_occ_extended)

      ! use the explicit loop in order to avoid temporary arrays.
      !
      ! The assignment statement : mo_occ_extended(1:nmo_scf) = mo_occ_scf(1:nmo_scf)
      ! implies temporary arrays as a compiler does not know in advance that the pointers
      ! on both sides of the statement point to non-overlapped memory regions
      DO iocc = 1, nmo_scf
         mo_occ_extended(iocc) = mo_occ_scf(iocc)
      END DO
      mo_occ_extended(nmo_scf+1:) = 0.0_dp

      ! ++ create temporary matrices
      NULLIFY (matrix_ks_fm, ortho_fm, work_fm)
      CALL cp_fm_create(matrix_ks_fm, ao_ao_fm_struct)
      CALL cp_fm_create(ortho_fm, ao_ao_fm_struct)
      CALL cp_fm_create(work_fm, ao_ao_fm_struct)

      CALL cp_fm_struct_release(ao_mo_virt_fm_struct)
      CALL cp_fm_struct_release(ao_mo_occ_fm_struct)
      CALL cp_fm_struct_release(ao_ao_fm_struct)

      ! some stuff from the subroutine general_eigenproblem()
      CALL copy_dbcsr_to_fm(matrix_s, ortho_fm)
      CALL copy_dbcsr_to_fm(matrix_ks, matrix_ks_fm)

      IF (cholesky_method == cholesky_dbcsr) THEN
         CPABORT('CHOLESKY DBCSR_INVERSE is not implemented in TDDFT.')
      ELSE IF (cholesky_method == cholesky_off) THEN
         CPABORT('CHOLESKY OFF is not implemented in TDDFT.')
      ELSE
         CALL cp_fm_cholesky_decompose(ortho_fm)
         IF (cholesky_method == cholesky_inverse) THEN
            CALL cp_fm_triangular_invert(ortho_fm)
         END IF

         ! need this temporary variable, as the subroutine eigensolver() is going to update it.
         cholesky_method_inout = cholesky_method
         CALL eigensolver(matrix_ks_fm=matrix_ks_fm, mo_set=mos_extended, ortho=ortho_fm, &
                          work=work_fm, cholesky_method=cholesky_method_inout, &
                          do_level_shift=.FALSE., level_shift=0.0_dp, matrix_u_fm=null(), use_jacobi=.FALSE.)
      END IF

      ! -- clean up needless matrices
      CALL cp_fm_release(ortho_fm)
      CALL cp_fm_release(work_fm)
      CALL cp_fm_release(matrix_ks_fm)

      ! return the requested occupied and virtual molecular orbitals and corresponding orbital energies
      CALL cp_fm_to_fm_submat(msource=mo_coeff_extended, mtarget=mos_occ, nrow=nao, ncol=nmo_occ_actual, &
                              s_firstrow=1, s_firstcol=homo-nmo_occ_actual+1, t_firstrow=1, t_firstcol=1)
      evals_occ(1:nmo_occ_actual) = mo_evals_extended(homo-nmo_occ_actual+1:homo)

      CALL cp_fm_to_fm_submat(msource=mo_coeff_extended, mtarget=mos_virt, nrow=nao, ncol=nmo_virt_actual, &
                              s_firstrow=1, s_firstcol=homo+1, t_firstrow=1, t_firstcol=1)
      evals_virt(1:nmo_virt_actual) = mo_evals_extended(homo+1:homo+nmo_virt_actual)

      CALL deallocate_mo_set(mos_extended)
      CALL timestop(handle)
   END SUBROUTINE tddfpt_lumos

! **************************************************************************************************
!> \brief Generate initial guess vectors for TDDFT and restrict the number of
!>        occupied orbitals in the active space.
!> \param nstates_active the number of excited states in the active space
!> \param evecs         trial vectors, on exit SIZE(evecs) >= nstates_active (allocated on exit)
!> \param evals         initial excitation energies; SIZE(evals) == the number of excited states
!>                      to converge (initialised on exit)
!> \param mos_occ       occupied molecular orbitals (reallocated on exit)
!> \param evals_occ     occupied orbital energies (reallocated on exit)
!> \param added_occ_mos the number of 'inactive' molecular orbitals to be added in the active space
!> \param mos_virt      virtual molecular orbitals
!> \param evals_virt    virtual orbital energies
!> \param threshold     energy threshold which controls when excited states
!>                      are considered to be degenerate
!> \param add_degenerate should we automatically include all degenerate excited states
!>                      into the active space
!> \param blacs_env     BLACS environment for newly allocated distributed matrices
!> \param para_env      parallel environment for newly allocated distributed matrices
!> \param log_unit      output unit
!> \par History
!>    * 05.2016 created [Sergey Chulkov]
!> \note Based on the subroutine co_initial_guess() which was originally created by
!>       Thomas Chassaing on 06.2003.
! **************************************************************************************************
   SUBROUTINE tddfpt_guess(nstates_active, evecs, evals, mos_occ, evals_occ, added_occ_mos, &
                           mos_virt, evals_virt, threshold, add_degenerate, blacs_env, para_env, log_unit)
      INTEGER, INTENT(in)                                :: nstates_active
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:), &
         INTENT(out)                                     :: evecs
      REAL(kind=dp), DIMENSION(:), INTENT(out)           :: evals
      TYPE(cp_fm_type), POINTER                          :: mos_occ
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(inout)                                   :: evals_occ
      INTEGER, INTENT(in)                                :: added_occ_mos
      TYPE(cp_fm_type), POINTER                          :: mos_virt
      REAL(kind=dp), DIMENSION(:), INTENT(in)            :: evals_virt
      REAL(kind=dp), INTENT(in)                          :: threshold
      LOGICAL, INTENT(in)                                :: add_degenerate
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(in)                                :: log_unit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_guess', routineP = moduleN//':'//routineN

      INTEGER :: handle, imo_occ, imo_virt, istate, nao, nmo_occ_active, nmo_occ_avail, &
         nmo_occ_selected, nmo_virt_avail, nmo_virt_selected, nstates, nstates_selected, &
         nstates_total
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: inds
      REAL(kind=dp)                                      :: e_occ
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: e_virt_minus_occ
      TYPE(cp_fm_struct_type), POINTER                   :: ao_mo_occ_fm_struct
      TYPE(cp_fm_type), POINTER                          :: mos_occ_new

      CALL timeset(routineN, handle)

      nstates = SIZE(evals)
      nstates_total = nstates_active
      nmo_occ_avail = SIZE(evals_occ)
      CALL cp_fm_get_info(mos_virt, nrow_global=nao, ncol_global=nmo_virt_avail)

      DO nmo_occ_selected = MIN(nstates_total, nmo_occ_avail), nmo_occ_avail-1
         IF (evals_occ(nmo_occ_avail-nmo_occ_selected+1)-evals_occ(nmo_occ_avail-nmo_occ_selected) > threshold) EXIT
      END DO

      DO nmo_virt_selected = MIN(nstates_total, nmo_virt_avail), nmo_virt_avail-1
         IF (evals_virt(nmo_virt_selected+1)-evals_virt(nmo_virt_selected) > threshold) EXIT
      END DO

      ! TO DO: the variable 'nstates_selected' should probably be declared as INTEGER(kind=int_8),
      !        however we need a special version of the subroutine sort() in order to do so
      nstates_selected = nmo_occ_selected*nmo_virt_selected

      ALLOCATE (inds(nstates_selected))
      ALLOCATE (e_virt_minus_occ(nstates_selected))

      istate = 0
      DO imo_occ = 1, nmo_occ_selected
         e_occ = evals_occ(nmo_occ_avail-imo_occ+1)
         DO imo_virt = 1, nmo_virt_selected
            istate = istate+1
            e_virt_minus_occ(istate) = evals_virt(imo_virt)-e_occ
         END DO
      END DO

      CALL sort(e_virt_minus_occ, nstates_selected, inds)

      nmo_occ_active = 0
      IF (log_unit > 0) THEN
         WRITE (log_unit, '(/,18X,A)') "TDDFT initial guess"
         WRITE (log_unit, '(1X,52("-"))')
         WRITE (log_unit, '(5X,A)') "State     Occupied  ->  Virtual      Excitation"
         WRITE (log_unit, '(5X,A)') "number     orbital      orbital      energy (eV)"
         WRITE (log_unit, '(1X,52("-"))')
      END IF
      DO istate = 1, nstates_total
         IF (log_unit > 0 .AND. istate == nstates+1) &
            WRITE (log_unit, '(1X,15("-"),1X,A,1X,16("-"))') "Other active states"

         imo_occ = (inds(istate)-1)/nmo_virt_selected+1
         imo_virt = MOD(inds(istate)-1, nmo_virt_selected)+1
         IF (log_unit > 0) &
            WRITE (log_unit, '(1X,I8,4X,I8,5X,I8,5X,F14.5)') istate, nmo_occ_avail-imo_occ+1, &
            imo_virt+nmo_occ_avail, e_virt_minus_occ(istate)*evolt

         IF (nmo_occ_active < imo_occ) &
            nmo_occ_active = imo_occ

         IF (istate <= nstates) &
            evals(istate) = e_virt_minus_occ(istate)
      END DO

      DO istate = nstates_total+1, nstates_selected
         IF (e_virt_minus_occ(istate)-e_virt_minus_occ(istate-1) > threshold) EXIT

         IF (log_unit > 0 .AND. istate == nstates_total+1) THEN
            IF (add_degenerate) THEN
               IF (nstates == nstates_total) &
                  WRITE (log_unit, '(1X,15("-"),1X,A,1X,16("-"))') "Other active states"
            ELSE
               WRITE (log_unit, '(1X,17("-"),1X,A,1X,18("-"))') "Inactive states"
            END IF
         END IF

         imo_occ = (inds(istate)-1)/nmo_virt_selected+1
         imo_virt = MOD(inds(istate)-1, nmo_virt_selected)+1

         IF (add_degenerate .AND. nmo_occ_active < imo_occ) &
            nmo_occ_active = imo_occ

         IF (log_unit > 0) &
            WRITE (log_unit, '(1X,I8,4X,I8,4X,I8,6X,F14.5)') istate, nmo_occ_avail-imo_occ+1, &
            imo_virt+nmo_occ_avail, e_virt_minus_occ(istate)*evolt
      END DO
      DEALLOCATE (e_virt_minus_occ)

      ! limit the number of occupied molecular orbitals in the active space and create guess vectors
      CPASSERT(nmo_occ_active > 0)

      IF (added_occ_mos < 0 .OR. nmo_occ_avail < nmo_occ_active+added_occ_mos) THEN
         nmo_occ_active = nmo_occ_avail
      ELSE
         nmo_occ_active = nmo_occ_active+added_occ_mos
      END IF

      IF (log_unit > 0) THEN
         WRITE (log_unit, '(/,1X,A,T72,I8)') 'Number of occupied molecular orbitals in the active space:', nmo_occ_active
      END IF

      IF (istate > nstates_total+1) THEN
         IF (add_degenerate) THEN
            nstates_total = istate-1
         ELSE
            CALL cp_warn(__LOCATION__, "It is advised to include all proposed"// &
                         " degenerate excited states into TDDFT active space")
         END IF
      END IF

      NULLIFY (ao_mo_occ_fm_struct)
      CALL cp_fm_struct_create(ao_mo_occ_fm_struct, nrow_global=nao, ncol_global=nmo_occ_active, &
                               para_env=para_env, context=blacs_env)

      IF (nmo_occ_active < nmo_occ_avail) THEN
         NULLIFY (mos_occ_new)
         CALL cp_fm_create(mos_occ_new, ao_mo_occ_fm_struct)
         CALL cp_fm_to_fm_submat(msource=mos_occ, mtarget=mos_occ_new, nrow=nao, ncol=nmo_occ_active, &
                                 s_firstrow=1, s_firstcol=nmo_occ_avail-nmo_occ_active+1, t_firstrow=1, t_firstcol=1)
         CALL cp_fm_release(mos_occ)
         mos_occ => mos_occ_new

         ALLOCATE (e_virt_minus_occ(nmo_occ_active))
         e_virt_minus_occ(1:nmo_occ_active) = evals_occ(nmo_occ_avail-nmo_occ_active+1:nmo_occ_avail)

         DEALLOCATE (evals_occ)
         ALLOCATE (evals_occ(nmo_occ_active))
         evals_occ(:) = e_virt_minus_occ(:)

         DEALLOCATE (e_virt_minus_occ)
      END IF

      ALLOCATE (evecs(nstates_total))
      DO istate = 1, nstates_total
         NULLIFY (evecs(istate)%matrix)
         CALL cp_fm_create(evecs(istate)%matrix, ao_mo_occ_fm_struct)
         CALL cp_fm_set_all(evecs(istate)%matrix, 0.0_dp)

         imo_occ = (inds(istate)-1)/nmo_virt_selected+1
         imo_virt = MOD(inds(istate)-1, nmo_virt_selected)+1

         CALL cp_fm_to_fm_submat(msource=mos_virt, mtarget=evecs(istate)%matrix, nrow=nao, ncol=1, &
                                 s_firstrow=1, s_firstcol=imo_virt, t_firstrow=1, t_firstcol=nmo_occ_active-imo_occ+1)
      END DO

      CALL cp_fm_struct_release(ao_mo_occ_fm_struct)
      DEALLOCATE (inds)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_guess

! **************************************************************************************************
!> \brief Initialise the action_matrix using the following expression:
!>        trace { C_{1,i}^T * ( KS * C_{1,j} - S * C_{1,j} * diag(evals_occ) ) } .
!> \param a_kv      action_matrix (initialised on exit)
!> \param evecs     TDDFPT trial vectors
!> \param evals_occ occupied orbital energies
!> \param matrix_ks Kohn-Sham matrix for a given spin
!> \param matrix_s  overlap matrix
!> \par History
!>    * 05.2016 initialise all matrix elements in one go [Sergey Chulkov]
!> \note Based on the subroutine p_op_l1() which was originally created by
!>       Thomas Chassaing on 08.2002.
! **************************************************************************************************
   SUBROUTINE tddfpt_init_by_energy_diff(a_kv, evecs, evals_occ, matrix_ks, matrix_s)
      REAL(kind=dp), DIMENSION(:, :), INTENT(out)        :: a_kv
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(in)       :: evecs
      REAL(kind=dp), DIMENSION(:), INTENT(in)            :: evals_occ
      TYPE(cp_dbcsr_type), POINTER                       :: matrix_ks, matrix_s

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_init_by_energy_diff', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, istate, jstate, nmo, nstates
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: ks_psi1_fm, s_psi1_fm

      CALL timeset(routineN, handle)

      nstates = SIZE(evecs)

      NULLIFY (fm_struct, s_psi1_fm)
      CALL cp_fm_get_info(matrix=evecs(1)%matrix, ncol_global=nmo, matrix_struct=fm_struct)
      CALL cp_fm_create(ks_psi1_fm, fm_struct)
      CALL cp_fm_create(s_psi1_fm, fm_struct)

      DO jstate = 1, nstates
         CALL cp_dbcsr_sm_fm_multiply(matrix_ks, evecs(jstate)%matrix, ks_psi1_fm, ncol=nmo, alpha=1.0_dp, beta=0.0_dp)
         CALL cp_dbcsr_sm_fm_multiply(matrix_s, evecs(jstate)%matrix, s_psi1_fm, ncol=nmo, alpha=1.0_dp, beta=0.0_dp)
         CALL cp_fm_column_scale(s_psi1_fm, evals_occ)
         !  KS * C_1 - S * C_1 * occupied_orbital_energies
         CALL cp_fm_scale_and_add(1.0_dp, ks_psi1_fm, -1.0_dp, s_psi1_fm)

         DO istate = 1, jstate
            CALL cp_fm_trace(evecs(istate)%matrix, ks_psi1_fm, a_kv(istate, jstate))
         END DO
      END DO

      ! on exit : a_kv(i, j) = \sum_{occ, virt} (C_virt * epsilon_virt - C_occ * epsilon_occ)
      CALL cp_fm_release(s_psi1_fm)
      CALL cp_fm_release(ks_psi1_fm)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_init_by_energy_diff

! **************************************************************************************************
!> \brief  Update the action matrix by adding coulomb terms
!>         a_kv(a,b) += alpha * trace { C_{1,a}^T * sum_j {(mu nu || j b) * C_{0,j} } } .
!> \param a_kv     action matrix (updated on exit)
!> \param alpha    scale factor
!> \param evecs    TDDFPT trial vectors
!> \param mos_occ  occupied molecular orbitals
!> \param qs_env   Quickstep environment
!> \param matrix_s overlap matrix
!> \par History
!>    * 05.2016 compute all coulomb terms in one go [Sergey Chulkov]
!> \note Based on the subroutine kpp1_calc_k_p_p1() which was originally created by
!>       Mohamed Fawzi on 10.2002.
! **************************************************************************************************
   SUBROUTINE tddfpt_apply_coulomb(a_kv, alpha, evecs, mos_occ, qs_env, matrix_s)
      REAL(kind=dp), DIMENSION(:, :), INTENT(out)        :: a_kv
      REAL(kind=dp), INTENT(in)                          :: alpha
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(in)       :: evecs
      TYPE(cp_fm_type), POINTER                          :: mos_occ
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_dbcsr_type), POINTER                       :: matrix_s

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_apply_coulomb', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, istate, jstate, nao, nmo, nstates
      REAL(kind=dp)                                      :: J_ia_jb
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_dbcsr_p_type)                              :: J_ia_munu
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: rho_ia_ao
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: J_ia_jnu, rho_ia_fm
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: v_gspace, v_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_ia_g
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_rho_type), POINTER                         :: rho_ia_struct

      CALL timeset(routineN, handle)

      NULLIFY (blacs_env, para_env, pw_env)
      CALL get_qs_env(qs_env, blacs_env=blacs_env, para_env=para_env, pw_env=pw_env)

      NULLIFY (auxbas_pw_pool, poisson_env)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, poisson_env=poisson_env)

      nstates = SIZE(evecs)
      CALL cp_fm_get_info(matrix=evecs(1)%matrix, nrow_global=nao, ncol_global=nmo)

      NULLIFY (rho_ia_ao)
      CALL cp_dbcsr_allocate_matrix_set(rho_ia_ao, 1)
      CALL cp_dbcsr_init_p(rho_ia_ao(1)%matrix)
      CALL cp_dbcsr_copy(rho_ia_ao(1)%matrix, matrix_s)

      NULLIFY (rho_ia_struct, rho_ia_g)
      CALL qs_rho_create(rho_ia_struct)
      CALL qs_rho_set(rho_ia_struct, rho_ao=rho_ia_ao)

      CALL qs_rho_rebuild(rho_ia_struct, qs_env, rebuild_ao=.FALSE., rebuild_grids=.TRUE.)
      CALL qs_rho_get(rho_ia_struct, rho_g=rho_ia_g)

      NULLIFY (J_ia_munu%matrix)
      CALL cp_dbcsr_init_p(J_ia_munu%matrix)
      ! TO DO: is there a better way to initialise a DBCSR matrix?
      ! The subroutine cp_dbcsr_create() alone does not work here,
      ! as block distribution of the matrix J_ia_munu should be identical with the matrix_s
      ! (the subroutine integrate_v_rspace() will fail otherwise) .
      CALL cp_dbcsr_copy(J_ia_munu%matrix, matrix_s)

      NULLIFY (fm_struct, rho_ia_fm)
      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, para_env=para_env, context=blacs_env)
      CALL cp_fm_create(rho_ia_fm, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      NULLIFY (fm_struct, J_ia_jnu)
      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nmo, para_env=para_env, context=blacs_env)
      CALL cp_fm_create(J_ia_jnu, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      NULLIFY (v_gspace%pw, v_rspace%pw)
      CALL pw_pool_create_pw(auxbas_pw_pool, v_gspace%pw, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, v_rspace%pw, use_data=REALDATA3D, in_space=REALSPACE)

      DO jstate = 1, nstates
         CALL cp_gemm('N', 'T', nao, nao, nmo, 0.5_dp, mos_occ, evecs(jstate)%matrix, 0.0_dp, rho_ia_fm)
         CALL cp_gemm('N', 'T', nao, nao, nmo, 0.5_dp, evecs(jstate)%matrix, mos_occ, 1.0_dp, rho_ia_fm)

         CALL copy_fm_to_dbcsr(rho_ia_fm, rho_ia_ao(1)%matrix, keep_sparsity=.TRUE.)
         CALL qs_rho_update_rho(rho_ia_struct, qs_env)

         ! here J_ia_jb is used as a temporary variable to store needless pair energy
         CALL pw_poisson_solve(poisson_env, rho_ia_g(1)%pw, J_ia_jb, v_gspace%pw)
         CALL pw_transfer(v_gspace%pw, v_rspace%pw)
         CALL pw_scale(v_rspace%pw, v_rspace%pw%pw_grid%dvol)

         CALL cp_dbcsr_set(J_ia_munu%matrix, 0.0_dp)
         CALL integrate_v_rspace(v_rspace=v_rspace, hmat=J_ia_munu, qs_env=qs_env, calculate_forces=.FALSE., gapw=.FALSE.)

         CALL cp_dbcsr_sm_fm_multiply(J_ia_munu%matrix, mos_occ, J_ia_jnu, ncol=nmo, alpha=1.0_dp, beta=0.0_dp)
         DO istate = 1, jstate
            CALL cp_fm_trace(evecs(istate)%matrix, J_ia_jnu, J_ia_jb)
            a_kv(istate, jstate) = a_kv(istate, jstate)+alpha*J_ia_jb
         END DO
      END DO

      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_rspace%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_gspace%pw)

      CALL cp_fm_release(J_ia_jnu)
      CALL cp_fm_release(rho_ia_fm)

      CALL cp_dbcsr_deallocate_matrix(J_ia_munu%matrix)
      CALL qs_rho_release(rho_ia_struct)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_apply_coulomb

! **************************************************************************************************
!> \brief Update the action matrix by adding kernel terms.
!> \param a_kv         action matrix (updated on exit)
!> \param alpha        scale factor
!> \param evecs        TDDFPT trial vectors
!> \param mos_occ      occupied molecular orbitals
!> \param xc_rho_set   some variable to compute the 2nd derivative of XC-functional
!> \param xc_deriv_set some variable to compute the 2nd derivative of XC-functional
!> \param qs_env       Quickstep environment
!> \param matrix_s     overlap matrix
!> \param xc_section   XC input section
!> \par History
!>    * 05.2016 compute all kernel terms in one go [Sergey Chulkov]
!> \note Based on the subroutine kpp1_calc_k_p_p1() which was originally created by
!>       Mohamed Fawzi on 10.2002.
! **************************************************************************************************
   SUBROUTINE tddfpt_apply_kernel(a_kv, alpha, evecs, mos_occ, xc_rho_set, xc_deriv_set, qs_env, matrix_s, xc_section)
      REAL(kind=dp), DIMENSION(:, :), INTENT(out)        :: a_kv
      REAL(kind=dp), INTENT(in)                          :: alpha
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(in)       :: evecs
      TYPE(cp_fm_type), POINTER                          :: mos_occ
      TYPE(xc_rho_set_type), POINTER                     :: xc_rho_set
      TYPE(xc_derivative_set_type), POINTER              :: xc_deriv_set
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_dbcsr_type), POINTER                       :: matrix_s
      TYPE(section_vals_type), POINTER                   :: xc_section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_apply_kernel', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, istate, jstate, nao, nmo, &
                                                            nstates, xc_deriv_method_id, &
                                                            xc_rho_smooth_id
      REAL(kind=dp)                                      :: F_ia_jb
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_dbcsr_p_type)                              :: F_ia_munu
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: rho_ia_ao
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: F_ia_jnu, rho_ia_fm
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_ia_g, rho_ia_r, v_xc
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_rho_type), POINTER                         :: rho_ia_struct
      TYPE(section_vals_type), POINTER                   :: xc_fun_section
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type), POINTER                     :: xc_rho1_set

      CALL timeset(routineN, handle)

      NULLIFY (blacs_env, para_env, pw_env)
      CALL get_qs_env(qs_env, blacs_env=blacs_env, para_env=para_env, pw_env=pw_env)

      NULLIFY (auxbas_pw_pool)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      nstates = SIZE(evecs)
      CALL cp_fm_get_info(matrix=evecs(1)%matrix, nrow_global=nao, ncol_global=nmo)

      NULLIFY (rho_ia_ao)
      CALL cp_dbcsr_allocate_matrix_set(rho_ia_ao, 1)
      CALL cp_dbcsr_init_p(rho_ia_ao(1)%matrix)
      CALL cp_dbcsr_copy(rho_ia_ao(1)%matrix, matrix_s)

      NULLIFY (rho_ia_struct, rho_ia_r, rho_ia_g)
      CALL qs_rho_create(rho_ia_struct)
      CALL qs_rho_set(rho_ia_struct, rho_ao=rho_ia_ao)

      CALL qs_rho_rebuild(rho_ia_struct, qs_env, rebuild_ao=.FALSE., rebuild_grids=.TRUE.)
      CALL qs_rho_get(rho_ia_struct, rho_r=rho_ia_r, rho_g=rho_ia_g)

      NULLIFY (F_ia_munu%matrix)
      CALL cp_dbcsr_init_p(F_ia_munu%matrix)
      CALL cp_dbcsr_copy(F_ia_munu%matrix, matrix_s)

      NULLIFY (fm_struct, rho_ia_fm)
      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, para_env=para_env, context=blacs_env)
      CALL cp_fm_create(rho_ia_fm, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      NULLIFY (fm_struct, F_ia_jnu)
      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nmo, para_env=para_env, context=blacs_env)
      CALL cp_fm_create(F_ia_jnu, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      ALLOCATE (v_xc(1))
      NULLIFY (v_xc(1)%pw)
      CALL pw_pool_create_pw(auxbas_pw_pool, v_xc(1)%pw, use_data=REALDATA3D, in_space=REALSPACE)

      NULLIFY (xc_rho1_set)
      CALL xc_rho_set_create(xc_rho1_set, rho_ia_r(1)%pw%pw_grid%bounds_local, &
                             rho_cutoff=section_get_rval(xc_section, "DENSITY_CUTOFF"), &
                             drho_cutoff=section_get_rval(xc_section, "GRADIENT_CUTOFF"), &
                             tau_cutoff=section_get_rval(xc_section, "TAU_CUTOFF"))

      xc_fun_section => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      xc_deriv_method_id = section_get_ival(xc_section, "XC_GRID%XC_DERIV")
      xc_rho_smooth_id = section_get_ival(xc_section, "XC_GRID%XC_SMOOTH_RHO")

      ! get arguments needed by the given functional
      needs = xc_functionals_get_needs(functionals=xc_fun_section, lsd=.FALSE., add_basic_components=.TRUE.)

      DO jstate = 1, nstates
         CALL cp_gemm('N', 'T', nao, nao, nmo, 0.5_dp, mos_occ, evecs(jstate)%matrix, 0.0_dp, rho_ia_fm)
         CALL cp_gemm('N', 'T', nao, nao, nmo, 0.5_dp, evecs(jstate)%matrix, mos_occ, 1.0_dp, rho_ia_fm)

         CALL copy_fm_to_dbcsr(rho_ia_fm, rho_ia_ao(1)%matrix, keep_sparsity=.TRUE.)
         CALL qs_rho_update_rho(rho_ia_struct, qs_env)

         CALL xc_rho_set_update(rho_set=xc_rho1_set, rho_r=rho_ia_r, rho_g=rho_ia_g, &
                                tau=null(), needs=needs, xc_deriv_method_id=xc_deriv_method_id, &
                                xc_rho_smooth_id=xc_rho_smooth_id, pw_pool=auxbas_pw_pool)

         CALL pw_zero(v_xc(1)%pw)
         CALL xc_calc_2nd_deriv(v_xc=v_xc, deriv_set=xc_deriv_set, rho_set=xc_rho_set, &
                                rho1_set=xc_rho1_set, pw_pool=auxbas_pw_pool, &
                                xc_section=xc_section, gapw=.FALSE., tddfpt_fac=0.0_dp)

         CALL pw_scale(v_xc(1)%pw, v_xc(1)%pw%pw_grid%dvol)

         CALL cp_dbcsr_set(F_ia_munu%matrix, 0.0_dp)
         CALL integrate_v_rspace(v_rspace=v_xc(1), hmat=F_ia_munu, qs_env=qs_env, calculate_forces=.FALSE., gapw=.FALSE.)

         CALL cp_dbcsr_sm_fm_multiply(F_ia_munu%matrix, mos_occ, F_ia_jnu, ncol=nmo, alpha=1.0_dp, beta=0.0_dp)

         DO istate = 1, jstate
            CALL cp_fm_trace(evecs(istate)%matrix, F_ia_jnu, F_ia_jb)
            a_kv(istate, jstate) = a_kv(istate, jstate)+alpha*F_ia_jb
         END DO
      END DO

      CALL xc_rho_set_release(xc_rho1_set)

      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_xc(1)%pw)
      DEALLOCATE (v_xc)

      CALL cp_fm_release(F_ia_jnu)
      CALL cp_fm_release(rho_ia_fm)

      CALL cp_dbcsr_deallocate_matrix(F_ia_munu%matrix)
      CALL qs_rho_release(rho_ia_struct)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_apply_kernel

! **************************************************************************************************
!> \brief Compute half-transformed four-centre integrals ( i j || mu nu) where indices
!>        i,j and mu,nu denote occupied molecular orbitals and atomic basis functions respectively.
!> \param K_ij_munu  a vector to store the integrals
!> \param mos_occ    occupied molecular orbitals
!> \param qs_env     Quickstep environment
!> \param matrix_s   overlap matrix
!> \par History
!>    * 05.2016 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_compute_hfx(K_ij_munu, mos_occ, qs_env, matrix_s)
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(in)       :: K_ij_munu
      TYPE(cp_fm_type), POINTER                          :: mos_occ
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_dbcsr_type), POINTER                       :: matrix_s

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_compute_hfx', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i_i, i_j, nmo
      REAL(kind=dp)                                      :: pair_energy
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_dbcsr_p_type)                              :: hfx_munu
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: psi0_i, psi0_j, rho_g, v_gspace
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      NULLIFY (atomic_kind_set, cell, dft_control, particle_set, pw_env, qs_kind_set)
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, cell=cell, &
                      dft_control=dft_control, particle_set=particle_set, &
                      pw_env=pw_env, qs_kind_set=qs_kind_set)

      NULLIFY (auxbas_pw_pool, poisson_env)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, poisson_env=poisson_env)

      CALL cp_fm_get_info(matrix=mos_occ, ncol_global=nmo)

      NULLIFY (psi0_i%pw, psi0_j%pw)
      CALL pw_pool_create_pw(auxbas_pw_pool, psi0_i%pw, use_data=REALDATA3D, in_space=REALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, psi0_j%pw, use_data=REALDATA3D, in_space=REALSPACE)

      NULLIFY (v_gspace%pw, rho_g%pw)
      CALL pw_pool_create_pw(auxbas_pw_pool, v_gspace%pw, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, rho_g%pw, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)

      NULLIFY (hfx_munu%matrix)
      CALL cp_dbcsr_init_p(hfx_munu%matrix)
      CALL cp_dbcsr_copy(hfx_munu%matrix, matrix_s)

      DO i_j = 1, nmo
         CALL calculate_wavefunction(mos_occ, i_j, psi0_j, rho_g, atomic_kind_set, &
                                     qs_kind_set, cell, dft_control, particle_set, pw_env)

         DO i_i = 1, i_j
            CALL calculate_wavefunction(mos_occ, i_i, psi0_i, rho_g, atomic_kind_set, &
                                        qs_kind_set, cell, dft_control, particle_set, pw_env)

            psi0_i%pw%cr3d = psi0_i%pw%cr3d*psi0_j%pw%cr3d

            CALL pw_transfer(psi0_i%pw, rho_g%pw)
            CALL pw_poisson_solve(poisson_env, rho_g%pw, pair_energy, v_gspace%pw)
            CALL pw_transfer(v_gspace%pw, psi0_i%pw)
            CALL pw_scale(psi0_i%pw, psi0_i%pw%pw_grid%dvol)

            CALL cp_dbcsr_set(hfx_munu%matrix, 0.0_dp)
            CALL integrate_v_rspace(v_rspace=psi0_i, hmat=hfx_munu, qs_env=qs_env, calculate_forces=.FALSE., gapw=.FALSE.)

            ! compute the index of half-transformed integrals on the fly using the arithmetic series formula,
            ! so the uppermost DO loop can be executed in parallel
            CALL copy_dbcsr_to_fm(hfx_munu%matrix, K_ij_munu((i_j-1)*i_j/2+i_i)%matrix)
         END DO
      END DO

      CALL cp_dbcsr_deallocate_matrix(hfx_munu%matrix)

      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_g%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_gspace%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, psi0_j%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, psi0_i%pw)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_compute_hfx

! **************************************************************************************************
!> \brief  Update the action matrix by adding exact-exchange terms
!>         a_kv(a, b) += alpha * sum_{i<=j} trace { C_{1,i->a}^T * (i j || mu nu) * C_{1,j->b} } .
!> \param a_kv      action matrix (updated on exit)
!> \param K_ij_munu half-transformed four-centre integrals
!> \param evecs     TDDFPT trial vectors
!> \param alpha     scale factor
!> \param qs_env    Quickstep environment
!> \par History
!>    * 05.2016 compute all exact-exchange terms in one go [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_apply_hfx(a_kv, K_ij_munu, evecs, alpha, qs_env)
      REAL(kind=dp), DIMENSION(:, :), INTENT(out)        :: a_kv
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(in)       :: K_ij_munu, evecs
      REAL(kind=dp), INTENT(in)                          :: alpha
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_apply_hfx', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i_a, i_b, i_i, i_j, ind, nao, &
                                                            nmo, nstates
      REAL(kind=dp)                                      :: K_ij_ab
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: hfx_psi0, mo_coeff_a, mo_coeff_b
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      NULLIFY (blacs_env, para_env)
      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env)

      CALL cp_fm_get_info(matrix=evecs(1)%matrix, nrow_global=nao, ncol_global=nmo)
      nstates = SIZE(evecs)

      NULLIFY (fm_struct, hfx_psi0)
      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nmo, para_env=para_env, context=blacs_env)
      CALL cp_fm_create(hfx_psi0, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      NULLIFY (fm_struct, mo_coeff_a, mo_coeff_b)
      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=1, para_env=para_env, context=blacs_env)
      CALL cp_fm_create(mo_coeff_a, fm_struct)
      CALL cp_fm_create(mo_coeff_b, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      ! TO DO: consider to combine the next two loops into one loop
      DO i_j = 1, nmo
         DO i_i = 1, i_j
            ! compute the index on the fly using the arithmetic series formula,
            ! so the uppermost DO loop can be executed in parallel
            ind = (i_j-1)*i_j/2+i_i

            DO i_a = 1, nstates
               CALL cp_fm_symm('L', 'U', nao, nmo, 1.0_dp, K_ij_munu(ind)%matrix, evecs(i_a)%matrix, 0.0_dp, hfx_psi0)
               CALL cp_fm_to_fm_submat(hfx_psi0, mo_coeff_a, nao, 1, 1, i_i, 1, 1)

               DO i_b = 1, i_a
                  CALL cp_fm_to_fm_submat(evecs(i_b)%matrix, mo_coeff_b, nao, 1, 1, i_j, 1, 1)
                  CALL cp_fm_trace(mo_coeff_b, mo_coeff_a, K_ij_ab)

                  a_kv(i_b, i_a) = a_kv(i_b, i_a)-alpha*K_ij_ab
               END DO
            END DO
         END DO
      END DO

      CALL cp_fm_release(mo_coeff_b)
      CALL cp_fm_release(mo_coeff_a)
      CALL cp_fm_release(hfx_psi0)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_apply_hfx

! **************************************************************************************************
!> \brief Make a TDDFPT trial vector orthogonal to all occupied molecular orbitals.
!> \param C1       trial vector (modified on exit)
!> \param C0_C0T_S matrix product C_0 * C_0^T * S, where C_0 is the ground state wave-function
!>                     expressed in atomic basis set, and S is the corresponding overlap matrix
!> \par History
!>    * 05.2016 created [Sergey Chulkov]
!> \note  Based on the subroutine p_preortho() which was created by Thomas Chassaing on 09.2002.
!>        Seems to be useless at the moment, as all virtual states are already orthogonal to
!>        occupied states by design.
! **************************************************************************************************
   SUBROUTINE tddfpt_orthogonalize_psi1_psi0(C1, C0_C0T_S)
      TYPE(cp_fm_type), POINTER                          :: C1, C0_C0T_S

      INTEGER                                            :: nao, nmo
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: C0_C0T_S_C1

      NULLIFY (fm_struct, C0_C0T_S_C1)
      CALL cp_fm_get_info(matrix=C1, nrow_global=nao, ncol_global=nmo, matrix_struct=fm_struct)
      CALL cp_fm_create(C0_C0T_S_C1, fm_struct)

      CALL cp_gemm('N', 'N', nao, nmo, nao, 1.0_dp, C0_C0T_S, C1, 0.0_dp, C0_C0T_S_C1)
      CALL cp_fm_scale_and_add(1.0_dp, C1, -1.0_dp, C0_C0T_S_C1)

      CALL cp_fm_release(C0_C0T_S_C1)
   END SUBROUTINE tddfpt_orthogonalize_psi1_psi0

! **************************************************************************************************
!> \brief Make a TDDFPT trial vector orthonormal to all previous trial vectors.
!> \param C1     trial vectors (the last vector is modified on exit)
!> \param matrix_s overlap matrix
!> \param n        number of trial vectors
!> \par History
!>    * 05.2016 created [Sergey Chulkov]
!> \note Based on the subroutines reorthogonalize() and normalize() which were originally created
!>             by Thomas Chassaing on 03.2003.
! **************************************************************************************************
   SUBROUTINE tddfpt_orthonormalize_psi1_psi1(C1, matrix_s, n)
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(in)       :: C1
      TYPE(cp_dbcsr_type), POINTER                       :: matrix_s
      INTEGER, INTENT(in)                                :: n

      INTEGER                                            :: i, nmo
      REAL(kind=dp)                                      :: trace
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: S_C1

      NULLIFY (fm_struct, S_C1)
      CALL cp_fm_get_info(matrix=C1(1)%matrix, ncol_global=nmo, matrix_struct=fm_struct)
      CALL cp_fm_create(S_C1, fm_struct)

      DO i = 1, n-1
         CALL cp_dbcsr_sm_fm_multiply(matrix_s, C1(i)%matrix, S_C1, ncol=nmo, alpha=1.0_dp, beta=0.0_dp)
         CALL cp_fm_trace(C1(n)%matrix, S_C1, trace)
         CALL cp_fm_scale_and_add(1.0_dp, C1(n)%matrix, -trace, C1(i)%matrix)
      END DO

      CALL cp_dbcsr_sm_fm_multiply(matrix_s, C1(n)%matrix, S_C1, ncol=nmo, alpha=1.0_dp, beta=0.0_dp)
      CALL cp_fm_trace(C1(n)%matrix, S_C1, trace)
      CALL cp_fm_scale((1.0_dp/SQRT(trace)), C1(n)%matrix)

      CALL cp_fm_release(S_C1)
   END SUBROUTINE tddfpt_orthonormalize_psi1_psi1

! **************************************************************************************************
!> \brief Print final TDDFPT excitation energies and analysis.
!> \param log_unit output unit
!> \param nstates  number of excited states to print
!> \param evecs    TDDFPT trial vectors
!> \param evals    TDDFPT eigenvalues
!> \param nmo_occ_total total number of occupied molecular orbitals
!> \param mos_occ occupied molecular orbitals
!> \param evals_occ  occupied orbital energies
!> \param mos_virt virtual molecular orbitals
!> \param evals_virt  virtual orbital energies
!> \param matrix_s overlap matrix and their firts derivatives over x, y, and z directions
!> \param min_amplitude the smallest excitation amplitude to print
!> \par History
!>    * 05.2016 created [Sergey Chulkov]
!>    * 06.2016 transition dipole moments and oscillator strengths [Sergey Chulkov]
!> \note \parblock
!>       Adapted version of the subroutine find_contributions() which was originally created
!>       by Thomas Chassaing on 02.2005.
!>
!>       To compute transition dipole integrals we use the fact that the ground state wave-function
!>       and all response functions (which are in fact a linear combination of virtual molecular
!>       orbitals) are eigenfunctions of the Hamiltonian operator. Using the identity:
!>       \f[\vec{r}\hat{H} - \hat{H}\vec{r} = [\vec{r},\hat{H}] = [\vec{r},-1/2 \nabla^2] = \nabla\f]
!>       we can derive a relationship between transition dipole integrals and momentum integrals:
!>       \f[<\psi_i|\nabla|\psi_a> = <\psi_i|\vec{r}|\hat{H}\psi_a> - <\hat{H}\psi_i|\vec{r}|\psi_a>=
!>       (\epsilon_a - \epsilon_i) <\psi_i|\vec{r}|\psi_a> .\f]
!>       \endparblock
! **************************************************************************************************
   SUBROUTINE tddfpt_print_summary(log_unit, nstates, evecs, evals, nmo_occ_total, &
                                   mos_occ, evals_occ, mos_virt, evals_virt, matrix_s, min_amplitude)
      INTEGER, INTENT(in)                                :: log_unit, nstates
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(in)       :: evecs
      REAL(kind=dp), DIMENSION(:), INTENT(in)            :: evals
      INTEGER, INTENT(in)                                :: nmo_occ_total
      TYPE(cp_fm_type), POINTER                          :: mos_occ
      REAL(kind=dp), DIMENSION(:), INTENT(in)            :: evals_occ
      TYPE(cp_fm_type), POINTER                          :: mos_virt
      REAL(kind=dp), DIMENSION(:), INTENT(in)            :: evals_virt
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s
      REAL(kind=dp), INTENT(in)                          :: min_amplitude

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_print_summary', &
         routineP = moduleN//':'//routineN
      INTEGER, PARAMETER                                 :: nderivs = 3

      INTEGER                                            :: handle, ideriv, iexc, imo_occ, imo_virt, &
                                                            istate, nao, nmo_occ, nmo_virt, &
                                                            nmo_virt_occ
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: inds
      REAL(kind=dp)                                      :: oscillator_strength
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: weights_abs_1d
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: weights_2d
      REAL(kind=dp), DIMENSION(nderivs)                  :: trans_dipole
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_p_type), DIMENSION(nderivs)             :: dS_mos_occ
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: ediff_inv, ediff_inv_weights, &
                                                            S_mos_virt, scaled_evect, weights_fm
      TYPE(cp_para_env_type), POINTER                    :: para_env

      EXTERNAL :: dcopy

      CALL timeset(routineN, handle)

      ! ensure we have the overlap matrix itself and its first derivatives over x, y and z
      CPASSERT(SIZE(matrix_s) >= nderivs+1)

      NULLIFY (fm_struct, blacs_env, para_env, S_mos_virt)
      CALL cp_fm_get_info(matrix=mos_virt, nrow_global=nao, ncol_global=nmo_virt, &
                          matrix_struct=fm_struct, context=blacs_env, para_env=para_env)

      CALL cp_fm_create(S_mos_virt, fm_struct)

      NULLIFY (fm_struct)
      CALL cp_fm_get_info(matrix=evecs(1)%matrix, ncol_global=nmo_occ, matrix_struct=fm_struct)
      DO ideriv = 1, nderivs
         NULLIFY (dS_mos_occ(ideriv)%matrix)
         CALL cp_fm_create(dS_mos_occ(ideriv)%matrix, fm_struct)
      END DO

      CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, mos_virt, S_mos_virt, ncol=nmo_virt, alpha=1.0_dp, beta=0.0_dp)
      DO ideriv = 1, nderivs
         CALL cp_dbcsr_sm_fm_multiply(matrix_s(ideriv+1)%matrix, mos_occ, dS_mos_occ(ideriv)%matrix, &
                                      ncol=nmo_occ, alpha=1.0_dp, beta=0.0_dp)
      END DO

      NULLIFY (fm_struct, weights_fm, ediff_inv, ediff_inv_weights)
      CALL cp_fm_struct_create(fm_struct, nrow_global=nmo_virt, ncol_global=nmo_occ, para_env=para_env, context=blacs_env)
      CALL cp_fm_create(weights_fm, fm_struct)
      CALL cp_fm_create(ediff_inv, fm_struct)
      CALL cp_fm_create(ediff_inv_weights, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      CALL cp_fm_get_info(matrix=evecs(1)%matrix, matrix_struct=fm_struct)
      NULLIFY (scaled_evect)
      CALL cp_fm_create(scaled_evect, fm_struct)

      ALLOCATE (weights_2d(nmo_virt, nmo_occ))
      DO imo_occ = 1, nmo_occ
         DO imo_virt = 1, nmo_virt
            weights_2d(imo_virt, imo_occ) = 1.0_dp/(evals_virt(imo_virt)-evals_occ(imo_occ))
         END DO
      END DO
      CALL cp_fm_set_submatrix(ediff_inv, weights_2d)

      ! summary information
      IF (log_unit > 0) THEN
         WRITE (log_unit, '(/,1X,A,/)') "TDDFT singlet states"

         WRITE (log_unit, '(3X,A,T14,A,T31,A,T63,A)') "State", "Excitation", "Transition dipole (a.u.)", "Oscillator"
         WRITE (log_unit, '(3X,A,T14,A,T33,A,T43,A,T53,A,T61,A)') "number", "energy (eV)", "x", "y", "z", "strength (a.u.)"
         WRITE (log_unit, '(1X,74("-"))')
      END IF

      DO istate = 1, nstates
         CALL cp_gemm('T', 'N', nmo_virt, nmo_occ, nao, 1.0_dp, S_mos_virt, evecs(istate)%matrix, 0.0_dp, weights_fm)

         CALL cp_fm_schur_product(weights_fm, ediff_inv, ediff_inv_weights)
         CALL cp_gemm('N', 'N', nao, nmo_occ, nmo_virt, 1.0_dp, mos_virt, ediff_inv_weights, 0.0_dp, scaled_evect)

         DO ideriv = 1, nderivs
            CALL cp_fm_trace(dS_mos_occ(ideriv)%matrix, scaled_evect, trans_dipole(ideriv))
         END DO
         oscillator_strength = 2.0_dp/3.0_dp*evals(istate)*DOT_PRODUCT(trans_dipole, trans_dipole)

         IF (log_unit > 0) &
            WRITE (log_unit, '(1X,I8,T14,F11.5,T28,3(1X,F9.4),T63,F10.5)') &
            istate, evals(istate)*evolt, trans_dipole(1:nderivs), oscillator_strength
      END DO

      CALL cp_fm_release(scaled_evect)
      CALL cp_fm_release(ediff_inv_weights)
      CALL cp_fm_release(ediff_inv)

      ! excitation analysis
      IF (min_amplitude < 1.0_dp) THEN
         nmo_virt_occ = nmo_virt*nmo_occ
         ALLOCATE (weights_abs_1d(nmo_virt_occ))
         ALLOCATE (inds(nmo_virt_occ))

         IF (log_unit > 0) THEN
            WRITE (log_unit, '(/,1X,A,/)') "Excitation analysis"

            WRITE (log_unit, '(3X,A,T15,A,T28,A,T41,A)') "State", "Occupied", "Virtual", "Excitation"
            WRITE (log_unit, '(3X,A,T16,A,T28,A,T41,A)') "number", "orbital", "orbital", "amplitude"
            WRITE (log_unit, '(1X,49("-"))')
         END IF

         DO istate = 1, nstates
            CALL cp_gemm('T', 'N', nmo_virt, nmo_occ, nao, 1.0_dp, S_mos_virt, evecs(istate)%matrix, 0.0_dp, weights_fm)
            CALL cp_fm_get_submatrix(weights_fm, weights_2d)
            CALL dcopy(nmo_virt_occ, weights_2d, 1, weights_abs_1d, 1)
            weights_abs_1d = ABS(weights_abs_1d)
            CALL sort(weights_abs_1d, nmo_virt_occ, inds)

            IF (log_unit > 0) &
               WRITE (log_unit, '(1X,I8)') istate

            DO iexc = nmo_virt_occ, 1, -1
               IF (weights_abs_1d(iexc) < min_amplitude) EXIT

               imo_occ = (inds(iexc)-1)/nmo_virt+1
               imo_virt = MOD(inds(iexc)-1, nmo_virt)+1

               IF (log_unit > 0) &
                  WRITE (log_unit, '(T14,I8,T26,I8,T42,F9.6)') nmo_occ_total-nmo_occ+imo_occ, &
                  nmo_occ_total+imo_virt, weights_2d(imo_virt, imo_occ)
            END DO
         END DO

         DEALLOCATE (inds)
         DEALLOCATE (weights_abs_1d)
      END IF

      DEALLOCATE (weights_2d)
      CALL cp_fm_release(weights_fm)
      DO ideriv = 1, nderivs
         CALL cp_fm_release(dS_mos_occ(ideriv)%matrix)
      END DO
      CALL cp_fm_release(S_mos_virt)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_print_summary

END MODULE qs_tddfpt2_methods
