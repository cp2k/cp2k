!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for a linear scaling quickstep SCF run based on the density
!>        matrix
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! **************************************************************************************************
MODULE dm_ls_scf
   USE arnoldi_api,                     ONLY: arnoldi_extremal
   USE bibliography,                    ONLY: Kolafa2004,&
                                              Kuhne2007,&
                                              cite_reference
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_external_control,             ONLY: external_control
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_binary_read, dbcsr_binary_write, dbcsr_checksum, dbcsr_copy, &
        dbcsr_create, dbcsr_distribution_type, dbcsr_filter, dbcsr_get_info, dbcsr_get_occupation, &
        dbcsr_multiply, dbcsr_p_type, dbcsr_release, dbcsr_scale, dbcsr_set, dbcsr_type, &
        dbcsr_type_no_symmetry
   USE dm_ls_chebyshev,                 ONLY: compute_chebyshev
   USE dm_ls_scf_curvy,                 ONLY: deallocate_curvy_data,&
                                              dm_ls_curvy_optimization
   USE dm_ls_scf_methods,               ONLY: apply_matrix_preconditioner,&
                                              compute_homo_lumo,&
                                              density_matrix_sign,&
                                              density_matrix_sign_fixed_mu,&
                                              density_matrix_tc2,&
                                              density_matrix_trs4,&
                                              ls_scf_init_matrix_S
   USE dm_ls_scf_qs,                    ONLY: ls_scf_dm_to_ks,&
                                              ls_scf_init_qs,&
                                              ls_scf_qs_atomic_guess,&
                                              matrix_ls_create,&
                                              matrix_ls_to_qs,&
                                              matrix_qs_to_ls,&
                                              rho_mixing_ls_init
   USE dm_ls_scf_types,                 ONLY: ls_scf_env_type
   USE input_constants,                 ONLY: ls_cluster_atomic,&
                                              ls_scf_pexsi,&
                                              ls_scf_sign,&
                                              ls_scf_tc2,&
                                              ls_scf_trs4,&
                                              transport_transmission
   USE input_section_types,             ONLY: section_vals_type
   USE iterate_matrix,                  ONLY: purify_mcweeny
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE mathlib,                         ONLY: binomial
   USE molecule_types,                  ONLY: molecule_type
   USE pao_main,                        ONLY: pao_optimization_end,&
                                              pao_optimization_start,&
                                              pao_post_scf,&
                                              pao_update
   USE pexsi_methods,                   ONLY: density_matrix_pexsi,&
                                              pexsi_finalize_scf,&
                                              pexsi_init_scf,&
                                              pexsi_set_convergence_tolerance,&
                                              pexsi_to_qs
   USE qs_diis,                         ONLY: qs_diis_b_clear_sparse,&
                                              qs_diis_b_create_sparse,&
                                              qs_diis_b_step_4lscf
   USE qs_diis_types,                   ONLY: qs_diis_b_release_sparse,&
                                              qs_diis_buffer_type_sparse
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_scf_post_gpw,                 ONLY: qs_scf_post_moments,&
                                              write_mo_free_results
   USE qs_scf_post_tb,                  ONLY: scf_post_calculation_tb
   USE transport,                       ONLY: external_scf_method,&
                                              transport_initialize
   USE transport_env_types,             ONLY: transport_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dm_ls_scf'

   PUBLIC :: ls_scf

CONTAINS

! **************************************************************************************************
!> \brief perform an linear scaling scf procedure: entry point
!>
!> \param qs_env ...
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! **************************************************************************************************
   SUBROUTINE ls_scf(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'ls_scf', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      LOGICAL                                            :: pao_is_done
      TYPE(ls_scf_env_type), POINTER                     :: ls_scf_env

      CALL timeset(routineN, handle)
      CALL get_qs_env(qs_env, ls_scf_env=ls_scf_env)
      CALL pao_optimization_start(qs_env, ls_scf_env)

      pao_is_done = .FALSE.
      DO WHILE (.NOT. pao_is_done)
         CALL ls_scf_init_scf(qs_env, ls_scf_env)
         CALL pao_update(qs_env, ls_scf_env, pao_is_done)
         CALL ls_scf_main(qs_env, ls_scf_env)
         CALL pao_post_scf(qs_env, ls_scf_env, pao_is_done)
         CALL ls_scf_post(qs_env, ls_scf_env)
      END DO

      CALL pao_optimization_end(ls_scf_env)

      CALL timestop(handle)

   END SUBROUTINE ls_scf

! **************************************************************************************************
!> \brief initialization needed for scf
!> \param qs_env ...
!> \param ls_scf_env ...
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! **************************************************************************************************
   SUBROUTINE ls_scf_init_scf(qs_env, ls_scf_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ls_scf_env_type)                              :: ls_scf_env

      CHARACTER(len=*), PARAMETER :: routineN = 'ls_scf_init_scf', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nspin, unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, matrix_w
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(section_vals_type), POINTER                   :: input

      CALL timeset(routineN, handle)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      ! get basic quantities from the qs_env
      CALL get_qs_env(qs_env, nelectron_total=ls_scf_env%nelectron_total, &
                      matrix_s=matrix_s, &
                      matrix_w=matrix_w, &
                      ks_env=ks_env, &
                      dft_control=dft_control, &
                      molecule_set=molecule_set, &
                      input=input, &
                      has_unit_metric=ls_scf_env%has_unit_metric, &
                      para_env=ls_scf_env%para_env, &
                      nelectron_spin=ls_scf_env%nelectron_spin)

      ! needs forces ? There might be a better way to flag this
      ls_scf_env%calculate_forces = ASSOCIATED(matrix_w)

      ! some basic initialization of the QS side of things
      CALL ls_scf_init_qs(qs_env)

      ! create the matrix template for use in the ls procedures
      CALL matrix_ls_create(matrix_ls=ls_scf_env%matrix_s, matrix_qs=matrix_s(1)%matrix, &
                            ls_mstruct=ls_scf_env%ls_mstruct)

      nspin = ls_scf_env%nspins
      IF (ALLOCATED(ls_scf_env%matrix_p)) THEN
         DO ispin = 1, SIZE(ls_scf_env%matrix_p)
            CALL dbcsr_release(ls_scf_env%matrix_p(ispin))
         ENDDO
      ELSE
         ALLOCATE (ls_scf_env%matrix_p(nspin))
      END IF

      DO ispin = 1, nspin
         CALL dbcsr_create(ls_scf_env%matrix_p(ispin), template=ls_scf_env%matrix_s, &
                           matrix_type=dbcsr_type_no_symmetry)
      ENDDO

      ALLOCATE (ls_scf_env%matrix_ks(nspin))
      DO ispin = 1, nspin
         CALL dbcsr_create(ls_scf_env%matrix_ks(ispin), template=ls_scf_env%matrix_s, &
                           matrix_type=dbcsr_type_no_symmetry)
      ENDDO

      ! set up matrix S, and needed functions of S
      CALL ls_scf_init_matrix_s(matrix_s(1)%matrix, ls_scf_env)

      ! get the initial guess for the SCF
      CALL ls_scf_initial_guess(qs_env, ls_scf_env)

      IF (ls_scf_env%do_rho_mixing) THEN
         CALL rho_mixing_ls_init(qs_env, ls_scf_env)
      ENDIF

      IF (ls_scf_env%do_pexsi) THEN
         CALL pexsi_init_scf(ks_env, ls_scf_env%pexsi, matrix_s(1)%matrix)
      ENDIF

      IF (qs_env%do_transport) THEN
         CALL transport_initialize(ks_env, qs_env%transport_env, matrix_s(1)%matrix)
      END IF

      CALL timestop(handle)

   END SUBROUTINE ls_scf_init_scf

! **************************************************************************************************
!> \brief deal with the scf initial guess
!> \param qs_env ...
!> \param ls_scf_env ...
!> \par History
!>       2012.11 created [Joost VandeVondele]
!> \author Joost VandeVondele
! **************************************************************************************************
   SUBROUTINE ls_scf_initial_guess(qs_env, ls_scf_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ls_scf_env_type)                              :: ls_scf_env

      CHARACTER(len=*), PARAMETER :: routineN = 'ls_scf_initial_guess', &
         routineP = moduleN//':'//routineN
      INTEGER, PARAMETER                                 :: aspc_guess = 2, atomic_guess = 1, &
                                                            restart_guess = 3

      CHARACTER(LEN=default_path_length)                 :: file_name, project_name
      INTEGER                                            :: handle, iaspc, initial_guess_type, &
                                                            ispin, istore, naspc, unit_nr
      REAL(KIND=dp)                                      :: alpha, cs_pos
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_distribution_type)                      :: dist
      TYPE(dbcsr_type)                                   :: matrix_tmp1

      IF (ls_scf_env%do_pao) RETURN ! pao has its own initial guess

      CALL timeset(routineN, handle)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      IF (unit_nr > 0) WRITE (unit_nr, '()')
      ! if there is no history go for the atomic guess, otherwise extrapolate the dm history
      IF (ls_scf_env%scf_history%istore == 0) THEN
         IF (ls_scf_env%restart_read) THEN
            initial_guess_type = restart_guess
         ELSE
            initial_guess_type = atomic_guess
         ENDIF
      ELSE
         initial_guess_type = aspc_guess
      ENDIF

      ! how to get the initial guess
      SELECT CASE (initial_guess_type)
      CASE (atomic_guess)
         CALL ls_scf_qs_atomic_guess(qs_env, ls_scf_env%energy_init)
      CASE (restart_guess)
         project_name = logger%iter_info%project_name
         DO ispin = 1, SIZE(ls_scf_env%matrix_p)
            WRITE (file_name, '(A,I0,A)') TRIM(project_name)//"_LS_DM_SPIN_", ispin, "_RESTART.dm"
            CALL dbcsr_get_info(ls_scf_env%matrix_p(1), distribution=dist)
            CALL dbcsr_binary_read(file_name, distribution=dist, matrix_new=ls_scf_env%matrix_p(ispin))
            cs_pos = dbcsr_checksum(ls_scf_env%matrix_p(ispin), pos=.TRUE.)
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(T2,A,E20.8)') "Read restart DM "//TRIM(file_name)//" with checksum: ", cs_pos
            ENDIF
         ENDDO

         ! directly go to computing the corresponding energy and ks matrix
         CALL ls_scf_dm_to_ks(qs_env, ls_scf_env, ls_scf_env%energy_init, iscf=0)
      CASE (aspc_guess)
         CALL cite_reference(Kolafa2004)
         CALL cite_reference(Kuhne2007)
         naspc = MIN(ls_scf_env%scf_history%istore, ls_scf_env%scf_history%nstore)
         DO ispin = 1, SIZE(ls_scf_env%matrix_p)
            ! actual extrapolation
            CALL dbcsr_set(ls_scf_env%matrix_p(ispin), 0.0_dp)
            DO iaspc = 1, naspc
               alpha = (-1.0_dp)**(iaspc + 1)*REAL(iaspc, KIND=dp)* &
                       binomial(2*naspc, naspc - iaspc)/binomial(2*naspc - 2, naspc - 1)
               istore = MOD(ls_scf_env%scf_history%istore - iaspc, ls_scf_env%scf_history%nstore) + 1
               CALL dbcsr_add(ls_scf_env%matrix_p(ispin), ls_scf_env%scf_history%matrix(ispin, istore), 1.0_dp, alpha)
            ENDDO
         ENDDO
      END SELECT

      ! which cases need getting purified and non-orthogonal ?
      SELECT CASE (initial_guess_type)
      CASE (atomic_guess, restart_guess)
         ! do nothing
      CASE (aspc_guess)
         ! purification can't be done on the pexsi matrix, which is not necessarily idempotent,
         ! and not stored in an ortho basis form
         IF (.NOT. (ls_scf_env%do_pexsi)) THEN
            DO ispin = 1, SIZE(ls_scf_env%matrix_p)
               ! linear combination of P's is not idempotent. A bit of McWeeny is needed to ensure it is again
               IF (SIZE(ls_scf_env%matrix_p) == 1) CALL dbcsr_scale(ls_scf_env%matrix_p(ispin), 0.5_dp)
               ! to ensure that noisy blocks do not build up during MD (in particular with curvy) filter that guess a bit more
               CALL dbcsr_filter(ls_scf_env%matrix_p(ispin), ls_scf_env%eps_filter**(2.0_dp/3.0_dp))
               CALL purify_mcweeny(ls_scf_env%matrix_p(ispin:ispin), ls_scf_env%eps_filter, 3)
               IF (SIZE(ls_scf_env%matrix_p) == 1) CALL dbcsr_scale(ls_scf_env%matrix_p(ispin), 2.0_dp)

               IF (ls_scf_env%use_s_sqrt) THEN
                  ! need to get P in the non-orthogonal basis if it was stored differently
                  CALL dbcsr_create(matrix_tmp1, template=ls_scf_env%matrix_s, &
                                    matrix_type=dbcsr_type_no_symmetry)
                  CALL dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_s_sqrt_inv, ls_scf_env%matrix_p(ispin), &
                                      0.0_dp, matrix_tmp1, filter_eps=ls_scf_env%eps_filter)
                  CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, ls_scf_env%matrix_s_sqrt_inv, &
                                      0.0_dp, ls_scf_env%matrix_p(ispin), &
                                      filter_eps=ls_scf_env%eps_filter)
                  CALL dbcsr_release(matrix_tmp1)

                  IF (ls_scf_env%has_s_preconditioner) THEN
                     CALL apply_matrix_preconditioner(ls_scf_env%matrix_p(ispin), "forward", &
                                                      ls_scf_env%matrix_bs_sqrt, ls_scf_env%matrix_bs_sqrt_inv)
                  ENDIF
               ENDIF
            ENDDO
         ENDIF

         ! compute corresponding energy and ks matrix
         CALL ls_scf_dm_to_ks(qs_env, ls_scf_env, ls_scf_env%energy_init, iscf=0)
      END SELECT

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A,F20.9)') "Energy with the initial guess:", ls_scf_env%energy_init
         WRITE (unit_nr, '()')
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE ls_scf_initial_guess

! **************************************************************************************************
!> \brief store a history of matrices for later use in ls_scf_initial_guess
!> \param ls_scf_env ...
!> \par History
!>       2012.11 created [Joost VandeVondele]
!> \author Joost VandeVondele
! **************************************************************************************************
   SUBROUTINE ls_scf_store_result(ls_scf_env)
      TYPE(ls_scf_env_type)                              :: ls_scf_env

      CHARACTER(len=*), PARAMETER :: routineN = 'ls_scf_store_result', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_path_length)                 :: file_name, project_name
      INTEGER                                            :: handle, ispin, istore, unit_nr
      REAL(KIND=dp)                                      :: cs_pos
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type)                                   :: matrix_tmp1

      CALL timeset(routineN, handle)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      IF (ls_scf_env%restart_write) THEN
         DO ispin = 1, SIZE(ls_scf_env%matrix_p)
            project_name = logger%iter_info%project_name
            WRITE (file_name, '(A,I0,A)') TRIM(project_name)//"_LS_DM_SPIN_", ispin, "_RESTART.dm"
            cs_pos = dbcsr_checksum(ls_scf_env%matrix_p(ispin), pos=.TRUE.)
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(T2,A,E20.8)') "Writing restart DM "//TRIM(file_name)//" with checksum: ", cs_pos
            ENDIF
            IF (ls_scf_env%do_transport .OR. ls_scf_env%do_pexsi) THEN
               IF (unit_nr > 0) THEN
                  WRITE (unit_nr, '(T6,A)') "The restart DM "//TRIM(file_name)//" has the sparsity of S, therefore,"
                  WRITE (unit_nr, '(T6,A)') "not compatible with methods that require a full DM! "
               END IF
            END IF
            CALL dbcsr_binary_write(ls_scf_env%matrix_p(ispin), file_name)
         ENDDO
      END IF

      IF (ls_scf_env%scf_history%nstore > 0) THEN
         ls_scf_env%scf_history%istore = ls_scf_env%scf_history%istore + 1
         DO ispin = 1, SIZE(ls_scf_env%matrix_p)
            istore = MOD(ls_scf_env%scf_history%istore - 1, ls_scf_env%scf_history%nstore) + 1
            CALL dbcsr_copy(ls_scf_env%scf_history%matrix(ispin, istore), ls_scf_env%matrix_p(ispin))

            ! if we have the sqrt around, we use it to go to the orthogonal basis
            IF (ls_scf_env%use_s_sqrt) THEN
               ! usually sqrt(S) * P * sqrt(S) should be available, or could be stored at least,
               ! so that the next multiplications could be saved.
               CALL dbcsr_create(matrix_tmp1, template=ls_scf_env%matrix_s, &
                                 matrix_type=dbcsr_type_no_symmetry)

               IF (ls_scf_env%has_s_preconditioner) THEN
                  CALL apply_matrix_preconditioner(ls_scf_env%scf_history%matrix(ispin, istore), "backward", &
                                                   ls_scf_env%matrix_bs_sqrt, ls_scf_env%matrix_bs_sqrt_inv)
               ENDIF
               CALL dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_s_sqrt, ls_scf_env%scf_history%matrix(ispin, istore), &
                                   0.0_dp, matrix_tmp1, filter_eps=ls_scf_env%eps_filter)
               CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, ls_scf_env%matrix_s_sqrt, &
                                   0.0_dp, ls_scf_env%scf_history%matrix(ispin, istore), &
                                   filter_eps=ls_scf_env%eps_filter)
               CALL dbcsr_release(matrix_tmp1)
            ENDIF

         ENDDO
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE ls_scf_store_result

! **************************************************************************************************
!> \brief Main SCF routine. Can we keep it clean ?
!> \param qs_env ...
!> \param ls_scf_env ...
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! **************************************************************************************************
   SUBROUTINE ls_scf_main(qs_env, ls_scf_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ls_scf_env_type)                              :: ls_scf_env

      CHARACTER(len=*), PARAMETER :: routineN = 'ls_scf_main', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iscf, ispin, &
                                                            nelectron_spin_real, nmixing, nspin, &
                                                            unit_nr
      LOGICAL :: check_convergence, diis_step, do_transport, extra_scf, maxscf_reached, &
         scf_converged, should_stop, transm_maxscf_reached, transm_scf_converged
      REAL(KIND=dp)                                      :: energy_diff, energy_new, energy_old, &
                                                            eps_diis, t1, t2
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:)        :: matrix_ks_deviation, matrix_mixing_old
      TYPE(qs_diis_buffer_type_sparse), POINTER          :: diis_buffer
      TYPE(transport_env_type), POINTER                  :: transport_env

      CALL timeset(routineN, handle)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      nspin = ls_scf_env%nspins

      ! old quantities, useful for mixing
      ALLOCATE (matrix_mixing_old(nspin), matrix_ks_deviation(nspin))
      DO ispin = 1, nspin
         CALL dbcsr_create(matrix_mixing_old(ispin), template=ls_scf_env%matrix_ks(ispin))

         CALL dbcsr_create(matrix_ks_deviation(ispin), template=ls_scf_env%matrix_ks(ispin))
         CALL dbcsr_set(matrix_ks_deviation(ispin), 0.0_dp)
      ENDDO
      ls_scf_env%homo_spin(:) = 0.0_dp
      ls_scf_env%lumo_spin(:) = 0.0_dp

      transm_scf_converged = .FALSE.
      transm_maxscf_reached = .FALSE.

      energy_old = 0.0_dp
      IF (ls_scf_env%scf_history%istore > 0) energy_old = ls_scf_env%energy_init
      check_convergence = .TRUE.
      iscf = 0
      IF (ls_scf_env%ls_diis) THEN
         diis_step = .FALSE.
         eps_diis = ls_scf_env%eps_diis
         nmixing = ls_scf_env%nmixing
         NULLIFY (diis_buffer)
         IF (.NOT. ASSOCIATED(diis_buffer)) THEN
            CALL qs_diis_b_create_sparse(diis_buffer, &
                                         nbuffer=ls_scf_env%max_diis)
         END IF
         CALL qs_diis_b_clear_sparse(diis_buffer)
         CALL get_qs_env(qs_env, matrix_s=matrix_s)
      END IF

      CALL get_qs_env(qs_env, transport_env=transport_env, do_transport=do_transport)

      ! the real SCF loop
      DO

         ! check on max SCF or timing/exit
         CALL external_control(should_stop, "SCF", start_time=qs_env%start_time, target_time=qs_env%target_time)
         IF (do_transport) THEN
            maxscf_reached = should_stop .OR. iscf >= ls_scf_env%max_scf
            ! one extra scf step for post-processing in transmission calculations
            IF (transport_env%params%method .EQ. transport_transmission) THEN
               IF (transm_maxscf_reached) THEN
                  IF (unit_nr > 0) WRITE (unit_nr, '(T2,A)') "SCF not converged! "
                  EXIT
               END IF
               transm_maxscf_reached = maxscf_reached
            ELSE
               IF (maxscf_reached) THEN
                  IF (unit_nr > 0) WRITE (unit_nr, '(T2,A)') "SCF not converged! "
                  EXIT
               ENDIF
            END IF
         ELSE
            IF (should_stop .OR. iscf >= ls_scf_env%max_scf) THEN
               IF (unit_nr > 0) WRITE (unit_nr, '(T2,A)') "SCF not converged! "
               EXIT
            ENDIF
         END IF

         t1 = m_walltime()
         iscf = iscf + 1

         ! first get a copy of the current KS matrix
         CALL get_qs_env(qs_env, matrix_ks=matrix_ks)
         DO ispin = 1, nspin
            CALL matrix_qs_to_ls(ls_scf_env%matrix_ks(ispin), matrix_ks(ispin)%matrix, &
                                 ls_scf_env%ls_mstruct, covariant=.TRUE.)
            IF (ls_scf_env%has_s_preconditioner) THEN
               CALL apply_matrix_preconditioner(ls_scf_env%matrix_ks(ispin), "forward", &
                                                ls_scf_env%matrix_bs_sqrt, ls_scf_env%matrix_bs_sqrt_inv)
            ENDIF
            CALL dbcsr_filter(ls_scf_env%matrix_ks(ispin), ls_scf_env%eps_filter)
         ENDDO
         ! run curvy steps if required. Needs an idempotent DM (either perification or restart)
         IF ((iscf > 1 .OR. ls_scf_env%scf_history%istore > 0) .AND. ls_scf_env%curvy_steps) THEN
            CALL dm_ls_curvy_optimization(ls_scf_env, energy_old, check_convergence)
         ELSE
            ! turn the KS matrix in a density matrix
            DO ispin = 1, nspin
               IF (ls_scf_env%do_rho_mixing) THEN
                  CALL dbcsr_copy(matrix_mixing_old(ispin), ls_scf_env%matrix_ks(ispin))
               ELSE
                  IF (iscf == 1) THEN
                     ! initialize the mixing matrix with the current state if needed
                     CALL dbcsr_copy(matrix_mixing_old(ispin), ls_scf_env%matrix_ks(ispin))
                  ELSE
                     IF (ls_scf_env%ls_diis) THEN ! ------- IF-DIIS+MIX--- START
                        IF (diis_step .AND. (iscf - 1) .GE. ls_scf_env%iter_ini_diis) THEN
                           IF (unit_nr > 0) THEN
                              WRITE (unit_nr, '(A61)') &
                                 '*************************************************************'
                              WRITE (unit_nr, '(A50,2(I3,A1),L1,A1)') &
                                 " Using DIIS mixed KS:  (iscf,INI_DIIS,DIIS_STEP)=(", &
                                 iscf, ",", ls_scf_env%iter_ini_diis, ",", diis_step, ")"
                              WRITE (unit_nr, '(A52)') &
                                 " KS_nw= DIIS-Linear-Combination-Previous KS matrices"
                              WRITE (unit_nr, '(61A)') &
                                 "*************************************************************"
                           ENDIF
                           CALL dbcsr_copy(matrix_mixing_old(ispin), & ! out
                                           ls_scf_env%matrix_ks(ispin)) ! in
                        ELSE
                           IF (unit_nr > 0) THEN
                              WRITE (unit_nr, '(A57)') &
                                 "*********************************************************"
                              WRITE (unit_nr, '(A23,F5.3,A25,I3)') &
                                 " Using MIXING_FRACTION=", ls_scf_env%mixing_fraction, &
                                 " to mix KS matrix:  iscf=", iscf
                              WRITE (unit_nr, '(A7,F5.3,A6,F5.3,A7)') &
                                 " KS_nw=", ls_scf_env%mixing_fraction, "*KS + ", &
                                 1.0_dp - ls_scf_env%mixing_fraction, "*KS_old"
                              WRITE (unit_nr, '(A57)') &
                                 "*********************************************************"
                           ENDIF
                           ! perform the mixing of ks matrices
                           CALL dbcsr_add(matrix_mixing_old(ispin), &
                                          ls_scf_env%matrix_ks(ispin), &
                                          1.0_dp - ls_scf_env%mixing_fraction, &
                                          ls_scf_env%mixing_fraction)
                        ENDIF
                     ELSE ! otherwise
                        IF (unit_nr > 0) THEN
                           WRITE (unit_nr, '(A57)') &
                              "*********************************************************"
                           WRITE (unit_nr, '(A23,F5.3,A25,I3)') &
                              " Using MIXING_FRACTION=", ls_scf_env%mixing_fraction, &
                              " to mix KS matrix:  iscf=", iscf
                           WRITE (unit_nr, '(A7,F5.3,A6,F5.3,A7)') &
                              " KS_nw=", ls_scf_env%mixing_fraction, "*KS + ", &
                              1.0_dp - ls_scf_env%mixing_fraction, "*KS_old"
                           WRITE (unit_nr, '(A57)') &
                              "*********************************************************"
                        ENDIF
                        ! perform the mixing of ks matrices
                        CALL dbcsr_add(matrix_mixing_old(ispin), &
                                       ls_scf_env%matrix_ks(ispin), &
                                       1.0_dp - ls_scf_env%mixing_fraction, &
                                       ls_scf_env%mixing_fraction)
                     ENDIF ! ------- IF-DIIS+MIX--- END
                  ENDIF
               ENDIF

               ! compute the density matrix that matches it
               ! we need the proper number of states
               nelectron_spin_real = ls_scf_env%nelectron_spin(ispin)
               IF (ls_scf_env%nspins == 1) nelectron_spin_real = nelectron_spin_real/2

               IF (do_transport) THEN
                  IF (ls_scf_env%has_s_preconditioner) &
                     CPABORT("NOT YET IMPLEMENTED with S preconditioner. ")
                  IF (ls_scf_env%ls_mstruct%cluster_type .NE. ls_cluster_atomic) &
                     CPABORT("NOT YET IMPLEMENTED with molecular clustering. ")

                  extra_scf = maxscf_reached .OR. scf_converged
                  ! get the current Kohn-Sham matrix (ks) and return matrix_p evaluated using an external C routine
                  CALL external_scf_method(transport_env, ls_scf_env%matrix_s, matrix_mixing_old(ispin), &
                                           ls_scf_env%matrix_p(ispin), nelectron_spin_real, ls_scf_env%natoms, &
                                           energy_diff, iscf, extra_scf)

               ELSE
                  SELECT CASE (ls_scf_env%purification_method)
                  CASE (ls_scf_sign)
                     CALL density_matrix_sign(ls_scf_env%matrix_p(ispin), ls_scf_env%mu_spin(ispin), ls_scf_env%fixed_mu, &
                                              ls_scf_env%sign_method, ls_scf_env%sign_order, matrix_mixing_old(ispin), &
                                              ls_scf_env%matrix_s, ls_scf_env%matrix_s_inv, nelectron_spin_real, &
                                              ls_scf_env%eps_filter, ls_scf_env%sign_symmetric, ls_scf_env%submatrix_sign_method, &
                                              ls_scf_env%matrix_s_sqrt_inv)
                  CASE (ls_scf_tc2)
                     CALL density_matrix_tc2(ls_scf_env%matrix_p(ispin), matrix_mixing_old(ispin), ls_scf_env%matrix_s_sqrt_inv, &
                                             nelectron_spin_real, ls_scf_env%eps_filter, ls_scf_env%homo_spin(ispin), &
                                             ls_scf_env%lumo_spin(ispin), non_monotonic=ls_scf_env%non_monotonic, &
                                             eps_lanczos=ls_scf_env%eps_lanczos, max_iter_lanczos=ls_scf_env%max_iter_lanczos)
                  CASE (ls_scf_trs4)
                     CALL density_matrix_trs4(ls_scf_env%matrix_p(ispin), matrix_mixing_old(ispin), ls_scf_env%matrix_s_sqrt_inv, &
                                              nelectron_spin_real, ls_scf_env%eps_filter, ls_scf_env%homo_spin(ispin), &
                                              ls_scf_env%lumo_spin(ispin), ls_scf_env%mu_spin(ispin), &
                                              dynamic_threshold=ls_scf_env%dynamic_threshold, &
                                              matrix_ks_deviation=matrix_ks_deviation(ispin), &
                                              eps_lanczos=ls_scf_env%eps_lanczos, max_iter_lanczos=ls_scf_env%max_iter_lanczos)
                  CASE (ls_scf_pexsi)
                     IF (ls_scf_env%has_s_preconditioner) &
                        CPABORT("S preconditioning not implemented in combination with the PEXSI library. ")
                     IF (ls_scf_env%ls_mstruct%cluster_type .NE. ls_cluster_atomic) &
                        CALL cp_abort(__LOCATION__, &
                                      "Molecular clustering not implemented in combination with the PEXSI library. ")
                     IF (ls_scf_env%ls_mstruct%single_precision) &
                        CPABORT("PEXSI library requires double precision datatype. ")
                     CALL density_matrix_pexsi(ls_scf_env%pexsi, ls_scf_env%matrix_p(ispin), ls_scf_env%pexsi%matrix_w(ispin), &
                                               ls_scf_env%pexsi%kTS(ispin), matrix_mixing_old(ispin), ls_scf_env%matrix_s, &
                                               nelectron_spin_real, ls_scf_env%mu_spin(ispin), iscf, ispin)
                  END SELECT
               END IF

               IF (ls_scf_env%has_s_preconditioner) THEN
                  CALL apply_matrix_preconditioner(ls_scf_env%matrix_p(ispin), "forward", &
                                                   ls_scf_env%matrix_bs_sqrt, ls_scf_env%matrix_bs_sqrt_inv)
               ENDIF
               CALL dbcsr_filter(ls_scf_env%matrix_p(ispin), ls_scf_env%eps_filter)

               IF (ls_scf_env%nspins == 1) CALL dbcsr_scale(ls_scf_env%matrix_p(ispin), 2.0_dp)

            ENDDO
         END IF

         ! compute the corresponding new energy KS matrix and new energy
         CALL ls_scf_dm_to_ks(qs_env, ls_scf_env, energy_new, iscf)

         IF (ls_scf_env%do_pexsi) THEN
            CALL pexsi_to_qs(ls_scf_env, qs_env, kTS=ls_scf_env%pexsi%kTS)
         ENDIF

         ! report current SCF loop
         energy_diff = energy_new - energy_old
         energy_old = energy_new

         t2 = m_walltime()
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, *)
            WRITE (unit_nr, '(T2,A,I6,F20.9,F20.9,F12.6)') "SCF", iscf, energy_new, energy_diff, t2 - t1
            WRITE (unit_nr, *)
            CALL m_flush(unit_nr)
         ENDIF

         IF (do_transport) THEN
            scf_converged = check_convergence .AND. ABS(energy_diff) < ls_scf_env%eps_scf*ls_scf_env%nelectron_total
            ! one extra scf step for post-processing in transmission calculations
            IF (transport_env%params%method .EQ. transport_transmission) THEN
               IF (transm_scf_converged) EXIT
               transm_scf_converged = scf_converged
            ELSE
               IF (scf_converged) EXIT
            END IF
         ELSE
            ! exit criterion on the energy only for the time being
            IF (check_convergence .AND. ABS(energy_diff) < ls_scf_env%eps_scf*ls_scf_env%nelectron_total) EXIT
         END IF

         IF (ls_scf_env%ls_diis) THEN
! diis_buffer, buffer with 1) Kohn-Sham history matrix,
!                          2) KS error history matrix (f=KPS-SPK),
!                          3) B matrix (for finding DIIS weighting coefficients)
            CALL qs_diis_b_step_4lscf(diis_buffer, qs_env, ls_scf_env, unit_nr, &
                                      iscf, diis_step, eps_diis, nmixing, matrix_s(1)%matrix, &
                                      ls_scf_env%eps_filter)
         ENDIF

         IF (ls_scf_env%do_pexsi) THEN
            CALL pexsi_set_convergence_tolerance(ls_scf_env%pexsi, energy_diff, &
                                                 ls_scf_env%eps_scf*ls_scf_env%nelectron_total, &
                                                 ! initialize in second scf step of first SCF cycle:
                                                 (iscf .EQ. 2) .AND. (ls_scf_env%scf_history%istore .EQ. 0), &
                                                 check_convergence)
         ENDIF

      ENDDO

      ! free storage
      IF (ls_scf_env%ls_diis) THEN
         CALL qs_diis_b_release_sparse(diis_buffer)
      ENDIF
      DO ispin = 1, nspin
         CALL dbcsr_release(matrix_mixing_old(ispin))
         CALL dbcsr_release(matrix_ks_deviation(ispin))
      ENDDO
      DEALLOCATE (matrix_mixing_old, matrix_ks_deviation)

      CALL timestop(handle)

   END SUBROUTINE ls_scf_main

! **************************************************************************************************
!> \brief after SCF we have a density matrix, and the self consistent KS matrix
!>        analyze its properties.
!> \param qs_env ...
!> \param ls_scf_env ...
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! **************************************************************************************************
   SUBROUTINE ls_scf_post(qs_env, ls_scf_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ls_scf_env_type)                              :: ls_scf_env

      CHARACTER(len=*), PARAMETER :: routineN = 'ls_scf_post', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, unit_nr
      REAL(KIND=dp)                                      :: occ
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_w
      TYPE(dft_control_type), POINTER                    :: dft_control

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, dft_control=dft_control)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      ! store the matrix for a next scf run
      IF (.NOT. ls_scf_env%do_pao) &
         CALL ls_scf_store_result(ls_scf_env)

      ! write homo and lumo energy and occupation (if not already part of the output)
      IF (ls_scf_env%curvy_steps) THEN
         CALL post_scf_homo_lumo(ls_scf_env)

         ! always report P occ
         IF (unit_nr > 0) WRITE (unit_nr, *) ""
         DO ispin = 1, ls_scf_env%nspins
            occ = dbcsr_get_occupation(ls_scf_env%matrix_p(ispin))
            IF (unit_nr > 0) WRITE (unit_nr, '(T2,A,F20.12)') "Density matrix (P) occupation ", occ
         ENDDO
      ENDIF

      ! compute the matrix_w if associated
      IF (ls_scf_env%calculate_forces) THEN
         CALL get_qs_env(qs_env, matrix_w=matrix_w)
         CPASSERT(ASSOCIATED(matrix_w))
         IF (ls_scf_env%do_pexsi) THEN
            CALL pexsi_to_qs(ls_scf_env, qs_env, matrix_w=ls_scf_env%pexsi%matrix_w)
         ELSE
            CALL calculate_w_matrix(matrix_w, ls_scf_env)
         ENDIF
      ENDIF

      ! compute properties

      IF (ls_scf_env%perform_mu_scan) CALL post_scf_mu_scan(ls_scf_env)

      IF (ls_scf_env%report_all_sparsities) CALL post_scf_sparsities(ls_scf_env)

      IF (dft_control%qs_control%dftb) THEN
         CALL scf_post_calculation_tb(qs_env, "DFTB", .TRUE.)
      ELSE IF (dft_control%qs_control%xtb) THEN
         CALL scf_post_calculation_tb(qs_env, "xTB", .TRUE.)
      ELSE
         CALL write_mo_free_results(qs_env)
      END IF

      IF (ls_scf_env%chebyshev%compute_chebyshev) CALL compute_chebyshev(qs_env, ls_scf_env)

      IF (.TRUE.) CALL post_scf_experiment()

      IF (dft_control%qs_control%dftb .OR. dft_control%qs_control%xtb) THEN
         !
      ELSE
         CALL qs_scf_post_moments(qs_env%input, logger, qs_env, unit_nr)
      END IF

      ! clean up used data

      CALL dbcsr_release(ls_scf_env%matrix_s)
      CALL deallocate_curvy_data(ls_scf_env%curvy_data)

      IF (ls_scf_env%has_s_preconditioner) THEN
         CALL dbcsr_release(ls_scf_env%matrix_bs_sqrt)
         CALL dbcsr_release(ls_scf_env%matrix_bs_sqrt_inv)
      ENDIF

      IF (ls_scf_env%needs_s_inv) THEN
         CALL dbcsr_release(ls_scf_env%matrix_s_inv)
      ENDIF

      IF (ls_scf_env%use_s_sqrt) THEN
         CALL dbcsr_release(ls_scf_env%matrix_s_sqrt)
         CALL dbcsr_release(ls_scf_env%matrix_s_sqrt_inv)
      ENDIF

      DO ispin = 1, SIZE(ls_scf_env%matrix_ks)
         CALL dbcsr_release(ls_scf_env%matrix_ks(ispin))
      ENDDO
      DEALLOCATE (ls_scf_env%matrix_ks)

      IF (ls_scf_env%do_pexsi) &
         CALL pexsi_finalize_scf(ls_scf_env%pexsi, ls_scf_env%mu_spin)

      CALL timestop(handle)

   END SUBROUTINE ls_scf_post

! **************************************************************************************************
!> \brief Compute the HOMO LUMO energies post SCF
!> \param ls_scf_env ...
!> \par History
!>       2013.06 created [Joost VandeVondele]
!> \author Joost VandeVondele
! **************************************************************************************************
   SUBROUTINE post_scf_homo_lumo(ls_scf_env)
      TYPE(ls_scf_env_type)                              :: ls_scf_env

      CHARACTER(len=*), PARAMETER :: routineN = 'post_scf_homo_lumo', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nspin, unit_nr
      LOGICAL                                            :: converged
      REAL(KIND=dp)                                      :: eps_max, eps_min, homo, lumo
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type)                                   :: matrix_k, matrix_p

      CALL timeset(routineN, handle)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      IF (unit_nr > 0) WRITE (unit_nr, '(T2,A)') ""

      ! TODO: remove these limitations
      CPASSERT(.NOT. ls_scf_env%has_s_preconditioner)
      CPASSERT(ls_scf_env%use_s_sqrt)

      nspin = ls_scf_env%nspins

      CALL dbcsr_create(matrix_p, template=ls_scf_env%matrix_p(1), matrix_type=dbcsr_type_no_symmetry)

      CALL dbcsr_create(matrix_k, template=ls_scf_env%matrix_p(1), matrix_type=dbcsr_type_no_symmetry)

      DO ispin = 1, nspin
         ! ortho basis ks
         CALL dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_s_sqrt_inv, ls_scf_env%matrix_ks(ispin), &
                             0.0_dp, matrix_k, filter_eps=ls_scf_env%eps_filter)
         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_k, ls_scf_env%matrix_s_sqrt_inv, &
                             0.0_dp, matrix_k, filter_eps=ls_scf_env%eps_filter)

         ! extremal eigenvalues ks
         CALL arnoldi_extremal(matrix_k, eps_max, eps_min, max_iter=ls_scf_env%max_iter_lanczos, &
                               threshold=ls_scf_env%eps_lanczos, converged=converged)

         ! ortho basis p
         CALL dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_s_sqrt, ls_scf_env%matrix_p(ispin), &
                             0.0_dp, matrix_p, filter_eps=ls_scf_env%eps_filter)
         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_p, ls_scf_env%matrix_s_sqrt, &
                             0.0_dp, matrix_p, filter_eps=ls_scf_env%eps_filter)
         IF (nspin == 1) CALL dbcsr_scale(matrix_p, 0.5_dp)

         ! go compute homo lumo
         CALL compute_homo_lumo(matrix_k, matrix_p, eps_min, eps_max, ls_scf_env%eps_filter, &
                                ls_scf_env%max_iter_lanczos, ls_scf_env%eps_lanczos, homo, lumo, unit_nr)

      ENDDO

      CALL dbcsr_release(matrix_p)
      CALL dbcsr_release(matrix_k)

      CALL timestop(handle)

   END SUBROUTINE post_scf_homo_lumo

! **************************************************************************************************
!> \brief Compute the density matrix for various values of the chemical potential
!> \param ls_scf_env ...
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! **************************************************************************************************
   SUBROUTINE post_scf_mu_scan(ls_scf_env)
      TYPE(ls_scf_env_type)                              :: ls_scf_env

      CHARACTER(len=*), PARAMETER :: routineN = 'post_scf_mu_scan', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, imu, ispin, nelectron_spin_real, &
                                                            nmu, nspin, unit_nr
      REAL(KIND=dp)                                      :: mu, t1, t2, trace
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type)                                   :: matrix_p

      CALL timeset(routineN, handle)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      nspin = ls_scf_env%nspins

      CALL dbcsr_create(matrix_p, template=ls_scf_env%matrix_p(1))

      nmu = 10
      DO imu = 0, nmu

         t1 = m_walltime()

         mu = -0.4_dp + imu*(0.1_dp + 0.4_dp)/nmu

         IF (unit_nr > 0) WRITE (unit_nr, *) "------- starting with mu ", mu

         DO ispin = 1, nspin
            ! we need the proper number of states
            nelectron_spin_real = ls_scf_env%nelectron_spin(ispin)
            IF (ls_scf_env%nspins == 1) nelectron_spin_real = nelectron_spin_real/2

            CALL density_matrix_sign_fixed_mu(matrix_p, trace, mu, ls_scf_env%sign_method, &
                                              ls_scf_env%sign_order, ls_scf_env%matrix_ks(ispin), &
                                              ls_scf_env%matrix_s, ls_scf_env%matrix_s_inv, &
                                              ls_scf_env%eps_filter, ls_scf_env%sign_symmetric, &
                                              ls_scf_env%submatrix_sign_method, ls_scf_env%matrix_s_sqrt_inv)
         ENDDO

         t2 = m_walltime()

         IF (unit_nr > 0) WRITE (unit_nr, *) " obtained ", mu, trace, t2 - t1

      ENDDO

      CALL dbcsr_release(matrix_p)

      CALL timestop(handle)

   END SUBROUTINE post_scf_mu_scan

! **************************************************************************************************
!> \brief Report on the sparsities of various interesting matrices.
!>
!> \param ls_scf_env ...
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! **************************************************************************************************
   SUBROUTINE post_scf_sparsities(ls_scf_env)
      TYPE(ls_scf_env_type)                              :: ls_scf_env

      CHARACTER(len=*), PARAMETER :: routineN = 'post_scf_sparsities', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_string_length)               :: title
      INTEGER                                            :: handle, ispin, nspin, unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type)                                   :: matrix_tmp1, matrix_tmp2

      CALL timeset(routineN, handle)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      nspin = ls_scf_env%nspins

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '()')
         WRITE (unit_nr, '(T2,A,E17.3)') "Sparsity reports for eps_filter: ", ls_scf_env%eps_filter
         WRITE (unit_nr, '()')
      ENDIF

      CALL report_matrix_sparsity(ls_scf_env%matrix_s, unit_nr, "overlap matrix (S)", &
                                  ls_scf_env%eps_filter)

      DO ispin = 1, nspin
         WRITE (title, '(A,I3)') "Kohn-Sham matrix (H) for spin ", ispin
         CALL report_matrix_sparsity(ls_scf_env%matrix_ks(ispin), unit_nr, title, &
                                     ls_scf_env%eps_filter)
      ENDDO

      CALL dbcsr_create(matrix_tmp1, template=ls_scf_env%matrix_s, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(matrix_tmp2, template=ls_scf_env%matrix_s, matrix_type=dbcsr_type_no_symmetry)

      DO ispin = 1, nspin
         WRITE (title, '(A,I3)') "Density matrix (P) for spin ", ispin
         CALL report_matrix_sparsity(ls_scf_env%matrix_p(ispin), unit_nr, title, &
                                     ls_scf_env%eps_filter)

         CALL dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_s, ls_scf_env%matrix_p(ispin), &
                             0.0_dp, matrix_tmp1, filter_eps=ls_scf_env%eps_filter)

         WRITE (title, '(A,I3,A)') "S * P(", ispin, ")"
         CALL report_matrix_sparsity(matrix_tmp1, unit_nr, title, ls_scf_env%eps_filter)

         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, ls_scf_env%matrix_s, &
                             0.0_dp, matrix_tmp2, filter_eps=ls_scf_env%eps_filter)
         WRITE (title, '(A,I3,A)') "S * P(", ispin, ") * S"
         CALL report_matrix_sparsity(matrix_tmp2, unit_nr, title, ls_scf_env%eps_filter)
      ENDDO

      IF (ls_scf_env%needs_s_inv) THEN
         CALL report_matrix_sparsity(ls_scf_env%matrix_s_inv, unit_nr, "inv(S)", &
                                     ls_scf_env%eps_filter)
         DO ispin = 1, nspin
            CALL dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_s_inv, ls_scf_env%matrix_ks(ispin), &
                                0.0_dp, matrix_tmp1, filter_eps=ls_scf_env%eps_filter)

            WRITE (title, '(A,I3,A)') "inv(S) * H(", ispin, ")"
            CALL report_matrix_sparsity(matrix_tmp1, unit_nr, title, ls_scf_env%eps_filter)
         ENDDO
      ENDIF

      IF (ls_scf_env%use_s_sqrt) THEN

         CALL report_matrix_sparsity(ls_scf_env%matrix_s_sqrt, unit_nr, "sqrt(S)", &
                                     ls_scf_env%eps_filter)
         CALL report_matrix_sparsity(ls_scf_env%matrix_s_sqrt_inv, unit_nr, "inv(sqrt(S))", &
                                     ls_scf_env%eps_filter)

         DO ispin = 1, nspin
            CALL dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_s_sqrt_inv, ls_scf_env%matrix_ks(ispin), &
                                0.0_dp, matrix_tmp1, filter_eps=ls_scf_env%eps_filter)
            CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, ls_scf_env%matrix_s_sqrt_inv, &
                                0.0_dp, matrix_tmp2, filter_eps=ls_scf_env%eps_filter)
            WRITE (title, '(A,I3,A)') "inv(sqrt(S)) * H(", ispin, ") * inv(sqrt(S))"
            CALL report_matrix_sparsity(matrix_tmp2, unit_nr, title, ls_scf_env%eps_filter)
         ENDDO

         DO ispin = 1, nspin
            CALL dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_s_sqrt, ls_scf_env%matrix_p(ispin), &
                                0.0_dp, matrix_tmp1, filter_eps=ls_scf_env%eps_filter)
            CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, ls_scf_env%matrix_s_sqrt, &
                                0.0_dp, matrix_tmp2, filter_eps=ls_scf_env%eps_filter)
            WRITE (title, '(A,I3,A)') "sqrt(S) * P(", ispin, ") * sqrt(S)"
            CALL report_matrix_sparsity(matrix_tmp2, unit_nr, title, ls_scf_env%eps_filter)
         ENDDO

      ENDIF

      CALL dbcsr_release(matrix_tmp1)
      CALL dbcsr_release(matrix_tmp2)

      CALL timestop(handle)

   END SUBROUTINE post_scf_sparsities

! **************************************************************************************************
!> \brief Helper routine to report on the sparsity of a single matrix,
!>        for several filtering values
!> \param matrix ...
!> \param unit_nr ...
!> \param title ...
!> \param eps ...
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! **************************************************************************************************
   SUBROUTINE report_matrix_sparsity(matrix, unit_nr, title, eps)
      TYPE(dbcsr_type)                                   :: matrix
      INTEGER                                            :: unit_nr
      CHARACTER(LEN=*)                                   :: title
      REAL(KIND=dp)                                      :: eps

      CHARACTER(len=*), PARAMETER :: routineN = 'report_matrix_sparsity', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: eps_local, occ
      TYPE(dbcsr_type)                                   :: matrix_tmp

      CALL timeset(routineN, handle)
      CALL dbcsr_create(matrix_tmp, template=matrix, name=TRIM(title))
      CALL dbcsr_copy(matrix_tmp, matrix, name=TRIM(title))

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A)') "Sparsity for : "//TRIM(title)
      ENDIF

      eps_local = MAX(eps, 10E-14_dp)
      DO
         IF (eps_local > 1.1_dp) EXIT
         CALL dbcsr_filter(matrix_tmp, eps_local)
         occ = dbcsr_get_occupation(matrix_tmp)
         IF (unit_nr > 0) WRITE (unit_nr, '(T2,F16.12,A3,F16.12)') eps_local, " : ", occ
         eps_local = eps_local*10
      ENDDO

      CALL dbcsr_release(matrix_tmp)

      CALL timestop(handle)

   END SUBROUTINE report_matrix_sparsity

! **************************************************************************************************
!> \brief Compute matrix_w as needed for the forces
!> \param matrix_w ...
!> \param ls_scf_env ...
!> \par History
!>       2010.11 created [Joost VandeVondele]
!> \author Joost VandeVondele
! **************************************************************************************************
   SUBROUTINE calculate_w_matrix(matrix_w, ls_scf_env)
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_w
      TYPE(ls_scf_env_type)                              :: ls_scf_env

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_w_matrix', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin
      REAL(KIND=dp)                                      :: scaling
      TYPE(dbcsr_type)                                   :: matrix_tmp1, matrix_tmp2, matrix_tmp3

      CALL timeset(routineN, handle)

      CALL dbcsr_create(matrix_tmp1, template=ls_scf_env%matrix_s, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(matrix_tmp2, template=ls_scf_env%matrix_s, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(matrix_tmp3, template=ls_scf_env%matrix_s, matrix_type=dbcsr_type_no_symmetry)

      IF (ls_scf_env%nspins == 1) THEN
         scaling = 0.5_dp
      ELSE
         scaling = 1.0_dp
      ENDIF

      DO ispin = 1, ls_scf_env%nspins

         CALL dbcsr_copy(matrix_tmp3, ls_scf_env%matrix_ks(ispin))
         IF (ls_scf_env%has_s_preconditioner) THEN
            CALL apply_matrix_preconditioner(matrix_tmp3, "backward", &
                                             ls_scf_env%matrix_bs_sqrt, ls_scf_env%matrix_bs_sqrt_inv)
         ENDIF
         CALL dbcsr_filter(matrix_tmp3, ls_scf_env%eps_filter)

         CALL dbcsr_multiply("N", "N", scaling, ls_scf_env%matrix_p(ispin), matrix_tmp3, &
                             0.0_dp, matrix_tmp1, filter_eps=ls_scf_env%eps_filter)
         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, ls_scf_env%matrix_p(ispin), &
                             0.0_dp, matrix_tmp2, filter_eps=ls_scf_env%eps_filter)
         CALL matrix_ls_to_qs(matrix_w(ispin)%matrix, matrix_tmp2, ls_scf_env%ls_mstruct, covariant=.FALSE.)
      ENDDO

      CALL dbcsr_release(matrix_tmp1)
      CALL dbcsr_release(matrix_tmp2)
      CALL dbcsr_release(matrix_tmp3)

      CALL timestop(handle)

   END SUBROUTINE calculate_w_matrix

! **************************************************************************************************
!> \brief a place for quick experiments
!> \par History
!>       2010.11 created [Joost VandeVondele]
!> \author Joost VandeVondele
! **************************************************************************************************
   SUBROUTINE post_scf_experiment()

      CHARACTER(len=*), PARAMETER :: routineN = 'post_scf_experiment', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      CALL timestop(handle)
   END SUBROUTINE post_scf_experiment

END MODULE dm_ls_scf
