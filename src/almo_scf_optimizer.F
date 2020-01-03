!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Optimization routines for all ALMO-based SCF methods
!> \par History
!>       2011.05 created [Rustam Z Khaliullin]
!>       2014.10 as a separate file [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
MODULE almo_scf_optimizer
   USE almo_scf_diis_types,             ONLY: almo_scf_diis_extrapolate,&
                                              almo_scf_diis_init,&
                                              almo_scf_diis_push,&
                                              almo_scf_diis_release,&
                                              almo_scf_diis_type
   USE almo_scf_methods,                ONLY: &
        almo_scf_ks_blk_to_tv_blk, almo_scf_ks_to_ks_blk, almo_scf_ks_to_ks_xx, &
        almo_scf_ks_xx_to_tv_xx, almo_scf_p_blk_to_t_blk, almo_scf_t_rescaling, &
        almo_scf_t_to_proj, apply_domain_operators, apply_projector, &
        construct_domain_preconditioner, construct_domain_r_down, construct_domain_s_inv, &
        construct_domain_s_sqrt, get_overlap, orthogonalize_mos, pseudo_invert_diagonal_blk, &
        xalmo_initial_guess
   USE almo_scf_qs,                     ONLY: almo_dm_to_almo_ks,&
                                              almo_dm_to_qs_env,&
                                              almo_scf_update_ks_energy
   USE almo_scf_types,                  ONLY: almo_scf_env_type,&
                                              optimizer_options_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_dbcsr_cholesky,               ONLY: cp_dbcsr_cholesky_decompose,&
                                              cp_dbcsr_cholesky_invert,&
                                              cp_dbcsr_cholesky_restore
   USE cp_external_control,             ONLY: external_control
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE ct_methods,                      ONLY: analytic_line_search,&
                                              ct_step_execute,&
                                              diagonalize_diagonal_blocks
   USE ct_types,                        ONLY: ct_step_env_clean,&
                                              ct_step_env_get,&
                                              ct_step_env_init,&
                                              ct_step_env_set,&
                                              ct_step_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_add_on_diag, dbcsr_copy, dbcsr_create, dbcsr_desymmetrize, &
        dbcsr_distribution_get, dbcsr_distribution_type, dbcsr_dot, dbcsr_filter, dbcsr_finalize, &
        dbcsr_frobenius_norm, dbcsr_func_dtanh, dbcsr_func_inverse, dbcsr_func_tanh, &
        dbcsr_function_of_elements, dbcsr_get_block_p, dbcsr_get_diag, dbcsr_get_info, &
        dbcsr_hadamard_product, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, &
        dbcsr_nblkcols_total, dbcsr_nblkrows_total, dbcsr_norm, dbcsr_norm_maxabsnorm, &
        dbcsr_p_type, dbcsr_print_block_sum, dbcsr_release, dbcsr_reserve_block2d, dbcsr_scale, &
        dbcsr_set, dbcsr_set_diag, dbcsr_triu, dbcsr_type, dbcsr_type_no_symmetry, &
        dbcsr_work_create
   USE domain_submatrix_methods,        ONLY: add_submatrices,&
                                              construct_submatrices,&
                                              copy_submatrices,&
                                              init_submatrices,&
                                              maxnorm_submatrices,&
                                              release_submatrices
   USE domain_submatrix_types,          ONLY: domain_map_type,&
                                              domain_submatrix_type,&
                                              select_row
   USE input_constants,                 ONLY: &
        almo_occ_vol_penalty_none, almo_scf_diag, almo_scf_dm_sign, cg_dai_yuan, cg_fletcher, &
        cg_fletcher_reeves, cg_hager_zhang, cg_hestenes_stiefel, cg_liu_storey, cg_polak_ribiere, &
        cg_zero, virt_full, xalmo_case_block_diag, xalmo_case_fully_deloc, xalmo_case_normal, &
        xalmo_prec_domain, xalmo_prec_full, xalmo_prec_zero
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE iterate_matrix,                  ONLY: determinant,&
                                              invert_Hotelling,&
                                              matrix_sqrt_Newton_Schulz
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'almo_scf_optimizer'

   PUBLIC :: almo_scf_block_diagonal, &
             almo_scf_xalmo_eigensolver, &
             almo_scf_xalmo_pcg

   LOGICAL, PARAMETER :: debug_mode = .FALSE.
   LOGICAL, PARAMETER :: safe_mode = .FALSE.
   LOGICAL, PARAMETER :: almo_mathematica = .FALSE.
   INTEGER, PARAMETER :: hessian_path_reuse = 1, &
                         hessian_path_assemble = 2

CONTAINS

! **************************************************************************************************
!> \brief An SCF procedure that optimizes block-diagonal ALMOs using DIIS
!> \param qs_env ...
!> \param almo_scf_env ...
!> \param optimizer ...
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!>       2018.09 smearing support [Ruben Staub]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE almo_scf_block_diagonal(qs_env, almo_scf_env, optimizer)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(almo_scf_env_type)                            :: almo_scf_env
      TYPE(optimizer_options_type)                       :: optimizer

      CHARACTER(len=*), PARAMETER :: routineN = 'almo_scf_block_diagonal', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iscf, ispin, nspin, unit_nr
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: local_nocc_of_domain
      LOGICAL                                            :: converged, prepare_to_exit, should_stop, &
                                                            use_diis, use_prev_as_guess
      REAL(KIND=dp) :: density_rec, energy_diff, energy_new, energy_old, error_norm, &
         error_norm_ispin, kTS_sum, prev_error_norm, t1, t2, true_mixing_fraction
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: local_mu
      TYPE(almo_scf_diis_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: almo_diis
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:)        :: matrix_mixing_old_blk
      TYPE(qs_energy_type), POINTER                      :: qs_energy

!TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks

      CALL timeset(routineN, handle)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      ! use DIIS, it's superior to simple mixing
      use_diis = .TRUE.
      use_prev_as_guess = .FALSE.

      nspin = almo_scf_env%nspins
      ALLOCATE (local_mu(almo_scf_env%ndomains))
      ALLOCATE (local_nocc_of_domain(almo_scf_env%ndomains))

      ! init mixing matrices
      ALLOCATE (matrix_mixing_old_blk(nspin))
      ALLOCATE (almo_diis(nspin))
      DO ispin = 1, nspin
         CALL dbcsr_create(matrix_mixing_old_blk(ispin), &
                           template=almo_scf_env%matrix_ks_blk(ispin))
         CALL almo_scf_diis_init(diis_env=almo_diis(ispin), &
                                 sample_err=almo_scf_env%matrix_ks_blk(ispin), &
                                 sample_var=almo_scf_env%matrix_s_blk(1), &
                                 error_type=1, &
                                 max_length=optimizer%ndiis)
      ENDDO

      CALL get_qs_env(qs_env, energy=qs_energy)
      energy_old = qs_energy%total

      iscf = 0
      prepare_to_exit = .FALSE.
      true_mixing_fraction = 0.0_dp
      error_norm = 1.0E+10_dp ! arbitrary big step

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A,A,A)') REPEAT("-", 20), &
            " Optimization of block-diagonal ALMOs ", REPEAT("-", 21)
         WRITE (unit_nr, *)
         WRITE (unit_nr, '(T2,A13,A6,A23,A14,A14,A9)') "Method", "Iter", &
            "Total Energy", "Change", "Convergence", "Time"
         WRITE (unit_nr, '(T2,A)') REPEAT("-", 79)
      ENDIF

      ! the real SCF loop
      t1 = m_walltime()
      DO

         iscf = iscf + 1

         ! obtain projected KS matrix and the DIIS-error vector
         CALL almo_scf_ks_to_ks_blk(almo_scf_env)

         ! inform the DIIS handler about the new KS matrix and its error vector
         IF (use_diis) THEN
            DO ispin = 1, nspin
               CALL almo_scf_diis_push(diis_env=almo_diis(ispin), &
                                       var=almo_scf_env%matrix_ks_blk(ispin), &
                                       err=almo_scf_env%matrix_err_blk(ispin))
            ENDDO
         ENDIF

         ! get error_norm: choose the largest of the two spins
         prev_error_norm = error_norm
         DO ispin = 1, nspin
            !error_norm=dbcsr_frobenius_norm(almo_scf_env%matrix_err_blk(ispin))
            CALL dbcsr_norm(almo_scf_env%matrix_err_blk(ispin), &
                            dbcsr_norm_maxabsnorm, &
                            norm_scalar=error_norm_ispin)
            IF (ispin .EQ. 1) error_norm = error_norm_ispin
            IF (ispin .GT. 1 .AND. error_norm_ispin .GT. error_norm) &
               error_norm = error_norm_ispin
         ENDDO

         IF (error_norm .LT. almo_scf_env%eps_prev_guess) THEN
            use_prev_as_guess = .TRUE.
         ELSE
            use_prev_as_guess = .FALSE.
         ENDIF

         ! check convergence
         converged = .TRUE.
         IF (error_norm .GT. optimizer%eps_error) converged = .FALSE.

         ! check other exit criteria: max SCF steps and timing
         CALL external_control(should_stop, "SCF", &
                               start_time=qs_env%start_time, &
                               target_time=qs_env%target_time)
         IF (should_stop .OR. iscf >= optimizer%max_iter .OR. converged) THEN
            prepare_to_exit = .TRUE.
            IF (iscf == 1) energy_new = energy_old
         ENDIF

         ! if early stopping is on do at least one iteration
         IF (optimizer%early_stopping_on .AND. iscf .EQ. 1) &
            prepare_to_exit = .FALSE.

         IF (.NOT. prepare_to_exit) THEN ! update the ALMOs and density matrix

            ! perform mixing of KS matrices
            IF (iscf .NE. 1) THEN
               IF (use_diis) THEN ! use diis instead of mixing
                  DO ispin = 1, nspin
                     CALL almo_scf_diis_extrapolate(diis_env=almo_diis(ispin), &
                                                    extr_var=almo_scf_env%matrix_ks_blk(ispin))
                  ENDDO
               ELSE ! use mixing
                  true_mixing_fraction = almo_scf_env%mixing_fraction
                  DO ispin = 1, nspin
                     CALL dbcsr_add(almo_scf_env%matrix_ks_blk(ispin), &
                                    matrix_mixing_old_blk(ispin), &
                                    true_mixing_fraction, &
                                    1.0_dp - true_mixing_fraction)
                  END DO
               ENDIF
            ENDIF
            ! save the new matrix for the future mixing
            DO ispin = 1, nspin
               CALL dbcsr_copy(matrix_mixing_old_blk(ispin), &
                               almo_scf_env%matrix_ks_blk(ispin))
            ENDDO

            ! obtain ALMOs from the new KS matrix
            SELECT CASE (almo_scf_env%almo_update_algorithm)
            CASE (almo_scf_diag)

               CALL almo_scf_ks_blk_to_tv_blk(almo_scf_env)

            CASE (almo_scf_dm_sign)

               ! update the density matrix
               DO ispin = 1, nspin

                  local_nocc_of_domain(:) = almo_scf_env%nocc_of_domain(:, ispin)
                  local_mu(:) = almo_scf_env%mu_of_domain(:, ispin)
                  ! RZK UPDATE! the update algorithm is removed because
                  ! RZK UPDATE! it requires updating core LS_SCF routines
                  ! RZK UPDATE! (the code exists in the CVS version)
                  CPABORT("Density_matrix_sign has not been tested yet")
                  ! RZK UPDATE!  CALL density_matrix_sign(almo_scf_env%matrix_p_blk(ispin),&
                  ! RZK UPDATE!          local_mu,&
                  ! RZK UPDATE!          almo_scf_env%fixed_mu,&
                  ! RZK UPDATE!          almo_scf_env%matrix_ks_blk(ispin),&
                  ! RZK UPDATE!          !matrix_mixing_old_blk(ispin),&
                  ! RZK UPDATE!          almo_scf_env%matrix_s_blk(1), &
                  ! RZK UPDATE!          almo_scf_env%matrix_s_blk_inv(1), &
                  ! RZK UPDATE!          local_nocc_of_domain,&
                  ! RZK UPDATE!          almo_scf_env%eps_filter,&
                  ! RZK UPDATE!          almo_scf_env%domain_index_of_ao)
                  ! RZK UPDATE!
                  almo_scf_env%mu_of_domain(:, ispin) = local_mu(:)

               ENDDO

               ! obtain ALMOs from matrix_p_blk: T_new = P_blk S_blk T_old
               CALL almo_scf_p_blk_to_t_blk(almo_scf_env, ionic=.FALSE.)

               DO ispin = 1, almo_scf_env%nspins

                  CALL orthogonalize_mos(ket=almo_scf_env%matrix_t_blk(ispin), &
                                         overlap=almo_scf_env%matrix_sigma_blk(ispin), &
                                         metric=almo_scf_env%matrix_s_blk(1), &
                                         retain_locality=.TRUE., &
                                         only_normalize=.FALSE., &
                                         nocc_of_domain=almo_scf_env%nocc_of_domain(:, ispin), &
                                         eps_filter=almo_scf_env%eps_filter, &
                                         order_lanczos=almo_scf_env%order_lanczos, &
                                         eps_lanczos=almo_scf_env%eps_lanczos, &
                                         max_iter_lanczos=almo_scf_env%max_iter_lanczos)

               ENDDO

            END SELECT

            ! obtain density matrix from ALMOs
            DO ispin = 1, almo_scf_env%nspins

               !! Application of an occupation-rescaling trick for smearing, if requested
               IF (almo_scf_env%smear) THEN
                  CALL almo_scf_t_rescaling(matrix_t=almo_scf_env%matrix_t_blk(ispin), &
                                            mo_energies=almo_scf_env%mo_energies(:, ispin), &
                                            mu_of_domain=almo_scf_env%mu_of_domain(:, ispin), &
                                            real_ne_of_domain=almo_scf_env%real_ne_of_domain(:, ispin), &
                                            spin_kTS=almo_scf_env%kTS(ispin), &
                                            smear_e_temp=almo_scf_env%smear_e_temp, &
                                            ndomains=almo_scf_env%ndomains, &
                                            nocc_of_domain=almo_scf_env%nocc_of_domain(:, ispin))
               END IF

               CALL almo_scf_t_to_proj(t=almo_scf_env%matrix_t_blk(ispin), &
                                       p=almo_scf_env%matrix_p(ispin), &
                                       eps_filter=almo_scf_env%eps_filter, &
                                       orthog_orbs=.FALSE., &
                                       nocc_of_domain=almo_scf_env%nocc_of_domain(:, ispin), &
                                       s=almo_scf_env%matrix_s(1), &
                                       sigma=almo_scf_env%matrix_sigma(ispin), &
                                       sigma_inv=almo_scf_env%matrix_sigma_inv(ispin), &
                                       use_guess=use_prev_as_guess, &
                                       smear=almo_scf_env%smear, &
                                       algorithm=almo_scf_env%sigma_inv_algorithm, &
                                       inverse_accelerator=almo_scf_env%order_lanczos, &
                                       inv_eps_factor=almo_scf_env%matrix_iter_eps_error_factor, &
                                       eps_lanczos=almo_scf_env%eps_lanczos, &
                                       max_iter_lanczos=almo_scf_env%max_iter_lanczos, &
                                       para_env=almo_scf_env%para_env, &
                                       blacs_env=almo_scf_env%blacs_env)

            ENDDO

            IF (almo_scf_env%nspins == 1) THEN
               CALL dbcsr_scale(almo_scf_env%matrix_p(1), 2.0_dp)
               !! Rescaling electronic entropy contribution by spin_factor
               IF (almo_scf_env%smear) THEN
                  almo_scf_env%kTS(1) = almo_scf_env%kTS(1)*2.0_dp
               END IF
            ENDIF

            IF (almo_scf_env%smear) THEN
               kTS_sum = SUM(almo_scf_env%kTS)
            ELSE
               kTS_sum = 0.0_dp
            ENDIF

            ! compute the new KS matrix and new energy
            CALL almo_dm_to_almo_ks(qs_env, &
                                    almo_scf_env%matrix_p, &
                                    almo_scf_env%matrix_ks, &
                                    energy_new, &
                                    almo_scf_env%eps_filter, &
                                    almo_scf_env%mat_distr_aos, &
                                    smear=almo_scf_env%smear, &
                                    kTS_sum=kTS_sum)

         ENDIF ! prepare_to_exit

         energy_diff = energy_new - energy_old
         energy_old = energy_new
         almo_scf_env%almo_scf_energy = energy_new

         t2 = m_walltime()
         ! brief report on the current SCF loop
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A13,I6,F23.10,E14.5,F14.9,F9.2)') "ALMO SCF DIIS", &
               iscf, &
               energy_new, energy_diff, error_norm, t2 - t1
         ENDIF
         t1 = m_walltime()

         IF (prepare_to_exit) EXIT

      ENDDO ! end scf cycle

      !! Print number of electrons recovered if smearing was requested
      IF (almo_scf_env%smear) THEN
         DO ispin = 1, nspin
            CALL dbcsr_dot(almo_scf_env%matrix_p(ispin), almo_scf_env%matrix_s(1), density_rec)
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(T2,A20,F23.10)') "Electrons recovered:", density_rec
            END IF
         END DO
      END IF

      IF (.NOT. converged .AND. (.NOT. optimizer%early_stopping_on)) THEN
         IF (unit_nr > 0) THEN
            CPABORT("SCF for block-diagonal ALMOs not converged!")
         ENDIF
      ENDIF

      DO ispin = 1, nspin
         CALL dbcsr_release(matrix_mixing_old_blk(ispin))
         CALL almo_scf_diis_release(diis_env=almo_diis(ispin))
      ENDDO
      DEALLOCATE (almo_diis)
      DEALLOCATE (matrix_mixing_old_blk)
      DEALLOCATE (local_mu)
      DEALLOCATE (local_nocc_of_domain)

      CALL timestop(handle)

   END SUBROUTINE almo_scf_block_diagonal

! **************************************************************************************************
!> \brief An eigensolver-based SCF to optimize extended ALMOs (i.e. ALMOs on
!>        overlapping domains)
!> \param qs_env ...
!> \param almo_scf_env ...
!> \param optimizer ...
!> \par History
!>       2013.03 created [Rustam Z Khaliullin]
!>       2018.09 smearing support [Ruben Staub]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE almo_scf_xalmo_eigensolver(qs_env, almo_scf_env, optimizer)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(almo_scf_env_type)                            :: almo_scf_env
      TYPE(optimizer_options_type)                       :: optimizer

      CHARACTER(len=*), PARAMETER :: routineN = 'almo_scf_xalmo_eigensolver', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iscf, ispin, nspin, unit_nr
      LOGICAL                                            :: converged, prepare_to_exit, should_stop
      REAL(KIND=dp) :: denergy_tot, density_rec, energy_diff, energy_new, energy_old, error_norm, &
         error_norm_0, kTS_sum, spin_factor, t1, t2
      REAL(KIND=dp), DIMENSION(2)                        :: denergy_spin
      TYPE(almo_scf_diis_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: almo_diis
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type)                                   :: matrix_p_almo_scf_converged
      TYPE(domain_submatrix_type), ALLOCATABLE, &
         DIMENSION(:, :)                                 :: submatrix_mixing_old_blk

      CALL timeset(routineN, handle)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      nspin = almo_scf_env%nspins
      IF (nspin == 1) THEN
         spin_factor = 2.0_dp
      ELSE
         spin_factor = 1.0_dp
      ENDIF

      ! RZK-warning domain_s_sqrt and domain_s_sqrt_inv do not have spin
      ! components yet (may be used later)
      ispin = 1
      CALL construct_domain_s_sqrt( &
         matrix_s=almo_scf_env%matrix_s(1), &
         subm_s_sqrt=almo_scf_env%domain_s_sqrt(:, ispin), &
         subm_s_sqrt_inv=almo_scf_env%domain_s_sqrt_inv(:, ispin), &
         dpattern=almo_scf_env%quench_t(ispin), &
         map=almo_scf_env%domain_map(ispin), &
         node_of_domain=almo_scf_env%cpu_of_domain)
      ! TRY: construct s_inv
      !CALL construct_domain_s_inv(&
      !       matrix_s=almo_scf_env%matrix_s(1),&
      !       subm_s_inv=almo_scf_env%domain_s_inv(:,ispin),&
      !       dpattern=almo_scf_env%quench_t(ispin),&
      !       map=almo_scf_env%domain_map(ispin),&
      !       node_of_domain=almo_scf_env%cpu_of_domain)

      ! construct the domain template for the occupied orbitals
      DO ispin = 1, nspin
         ! RZK-warning we need only the matrix structure, not data
         ! replace construct_submatrices with lighter procedure with
         ! no heavy communications
         CALL construct_submatrices( &
            matrix=almo_scf_env%quench_t(ispin), &
            submatrix=almo_scf_env%domain_t(:, ispin), &
            distr_pattern=almo_scf_env%quench_t(ispin), &
            domain_map=almo_scf_env%domain_map(ispin), &
            node_of_domain=almo_scf_env%cpu_of_domain, &
            job_type=select_row)
      ENDDO

      ! init mixing matrices
      ALLOCATE (submatrix_mixing_old_blk(almo_scf_env%ndomains, nspin))
      CALL init_submatrices(submatrix_mixing_old_blk)
      ALLOCATE (almo_diis(nspin))

      ! TRY: construct block-projector
      !ALLOCATE(submatrix_tmp(almo_scf_env%ndomains))
      !DO ispin=1,nspin
      !   CALL init_submatrices(submatrix_tmp)
      !   CALL construct_domain_r_down(&
      !           matrix_t=almo_scf_env%matrix_t_blk(ispin),&
      !           matrix_sigma_inv=almo_scf_env%matrix_sigma_inv(ispin),&
      !           matrix_s=almo_scf_env%matrix_s(1),&
      !           subm_r_down=submatrix_tmp(:),&
      !           dpattern=almo_scf_env%quench_t(ispin),&
      !           map=almo_scf_env%domain_map(ispin),&
      !           node_of_domain=almo_scf_env%cpu_of_domain,&
      !           filter_eps=almo_scf_env%eps_filter)
      !   CALL multiply_submatrices('N','N',1.0_dp,&
      !           submatrix_tmp(:),&
      !           almo_scf_env%domain_s_inv(:,1),0.0_dp,&
      !           almo_scf_env%domain_r_down_up(:,ispin))
      !   CALL release_submatrices(submatrix_tmp)
      !ENDDO
      !DEALLOCATE(submatrix_tmp)

      DO ispin = 1, nspin
         ! use s_sqrt since they are already properly constructed
         ! and have the same distributions as domain_err and domain_ks_xx
         CALL almo_scf_diis_init(diis_env=almo_diis(ispin), &
                                 sample_err=almo_scf_env%domain_s_sqrt(:, ispin), &
                                 error_type=1, &
                                 max_length=optimizer%ndiis)
      ENDDO

      denergy_tot = 0.0_dp
      energy_old = 0.0_dp
      iscf = 0
      prepare_to_exit = .FALSE.

      ! the SCF loop
      t1 = m_walltime()
      DO

         iscf = iscf + 1

         ! obtain projected KS matrix and the DIIS-error vector
         CALL almo_scf_ks_to_ks_xx(almo_scf_env)

         ! inform the DIIS handler about the new KS matrix and its error vector
         DO ispin = 1, nspin
            CALL almo_scf_diis_push(diis_env=almo_diis(ispin), &
                                    d_var=almo_scf_env%domain_ks_xx(:, ispin), &
                                    d_err=almo_scf_env%domain_err(:, ispin))
         ENDDO

         ! check convergence
         converged = .TRUE.
         DO ispin = 1, nspin
            !error_norm=dbcsr_frobenius_norm(almo_scf_env%matrix_err_blk(ispin))
            CALL dbcsr_norm(almo_scf_env%matrix_err_xx(ispin), &
                            dbcsr_norm_maxabsnorm, &
                            norm_scalar=error_norm)
            CALL maxnorm_submatrices(almo_scf_env%domain_err(:, ispin), &
                                     norm=error_norm_0)
            IF (error_norm .GT. optimizer%eps_error) THEN
               converged = .FALSE.
               EXIT ! no need to check the other spin
            ENDIF
         ENDDO
         ! check other exit criteria: max SCF steps and timing
         CALL external_control(should_stop, "SCF", &
                               start_time=qs_env%start_time, &
                               target_time=qs_env%target_time)
         IF (should_stop .OR. iscf >= optimizer%max_iter .OR. converged) THEN
            prepare_to_exit = .TRUE.
         ENDIF

         ! if early stopping is on do at least one iteration
         IF (optimizer%early_stopping_on .AND. iscf .EQ. 1) &
            prepare_to_exit = .FALSE.

         IF (.NOT. prepare_to_exit) THEN ! update the ALMOs and density matrix

            ! perform mixing of KS matrices
            IF (iscf .NE. 1) THEN
               IF (.FALSE.) THEN ! use diis instead of mixing
                  DO ispin = 1, nspin
                     CALL add_submatrices( &
                        almo_scf_env%mixing_fraction, &
                        almo_scf_env%domain_ks_xx(:, ispin), &
                        1.0_dp - almo_scf_env%mixing_fraction, &
                        submatrix_mixing_old_blk(:, ispin), &
                        'N')
                  END DO
               ELSE
                  DO ispin = 1, nspin
                     CALL almo_scf_diis_extrapolate(diis_env=almo_diis(ispin), &
                                                    d_extr_var=almo_scf_env%domain_ks_xx(:, ispin))
                  ENDDO
               ENDIF
            ENDIF
            ! save the new matrix for the future mixing
            DO ispin = 1, nspin
               CALL copy_submatrices( &
                  almo_scf_env%domain_ks_xx(:, ispin), &
                  submatrix_mixing_old_blk(:, ispin), &
                  copy_data=.TRUE.)
            ENDDO

            ! obtain a new set of ALMOs from the updated KS matrix
            CALL almo_scf_ks_xx_to_tv_xx(almo_scf_env)

            ! update the density matrix
            DO ispin = 1, nspin

               ! save the initial density matrix (to get the perturbative energy lowering)
               IF (iscf .EQ. 1) THEN
                  CALL dbcsr_create(matrix_p_almo_scf_converged, &
                                    template=almo_scf_env%matrix_p(ispin))
                  CALL dbcsr_copy(matrix_p_almo_scf_converged, &
                                  almo_scf_env%matrix_p(ispin))
               ENDIF

               !! Application of an occupation-rescaling trick for smearing, if requested
               IF (almo_scf_env%smear) THEN
                  CALL almo_scf_t_rescaling(matrix_t=almo_scf_env%matrix_t_blk(ispin), &
                                            mo_energies=almo_scf_env%mo_energies(:, ispin), &
                                            mu_of_domain=almo_scf_env%mu_of_domain(:, ispin), &
                                            real_ne_of_domain=almo_scf_env%real_ne_of_domain(:, ispin), &
                                            spin_kTS=almo_scf_env%kTS(ispin), &
                                            smear_e_temp=almo_scf_env%smear_e_temp, &
                                            ndomains=almo_scf_env%ndomains, &
                                            nocc_of_domain=almo_scf_env%nocc_of_domain(:, ispin))
               END IF

               ! update now
               CALL almo_scf_t_to_proj( &
                  t=almo_scf_env%matrix_t(ispin), &
                  p=almo_scf_env%matrix_p(ispin), &
                  eps_filter=almo_scf_env%eps_filter, &
                  orthog_orbs=.FALSE., &
                  nocc_of_domain=almo_scf_env%nocc_of_domain(:, ispin), &
                  s=almo_scf_env%matrix_s(1), &
                  sigma=almo_scf_env%matrix_sigma(ispin), &
                  sigma_inv=almo_scf_env%matrix_sigma_inv(ispin), &
                  use_guess=.TRUE., &
                  smear=almo_scf_env%smear, &
                  algorithm=almo_scf_env%sigma_inv_algorithm, &
                  inverse_accelerator=almo_scf_env%order_lanczos, &
                  inv_eps_factor=almo_scf_env%matrix_iter_eps_error_factor, &
                  eps_lanczos=almo_scf_env%eps_lanczos, &
                  max_iter_lanczos=almo_scf_env%max_iter_lanczos, &
                  para_env=almo_scf_env%para_env, &
                  blacs_env=almo_scf_env%blacs_env)
               CALL dbcsr_scale(almo_scf_env%matrix_p(ispin), spin_factor)
               !! Rescaling electronic entropy contribution by spin_factor
               IF (almo_scf_env%smear) THEN
                  almo_scf_env%kTS(ispin) = almo_scf_env%kTS(ispin)*spin_factor
               END IF

               ! obtain perturbative estimate (at no additional cost)
               ! of the energy lowering relative to the block-diagonal ALMOs
               IF (iscf .EQ. 1) THEN

                  CALL dbcsr_add(matrix_p_almo_scf_converged, &
                                 almo_scf_env%matrix_p(ispin), -1.0_dp, 1.0_dp)
                  CALL dbcsr_dot(almo_scf_env%matrix_ks_0deloc(ispin), &
                                 matrix_p_almo_scf_converged, &
                                 denergy_spin(ispin))

                  CALL dbcsr_release(matrix_p_almo_scf_converged)

                  !! RS-WARNING: If smearing ALMO is requested, electronic entropy contribution should probably be included here

                  denergy_tot = denergy_tot + denergy_spin(ispin)

                  ! RZK-warning Energy correction can be evaluated using matrix_x
                  ! as shown in the attempt below and in the PCG procedure.
                  ! Using matrix_x allows immediate decomposition of the energy
                  ! lowering into 2-body components for EDA. However, it does not
                  ! work here because the diagonalization routine does not necessarily
                  ! produce orbitals with the same sign as the block-diagonal ALMOs
                  ! Any fixes?!

                  !CALL dbcsr_init(matrix_x)
                  !CALL dbcsr_create(matrix_x,&
                  !        template=almo_scf_env%matrix_t(ispin))
                  !
                  !CALL dbcsr_init(matrix_tmp_no)
                  !CALL dbcsr_create(matrix_tmp_no,&
                  !        template=almo_scf_env%matrix_t(ispin))
                  !
                  !CALL dbcsr_copy(matrix_x,&
                  !        almo_scf_env%matrix_t_blk(ispin))
                  !CALL dbcsr_add(matrix_x,almo_scf_env%matrix_t(ispin),&
                  !        -1.0_dp,1.0_dp)

                  !CALL dbcsr_dot(matrix_x, almo_scf_env%matrix_err_xx(ispin),denergy)

                  !denergy=denergy*spin_factor

                  !IF (unit_nr>0) THEN
                  !   WRITE(unit_nr,*) "_ENERGY-0: ", almo_scf_env%almo_scf_energy
                  !   WRITE(unit_nr,*) "_ENERGY-D: ", denergy
                  !   WRITE(unit_nr,*) "_ENERGY-F: ", almo_scf_env%almo_scf_energy+denergy
                  !ENDIF
                  !! RZK-warning update will not work since the energy is overwritten almost immediately
                  !!CALL almo_scf_update_ks_energy(qs_env,&
                  !!        almo_scf_env%almo_scf_energy+denergy)
                  !!

                  !! print out the results of the decomposition analysis
                  !CALL dbcsr_hadamard_product(matrix_x,&
                  !        almo_scf_env%matrix_err_xx(ispin),&
                  !        matrix_tmp_no)
                  !CALL dbcsr_scale(matrix_tmp_no,spin_factor)
                  !CALL dbcsr_filter(matrix_tmp_no,almo_scf_env%eps_filter)
                  !
                  !IF (unit_nr>0) THEN
                  !   WRITE(unit_nr,*)
                  !   WRITE(unit_nr,'(T2,A)') "DECOMPOSITION OF THE DELOCALIZATION ENERGY"
                  !ENDIF

                  !mynode=dbcsr_mp_mynode(dbcsr_distribution_mp(&
                  !   dbcsr_distribution(matrix_tmp_no)))
                  !WRITE(mynodestr,'(I6.6)') mynode
                  !mylogfile='EDA.'//TRIM(ADJUSTL(mynodestr))
                  !OPEN (iunit,file=mylogfile,status='REPLACE')
                  !CALL dbcsr_print_block_sum(matrix_tmp_no,iunit)
                  !CLOSE(iunit)
                  !
                  !CALL dbcsr_release(matrix_tmp_no)
                  !CALL dbcsr_release(matrix_x)

               ENDIF ! iscf.eq.1

            ENDDO

            ! print out the energy lowering
            IF (iscf .EQ. 1) THEN
               IF (unit_nr > 0) THEN
                  WRITE (unit_nr, *)
                  WRITE (unit_nr, '(T2,A35,F25.10)') "ENERGY OF BLOCK-DIAGONAL ALMOs:", &
                     almo_scf_env%almo_scf_energy
                  WRITE (unit_nr, '(T2,A35,F25.10)') "ENERGY LOWERING:", &
                     denergy_tot
                  WRITE (unit_nr, '(T2,A35,F25.10)') "CORRECTED ENERGY:", &
                     almo_scf_env%almo_scf_energy + denergy_tot
                  WRITE (unit_nr, *)
               ENDIF
               CALL almo_scf_update_ks_energy(qs_env, &
                                              energy=almo_scf_env%almo_scf_energy, &
                                              energy_singles_corr=denergy_tot)
            ENDIF

            ! compute the new KS matrix and new energy
            IF (.NOT. almo_scf_env%perturbative_delocalization) THEN

               IF (almo_scf_env%smear) THEN
                  kTS_sum = SUM(almo_scf_env%kTS)
               ELSE
                  kTS_sum = 0.0_dp
               ENDIF

               CALL almo_dm_to_almo_ks(qs_env, &
                                       almo_scf_env%matrix_p, &
                                       almo_scf_env%matrix_ks, &
                                       energy_new, &
                                       almo_scf_env%eps_filter, &
                                       almo_scf_env%mat_distr_aos, &
                                       smear=almo_scf_env%smear, &
                                       kTS_sum=kTS_sum)
            ENDIF

         ENDIF ! prepare_to_exit

         IF (almo_scf_env%perturbative_delocalization) THEN

            ! exit after the first step if we do not need the SCF procedure
            CALL almo_dm_to_qs_env(qs_env, almo_scf_env%matrix_p, almo_scf_env%mat_distr_aos)
            converged = .TRUE.
            prepare_to_exit = .TRUE.

         ELSE ! not a perturbative treatment

            energy_diff = energy_new - energy_old
            energy_old = energy_new
            almo_scf_env%almo_scf_energy = energy_new

            t2 = m_walltime()
            ! brief report on the current SCF loop
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(T2,A,I6,F20.9,E11.3,E11.3,E11.3,F8.2)') "ALMO SCF", &
                  iscf, &
                  energy_new, energy_diff, error_norm, error_norm_0, t2 - t1
            ENDIF
            t1 = m_walltime()

         ENDIF

         IF (prepare_to_exit) EXIT

      ENDDO ! end scf cycle

      !! Print number of electrons recovered if smearing was requested
      IF (almo_scf_env%smear) THEN
         DO ispin = 1, nspin
            CALL dbcsr_dot(almo_scf_env%matrix_p(ispin), almo_scf_env%matrix_s(1), density_rec)
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(T2,A20,F23.10)') "Electrons recovered:", density_rec
            END IF
         END DO
      END IF

      IF (.NOT. converged .AND. .NOT. optimizer%early_stopping_on) THEN
         CPABORT("SCF for ALMOs on overlapping domains not converged! ")
      ENDIF

      DO ispin = 1, nspin
         CALL release_submatrices(submatrix_mixing_old_blk(:, ispin))
         CALL almo_scf_diis_release(diis_env=almo_diis(ispin))
      ENDDO
      DEALLOCATE (almo_diis)
      DEALLOCATE (submatrix_mixing_old_blk)

      CALL timestop(handle)

   END SUBROUTINE almo_scf_xalmo_eigensolver

! **************************************************************************************************
!> \brief Optimization of ALMOs using PCG-like minimizers
!> \param qs_env ...
!> \param almo_scf_env ...
!> \param optimizer   controls the optimization algorithm
!> \param quench_t ...
!> \param matrix_t_in ...
!> \param matrix_t_out ...
!> \param assume_t0_q0x - since it is extrememly difficult to converge the iterative
!>                        procedure using T as an optimized variable, assume
!>                        T = T_0 + (1-R_0)*X and optimize X
!>                        T_0 is assumed to be the zero-delocalization reference
!> \param perturbation_only - perturbative (do not update Hamiltonian)
!> \param special_case   to reduce the overhead special cases are implemented:
!>                       xalmo_case_normal - no special case (i.e. xALMOs)
!>                       xalmo_case_block_diag
!>                       xalmo_case_fully_deloc
!> \par History
!>       2011.11 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE almo_scf_xalmo_pcg(qs_env, almo_scf_env, optimizer, quench_t, &
                                 matrix_t_in, matrix_t_out, assume_t0_q0x, perturbation_only, &
                                 special_case)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(almo_scf_env_type)                            :: almo_scf_env
      TYPE(optimizer_options_type)                       :: optimizer
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:)        :: quench_t, matrix_t_in, matrix_t_out
      LOGICAL, INTENT(IN)                                :: assume_t0_q0x, perturbation_only
      INTEGER, INTENT(IN), OPTIONAL                      :: special_case

      CHARACTER(len=*), PARAMETER :: routineN = 'almo_scf_xalmo_pcg', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=20)                                  :: iter_type
      INTEGER :: cg_iteration, eda_unit, fixed_line_search_niter, handle, ispin, iteration, &
         line_search_iteration, max_iter, my_special_case, ndomains, nspins, outer_iteration, &
         outer_max_iter, prec_type, unit_nr
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nocc
      LOGICAL :: blissful_neglect, converged, just_started, line_search, normalize_orbitals, &
         optimize_theta, outer_prepare_to_exit, penalty_occ_vol, prepare_to_exit, &
         reset_conjugator, skip_grad, use_guess
      REAL(kind=dp) :: appr_sec_der, beta, denom, denom2, det1, e0, e1, energy_diff, energy_ispin, &
         energy_new, energy_old, eps_skip_gradients, g0, g1, grad_norm, grad_norm_frob, &
         line_search_error, next_step_size_guess, penalty_amplitude, penalty_func_new, &
         spin_factor, step_size, t1, t2, tempreal
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: grad_norm_spin, &
                                                            penalty_occ_vol_g_prefactor, &
                                                            penalty_occ_vol_h_prefactor
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:) :: FTsiginv, grad, m_sig_sqrti_ii, m_t_in_local, &
         m_theta, prec_vv, prev_grad, prev_minus_prec_grad, prev_step, siginvTFTsiginv, ST, step, &
         STsiginv_0
      TYPE(domain_submatrix_type), ALLOCATABLE, &
         DIMENSION(:, :)                                 :: bad_modes_projector_down, domain_r_down
      TYPE(section_vals_type), POINTER                   :: almo_print_section, input

      CALL timeset(routineN, handle)

      my_special_case = xalmo_case_normal
      IF (PRESENT(special_case)) my_special_case = special_case

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      nspins = almo_scf_env%nspins

      ! if unprojected XALMOs are optimized
      ! then we must use the "blissful_neglect" procedure
      blissful_neglect = .FALSE.
      IF (my_special_case .EQ. xalmo_case_normal .AND. .NOT. assume_t0_q0x) THEN
         blissful_neglect = .TRUE.
      ENDIF

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *)
         SELECT CASE (my_special_case)
         CASE (xalmo_case_block_diag)
            WRITE (unit_nr, '(T2,A,A,A)') REPEAT("-", 20), &
               " Optimization of block-diagonal ALMOs ", REPEAT("-", 21)
         CASE (xalmo_case_fully_deloc)
            WRITE (unit_nr, '(T2,A,A,A)') REPEAT("-", 20), &
               " Optimization of fully delocalized MOs ", REPEAT("-", 20)
         CASE (xalmo_case_normal)
            IF (blissful_neglect) THEN
               WRITE (unit_nr, '(T2,A,A,A)') REPEAT("-", 25), &
                  " LCP optimization of XALMOs ", REPEAT("-", 26)
            ELSE
               WRITE (unit_nr, '(T2,A,A,A)') REPEAT("-", 27), &
                  " Optimization of XALMOs ", REPEAT("-", 28)
            ENDIF
         END SELECT
         WRITE (unit_nr, *)
         WRITE (unit_nr, '(T2,A13,A6,A23,A14,A14,A9)') "Method", "Iter", &
            "Objective Function", "Change", "Convergence", "Time"
         WRITE (unit_nr, '(T2,A)') REPEAT("-", 79)
      ENDIF

      ! set local parameters using developer's keywords
      ! RZK-warning: change to normal keywords later
      optimize_theta = almo_scf_env%logical05
      eps_skip_gradients = almo_scf_env%real01

      ! penalty amplitude adjusts the strenght of volume conservation
      ! the following guidelines are useful
      ! A = T for n = 2
      ! A = 2T for n = 4
      ! A = (32/6)T for n = 6
      penalty_occ_vol = (almo_scf_env%penalty%occ_vol_method .NE. almo_occ_vol_penalty_none .AND. &
                         my_special_case .EQ. xalmo_case_fully_deloc)
      normalize_orbitals = penalty_occ_vol
      ! not used with lndet: penalty_order=2
      penalty_amplitude = almo_scf_env%penalty%occ_vol_coeff
      ALLOCATE (penalty_occ_vol_g_prefactor(nspins))
      ALLOCATE (penalty_occ_vol_h_prefactor(nspins))
      penalty_occ_vol_g_prefactor(:) = 0.0_dp
      penalty_occ_vol_h_prefactor(:) = 0.0_dp
      penalty_func_new = 0.0_dp

      ! preconditioner control
      prec_type = optimizer%preconditioner

      ! control of the line search
      fixed_line_search_niter = 0 ! init to zero, change when eps is small enough

      IF (nspins == 1) THEN
         spin_factor = 2.0_dp
      ELSE
         spin_factor = 1.0_dp
      ENDIF

      ALLOCATE (grad_norm_spin(nspins))
      ALLOCATE (nocc(nspins))

      ! create a local copy of matrix_t_in because
      ! matrix_t_in and matrix_t_out can be the same matrix
      ! we need to make sure data in matrix_t_in is intact
      ! after we start writing to matrix_t_out
      ALLOCATE (m_t_in_local(nspins))
      DO ispin = 1, nspins
         CALL dbcsr_create(m_t_in_local(ispin), &
                           template=matrix_t_in(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_copy(m_t_in_local(ispin), matrix_t_in(ispin))
      ENDDO

      ! m_theta contains a set of variational parameters
      ! that define one-electron orbitals (simple, projected, etc.)
      ALLOCATE (m_theta(nspins))
      DO ispin = 1, nspins
         CALL dbcsr_create(m_theta(ispin), &
                           template=matrix_t_out(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
      ENDDO

      ! create initial guess from the initial orbitals
      CALL xalmo_initial_guess(m_guess=m_theta, &
                               m_t_in=m_t_in_local, &
                               m_t0=almo_scf_env%matrix_t_blk, &
                               m_quench_t=quench_t, &
                               m_overlap=almo_scf_env%matrix_s(1), &
                               m_sigma_tmpl=almo_scf_env%matrix_sigma_inv, &
                               nspins=nspins, &
                               xalmo_history=almo_scf_env%xalmo_history, &
                               assume_t0_q0x=assume_t0_q0x, &
                               optimize_theta=optimize_theta, &
                               envelope_amplitude=almo_scf_env%envelope_amplitude, &
                               eps_filter=almo_scf_env%eps_filter, &
                               order_lanczos=almo_scf_env%order_lanczos, &
                               eps_lanczos=almo_scf_env%eps_lanczos, &
                               max_iter_lanczos=almo_scf_env%max_iter_lanczos, &
                               nocc_of_domain=almo_scf_env%nocc_of_domain)

      ndomains = almo_scf_env%ndomains
      ALLOCATE (domain_r_down(ndomains, nspins))
      CALL init_submatrices(domain_r_down)
      ALLOCATE (bad_modes_projector_down(ndomains, nspins))
      CALL init_submatrices(bad_modes_projector_down)

      ALLOCATE (prec_vv(nspins))
      ALLOCATE (siginvTFTsiginv(nspins))
      ALLOCATE (STsiginv_0(nspins))
      ALLOCATE (FTsiginv(nspins))
      ALLOCATE (ST(nspins))
      ALLOCATE (prev_grad(nspins))
      ALLOCATE (grad(nspins))
      ALLOCATE (prev_step(nspins))
      ALLOCATE (step(nspins))
      ALLOCATE (prev_minus_prec_grad(nspins))
      ALLOCATE (m_sig_sqrti_ii(nspins))
      DO ispin = 1, nspins

         ! init temporary storage
         CALL dbcsr_create(prec_vv(ispin), &
                           template=almo_scf_env%matrix_ks(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(siginvTFTsiginv(ispin), &
                           template=almo_scf_env%matrix_sigma(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(STsiginv_0(ispin), &
                           template=matrix_t_out(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(FTsiginv(ispin), &
                           template=matrix_t_out(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(ST(ispin), &
                           template=matrix_t_out(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(prev_grad(ispin), &
                           template=matrix_t_out(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(grad(ispin), &
                           template=matrix_t_out(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(prev_step(ispin), &
                           template=matrix_t_out(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(step(ispin), &
                           template=matrix_t_out(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(prev_minus_prec_grad(ispin), &
                           template=matrix_t_out(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(m_sig_sqrti_ii(ispin), &
                           template=almo_scf_env%matrix_sigma_inv(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)

         CALL dbcsr_set(step(ispin), 0.0_dp)
         CALL dbcsr_set(prev_step(ispin), 0.0_dp)

         CALL dbcsr_get_info(almo_scf_env%matrix_sigma_inv(ispin), &
                             nfullrows_total=nocc(ispin))

         ! invert S domains if necessary
         ! Note: domains for alpha and beta electrons might be different
         ! that is why the inversion of the AO overlap is inside the spin loop
         IF (my_special_case .EQ. xalmo_case_normal) THEN
            CALL construct_domain_s_inv( &
               matrix_s=almo_scf_env%matrix_s(1), &
               subm_s_inv=almo_scf_env%domain_s_inv(:, ispin), &
               dpattern=quench_t(ispin), &
               map=almo_scf_env%domain_map(ispin), &
               node_of_domain=almo_scf_env%cpu_of_domain)

            CALL construct_domain_s_sqrt( &
               matrix_s=almo_scf_env%matrix_s(1), &
               subm_s_sqrt=almo_scf_env%domain_s_sqrt(:, ispin), &
               subm_s_sqrt_inv=almo_scf_env%domain_s_sqrt_inv(:, ispin), &
               dpattern=almo_scf_env%quench_t(ispin), &
               map=almo_scf_env%domain_map(ispin), &
               node_of_domain=almo_scf_env%cpu_of_domain)

         ENDIF

         IF (assume_t0_q0x) THEN

            ! save S.T_0.siginv_0
            IF (my_special_case .EQ. xalmo_case_fully_deloc) THEN
               CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                   almo_scf_env%matrix_s(1), &
                                   almo_scf_env%matrix_t_blk(ispin), &
                                   0.0_dp, ST(ispin), &
                                   filter_eps=almo_scf_env%eps_filter)
               CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                   ST(ispin), &
                                   almo_scf_env%matrix_sigma_inv_0deloc(ispin), &
                                   0.0_dp, STsiginv_0(ispin), &
                                   filter_eps=almo_scf_env%eps_filter)
            ENDIF

            ! construct domain-projector
            IF (my_special_case .EQ. xalmo_case_normal) THEN
               CALL construct_domain_r_down( &
                  matrix_t=almo_scf_env%matrix_t_blk(ispin), &
                  matrix_sigma_inv=almo_scf_env%matrix_sigma_inv(ispin), &
                  matrix_s=almo_scf_env%matrix_s(1), &
                  subm_r_down=domain_r_down(:, ispin), &
                  dpattern=quench_t(ispin), &
                  map=almo_scf_env%domain_map(ispin), &
                  node_of_domain=almo_scf_env%cpu_of_domain, &
                  filter_eps=almo_scf_env%eps_filter)
            ENDIF

         ENDIF ! assume_t0_q0x

      ENDDO ! ispin

      ! start the outer SCF loop
      outer_max_iter = optimizer%max_iter_outer_loop
      outer_prepare_to_exit = .FALSE.
      outer_iteration = 0
      grad_norm = 0.0_dp
      grad_norm_frob = 0.0_dp
      use_guess = .FALSE.

      DO

         ! start the inner SCF loop
         max_iter = optimizer%max_iter
         prepare_to_exit = .FALSE.
         line_search = .FALSE.
         converged = .FALSE.
         iteration = 0
         cg_iteration = 0
         line_search_iteration = 0
         energy_new = 0.0_dp
         energy_old = 0.0_dp
         line_search_error = 0.0_dp

         t1 = m_walltime()

         DO

            just_started = (iteration .EQ. 0) .AND. (outer_iteration .EQ. 0)

            DO ispin = 1, nspins

               ! compute MO coefficients from the main variable
               CALL compute_xalmos_from_main_var( &
                  m_var_in=m_theta(ispin), &
                  m_t_out=matrix_t_out(ispin), &
                  m_quench_t=quench_t(ispin), &
                  m_t0=almo_scf_env%matrix_t_blk(ispin), &
                  m_siginv=almo_scf_env%matrix_sigma_inv(ispin), &
                  m_STsiginv0=STsiginv_0(ispin), &
                  m_s=almo_scf_env%matrix_s(1), &
                  m_sig_sqrti_ii_out=m_sig_sqrti_ii(ispin), &
                  domain_r_down=domain_r_down(:, ispin), &
                  domain_s_inv=almo_scf_env%domain_s_inv(:, ispin), &
                  domain_map=almo_scf_env%domain_map(ispin), &
                  cpu_of_domain=almo_scf_env%cpu_of_domain, &
                  assume_t0_q0x=assume_t0_q0x, &
                  just_started=just_started, &
                  optimize_theta=optimize_theta, &
                  normalize_orbitals=normalize_orbitals, &
                  envelope_amplitude=almo_scf_env%envelope_amplitude, &
                  eps_filter=almo_scf_env%eps_filter, &
                  special_case=my_special_case, &
                  nocc_of_domain=almo_scf_env%nocc_of_domain(:, ispin), &
                  order_lanczos=almo_scf_env%order_lanczos, &
                  eps_lanczos=almo_scf_env%eps_lanczos, &
                  max_iter_lanczos=almo_scf_env%max_iter_lanczos)

               ! compute the global projectors (for the density matrix)
               CALL almo_scf_t_to_proj( &
                  t=matrix_t_out(ispin), &
                  p=almo_scf_env%matrix_p(ispin), &
                  eps_filter=almo_scf_env%eps_filter, &
                  orthog_orbs=.FALSE., &
                  nocc_of_domain=almo_scf_env%nocc_of_domain(:, ispin), &
                  s=almo_scf_env%matrix_s(1), &
                  sigma=almo_scf_env%matrix_sigma(ispin), &
                  sigma_inv=almo_scf_env%matrix_sigma_inv(ispin), &
                  use_guess=use_guess, &
                  algorithm=almo_scf_env%sigma_inv_algorithm, &
                  inv_eps_factor=almo_scf_env%matrix_iter_eps_error_factor, &
                  inverse_accelerator=almo_scf_env%order_lanczos, &
                  eps_lanczos=almo_scf_env%eps_lanczos, &
                  max_iter_lanczos=almo_scf_env%max_iter_lanczos, &
                  para_env=almo_scf_env%para_env, &
                  blacs_env=almo_scf_env%blacs_env)

               ! compute dm from the projector(s)
               CALL dbcsr_scale(almo_scf_env%matrix_p(ispin), &
                                spin_factor)

            ENDDO ! ispin

            ! update the KS matrix and energy if necessary
            IF (perturbation_only) THEN
               ! note: do not combine with the pert_only statement
               IF (just_started) THEN
                  DO ispin = 1, nspins
                     CALL dbcsr_copy(almo_scf_env%matrix_ks(ispin), &
                                     almo_scf_env%matrix_ks_0deloc(ispin))
                  ENDDO
               ENDIF
            ELSE
               ! the KS matrix is updated outside the spin loop
               CALL almo_dm_to_almo_ks(qs_env, &
                                       almo_scf_env%matrix_p, &
                                       almo_scf_env%matrix_ks, &
                                       energy_new, &
                                       almo_scf_env%eps_filter, &
                                       almo_scf_env%mat_distr_aos)
            ENDIF

            !!!! RZK-warning: put the following into a subr - objective function
            DO ispin = 1, nspins

               CALL compute_frequently_used_matrices( &
                  filter_eps=almo_scf_env%eps_filter, &
                  m_T_in=matrix_t_out(ispin), &
                  m_siginv_in=almo_scf_env%matrix_sigma_inv(ispin), &
                  m_S_in=almo_scf_env%matrix_s(1), &
                  m_F_in=almo_scf_env%matrix_ks(ispin), &
                  m_FTsiginv_out=FTsiginv(ispin), &
                  m_siginvTFTsiginv_out=siginvTFTsiginv(ispin), &
                  m_ST_out=ST(ispin))

               IF (perturbation_only) THEN
                  ! calculate objective function Tr(F_0 R)
                  IF (ispin .EQ. 1) energy_new = 0.0_dp
                  CALL dbcsr_dot(matrix_t_out(ispin), FTsiginv(ispin), energy_ispin)
                  energy_new = energy_new + energy_ispin*spin_factor
               ENDIF

               IF (penalty_occ_vol) THEN
!CALL dbcsr_print(matrix_t_out(ispin))
!CALL dbcsr_print(almo_scf_env%matrix_sigma(ispin))
                  CALL determinant(almo_scf_env%matrix_sigma(ispin), det1, &
                                   almo_scf_env%eps_filter)
                  penalty_func_new = &
                     -penalty_amplitude*spin_factor*nocc(ispin)* &
                     LOG(det1)
                  penalty_occ_vol_g_prefactor(ispin) = &
                     -2.0_dp*penalty_amplitude*spin_factor*nocc(ispin)
                  penalty_occ_vol_h_prefactor(ispin) = 0.0_dp
                  !!!! different penalty functions (lead to non-quadratic searches)
                  !   penalty_func_new = &
                  !           penalty_amplitude * spin_factor * nocc * &
                  !           (LOG(det1))**penalty_order
                  !   penalty_occ_vol_g_prefactor(ispin) = &
                  !           2.0_dp * &
                  !           penalty_amplitude * spin_factor * nocc * &
                  !           penalty_order * &
                  !           (LOG(det1))**(penalty_order-1)
                  !   penalty_occ_vol_h_prefactor(ispin) = &
                  !           2.0_dp * (penalty_order-1) / LOG(det1)
                  !
                  !   penalty_func_new = &
                  !           penalty_amplitude * spin_factor * nocc * &
                  !           (1.0_dp-1.0_dp/det1)**penalty_order
                  !   penalty_occ_vol_g_prefactor(ispin) = &
                  !           2.0_dp * &
                  !           penalty_amplitude * spin_factor * nocc * &
                  !           penalty_order * &
                  !           (1.0_dp/det1) * &
                  !           (1.0_dp-1.0_dp/det1)**(penalty_order-1)
                  !   penalty_occ_vol_h_prefactor(ispin) = &
                  !           2.0_dp * (penalty_order-det1) / (det1 - 1.0_dp)
                  !
                  !   penalty_func_new = &
                  !           penalty_amplitude * spin_factor * nocc * &
                  !           (1.0_dp-det1)**penalty_order
                  !   penalty_occ_vol_g_prefactor(ispin) = &
                  !           - 2.0_dp * &
                  !           penalty_amplitude * spin_factor * nocc * &
                  !           penalty_order * &
                  !           det1 * &
                  !           (1.0_dp-det1)**(penalty_order-1)
                  !   penalty_occ_vol_h_prefactor = &
                  !           2.0_dp * ( 1.0_dp - (penalty_order-1) * det1 / (1.0_dp-det1) )
                  IF (unit_nr > 0) THEN
                     WRITE (unit_nr, *) "penalty c0: ", penalty_occ_vol_g_prefactor(ispin)
                     WRITE (unit_nr, *) "penalty c1: ", penalty_occ_vol_h_prefactor(ispin)
                     WRITE (unit_nr, *) "penalty c0*c1: ", penalty_occ_vol_g_prefactor(ispin)*penalty_occ_vol_h_prefactor(ispin)
                     WRITE (unit_nr, *) "energy, penalty: ", energy_new, penalty_func_new
                  ENDIF
                  ! this is not pure energy anymore
                  energy_new = energy_new + penalty_func_new
               ENDIF

            ENDDO ! ispin
            !!!! -- end objective function

            DO ispin = 1, nspins

               IF (just_started .AND. almo_mathematica) THEN
                  IF (ispin .GT. 1) CPWARN("Mathematica files will be overwritten")
                  CALL print_mathematica_matrix(almo_scf_env%matrix_s(1), "matrixS.dat")
                  CALL print_mathematica_matrix(almo_scf_env%matrix_ks(ispin), "matrixF.dat")
                  CALL print_mathematica_matrix(matrix_t_out(ispin), "matrixT.dat")
                  CALL print_mathematica_matrix(quench_t(ispin), "matrixQ.dat")
               ENDIF

               ! save the previous gradient to compute beta
               ! do it only if the previous grad was computed
               ! for .NOT.line_search
               IF (line_search_iteration .EQ. 0 .AND. iteration .NE. 0) &
                  CALL dbcsr_copy(prev_grad(ispin), grad(ispin))

            ENDDO ! ispin

            ! compute the energy gradient if necessary
            skip_grad = (iteration .GT. 0 .AND. &
                         fixed_line_search_niter .NE. 0 .AND. &
                         line_search_iteration .NE. fixed_line_search_niter)

            IF (.NOT. skip_grad) THEN

               DO ispin = 1, nspins

                  CALL compute_gradient( &
                     m_grad_out=grad(ispin), &
                     m_ks=almo_scf_env%matrix_ks(ispin), &
                     m_s=almo_scf_env%matrix_s(1), &
                     m_t=matrix_t_out(ispin), &
                     m_t0=almo_scf_env%matrix_t_blk(ispin), &
                     m_siginv=almo_scf_env%matrix_sigma_inv(ispin), &
                     m_quench_t=quench_t(ispin), &
                     m_FTsiginv=FTsiginv(ispin), &
                     m_siginvTFTsiginv=siginvTFTsiginv(ispin), &
                     m_ST=ST(ispin), &
                     m_STsiginv0=STsiginv_0(ispin), &
                     m_theta=m_theta(ispin), &
                     m_sig_sqrti_ii=m_sig_sqrti_ii(ispin), &
                     domain_s_inv=almo_scf_env%domain_s_inv(:, ispin), &
                     domain_r_down=domain_r_down(:, ispin), &
                     cpu_of_domain=almo_scf_env%cpu_of_domain, &
                     domain_map=almo_scf_env%domain_map(ispin), &
                     assume_t0_q0x=assume_t0_q0x, &
                     optimize_theta=optimize_theta, &
                     normalize_orbitals=normalize_orbitals, &
                     penalty_occ_vol=penalty_occ_vol, &
                     penalty_occ_vol_prefactor=penalty_occ_vol_g_prefactor(ispin), &
                     envelope_amplitude=almo_scf_env%envelope_amplitude, &
                     eps_filter=almo_scf_env%eps_filter, &
                     spin_factor=spin_factor, &
                     special_case=my_special_case)

               ENDDO ! ispin

            ENDIF ! skip_grad

            ! if unprojected XALMOs are optimized then compute both
            ! HessianInv/preconditioner and the "bad-mode" projector

            IF (blissful_neglect) THEN
               DO ispin = 1, nspins
                  !compute the prec only for the first step,
                  !but project the gradient every step
                  IF (iteration .EQ. 0) THEN
                     CALL compute_preconditioner( &
                        domain_prec_out=almo_scf_env%domain_preconditioner(:, ispin), &
                        bad_modes_projector_down_out=bad_modes_projector_down(:, ispin), &
                        m_prec_out=prec_vv(ispin), &
                        m_ks=almo_scf_env%matrix_ks(ispin), &
                        m_s=almo_scf_env%matrix_s(1), &
                        m_siginv=almo_scf_env%matrix_sigma_inv(ispin), &
                        m_quench_t=quench_t(ispin), &
                        m_FTsiginv=FTsiginv(ispin), &
                        m_siginvTFTsiginv=siginvTFTsiginv(ispin), &
                        m_ST=ST(ispin), &
                        para_env=almo_scf_env%para_env, &
                        blacs_env=almo_scf_env%blacs_env, &
                        nocc_of_domain=almo_scf_env%nocc_of_domain(:, ispin), &
                        domain_s_inv=almo_scf_env%domain_s_inv(:, ispin), &
                        domain_s_inv_half=almo_scf_env%domain_s_sqrt_inv(:, ispin), &
                        domain_s_half=almo_scf_env%domain_s_sqrt(:, ispin), &
                        domain_r_down=domain_r_down(:, ispin), &
                        cpu_of_domain=almo_scf_env%cpu_of_domain, &
                        domain_map=almo_scf_env%domain_map(ispin), &
                        assume_t0_q0x=assume_t0_q0x, &
                        penalty_occ_vol=penalty_occ_vol, &
                        penalty_occ_vol_prefactor=penalty_occ_vol_g_prefactor(ispin), &
                        eps_filter=almo_scf_env%eps_filter, &
                        neg_thr=optimizer%neglect_threshold, &
                        spin_factor=spin_factor, &
                        special_case=my_special_case)
                  ENDIF
                  ! remove bad modes from the gradient
                  CALL apply_domain_operators( &
                     matrix_in=grad(ispin), &
                     matrix_out=grad(ispin), &
                     operator1=almo_scf_env%domain_s_inv(:, ispin), &
                     operator2=bad_modes_projector_down(:, ispin), &
                     dpattern=quench_t(ispin), &
                     map=almo_scf_env%domain_map(ispin), &
                     node_of_domain=almo_scf_env%cpu_of_domain, &
                     my_action=1, &
                     filter_eps=almo_scf_env%eps_filter)

               ENDDO ! ispin

            ENDIF ! blissful neglect

            ! check convergence and other exit criteria
            DO ispin = 1, nspins
               CALL dbcsr_norm(grad(ispin), dbcsr_norm_maxabsnorm, &
                               norm_scalar=grad_norm_spin(ispin))
               !grad_norm_frob = dbcsr_frobenius_norm(grad(ispin)) / &
               !                 dbcsr_frobenius_norm(quench_t(ispin))
               !IF (unit_nr > 0 ) WRITE(*,*) "Gradient RMS norm: ", grad_norm_frob
            ENDDO ! ispin
            grad_norm = MAXVAL(grad_norm_spin)

            converged = (grad_norm .LE. optimizer%eps_error)
            IF (converged .OR. (iteration .GE. max_iter)) THEN
               prepare_to_exit = .TRUE.
            ENDIF
            ! if early stopping is on do at least one iteration
            IF (optimizer%early_stopping_on .AND. just_started) &
               prepare_to_exit = .FALSE.

            IF (grad_norm .LT. almo_scf_env%eps_prev_guess) &
               use_guess = .TRUE.

            ! it is not time to exit just yet
            IF (.NOT. prepare_to_exit) THEN

               ! check the gradient along the step direction
               ! and decide whether to switch to the line-search mode
               ! do not do this in the first iteration
               IF (iteration .NE. 0) THEN

                  IF (fixed_line_search_niter .EQ. 0) THEN

                     ! enforce at least one line search
                     ! without even checking the error
                     IF (.NOT. line_search) THEN

                        line_search = .TRUE.
                        line_search_iteration = line_search_iteration + 1

                     ELSE

                        ! check the line-search error and decide whether to
                        ! change the direction
                        line_search_error = 0.0_dp
                        denom = 0.0_dp
                        denom2 = 0.0_dp

                        DO ispin = 1, nspins

                           CALL dbcsr_dot(grad(ispin), step(ispin), tempreal)
                           line_search_error = line_search_error + tempreal
                           CALL dbcsr_dot(grad(ispin), grad(ispin), tempreal)
                           denom = denom + tempreal
                           CALL dbcsr_dot(step(ispin), step(ispin), tempreal)
                           denom2 = denom2 + tempreal

                        ENDDO ! ispin

                        ! cosine of the angle between the step and grad
                        ! (must be close to zero at convergence)
                        line_search_error = line_search_error/SQRT(denom)/SQRT(denom2)

                        IF (ABS(line_search_error) .GT. optimizer%lin_search_eps_error) THEN
                           line_search = .TRUE.
                           line_search_iteration = line_search_iteration + 1
                        ELSE
                           line_search = .FALSE.
                           line_search_iteration = 0
                           IF (grad_norm .LT. eps_skip_gradients) THEN
                              fixed_line_search_niter = ABS(almo_scf_env%integer04)
                           ENDIF
                        ENDIF

                     ENDIF

                  ELSE ! decision for fixed_line_search_niter

                     IF (.NOT. line_search) THEN
                        line_search = .TRUE.
                        line_search_iteration = line_search_iteration + 1
                     ELSE
                        IF (line_search_iteration .EQ. fixed_line_search_niter) THEN
                           line_search = .FALSE.
                           line_search_iteration = 0
                           line_search_iteration = line_search_iteration + 1
                        ENDIF
                     ENDIF

                  ENDIF ! fixed_line_search_niter fork

               ENDIF ! iteration.ne.0

               IF (line_search) THEN
                  energy_diff = 0.0_dp
               ELSE
                  energy_diff = energy_new - energy_old
                  energy_old = energy_new
               ENDIF

               ! update the step direction
               IF (.NOT. line_search) THEN

                  !IF (unit_nr>0) THEN
                  !   WRITE(unit_nr,*) "....updating step direction...."
                  !ENDIF

                  cg_iteration = cg_iteration + 1

                  ! save the previous step
                  DO ispin = 1, nspins
                     CALL dbcsr_copy(prev_step(ispin), step(ispin))
                  ENDDO ! ispin

                  ! compute the new step (apply preconditioner if available)
                  SELECT CASE (prec_type)
                  CASE (xalmo_prec_full)

                     ! solving approximate Newton eq in the full (linearized) space
                     CALL newton_grad_to_step( &
                        optimizer=almo_scf_env%opt_xalmo_newton_pcg_solver, &
                        m_grad=grad(:), &
                        m_delta=step(:), &
                        m_s=almo_scf_env%matrix_s(:), &
                        m_ks=almo_scf_env%matrix_ks(:), &
                        m_siginv=almo_scf_env%matrix_sigma_inv(:), &
                        m_quench_t=quench_t(:), &
                        m_FTsiginv=FTsiginv(:), &
                        m_siginvTFTsiginv=siginvTFTsiginv(:), &
                        m_ST=ST(:), &
                        m_t=matrix_t_out(:), &
                        m_sig_sqrti_ii=m_sig_sqrti_ii(:), &
                        domain_s_inv=almo_scf_env%domain_s_inv(:, :), &
                        domain_r_down=domain_r_down(:, :), &
                        domain_map=almo_scf_env%domain_map(:), &
                        cpu_of_domain=almo_scf_env%cpu_of_domain, &
                        nocc_of_domain=almo_scf_env%nocc_of_domain(:, :), &
                        para_env=almo_scf_env%para_env, &
                        blacs_env=almo_scf_env%blacs_env, &
                        eps_filter=almo_scf_env%eps_filter, &
                        optimize_theta=optimize_theta, &
                        penalty_occ_vol=penalty_occ_vol, &
                        normalize_orbitals=normalize_orbitals, &
                        penalty_occ_vol_prefactor=penalty_occ_vol_g_prefactor(:), &
                        penalty_occ_vol_pf2=penalty_occ_vol_h_prefactor(:), &
                        special_case=my_special_case &
                        )

                  CASE (xalmo_prec_domain)

                     ! compute and invert preconditioner?
                     IF (.NOT. blissful_neglect .AND. &
                         ((just_started .AND. perturbation_only) .OR. &
                          (iteration .EQ. 0 .AND. (.NOT. perturbation_only))) &
                         ) THEN

                        ! computing preconditioner
                        DO ispin = 1, nspins
                           CALL compute_preconditioner( &
                              domain_prec_out=almo_scf_env%domain_preconditioner(:, ispin), &
                              m_prec_out=prec_vv(ispin), &
                              m_ks=almo_scf_env%matrix_ks(ispin), &
                              m_s=almo_scf_env%matrix_s(1), &
                              m_siginv=almo_scf_env%matrix_sigma_inv(ispin), &
                              m_quench_t=quench_t(ispin), &
                              m_FTsiginv=FTsiginv(ispin), &
                              m_siginvTFTsiginv=siginvTFTsiginv(ispin), &
                              m_ST=ST(ispin), &
                              para_env=almo_scf_env%para_env, &
                              blacs_env=almo_scf_env%blacs_env, &
                              nocc_of_domain=almo_scf_env%nocc_of_domain(:, ispin), &
                              domain_s_inv=almo_scf_env%domain_s_inv(:, ispin), &
                              domain_r_down=domain_r_down(:, ispin), &
                              cpu_of_domain=almo_scf_env%cpu_of_domain, &
                              domain_map=almo_scf_env%domain_map(ispin), &
                              assume_t0_q0x=assume_t0_q0x, &
                              penalty_occ_vol=penalty_occ_vol, &
                              penalty_occ_vol_prefactor=penalty_occ_vol_g_prefactor(ispin), &
                              eps_filter=almo_scf_env%eps_filter, &
                              neg_thr=0.5_dp, &
                              spin_factor=spin_factor, &
                              special_case=my_special_case)
                        ENDDO ! ispin
                     ENDIF ! compute_prec

                     !IF (unit_nr>0) THEN
                     !   WRITE(unit_nr,*) "....applying precomputed preconditioner...."
                     !ENDIF

                     IF (my_special_case .EQ. xalmo_case_block_diag .OR. &
                         my_special_case .EQ. xalmo_case_fully_deloc) THEN

                        DO ispin = 1, nspins

                           CALL dbcsr_multiply("N", "N", -1.0_dp, &
                                               prec_vv(ispin), &
                                               grad(ispin), &
                                               0.0_dp, step(ispin), &
                                               filter_eps=almo_scf_env%eps_filter)

                        ENDDO ! ispin

                     ELSE

                        !!! RZK-warning Currently for non-theta only
                        IF (optimize_theta) THEN
                           CPABORT("theta is NYI")
                        ENDIF

                        DO ispin = 1, nspins

                           CALL apply_domain_operators( &
                              matrix_in=grad(ispin), &
                              matrix_out=step(ispin), &
                              operator1=almo_scf_env%domain_preconditioner(:, ispin), &
                              dpattern=quench_t(ispin), &
                              map=almo_scf_env%domain_map(ispin), &
                              node_of_domain=almo_scf_env%cpu_of_domain, &
                              my_action=0, &
                              filter_eps=almo_scf_env%eps_filter)
                           CALL dbcsr_scale(step(ispin), -1.0_dp)

                           !CALL dbcsr_copy(m_tmp_no_3,&
                           !        quench_t(ispin))
                           !CALL dbcsr_function_of_elements(m_tmp_no_3,&
                           !        func=dbcsr_func_inverse,&
                           !        a0=0.0_dp,&
                           !        a1=1.0_dp)
                           !CALL dbcsr_copy(m_tmp_no_2,step)
                           !CALL dbcsr_hadamard_product(&
                           !        m_tmp_no_2,&
                           !        m_tmp_no_3,&
                           !        step)
                           !CALL dbcsr_copy(m_tmp_no_3,quench_t(ispin))

                        ENDDO ! ispin

                     ENDIF ! special case

                  CASE (xalmo_prec_zero)

                     ! no preconditioner
                     DO ispin = 1, nspins

                        CALL dbcsr_copy(step(ispin), grad(ispin))
                        CALL dbcsr_scale(step(ispin), -1.0_dp)

                     ENDDO ! ispin

                  END SELECT ! preconditioner type fork

                  ! check whether we need to reset conjugate directions
                  IF (iteration .EQ. 0) THEN
                     reset_conjugator = .TRUE.
                  ENDIF

                  ! compute the conjugation coefficient - beta
                  IF (.NOT. reset_conjugator) THEN

                     CALL compute_cg_beta( &
                        beta=beta, &
                        reset_conjugator=reset_conjugator, &
                        conjugator=optimizer%conjugator, &
                        grad=grad(:), &
                        prev_grad=prev_grad(:), &
                        step=step(:), &
                        prev_step=prev_step(:), &
                        prev_minus_prec_grad=prev_minus_prec_grad(:) &
                        )

                  ENDIF

                  IF (reset_conjugator) THEN

                     beta = 0.0_dp
                     IF (unit_nr > 0 .AND. (.NOT. just_started)) THEN
                        WRITE (unit_nr, '(T2,A35)') "Re-setting conjugator to zero"
                     ENDIF
                     reset_conjugator = .FALSE.

                  ENDIF

                  ! save the preconditioned gradient (useful for beta)
                  DO ispin = 1, nspins

                     CALL dbcsr_copy(prev_minus_prec_grad(ispin), step(ispin))

                     !IF (unit_nr>0) THEN
                     !   WRITE(unit_nr,*) "....final beta....", beta
                     !ENDIF

                     ! conjugate the step direction
                     CALL dbcsr_add(step(ispin), prev_step(ispin), 1.0_dp, beta)

                  ENDDO ! ispin

               ENDIF ! update the step direction

               ! estimate the step size
               IF (.NOT. line_search) THEN
                  ! we just changed the direction and
                  ! we have only E and grad from the current step
                  ! it is not enouhg to compute step_size - just guess it
                  e0 = energy_new
                  g0 = 0.0_dp
                  DO ispin = 1, nspins
                     CALL dbcsr_dot(grad(ispin), step(ispin), tempreal)
                     g0 = g0 + tempreal
                  ENDDO ! ispin
                  IF (iteration .EQ. 0) THEN
                     step_size = optimizer%lin_search_step_size_guess
                  ELSE
                     IF (next_step_size_guess .LE. 0.0_dp) THEN
                        step_size = optimizer%lin_search_step_size_guess
                     ELSE
                        ! take the last value
                        step_size = next_step_size_guess*1.05_dp
                     ENDIF
                  ENDIF
                  !IF (unit_nr > 0) THEN
                  !   WRITE (unit_nr, '(A2,3F12.5)') &
                  !      "EG", e0, g0, step_size
                  !ENDIF
                  next_step_size_guess = step_size
               ELSE
                  IF (fixed_line_search_niter .EQ. 0) THEN
                     e1 = energy_new
                     g1 = 0.0_dp
                     DO ispin = 1, nspins
                        CALL dbcsr_dot(grad(ispin), step(ispin), tempreal)
                        g1 = g1 + tempreal
                     ENDDO ! ispin
                     ! we have accumulated some points along this direction
                     ! use only the most recent g0 (quadratic approximation)
                     appr_sec_der = (g1 - g0)/step_size
                     !IF (unit_nr > 0) THEN
                     !   WRITE (unit_nr, '(A2,7F12.5)') &
                     !      "EG", e0, e1, g0, g1, appr_sec_der, step_size, -g1/appr_sec_der
                     !ENDIF
                     step_size = -g1/appr_sec_der
                     e0 = e1
                     g0 = g1
                  ELSE
                     ! use e0, g0 and e1 to compute g1 and make a step
                     ! if the next iteration is also line_search
                     ! use e1 and the calculated g1 as e0 and g0
                     e1 = energy_new
                     appr_sec_der = 2.0*((e1 - e0)/step_size - g0)/step_size
                     g1 = appr_sec_der*step_size + g0
                     !IF (unit_nr > 0) THEN
                     !   WRITE (unit_nr, '(A2,7F12.5)') &
                     !      "EG", e0, e1, g0, g1, appr_sec_der, step_size, -g1/appr_sec_der
                     !ENDIF
                     !appr_sec_der=(g1-g0)/step_size
                     step_size = -g1/appr_sec_der
                     e0 = e1
                     g0 = g1
                  ENDIF
                  next_step_size_guess = next_step_size_guess + step_size
               ENDIF

               ! update theta
               DO ispin = 1, nspins
                  CALL dbcsr_add(m_theta(ispin), step(ispin), 1.0_dp, step_size)
               ENDDO ! ispin

            ENDIF ! not.prepare_to_exit

            IF (line_search) THEN
               iter_type = "LS"
            ELSE
               iter_type = "CG"
            ENDIF

            t2 = m_walltime()
            IF (unit_nr > 0) THEN
               iter_type = TRIM("ALMO SCF "//iter_type)
               WRITE (unit_nr, '(T2,A13,I6,F23.10,E14.5,F14.9,F9.2)') &
                  iter_type, iteration, &
                  energy_new, energy_diff, grad_norm, &
                  t2 - t1
               IF (penalty_occ_vol) THEN
                  WRITE (unit_nr, '(T2,A19,F23.10)') &
                     "Energy component:", energy_new - penalty_func_new
                  WRITE (unit_nr, '(T2,A19,F23.10)') &
                     "Penalty component:", penalty_func_new
               ENDIF
            ENDIF

            IF (my_special_case .EQ. xalmo_case_block_diag) THEN
               IF (penalty_occ_vol) THEN
                  almo_scf_env%almo_scf_energy = energy_new - penalty_func_new
               ELSE
                  almo_scf_env%almo_scf_energy = energy_new
               ENDIF
            ENDIF

            t1 = m_walltime()

            iteration = iteration + 1
            IF (prepare_to_exit) EXIT

         ENDDO ! inner SCF loop

         IF (converged .OR. (outer_iteration .GE. outer_max_iter)) THEN
            outer_prepare_to_exit = .TRUE.
         ENDIF

         outer_iteration = outer_iteration + 1
         IF (outer_prepare_to_exit) EXIT

      ENDDO ! outer SCF loop

      DO ispin = 1, nspins
         IF (converged .AND. almo_mathematica) THEN
            IF (ispin .GT. 1) CPWARN("Mathematica files will be overwritten")
            CALL print_mathematica_matrix(matrix_t_out(ispin), "matrixTf.dat")
         ENDIF
      ENDDO ! ispin

      ! post SCF-loop calculations
      IF (converged) THEN

         ! RZK-warning: must obtain MO coefficients from final theta

         IF (perturbation_only) THEN

            ! return perturbed density to qs_env
            CALL almo_dm_to_qs_env(qs_env, almo_scf_env%matrix_p, &
                                   almo_scf_env%mat_distr_aos)

            ! compute energy correction and perform
            ! detailed decomposition analysis (if requested)
            ! reuse step and grad matrices to store decomposition results
            CALL xalmo_analysis( &
               detailed_analysis=almo_scf_env%almo_analysis%do_analysis, &
               eps_filter=almo_scf_env%eps_filter, &
               m_T_in=matrix_t_out(:), &
               m_T0_in=almo_scf_env%matrix_t_blk(:), &
               m_siginv_in=almo_scf_env%matrix_sigma_inv(:), &
               m_siginv0_in=almo_scf_env%matrix_sigma_inv_0deloc(:), &
               m_S_in=almo_scf_env%matrix_s(:), &
               m_KS0_in=almo_scf_env%matrix_ks_0deloc(:), &
               m_quench_t_in=quench_t(:), &
               energy_out=energy_new, &
               m_eda_out=step(:), &
               m_cta_out=grad(:) &
               )

            IF (almo_scf_env%almo_analysis%do_analysis) THEN

               DO ispin = 1, nspins

                  ! energy decomposition analysis (EDA)
                  IF (unit_nr > 0) THEN
                     WRITE (unit_nr, '(T2,A)') "DECOMPOSITION OF THE DELOCALIZATION ENERGY"
                  ENDIF

                  ! open the output file, print and close
                  CALL get_qs_env(qs_env, input=input)
                  almo_print_section => section_vals_get_subs_vals(input, "DFT%ALMO_SCF%ANALYSIS%PRINT")
                  eda_unit = cp_print_key_unit_nr(logger, almo_print_section, &
                                                  "ALMO_EDA_CT", extension=".dat", local=.TRUE.)
                  CALL dbcsr_print_block_sum(step(ispin), eda_unit)
                  CALL cp_print_key_finished_output(eda_unit, logger, almo_print_section, &
                                                    "ALMO_EDA_CT", local=.TRUE.)

                  ! charge transfer analysis (CTA)
                  IF (unit_nr > 0) THEN
                     WRITE (unit_nr, '(T2,A)') "DECOMPOSITION OF CHARGE TRANSFER TERMS"
                  ENDIF

                  eda_unit = cp_print_key_unit_nr(logger, almo_print_section, &
                                                  "ALMO_CTA", extension=".dat", local=.TRUE.)
                  CALL dbcsr_print_block_sum(grad(ispin), eda_unit)
                  CALL cp_print_key_finished_output(eda_unit, logger, almo_print_section, &
                                                    "ALMO_CTA", local=.TRUE.)

               ENDDO ! ispin

            ENDIF ! do ALMO EDA/CTA

            ! print out the energy lowering
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, *)
               WRITE (unit_nr, '(T2,A35,F25.10)') "ENERGY OF BLOCK-DIAGONAL ALMOs:", &
                  almo_scf_env%almo_scf_energy
               WRITE (unit_nr, '(T2,A35,F25.10)') "ENERGY LOWERING:", &
                  energy_new
               WRITE (unit_nr, '(T2,A35,F25.10)') "CORRECTED ENERGY:", &
                  almo_scf_env%almo_scf_energy + energy_new
               WRITE (unit_nr, *)
            ENDIF
            CALL almo_scf_update_ks_energy(qs_env, &
                                           energy=almo_scf_env%almo_scf_energy, &
                                           energy_singles_corr=energy_new)

         ELSE ! non-perturbative

            CALL almo_scf_update_ks_energy(qs_env, &
                                           energy=energy_new)

         ENDIF ! if perturbation only

      ENDIF ! if converged

      DO ispin = 1, nspins
         CALL dbcsr_release(prec_vv(ispin))
         CALL dbcsr_release(STsiginv_0(ispin))
         CALL dbcsr_release(ST(ispin))
         CALL dbcsr_release(FTsiginv(ispin))
         CALL dbcsr_release(siginvTFTsiginv(ispin))
         CALL dbcsr_release(prev_grad(ispin))
         CALL dbcsr_release(prev_step(ispin))
         CALL dbcsr_release(grad(ispin))
         CALL dbcsr_release(step(ispin))
         CALL dbcsr_release(prev_minus_prec_grad(ispin))
         CALL dbcsr_release(m_theta(ispin))
         CALL dbcsr_release(m_t_in_local(ispin))
         CALL dbcsr_release(m_sig_sqrti_ii(ispin))
         CALL release_submatrices(domain_r_down(:, ispin))
         CALL release_submatrices(bad_modes_projector_down(:, ispin))
      ENDDO ! ispin

      DEALLOCATE (prec_vv)
      DEALLOCATE (siginvTFTsiginv)
      DEALLOCATE (STsiginv_0)
      DEALLOCATE (FTsiginv)
      DEALLOCATE (ST)
      DEALLOCATE (prev_grad)
      DEALLOCATE (grad)
      DEALLOCATE (prev_step)
      DEALLOCATE (step)
      DEALLOCATE (prev_minus_prec_grad)
      DEALLOCATE (m_sig_sqrti_ii)

      DEALLOCATE (domain_r_down)
      DEALLOCATE (bad_modes_projector_down)

      DEALLOCATE (penalty_occ_vol_g_prefactor)
      DEALLOCATE (penalty_occ_vol_h_prefactor)
      DEALLOCATE (grad_norm_spin)
      DEALLOCATE (nocc)

      DEALLOCATE (m_theta, m_t_in_local)

      IF (.NOT. converged .AND. .NOT. optimizer%early_stopping_on) THEN
         CPABORT("Optimization not converged! ")
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE almo_scf_xalmo_pcg

! **************************************************************************************************
!> \brief Analysis of the orbitals
!> \param detailed_analysis ...
!> \param eps_filter ...
!> \param m_T_in ...
!> \param m_T0_in ...
!> \param m_siginv_in ...
!> \param m_siginv0_in ...
!> \param m_S_in ...
!> \param m_KS0_in ...
!> \param m_quench_t_in ...
!> \param energy_out ...
!> \param m_eda_out ...
!> \param m_cta_out ...
!> \par History
!>       2017.07 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE xalmo_analysis(detailed_analysis, eps_filter, m_T_in, m_T0_in, &
                             m_siginv_in, m_siginv0_in, m_S_in, m_KS0_in, m_quench_t_in, energy_out, &
                             m_eda_out, m_cta_out)

      LOGICAL, INTENT(IN)                                :: detailed_analysis
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter
      TYPE(dbcsr_type), DIMENSION(:), INTENT(IN)         :: m_T_in, m_T0_in, m_siginv_in, &
                                                            m_siginv0_in, m_S_in, m_KS0_in, &
                                                            m_quench_t_in
      REAL(KIND=dp), INTENT(INOUT)                       :: energy_out
      TYPE(dbcsr_type), DIMENSION(:), INTENT(INOUT)      :: m_eda_out, m_cta_out

      CHARACTER(len=*), PARAMETER :: routineN = 'xalmo_analysis', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nspins
      REAL(KIND=dp)                                      :: energy_ispin, spin_factor
      TYPE(dbcsr_type)                                   :: FTsiginv0, Fvo0, m_X, siginvTFTsiginv0, &
                                                            ST0

      CALL timeset(routineN, handle)

      nspins = SIZE(m_T_in)

      IF (nspins == 1) THEN
         spin_factor = 2.0_dp
      ELSE
         spin_factor = 1.0_dp
      ENDIF

      energy_out = 0.0_dp
      DO ispin = 1, nspins

         ! create temporary matrices
         CALL dbcsr_create(Fvo0, &
                           template=m_T_in(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(FTsiginv0, &
                           template=m_T_in(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(ST0, &
                           template=m_T_in(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(m_X, &
                           template=m_T_in(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(siginvTFTsiginv0, &
                           template=m_siginv0_in(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)

         ! compute F_{virt,occ} for the zero-delocalization state
         CALL compute_frequently_used_matrices( &
            filter_eps=eps_filter, &
            m_T_in=m_T0_in(ispin), &
            m_siginv_in=m_siginv0_in(ispin), &
            m_S_in=m_S_in(1), &
            m_F_in=m_KS0_in(ispin), &
            m_FTsiginv_out=FTsiginv0, &
            m_siginvTFTsiginv_out=siginvTFTsiginv0, &
            m_ST_out=ST0)
         CALL dbcsr_copy(Fvo0, m_quench_t_in(ispin))
         CALL dbcsr_copy(Fvo0, FTsiginv0, keep_sparsity=.TRUE.)
         CALL dbcsr_multiply("N", "N", -1.0_dp, &
                             ST0, &
                             siginvTFTsiginv0, &
                             1.0_dp, Fvo0, &
                             retain_sparsity=.TRUE.)

         ! get single excitation amplitudes
         CALL dbcsr_copy(m_X, m_T0_in(ispin))
         CALL dbcsr_add(m_X, m_T_in(ispin), -1.0_dp, 1.0_dp)

         CALL dbcsr_dot(m_X, Fvo0, energy_ispin)
         energy_out = energy_out + energy_ispin*spin_factor

         IF (detailed_analysis) THEN

            CALL dbcsr_hadamard_product(m_X, Fvo0, m_eda_out(ispin))
            CALL dbcsr_scale(m_eda_out(ispin), spin_factor)
            CALL dbcsr_filter(m_eda_out(ispin), eps_filter)

            ! first, compute [QR'R]_mu^i = [(S-SRS).X.siginv']_mu^i
            ! a. FTsiginv0 = S.T0*siginv0
            CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                ST0, &
                                m_siginv0_in(ispin), &
                                0.0_dp, FTsiginv0, &
                                filter_eps=eps_filter)
            ! c. tmp1(use ST0) = S.X
            CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                m_S_in(1), &
                                m_X, &
                                0.0_dp, ST0, &
                                filter_eps=eps_filter)
            ! d. tmp2 = tr(T0).tmp1 = tr(T0).S.X
            CALL dbcsr_multiply("T", "N", 1.0_dp, &
                                m_T0_in(ispin), &
                                ST0, &
                                0.0_dp, siginvTFTsiginv0, &
                                filter_eps=eps_filter)
            ! e. tmp1 = tmp1 - tmp3.tmp2 = S.X - S.T0.siginv0*tr(T0).S.X
            !         = (1-S.R0).S.X
            CALL dbcsr_multiply("N", "N", -1.0_dp, &
                                FTsiginv0, &
                                siginvTFTsiginv0, &
                                1.0_dp, ST0, &
                                filter_eps=eps_filter)
            ! f. tmp2(use FTsiginv0) = tmp1*siginv
            CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                ST0, &
                                m_siginv_in(ispin), &
                                0.0_dp, FTsiginv0, &
                                filter_eps=eps_filter)
            ! second, compute traces of blocks [RR'Q]^x_y * [X]^y_x
            CALL dbcsr_hadamard_product(m_X, &
                                        FTsiginv0, m_cta_out(ispin))
            CALL dbcsr_scale(m_cta_out(ispin), spin_factor)
            CALL dbcsr_filter(m_cta_out(ispin), eps_filter)

         ENDIF ! do ALMO EDA/CTA

         CALL dbcsr_release(Fvo0)
         CALL dbcsr_release(FTsiginv0)
         CALL dbcsr_release(ST0)
         CALL dbcsr_release(m_X)
         CALL dbcsr_release(siginvTFTsiginv0)

      ENDDO ! ispin

      CALL timestop(handle)

   END SUBROUTINE xalmo_analysis

! **************************************************************************************************
!> \brief Compute matrices that are used often in various parts of the
!>        optimization procedure
!> \param filter_eps ...
!> \param m_T_in ...
!> \param m_siginv_in ...
!> \param m_S_in ...
!> \param m_F_in ...
!> \param m_FTsiginv_out ...
!> \param m_siginvTFTsiginv_out ...
!> \param m_ST_out ...
!> \par History
!>       2016.12 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE compute_frequently_used_matrices(filter_eps, &
                                               m_T_in, m_siginv_in, m_S_in, m_F_in, m_FTsiginv_out, &
                                               m_siginvTFTsiginv_out, m_ST_out)

      REAL(KIND=dp), INTENT(IN)                          :: filter_eps
      TYPE(dbcsr_type), INTENT(IN)                       :: m_T_in, m_siginv_in, m_S_in, m_F_in
      TYPE(dbcsr_type), INTENT(INOUT)                    :: m_FTsiginv_out, m_siginvTFTsiginv_out, &
                                                            m_ST_out

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_frequently_used_matrices', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      TYPE(dbcsr_type)                                   :: m_tmp_no_1, m_tmp_oo_1

      CALL timeset(routineN, handle)

      CALL dbcsr_create(m_tmp_no_1, &
                        template=m_T_in, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(m_tmp_oo_1, &
                        template=m_siginv_in, &
                        matrix_type=dbcsr_type_no_symmetry)

      CALL dbcsr_multiply("N", "N", 1.0_dp, &
                          m_F_in, &
                          m_T_in, &
                          0.0_dp, m_tmp_no_1, &
                          filter_eps=filter_eps)

      CALL dbcsr_multiply("N", "N", 1.0_dp, &
                          m_tmp_no_1, &
                          m_siginv_in, &
                          0.0_dp, m_FTsiginv_out, &
                          filter_eps=filter_eps)

      CALL dbcsr_multiply("T", "N", 1.0_dp, &
                          m_T_in, &
                          m_FTsiginv_out, &
                          0.0_dp, m_tmp_oo_1, &
                          filter_eps=filter_eps)

      CALL dbcsr_multiply("N", "N", 1.0_dp, &
                          m_siginv_in, &
                          m_tmp_oo_1, &
                          0.0_dp, m_siginvTFTsiginv_out, &
                          filter_eps=filter_eps)

      CALL dbcsr_multiply("N", "N", 1.0_dp, &
                          m_S_in, &
                          m_T_in, &
                          0.0_dp, m_ST_out, &
                          filter_eps=filter_eps)

      CALL dbcsr_release(m_tmp_no_1)
      CALL dbcsr_release(m_tmp_oo_1)

      CALL timestop(handle)

   END SUBROUTINE compute_frequently_used_matrices

! **************************************************************************************************
!> \brief Split the matrix of virtual orbitals into two:
!>        retained orbs and discarded
!> \param almo_scf_env ...
!> \par History
!>       2011.09 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE split_v_blk(almo_scf_env)

      TYPE(almo_scf_env_type)                            :: almo_scf_env

      CHARACTER(len=*), PARAMETER :: routineN = 'split_v_blk', routineP = moduleN//':'//routineN

      INTEGER                                            :: discarded_v, handle, iblock_col, &
                                                            iblock_col_size, iblock_row, &
                                                            iblock_row_size, ispin, retained_v
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: data_p, p_new_block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      DO ispin = 1, almo_scf_env%nspins

         CALL dbcsr_work_create(almo_scf_env%matrix_v_blk(ispin), &
                                work_mutable=.TRUE.)
         CALL dbcsr_work_create(almo_scf_env%matrix_v_disc_blk(ispin), &
                                work_mutable=.TRUE.)

         CALL dbcsr_iterator_start(iter, almo_scf_env%matrix_v_full_blk(ispin))

         DO WHILE (dbcsr_iterator_blocks_left(iter))

            CALL dbcsr_iterator_next_block(iter, iblock_row, iblock_col, data_p, &
                                           row_size=iblock_row_size, col_size=iblock_col_size)

            IF (iblock_row .NE. iblock_col) THEN
               CPABORT("off-diagonal block found")
            ENDIF

            retained_v = almo_scf_env%nvirt_of_domain(iblock_col, ispin)
            discarded_v = almo_scf_env%nvirt_disc_of_domain(iblock_col, ispin)
            CPASSERT(retained_v .GT. 0)
            CPASSERT(discarded_v .GT. 0)

            NULLIFY (p_new_block)
            CALL dbcsr_reserve_block2d(almo_scf_env%matrix_v_disc_blk(ispin), &
                                       iblock_row, iblock_col, p_new_block)
            CPASSERT(ASSOCIATED(p_new_block))
            CPASSERT(retained_v + discarded_v .EQ. iblock_col_size)
            p_new_block(:, :) = data_p(:, (retained_v + 1):iblock_col_size)

            NULLIFY (p_new_block)
            CALL dbcsr_reserve_block2d(almo_scf_env%matrix_v_blk(ispin), &
                                       iblock_row, iblock_col, p_new_block)
            CPASSERT(ASSOCIATED(p_new_block))
            p_new_block(:, :) = data_p(:, 1:retained_v)

         ENDDO ! iterator
         CALL dbcsr_iterator_stop(iter)

         CALL dbcsr_finalize(almo_scf_env%matrix_v_blk(ispin))
         CALL dbcsr_finalize(almo_scf_env%matrix_v_disc_blk(ispin))

      ENDDO ! ispin

      CALL timestop(handle)

   END SUBROUTINE split_v_blk

! **************************************************************************************************
!> \brief various methods for calculating the Harris-Foulkes correction
!> \param almo_scf_env ...
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE harris_foulkes_correction(almo_scf_env)

      TYPE(almo_scf_env_type)                            :: almo_scf_env

      CHARACTER(len=*), PARAMETER :: routineN = 'harris_foulkes_correction', &
         routineP = moduleN//':'//routineN
      INTEGER, PARAMETER                                 :: cayley_transform = 1, dm_ls_step = 2

      INTEGER :: algorithm_id, handle, handle1, handle2, handle3, handle4, handle5, handle6, &
         handle7, handle8, ispin, iteration, n, nmins, nspin, opt_k_max_iter, &
         outer_opt_k_iteration, outer_opt_k_max_iter, unit_nr
      INTEGER, DIMENSION(1)                              :: fake, nelectron_spin_real
      LOGICAL :: converged, line_search, md_in_k_space, outer_opt_k_prepare_to_exit, &
         prepare_to_exit, reset_conjugator, reset_step_size, use_cubic_approximation, &
         use_quadratic_approximation
      REAL(KIND=dp) :: aa, bb, beta, conjugacy_error, conjugacy_error_threshold, &
         delta_obj_function, denom, energy_correction_final, frob_matrix, frob_matrix_base, fun0, &
         fun1, gfun0, gfun1, grad_norm, grad_norm_frob, kappa, kin_energy, line_search_error, &
         line_search_error_threshold, num_threshold, numer, obj_function, quadratic_approx_error, &
         quadratic_approx_error_threshold, safety_multiplier, spin_factor, step_size, &
         step_size_quadratic_approx, step_size_quadratic_approx2, t1, t1a, t1cholesky, t2, t2a, &
         t2cholesky, tau, time_step, x_opt_eps_adaptive, x_opt_eps_adaptive_factor
      REAL(KIND=dp), DIMENSION(1)                        :: local_mu
      REAL(KIND=dp), DIMENSION(2)                        :: energy_correction
      REAL(KIND=dp), DIMENSION(3)                        :: minima
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(ct_step_env_type)                             :: ct_step_env
      TYPE(dbcsr_type) :: grad, k_vd_index_down, k_vr_index_down, matrix_k_central, matrix_tmp1, &
         matrix_tmp2, prec, prev_grad, prev_minus_prec_grad, prev_step, sigma_oo_curr, &
         sigma_oo_curr_inv, sigma_vv_sqrt, sigma_vv_sqrt_guess, sigma_vv_sqrt_inv, &
         sigma_vv_sqrt_inv_guess, step, t_curr, tmp1_n_vr, tmp2_n_o, tmp3_vd_vr, tmp4_o_vr, &
         tmp_k_blk, vd_fixed, vd_index_sqrt, vd_index_sqrt_inv, velocity, vr_fixed, vr_index_sqrt, &
         vr_index_sqrt_inv
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:)        :: matrix_p_almo_scf_converged

!TYPE(dbcsr_iterator_type)                  :: iter
!TYPE(dbcsr_type) :: tmp11,tmp22,tmp33
!REAL(kind=dp)                            :: k_var1, k_var2
!TYPE(dbcsr_type)                      :: fake_step
!
!TYPE(dbcsr_type)                      :: sigma_dr, sigma_dr2, sigma_rr, sigma_rr2
!TYPE(dbcsr_type)                      :: fake_a,fake_b,fake_k0
!INTEGER                                  :: retained_v,discarded_v,i_row,j_col
!TYPE(dbcsr_type) :: matrix_rst0, matrix_rst1, matrix_rst2, ss_vv
!REAL(KIND=dp)       :: filter_memorize, init_filter, occ_vv
!INTEGER             :: ppp

      CALL timeset(routineN, handle)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      nspin = almo_scf_env%nspins
      energy_correction_final = 0.0_dp
      IF (nspin .EQ. 1) THEN
         spin_factor = 2.0_dp
      ELSE
         spin_factor = 1.0_dp
      ENDIF

      IF (almo_scf_env%deloc_use_occ_orbs) THEN
         algorithm_id = cayley_transform
      ELSE
         algorithm_id = dm_ls_step
      ENDIF

      t1 = m_walltime()

      SELECT CASE (algorithm_id)
      CASE (cayley_transform)

         ! rescale density matrix by spin factor
         ! so the orbitals and density are consistent with each other
         IF (almo_scf_env%nspins == 1) THEN
            CALL dbcsr_scale(almo_scf_env%matrix_p(1), 1.0_dp/spin_factor)
         ENDIF

         ! transform matrix_t not matrix_t_blk (we might need ALMOs later)
         DO ispin = 1, nspin

            CALL dbcsr_copy(almo_scf_env%matrix_t(ispin), &
                            almo_scf_env%matrix_t_blk(ispin))

            ! obtain orthogonalization matrices for ALMOs
            ! RZK-warning - remove this sqrt(sigma) and inv(sqrt(sigma))
            ! ideally ALMO scf should use sigma and sigma_inv in
            ! the tensor_up_down representation

            IF (unit_nr > 0) THEN
               WRITE (unit_nr, *) "sqrt and inv(sqrt) of MO overlap matrix"
            ENDIF
            CALL dbcsr_create(almo_scf_env%matrix_sigma_sqrt(ispin), &
                              template=almo_scf_env%matrix_sigma(ispin), &
                              matrix_type=dbcsr_type_no_symmetry)
            CALL dbcsr_create(almo_scf_env%matrix_sigma_sqrt_inv(ispin), &
                              template=almo_scf_env%matrix_sigma(ispin), &
                              matrix_type=dbcsr_type_no_symmetry)

            CALL matrix_sqrt_Newton_Schulz(almo_scf_env%matrix_sigma_sqrt(ispin), &
                                           almo_scf_env%matrix_sigma_sqrt_inv(ispin), &
                                           almo_scf_env%matrix_sigma(ispin), &
                                           threshold=almo_scf_env%eps_filter, &
                                           order=almo_scf_env%order_lanczos, &
                                           eps_lanczos=almo_scf_env%eps_lanczos, &
                                           max_iter_lanczos=almo_scf_env%max_iter_lanczos)

            IF (safe_mode) THEN
               CALL dbcsr_create(matrix_tmp1, template=almo_scf_env%matrix_sigma(ispin), &
                                 matrix_type=dbcsr_type_no_symmetry)
               CALL dbcsr_create(matrix_tmp2, template=almo_scf_env%matrix_sigma(ispin), &
                                 matrix_type=dbcsr_type_no_symmetry)

               CALL dbcsr_multiply("N", "N", 1.0_dp, almo_scf_env%matrix_sigma_sqrt_inv(ispin), &
                                   almo_scf_env%matrix_sigma(ispin), &
                                   0.0_dp, matrix_tmp1, filter_eps=almo_scf_env%eps_filter)
               CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, &
                                   almo_scf_env%matrix_sigma_sqrt_inv(ispin), &
                                   0.0_dp, matrix_tmp2, filter_eps=almo_scf_env%eps_filter)

               frob_matrix_base = dbcsr_frobenius_norm(matrix_tmp2)
               CALL dbcsr_add_on_diag(matrix_tmp2, -1.0_dp)
               frob_matrix = dbcsr_frobenius_norm(matrix_tmp2)
               IF (unit_nr > 0) THEN
                  WRITE (unit_nr, *) "Error for (inv(sqrt(SIG))*SIG*inv(sqrt(SIG))-I)", frob_matrix/frob_matrix_base
               ENDIF

               CALL dbcsr_release(matrix_tmp1)
               CALL dbcsr_release(matrix_tmp2)
            ENDIF
         ENDDO

         IF (almo_scf_env%almo_update_algorithm .EQ. almo_scf_diag) THEN

            DO ispin = 1, nspin

               t1a = m_walltime()

               line_search_error_threshold = almo_scf_env%real01
               conjugacy_error_threshold = almo_scf_env%real02
               quadratic_approx_error_threshold = almo_scf_env%real03
               x_opt_eps_adaptive_factor = almo_scf_env%real04

               !! the outer loop for k optimization
               outer_opt_k_max_iter = almo_scf_env%opt_k_outer_max_iter
               outer_opt_k_prepare_to_exit = .FALSE.
               outer_opt_k_iteration = 0
               grad_norm = 0.0_dp
               grad_norm_frob = 0.0_dp
               CALL dbcsr_set(almo_scf_env%matrix_x(ispin), 0.0_dp)
               IF (almo_scf_env%deloc_truncate_virt .EQ. virt_full) outer_opt_k_max_iter = 0

               DO

                  ! obtain proper retained virtuals (1-R)|ALMO_vr>
                  CALL apply_projector(psi_in=almo_scf_env%matrix_v_blk(ispin), &
                                       psi_out=almo_scf_env%matrix_v(ispin), &
                                       psi_projector=almo_scf_env%matrix_t_blk(ispin), &
                                       metric=almo_scf_env%matrix_s(1), &
                                       project_out=.TRUE., &
                                       psi_projector_orthogonal=.FALSE., &
                                       proj_in_template=almo_scf_env%matrix_ov(ispin), &
                                       eps_filter=almo_scf_env%eps_filter, &
                                       sig_inv_projector=almo_scf_env%matrix_sigma_inv(ispin))
                  !sig_inv_template=almo_scf_env%matrix_sigma_inv(ispin),&

                  ! save initial retained virtuals
                  CALL dbcsr_create(vr_fixed, &
                                    template=almo_scf_env%matrix_v(ispin))
                  CALL dbcsr_copy(vr_fixed, almo_scf_env%matrix_v(ispin))

                  ! init matrices common for optimized and non-optimized virts
                  CALL dbcsr_create(sigma_vv_sqrt, &
                                    template=almo_scf_env%matrix_sigma_vv(ispin), &
                                    matrix_type=dbcsr_type_no_symmetry)
                  CALL dbcsr_create(sigma_vv_sqrt_inv, &
                                    template=almo_scf_env%matrix_sigma_vv(ispin), &
                                    matrix_type=dbcsr_type_no_symmetry)
                  CALL dbcsr_create(sigma_vv_sqrt_inv_guess, &
                                    template=almo_scf_env%matrix_sigma_vv(ispin), &
                                    matrix_type=dbcsr_type_no_symmetry)
                  CALL dbcsr_create(sigma_vv_sqrt_guess, &
                                    template=almo_scf_env%matrix_sigma_vv(ispin), &
                                    matrix_type=dbcsr_type_no_symmetry)
                  CALL dbcsr_set(sigma_vv_sqrt_guess, 0.0_dp)
                  CALL dbcsr_add_on_diag(sigma_vv_sqrt_guess, 1.0_dp)
                  CALL dbcsr_filter(sigma_vv_sqrt_guess, almo_scf_env%eps_filter)
                  CALL dbcsr_set(sigma_vv_sqrt_inv_guess, 0.0_dp)
                  CALL dbcsr_add_on_diag(sigma_vv_sqrt_inv_guess, 1.0_dp)
                  CALL dbcsr_filter(sigma_vv_sqrt_inv_guess, almo_scf_env%eps_filter)

                  ! do things required to optimize virtuals
                  IF (almo_scf_env%deloc_truncate_virt .NE. virt_full) THEN

                     ! project retained virtuals out of discarded block-by-block
                     ! (1-Q^VR_ALMO)|ALMO_vd>
                     ! this is probably not necessary, do it just to be safe
                     !CALL apply_projector(psi_in=almo_scf_env%matrix_v_disc_blk(ispin),&
                     !        psi_out=almo_scf_env%matrix_v_disc(ispin),&
                     !        psi_projector=almo_scf_env%matrix_v_blk(ispin),&
                     !        metric=almo_scf_env%matrix_s_blk(1),&
                     !        project_out=.TRUE.,&
                     !        psi_projector_orthogonal=.FALSE.,&
                     !        proj_in_template=almo_scf_env%matrix_k_tr(ispin),&
                     !        eps_filter=almo_scf_env%eps_filter,&
                     !        sig_inv_template=almo_scf_env%matrix_sigma_vv(ispin))
                     !CALL dbcsr_copy(almo_scf_env%matrix_v_disc_blk(ispin),&
                     !        almo_scf_env%matrix_v_disc(ispin))

                     ! construct discarded virtuals (1-R)|ALMO_vd>
                     CALL apply_projector(psi_in=almo_scf_env%matrix_v_disc_blk(ispin), &
                                          psi_out=almo_scf_env%matrix_v_disc(ispin), &
                                          psi_projector=almo_scf_env%matrix_t_blk(ispin), &
                                          metric=almo_scf_env%matrix_s(1), &
                                          project_out=.TRUE., &
                                          psi_projector_orthogonal=.FALSE., &
                                          proj_in_template=almo_scf_env%matrix_ov_disc(ispin), &
                                          eps_filter=almo_scf_env%eps_filter, &
                                          sig_inv_projector=almo_scf_env%matrix_sigma_inv(ispin))
                     !sig_inv_template=almo_scf_env%matrix_sigma_inv(ispin),&

                     ! save initial discarded
                     CALL dbcsr_create(vd_fixed, &
                                       template=almo_scf_env%matrix_v_disc(ispin))
                     CALL dbcsr_copy(vd_fixed, almo_scf_env%matrix_v_disc(ispin))

                     !! create the down metric in the retained k-subspace
                     CALL dbcsr_create(k_vr_index_down, &
                                       template=almo_scf_env%matrix_sigma_vv_blk(ispin), &
                                       matrix_type=dbcsr_type_no_symmetry)
                     !CALL dbcsr_copy(k_vr_index_down,&
                     !        almo_scf_env%matrix_sigma_vv_blk(ispin))

                     !CALL get_overlap(bra=almo_scf_env%matrix_v_blk(ispin),&
                     !        ket=almo_scf_env%matrix_v_blk(ispin),&
                     !        overlap=k_vr_index_down,&
                     !        metric=almo_scf_env%matrix_s_blk(1),&
                     !        retain_overlap_sparsity=.FALSE.,&
                     !        eps_filter=almo_scf_env%eps_filter)

                     !! create the up metric in the discarded k-subspace
                     CALL dbcsr_create(k_vd_index_down, &
                                       template=almo_scf_env%matrix_vv_disc_blk(ispin), &
                                       matrix_type=dbcsr_type_no_symmetry)
                     !CALL dbcsr_init(k_vd_index_up)
                     !CALL dbcsr_create(k_vd_index_up,&
                     !        template=almo_scf_env%matrix_vv_disc_blk(ispin),&
                     !        matrix_type=dbcsr_type_no_symmetry)
                     !CALL dbcsr_copy(k_vd_index_down,&
                     !        almo_scf_env%matrix_vv_disc_blk(ispin))

                     !CALL get_overlap(bra=almo_scf_env%matrix_v_disc_blk(ispin),&
                     !        ket=almo_scf_env%matrix_v_disc_blk(ispin),&
                     !        overlap=k_vd_index_down,&
                     !        metric=almo_scf_env%matrix_s_blk(1),&
                     !        retain_overlap_sparsity=.FALSE.,&
                     !        eps_filter=almo_scf_env%eps_filter)

                     !IF (unit_nr>0) THEN
                     !   WRITE(unit_nr,*) "Inverting blocked overlap matrix of discarded virtuals"
                     !ENDIF
                     !CALL invert_Hotelling(k_vd_index_up,&
                     !        k_vd_index_down,&
                     !        almo_scf_env%eps_filter)
                     !IF (safe_mode) THEN
                     !   CALL dbcsr_init(matrix_tmp1)
                     !   CALL dbcsr_create(matrix_tmp1,template=k_vd_index_down,&
                     !                        matrix_type=dbcsr_type_no_symmetry)
                     !   CALL dbcsr_multiply("N","N",1.0_dp,k_vd_index_up,&
                     !                          k_vd_index_down,&
                     !                          0.0_dp, matrix_tmp1,&
                     !                          filter_eps=almo_scf_env%eps_filter)
                     !   frob_matrix_base=dbcsr_frobenius_norm(matrix_tmp1)
                     !   CALL dbcsr_add_on_diag(matrix_tmp1,-1.0_dp)
                     !   frob_matrix=dbcsr_frobenius_norm(matrix_tmp1)
                     !   IF (unit_nr>0) THEN
                     !      WRITE(unit_nr,*) "Error for (inv(SIG)*SIG-I)",&
                     !            frob_matrix/frob_matrix_base
                     !   ENDIF
                     !   CALL dbcsr_release(matrix_tmp1)
                     !ENDIF

                     ! init matrices necessary for optimization of truncated virts
                     ! init blocked gradient before setting K to zero
                     ! otherwise the block structure might be lost
                     CALL dbcsr_create(grad, &
                                       template=almo_scf_env%matrix_k_blk(ispin))
                     CALL dbcsr_copy(grad, almo_scf_env%matrix_k_blk(ispin))

                     ! init MD in the k-space
                     md_in_k_space = almo_scf_env%logical01
                     IF (md_in_k_space) THEN
                        CALL dbcsr_create(velocity, &
                                          template=almo_scf_env%matrix_k_blk(ispin))
                        CALL dbcsr_copy(velocity, almo_scf_env%matrix_k_blk(ispin))
                        CALL dbcsr_set(velocity, 0.0_dp)
                        time_step = almo_scf_env%opt_k_trial_step_size
                     ENDIF

                     CALL dbcsr_create(prev_step, &
                                       template=almo_scf_env%matrix_k_blk(ispin))

                     CALL dbcsr_create(prev_minus_prec_grad, &
                                       template=almo_scf_env%matrix_k_blk(ispin))

                     ! initialize diagonal blocks of the preconditioner to 1.0_dp
                     CALL dbcsr_create(prec, &
                                       template=almo_scf_env%matrix_k_blk(ispin))
                     CALL dbcsr_copy(prec, almo_scf_env%matrix_k_blk(ispin))
                     CALL dbcsr_set(prec, 1.0_dp)

                     ! generate initial K (extrapolate if previous values are available)
                     CALL dbcsr_set(almo_scf_env%matrix_k_blk(ispin), 0.0_dp)
                     ! matrix_k_central stores current k because matrix_k_blk is updated
                     ! during linear search
                     CALL dbcsr_create(matrix_k_central, &
                                       template=almo_scf_env%matrix_k_blk(ispin))
                     CALL dbcsr_copy(matrix_k_central, &
                                     almo_scf_env%matrix_k_blk(ispin))
                     CALL dbcsr_create(tmp_k_blk, &
                                       template=almo_scf_env%matrix_k_blk(ispin))
                     CALL dbcsr_create(step, &
                                       template=almo_scf_env%matrix_k_blk(ispin))
                     CALL dbcsr_set(step, 0.0_dp)
                     CALL dbcsr_create(t_curr, &
                                       template=almo_scf_env%matrix_t(ispin))
                     CALL dbcsr_create(sigma_oo_curr, &
                                       template=almo_scf_env%matrix_sigma(ispin), &
                                       matrix_type=dbcsr_type_no_symmetry)
                     CALL dbcsr_create(sigma_oo_curr_inv, &
                                       template=almo_scf_env%matrix_sigma(ispin), &
                                       matrix_type=dbcsr_type_no_symmetry)
                     CALL dbcsr_create(tmp1_n_vr, &
                                       template=almo_scf_env%matrix_v(ispin))
                     CALL dbcsr_create(tmp3_vd_vr, &
                                       template=almo_scf_env%matrix_k_blk(ispin))
                     CALL dbcsr_create(tmp2_n_o, &
                                       template=almo_scf_env%matrix_t(ispin))
                     CALL dbcsr_create(tmp4_o_vr, &
                                       template=almo_scf_env%matrix_ov(ispin))
                     CALL dbcsr_create(prev_grad, &
                                       template=almo_scf_env%matrix_k_blk(ispin))
                     CALL dbcsr_set(prev_grad, 0.0_dp)

                     !CALL dbcsr_init(sigma_oo_guess)
                     !CALL dbcsr_create(sigma_oo_guess,&
                     !        template=almo_scf_env%matrix_sigma(ispin),&
                     !        matrix_type=dbcsr_type_no_symmetry)
                     !CALL dbcsr_set(sigma_oo_guess,0.0_dp)
                     !CALL dbcsr_add_on_diag(sigma_oo_guess,1.0_dp)
                     !CALL dbcsr_filter(sigma_oo_guess,almo_scf_env%eps_filter)
                     !CALL dbcsr_print(sigma_oo_guess)

                  ENDIF ! done constructing discarded virtuals

                  ! init variables
                  opt_k_max_iter = almo_scf_env%opt_k_max_iter
                  iteration = 0
                  converged = .FALSE.
                  prepare_to_exit = .FALSE.
                  beta = 0.0_dp
                  line_search = .FALSE.
                  obj_function = 0.0_dp
                  conjugacy_error = 0.0_dp
                  line_search_error = 0.0_dp
                  fun0 = 0.0_dp
                  fun1 = 0.0_dp
                  gfun0 = 0.0_dp
                  gfun1 = 0.0_dp
                  step_size_quadratic_approx = 0.0_dp
                  reset_step_size = .TRUE.
                  IF (almo_scf_env%deloc_truncate_virt .EQ. virt_full) opt_k_max_iter = 0

                  ! start cg iterations to optimize matrix_k_blk
                  DO

                     CALL timeset('k_opt_vr', handle1)

                     IF (almo_scf_env%deloc_truncate_virt .NE. virt_full) THEN

                        ! construct k-excited virtuals
                        CALL dbcsr_multiply("N", "N", 1.0_dp, vd_fixed, &
                                            almo_scf_env%matrix_k_blk(ispin), &
                                            0.0_dp, almo_scf_env%matrix_v(ispin), &
                                            filter_eps=almo_scf_env%eps_filter)
                        CALL dbcsr_add(almo_scf_env%matrix_v(ispin), vr_fixed, &
                                       +1.0_dp, +1.0_dp)
                     ENDIF

                     ! decompose the overlap matrix of the current retained orbitals
                     !IF (unit_nr>0) THEN
                     !   WRITE(unit_nr,*) "decompose the active VV overlap matrix"
                     !ENDIF
                     CALL get_overlap(bra=almo_scf_env%matrix_v(ispin), &
                                      ket=almo_scf_env%matrix_v(ispin), &
                                      overlap=almo_scf_env%matrix_sigma_vv(ispin), &
                                      metric=almo_scf_env%matrix_s(1), &
                                      retain_overlap_sparsity=.FALSE., &
                                      eps_filter=almo_scf_env%eps_filter)
                     ! use either cholesky or sqrt
                     !! RZK-warning: strangely, cholesky does not work with k-optimization
                     IF (almo_scf_env%deloc_truncate_virt .EQ. virt_full) THEN
                        CALL timeset('cholesky', handle2)
                        t1cholesky = m_walltime()

                        ! re-create sigma_vv_sqrt because desymmetrize is buggy -
                        ! it will create multiple copies of blocks
                        CALL dbcsr_create(sigma_vv_sqrt, &
                                          template=almo_scf_env%matrix_sigma_vv(ispin), &
                                          matrix_type=dbcsr_type_no_symmetry)
                        CALL dbcsr_desymmetrize(almo_scf_env%matrix_sigma_vv(ispin), &
                                                sigma_vv_sqrt)
                        CALL cp_dbcsr_cholesky_decompose(sigma_vv_sqrt, &
                                                         para_env=almo_scf_env%para_env, &
                                                         blacs_env=almo_scf_env%blacs_env)
                        CALL dbcsr_triu(sigma_vv_sqrt)
                        CALL dbcsr_filter(sigma_vv_sqrt, almo_scf_env%eps_filter)
                        ! apply SOLVE to compute U^(-1) : U*U^(-1)=I
                        CALL dbcsr_get_info(sigma_vv_sqrt, nfullrows_total=n)
                        CALL dbcsr_create(matrix_tmp1, template=almo_scf_env%matrix_sigma_vv(ispin), &
                                          matrix_type=dbcsr_type_no_symmetry)
                        CALL dbcsr_set(matrix_tmp1, 0.0_dp)
                        CALL dbcsr_add_on_diag(matrix_tmp1, 1.0_dp)
                        CALL cp_dbcsr_cholesky_restore(matrix_tmp1, n, sigma_vv_sqrt, &
                                                       sigma_vv_sqrt_inv, op="SOLVE", pos="RIGHT", &
                                                       para_env=almo_scf_env%para_env, &
                                                       blacs_env=almo_scf_env%blacs_env)
                        CALL dbcsr_filter(sigma_vv_sqrt_inv, almo_scf_env%eps_filter)
                        CALL dbcsr_release(matrix_tmp1)
                        IF (safe_mode) THEN
                           CALL dbcsr_create(matrix_tmp1, template=almo_scf_env%matrix_sigma_vv(ispin), &
                                             matrix_type=dbcsr_type_no_symmetry)
                           CALL dbcsr_desymmetrize(almo_scf_env%matrix_sigma_vv(ispin), &
                                                   matrix_tmp1)
                           CALL dbcsr_multiply("T", "N", 1.0_dp, sigma_vv_sqrt, &
                                               sigma_vv_sqrt, &
                                               -1.0_dp, matrix_tmp1, filter_eps=almo_scf_env%eps_filter)
                           frob_matrix = dbcsr_frobenius_norm(matrix_tmp1)
                           CALL dbcsr_add_on_diag(matrix_tmp1, 1.0_dp)
                           frob_matrix_base = dbcsr_frobenius_norm(matrix_tmp1)
                           IF (unit_nr > 0) THEN
                              WRITE (unit_nr, *) "Error for ( U^T * U - Sig )", &
                                 frob_matrix/frob_matrix_base
                           ENDIF
                           CALL dbcsr_multiply("N", "N", 1.0_dp, sigma_vv_sqrt_inv, &
                                               sigma_vv_sqrt, &
                                               0.0_dp, matrix_tmp1, filter_eps=almo_scf_env%eps_filter)
                           frob_matrix_base = dbcsr_frobenius_norm(matrix_tmp1)
                           CALL dbcsr_add_on_diag(matrix_tmp1, -1.0_dp)
                           frob_matrix = dbcsr_frobenius_norm(matrix_tmp1)
                           IF (unit_nr > 0) THEN
                              WRITE (unit_nr, *) "Error for ( inv(U) * U - I )", &
                                 frob_matrix/frob_matrix_base
                           ENDIF
                           CALL dbcsr_release(matrix_tmp1)
                        ENDIF ! safe_mode
                        t2cholesky = m_walltime()
                        IF (unit_nr > 0) THEN
                           WRITE (unit_nr, *) "Cholesky+inverse wall-time: ", t2cholesky - t1cholesky
                        ENDIF
                        CALL timestop(handle2)
                     ELSE
                        CALL matrix_sqrt_Newton_Schulz(sigma_vv_sqrt, &
                                                       sigma_vv_sqrt_inv, &
                                                       almo_scf_env%matrix_sigma_vv(ispin), &
                                                       !matrix_sqrt_inv_guess=sigma_vv_sqrt_inv_guess,&
                                                       !matrix_sqrt_guess=sigma_vv_sqrt_guess,&
                                                       threshold=almo_scf_env%eps_filter, &
                                                       order=almo_scf_env%order_lanczos, &
                                                       eps_lanczos=almo_scf_env%eps_lanczos, &
                                                       max_iter_lanczos=almo_scf_env%max_iter_lanczos)
                        CALL dbcsr_copy(sigma_vv_sqrt_inv_guess, sigma_vv_sqrt_inv)
                        CALL dbcsr_copy(sigma_vv_sqrt_guess, sigma_vv_sqrt)
                        IF (safe_mode) THEN
                           CALL dbcsr_create(matrix_tmp1, template=almo_scf_env%matrix_sigma_vv(ispin), &
                                             matrix_type=dbcsr_type_no_symmetry)
                           CALL dbcsr_create(matrix_tmp2, template=almo_scf_env%matrix_sigma_vv(ispin), &
                                             matrix_type=dbcsr_type_no_symmetry)

                           CALL dbcsr_multiply("N", "N", 1.0_dp, sigma_vv_sqrt_inv, &
                                               almo_scf_env%matrix_sigma_vv(ispin), &
                                               0.0_dp, matrix_tmp1, filter_eps=almo_scf_env%eps_filter)
                           CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, &
                                               sigma_vv_sqrt_inv, &
                                               0.0_dp, matrix_tmp2, filter_eps=almo_scf_env%eps_filter)

                           frob_matrix_base = dbcsr_frobenius_norm(matrix_tmp2)
                           CALL dbcsr_add_on_diag(matrix_tmp2, -1.0_dp)
                           frob_matrix = dbcsr_frobenius_norm(matrix_tmp2)
                           IF (unit_nr > 0) THEN
                              WRITE (unit_nr, *) "Error for (inv(sqrt(SIGVV))*SIGVV*inv(sqrt(SIGVV))-I)", &
                                 frob_matrix/frob_matrix_base
                           ENDIF

                           CALL dbcsr_release(matrix_tmp1)
                           CALL dbcsr_release(matrix_tmp2)
                        ENDIF
                     ENDIF
                     CALL timestop(handle1)

                     ! compute excitation amplitudes (to the current set of retained virtuals)
                     ! set convergence criterion for x-optimization
                     IF ((iteration .EQ. 0) .AND. (.NOT. line_search) .AND. &
                         (outer_opt_k_iteration .EQ. 0)) THEN
                        x_opt_eps_adaptive = &
                           almo_scf_env%deloc_cayley_eps_convergence
                     ELSE
                        x_opt_eps_adaptive = &
                           MAX(ABS(almo_scf_env%deloc_cayley_eps_convergence), &
                               ABS(x_opt_eps_adaptive_factor*grad_norm))
                     ENDIF
                     CALL ct_step_env_init(ct_step_env)
                     CALL ct_step_env_set(ct_step_env, &
                                          para_env=almo_scf_env%para_env, &
                                          blacs_env=almo_scf_env%blacs_env, &
                                          use_occ_orbs=.TRUE., &
                                          use_virt_orbs=.TRUE., &
                                          occ_orbs_orthogonal=.FALSE., &
                                          virt_orbs_orthogonal=.FALSE., &
                                          pp_preconditioner_full=almo_scf_env%deloc_cayley_occ_precond, &
                                          qq_preconditioner_full=almo_scf_env%deloc_cayley_vir_precond, &
                                          tensor_type=almo_scf_env%deloc_cayley_tensor_type, &
                                          neglect_quadratic_term=almo_scf_env%deloc_cayley_linear, &
                                          conjugator=almo_scf_env%deloc_cayley_conjugator, &
                                          max_iter=almo_scf_env%deloc_cayley_max_iter, &
                                          calculate_energy_corr=.TRUE., &
                                          update_p=.FALSE., &
                                          update_q=.FALSE., &
                                          eps_convergence=x_opt_eps_adaptive, &
                                          eps_filter=almo_scf_env%eps_filter, &
                                          !nspins=1,&
                                          q_index_up=sigma_vv_sqrt_inv, &
                                          q_index_down=sigma_vv_sqrt, &
                                          p_index_up=almo_scf_env%matrix_sigma_sqrt_inv(ispin), &
                                          p_index_down=almo_scf_env%matrix_sigma_sqrt(ispin), &
                                          matrix_ks=almo_scf_env%matrix_ks_0deloc(ispin), &
                                          matrix_t=almo_scf_env%matrix_t(ispin), &
                                          matrix_qp_template=almo_scf_env%matrix_vo(ispin), &
                                          matrix_pq_template=almo_scf_env%matrix_ov(ispin), &
                                          matrix_v=almo_scf_env%matrix_v(ispin), &
                                          matrix_x_guess=almo_scf_env%matrix_x(ispin))
                     ! perform calculations
                     CALL ct_step_execute(ct_step_env)
                     ! get the energy correction
                     CALL ct_step_env_get(ct_step_env, &
                                          energy_correction=energy_correction(ispin), &
                                          copy_matrix_x=almo_scf_env%matrix_x(ispin))
                     CALL ct_step_env_clean(ct_step_env)
                     ! RZK-warning matrix_x is being transformed
                     ! back and forth between orth and up_down representations
                     energy_correction(1) = energy_correction(1)*spin_factor

                     IF (opt_k_max_iter .NE. 0) THEN

                        CALL timeset('k_opt_t_curr', handle3)

                        ! construct current occupied orbitals T_blk + V_r*X
                        CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                            almo_scf_env%matrix_v(ispin), &
                                            almo_scf_env%matrix_x(ispin), &
                                            0.0_dp, t_curr, &
                                            filter_eps=almo_scf_env%eps_filter)
                        CALL dbcsr_add(t_curr, almo_scf_env%matrix_t_blk(ispin), &
                                       +1.0_dp, +1.0_dp)

                        ! calculate current occupied overlap
                        !IF (unit_nr>0) THEN
                        !   WRITE(unit_nr,*) "Inverting current occ overlap matrix"
                        !ENDIF
                        CALL get_overlap(bra=t_curr, &
                                         ket=t_curr, &
                                         overlap=sigma_oo_curr, &
                                         metric=almo_scf_env%matrix_s(1), &
                                         retain_overlap_sparsity=.FALSE., &
                                         eps_filter=almo_scf_env%eps_filter)
                        IF (iteration .EQ. 0) THEN
                           CALL invert_Hotelling(sigma_oo_curr_inv, &
                                                 sigma_oo_curr, &
                                                 threshold=almo_scf_env%eps_filter, &
                                                 use_inv_as_guess=.FALSE.)
                        ELSE
                           CALL invert_Hotelling(sigma_oo_curr_inv, &
                                                 sigma_oo_curr, &
                                                 threshold=almo_scf_env%eps_filter, &
                                                 use_inv_as_guess=.TRUE.)
                           !CALL dbcsr_copy(sigma_oo_guess,sigma_oo_curr_inv)
                        ENDIF
                        IF (safe_mode) THEN
                           CALL dbcsr_create(matrix_tmp1, template=sigma_oo_curr, &
                                             matrix_type=dbcsr_type_no_symmetry)
                           CALL dbcsr_multiply("N", "N", 1.0_dp, sigma_oo_curr, &
                                               sigma_oo_curr_inv, &
                                               0.0_dp, matrix_tmp1, &
                                               filter_eps=almo_scf_env%eps_filter)
                           frob_matrix_base = dbcsr_frobenius_norm(matrix_tmp1)
                           CALL dbcsr_add_on_diag(matrix_tmp1, -1.0_dp)
                           frob_matrix = dbcsr_frobenius_norm(matrix_tmp1)
                           !CALL dbcsr_filter(matrix_tmp1,almo_scf_env%eps_filter)
                           !CALL dbcsr_print(matrix_tmp1)
                           IF (unit_nr > 0) THEN
                              WRITE (unit_nr, *) "Error for (SIG*inv(SIG)-I)", &
                                 frob_matrix/frob_matrix_base, frob_matrix_base
                           ENDIF
                           CALL dbcsr_release(matrix_tmp1)
                        ENDIF
                        IF (safe_mode) THEN
                           CALL dbcsr_create(matrix_tmp1, template=sigma_oo_curr, &
                                             matrix_type=dbcsr_type_no_symmetry)
                           CALL dbcsr_multiply("N", "N", 1.0_dp, sigma_oo_curr_inv, &
                                               sigma_oo_curr, &
                                               0.0_dp, matrix_tmp1, &
                                               filter_eps=almo_scf_env%eps_filter)
                           frob_matrix_base = dbcsr_frobenius_norm(matrix_tmp1)
                           CALL dbcsr_add_on_diag(matrix_tmp1, -1.0_dp)
                           frob_matrix = dbcsr_frobenius_norm(matrix_tmp1)
                           !CALL dbcsr_filter(matrix_tmp1,almo_scf_env%eps_filter)
                           !CALL dbcsr_print(matrix_tmp1)
                           IF (unit_nr > 0) THEN
                              WRITE (unit_nr, *) "Error for (inv(SIG)*SIG-I)", &
                                 frob_matrix/frob_matrix_base, frob_matrix_base
                           ENDIF
                           CALL dbcsr_release(matrix_tmp1)
                        ENDIF

                        CALL timestop(handle3)
                        CALL timeset('k_opt_vd', handle4)

                        ! construct current discarded virtuals:
                        ! (1-R_curr)(1-Q^VR_curr)|ALMO_vd_basis> =
                        ! = (1-Q^VR_curr)|ALMO_vd_basis>
                        ! use sigma_vv_sqrt to store the inverse of the overlap
                        ! sigma_vv_inv is computed from sqrt/cholesky
                        CALL dbcsr_multiply("N", "T", 1.0_dp, &
                                            sigma_vv_sqrt_inv, &
                                            sigma_vv_sqrt_inv, &
                                            0.0_dp, sigma_vv_sqrt, &
                                            filter_eps=almo_scf_env%eps_filter)
                        CALL apply_projector(psi_in=almo_scf_env%matrix_v_disc_blk(ispin), &
                                             psi_out=almo_scf_env%matrix_v_disc(ispin), &
                                             psi_projector=almo_scf_env%matrix_v(ispin), &
                                             metric=almo_scf_env%matrix_s(1), &
                                             project_out=.FALSE., &
                                             psi_projector_orthogonal=.FALSE., &
                                             proj_in_template=almo_scf_env%matrix_k_tr(ispin), &
                                             eps_filter=almo_scf_env%eps_filter, &
                                             sig_inv_projector=sigma_vv_sqrt)
                        !sig_inv_template=almo_scf_env%matrix_sigma_vv(ispin),&
                        CALL dbcsr_add(almo_scf_env%matrix_v_disc(ispin), &
                                       vd_fixed, -1.0_dp, +1.0_dp)

                        CALL timestop(handle4)
                        CALL timeset('k_opt_grad', handle5)

                        ! evaluate the gradient from the assembled components
                        ! grad_xx = c0 [ (Vd_curr^tr)*F*T_curr*sigma_oo_curr_inv*(X^tr)]_xx
                        ! save previous gradient to calculate conjugation coef
                        IF (line_search) THEN
                           CALL dbcsr_copy(prev_grad, grad)
                        ENDIF
                        CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                            almo_scf_env%matrix_ks_0deloc(ispin), &
                                            t_curr, &
                                            0.0_dp, tmp2_n_o, &
                                            filter_eps=almo_scf_env%eps_filter)
                        CALL dbcsr_multiply("N", "T", 1.0_dp, &
                                            sigma_oo_curr_inv, &
                                            almo_scf_env%matrix_x(ispin), &
                                            0.0_dp, tmp4_o_vr, &
                                            filter_eps=almo_scf_env%eps_filter)
                        CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                            tmp2_n_o, &
                                            tmp4_o_vr, &
                                            0.0_dp, tmp1_n_vr, &
                                            filter_eps=almo_scf_env%eps_filter)
                        CALL dbcsr_multiply("T", "N", 2.0_dp*spin_factor, &
                                            almo_scf_env%matrix_v_disc(ispin), &
                                            tmp1_n_vr, &
                                            0.0_dp, grad, &
                                            retain_sparsity=.TRUE.)
                        !filter_eps=almo_scf_env%eps_filter,&
                        ! keep tmp2_n_o for the next step
                        ! keep tmp4_o_vr for the preconditioner

                        ! check convergence and other exit criteria
                        grad_norm_frob = dbcsr_frobenius_norm(grad)
                        CALL dbcsr_norm(grad, dbcsr_norm_maxabsnorm, norm_scalar=grad_norm)
                        converged = (grad_norm .LT. almo_scf_env%opt_k_eps_convergence)
                        IF (converged .OR. (iteration .GE. opt_k_max_iter)) THEN
                           prepare_to_exit = .TRUE.
                        ENDIF
                        CALL timestop(handle5)

                        IF (.NOT. prepare_to_exit) THEN

                           CALL timeset('k_opt_energy', handle6)

                           ! compute "energy" c0*Tr[sig_inv_oo*t*F*t]
                           CALL dbcsr_multiply("T", "N", spin_factor, &
                                               t_curr, &
                                               tmp2_n_o, &
                                               0.0_dp, sigma_oo_curr, &
                                               filter_eps=almo_scf_env%eps_filter)
                           delta_obj_function = fun0
                           CALL dbcsr_dot(sigma_oo_curr_inv, sigma_oo_curr, obj_function)
                           delta_obj_function = obj_function - delta_obj_function
                           IF (line_search) THEN
                              fun1 = obj_function
                           ELSE
                              fun0 = obj_function
                           ENDIF

                           CALL timestop(handle6)

                           ! update the step direction
                           IF (.NOT. line_search) THEN

                              CALL timeset('k_opt_step', handle7)

                              IF ((.NOT. md_in_k_space) .AND. &
                                  (iteration .GE. MAX(0, almo_scf_env%opt_k_prec_iter_start) .AND. &
                                   MOD(iteration - almo_scf_env%opt_k_prec_iter_start, &
                                       almo_scf_env%opt_k_prec_iter_freq) .EQ. 0)) THEN

                                 !IF ((iteration.eq.0).AND.(.NOT.md_in_k_space)) THEN

                                 ! compute the preconditioner
                                 IF (unit_nr > 0) THEN
                                    WRITE (unit_nr, *) "Computing preconditioner"
                                 ENDIF
                                 !CALL opt_k_create_preconditioner(prec,&
                                 !        almo_scf_env%matrix_v_disc(ispin),&
                                 !        almo_scf_env%matrix_ks_0deloc(ispin),&
                                 !        almo_scf_env%matrix_x(ispin),&
                                 !        tmp4_o_vr,&
                                 !        almo_scf_env%matrix_s(1),&
                                 !        grad,&
                                 !        !almo_scf_env%matrix_v_disc_blk(ispin),&
                                 !        vd_fixed,&
                                 !        t_curr,&
                                 !        k_vd_index_up,&
                                 !        k_vr_index_down,&
                                 !        tmp1_n_vr,&
                                 !        spin_factor,&
                                 !        almo_scf_env%eps_filter)
                                 CALL opt_k_create_preconditioner_blk(almo_scf_env, &
                                                                      almo_scf_env%matrix_v_disc(ispin), &
                                                                      tmp4_o_vr, &
                                                                      t_curr, &
                                                                      ispin, &
                                                                      spin_factor)

                              ENDIF

                              ! save the previous step
                              CALL dbcsr_copy(prev_step, step)

                              ! compute the new step
                              CALL opt_k_apply_preconditioner_blk(almo_scf_env, &
                                                                  step, grad, ispin)
                              !CALL dbcsr_hadamard_product(prec,grad,step)
                              CALL dbcsr_scale(step, -1.0_dp)

                              ! check whether we need to reset conjugate directions
                              reset_conjugator = .FALSE.
                              ! first check if manual reset is active
                              IF (iteration .LT. MAX(almo_scf_env%opt_k_conj_iter_start, 1) .OR. &
                                  MOD(iteration - almo_scf_env%opt_k_conj_iter_start, &
                                      almo_scf_env%opt_k_conj_iter_freq) .EQ. 0) THEN

                                 reset_conjugator = .TRUE.

                              ELSE

                                 ! check for the errors in the cg algorithm
                                 !CALL dbcsr_hadamard_product(prec,prev_grad,tmp_k_blk)
                                 !CALL dbcsr_dot(grad,tmp_k_blk,numer)
                                 !CALL dbcsr_dot(prev_grad,tmp_k_blk,denom)
                                 CALL dbcsr_dot(grad, prev_minus_prec_grad, numer)
                                 CALL dbcsr_dot(prev_grad, prev_minus_prec_grad, denom)
                                 conjugacy_error = numer/denom

                                 IF (conjugacy_error .GT. MIN(0.5_dp, conjugacy_error_threshold)) THEN
                                    reset_conjugator = .TRUE.
                                    IF (unit_nr > 0) THEN
                                       WRITE (unit_nr, *) "Lack of progress, conjugacy error is ", conjugacy_error
                                    ENDIF
                                 ENDIF

                                 ! check the gradient along the previous direction
                                 IF ((iteration .NE. 0) .AND. (.NOT. reset_conjugator)) THEN
                                    CALL dbcsr_dot(grad, prev_step, numer)
                                    CALL dbcsr_dot(prev_grad, prev_step, denom)
                                    line_search_error = numer/denom
                                    IF (line_search_error .GT. line_search_error_threshold) THEN
                                       reset_conjugator = .TRUE.
                                       IF (unit_nr > 0) THEN
                                          WRITE (unit_nr, *) "Bad line search, line search error is ", line_search_error
                                       ENDIF
                                    ENDIF
                                 ENDIF

                              ENDIF

                              ! compute the conjugation coefficient - beta
                              IF (.NOT. reset_conjugator) THEN

                                 SELECT CASE (almo_scf_env%opt_k_conjugator)
                                 CASE (cg_hestenes_stiefel)
                                    CALL dbcsr_copy(tmp_k_blk, grad)
                                    CALL dbcsr_add(tmp_k_blk, prev_grad, 1.0_dp, -1.0_dp)
                                    CALL dbcsr_dot(tmp_k_blk, step, numer)
                                    CALL dbcsr_dot(tmp_k_blk, prev_step, denom)
                                    beta = -1.0_dp*numer/denom
                                 CASE (cg_fletcher_reeves)
                                    !CALL dbcsr_hadamard_product(prec,prev_grad,tmp_k_blk)
                                    !CALL dbcsr_dot(prev_grad,tmp_k_blk,denom)
                                    !CALL dbcsr_hadamard_product(prec,grad,tmp_k_blk)
                                    !CALL dbcsr_dot(grad,tmp_k_blk,numer)
                                    !beta=numer/denom
                                    CALL dbcsr_dot(grad, step, numer)
                                    CALL dbcsr_dot(prev_grad, prev_minus_prec_grad, denom)
                                    beta = numer/denom
                                 CASE (cg_polak_ribiere)
                                    !CALL dbcsr_hadamard_product(prec,prev_grad,tmp_k_blk)
                                    !CALL dbcsr_dot(prev_grad,tmp_k_blk,denom)
                                    !CALL dbcsr_add(prev_grad,grad,-1.0_dp,1.0_dp)
                                    !CALL dbcsr_hadamard_product(prec,prev_grad,tmp_k_blk)
                                    !CALL dbcsr_dot(tmp_k_blk,grad,numer)
                                    CALL dbcsr_dot(prev_grad, prev_minus_prec_grad, denom)
                                    CALL dbcsr_copy(tmp_k_blk, grad)
                                    CALL dbcsr_add(tmp_k_blk, prev_grad, 1.0_dp, -1.0_dp)
                                    CALL dbcsr_dot(tmp_k_blk, step, numer)
                                    beta = numer/denom
                                 CASE (cg_fletcher)
                                    !CALL dbcsr_hadamard_product(prec,grad,tmp_k_blk)
                                    !CALL dbcsr_dot(grad,tmp_k_blk,numer)
                                    !CALL dbcsr_dot(prev_grad,prev_step,denom)
                                    !beta=-1.0_dp*numer/denom
                                    CALL dbcsr_dot(grad, step, numer)
                                    CALL dbcsr_dot(prev_grad, prev_step, denom)
                                    beta = numer/denom
                                 CASE (cg_liu_storey)
                                    CALL dbcsr_dot(prev_grad, prev_step, denom)
                                    !CALL dbcsr_add(prev_grad,grad,-1.0_dp,1.0_dp)
                                    !CALL dbcsr_hadamard_product(prec,prev_grad,tmp_k_blk)
                                    !CALL dbcsr_dot(tmp_k_blk,grad,numer)
                                    CALL dbcsr_copy(tmp_k_blk, grad)
                                    CALL dbcsr_add(tmp_k_blk, prev_grad, 1.0_dp, -1.0_dp)
                                    CALL dbcsr_dot(tmp_k_blk, step, numer)
                                    beta = numer/denom
                                 CASE (cg_dai_yuan)
                                    !CALL dbcsr_hadamard_product(prec,grad,tmp_k_blk)
                                    !CALL dbcsr_dot(grad,tmp_k_blk,numer)
                                    !CALL dbcsr_add(prev_grad,grad,-1.0_dp,1.0_dp)
                                    !CALL dbcsr_dot(prev_grad,prev_step,denom)
                                    !beta=numer/denom
                                    CALL dbcsr_dot(grad, step, numer)
                                    CALL dbcsr_copy(tmp_k_blk, grad)
                                    CALL dbcsr_add(tmp_k_blk, prev_grad, 1.0_dp, -1.0_dp)
                                    CALL dbcsr_dot(tmp_k_blk, prev_step, denom)
                                    beta = -1.0_dp*numer/denom
                                 CASE (cg_hager_zhang)
                                    !CALL dbcsr_add(prev_grad,grad,-1.0_dp,1.0_dp)
                                    !CALL dbcsr_dot(prev_grad,prev_step,denom)
                                    !CALL dbcsr_hadamard_product(prec,prev_grad,tmp_k_blk)
                                    !CALL dbcsr_dot(tmp_k_blk,prev_grad,numer)
                                    !kappa=2.0_dp*numer/denom
                                    !CALL dbcsr_dot(tmp_k_blk,grad,numer)
                                    !tau=numer/denom
                                    !CALL dbcsr_dot(prev_step,grad,numer)
                                    !beta=tau-kappa*numer/denom
                                    CALL dbcsr_copy(tmp_k_blk, grad)
                                    CALL dbcsr_add(tmp_k_blk, prev_grad, 1.0_dp, -1.0_dp)
                                    CALL dbcsr_dot(tmp_k_blk, prev_step, denom)
                                    CALL dbcsr_dot(tmp_k_blk, prev_minus_prec_grad, numer)
                                    kappa = -2.0_dp*numer/denom
                                    CALL dbcsr_dot(tmp_k_blk, step, numer)
                                    tau = -1.0_dp*numer/denom
                                    CALL dbcsr_dot(prev_step, grad, numer)
                                    beta = tau - kappa*numer/denom
                                 CASE (cg_zero)
                                    beta = 0.0_dp
                                 CASE DEFAULT
                                    CPABORT("illegal conjugator")
                                 END SELECT

                                 IF (beta .LT. 0.0_dp) THEN
                                    IF (unit_nr > 0) THEN
                                       WRITE (unit_nr, *) "Beta is negative, ", beta
                                    ENDIF
                                    reset_conjugator = .TRUE.
                                 ENDIF

                              ENDIF

                              IF (md_in_k_space) THEN
                                 reset_conjugator = .TRUE.
                              ENDIF

                              IF (reset_conjugator) THEN

                                 beta = 0.0_dp
                                 !reset_step_size=.TRUE.

                                 IF (unit_nr > 0) THEN
                                    WRITE (unit_nr, *) "(Re)-setting conjugator to zero"
                                 ENDIF

                              ENDIF

                              ! save the preconditioned gradient
                              CALL dbcsr_copy(prev_minus_prec_grad, step)

                              ! conjugate the step direction
                              CALL dbcsr_add(step, prev_step, 1.0_dp, beta)

                              CALL timestop(handle7)

                              ! update the step direction
                           ELSE ! step update
                              conjugacy_error = 0.0_dp
                           ENDIF

                           ! compute the gradient with respect to the step size in the curr direction
                           IF (line_search) THEN
                              CALL dbcsr_dot(grad, step, gfun1)
                              line_search_error = gfun1/gfun0
                           ELSE
                              CALL dbcsr_dot(grad, step, gfun0)
                           ENDIF

                           ! make a step - update k
                           IF (line_search) THEN

                              ! check if the trial step provides enough numerical accuracy
                              safety_multiplier = 1.0E+1_dp ! must be more than one
                              num_threshold = MAX(EPSILON(1.0_dp), &
                                                  safety_multiplier*(almo_scf_env%eps_filter**2)*almo_scf_env%ndomains)
                              IF (ABS(fun1 - fun0 - gfun0*step_size) .LT. num_threshold) THEN
                                 IF (unit_nr > 0) THEN
                                    WRITE (unit_nr, '(T3,A,1X,E17.7)') &
                                       "Numerical accuracy is too low to observe non-linear behavior", &
                                       ABS(fun1 - fun0 - gfun0*step_size)
                                    WRITE (unit_nr, '(T3,A,1X,E17.7,A,1X,E12.3)') "Error computing ", &
                                       ABS(gfun0), &
                                       " is smaller than the threshold", num_threshold
                                 ENDIF
                                 CPABORT("")
                              ENDIF
                              IF (ABS(gfun0) .LT. num_threshold) THEN
                                 IF (unit_nr > 0) THEN
                                    WRITE (unit_nr, '(T3,A,1X,E17.7,A,1X,E12.3)') "Linear gradient", &
                                       ABS(gfun0), &
                                       " is smaller than the threshold", num_threshold
                                 ENDIF
                                 CPABORT("")
                              ENDIF

                              use_quadratic_approximation = .TRUE.
                              use_cubic_approximation = .FALSE.

                              ! find the minimum assuming quadratic form
                              ! use f0, f1, g0
                              step_size_quadratic_approx = -(gfun0*step_size*step_size)/(2.0_dp*(fun1 - fun0 - gfun0*step_size))
                              ! use f0, f1, g1
                             step_size_quadratic_approx2 = -(fun1 - fun0 - step_size*gfun1/2.0_dp)/(gfun1 - (fun1 - fun0)/step_size)

                              IF ((step_size_quadratic_approx .LT. 0.0_dp) .AND. &
                                  (step_size_quadratic_approx2 .LT. 0.0_dp)) THEN
                                 IF (unit_nr > 0) THEN
                                    WRITE (unit_nr, '(T3,A,1X,E17.7,1X,E17.7,1X,A)') &
                                       "Quadratic approximation gives negative steps", &
                                       step_size_quadratic_approx, step_size_quadratic_approx2, &
                                       "trying cubic..."
                                 ENDIF
                                 use_cubic_approximation = .TRUE.
                                 use_quadratic_approximation = .FALSE.
                              ELSE
                                 IF (step_size_quadratic_approx .LT. 0.0_dp) THEN
                                    step_size_quadratic_approx = step_size_quadratic_approx2
                                 ENDIF
                                 IF (step_size_quadratic_approx2 .LT. 0.0_dp) THEN
                                    step_size_quadratic_approx2 = step_size_quadratic_approx
                                 ENDIF
                              ENDIF

                              ! check accuracy of the quadratic approximation
                              IF (use_quadratic_approximation) THEN
                                 quadratic_approx_error = ABS(step_size_quadratic_approx - &
                                                              step_size_quadratic_approx2)/step_size_quadratic_approx
                                 IF (quadratic_approx_error .GT. quadratic_approx_error_threshold) THEN
                                    IF (unit_nr > 0) THEN
                                       WRITE (unit_nr, '(T3,A,1X,E17.7,1X,E17.7,1X,A)') "Quadratic approximation is poor", &
                                          step_size_quadratic_approx, step_size_quadratic_approx2, &
                                          "Try cubic approximation"
                                    ENDIF
                                    use_cubic_approximation = .TRUE.
                                    use_quadratic_approximation = .FALSE.
                                 ENDIF
                              ENDIF

                              ! check if numerics is fine enough to capture the cubic form
                              IF (use_cubic_approximation) THEN

                                 ! if quadratic approximation is not accurate enough
                                 ! try to find the minimum assuming cubic form
                                 ! aa*x**3 + bb*x**2 + cc*x + dd = f(x)
                                 bb = (-step_size*gfun1 + 3.0_dp*(fun1 - fun0) - 2.0_dp*step_size*gfun0)/(step_size*step_size)
                                 aa = (gfun1 - 2.0_dp*step_size*bb - gfun0)/(3.0_dp*step_size*step_size)

                                 IF (ABS(gfun1 - 2.0_dp*step_size*bb - gfun0) .LT. num_threshold) THEN
                                    IF (unit_nr > 0) THEN
                                       WRITE (unit_nr, '(T3,A,1X,E17.7)') &
                                          "Numerical accuracy is too low to observe cubic behavior", &
                                          ABS(gfun1 - 2.0_dp*step_size*bb - gfun0)
                                    ENDIF
                                    use_cubic_approximation = .FALSE.
                                    use_quadratic_approximation = .TRUE.
                                 ENDIF
                                 IF (ABS(gfun1) .LT. num_threshold) THEN
                                    IF (unit_nr > 0) THEN
                                       WRITE (unit_nr, '(T3,A,1X,E17.7,A,1X,E12.3)') "Linear gradient", &
                                          ABS(gfun1), &
                                          " is smaller than the threshold", num_threshold
                                    ENDIF
                                    use_cubic_approximation = .FALSE.
                                    use_quadratic_approximation = .TRUE.
                                 ENDIF
                              ENDIF

                              ! find the step assuming cubic approximation
                              IF (use_cubic_approximation) THEN
                                 ! to obtain the minimum of the cubic function solve the quadratic equation
                                 ! 0.0*x**3 + 3.0*aa*x**2 + 2.0*bb*x + cc = 0
                                 CALL analytic_line_search(0.0_dp, 3.0_dp*aa, 2.0_dp*bb, gfun0, minima, nmins)
                                 IF (nmins .LT. 1) THEN
                                    IF (unit_nr > 0) THEN
                                       WRITE (unit_nr, '(T3,A)') &
                                          "Cubic approximation gives zero soultions! Use quadratic approximation"
                                    ENDIF
                                    use_quadratic_approximation = .TRUE.
                                    use_cubic_approximation = .TRUE.
                                 ELSE
                                    step_size = minima(1)
                                    IF (nmins .GT. 1) THEN
                                       IF (unit_nr > 0) THEN
                                          WRITE (unit_nr, '(T3,A)') &
                                             "More than one solution found! Use quadratic approximation"
                                       ENDIF
                                       use_quadratic_approximation = .TRUE.
                                       use_cubic_approximation = .TRUE.
                                    ENDIF
                                 ENDIF
                              ENDIF

                              IF (use_quadratic_approximation) THEN ! use quadratic approximation
                                 IF (unit_nr > 0) THEN
                                    WRITE (unit_nr, '(T3,A)') "Use quadratic approximation"
                                 ENDIF
                                 step_size = (step_size_quadratic_approx + step_size_quadratic_approx2)*0.5_dp
                              ENDIF

                              ! one more check on the step size
                              IF (step_size .LT. 0.0_dp) THEN
                                 CPABORT("Negative step proposed")
                              ENDIF

                              CALL dbcsr_copy(almo_scf_env%matrix_k_blk(ispin), &
                                              matrix_k_central)
                              CALL dbcsr_add(almo_scf_env%matrix_k_blk(ispin), &
                                             step, 1.0_dp, step_size)
                              CALL dbcsr_copy(matrix_k_central, &
                                              almo_scf_env%matrix_k_blk(ispin))
                              line_search = .FALSE.

                           ELSE

                              IF (md_in_k_space) THEN

                                 ! update velocities v(i) = v(i-1) + 0.5*dT*(a(i-1) + a(i))
                                 IF (iteration .NE. 0) THEN
                                    CALL dbcsr_add(velocity, &
                                                   step, 1.0_dp, 0.5_dp*time_step)
                                    CALL dbcsr_add(velocity, &
                                                   prev_step, 1.0_dp, 0.5_dp*time_step)
                                 ENDIF
                                 kin_energy = dbcsr_frobenius_norm(velocity)
                                 kin_energy = 0.5_dp*kin_energy*kin_energy

                                 ! update positions k(i) = k(i-1) + dT*v(i-1) + 0.5*dT*dT*a(i-1)
                                 CALL dbcsr_add(almo_scf_env%matrix_k_blk(ispin), &
                                                velocity, 1.0_dp, time_step)
                                 CALL dbcsr_add(almo_scf_env%matrix_k_blk(ispin), &
                                                step, 1.0_dp, 0.5_dp*time_step*time_step)

                              ELSE

                                 IF (reset_step_size) THEN
                                    step_size = almo_scf_env%opt_k_trial_step_size
                                    reset_step_size = .FALSE.
                                 ELSE
                                    step_size = step_size*almo_scf_env%opt_k_trial_step_size_multiplier
                                 ENDIF
                                 CALL dbcsr_copy(almo_scf_env%matrix_k_blk(ispin), &
                                                 matrix_k_central)
                                 CALL dbcsr_add(almo_scf_env%matrix_k_blk(ispin), &
                                                step, 1.0_dp, step_size)
                                 line_search = .TRUE.
                              ENDIF

                           ENDIF

                        ENDIF ! .NOT.prepare_to_exit

                        ! print the status of the optimization
                        t2a = m_walltime()
                        IF (unit_nr > 0) THEN
                           IF (md_in_k_space) THEN
                              WRITE (unit_nr, '(T6,A,1X,I5,1X,E12.3,E16.7,F15.9,F15.9,F15.9,E12.3,F15.9,F15.9,F8.3)') &
                                 "K iter CG", iteration, time_step, time_step*iteration, &
                                 energy_correction(ispin), obj_function, delta_obj_function, grad_norm, &
                                 kin_energy, kin_energy + obj_function, beta
                           ELSE
                              IF (line_search .OR. prepare_to_exit) THEN
                                 WRITE (unit_nr, '(T6,A,1X,I3,1X,E12.3,F16.10,F16.10,E12.3,E12.3,E12.3,F8.3,F8.3,F10.3)') &
                                    "K iter CG", iteration, step_size, &
                                    energy_correction(ispin), delta_obj_function, grad_norm, &
                                    gfun0, line_search_error, beta, conjugacy_error, t2a - t1a
                                 !(flop1+flop2)/(1.0E6_dp*(t2-t1))
                              ELSE
                                 WRITE (unit_nr, '(T6,A,1X,I3,1X,E12.3,F16.10,F16.10,E12.3,E12.3,E12.3,F8.3,F8.3,F10.3)') &
                                    "K iter LS", iteration, step_size, &
                                    energy_correction(ispin), delta_obj_function, grad_norm, &
                                    gfun1, line_search_error, beta, conjugacy_error, t2a - t1a
                                 !(flop1+flop2)/(1.0E6_dp*(t2-t1))
                              ENDIF
                           ENDIF
                           CALL m_flush(unit_nr)
                        ENDIF
                        t1a = m_walltime()

                     ELSE ! opt_k_max_iter .eq. 0
                        prepare_to_exit = .TRUE.
                     ENDIF ! opt_k_max_iter .ne. 0

                     IF (.NOT. line_search) iteration = iteration + 1

                     IF (prepare_to_exit) EXIT

                  ENDDO ! end iterations on K

                  IF (converged .OR. (outer_opt_k_iteration .GE. outer_opt_k_max_iter)) THEN
                     outer_opt_k_prepare_to_exit = .TRUE.
                  ENDIF

                  IF (almo_scf_env%deloc_truncate_virt .NE. virt_full) THEN

                     IF (unit_nr > 0) THEN
                        WRITE (unit_nr, *) "Updating ALMO virtuals"
                     ENDIF

                     CALL timeset('k_opt_v0_update', handle8)

                     ! update retained ALMO virtuals to restart the cg iterations
                     CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                         almo_scf_env%matrix_v_disc_blk(ispin), &
                                         almo_scf_env%matrix_k_blk(ispin), &
                                         0.0_dp, vr_fixed, &
                                         filter_eps=almo_scf_env%eps_filter)
                     CALL dbcsr_add(vr_fixed, almo_scf_env%matrix_v_blk(ispin), &
                                    +1.0_dp, +1.0_dp)

                     ! update discarded ALMO virtuals to restart the cg iterations
                     CALL dbcsr_multiply("N", "T", 1.0_dp, &
                                         almo_scf_env%matrix_v_blk(ispin), &
                                         almo_scf_env%matrix_k_blk(ispin), &
                                         0.0_dp, vd_fixed, &
                                         filter_eps=almo_scf_env%eps_filter)
                     CALL dbcsr_add(vd_fixed, almo_scf_env%matrix_v_disc_blk(ispin), &
                                    -1.0_dp, +1.0_dp)

                     ! orthogonalize new orbitals on fragments
                     CALL get_overlap(bra=vr_fixed, &
                                      ket=vr_fixed, &
                                      overlap=k_vr_index_down, &
                                      metric=almo_scf_env%matrix_s_blk(1), &
                                      retain_overlap_sparsity=.FALSE., &
                                      eps_filter=almo_scf_env%eps_filter)
                     CALL dbcsr_create(vr_index_sqrt_inv, template=k_vr_index_down, &
                                       matrix_type=dbcsr_type_no_symmetry)
                     CALL dbcsr_create(vr_index_sqrt, template=k_vr_index_down, &
                                       matrix_type=dbcsr_type_no_symmetry)
                     CALL matrix_sqrt_Newton_Schulz(vr_index_sqrt, &
                                                    vr_index_sqrt_inv, &
                                                    k_vr_index_down, &
                                                    threshold=almo_scf_env%eps_filter, &
                                                    order=almo_scf_env%order_lanczos, &
                                                    eps_lanczos=almo_scf_env%eps_lanczos, &
                                                    max_iter_lanczos=almo_scf_env%max_iter_lanczos)
                     IF (safe_mode) THEN
                        CALL dbcsr_create(matrix_tmp1, template=k_vr_index_down, &
                                          matrix_type=dbcsr_type_no_symmetry)
                        CALL dbcsr_create(matrix_tmp2, template=k_vr_index_down, &
                                          matrix_type=dbcsr_type_no_symmetry)

                        CALL dbcsr_multiply("N", "N", 1.0_dp, vr_index_sqrt_inv, &
                                            k_vr_index_down, &
                                            0.0_dp, matrix_tmp1, filter_eps=almo_scf_env%eps_filter)
                        CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, &
                                            vr_index_sqrt_inv, &
                                            0.0_dp, matrix_tmp2, filter_eps=almo_scf_env%eps_filter)

                        frob_matrix_base = dbcsr_frobenius_norm(matrix_tmp2)
                        CALL dbcsr_add_on_diag(matrix_tmp2, -1.0_dp)
                        frob_matrix = dbcsr_frobenius_norm(matrix_tmp2)
                        IF (unit_nr > 0) THEN
                           WRITE (unit_nr, *) "Error for (inv(sqrt(SIGVV))*SIGVV*inv(sqrt(SIGVV))-I)", &
                              frob_matrix/frob_matrix_base
                        ENDIF

                        CALL dbcsr_release(matrix_tmp1)
                        CALL dbcsr_release(matrix_tmp2)
                     ENDIF
                     CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                         vr_fixed, &
                                         vr_index_sqrt_inv, &
                                         0.0_dp, almo_scf_env%matrix_v_blk(ispin), &
                                         filter_eps=almo_scf_env%eps_filter)

                     CALL get_overlap(bra=vd_fixed, &
                                      ket=vd_fixed, &
                                      overlap=k_vd_index_down, &
                                      metric=almo_scf_env%matrix_s_blk(1), &
                                      retain_overlap_sparsity=.FALSE., &
                                      eps_filter=almo_scf_env%eps_filter)
                     CALL dbcsr_create(vd_index_sqrt_inv, template=k_vd_index_down, &
                                       matrix_type=dbcsr_type_no_symmetry)
                     CALL dbcsr_create(vd_index_sqrt, template=k_vd_index_down, &
                                       matrix_type=dbcsr_type_no_symmetry)
                     CALL matrix_sqrt_Newton_Schulz(vd_index_sqrt, &
                                                    vd_index_sqrt_inv, &
                                                    k_vd_index_down, &
                                                    threshold=almo_scf_env%eps_filter, &
                                                    order=almo_scf_env%order_lanczos, &
                                                    eps_lanczos=almo_scf_env%eps_lanczos, &
                                                    max_iter_lanczos=almo_scf_env%max_iter_lanczos)
                     IF (safe_mode) THEN
                        CALL dbcsr_create(matrix_tmp1, template=k_vd_index_down, &
                                          matrix_type=dbcsr_type_no_symmetry)
                        CALL dbcsr_create(matrix_tmp2, template=k_vd_index_down, &
                                          matrix_type=dbcsr_type_no_symmetry)

                        CALL dbcsr_multiply("N", "N", 1.0_dp, vd_index_sqrt_inv, &
                                            k_vd_index_down, &
                                            0.0_dp, matrix_tmp1, filter_eps=almo_scf_env%eps_filter)
                        CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, &
                                            vd_index_sqrt_inv, &
                                            0.0_dp, matrix_tmp2, filter_eps=almo_scf_env%eps_filter)

                        frob_matrix_base = dbcsr_frobenius_norm(matrix_tmp2)
                        CALL dbcsr_add_on_diag(matrix_tmp2, -1.0_dp)
                        frob_matrix = dbcsr_frobenius_norm(matrix_tmp2)
                        IF (unit_nr > 0) THEN
                           WRITE (unit_nr, *) "Error for (inv(sqrt(SIGVV))*SIGVV*inv(sqrt(SIGVV))-I)", &
                              frob_matrix/frob_matrix_base
                        ENDIF

                        CALL dbcsr_release(matrix_tmp1)
                        CALL dbcsr_release(matrix_tmp2)
                     ENDIF
                     CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                         vd_fixed, &
                                         vd_index_sqrt_inv, &
                                         0.0_dp, almo_scf_env%matrix_v_disc_blk(ispin), &
                                         filter_eps=almo_scf_env%eps_filter)

                     CALL dbcsr_release(vr_index_sqrt_inv)
                     CALL dbcsr_release(vr_index_sqrt)
                     CALL dbcsr_release(vd_index_sqrt_inv)
                     CALL dbcsr_release(vd_index_sqrt)

                     CALL timestop(handle8)

                  ENDIF ! ne.virt_full

                  ! RZK-warning released outside the outer loop
                  CALL dbcsr_release(sigma_vv_sqrt)
                  CALL dbcsr_release(sigma_vv_sqrt_inv)
                  IF (almo_scf_env%deloc_truncate_virt .NE. virt_full) THEN
                     CALL dbcsr_release(k_vr_index_down)
                     CALL dbcsr_release(k_vd_index_down)
                     !CALL dbcsr_release(k_vd_index_up)
                     CALL dbcsr_release(matrix_k_central)
                     CALL dbcsr_release(vr_fixed)
                     CALL dbcsr_release(vd_fixed)
                     CALL dbcsr_release(grad)
                     CALL dbcsr_release(prec)
                     CALL dbcsr_release(prev_grad)
                     CALL dbcsr_release(tmp3_vd_vr)
                     CALL dbcsr_release(tmp1_n_vr)
                     CALL dbcsr_release(tmp_k_blk)
                     CALL dbcsr_release(t_curr)
                     CALL dbcsr_release(sigma_oo_curr)
                     CALL dbcsr_release(sigma_oo_curr_inv)
                     CALL dbcsr_release(step)
                     CALL dbcsr_release(tmp2_n_o)
                     CALL dbcsr_release(tmp4_o_vr)
                     CALL dbcsr_release(prev_step)
                     CALL dbcsr_release(prev_minus_prec_grad)
                     IF (md_in_k_space) THEN
                        CALL dbcsr_release(velocity)
                     ENDIF

                  ENDIF

                  outer_opt_k_iteration = outer_opt_k_iteration + 1
                  IF (outer_opt_k_prepare_to_exit) EXIT

               ENDDO ! outer loop for k

            ENDDO ! ispin

            ! RZK-warning update mo orbitals

         ELSE ! virtual orbitals might not be available use projected AOs

            ! compute sqrt(S) and inv(sqrt(S))
            ! RZK-warning - remove this sqrt(S) and inv(sqrt(S))
            ! ideally ALMO scf should use sigma and sigma_inv in
            ! the tensor_up_down representation
            IF (.NOT. almo_scf_env%s_sqrt_done) THEN

               IF (unit_nr > 0) THEN
                  WRITE (unit_nr, *) "sqrt and inv(sqrt) of AO overlap matrix"
               ENDIF
               CALL dbcsr_create(almo_scf_env%matrix_s_sqrt(1), &
                                 template=almo_scf_env%matrix_s(1), &
                                 matrix_type=dbcsr_type_no_symmetry)
               CALL dbcsr_create(almo_scf_env%matrix_s_sqrt_inv(1), &
                                 template=almo_scf_env%matrix_s(1), &
                                 matrix_type=dbcsr_type_no_symmetry)

               CALL matrix_sqrt_Newton_Schulz(almo_scf_env%matrix_s_sqrt(1), &
                                              almo_scf_env%matrix_s_sqrt_inv(1), &
                                              almo_scf_env%matrix_s(1), &
                                              threshold=almo_scf_env%eps_filter, &
                                              order=almo_scf_env%order_lanczos, &
                                              eps_lanczos=almo_scf_env%eps_lanczos, &
                                              max_iter_lanczos=almo_scf_env%max_iter_lanczos)

               IF (safe_mode) THEN
                  CALL dbcsr_create(matrix_tmp1, template=almo_scf_env%matrix_s(1), &
                                    matrix_type=dbcsr_type_no_symmetry)
                  CALL dbcsr_create(matrix_tmp2, template=almo_scf_env%matrix_s(1), &
                                    matrix_type=dbcsr_type_no_symmetry)

                  CALL dbcsr_multiply("N", "N", 1.0_dp, almo_scf_env%matrix_s_sqrt_inv(1), &
                                      almo_scf_env%matrix_s(1), &
                                      0.0_dp, matrix_tmp1, filter_eps=almo_scf_env%eps_filter)
                  CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, almo_scf_env%matrix_s_sqrt_inv(1), &
                                      0.0_dp, matrix_tmp2, filter_eps=almo_scf_env%eps_filter)

                  frob_matrix_base = dbcsr_frobenius_norm(matrix_tmp2)
                  CALL dbcsr_add_on_diag(matrix_tmp2, -1.0_dp)
                  frob_matrix = dbcsr_frobenius_norm(matrix_tmp2)
                  IF (unit_nr > 0) THEN
                     WRITE (unit_nr, *) "Error for (inv(sqrt(S))*S*inv(sqrt(S))-I)", frob_matrix/frob_matrix_base
                  ENDIF

                  CALL dbcsr_release(matrix_tmp1)
                  CALL dbcsr_release(matrix_tmp2)
               ENDIF

               almo_scf_env%s_sqrt_done = .TRUE.

            ENDIF

            DO ispin = 1, nspin

               CALL ct_step_env_init(ct_step_env)
               CALL ct_step_env_set(ct_step_env, &
                                    para_env=almo_scf_env%para_env, &
                                    blacs_env=almo_scf_env%blacs_env, &
                                    use_occ_orbs=.TRUE., &
                                    use_virt_orbs=almo_scf_env%deloc_cayley_use_virt_orbs, &
                                    occ_orbs_orthogonal=.FALSE., &
                                    virt_orbs_orthogonal=almo_scf_env%orthogonal_basis, &
                                    tensor_type=almo_scf_env%deloc_cayley_tensor_type, &
                                    neglect_quadratic_term=almo_scf_env%deloc_cayley_linear, &
                                    calculate_energy_corr=.TRUE., &
                                    update_p=.TRUE., &
                                    update_q=.FALSE., &
                                    pp_preconditioner_full=almo_scf_env%deloc_cayley_occ_precond, &
                                    qq_preconditioner_full=almo_scf_env%deloc_cayley_vir_precond, &
                                    eps_convergence=almo_scf_env%deloc_cayley_eps_convergence, &
                                    eps_filter=almo_scf_env%eps_filter, &
                                    !nspins=almo_scf_env%nspins,&
                                    q_index_up=almo_scf_env%matrix_s_sqrt_inv(1), &
                                    q_index_down=almo_scf_env%matrix_s_sqrt(1), &
                                    p_index_up=almo_scf_env%matrix_sigma_sqrt_inv(ispin), &
                                    p_index_down=almo_scf_env%matrix_sigma_sqrt(ispin), &
                                    matrix_ks=almo_scf_env%matrix_ks_0deloc(ispin), &
                                    matrix_p=almo_scf_env%matrix_p(ispin), &
                                    matrix_qp_template=almo_scf_env%matrix_t(ispin), &
                                    matrix_pq_template=almo_scf_env%matrix_t_tr(ispin), &
                                    matrix_t=almo_scf_env%matrix_t(ispin), &
                                    conjugator=almo_scf_env%deloc_cayley_conjugator, &
                                    max_iter=almo_scf_env%deloc_cayley_max_iter)

               ! perform calculations
               CALL ct_step_execute(ct_step_env)

               ! for now we do not need the new set of orbitals
               ! just get the energy correction
               CALL ct_step_env_get(ct_step_env, &
                                    energy_correction=energy_correction(ispin))
               !copy_da_energy_matrix=matrix_eda(ispin),&
               !copy_da_charge_matrix=matrix_cta(ispin),&

               CALL ct_step_env_clean(ct_step_env)

            ENDDO

            energy_correction(1) = energy_correction(1)*spin_factor

         ENDIF

         ! print the energy correction and exit
         DO ispin = 1, nspin

            IF (unit_nr > 0) THEN
               WRITE (unit_nr, *)
               WRITE (unit_nr, '(T2,A,I6,F20.9)') "ECORR", ispin, &
                  energy_correction(ispin)
               WRITE (unit_nr, *)
            ENDIF
            energy_correction_final = energy_correction_final + energy_correction(ispin)

            !!! print out the results of decomposition analysis
            !!IF (unit_nr>0) THEN
            !!   WRITE(unit_nr,*)
            !!   WRITE(unit_nr,'(T2,A)') "ENERGY DECOMPOSITION"
            !!ENDIF
            !!CALL dbcsr_print_block_sum(eda_matrix(ispin))
            !!IF (unit_nr>0) THEN
            !!   WRITE(unit_nr,*)
            !!   WRITE(unit_nr,'(T2,A)') "CHARGE DECOMPOSITION"
            !!ENDIF
            !!CALL dbcsr_print_block_sum(cta_matrix(ispin))

            ! obtain density matrix from updated MOs
            ! RZK-later sigma and sigma_inv are lost here
            CALL almo_scf_t_to_proj(t=almo_scf_env%matrix_t(ispin), &
                                    p=almo_scf_env%matrix_p(ispin), &
                                    eps_filter=almo_scf_env%eps_filter, &
                                    orthog_orbs=.FALSE., &
                                    nocc_of_domain=almo_scf_env%nocc_of_domain(:, ispin), &
                                    s=almo_scf_env%matrix_s(1), &
                                    sigma=almo_scf_env%matrix_sigma(ispin), &
                                    sigma_inv=almo_scf_env%matrix_sigma_inv(ispin), &
                                    !use_guess=use_guess, &
                                    algorithm=almo_scf_env%sigma_inv_algorithm, &
                                    inverse_accelerator=almo_scf_env%order_lanczos, &
                                    inv_eps_factor=almo_scf_env%matrix_iter_eps_error_factor, &
                                    eps_lanczos=almo_scf_env%eps_lanczos, &
                                    max_iter_lanczos=almo_scf_env%max_iter_lanczos, &
                                    para_env=almo_scf_env%para_env, &
                                    blacs_env=almo_scf_env%blacs_env)

            IF (almo_scf_env%nspins == 1) &
               CALL dbcsr_scale(almo_scf_env%matrix_p(ispin), &
                                spin_factor)

         ENDDO

      CASE (dm_ls_step)

         ! compute the inverse of S
         IF (.NOT. almo_scf_env%s_inv_done) THEN
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, *) "Inverting AO overlap matrix"
            ENDIF
            CALL dbcsr_create(almo_scf_env%matrix_s_inv(1), &
                              template=almo_scf_env%matrix_s(1), &
                              matrix_type=dbcsr_type_no_symmetry)
            IF (.NOT. almo_scf_env%s_sqrt_done) THEN
               CALL invert_Hotelling(almo_scf_env%matrix_s_inv(1), &
                                     almo_scf_env%matrix_s(1), &
                                     threshold=almo_scf_env%eps_filter)
            ELSE
               CALL dbcsr_multiply("N", "N", 1.0_dp, almo_scf_env%matrix_s_sqrt_inv(1), &
                                   almo_scf_env%matrix_s_sqrt_inv(1), &
                                   0.0_dp, almo_scf_env%matrix_s_inv(1), &
                                   filter_eps=almo_scf_env%eps_filter)
            ENDIF

            IF (safe_mode) THEN
               CALL dbcsr_create(matrix_tmp1, template=almo_scf_env%matrix_s(1), &
                                 matrix_type=dbcsr_type_no_symmetry)
               CALL dbcsr_multiply("N", "N", 1.0_dp, almo_scf_env%matrix_s_inv(1), &
                                   almo_scf_env%matrix_s(1), &
                                   0.0_dp, matrix_tmp1, &
                                   filter_eps=almo_scf_env%eps_filter)
               frob_matrix_base = dbcsr_frobenius_norm(matrix_tmp1)
               CALL dbcsr_add_on_diag(matrix_tmp1, -1.0_dp)
               frob_matrix = dbcsr_frobenius_norm(matrix_tmp1)
               IF (unit_nr > 0) THEN
                  WRITE (unit_nr, *) "Error for (inv(S)*S-I)", &
                     frob_matrix/frob_matrix_base
               ENDIF
               CALL dbcsr_release(matrix_tmp1)
            ENDIF

            almo_scf_env%s_inv_done = .TRUE.

         ENDIF

         DO ispin = 1, nspin
            ! RZK-warning the preconditioner is very important
            !       IF (.FALSE.) THEN
            !           CALL apply_matrix_preconditioner(almo_scf_env%matrix_ks(ispin),&
            !                   "forward",almo_scf_env%matrix_s_blk_sqrt(1),&
            !                   almo_scf_env%matrix_s_blk_sqrt_inv(1))
            !       ENDIF
            !CALL dbcsr_filter(almo_scf_env%matrix_ks(ispin),&
            !         almo_scf_env%eps_filter)
         ENDDO

         ALLOCATE (matrix_p_almo_scf_converged(nspin))
         DO ispin = 1, nspin
            CALL dbcsr_create(matrix_p_almo_scf_converged(ispin), &
                              template=almo_scf_env%matrix_p(ispin))
            CALL dbcsr_copy(matrix_p_almo_scf_converged(ispin), &
                            almo_scf_env%matrix_p(ispin))
         ENDDO

         ! update the density matrix
         DO ispin = 1, nspin

            nelectron_spin_real(1) = almo_scf_env%nelectrons_spin(ispin)
            IF (almo_scf_env%nspins == 1) &
               nelectron_spin_real(1) = nelectron_spin_real(1)/2

            local_mu(1) = SUM(almo_scf_env%mu_of_domain(:, ispin))/almo_scf_env%ndomains
            fake(1) = 123523

            ! RZK UPDATE! the update algorithm is removed because
            ! RZK UPDATE! it requires updating core LS_SCF routines
            ! RZK UPDATE! (the code exists in the CVS version)
            CPABORT("CVS only: density_matrix_sign has not been updated in SVN")
            ! RZK UPDATE!CALL density_matrix_sign(almo_scf_env%matrix_p(ispin),&
            ! RZK UPDATE!                     local_mu,&
            ! RZK UPDATE!                     almo_scf_env%fixed_mu,&
            ! RZK UPDATE!                     almo_scf_env%matrix_ks_0deloc(ispin),&
            ! RZK UPDATE!                     almo_scf_env%matrix_s(1), &
            ! RZK UPDATE!                     almo_scf_env%matrix_s_inv(1), &
            ! RZK UPDATE!                     nelectron_spin_real,&
            ! RZK UPDATE!                     almo_scf_env%eps_filter,&
            ! RZK UPDATE!                     fake)
            ! RZK UPDATE!
            almo_scf_env%mu = local_mu(1)

            !IF (almo_scf_env%has_s_preconditioner) THEN
            !    CALL apply_matrix_preconditioner(&
            !             almo_scf_env%matrix_p_blk(ispin),&
            !             "forward",almo_scf_env%matrix_s_blk_sqrt(1),&
            !             almo_scf_env%matrix_s_blk_sqrt_inv(1))
            !ENDIF
            !CALL dbcsr_filter(almo_scf_env%matrix_p(ispin),&
            !        almo_scf_env%eps_filter)

            IF (almo_scf_env%nspins == 1) &
               CALL dbcsr_scale(almo_scf_env%matrix_p(ispin), &
                                spin_factor)

            !CALL dbcsr_dot(almo_scf_env%matrix_ks_0deloc(ispin),&
            !  almo_scf_env%matrix_p(ispin),&
            !  energy_correction(ispin))
            !IF (unit_nr>0) THEN
            !   WRITE(unit_nr,*)
            !   WRITE(unit_nr,'(T2,A,I6,F20.9)') "EFAKE",ispin,&
            !           energy_correction(ispin)
            !   WRITE(unit_nr,*)
            !ENDIF
            CALL dbcsr_add(matrix_p_almo_scf_converged(ispin), &
                           almo_scf_env%matrix_p(ispin), -1.0_dp, 1.0_dp)
            CALL dbcsr_dot(almo_scf_env%matrix_ks_0deloc(ispin), &
                           matrix_p_almo_scf_converged(ispin), &
                           energy_correction(ispin))

            energy_correction_final = energy_correction_final + energy_correction(ispin)

            IF (unit_nr > 0) THEN
               WRITE (unit_nr, *)
               WRITE (unit_nr, '(T2,A,I6,F20.9)') "ECORR", ispin, &
                  energy_correction(ispin)
               WRITE (unit_nr, *)
            ENDIF

         ENDDO

         DO ispin = 1, nspin
            CALL dbcsr_release(matrix_p_almo_scf_converged(ispin))
         ENDDO
         DEALLOCATE (matrix_p_almo_scf_converged)

      END SELECT ! algorithm selection

      t2 = m_walltime()

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *)
         WRITE (unit_nr, '(T2,A,F18.9,F18.9,F18.9,F12.6)') "ETOT", &
            almo_scf_env%almo_scf_energy, &
            energy_correction_final, &
            almo_scf_env%almo_scf_energy + energy_correction_final, &
            t2 - t1
         WRITE (unit_nr, *)
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE harris_foulkes_correction

! **************************************************************************************************
!> \brief Computes a diagonal preconditioner for the cg optimization of k matrix
!> \param prec ...
!> \param vd_prop ...
!> \param f ...
!> \param x ...
!> \param oo_inv_x_tr ...
!> \param s ...
!> \param grad ...
!> \param vd_blk ...
!> \param t ...
!> \param template_vd_vd_blk ...
!> \param template_vr_vr_blk ...
!> \param template_n_vr ...
!> \param spin_factor ...
!> \param eps_filter ...
!> \par History
!>       2011.09 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE opt_k_create_preconditioner(prec, vd_prop, f, x, oo_inv_x_tr, s, grad, &
                                          vd_blk, t, template_vd_vd_blk, template_vr_vr_blk, template_n_vr, &
                                          spin_factor, eps_filter)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: prec
      TYPE(dbcsr_type), INTENT(IN)                       :: vd_prop, f, x, oo_inv_x_tr, s, grad, &
                                                            vd_blk, t, template_vd_vd_blk, &
                                                            template_vr_vr_blk, template_n_vr
      REAL(KIND=dp), INTENT(IN)                          :: spin_factor, eps_filter

      CHARACTER(len=*), PARAMETER :: routineN = 'opt_k_create_preconditioner', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, p_nrows, q_nrows
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: p_diagonal, q_diagonal
      TYPE(dbcsr_type)                                   :: pp_diag, qq_diag, t1, t2, tmp, &
                                                            tmp1_n_vr, tmp2_n_vr, tmp_n_vd, &
                                                            tmp_vd_vd_blk, tmp_vr_vr_blk

! init diag blocks outside
! init diag blocks otside
!INTEGER                                  :: iblock_row, iblock_col,&
!                                            nblkrows_tot, nblkcols_tot
!REAL(KIND=dp), DIMENSION(:, :), POINTER  :: p_new_block
!INTEGER                                  :: mynode, hold, row, col

      CALL timeset(routineN, handle)

      ! initialize a matrix to 1.0
      CALL dbcsr_create(tmp, template=prec)
      ! in order to use dbcsr_set matrix blocks must exist
      CALL dbcsr_copy(tmp, prec)
      CALL dbcsr_set(tmp, 1.0_dp)

      ! compute qq = (Vd^tr)*F*Vd
      CALL dbcsr_create(tmp_n_vd, template=vd_prop)
      CALL dbcsr_multiply("N", "N", 1.0_dp, f, vd_prop, &
                          0.0_dp, tmp_n_vd, filter_eps=eps_filter)
      CALL dbcsr_create(tmp_vd_vd_blk, &
                        template=template_vd_vd_blk)
      CALL dbcsr_copy(tmp_vd_vd_blk, template_vd_vd_blk)
      CALL dbcsr_multiply("T", "N", 1.0_dp, vd_prop, tmp_n_vd, &
                          0.0_dp, tmp_vd_vd_blk, &
                          retain_sparsity=.TRUE., &
                          filter_eps=eps_filter)
      ! copy diagonal elements of the result into rows of a matrix
      CALL dbcsr_get_info(tmp_vd_vd_blk, nfullrows_total=q_nrows)
      ALLOCATE (q_diagonal(q_nrows))
      CALL dbcsr_get_diag(tmp_vd_vd_blk, q_diagonal)
      CALL dbcsr_create(qq_diag, &
                        template=template_vd_vd_blk)
      CALL dbcsr_add_on_diag(qq_diag, 1.0_dp)
      CALL dbcsr_set_diag(qq_diag, q_diagonal)
      CALL dbcsr_create(t1, template=prec)
      CALL dbcsr_multiply("N", "N", 1.0_dp, qq_diag, tmp, &
                          0.0_dp, t1, filter_eps=eps_filter)

      ! compute pp = X*sigma_oo_inv*X^tr
      CALL dbcsr_create(tmp_vr_vr_blk, template=template_vr_vr_blk)
      CALL dbcsr_copy(tmp_vr_vr_blk, template_vr_vr_blk)
      CALL dbcsr_multiply("N", "N", 1.0_dp, x, oo_inv_x_tr, &
                          0.0_dp, tmp_vr_vr_blk, &
                          retain_sparsity=.TRUE., &
                          filter_eps=eps_filter)
      ! copy diagonal elements of the result into cols of a matrix
      CALL dbcsr_get_info(tmp_vr_vr_blk, nfullrows_total=p_nrows)
      ALLOCATE (p_diagonal(p_nrows))
      CALL dbcsr_get_diag(tmp_vr_vr_blk, p_diagonal)
      CALL dbcsr_create(pp_diag, template=template_vr_vr_blk)
      CALL dbcsr_add_on_diag(pp_diag, 1.0_dp)
      CALL dbcsr_set_diag(pp_diag, p_diagonal)
      CALL dbcsr_set(tmp, 1.0_dp)
      CALL dbcsr_create(t2, template=prec)
      CALL dbcsr_multiply("N", "N", 1.0_dp, tmp, pp_diag, &
                          0.0_dp, t2, filter_eps=eps_filter)

      CALL dbcsr_hadamard_product(t1, t2, prec)

      ! compute qq = (Vd^tr)*S*Vd
      CALL dbcsr_multiply("N", "N", 1.0_dp, s, vd_prop, &
                          0.0_dp, tmp_n_vd, filter_eps=eps_filter)
      CALL dbcsr_multiply("T", "N", 1.0_dp, vd_prop, tmp_n_vd, &
                          0.0_dp, tmp_vd_vd_blk, &
                          retain_sparsity=.TRUE., &
                          filter_eps=eps_filter)
      ! copy diagonal elements of the result into rows of a matrix
      CALL dbcsr_get_diag(tmp_vd_vd_blk, q_diagonal)
      CALL dbcsr_add_on_diag(qq_diag, 1.0_dp)
      CALL dbcsr_set_diag(qq_diag, q_diagonal)
      CALL dbcsr_set(tmp, 1.0_dp)
      CALL dbcsr_multiply("N", "N", 1.0_dp, qq_diag, tmp, &
                          0.0_dp, t1, filter_eps=eps_filter)

      ! compute pp = X*sig_oo_inv*(T^tr)*F*T*sig_oo_inv*(X^tr)
      CALL dbcsr_create(tmp1_n_vr, template=template_n_vr)
      CALL dbcsr_create(tmp2_n_vr, template=template_n_vr)
      CALL dbcsr_multiply("N", "N", 1.0_dp, t, oo_inv_x_tr, &
                          0.0_dp, tmp1_n_vr, filter_eps=eps_filter)
      CALL dbcsr_multiply("N", "N", 1.0_dp, f, tmp1_n_vr, &
                          0.0_dp, tmp2_n_vr, filter_eps=eps_filter)
      CALL dbcsr_multiply("T", "N", 1.0_dp, tmp1_n_vr, tmp2_n_vr, &
                          0.0_dp, tmp_vr_vr_blk, &
                          retain_sparsity=.TRUE., &
                          filter_eps=eps_filter)
      ! copy diagonal elements of the result into cols of a matrix
      CALL dbcsr_get_diag(tmp_vr_vr_blk, p_diagonal)
      CALL dbcsr_add_on_diag(pp_diag, 1.0_dp)
      CALL dbcsr_set_diag(pp_diag, p_diagonal)
      CALL dbcsr_set(tmp, 1.0_dp)
      CALL dbcsr_multiply("N", "N", 1.0_dp, tmp, pp_diag, &
                          0.0_dp, t2, filter_eps=eps_filter)

      CALL dbcsr_hadamard_product(t1, t2, tmp)
      CALL dbcsr_add(prec, tmp, 1.0_dp, -1.0_dp)
      CALL dbcsr_scale(prec, 2.0_dp*spin_factor)

      ! compute qp = X*sig_oo_inv*(T^tr)*S*Vd
      CALL dbcsr_multiply("N", "N", 1.0_dp, s, vd_blk, &
                          0.0_dp, tmp_n_vd, filter_eps=eps_filter)
      CALL dbcsr_multiply("T", "N", 1.0_dp, tmp_n_vd, tmp1_n_vr, &
                          0.0_dp, tmp, retain_sparsity=.TRUE., &
                          filter_eps=eps_filter)
      CALL dbcsr_hadamard_product(grad, tmp, t1)
      ! gradient already contains 2.0*spin_factor
      CALL dbcsr_scale(t1, -2.0_dp)

      CALL dbcsr_add(prec, t1, 1.0_dp, 1.0_dp)

      CALL dbcsr_function_of_elements(prec, dbcsr_func_inverse)
      CALL dbcsr_filter(prec, eps_filter)

      DEALLOCATE (q_diagonal)
      DEALLOCATE (p_diagonal)
      CALL dbcsr_release(tmp)
      CALL dbcsr_release(qq_diag)
      CALL dbcsr_release(t1)
      CALL dbcsr_release(pp_diag)
      CALL dbcsr_release(t2)
      CALL dbcsr_release(tmp_n_vd)
      CALL dbcsr_release(tmp_vd_vd_blk)
      CALL dbcsr_release(tmp_vr_vr_blk)
      CALL dbcsr_release(tmp1_n_vr)
      CALL dbcsr_release(tmp2_n_vr)

      CALL timestop(handle)

   END SUBROUTINE opt_k_create_preconditioner

! **************************************************************************************************
!> \brief Computes a block-diagonal preconditioner for the optimization of
!>        k matrix
!> \param almo_scf_env ...
!> \param vd_prop ...
!> \param oo_inv_x_tr ...
!> \param t_curr ...
!> \param ispin ...
!> \param spin_factor ...
!> \par History
!>       2011.10 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE opt_k_create_preconditioner_blk(almo_scf_env, vd_prop, oo_inv_x_tr, &
                                              t_curr, ispin, spin_factor)

      TYPE(almo_scf_env_type), INTENT(INOUT)             :: almo_scf_env
      TYPE(dbcsr_type), INTENT(IN)                       :: vd_prop, oo_inv_x_tr, t_curr
      INTEGER, INTENT(IN)                                :: ispin
      REAL(KIND=dp), INTENT(IN)                          :: spin_factor

      CHARACTER(len=*), PARAMETER :: routineN = 'opt_k_create_preconditioner_blk', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: eps_filter
      TYPE(dbcsr_type)                                   :: opt_k_e_dd, opt_k_e_rr, s_dd_sqrt, &
                                                            s_rr_sqrt, t1, tmp, tmp1_n_vr, &
                                                            tmp2_n_vr, tmp_n_vd, tmp_vd_vd_blk, &
                                                            tmp_vr_vr_blk

! matrices that has been computed outside the routine already

      CALL timeset(routineN, handle)

      eps_filter = almo_scf_env%eps_filter

      ! compute S_qq = (Vd^tr)*S*Vd
      CALL dbcsr_create(tmp_n_vd, template=almo_scf_env%matrix_v_disc(ispin))
      CALL dbcsr_create(tmp_vd_vd_blk, &
                        template=almo_scf_env%matrix_vv_disc_blk(ispin), &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_multiply("N", "N", 1.0_dp, &
                          almo_scf_env%matrix_s(1), &
                          vd_prop, &
                          0.0_dp, tmp_n_vd, filter_eps=eps_filter)
      CALL dbcsr_copy(tmp_vd_vd_blk, &
                      almo_scf_env%matrix_vv_disc_blk(ispin))
      CALL dbcsr_multiply("T", "N", 1.0_dp, vd_prop, tmp_n_vd, &
                          0.0_dp, tmp_vd_vd_blk, &
                          retain_sparsity=.TRUE.)

      CALL dbcsr_create(s_dd_sqrt, &
                        template=almo_scf_env%matrix_vv_disc_blk(ispin), &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL matrix_sqrt_Newton_Schulz(s_dd_sqrt, &
                                     almo_scf_env%opt_k_t_dd(ispin), &
                                     tmp_vd_vd_blk, &
                                     threshold=eps_filter, &
                                     order=almo_scf_env%order_lanczos, &
                                     eps_lanczos=almo_scf_env%eps_lanczos, &
                                     max_iter_lanczos=almo_scf_env%max_iter_lanczos)

      ! compute F_qq = (Vd^tr)*F*Vd
      CALL dbcsr_multiply("N", "N", 1.0_dp, &
                          almo_scf_env%matrix_ks_0deloc(ispin), &
                          vd_prop, &
                          0.0_dp, tmp_n_vd, filter_eps=eps_filter)
      CALL dbcsr_copy(tmp_vd_vd_blk, &
                      almo_scf_env%matrix_vv_disc_blk(ispin))
      CALL dbcsr_multiply("T", "N", 1.0_dp, vd_prop, tmp_n_vd, &
                          0.0_dp, tmp_vd_vd_blk, &
                          retain_sparsity=.TRUE.)
      CALL dbcsr_release(tmp_n_vd)

      ! bring to the blocked-orthogonalized basis
      CALL dbcsr_multiply("N", "N", 1.0_dp, &
                          tmp_vd_vd_blk, &
                          almo_scf_env%opt_k_t_dd(ispin), &
                          0.0_dp, s_dd_sqrt, filter_eps=eps_filter)
      CALL dbcsr_multiply("N", "N", 1.0_dp, &
                          almo_scf_env%opt_k_t_dd(ispin), &
                          s_dd_sqrt, &
                          0.0_dp, tmp_vd_vd_blk, filter_eps=eps_filter)

      ! diagonalize the matrix
      CALL dbcsr_create(opt_k_e_dd, &
                        template=almo_scf_env%matrix_vv_disc_blk(ispin))
      CALL dbcsr_release(s_dd_sqrt)
      CALL dbcsr_create(s_dd_sqrt, &
                        template=almo_scf_env%matrix_vv_disc_blk(ispin), &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL diagonalize_diagonal_blocks(tmp_vd_vd_blk, &
                                       s_dd_sqrt, &
                                       opt_k_e_dd)

      ! obtain the transformation matrix in the discarded subspace
      ! T = S^{-1/2}.U
      CALL dbcsr_copy(tmp_vd_vd_blk, &
                      almo_scf_env%opt_k_t_dd(ispin))
      CALL dbcsr_multiply("N", "N", 1.0_dp, &
                          tmp_vd_vd_blk, &
                          s_dd_sqrt, &
                          0.0_dp, almo_scf_env%opt_k_t_dd(ispin), &
                          filter_eps=eps_filter)
      CALL dbcsr_release(s_dd_sqrt)
      CALL dbcsr_release(tmp_vd_vd_blk)

      ! copy diagonal elements of the result into rows of a matrix
      CALL dbcsr_create(tmp, &
                        template=almo_scf_env%matrix_k_blk_ones(ispin))
      CALL dbcsr_copy(tmp, &
                      almo_scf_env%matrix_k_blk_ones(ispin))
      CALL dbcsr_create(t1, &
                        template=almo_scf_env%matrix_k_blk_ones(ispin))
      CALL dbcsr_multiply("N", "N", 1.0_dp, &
                          opt_k_e_dd, tmp, &
                          0.0_dp, t1, filter_eps=eps_filter)
      CALL dbcsr_release(opt_k_e_dd)

      ! compute S_pp = X*sigma_oo_inv*X^tr
      CALL dbcsr_create(tmp_vr_vr_blk, &
                        template=almo_scf_env%matrix_sigma_vv_blk(ispin), &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_copy(tmp_vr_vr_blk, &
                      almo_scf_env%matrix_sigma_vv_blk(ispin))
      CALL dbcsr_multiply("N", "N", 1.0_dp, &
                          almo_scf_env%matrix_x(ispin), &
                          oo_inv_x_tr, &
                          0.0_dp, tmp_vr_vr_blk, &
                          retain_sparsity=.TRUE.)

      ! obtain the orthogonalization matrix
      CALL dbcsr_create(s_rr_sqrt, &
                        template=almo_scf_env%matrix_sigma_vv_blk(ispin), &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL matrix_sqrt_Newton_Schulz(s_rr_sqrt, &
                                     almo_scf_env%opt_k_t_rr(ispin), &
                                     tmp_vr_vr_blk, &
                                     threshold=eps_filter, &
                                     order=almo_scf_env%order_lanczos, &
                                     eps_lanczos=almo_scf_env%eps_lanczos, &
                                     max_iter_lanczos=almo_scf_env%max_iter_lanczos)

      ! compute F_pp = X*sig_oo_inv*(T^tr)*F*T*sig_oo_inv*(X^tr)
      CALL dbcsr_create(tmp1_n_vr, &
                        template=almo_scf_env%matrix_v(ispin))
      CALL dbcsr_create(tmp2_n_vr, &
                        template=almo_scf_env%matrix_v(ispin))
      CALL dbcsr_multiply("N", "N", 1.0_dp, t_curr, oo_inv_x_tr, &
                          0.0_dp, tmp1_n_vr, filter_eps=eps_filter)
      CALL dbcsr_multiply("N", "N", 1.0_dp, &
                          almo_scf_env%matrix_ks_0deloc(ispin), &
                          tmp1_n_vr, &
                          0.0_dp, tmp2_n_vr, filter_eps=eps_filter)
      CALL dbcsr_multiply("T", "N", 1.0_dp, tmp1_n_vr, tmp2_n_vr, &
                          0.0_dp, tmp_vr_vr_blk, &
                          retain_sparsity=.TRUE.)
      CALL dbcsr_release(tmp1_n_vr)
      CALL dbcsr_release(tmp2_n_vr)

      ! bring to the blocked-orthogonalized basis
      CALL dbcsr_multiply("N", "N", 1.0_dp, &
                          tmp_vr_vr_blk, &
                          almo_scf_env%opt_k_t_rr(ispin), &
                          0.0_dp, s_rr_sqrt, filter_eps=eps_filter)
      CALL dbcsr_multiply("N", "N", 1.0_dp, &
                          almo_scf_env%opt_k_t_rr(ispin), &
                          s_rr_sqrt, &
                          0.0_dp, tmp_vr_vr_blk, filter_eps=eps_filter)

      ! diagonalize the matrix
      CALL dbcsr_create(opt_k_e_rr, &
                        template=almo_scf_env%matrix_sigma_vv_blk(ispin))
      CALL dbcsr_release(s_rr_sqrt)
      CALL dbcsr_create(s_rr_sqrt, &
                        template=almo_scf_env%matrix_sigma_vv_blk(ispin), &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL diagonalize_diagonal_blocks(tmp_vr_vr_blk, &
                                       s_rr_sqrt, &
                                       opt_k_e_rr)

      ! obtain the transformation matrix in the retained subspace
      ! T = S^{-1/2}.U
      CALL dbcsr_copy(tmp_vr_vr_blk, &
                      almo_scf_env%opt_k_t_rr(ispin))
      CALL dbcsr_multiply("N", "N", 1.0_dp, &
                          tmp_vr_vr_blk, &
                          s_rr_sqrt, &
                          0.0_dp, almo_scf_env%opt_k_t_rr(ispin), &
                          filter_eps=eps_filter)
      CALL dbcsr_release(s_rr_sqrt)
      CALL dbcsr_release(tmp_vr_vr_blk)

      ! copy diagonal elements of the result into cols of a matrix
      CALL dbcsr_multiply("N", "N", 1.0_dp, &
                          tmp, opt_k_e_rr, &
                          0.0_dp, almo_scf_env%opt_k_denom(ispin), &
                          filter_eps=eps_filter)
      CALL dbcsr_release(opt_k_e_rr)
      CALL dbcsr_release(tmp)

      ! form the denominator matrix
      CALL dbcsr_add(almo_scf_env%opt_k_denom(ispin), t1, &
                     -1.0_dp, 1.0_dp)
      CALL dbcsr_release(t1)
      CALL dbcsr_scale(almo_scf_env%opt_k_denom(ispin), &
                       2.0_dp*spin_factor)

      CALL dbcsr_function_of_elements(almo_scf_env%opt_k_denom(ispin), &
                                      dbcsr_func_inverse)
      CALL dbcsr_filter(almo_scf_env%opt_k_denom(ispin), &
                        eps_filter)

      CALL timestop(handle)

   END SUBROUTINE opt_k_create_preconditioner_blk

! **************************************************************************************************
!> \brief Applies a block-diagonal preconditioner for the optimization of
!>        k matrix (preconditioner matrices must be calculated and stored
!>        beforehand)
!> \param almo_scf_env ...
!> \param step ...
!> \param grad ...
!> \param ispin ...
!> \par History
!>       2011.10 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE opt_k_apply_preconditioner_blk(almo_scf_env, step, grad, ispin)

      TYPE(almo_scf_env_type), INTENT(INOUT)             :: almo_scf_env
      TYPE(dbcsr_type), INTENT(OUT)                      :: step
      TYPE(dbcsr_type), INTENT(IN)                       :: grad
      INTEGER, INTENT(IN)                                :: ispin

      CHARACTER(len=*), PARAMETER :: routineN = 'opt_k_apply_preconditioner_blk', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: eps_filter
      TYPE(dbcsr_type)                                   :: tmp_k

      CALL timeset(routineN, handle)

      eps_filter = almo_scf_env%eps_filter

      CALL dbcsr_create(tmp_k, template=almo_scf_env%matrix_k_blk(ispin))

      ! transform gradient to the correct "diagonal" basis
      CALL dbcsr_multiply("N", "N", 1.0_dp, &
                          grad, almo_scf_env%opt_k_t_rr(ispin), &
                          0.0_dp, tmp_k, filter_eps=eps_filter)
      CALL dbcsr_multiply("T", "N", 1.0_dp, &
                          almo_scf_env%opt_k_t_dd(ispin), tmp_k, &
                          0.0_dp, step, filter_eps=eps_filter)

      ! apply diagonal preconditioner
      CALL dbcsr_hadamard_product(step, &
                                  almo_scf_env%opt_k_denom(ispin), tmp_k)

      ! back-transform the result to the initial basis
      CALL dbcsr_multiply("N", "N", 1.0_dp, &
                          almo_scf_env%opt_k_t_dd(ispin), tmp_k, &
                          0.0_dp, step, filter_eps=eps_filter)
      CALL dbcsr_multiply("N", "T", 1.0_dp, &
                          step, almo_scf_env%opt_k_t_rr(ispin), &
                          0.0_dp, tmp_k, filter_eps=eps_filter)

      CALL dbcsr_copy(step, tmp_k)

      CALL dbcsr_release(tmp_k)

      CALL timestop(handle)

   END SUBROUTINE opt_k_apply_preconditioner_blk

!! **************************************************************************************************
!!> \brief Reduce the number of virtual orbitals by rotating them within
!!>        a domain. The rotation is such that minimizes the frobenius norm of
!!>        the Fov domain-blocks of the discarded virtuals
!!> \par History
!!>       2011.08 created [Rustam Z Khaliullin]
!!> \author Rustam Z Khaliullin
!! **************************************************************************************************
!  SUBROUTINE truncate_subspace_v_blk(qs_env,almo_scf_env)
!
!    TYPE(qs_environment_type), POINTER       :: qs_env
!    TYPE(almo_scf_env_type)                  :: almo_scf_env
!
!    CHARACTER(len=*), PARAMETER :: routineN = 'truncate_subspace_v_blk', &
!      routineP = moduleN//':'//routineN
!
!    INTEGER                                  :: handle, ispin, iblock_row, &
!                                                iblock_col, iblock_row_size, &
!                                                iblock_col_size, retained_v, &
!                                                iteration, line_search_step, &
!                                                unit_nr, line_search_step_last
!    REAL(KIND=dp)                            :: t1, obj_function, grad_norm,&
!                                                c0, b0, a0, obj_function_new,&
!                                                t2, alpha, ff1, ff2, step1,&
!                                                step2,&
!                                                frob_matrix_base,&
!                                                frob_matrix
!    LOGICAL                                  :: safe_mode, converged, &
!                                                prepare_to_exit, failure
!    TYPE(cp_logger_type), POINTER            :: logger
!    TYPE(dbcsr_type)                      :: Fon, Fov, Fov_filtered, &
!                                                temp1_oo, temp2_oo, Fov_original, &
!                                                temp0_ov, U_blk_tot, U_blk, &
!                                                grad_blk, step_blk, matrix_filter, &
!                                                v_full_new,v_full_tmp,&
!                                                matrix_sigma_vv_full,&
!                                                matrix_sigma_vv_full_sqrt,&
!                                                matrix_sigma_vv_full_sqrt_inv,&
!                                                matrix_tmp1,&
!                                                matrix_tmp2
!
!    REAL(kind=dp), DIMENSION(:, :), POINTER  :: data_p, p_new_block
!    TYPE(dbcsr_iterator_type)                  :: iter
!
!
!REAL(kind=dp), DIMENSION(:), ALLOCATABLE     :: eigenvalues, WORK
!REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE   :: data_copy, left_vectors, right_vectors
!INTEGER                                      :: LWORK, INFO
!TYPE(dbcsr_type)                          :: temp_u_v_full_blk
!
!    CALL timeset(routineN,handle)
!
!    safe_mode=.TRUE.
!
!    ! get a useful output_unit
!    logger => cp_get_default_logger()
!    IF (logger%para_env%mepos==logger%para_env%source) THEN
!       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
!    ELSE
!       unit_nr=-1
!    ENDIF
!
!    DO ispin=1,almo_scf_env%nspins
!
!       t1 = m_walltime()
!
!       !!!!!!!!!!!!!!!!!
!       ! 0. Orthogonalize virtuals
!       !    Unfortunately, we have to do it in the FULL V subspace :(
!
!       CALL dbcsr_init(v_full_new)
!       CALL dbcsr_create(v_full_new,&
!               template=almo_scf_env%matrix_v_full_blk(ispin),&
!               matrix_type=dbcsr_type_no_symmetry)
!
!       ! project the occupied subspace out
!       CALL almo_scf_p_out_from_v(almo_scf_env%matrix_v_full_blk(ispin),&
!              v_full_new,almo_scf_env%matrix_ov_full(ispin),&
!              ispin,almo_scf_env)
!
!       ! init overlap and its functions
!       CALL dbcsr_init(matrix_sigma_vv_full)
!       CALL dbcsr_init(matrix_sigma_vv_full_sqrt)
!       CALL dbcsr_init(matrix_sigma_vv_full_sqrt_inv)
!       CALL dbcsr_create(matrix_sigma_vv_full,&
!               template=almo_scf_env%matrix_vv_full_blk(ispin),&
!               matrix_type=dbcsr_type_no_symmetry)
!       CALL dbcsr_create(matrix_sigma_vv_full_sqrt,&
!               template=almo_scf_env%matrix_vv_full_blk(ispin),&
!               matrix_type=dbcsr_type_no_symmetry)
!       CALL dbcsr_create(matrix_sigma_vv_full_sqrt_inv,&
!               template=almo_scf_env%matrix_vv_full_blk(ispin),&
!               matrix_type=dbcsr_type_no_symmetry)
!
!       ! construct VV overlap
!       CALL almo_scf_mo_to_sigma(v_full_new,&
!               matrix_sigma_vv_full,&
!               almo_scf_env%matrix_s(1),&
!               almo_scf_env%eps_filter)
!
!       IF (unit_nr>0) THEN
!          WRITE(unit_nr,*) "sqrt and inv(sqrt) of the FULL virtual MO overlap"
!       ENDIF
!
!       ! construct orthogonalization matrices
!       CALL matrix_sqrt_Newton_Schulz(matrix_sigma_vv_full_sqrt,&
!                                      matrix_sigma_vv_full_sqrt_inv,&
!                                      matrix_sigma_vv_full,&
!                                      threshold=almo_scf_env%eps_filter,&
!                                      order=almo_scf_env%order_lanczos,&
!                                      eps_lanczos=almo_scf_env%eps_lanczos,&
!                                      max_iter_lanczos=almo_scf_env%max_iter_lanczos)
!       IF (safe_mode) THEN
!          CALL dbcsr_init(matrix_tmp1)
!          CALL dbcsr_create(matrix_tmp1,template=matrix_sigma_vv_full,&
!                               matrix_type=dbcsr_type_no_symmetry)
!          CALL dbcsr_init(matrix_tmp2)
!          CALL dbcsr_create(matrix_tmp2,template=matrix_sigma_vv_full,&
!                               matrix_type=dbcsr_type_no_symmetry)
!
!          CALL dbcsr_multiply("N","N",1.0_dp,matrix_sigma_vv_full_sqrt_inv,&
!                                 matrix_sigma_vv_full,&
!                                 0.0_dp,matrix_tmp1,filter_eps=almo_scf_env%eps_filter)
!          CALL dbcsr_multiply("N","N",1.0_dp,matrix_tmp1,&
!                                 matrix_sigma_vv_full_sqrt_inv,&
!                                 0.0_dp,matrix_tmp2,filter_eps=almo_scf_env%eps_filter)
!
!          frob_matrix_base=dbcsr_frobenius_norm(matrix_tmp2)
!          CALL dbcsr_add_on_diag(matrix_tmp2,-1.0_dp)
!          frob_matrix=dbcsr_frobenius_norm(matrix_tmp2)
!          IF (unit_nr>0) THEN
!             WRITE(unit_nr,*) "Error for (inv(sqrt(SIGVV))*SIGVV*inv(sqrt(SIGVV))-I)",frob_matrix/frob_matrix_base
!          ENDIF
!
!          CALL dbcsr_release(matrix_tmp1)
!          CALL dbcsr_release(matrix_tmp2)
!       ENDIF
!
!       ! discard unnecessary overlap functions
!       CALL dbcsr_release(matrix_sigma_vv_full)
!       CALL dbcsr_release(matrix_sigma_vv_full_sqrt)
!
!! this can be re-written because we have (1-P)|v>
!
!       !!!!!!!!!!!!!!!!!!!
!       ! 1. Compute F_ov
!       CALL dbcsr_init(Fon)
!       CALL dbcsr_create(Fon,&
!               template=almo_scf_env%matrix_v_full_blk(ispin))
!       CALL dbcsr_init(Fov)
!       CALL dbcsr_create(Fov,&
!               template=almo_scf_env%matrix_ov_full(ispin))
!       CALL dbcsr_init(Fov_filtered)
!       CALL dbcsr_create(Fov_filtered,&
!               template=almo_scf_env%matrix_ov_full(ispin))
!       CALL dbcsr_init(temp1_oo)
!       CALL dbcsr_create(temp1_oo,&
!               template=almo_scf_env%matrix_sigma(ispin),&
!               !matrix_type=dbcsr_type_no_symmetry)
!       CALL dbcsr_init(temp2_oo)
!       CALL dbcsr_create(temp2_oo,&
!               template=almo_scf_env%matrix_sigma(ispin),&
!               matrix_type=dbcsr_type_no_symmetry)
!
!       CALL dbcsr_multiply("T","N",1.0_dp,almo_scf_env%matrix_t_blk(ispin),&
!               almo_scf_env%matrix_ks_0deloc(ispin),&
!               0.0_dp,Fon,filter_eps=almo_scf_env%eps_filter)
!
!       CALL dbcsr_multiply("N","N",1.0_dp,Fon,&
!               almo_scf_env%matrix_v_full_blk(ispin),&
!               0.0_dp,Fov,filter_eps=almo_scf_env%eps_filter)
!
!       CALL dbcsr_multiply("N","N",1.0_dp,Fon,&
!               almo_scf_env%matrix_t_blk(ispin),&
!               0.0_dp,temp1_oo,filter_eps=almo_scf_env%eps_filter)
!
!       CALL dbcsr_multiply("N","N",1.0_dp,temp1_oo,&
!               almo_scf_env%matrix_sigma_inv(ispin),&
!               0.0_dp,temp2_oo,filter_eps=almo_scf_env%eps_filter)
!       CALL dbcsr_release(temp1_oo)
!
!       CALL dbcsr_multiply("T","N",1.0_dp,almo_scf_env%matrix_t_blk(ispin),&
!               almo_scf_env%matrix_s(1),&
!               0.0_dp,Fon,filter_eps=almo_scf_env%eps_filter)
!
!       CALL dbcsr_multiply("N","N",1.0_dp,Fon,&
!               almo_scf_env%matrix_v_full_blk(ispin),&
!               0.0_dp,Fov_filtered,filter_eps=almo_scf_env%eps_filter)
!       CALL dbcsr_release(Fon)
!
!       CALL dbcsr_multiply("N","N",-1.0_dp,temp2_oo,&
!               Fov_filtered,&
!               1.0_dp,Fov,filter_eps=almo_scf_env%eps_filter)
!       CALL dbcsr_release(temp2_oo)
!
!       CALL dbcsr_multiply("N","N",1.0_dp,almo_scf_env%matrix_sigma_inv(ispin),&
!               Fov,0.0_dp,Fov_filtered,filter_eps=almo_scf_env%eps_filter)
!
!       CALL dbcsr_multiply("N","N",1.0_dp,Fov_filtered,&
!               matrix_sigma_vv_full_sqrt_inv,&
!               0.0_dp,Fov,filter_eps=almo_scf_env%eps_filter)
!       !CALL dbcsr_copy(Fov,Fov_filtered)
!CALL dbcsr_print(Fov)
!
!       IF (safe_mode) THEN
!          CALL dbcsr_init(Fov_original)
!          CALL dbcsr_create(Fov_original,template=Fov)
!          CALL dbcsr_copy(Fov_original,Fov)
!       ENDIF
!
!!! remove diagonal blocks
!!CALL dbcsr_iterator_start(iter,Fov)
!!DO WHILE (dbcsr_iterator_blocks_left(iter))
!!
!!   CALL dbcsr_iterator_next_block(iter,iblock_row,iblock_col,data_p,&
!!           row_size=iblock_row_size,col_size=iblock_col_size)
!!
!!   IF (iblock_row.eq.iblock_col) data_p(:,:)=0.0_dp
!!
!!ENDDO
!!CALL dbcsr_iterator_stop(iter)
!!CALL dbcsr_finalize(Fov)
!
!!! perform svd of blocks
!!!!! THIS ROUTINE WORKS ONLY ON ONE CPU AND ONLY FOR 2 MOLECULES !!!
!!CALL dbcsr_init(temp_u_v_full_blk)
!!CALL dbcsr_create(temp_u_v_full_blk,&
!!        template=almo_scf_env%matrix_vv_full_blk(ispin),&
!!        matrix_type=dbcsr_type_no_symmetry)
!!
!!CALL dbcsr_work_create(temp_u_v_full_blk,&
!!        work_mutable=.TRUE.)
!!CALL dbcsr_iterator_start(iter,Fov)
!!DO WHILE (dbcsr_iterator_blocks_left(iter))
!!
!!   CALL dbcsr_iterator_next_block(iter,iblock_row,iblock_col,data_p,&
!!           row_size=iblock_row_size,col_size=iblock_col_size)
!!
!!   IF (iblock_row.ne.iblock_col) THEN
!!
!!      ! Prepare data
!!      allocate(eigenvalues(min(iblock_row_size,iblock_col_size)))
!!      allocate(data_copy(iblock_row_size,iblock_col_size))
!!      allocate(left_vectors(iblock_row_size,iblock_row_size))
!!      allocate(right_vectors(iblock_col_size,iblock_col_size))
!!      data_copy(:,:)=data_p(:,:)
!!
!!      ! Query the optimal workspace for dgesvd
!!      LWORK = -1
!!      allocate(WORK(MAX(1,LWORK)))
!!      CALL DGESVD('N','A',iblock_row_size,iblock_col_size,data_copy,&
!!              iblock_row_size,eigenvalues,left_vectors,iblock_row_size,&
!!              right_vectors,iblock_col_size,WORK,LWORK,INFO)
!!      LWORK = INT(WORK( 1 ))
!!      deallocate(WORK)
!!
!!      ! Allocate the workspace and perform svd
!!      allocate(WORK(MAX(1,LWORK)))
!!      CALL DGESVD('N','A',iblock_row_size,iblock_col_size,data_copy,&
!!              iblock_row_size,eigenvalues,left_vectors,iblock_row_size,&
!!              right_vectors,iblock_col_size,WORK,LWORK,INFO)
!!      deallocate(WORK)
!!      IF( INFO.NE.0 ) THEN
!!         CPABORT("DGESVD failed")
!!      END IF
!!
!!      ! copy right singular vectors into a unitary matrix
!!      NULLIFY (p_new_block)
!!      CALL dbcsr_reserve_block2d(temp_u_v_full_blk,iblock_col,iblock_col,p_new_block)
!!      CPASSERT(ASSOCIATED(p_new_block))
!!      p_new_block(:,:) = right_vectors(:,:)
!!
!!      deallocate(eigenvalues)
!!      deallocate(data_copy)
!!      deallocate(left_vectors)
!!      deallocate(right_vectors)
!!
!!   ENDIF
!!ENDDO
!!CALL dbcsr_iterator_stop(iter)
!!CALL dbcsr_finalize(temp_u_v_full_blk)
!!!CALL dbcsr_print(temp_u_v_full_blk)
!!CALL dbcsr_multiply("N","T",1.0_dp,Fov,temp_u_v_full_blk,&
!!        0.0_dp,Fov_filtered,filter_eps=almo_scf_env%eps_filter)
!!
!!CALL dbcsr_copy(Fov,Fov_filtered)
!!CALL dbcsr_print(Fov)
!
!       !!!!!!!!!!!!!!!!!!!
!       ! 2. Initialize variables
!
!       ! temp space
!       CALL dbcsr_init(temp0_ov)
!       CALL dbcsr_create(temp0_ov,&
!               template=almo_scf_env%matrix_ov_full(ispin))
!
!       ! current unitary matrix
!       CALL dbcsr_init(U_blk)
!       CALL dbcsr_create(U_blk,&
!               template=almo_scf_env%matrix_vv_full_blk(ispin),&
!               matrix_type=dbcsr_type_no_symmetry)
!
!       ! unitary matrix accumulator
!       CALL dbcsr_init(U_blk_tot)
!       CALL dbcsr_create(U_blk_tot,&
!               template=almo_scf_env%matrix_vv_full_blk(ispin),&
!               matrix_type=dbcsr_type_no_symmetry)
!       CALL dbcsr_add_on_diag(U_blk_tot,1.0_dp)
!
!!CALL dbcsr_add_on_diag(U_blk,1.0_dp)
!!CALL dbcsr_multiply("N","T",1.0_dp,U_blk,temp_u_v_full_blk,&
!!        0.0_dp,U_blk_tot,filter_eps=almo_scf_env%eps_filter)
!!
!!CALL dbcsr_release(temp_u_v_full_blk)
!
!       ! init gradient
!       CALL dbcsr_init(grad_blk)
!       CALL dbcsr_create(grad_blk,&
!               template=almo_scf_env%matrix_vv_full_blk(ispin),&
!               matrix_type=dbcsr_type_no_symmetry)
!
!       ! init step matrix
!       CALL dbcsr_init(step_blk)
!       CALL dbcsr_create(step_blk,&
!               template=almo_scf_env%matrix_vv_full_blk(ispin),&
!               matrix_type=dbcsr_type_no_symmetry)
!
!       ! "retain discarded" filter (0.0 - retain, 1.0 - discard)
!       CALL dbcsr_init(matrix_filter)
!       CALL dbcsr_create(matrix_filter,&
!               template=almo_scf_env%matrix_ov_full(ispin))
!       ! copy Fov into the filter matrix temporarily
!       ! so we know which blocks contain significant elements
!       CALL dbcsr_copy(matrix_filter,Fov)
!
!       ! fill out filter elements block-by-block
!       CALL dbcsr_iterator_start(iter,matrix_filter)
!       DO WHILE (dbcsr_iterator_blocks_left(iter))
!
!          CALL dbcsr_iterator_next_block(iter,iblock_row,iblock_col,data_p,&
!                  row_size=iblock_row_size,col_size=iblock_col_size)
!
!          retained_v=almo_scf_env%nvirt_of_domain(iblock_col,ispin)
!
!          data_p(:,1:retained_v)=0.0_dp
!          data_p(:,(retained_v+1):iblock_col_size)=1.0_dp
!
!       ENDDO
!       CALL dbcsr_iterator_stop(iter)
!       CALL dbcsr_finalize(matrix_filter)
!
!       ! apply the filter
!       CALL dbcsr_hadamard_product(Fov,matrix_filter,Fov_filtered)
!
!       !!!!!!!!!!!!!!!!!!!!!
!       ! 3. start iterative minimization of the elements to be discarded
!       iteration=0
!       converged=.FALSE.
!       prepare_to_exit=.FALSE.
!       DO
!
!          iteration=iteration+1
!
!          !!!!!!!!!!!!!!!!!!!!!!!!!
!          ! 4. compute the gradient
!          CALL dbcsr_set(grad_blk,0.0_dp)
!          ! create the diagonal blocks only
!          CALL dbcsr_add_on_diag(grad_blk,1.0_dp)
!
!          CALL dbcsr_multiply("T","N",2.0_dp,Fov_filtered,Fov,&
!                  0.0_dp,grad_blk,retain_sparsity=.TRUE.,&
!                  filter_eps=almo_scf_env%eps_filter)
!          CALL dbcsr_multiply("T","N",-2.0_dp,Fov,Fov_filtered,&
!                  1.0_dp,grad_blk,retain_sparsity=.TRUE.,&
!                  filter_eps=almo_scf_env%eps_filter)
!
!          !!!!!!!!!!!!!!!!!!!!!!!
!          ! 5. check convergence
!          obj_function = 0.5_dp*(dbcsr_frobenius_norm(Fov_filtered))**2
!          grad_norm = dbcsr_frobenius_norm(grad_blk)
!          converged=(grad_norm.lt.almo_scf_env%truncate_v_eps_convergence)
!          IF (converged.OR.(iteration.ge.almo_scf_env%truncate_v_max_iter)) THEN
!             prepare_to_exit=.TRUE.
!          ENDIF
!
!          IF (.NOT.prepare_to_exit) THEN
!
!             !!!!!!!!!!!!!!!!!!!!!!!
!             ! 6. perform steps in the direction of the gradient
!             !    a. first, perform a trial step to "see" the parameters
!             !       of the parabola along the gradient:
!             !       a0 * x^2 + b0 * x + c0
!             !    b. then perform the step to the bottom of the parabola
!
!             ! get c0
!             c0 = obj_function
!             ! get b0 <= d_f/d_alpha along grad
!             !!!CALL dbcsr_multiply("N","N",4.0_dp,Fov,grad_blk,&
!             !!!        0.0_dp,temp0_ov,&
!             !!!        filter_eps=almo_scf_env%eps_filter)
!             !!!CALL dbcsr_dot(Fov_filtered,temp0_ov,b0)
!
!             alpha=almo_scf_env%truncate_v_trial_step_size
!
!             line_search_step_last=3
!             DO line_search_step=1,line_search_step_last
!                CALL dbcsr_copy(step_blk,grad_blk)
!                CALL dbcsr_scale(step_blk,-1.0_dp*alpha)
!                CALL generator_to_unitary(step_blk,U_blk,&
!                        almo_scf_env%eps_filter)
!                CALL dbcsr_multiply("N","N",1.0_dp,Fov,U_blk,0.0_dp,temp0_ov,&
!                        filter_eps=almo_scf_env%eps_filter)
!                CALL dbcsr_hadamard_product(temp0_ov,matrix_filter,&
!                        Fov_filtered)
!
!                obj_function_new = 0.5_dp*(dbcsr_frobenius_norm(Fov_filtered))**2
!                IF (line_search_step.eq.1) THEN
!                   ff1 = obj_function_new
!                   step1 = alpha
!                ELSE IF (line_search_step.eq.2) THEN
!                   ff2 = obj_function_new
!                   step2 = alpha
!                ENDIF
!
!                IF (unit_nr>0.AND.(line_search_step.ne.line_search_step_last)) THEN
!                   WRITE(unit_nr,'(T6,A,1X,I3,1X,F10.3,E12.3,E12.3,E12.3)') &
!                         "JOINT_SVD_lin",&
!                         iteration,&
!                         alpha,&
!                         obj_function,&
!                         obj_function_new,&
!                         obj_function_new-obj_function
!                ENDIF
!
!                IF (line_search_step.eq.1) THEN
!                   alpha=2.0_dp*alpha
!                ENDIF
!                IF (line_search_step.eq.2) THEN
!                   a0 = ((ff1-c0)/step1 - (ff2-c0)/step2) / (step1 - step2)
!                   b0 = (ff1-c0)/step1 - a0*step1
!                   ! step size in to the bottom of "the parabola"
!                   alpha=-b0/(2.0_dp*a0)
!                   ! update the default step size
!                   almo_scf_env%truncate_v_trial_step_size=alpha
!                ENDIF
!                !!!IF (line_search_step.eq.1) THEN
!                !!!   a0 = (obj_function_new - b0 * alpha - c0) / (alpha*alpha)
!                !!!   ! step size in to the bottom of "the parabola"
!                !!!   alpha=-b0/(2.0_dp*a0)
!                !!!   !IF (alpha.gt.10.0_dp) alpha=10.0_dp
!                !!!ENDIF
!
!             ENDDO
!
!             ! update Fov and U_blk_tot (use grad_blk as tmp storage)
!             CALL dbcsr_copy(Fov,temp0_ov)
!             CALL dbcsr_multiply("N","N",1.0_dp,U_blk_tot,U_blk,&
!                     0.0_dp,grad_blk,&
!                     filter_eps=almo_scf_env%eps_filter)
!             CALL dbcsr_copy(U_blk_tot,grad_blk)
!
!          ENDIF
!
!          t2 = m_walltime()
!
!          IF (unit_nr>0) THEN
!             WRITE(unit_nr,'(T6,A,1X,I3,1X,F10.3,E12.3,E12.3,E12.3,E12.3,F10.3)') &
!                   "JOINT_SVD_itr",&
!                   iteration,&
!                   alpha,&
!                   obj_function,&
!                   obj_function_new,&
!                   obj_function_new-obj_function,&
!                   grad_norm,&
!                   t2-t1
!                   !(flop1+flop2)/(1.0E6_dp*(t2-t1))
!             CALL m_flush(unit_nr)
!          ENDIF
!
!          t1 = m_walltime()
!
!          IF (prepare_to_exit) EXIT
!
!       ENDDO ! stop iterations
!
!       IF (safe_mode) THEN
!          CALL dbcsr_multiply("N","N",1.0_dp,Fov_original,&
!                  U_blk_tot,0.0_dp,temp0_ov,&
!                  filter_eps=almo_scf_env%eps_filter)
!CALL dbcsr_print(temp0_ov)
!          CALL dbcsr_hadamard_product(temp0_ov,matrix_filter,&
!                  Fov_filtered)
!          obj_function_new = 0.5_dp*(dbcsr_frobenius_norm(Fov_filtered))**2
!
!          IF (unit_nr>0) THEN
!             WRITE(unit_nr,'(T6,A,1X,E12.3)') &
!                   "SANITY CHECK:",&
!                   obj_function_new
!             CALL m_flush(unit_nr)
!          ENDIF
!
!          CALL dbcsr_release(Fov_original)
!       ENDIF
!
!       CALL dbcsr_release(temp0_ov)
!       CALL dbcsr_release(U_blk)
!       CALL dbcsr_release(grad_blk)
!       CALL dbcsr_release(step_blk)
!       CALL dbcsr_release(matrix_filter)
!       CALL dbcsr_release(Fov)
!       CALL dbcsr_release(Fov_filtered)
!
!       ! compute rotated virtual orbitals
!       CALL dbcsr_init(v_full_tmp)
!       CALL dbcsr_create(v_full_tmp,&
!               template=almo_scf_env%matrix_v_full_blk(ispin),&
!               matrix_type=dbcsr_type_no_symmetry)
!       CALL dbcsr_multiply("N","N",1.0_dp,&
!               v_full_new,&
!               matrix_sigma_vv_full_sqrt_inv,0.0_dp,v_full_tmp,&
!               filter_eps=almo_scf_env%eps_filter)
!       CALL dbcsr_multiply("N","N",1.0_dp,&
!               v_full_tmp,&
!               U_blk_tot,0.0_dp,v_full_new,&
!               filter_eps=almo_scf_env%eps_filter)
!
!       CALL dbcsr_release(matrix_sigma_vv_full_sqrt_inv)
!       CALL dbcsr_release(v_full_tmp)
!       CALL dbcsr_release(U_blk_tot)
!
!!!!! orthogonalized virtuals are not blocked
!       ! copy new virtuals into the truncated matrix
!       !CALL dbcsr_work_create(almo_scf_env%matrix_v_blk(ispin),&
!       CALL dbcsr_work_create(almo_scf_env%matrix_v(ispin),&
!               work_mutable=.TRUE.)
!       CALL dbcsr_iterator_start(iter,v_full_new)
!       DO WHILE (dbcsr_iterator_blocks_left(iter))
!
!          CALL dbcsr_iterator_next_block(iter,iblock_row,iblock_col,data_p,&
!                  row_size=iblock_row_size,col_size=iblock_col_size)
!
!          retained_v=almo_scf_env%nvirt_of_domain(iblock_col,ispin)
!
!          NULLIFY (p_new_block)
!          !CALL dbcsr_reserve_block2d(almo_scf_env%matrix_v_blk(ispin),&
!          CALL dbcsr_reserve_block2d(almo_scf_env%matrix_v(ispin),&
!                  iblock_row,iblock_col,p_new_block)
!          CPASSERT(ASSOCIATED(p_new_block))
!          CPASSERT(retained_v.gt.0)
!          p_new_block(:,:) = data_p(:,1:retained_v)
!
!       ENDDO ! iterator
!       CALL dbcsr_iterator_stop(iter)
!       !!CALL dbcsr_finalize(almo_scf_env%matrix_v_blk(ispin))
!       CALL dbcsr_finalize(almo_scf_env%matrix_v(ispin))
!
!       CALL dbcsr_release(v_full_new)
!
!    ENDDO ! ispin
!
!    CALL timestop(handle)
!
!  END SUBROUTINE truncate_subspace_v_blk

! *****************************************************************************
!> \brief Compute the gradient wrt the main variable (e.g. Theta, X)
!> \param m_grad_out ...
!> \param m_ks ...
!> \param m_s ...
!> \param m_t ...
!> \param m_t0 ...
!> \param m_siginv ...
!> \param m_quench_t ...
!> \param m_FTsiginv ...
!> \param m_siginvTFTsiginv ...
!> \param m_ST ...
!> \param m_STsiginv0 ...
!> \param m_theta ...
!> \param domain_s_inv ...
!> \param domain_r_down ...
!> \param cpu_of_domain ...
!> \param domain_map ...
!> \param assume_t0_q0x ...
!> \param optimize_theta ...
!> \param normalize_orbitals ...
!> \param penalty_occ_vol ...
!> \param penalty_occ_vol_prefactor ...
!> \param envelope_amplitude ...
!> \param eps_filter ...
!> \param spin_factor ...
!> \param special_case ...
!> \param m_sig_sqrti_ii ...
!> \par History
!>       2015.03 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE compute_gradient(m_grad_out, m_ks, m_s, m_t, m_t0, &
                               m_siginv, m_quench_t, m_FTsiginv, m_siginvTFTsiginv, m_ST, m_STsiginv0, &
                               m_theta, domain_s_inv, domain_r_down, &
                               cpu_of_domain, domain_map, assume_t0_q0x, optimize_theta, &
                               normalize_orbitals, penalty_occ_vol, &
                               penalty_occ_vol_prefactor, envelope_amplitude, eps_filter, spin_factor, &
                               special_case, m_sig_sqrti_ii)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: m_grad_out
      TYPE(dbcsr_type), INTENT(IN)                       :: m_ks, m_s, m_t, m_t0, m_siginv, &
                                                            m_quench_t, m_FTsiginv, &
                                                            m_siginvTFTsiginv, m_ST, m_STsiginv0, &
                                                            m_theta
      TYPE(domain_submatrix_type), DIMENSION(:), &
         INTENT(IN)                                      :: domain_s_inv, domain_r_down
      INTEGER, DIMENSION(:), INTENT(IN)                  :: cpu_of_domain
      TYPE(domain_map_type), INTENT(IN)                  :: domain_map
      LOGICAL, INTENT(IN)                                :: assume_t0_q0x, optimize_theta, &
                                                            normalize_orbitals, penalty_occ_vol
      REAL(KIND=dp), INTENT(IN)                          :: penalty_occ_vol_prefactor, &
                                                            envelope_amplitude, eps_filter, &
                                                            spin_factor
      INTEGER, INTENT(IN)                                :: special_case
      TYPE(dbcsr_type), INTENT(IN), OPTIONAL             :: m_sig_sqrti_ii

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_gradient', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, nao
      REAL(KIND=dp)                                      :: energy_g_norm, penalty_occ_vol_g_norm
      TYPE(dbcsr_type)                                   :: m_tmp_no_1, m_tmp_no_2, m_tmp_no_3, &
                                                            m_tmp_oo_1, m_tmp_oo_2

!TYPE(dbcsr_type), INTENT(INOUT)          :: m_tmp_no_2

      CALL timeset(routineN, handle)

      IF (normalize_orbitals .AND. (.NOT. PRESENT(m_sig_sqrti_ii))) THEN
         CPABORT("Normalization matrix is required")
      ENDIF

      ! use this otherways unused variables
      CALL dbcsr_get_info(matrix=m_ks, nfullrows_total=nao)
      CALL dbcsr_get_info(matrix=m_s, nfullrows_total=nao)
      CALL dbcsr_get_info(matrix=m_t, nfullrows_total=nao)

      CALL dbcsr_create(m_tmp_no_1, &
                        template=m_quench_t, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(m_tmp_no_2, &
                        template=m_quench_t, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(m_tmp_no_3, &
                        template=m_quench_t, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(m_tmp_oo_1, &
                        template=m_siginv, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(m_tmp_oo_2, &
                        template=m_siginv, &
                        matrix_type=dbcsr_type_no_symmetry)

      ! do d_E/d_T first
      !IF (.NOT.PRESENT(m_FTsiginv)) THEN
      !   CALL dbcsr_multiply("N","N",1.0_dp,&
      !           m_ks,&
      !           m_t,&
      !           0.0_dp,m_tmp_no_1,&
      !           filter_eps=eps_filter)
      !   CALL dbcsr_multiply("N","N",1.0_dp,&
      !           m_tmp_no_1,&
      !           m_siginv,&
      !           0.0_dp,m_FTsiginv,&
      !           filter_eps=eps_filter)
      !ENDIF

      CALL dbcsr_copy(m_tmp_no_2, m_quench_t)
      CALL dbcsr_copy(m_tmp_no_2, m_FTsiginv, keep_sparsity=.TRUE.)

      !IF (.NOT.PRESENT(m_siginvTFTsiginv)) THEN
      !   CALL dbcsr_multiply("T","N",1.0_dp,&
      !           m_t,&
      !           m_FTsiginv,&
      !           0.0_dp,m_tmp_oo_1,&
      !           filter_eps=eps_filter)
      !   CALL dbcsr_multiply("N","N",1.0_dp,&
      !           m_siginv,&
      !           m_tmp_oo_1,&
      !           0.0_dp,m_siginvTFTsiginv,&
      !           filter_eps=eps_filter)
      !ENDIF

      !IF (.NOT.PRESENT(m_ST)) THEN
      !   CALL dbcsr_multiply("N","N",1.0_dp,&
      !           m_s,&
      !           m_t,&
      !           0.0_dp,m_ST,&
      !           filter_eps=eps_filter)
      !ENDIF

      CALL dbcsr_multiply("N", "N", -1.0_dp, &
                          m_ST, &
                          m_siginvTFTsiginv, &
                          1.0_dp, m_tmp_no_2, &
                          retain_sparsity=.TRUE.)
      CALL dbcsr_scale(m_tmp_no_2, 2.0_dp*spin_factor)

      ! add penalty on the occupied volume: det(sigma)
      IF (penalty_occ_vol) THEN
         !RZK-warning CALL dbcsr_multiply("N","N",&
         !RZK-warning         penalty_occ_vol_prefactor,&
         !RZK-warning         m_ST,&
         !RZK-warning         m_siginv,&
         !RZK-warning         1.0_dp,m_tmp_no_2,&
         !RZK-warning         retain_sparsity=.TRUE.,&
         !RZK-warning         )
         CALL dbcsr_copy(m_tmp_no_1, m_quench_t)
         CALL dbcsr_multiply("N", "N", &
                             penalty_occ_vol_prefactor, &
                             m_ST, &
                             m_siginv, &
                             0.0_dp, m_tmp_no_1, &
                             retain_sparsity=.TRUE.)
! this norm does not contain the normalization factors
         CALL dbcsr_norm(m_tmp_no_1, dbcsr_norm_maxabsnorm, &
                         norm_scalar=penalty_occ_vol_g_norm)
         CALL dbcsr_norm(m_tmp_no_2, dbcsr_norm_maxabsnorm, &
                         norm_scalar=energy_g_norm)
         !WRITE (*, "(A30,2F20.10)") "Energy/penalty g norms (no norm): ", energy_g_norm, penalty_occ_vol_g_norm
         CALL dbcsr_add(m_tmp_no_2, m_tmp_no_1, 1.0_dp, 1.0_dp)
      ENDIF

      ! take into account the factor from the normalization constraint
      IF (normalize_orbitals) THEN

         ! G = ( G - ST.[tr(T).G]_ii ) . [sig_sqrti]_ii
         ! this expression can be simplified to
         ! G = ( G - c0*ST ) . [sig_sqrti]_ii
         ! where c0 = penalty_occ_vol_prefactor
         ! This is because tr(T).G_Energy = 0 and
         !                 tr(T).G_Penalty = c0*I
         CALL dbcsr_copy(m_tmp_no_1, m_quench_t)
         CALL dbcsr_copy(m_tmp_no_1, m_ST, keep_sparsity=.TRUE.)
         CALL dbcsr_add(m_tmp_no_2, m_tmp_no_1, 1.0_dp, -penalty_occ_vol_prefactor)
         CALL dbcsr_copy(m_tmp_no_1, m_quench_t)
         CALL dbcsr_multiply("N", "N", 1.0_dp, &
                             m_tmp_no_2, &
                             m_sig_sqrti_ii, &
                             0.0_dp, m_tmp_no_1, &
                             retain_sparsity=.TRUE.)

         !!! slower way of taking the norm into account
         !!CALL dbcsr_copy(m_tmp_no_1,m_tmp_no_2)
         !!CALL dbcsr_multiply("N","N",1.0_dp,&
         !!        m_tmp_no_2,&
         !!        m_sig_sqrti_ii,&
         !!        0.0_dp,m_tmp_no_1,&
         !!        retain_sparsity=.TRUE.,&
         !!        )
         !!
         !!! get [tr(T).G]_ii
         !!CALL dbcsr_copy(m_tmp_oo_1,m_sig_sqrti_ii)
         !!CALL dbcsr_multiply("T","N",1.0_dp,&
         !!        m_t,&
         !!        m_tmp_no_2,&
         !!        0.0_dp,m_tmp_oo_1,&
         !!        retain_sparsity=.TRUE.,&
         !!        )
         !!CALL dbcsr_get_info(m_sig_sqrti_ii, nfullrows_total=dim0 )
         !!ALLOCATE(tg_diagonal(dim0))
         !!CALL dbcsr_get_diag(m_tmp_oo_1,tg_diagonal)
         !!CALL dbcsr_set(m_tmp_oo_1,0.0_dp)
         !!CALL dbcsr_set_diag(m_tmp_oo_1,tg_diagonal)
         !!DEALLOCATE(tg_diagonal)
         !!
         !!CALL dbcsr_multiply("N","N",1.0_dp,&
         !!        m_sig_sqrti_ii,&
         !!        m_tmp_oo_1,&
         !!        0.0_dp,m_tmp_oo_2,&
         !!        filter_eps=eps_filter,&
         !!        )
         !!CALL dbcsr_multiply("N","N",-1.0_dp,&
         !!        m_ST,&
         !!        m_tmp_oo_2,&
         !!        1.0_dp,m_tmp_no_1,&
         !!        retain_sparsity=.TRUE.,&
         !!        )

      ELSE

         CALL dbcsr_copy(m_tmp_no_1, m_tmp_no_2)

      ENDIF ! normalize_orbitals

      ! project out the occupied space from the gradient
      IF (assume_t0_q0x) THEN
         IF (special_case .EQ. xalmo_case_fully_deloc) THEN
            CALL dbcsr_copy(m_grad_out, m_tmp_no_1)
            CALL dbcsr_multiply("T", "N", 1.0_dp, &
                                m_t0, &
                                m_grad_out, &
                                0.0_dp, m_tmp_oo_1, &
                                filter_eps=eps_filter)
            CALL dbcsr_multiply("N", "N", -1.0_dp, &
                                m_STsiginv0, &
                                m_tmp_oo_1, &
                                1.0_dp, m_grad_out, &
                                filter_eps=eps_filter)
         ELSE IF (special_case .EQ. xalmo_case_block_diag) THEN
            CPABORT("Cannot project the zero-order space from itself")
         ELSE
            ! no special case: normal xALMOs
            CALL apply_domain_operators( &
               matrix_in=m_tmp_no_1, &
               matrix_out=m_grad_out, &
               operator2=domain_r_down(:), &
               operator1=domain_s_inv(:), &
               dpattern=m_quench_t, &
               map=domain_map, &
               node_of_domain=cpu_of_domain, &
               my_action=1, &
               filter_eps=eps_filter, &
               !matrix_trimmer=,&
               use_trimmer=.FALSE.)
         ENDIF ! my_special_case
         CALL dbcsr_copy(m_tmp_no_1, m_grad_out)
      ENDIF

      !! check whether the gradient lies entirely in R or Q
      !CALL dbcsr_multiply("T","N",1.0_dp,&
      !        m_t,&
      !        m_tmp_no_1,&
      !        0.0_dp,m_tmp_oo_1,&
      !        filter_eps=eps_filter,&
      !        )
      !CALL dbcsr_multiply("N","N",1.0_dp,&
      !        m_siginv,&
      !        m_tmp_oo_1,&
      !        0.0_dp,m_tmp_oo_2,&
      !        filter_eps=eps_filter,&
      !        )
      !CALL dbcsr_copy(m_tmp_no_2,m_tmp_no_1)
      !CALL dbcsr_multiply("N","N",-1.0_dp,&
      !        m_ST,&
      !        m_tmp_oo_2,&
      !        1.0_dp,m_tmp_no_2,&
      !        retain_sparsity=.TRUE.,&
      !        )
      !CALL dbcsr_norm(m_tmp_no_2, dbcsr_norm_maxabsnorm,&
      !        norm_scalar=penalty_occ_vol_g_norm, )
      !WRITE(*,"(A50,2F20.10)") "Virtual-space projection of the gradient", penalty_occ_vol_g_norm
      !CALL dbcsr_add(m_tmp_no_2,m_tmp_no_1,1.0_dp,-1.0_dp)
      !CALL dbcsr_norm(m_tmp_no_2, dbcsr_norm_maxabsnorm,&
      !        norm_scalar=penalty_occ_vol_g_norm, )
      !WRITE(*,"(A50,2F20.10)") "Occupied-space projection of the gradient", penalty_occ_vol_g_norm
      !CALL dbcsr_norm(m_tmp_no_1, dbcsr_norm_maxabsnorm,&
      !        norm_scalar=penalty_occ_vol_g_norm, )
      !WRITE(*,"(A50,2F20.10)") "Full gradient", penalty_occ_vol_g_norm

      ! transform d_E/d_T to d_E/d_theta
      IF (optimize_theta) THEN
         CALL dbcsr_copy(m_tmp_no_2, m_theta)
         CALL dbcsr_function_of_elements(m_tmp_no_2, &
                                         !func=dbcsr_func_cos,&
                                         func=dbcsr_func_dtanh, &
                                         a0=0.0_dp, &
                                         a1=1.0_dp/envelope_amplitude)
         CALL dbcsr_scale(m_tmp_no_2, &
                          envelope_amplitude)
         CALL dbcsr_set(m_tmp_no_3, 0.0_dp)
         CALL dbcsr_filter(m_tmp_no_3, eps=eps_filter)
         CALL dbcsr_hadamard_product(m_tmp_no_1, &
                                     m_tmp_no_2, &
                                     m_tmp_no_3, &
                                     b_assume_value=1.0_dp)
         CALL dbcsr_hadamard_product(m_tmp_no_3, &
                                     m_quench_t, &
                                     m_grad_out)
      ELSE ! simply copy
         CALL dbcsr_hadamard_product(m_tmp_no_1, &
                                     m_quench_t, &
                                     m_grad_out)
      ENDIF
      CALL dbcsr_filter(m_grad_out, eps=eps_filter)

      CALL dbcsr_release(m_tmp_no_1)
      CALL dbcsr_release(m_tmp_no_2)
      CALL dbcsr_release(m_tmp_no_3)
      CALL dbcsr_release(m_tmp_oo_1)
      CALL dbcsr_release(m_tmp_oo_2)

      CALL timestop(handle)

   END SUBROUTINE compute_gradient

! *****************************************************************************
!> \brief Serial code that prints matrices readable by Mathematica
!> \param matrix - matrix to print
!> \param filename ...
!> \par History
!>       2015.05 created [Rustam Z. Khaliullin]
!> \author Rustam Z. Khaliullin
! **************************************************************************************************
   SUBROUTINE print_mathematica_matrix(matrix, filename)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      CHARACTER(len=*)                                   :: filename

      CHARACTER(len=*), PARAMETER :: routineN = 'print_mathematica_matrix', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=20)                                  :: formatstr, Scols
      INTEGER                                            :: col, fiunit, handle, hori_offset, jj, &
                                                            nblkcols_tot, nblkrows_tot, Ncols, &
                                                            ncores, Nrows, row, unit_nr, &
                                                            vert_offset
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ao_block_sizes, mo_block_sizes
      INTEGER, DIMENSION(:), POINTER                     :: ao_blk_sizes, mo_blk_sizes
      LOGICAL                                            :: found
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: H
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block_p
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_distribution_type)                      :: dist
      TYPE(dbcsr_type)                                   :: matrix_asym

      CALL timeset(routineN, handle)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      ! serial code only
      CALL dbcsr_get_info(matrix, distribution=dist)
      CALL dbcsr_distribution_get(dist, numnodes=ncores)
      IF (ncores .GT. 1) THEN
         CPABORT("mathematica files: serial code only")
      ENDIF

      nblkrows_tot = dbcsr_nblkrows_total(matrix)
      nblkcols_tot = dbcsr_nblkcols_total(matrix)
      CPASSERT(nblkrows_tot == nblkcols_tot)
      CALL dbcsr_get_info(matrix, row_blk_size=ao_blk_sizes)
      CALL dbcsr_get_info(matrix, col_blk_size=mo_blk_sizes)
      ALLOCATE (mo_block_sizes(nblkcols_tot), ao_block_sizes(nblkcols_tot))
      mo_block_sizes(:) = mo_blk_sizes(:)
      ao_block_sizes(:) = ao_blk_sizes(:)

      CALL dbcsr_create(matrix_asym, &
                        template=matrix, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_desymmetrize(matrix, matrix_asym)

      Ncols = SUM(mo_block_sizes)
      Nrows = SUM(ao_block_sizes)
      ALLOCATE (H(Nrows, Ncols))
      H(:, :) = 0.0_dp

      hori_offset = 0
      DO col = 1, nblkcols_tot

         vert_offset = 0
         DO row = 1, nblkrows_tot

            CALL dbcsr_get_block_p(matrix_asym, row, col, block_p, found)
            IF (found) THEN

               H(vert_offset + 1:vert_offset + ao_block_sizes(row), &
                 hori_offset + 1:hori_offset + mo_block_sizes(col)) &
                  = block_p(:, :)

            ENDIF

            vert_offset = vert_offset + ao_block_sizes(row)

         ENDDO

         hori_offset = hori_offset + mo_block_sizes(col)

      ENDDO ! loop over electron blocks

      CALL dbcsr_release(matrix_asym)

      IF (unit_nr > 0) THEN
         CALL open_file(filename, unit_number=fiunit, file_status='REPLACE')
         WRITE (Scols, "(I10)") Ncols
         formatstr = "("//TRIM(Scols)//"E27.17)"
         DO jj = 1, Nrows
            WRITE (fiunit, formatstr) H(jj, :)
         ENDDO
         CALL close_file(fiunit)
      ENDIF

      DEALLOCATE (mo_block_sizes)
      DEALLOCATE (ao_block_sizes)
      DEALLOCATE (H)

      CALL timestop(handle)

   END SUBROUTINE print_mathematica_matrix

! *****************************************************************************
!> \brief Compute MO coeffs from the main optimized variable (e.g. Theta, X)
!> \param m_var_in ...
!> \param m_t_out ...
!> \param m_quench_t ...
!> \param m_t0 ...
!> \param m_siginv ...
!> \param m_STsiginv0 ...
!> \param m_s ...
!> \param m_sig_sqrti_ii_out ...
!> \param domain_r_down ...
!> \param domain_s_inv ...
!> \param domain_map ...
!> \param cpu_of_domain ...
!> \param assume_t0_q0x ...
!> \param just_started ...
!> \param optimize_theta ...
!> \param normalize_orbitals ...
!> \param envelope_amplitude ...
!> \param eps_filter ...
!> \param special_case ...
!> \param nocc_of_domain ...
!> \param order_lanczos ...
!> \param eps_lanczos ...
!> \param max_iter_lanczos ...
!> \par History
!>       2015.03 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE compute_xalmos_from_main_var(m_var_in, m_t_out, m_quench_t, &
                                           m_t0, m_siginv, m_STsiginv0, m_s, m_sig_sqrti_ii_out, domain_r_down, &
                                           domain_s_inv, domain_map, cpu_of_domain, assume_t0_q0x, just_started, &
                                           optimize_theta, normalize_orbitals, envelope_amplitude, eps_filter, &
                                           special_case, nocc_of_domain, order_lanczos, eps_lanczos, max_iter_lanczos)

      TYPE(dbcsr_type), INTENT(IN)                       :: m_var_in
      TYPE(dbcsr_type), INTENT(INOUT)                    :: m_t_out, m_quench_t
      TYPE(dbcsr_type), INTENT(IN)                       :: m_t0
      TYPE(dbcsr_type), INTENT(INOUT)                    :: m_siginv
      TYPE(dbcsr_type), INTENT(IN)                       :: m_STsiginv0, m_s
      TYPE(dbcsr_type), INTENT(INOUT)                    :: m_sig_sqrti_ii_out
      TYPE(domain_submatrix_type), DIMENSION(:), &
         INTENT(IN)                                      :: domain_r_down, domain_s_inv
      TYPE(domain_map_type), INTENT(IN)                  :: domain_map
      INTEGER, DIMENSION(:), INTENT(IN)                  :: cpu_of_domain
      LOGICAL, INTENT(IN)                                :: assume_t0_q0x, just_started, &
                                                            optimize_theta, normalize_orbitals
      REAL(KIND=dp), INTENT(IN)                          :: envelope_amplitude, eps_filter
      INTEGER, INTENT(IN)                                :: special_case
      INTEGER, DIMENSION(:), INTENT(IN)                  :: nocc_of_domain
      INTEGER, INTENT(IN)                                :: order_lanczos
      REAL(KIND=dp), INTENT(IN)                          :: eps_lanczos
      INTEGER, INTENT(IN)                                :: max_iter_lanczos

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_xalmos_from_main_var', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, unit_nr
      REAL(KIND=dp)                                      :: t_norm
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type)                                   :: m_tmp_no_1, m_tmp_oo_1

      CALL timeset(routineN, handle)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      CALL dbcsr_create(m_tmp_no_1, &
                        template=m_quench_t, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(m_tmp_oo_1, &
                        template=m_siginv, &
                        matrix_type=dbcsr_type_no_symmetry)

      CALL dbcsr_copy(m_tmp_no_1, m_var_in)
      IF (optimize_theta) THEN
         ! check that all MO coefficients of the guess are less
         ! than the maximum allowed amplitude
         CALL dbcsr_norm(m_tmp_no_1, &
                         dbcsr_norm_maxabsnorm, norm_scalar=t_norm)
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, *) "Maximum norm of the initial guess: ", t_norm
            WRITE (unit_nr, *) "Maximum allowed amplitude: ", &
               envelope_amplitude
         ENDIF
         IF (t_norm .GT. envelope_amplitude .AND. just_started) THEN
            CPABORT("Max norm of the initial guess is too large")
         ENDIF
         ! use artanh to tame MOs
         CALL dbcsr_function_of_elements(m_tmp_no_1, &
                                         func=dbcsr_func_tanh, &
                                         a0=0.0_dp, &
                                         a1=1.0_dp/envelope_amplitude)
         CALL dbcsr_scale(m_tmp_no_1, &
                          envelope_amplitude)
      ENDIF
      CALL dbcsr_hadamard_product(m_tmp_no_1, m_quench_t, &
                                  m_t_out)

      ! project out R_0
      IF (assume_t0_q0x) THEN
         IF (special_case .EQ. xalmo_case_fully_deloc) THEN
            CALL dbcsr_multiply("T", "N", 1.0_dp, &
                                m_STsiginv0, &
                                m_t_out, &
                                0.0_dp, m_tmp_oo_1, &
                                filter_eps=eps_filter)
            CALL dbcsr_multiply("N", "N", -1.0_dp, &
                                m_t0, &
                                m_tmp_oo_1, &
                                1.0_dp, m_t_out, &
                                filter_eps=eps_filter)
         ELSE IF (special_case .EQ. xalmo_case_block_diag) THEN
            CPABORT("cannot use projector with block-daigonal ALMOs")
         ELSE
            ! no special case
            CALL apply_domain_operators( &
               matrix_in=m_t_out, &
               matrix_out=m_tmp_no_1, &
               operator1=domain_r_down, &
               operator2=domain_s_inv, &
               dpattern=m_quench_t, &
               map=domain_map, &
               node_of_domain=cpu_of_domain, &
               my_action=1, &
               filter_eps=eps_filter, &
               use_trimmer=.FALSE.)
            CALL dbcsr_copy(m_t_out, &
                            m_tmp_no_1)
         ENDIF ! special case
         CALL dbcsr_add(m_t_out, &
                        m_t0, 1.0_dp, 1.0_dp)
      ENDIF

      IF (normalize_orbitals) THEN
         CALL orthogonalize_mos( &
            ket=m_t_out, &
            overlap=m_tmp_oo_1, &
            metric=m_s, &
            retain_locality=.TRUE., &
            only_normalize=.TRUE., &
            nocc_of_domain=nocc_of_domain(:), &
            eps_filter=eps_filter, &
            order_lanczos=order_lanczos, &
            eps_lanczos=eps_lanczos, &
            max_iter_lanczos=max_iter_lanczos, &
            overlap_sqrti=m_sig_sqrti_ii_out)
      ENDIF

      CALL dbcsr_filter(m_t_out, eps=eps_filter)

      CALL dbcsr_release(m_tmp_no_1)
      CALL dbcsr_release(m_tmp_oo_1)

      CALL timestop(handle)

   END SUBROUTINE compute_xalmos_from_main_var

! *****************************************************************************
!> \brief Compute the preconditioner matrices and invert them if necessary
!> \param domain_prec_out ...
!> \param m_prec_out ...
!> \param m_ks ...
!> \param m_s ...
!> \param m_siginv ...
!> \param m_quench_t ...
!> \param m_FTsiginv ...
!> \param m_siginvTFTsiginv ...
!> \param m_ST ...
!> \param m_STsiginv_out ...
!> \param m_s_vv_out ...
!> \param m_f_vv_out ...
!> \param para_env ...
!> \param blacs_env ...
!> \param nocc_of_domain ...
!> \param domain_s_inv ...
!> \param domain_s_inv_half ...
!> \param domain_s_half ...
!> \param domain_r_down ...
!> \param cpu_of_domain ...
!> \param domain_map ...
!> \param assume_t0_q0x ...
!> \param penalty_occ_vol ...
!> \param penalty_occ_vol_prefactor ...
!> \param eps_filter ...
!> \param neg_thr ...
!> \param spin_factor ...
!> \param special_case ...
!> \param bad_modes_projector_down_out ...
!> \par History
!>       2015.03 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE compute_preconditioner(domain_prec_out, m_prec_out, m_ks, m_s, &
                                     m_siginv, m_quench_t, m_FTsiginv, m_siginvTFTsiginv, m_ST, &
                                     m_STsiginv_out, m_s_vv_out, m_f_vv_out, para_env, &
                                     blacs_env, nocc_of_domain, domain_s_inv, domain_s_inv_half, domain_s_half, &
                                     domain_r_down, cpu_of_domain, &
                                     domain_map, assume_t0_q0x, penalty_occ_vol, penalty_occ_vol_prefactor, &
                                     eps_filter, neg_thr, spin_factor, special_case, bad_modes_projector_down_out)

      TYPE(domain_submatrix_type), DIMENSION(:), &
         INTENT(INOUT)                                   :: domain_prec_out
      TYPE(dbcsr_type), INTENT(INOUT)                    :: m_prec_out, m_ks, m_s
      TYPE(dbcsr_type), INTENT(IN)                       :: m_siginv
      TYPE(dbcsr_type), INTENT(INOUT)                    :: m_quench_t
      TYPE(dbcsr_type), INTENT(IN)                       :: m_FTsiginv, m_siginvTFTsiginv, m_ST
      TYPE(dbcsr_type), INTENT(INOUT), OPTIONAL          :: m_STsiginv_out, m_s_vv_out, m_f_vv_out
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      INTEGER, DIMENSION(:), INTENT(IN)                  :: nocc_of_domain
      TYPE(domain_submatrix_type), DIMENSION(:), &
         INTENT(IN)                                      :: domain_s_inv
      TYPE(domain_submatrix_type), DIMENSION(:), &
         INTENT(IN), OPTIONAL                            :: domain_s_inv_half, domain_s_half
      TYPE(domain_submatrix_type), DIMENSION(:), &
         INTENT(IN)                                      :: domain_r_down
      INTEGER, DIMENSION(:), INTENT(IN)                  :: cpu_of_domain
      TYPE(domain_map_type), INTENT(IN)                  :: domain_map
      LOGICAL, INTENT(IN)                                :: assume_t0_q0x, penalty_occ_vol
      REAL(KIND=dp), INTENT(IN)                          :: penalty_occ_vol_prefactor, eps_filter, &
                                                            neg_thr, spin_factor
      INTEGER, INTENT(IN)                                :: special_case
      TYPE(domain_submatrix_type), DIMENSION(:), &
         INTENT(INOUT), OPTIONAL                         :: bad_modes_projector_down_out

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_preconditioner', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, precond_domain_projector
      TYPE(dbcsr_type)                                   :: m_tmp_nn_1, m_tmp_no_3

      CALL timeset(routineN, handle)

      CALL dbcsr_create(m_tmp_nn_1, &
                        template=m_s, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(m_tmp_no_3, &
                        template=m_quench_t, &
                        matrix_type=dbcsr_type_no_symmetry)

      ! calculate (1-R)F(1-R) and S-SRS
      ! RZK-warning take advantage: some elements will be removed by the quencher
      ! RZK-warning S operations can be performed outside the spin loop to save time
      ! IT IS REQUIRED THAT PRECONDITIONER DOES NOT BREAK THE LOCALITY!!!!
      ! RZK-warning: further optimization is ABSOLUTELY NECESSARY

      ! First S-SRS
      CALL dbcsr_multiply("N", "N", 1.0_dp, &
                          m_ST, &
                          m_siginv, &
                          0.0_dp, m_tmp_no_3, &
                          filter_eps=eps_filter)
      CALL dbcsr_desymmetrize(m_s, m_tmp_nn_1)
      ! return STsiginv if necessary
      IF (PRESENT(m_STsiginv_out)) THEN
         CALL dbcsr_copy(m_STsiginv_out, m_tmp_no_3)
      ENDIF
      IF (special_case .EQ. xalmo_case_fully_deloc) THEN
         ! use S instead of S-SRS
      ELSE
         CALL dbcsr_multiply("N", "T", -1.0_dp, &
                             m_ST, &
                             m_tmp_no_3, &
                             1.0_dp, m_tmp_nn_1, &
                             filter_eps=eps_filter)
      ENDIF
      ! return S_vv = (S or S-SRS) if necessary
      IF (PRESENT(m_s_vv_out)) THEN
         CALL dbcsr_copy(m_s_vv_out, m_tmp_nn_1)
      ENDIF

      ! Second (1-R)F(1-R)
      ! re-create matrix because desymmetrize is buggy -
      ! it will create multiple copies of blocks
      CALL dbcsr_desymmetrize(m_ks, m_prec_out)
      CALL dbcsr_multiply("N", "T", -1.0_dp, &
                          m_FTsiginv, &
                          m_ST, &
                          1.0_dp, m_prec_out, &
                          filter_eps=eps_filter)
      CALL dbcsr_multiply("N", "T", -1.0_dp, &
                          m_ST, &
                          m_FTsiginv, &
                          1.0_dp, m_prec_out, &
                          filter_eps=eps_filter)
      CALL dbcsr_multiply("N", "N", 1.0_dp, &
                          m_ST, &
                          m_siginvTFTsiginv, &
                          0.0_dp, m_tmp_no_3, &
                          filter_eps=eps_filter)
      CALL dbcsr_multiply("N", "T", 1.0_dp, &
                          m_tmp_no_3, &
                          m_ST, &
                          1.0_dp, m_prec_out, &
                          filter_eps=eps_filter)
      ! return F_vv = (I-SR)F(I-RS) if necessary
      IF (PRESENT(m_f_vv_out)) THEN
         CALL dbcsr_copy(m_f_vv_out, m_prec_out)
      ENDIF

#if 0
!penalty_only=.TRUE.
      WRITE (unit_nr, *) "prefactor0:", penalty_occ_vol_prefactor
      !IF (penalty_occ_vol) THEN
      CALL dbcsr_desymmetrize(m_s, &
                              m_prec_out)
      !CALL dbcsr_scale(m_prec_out,-penalty_occ_vol_prefactor)
      !ENDIF
#else
      ! sum up the F_vv and S_vv terms
      CALL dbcsr_add(m_prec_out, m_tmp_nn_1, &
                     1.0_dp, 1.0_dp)
      ! Scale to obtain unit step length
      CALL dbcsr_scale(m_prec_out, 2.0_dp*spin_factor)

      ! add the contribution from the penalty on the occupied volume
      IF (penalty_occ_vol) THEN
         CALL dbcsr_add(m_prec_out, m_tmp_nn_1, &
                        1.0_dp, penalty_occ_vol_prefactor)
      ENDIF
#endif

      CALL dbcsr_copy(m_tmp_nn_1, m_prec_out)

      ! invert using various algorithms
      IF (special_case .EQ. xalmo_case_block_diag) THEN ! non-overlapping diagonal blocks

         CALL pseudo_invert_diagonal_blk( &
            matrix_in=m_tmp_nn_1, &
            matrix_out=m_prec_out, &
            nocc=nocc_of_domain(:) &
            )

      ELSE IF (special_case .EQ. xalmo_case_fully_deloc) THEN ! the entire system is a block

         ! invert using cholesky (works with S matrix, will not work with S-SRS matrix)
         CALL cp_dbcsr_cholesky_decompose(m_prec_out, &
                                          para_env=para_env, &
                                          blacs_env=blacs_env)
         CALL cp_dbcsr_cholesky_invert(m_prec_out, &
                                       para_env=para_env, &
                                       blacs_env=blacs_env, &
                                       upper_to_full=.TRUE.)
         CALL dbcsr_filter(m_prec_out, &
                           eps=eps_filter)

      ELSE

         !!! use a true domain preconditioner with overlapping domains
         IF (assume_t0_q0x) THEN
            precond_domain_projector = -1
         ELSE
            precond_domain_projector = 0
         ENDIF
         !! RZK-warning: use PRESENT to make two nearly-identical calls
         !! this is done because intel compiler does not seem to conform
         !! to the FORTRAN standard for passing through optional arguments
         IF (PRESENT(bad_modes_projector_down_out)) THEN
            CALL construct_domain_preconditioner( &
               matrix_main=m_tmp_nn_1, &
               subm_s_inv=domain_s_inv(:), &
               subm_s_inv_half=domain_s_inv_half(:), &
               subm_s_half=domain_s_half(:), &
               subm_r_down=domain_r_down(:), &
               matrix_trimmer=m_quench_t, &
               dpattern=m_quench_t, &
               map=domain_map, &
               node_of_domain=cpu_of_domain, &
               preconditioner=domain_prec_out(:), &
               use_trimmer=.FALSE., &
               bad_modes_projector_down=bad_modes_projector_down_out(:), &
               eps_zero_eigenvalues=neg_thr, &
               my_action=precond_domain_projector &
               )
         ELSE
            CALL construct_domain_preconditioner( &
               matrix_main=m_tmp_nn_1, &
               subm_s_inv=domain_s_inv(:), &
               subm_r_down=domain_r_down(:), &
               matrix_trimmer=m_quench_t, &
               dpattern=m_quench_t, &
               map=domain_map, &
               node_of_domain=cpu_of_domain, &
               preconditioner=domain_prec_out(:), &
               use_trimmer=.FALSE., &
               !eps_zero_eigenvalues=neg_thr,&
               my_action=precond_domain_projector &
               )
         ENDIF
      ENDIF

      ! invert using cholesky (works with S matrix, will not work with S-SRS matrix)
      !!!CALL cp_dbcsr_cholesky_decompose(prec_vv,&
      !!!        para_env=almo_scf_env%para_env,&
      !!!        blacs_env=almo_scf_env%blacs_env)
      !!!CALL cp_dbcsr_cholesky_invert(prec_vv,&
      !!!        para_env=almo_scf_env%para_env,&
      !!!        blacs_env=almo_scf_env%blacs_env,&
      !!!        upper_to_full=.TRUE.)
      !!!CALL dbcsr_filter(prec_vv,&
      !!!        eps=almo_scf_env%eps_filter)
      !!!

      ! re-create the matrix because desymmetrize is buggy -
      ! it will create multiple copies of blocks
      !!!DESYM!CALL dbcsr_create(prec_vv,&
      !!!DESYM!        template=almo_scf_env%matrix_s(1),&
      !!!DESYM!        matrix_type=dbcsr_type_no_symmetry)
      !!!DESYM!CALL dbcsr_desymmetrize(almo_scf_env%matrix_s(1),&
      !!!DESYM!        prec_vv)
      !CALL dbcsr_multiply("N","N",1.0_dp,&
      !        almo_scf_env%matrix_s(1),&
      !        matrix_t_out(ispin),&
      !        0.0_dp,m_tmp_no_1,&
      !        filter_eps=almo_scf_env%eps_filter)
      !CALL dbcsr_multiply("N","N",1.0_dp,&
      !        m_tmp_no_1,&
      !        almo_scf_env%matrix_sigma_inv(ispin),&
      !        0.0_dp,m_tmp_no_3,&
      !        filter_eps=almo_scf_env%eps_filter)
      !CALL dbcsr_multiply("N","T",-1.0_dp,&
      !        m_tmp_no_3,&
      !        m_tmp_no_1,&
      !        1.0_dp,prec_vv,&
      !        filter_eps=almo_scf_env%eps_filter)
      !CALL dbcsr_add_on_diag(prec_vv,&
      !        prec_sf_mixing_s)

      !CALL dbcsr_create(prec_oo,&
      !        template=almo_scf_env%matrix_sigma(ispin),&
      !        matrix_type=dbcsr_type_no_symmetry)
      !CALL dbcsr_desymmetrize(almo_scf_env%matrix_sigma(ispin),&
      !        matrix_type=dbcsr_type_no_symmetry)
      !CALL dbcsr_desymmetrize(almo_scf_env%matrix_sigma(ispin),&
      !        prec_oo)
      !CALL dbcsr_filter(prec_oo,&
      !        eps=almo_scf_env%eps_filter)

      !! invert using cholesky
      !CALL dbcsr_create(prec_oo_inv,&
      !        template=prec_oo,&
      !        matrix_type=dbcsr_type_no_symmetry)
      !CALL dbcsr_desymmetrize(prec_oo,&
      !        prec_oo_inv)
      !CALL cp_dbcsr_cholesky_decompose(prec_oo_inv,&
      !        para_env=almo_scf_env%para_env,&
      !        blacs_env=almo_scf_env%blacs_env)
      !CALL cp_dbcsr_cholesky_invert(prec_oo_inv,&
      !        para_env=almo_scf_env%para_env,&
      !        blacs_env=almo_scf_env%blacs_env,&
      !        upper_to_full=.TRUE.)

      CALL dbcsr_release(m_tmp_nn_1)
      CALL dbcsr_release(m_tmp_no_3)

      CALL timestop(handle)

   END SUBROUTINE compute_preconditioner

! *****************************************************************************
!> \brief Compute beta for conjugate gradient algorithms
!> \param beta ...
!> \param numer ...
!> \param denom ...
!> \param reset_conjugator ...
!> \param conjugator ...
!> \param grad ...
!> \param prev_grad ...
!> \param step ...
!> \param prev_step ...
!> \param prev_minus_prec_grad ...
!> \par History
!>       2015.04 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE compute_cg_beta(beta, numer, denom, reset_conjugator, conjugator, &
                              grad, prev_grad, step, prev_step, prev_minus_prec_grad)

      REAL(KIND=dp), INTENT(INOUT)                       :: beta
      REAL(KIND=dp), INTENT(INOUT), OPTIONAL             :: numer, denom
      LOGICAL, INTENT(INOUT)                             :: reset_conjugator
      INTEGER, INTENT(IN)                                :: conjugator
      TYPE(dbcsr_type), DIMENSION(:), INTENT(INOUT)      :: grad, prev_grad, step, prev_step
      TYPE(dbcsr_type), DIMENSION(:), INTENT(INOUT), &
         OPTIONAL                                        :: prev_minus_prec_grad

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_cg_beta', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, nsize, unit_nr
      REAL(KIND=dp)                                      :: den, kappa, my_denom, my_numer, &
                                                            my_numer2, my_numer3, num, num2, num3, &
                                                            tau
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type)                                   :: m_tmp_no_1

      CALL timeset(routineN, handle)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      IF (.NOT. PRESENT(prev_minus_prec_grad)) THEN
         IF (conjugator .EQ. cg_fletcher_reeves .OR. &
             conjugator .EQ. cg_polak_ribiere .OR. &
             conjugator .EQ. cg_hager_zhang) THEN
            CPABORT("conjugator needs more input")
         ENDIF
      ENDIF

      ! return num denom so beta can be calculated spin-by-spin
      IF (PRESENT(numer) .OR. PRESENT(denom)) THEN
         IF (conjugator .EQ. cg_hestenes_stiefel .OR. &
             conjugator .EQ. cg_dai_yuan .OR. &
             conjugator .EQ. cg_hager_zhang) THEN
            CPABORT("cannot return numer/denom")
         ENDIF
      ENDIF

      nsize = SIZE(grad)

      my_numer = 0.0_dp
      my_numer2 = 0.0_dp
      my_numer3 = 0.0_dp
      my_denom = 0.0_dp

      DO i = 1, nsize

         CALL dbcsr_create(m_tmp_no_1, &
                           template=grad(i), &
                           matrix_type=dbcsr_type_no_symmetry)

         SELECT CASE (conjugator)
         CASE (cg_hestenes_stiefel)
            CALL dbcsr_copy(m_tmp_no_1, grad(i))
            CALL dbcsr_add(m_tmp_no_1, prev_grad(i), &
                           1.0_dp, -1.0_dp)
            CALL dbcsr_dot(m_tmp_no_1, step(i), num)
            CALL dbcsr_dot(m_tmp_no_1, prev_step(i), den)
         CASE (cg_fletcher_reeves)
            CALL dbcsr_dot(grad(i), step(i), num)
            CALL dbcsr_dot(prev_grad(i), prev_minus_prec_grad(i), den)
         CASE (cg_polak_ribiere)
            CALL dbcsr_dot(prev_grad(i), prev_minus_prec_grad(i), den)
            CALL dbcsr_copy(m_tmp_no_1, grad(i))
            CALL dbcsr_add(m_tmp_no_1, prev_grad(i), 1.0_dp, -1.0_dp)
            CALL dbcsr_dot(m_tmp_no_1, step(i), num)
         CASE (cg_fletcher)
            CALL dbcsr_dot(grad(i), step(i), num)
            CALL dbcsr_dot(prev_grad(i), prev_step(i), den)
         CASE (cg_liu_storey)
            CALL dbcsr_dot(prev_grad(i), prev_step(i), den)
            CALL dbcsr_copy(m_tmp_no_1, grad(i))
            CALL dbcsr_add(m_tmp_no_1, prev_grad(i), 1.0_dp, -1.0_dp)
            CALL dbcsr_dot(m_tmp_no_1, step(i), num)
         CASE (cg_dai_yuan)
            CALL dbcsr_dot(grad(i), step(i), num)
            CALL dbcsr_copy(m_tmp_no_1, grad(i))
            CALL dbcsr_add(m_tmp_no_1, prev_grad(i), 1.0_dp, -1.0_dp)
            CALL dbcsr_dot(m_tmp_no_1, prev_step(i), den)
         CASE (cg_hager_zhang)
            CALL dbcsr_copy(m_tmp_no_1, grad(i))
            CALL dbcsr_add(m_tmp_no_1, prev_grad(i), 1.0_dp, -1.0_dp)
            CALL dbcsr_dot(m_tmp_no_1, prev_step(i), den)
            CALL dbcsr_dot(m_tmp_no_1, prev_minus_prec_grad(i), num)
            CALL dbcsr_dot(m_tmp_no_1, step(i), num2)
            CALL dbcsr_dot(prev_step(i), grad(i), num3)
            my_numer2 = my_numer2 + num2
            my_numer3 = my_numer3 + num3
         CASE (cg_zero)
            num = 0.0_dp
            den = 1.0_dp
         CASE DEFAULT
            CPABORT("illegal conjugator")
         END SELECT
         my_numer = my_numer + num
         my_denom = my_denom + den

         CALL dbcsr_release(m_tmp_no_1)

      ENDDO ! i - nsize

      DO i = 1, nsize

         SELECT CASE (conjugator)
         CASE (cg_hestenes_stiefel, cg_dai_yuan)
            beta = -1.0_dp*my_numer/my_denom
         CASE (cg_fletcher_reeves, cg_polak_ribiere, cg_fletcher, cg_liu_storey)
            beta = my_numer/my_denom
         CASE (cg_hager_zhang)
            kappa = -2.0_dp*my_numer/my_denom
            tau = -1.0_dp*my_numer2/my_denom
            beta = tau - kappa*my_numer3/my_denom
         CASE (cg_zero)
            beta = 0.0_dp
         CASE DEFAULT
            CPABORT("illegal conjugator")
         END SELECT

      ENDDO ! i - nsize

      IF (beta .LT. 0.0_dp) THEN
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, *) " Resetting conjugator because beta is negative: ", beta
         ENDIF
         reset_conjugator = .TRUE.
      ENDIF

      IF (PRESENT(numer)) THEN
         numer = my_numer
      ENDIF
      IF (PRESENT(denom)) THEN
         denom = my_denom
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE compute_cg_beta

! *****************************************************************************
!> \brief computes the step matrix from the gradient and Hessian using
!>         the Newton-Raphson method
!> \param optimizer ...
!> \param m_grad ...
!> \param m_delta ...
!> \param m_s ...
!> \param m_ks ...
!> \param m_siginv ...
!> \param m_quench_t ...
!> \param m_FTsiginv ...
!> \param m_siginvTFTsiginv ...
!> \param m_ST ...
!> \param m_t ...
!> \param m_sig_sqrti_ii ...
!> \param domain_s_inv ...
!> \param domain_r_down ...
!> \param domain_map ...
!> \param cpu_of_domain ...
!> \param nocc_of_domain ...
!> \param para_env ...
!> \param blacs_env ...
!> \param eps_filter ...
!> \param optimize_theta ...
!> \param penalty_occ_vol ...
!> \param normalize_orbitals ...
!> \param penalty_occ_vol_prefactor ...
!> \param penalty_occ_vol_pf2 ...
!> \param special_case ...
!> \par History
!>       2015.04 created [Rustam Z. Khaliullin]
!> \author Rustam Z. Khaliullin
! **************************************************************************************************
   SUBROUTINE newton_grad_to_step(optimizer, m_grad, m_delta, m_s, m_ks, &
                                  m_siginv, m_quench_t, m_FTsiginv, m_siginvTFTsiginv, m_ST, m_t, &
                                  m_sig_sqrti_ii, domain_s_inv, domain_r_down, domain_map, cpu_of_domain, &
                                  nocc_of_domain, para_env, blacs_env, eps_filter, optimize_theta, &
                                  penalty_occ_vol, normalize_orbitals, penalty_occ_vol_prefactor, &
                                  penalty_occ_vol_pf2, special_case)

      TYPE(optimizer_options_type)                       :: optimizer
      TYPE(dbcsr_type), DIMENSION(:), INTENT(IN)         :: m_grad
      TYPE(dbcsr_type), DIMENSION(:), INTENT(INOUT)      :: m_delta, m_s, m_ks, m_siginv, m_quench_t
      TYPE(dbcsr_type), DIMENSION(:), INTENT(IN)         :: m_FTsiginv, m_siginvTFTsiginv, m_ST, &
                                                            m_t, m_sig_sqrti_ii
      TYPE(domain_submatrix_type), DIMENSION(:, :), &
         INTENT(IN)                                      :: domain_s_inv, domain_r_down
      TYPE(domain_map_type), DIMENSION(:), INTENT(IN)    :: domain_map
      INTEGER, DIMENSION(:), INTENT(IN)                  :: cpu_of_domain
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: nocc_of_domain
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter
      LOGICAL, INTENT(IN)                                :: optimize_theta, penalty_occ_vol, &
                                                            normalize_orbitals
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: penalty_occ_vol_prefactor, &
                                                            penalty_occ_vol_pf2
      INTEGER, INTENT(IN)                                :: special_case

      CHARACTER(len=*), PARAMETER :: routineN = 'newton_grad_to_step', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=20)                                  :: iter_type
      INTEGER                                            :: handle, ispin, iteration, max_iter, &
                                                            ndomains, nspins, outer_iteration, &
                                                            outer_max_iter, unit_nr
      LOGICAL :: converged, do_exact_inversion, outer_prepare_to_exit, prepare_to_exit, &
         reset_conjugator, use_preconditioner
      REAL(KIND=dp)                                      :: alpha, beta, denom, denom_ispin, &
                                                            eps_error_target, numer, numer_ispin, &
                                                            residue_norm, spin_factor, t1, t2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: residue_max_norm
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type)                                   :: m_tmp_oo_1, m_tmp_oo_2
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:)        :: m_f_vo, m_f_vv, m_Hstep, m_prec, &
                                                            m_residue, m_residue_prev, m_s_vv, &
                                                            m_step, m_STsiginv, m_zet, m_zet_prev
      TYPE(domain_submatrix_type), ALLOCATABLE, &
         DIMENSION(:, :)                                 :: domain_prec

      CALL timeset(routineN, handle)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      !!! Currently for non-theta only
      IF (optimize_theta) THEN
         CPABORT("theta is NYI")
      ENDIF

      ! set optimizer options
      use_preconditioner = (optimizer%preconditioner .NE. xalmo_prec_zero)
      outer_max_iter = optimizer%max_iter_outer_loop
      max_iter = optimizer%max_iter
      eps_error_target = optimizer%eps_error

      ! set key dimensions
      nspins = SIZE(m_ks)
      ndomains = SIZE(domain_s_inv, 1)

      IF (nspins == 1) THEN
         spin_factor = 2.0_dp
      ELSE
         spin_factor = 1.0_dp
      ENDIF

      ALLOCATE (domain_prec(ndomains, nspins))
      CALL init_submatrices(domain_prec)

      ! allocate matrices
      ALLOCATE (m_residue(nspins))
      ALLOCATE (m_residue_prev(nspins))
      ALLOCATE (m_step(nspins))
      ALLOCATE (m_zet(nspins))
      ALLOCATE (m_zet_prev(nspins))
      ALLOCATE (m_Hstep(nspins))
      ALLOCATE (m_prec(nspins))
      ALLOCATE (m_s_vv(nspins))
      ALLOCATE (m_f_vv(nspins))
      ALLOCATE (m_f_vo(nspins))
      ALLOCATE (m_STsiginv(nspins))

      ALLOCATE (residue_max_norm(nspins))

      ! initiate objects before iterations
      DO ispin = 1, nspins

         ! init matrices
         CALL dbcsr_create(m_residue(ispin), &
                           template=m_quench_t(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(m_residue_prev(ispin), &
                           template=m_quench_t(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(m_step(ispin), &
                           template=m_quench_t(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(m_zet_prev(ispin), &
                           template=m_quench_t(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(m_zet(ispin), &
                           template=m_quench_t(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(m_Hstep(ispin), &
                           template=m_quench_t(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(m_f_vo(ispin), &
                           template=m_quench_t(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(m_STsiginv(ispin), &
                           template=m_quench_t(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(m_f_vv(ispin), &
                           template=m_ks(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(m_s_vv(ispin), &
                           template=m_s(1), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(m_prec(ispin), &
                           template=m_ks(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)

         ! compute the full "gradient" - it is necessary to
         ! evaluate Hessian.X
         CALL dbcsr_copy(m_f_vo(ispin), m_FTsiginv(ispin))
         CALL dbcsr_multiply("N", "N", -1.0_dp, &
                             m_ST(ispin), &
                             m_siginvTFTsiginv(ispin), &
                             1.0_dp, m_f_vo(ispin), &
                             filter_eps=eps_filter)

! RZK-warning
! compute preconditioner even if we do not use it
! this is for debugging because compute_preconditioner includes
! computing F_vv and S_vv necessary for
!       IF ( use_preconditioner ) THEN

! domain_s_inv and domain_r_down are never used with assume_t0_q0x=FALSE
         CALL compute_preconditioner( &
            domain_prec_out=domain_prec(:, ispin), &
            m_prec_out=m_prec(ispin), &
            m_ks=m_ks(ispin), &
            m_s=m_s(1), &
            m_siginv=m_siginv(ispin), &
            m_quench_t=m_quench_t(ispin), &
            m_FTsiginv=m_FTsiginv(ispin), &
            m_siginvTFTsiginv=m_siginvTFTsiginv(ispin), &
            m_ST=m_ST(ispin), &
            m_STsiginv_out=m_STsiginv(ispin), &
            m_s_vv_out=m_s_vv(ispin), &
            m_f_vv_out=m_f_vv(ispin), &
            para_env=para_env, &
            blacs_env=blacs_env, &
            nocc_of_domain=nocc_of_domain(:, ispin), &
            domain_s_inv=domain_s_inv(:, ispin), &
            domain_r_down=domain_r_down(:, ispin), &
            cpu_of_domain=cpu_of_domain(:), &
            domain_map=domain_map(ispin), &
            assume_t0_q0x=.FALSE., &
            penalty_occ_vol=penalty_occ_vol, &
            penalty_occ_vol_prefactor=penalty_occ_vol_prefactor(ispin), &
            eps_filter=eps_filter, &
            neg_thr=0.5_dp, &
            spin_factor=spin_factor, &
            special_case=special_case &
            )

!       ENDIF ! use_preconditioner

         ! initial guess
         CALL dbcsr_copy(m_delta(ispin), m_quench_t(ispin))
         ! in order to use dbcsr_set matrix blocks must exist
         CALL dbcsr_set(m_delta(ispin), 0.0_dp)
         CALL dbcsr_copy(m_residue(ispin), m_grad(ispin))
         CALL dbcsr_scale(m_residue(ispin), -1.0_dp)

         do_exact_inversion = .FALSE.
         IF (do_exact_inversion) THEN

            ! copy grad to m_step temporarily
            ! use m_step as input to the inversion routine
            CALL dbcsr_copy(m_step(ispin), m_grad(ispin))

            ! expensive "exact" inversion of the "nearly-exact" Hessian
            ! hopefully returns Z=-H^(-1).G
            CALL hessian_diag_apply( &
               matrix_grad=m_step(ispin), &
               matrix_step=m_zet(ispin), &
               matrix_S_ao=m_s_vv(ispin), &
               matrix_F_ao=m_f_vv(ispin), &
               !matrix_S_ao=m_s(ispin),&
               !matrix_F_ao=m_ks(ispin),&
               matrix_S_mo=m_siginv(ispin), &
               matrix_F_mo=m_siginvTFTsiginv(ispin), &
               matrix_S_vo=m_STsiginv(ispin), &
               matrix_F_vo=m_f_vo(ispin), &
               quench_t=m_quench_t(ispin), &
               spin_factor=spin_factor, &
               eps_zero=eps_filter*10.0_dp, &
               penalty_occ_vol=penalty_occ_vol, &
               penalty_occ_vol_prefactor=penalty_occ_vol_prefactor(ispin), &
               penalty_occ_vol_pf2=penalty_occ_vol_pf2(ispin), &
               m_s=m_s(1), &
               para_env=para_env, &
               blacs_env=blacs_env &
               )
            ! correct solution by the spin factor
            !CALL dbcsr_scale(m_zet(ispin),1.0_dp/(2.0_dp*spin_factor))

         ELSE ! use PCG to solve H.D=-G

            IF (use_preconditioner) THEN

               IF (special_case .EQ. xalmo_case_block_diag .OR. &
                   special_case .EQ. xalmo_case_fully_deloc) THEN

                  CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                      m_prec(ispin), &
                                      m_residue(ispin), &
                                      0.0_dp, m_zet(ispin), &
                                      filter_eps=eps_filter)

               ELSE

                  CALL apply_domain_operators( &
                     matrix_in=m_residue(ispin), &
                     matrix_out=m_zet(ispin), &
                     operator1=domain_prec(:, ispin), &
                     dpattern=m_quench_t(ispin), &
                     map=domain_map(ispin), &
                     node_of_domain=cpu_of_domain(:), &
                     my_action=0, &
                     filter_eps=eps_filter &
                     !matrix_trimmer=,&
                     !use_trimmer=.FALSE.,&
                     )

               ENDIF ! special_case

            ELSE ! do not use preconditioner

               CALL dbcsr_copy(m_zet(ispin), m_residue(ispin))

            ENDIF ! use_preconditioner

         ENDIF ! do_exact_inversion

         CALL dbcsr_copy(m_step(ispin), m_zet(ispin))

      ENDDO !ispin

      ! start the outer SCF loop
      outer_prepare_to_exit = .FALSE.
      outer_iteration = 0
      residue_norm = 0.0_dp

      DO

         ! start the inner SCF loop
         prepare_to_exit = .FALSE.
         converged = .FALSE.
         iteration = 0
         t1 = m_walltime()

         DO

            ! apply hessian to the step matrix
            CALL apply_hessian( &
               m_x_in=m_step, &
               m_x_out=m_Hstep, &
               m_ks=m_ks, &
               m_s=m_s, &
               m_siginv=m_siginv, &
               m_quench_t=m_quench_t, &
               m_FTsiginv=m_FTsiginv, &
               m_siginvTFTsiginv=m_siginvTFTsiginv, &
               m_ST=m_ST, &
               m_STsiginv=m_STsiginv, &
               m_s_vv=m_s_vv, &
               m_ks_vv=m_f_vv, &
               !m_s_vv=m_s,&
               !m_ks_vv=m_ks,&
               m_g_full=m_f_vo, &
               m_t=m_t, &
               m_sig_sqrti_ii=m_sig_sqrti_ii, &
               penalty_occ_vol=penalty_occ_vol, &
               normalize_orbitals=normalize_orbitals, &
               penalty_occ_vol_prefactor=penalty_occ_vol_prefactor, &
               eps_filter=eps_filter, &
               path_num=hessian_path_reuse &
               )

            ! alpha is computed outside the spin loop
            numer = 0.0_dp
            denom = 0.0_dp
            DO ispin = 1, nspins

               CALL dbcsr_dot(m_residue(ispin), m_zet(ispin), numer_ispin)
               CALL dbcsr_dot(m_step(ispin), m_Hstep(ispin), denom_ispin)

               numer = numer + numer_ispin
               denom = denom + denom_ispin

            ENDDO !ispin

            alpha = numer/denom

            DO ispin = 1, nspins

               ! update the variable
               CALL dbcsr_add(m_delta(ispin), m_step(ispin), 1.0_dp, alpha)
               CALL dbcsr_copy(m_residue_prev(ispin), m_residue(ispin))
               CALL dbcsr_add(m_residue(ispin), m_Hstep(ispin), &
                              1.0_dp, -1.0_dp*alpha)
               CALL dbcsr_norm(m_residue(ispin), dbcsr_norm_maxabsnorm, &
                               norm_scalar=residue_max_norm(ispin))

            ENDDO ! ispin

            ! check convergence and other exit criteria
            residue_norm = MAXVAL(residue_max_norm)
            converged = (residue_norm .LT. eps_error_target)
            IF (converged .OR. (iteration .GE. max_iter)) THEN
               prepare_to_exit = .TRUE.
            ENDIF

            IF (.NOT. prepare_to_exit) THEN

               DO ispin = 1, nspins

                  ! save current z before the update
                  CALL dbcsr_copy(m_zet_prev(ispin), m_zet(ispin))

                  ! compute the new step (apply preconditioner if available)
                  IF (use_preconditioner) THEN

                     !IF (unit_nr>0) THEN
                     !   WRITE(unit_nr,*) "....applying preconditioner...."
                     !ENDIF

                     IF (special_case .EQ. xalmo_case_block_diag .OR. &
                         special_case .EQ. xalmo_case_fully_deloc) THEN

                        CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                            m_prec(ispin), &
                                            m_residue(ispin), &
                                            0.0_dp, m_zet(ispin), &
                                            filter_eps=eps_filter)

                     ELSE

                        CALL apply_domain_operators( &
                           matrix_in=m_residue(ispin), &
                           matrix_out=m_zet(ispin), &
                           operator1=domain_prec(:, ispin), &
                           dpattern=m_quench_t(ispin), &
                           map=domain_map(ispin), &
                           node_of_domain=cpu_of_domain(:), &
                           my_action=0, &
                           filter_eps=eps_filter &
                           !matrix_trimmer=,&
                           !use_trimmer=.FALSE.,&
                           )

                     ENDIF ! special case

                  ELSE

                     CALL dbcsr_copy(m_zet(ispin), m_residue(ispin))

                  ENDIF

               ENDDO !ispin

               ! compute the conjugation coefficient - beta
               CALL compute_cg_beta( &
                  beta=beta, &
                  reset_conjugator=reset_conjugator, &
                  conjugator=cg_fletcher, &
                  grad=m_residue, &
                  prev_grad=m_residue_prev, &
                  step=m_zet, &
                  prev_step=m_zet_prev)

               DO ispin = 1, nspins

                  ! conjugate the step direction
                  CALL dbcsr_add(m_step(ispin), m_zet(ispin), beta, 1.0_dp)

               ENDDO !ispin

            ENDIF ! not.prepare_to_exit

            t2 = m_walltime()
            IF (unit_nr > 0) THEN
               !iter_type=TRIM("ALMO SCF "//iter_type)
               iter_type = TRIM("NR STEP")
               WRITE (unit_nr, '(T6,A9,I6,F14.5,F14.5,F15.10,F9.2)') &
                  iter_type, iteration, &
                  alpha, beta, residue_norm, &
                  t2 - t1
            ENDIF
            t1 = m_walltime()

            iteration = iteration + 1
            IF (prepare_to_exit) EXIT

         ENDDO ! inner loop

         IF (converged .OR. (outer_iteration .GE. outer_max_iter)) THEN
            outer_prepare_to_exit = .TRUE.
         ENDIF

         outer_iteration = outer_iteration + 1
         IF (outer_prepare_to_exit) EXIT

      ENDDO ! outer loop

! is not necessary if penalty_occ_vol_pf2=0.0
#if 0

      IF (penalty_occ_vol) THEN

         DO ispin = 1, nspins

            CALL dbcsr_copy(m_zet(ispin), m_grad(ispin))
            CALL dbcsr_dot(m_delta(ispin), m_zet(ispin), alpha)
            WRITE (unit_nr, *) "trace(grad.delta): ", alpha
            alpha = -1.0_dp/(penalty_occ_vol_pf2(ispin)*alpha - 1.0_dp)
            WRITE (unit_nr, *) "correction alpha: ", alpha
            CALL dbcsr_scale(m_delta(ispin), alpha)

         ENDDO

      ENDIF

#endif

      DO ispin = 1, nspins

         ! check whether the step lies entirely in R or Q
         CALL dbcsr_create(m_tmp_oo_1, &
                           template=m_siginv(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(m_tmp_oo_2, &
                           template=m_siginv(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_multiply("T", "N", 1.0_dp, &
                             m_ST(ispin), &
                             m_delta(ispin), &
                             0.0_dp, m_tmp_oo_1, &
                             filter_eps=eps_filter)
         CALL dbcsr_multiply("N", "N", 1.0_dp, &
                             m_siginv(ispin), &
                             m_tmp_oo_1, &
                             0.0_dp, m_tmp_oo_2, &
                             filter_eps=eps_filter)
         CALL dbcsr_copy(m_zet(ispin), m_quench_t(ispin))
         CALL dbcsr_multiply("N", "N", 1.0_dp, &
                             m_t(ispin), &
                             m_tmp_oo_2, &
                             0.0_dp, m_zet(ispin), &
                             retain_sparsity=.TRUE.)
         CALL dbcsr_norm(m_zet(ispin), dbcsr_norm_maxabsnorm, &
                         norm_scalar=alpha)
         WRITE (unit_nr, "(A50,2F20.10)") "Occupied-space projection of the step", alpha
         CALL dbcsr_add(m_zet(ispin), m_delta(ispin), -1.0_dp, 1.0_dp)
         CALL dbcsr_norm(m_zet(ispin), dbcsr_norm_maxabsnorm, &
                         norm_scalar=alpha)
         WRITE (unit_nr, "(A50,2F20.10)") "Virtual-space projection of the step", alpha
         CALL dbcsr_norm(m_delta(ispin), dbcsr_norm_maxabsnorm, &
                         norm_scalar=alpha)
         WRITE (unit_nr, "(A50,2F20.10)") "Full step", alpha
         CALL dbcsr_release(m_tmp_oo_1)
         CALL dbcsr_release(m_tmp_oo_2)

      ENDDO

      ! clean up
      DO ispin = 1, nspins
         CALL release_submatrices(domain_prec(:, ispin))
         CALL dbcsr_release(m_residue(ispin))
         CALL dbcsr_release(m_residue_prev(ispin))
         CALL dbcsr_release(m_step(ispin))
         CALL dbcsr_release(m_zet(ispin))
         CALL dbcsr_release(m_zet_prev(ispin))
         CALL dbcsr_release(m_Hstep(ispin))
         CALL dbcsr_release(m_f_vo(ispin))
         CALL dbcsr_release(m_f_vv(ispin))
         CALL dbcsr_release(m_s_vv(ispin))
         CALL dbcsr_release(m_prec(ispin))
         CALL dbcsr_release(m_STsiginv(ispin))
      ENDDO !ispin
      DEALLOCATE (domain_prec)
      DEALLOCATE (m_residue)
      DEALLOCATE (m_residue_prev)
      DEALLOCATE (m_step)
      DEALLOCATE (m_zet)
      DEALLOCATE (m_zet_prev)
      DEALLOCATE (m_prec)
      DEALLOCATE (m_Hstep)
      DEALLOCATE (m_s_vv)
      DEALLOCATE (m_f_vv)
      DEALLOCATE (m_f_vo)
      DEALLOCATE (m_STsiginv)
      DEALLOCATE (residue_max_norm)

      IF (.NOT. converged) THEN
         CPABORT("Optimization not converged!")
      ENDIF

      ! check that the step satisfies H.step=-grad

      CALL timestop(handle)

   END SUBROUTINE newton_grad_to_step

! *****************************************************************************
!> \brief Computes Hessian.X
!> \param m_x_in ...
!> \param m_x_out ...
!> \param m_ks ...
!> \param m_s ...
!> \param m_siginv ...
!> \param m_quench_t ...
!> \param m_FTsiginv ...
!> \param m_siginvTFTsiginv ...
!> \param m_ST ...
!> \param m_STsiginv ...
!> \param m_s_vv ...
!> \param m_ks_vv ...
!> \param m_g_full ...
!> \param m_t ...
!> \param m_sig_sqrti_ii ...
!> \param penalty_occ_vol ...
!> \param normalize_orbitals ...
!> \param penalty_occ_vol_prefactor ...
!> \param eps_filter ...
!> \param path_num ...
!> \par History
!>       2015.04 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE apply_hessian(m_x_in, m_x_out, m_ks, m_s, m_siginv, &
                            m_quench_t, m_FTsiginv, m_siginvTFTsiginv, m_ST, m_STsiginv, m_s_vv, &
                            m_ks_vv, m_g_full, m_t, m_sig_sqrti_ii, penalty_occ_vol, &
                            normalize_orbitals, penalty_occ_vol_prefactor, eps_filter, path_num)

      TYPE(dbcsr_type), DIMENSION(:), INTENT(INOUT)      :: m_x_in, m_x_out, m_ks, m_s
      TYPE(dbcsr_type), DIMENSION(:), INTENT(IN)         :: m_siginv
      TYPE(dbcsr_type), DIMENSION(:), INTENT(INOUT)      :: m_quench_t
      TYPE(dbcsr_type), DIMENSION(:), INTENT(IN)         :: m_FTsiginv, m_siginvTFTsiginv, m_ST, &
                                                            m_STsiginv
      TYPE(dbcsr_type), DIMENSION(:), INTENT(INOUT)      :: m_s_vv, m_ks_vv, m_g_full
      TYPE(dbcsr_type), DIMENSION(:), INTENT(IN)         :: m_t, m_sig_sqrti_ii
      LOGICAL, INTENT(IN)                                :: penalty_occ_vol, normalize_orbitals
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: penalty_occ_vol_prefactor
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter
      INTEGER, INTENT(IN)                                :: path_num

      CHARACTER(len=*), PARAMETER :: routineN = 'apply_hessian', routineP = moduleN//':'//routineN

      INTEGER                                            :: dim0, handle, ispin, nspins
      REAL(KIND=dp)                                      :: penalty_prefactor_local, spin_factor
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: tg_diagonal
      TYPE(dbcsr_type)                                   :: m_tmp_no_1, m_tmp_no_2, m_tmp_oo_1, &
                                                            m_tmp_x_in

      CALL timeset(routineN, handle)

      !JHU: test and use for unused debug variables
      IF (penalty_occ_vol) penalty_prefactor_local = 1._dp
      CPASSERT(SIZE(m_STsiginv) >= 0)
      CPASSERT(SIZE(m_siginvTFTsiginv) >= 0)
      CPASSERT(SIZE(m_s) >= 0)
      CPASSERT(SIZE(m_g_full) >= 0)
      CPASSERT(SIZE(m_FTsiginv) >= 0)

      nspins = SIZE(m_ks)

      IF (nspins .EQ. 1) THEN
         spin_factor = 2.0_dp
      ELSE
         spin_factor = 1.0_dp
      ENDIF

      DO ispin = 1, nspins

         penalty_prefactor_local = penalty_occ_vol_prefactor(ispin)/(2.0_dp*spin_factor)

         CALL dbcsr_create(m_tmp_oo_1, &
                           template=m_siginv(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(m_tmp_no_1, &
                           template=m_quench_t(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(m_tmp_no_2, &
                           template=m_quench_t(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(m_tmp_x_in, &
                           template=m_quench_t(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)

         ! transform the input X to take into account the normalization constraint
         IF (normalize_orbitals) THEN

            ! H.D = ( (H.D) - ST.[tr(T).(H.D)]_ii ) . [sig_sqrti]_ii

            ! get [tr(T).HD]_ii
            CALL dbcsr_copy(m_tmp_oo_1, m_sig_sqrti_ii(ispin))
            CALL dbcsr_multiply("T", "N", 1.0_dp, &
                                m_x_in(ispin), &
                                m_ST(ispin), &
                                0.0_dp, m_tmp_oo_1, &
                                retain_sparsity=.TRUE.)
            CALL dbcsr_get_info(m_sig_sqrti_ii(ispin), nfullrows_total=dim0)
            ALLOCATE (tg_diagonal(dim0))
            CALL dbcsr_get_diag(m_tmp_oo_1, tg_diagonal)
            CALL dbcsr_set(m_tmp_oo_1, 0.0_dp)
            CALL dbcsr_set_diag(m_tmp_oo_1, tg_diagonal)
            DEALLOCATE (tg_diagonal)

            CALL dbcsr_copy(m_tmp_no_1, m_x_in(ispin))
            CALL dbcsr_multiply("N", "N", -1.0_dp, &
                                m_t(ispin), &
                                m_tmp_oo_1, &
                                1.0_dp, m_tmp_no_1, &
                                filter_eps=eps_filter)
            CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                m_tmp_no_1, &
                                m_sig_sqrti_ii(ispin), &
                                0.0_dp, m_tmp_x_in, &
                                filter_eps=eps_filter)

         ELSE

            CALL dbcsr_copy(m_tmp_x_in, m_x_in(ispin))

         ENDIF ! normalize_orbitals

         IF (path_num .EQ. hessian_path_reuse) THEN

            ! apply pre-computed F_vv and S_vv to X

#if 0
! RZK-warning: negative sign at penalty_prefactor_local is that
! magical fix for the negative definite problem
! (since penalty_prefactor_local<0 the coeff before S_vv must
! be multiplied by -1 to take the step in the right direction)
!CALL dbcsr_multiply("N","N",-4.0_dp*penalty_prefactor_local,&
!        m_s_vv(ispin),&
!        m_tmp_x_in,&
!        0.0_dp,m_tmp_no_1,&
!        filter_eps=eps_filter)
!CALL dbcsr_copy(m_x_out(ispin),m_quench_t(ispin))
!CALL dbcsr_multiply("N","N",1.0_dp,&
!        m_tmp_no_1,&
!        m_siginv(ispin),&
!        0.0_dp,m_x_out(ispin),&
!        retain_sparsity=.TRUE.)

            CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                m_s(1), &
                                m_tmp_x_in, &
                                0.0_dp, m_tmp_no_1, &
                                filter_eps=eps_filter)
            CALL dbcsr_copy(m_x_out(ispin), m_quench_t(ispin))
            CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                m_tmp_no_1, &
                                m_siginv(ispin), &
                                0.0_dp, m_x_out(ispin), &
                                retain_sparsity=.TRUE.)

!CALL dbcsr_copy(m_x_out(ispin),m_quench_t(ispin))
!CALL dbcsr_multiply("N","N",1.0_dp,&
!        m_s(1),&
!        m_tmp_x_in,&
!        0.0_dp,m_x_out(ispin),&
!        retain_sparsity=.TRUE.)

#else

            ! debugging: only vv matrices, oo matrices are kronecker
            CALL dbcsr_copy(m_x_out(ispin), m_quench_t(ispin))
            CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                m_ks_vv(ispin), &
                                m_tmp_x_in, &
                                0.0_dp, m_x_out(ispin), &
                                retain_sparsity=.TRUE.)

            CALL dbcsr_copy(m_tmp_no_2, m_quench_t(ispin))
            CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                m_s_vv(ispin), &
                                m_tmp_x_in, &
                                0.0_dp, m_tmp_no_2, &
                                retain_sparsity=.TRUE.)
            CALL dbcsr_add(m_x_out(ispin), m_tmp_no_2, &
                           1.0_dp, -4.0_dp*penalty_prefactor_local + 1.0_dp)
#endif

!          ! F_vv.X.S_oo
!          CALL dbcsr_multiply("N","N",1.0_dp,&
!                  m_ks_vv(ispin),&
!                  m_tmp_x_in,&
!                  0.0_dp,m_tmp_no_1,&
!                  filter_eps=eps_filter,&
!                  )
!          CALL dbcsr_copy(m_x_out(ispin),m_quench_t(ispin))
!          CALL dbcsr_multiply("N","N",1.0_dp,&
!                  m_tmp_no_1,&
!                  m_siginv(ispin),&
!                  0.0_dp,m_x_out(ispin),&
!                  retain_sparsity=.TRUE.,&
!                  )
!
!          ! S_vv.X.F_oo
!          CALL dbcsr_multiply("N","N",1.0_dp,&
!                  m_s_vv(ispin),&
!                  m_tmp_x_in,&
!                  0.0_dp,m_tmp_no_1,&
!                  filter_eps=eps_filter,&
!                  )
!          CALL dbcsr_copy(m_tmp_no_2,m_quench_t(ispin))
!          CALL dbcsr_multiply("N","N",1.0_dp,&
!                  m_tmp_no_1,&
!                  m_siginvTFTsiginv(ispin),&
!                  0.0_dp,m_tmp_no_2,&
!                  retain_sparsity=.TRUE.,&
!                  )
!          CALL dbcsr_add(m_x_out(ispin),m_tmp_no_2,&
!               1.0_dp,-1.0_dp)
!! we have to add occ voll penalty here (the Svv termi (i.e. both Svv.D.Soo)
!!  and STsiginv terms)
!
!         ! S_vo.X^t.F_vo
!          CALL dbcsr_multiply("T","N",1.0_dp,&
!                  m_tmp_x_in,&
!                  m_g_full(ispin),&
!                  0.0_dp,m_tmp_oo_1,&
!                  filter_eps=eps_filter,&
!                  )
!          CALL dbcsr_copy(m_tmp_no_2,m_quench_t(ispin))
!          CALL dbcsr_multiply("N","N",1.0_dp,&
!                  m_STsiginv(ispin),&
!                  m_tmp_oo_1,&
!                  0.0_dp,m_tmp_no_2,&
!                  retain_sparsity=.TRUE.,&
!                  )
!          CALL dbcsr_add(m_x_out(ispin),m_tmp_no_2,&
!                  1.0_dp,-1.0_dp)
!
!          ! S_vo.X^t.F_vo
!          CALL dbcsr_multiply("T","N",1.0_dp,&
!                  m_tmp_x_in,&
!                  m_STsiginv(ispin),&
!                  0.0_dp,m_tmp_oo_1,&
!                  filter_eps=eps_filter,&
!                  )
!          CALL dbcsr_copy(m_tmp_no_2,m_quench_t(ispin))
!          CALL dbcsr_multiply("N","N",1.0_dp,&
!                  m_g_full(ispin),&
!                  m_tmp_oo_1,&
!                  0.0_dp,m_tmp_no_2,&
!                  retain_sparsity=.TRUE.,&
!                  )
!          CALL dbcsr_add(m_x_out(ispin),m_tmp_no_2,&
!                  1.0_dp,-1.0_dp)

         ELSE IF (path_num .EQ. hessian_path_assemble) THEN

            ! compute F_vv.X and S_vv.X directly
            ! this path will be advantageous if the number
            ! of PCG iterations is small
            CPABORT("path is NYI")

         ELSE
            CPABORT("illegal path")
         ENDIF ! path

         ! transform the output to take into account the normalization constraint
         IF (normalize_orbitals) THEN

            ! H.D = ( (H.D) - ST.[tr(T).(H.D)]_ii ) . [sig_sqrti]_ii

            ! get [tr(T).HD]_ii
            CALL dbcsr_copy(m_tmp_oo_1, m_sig_sqrti_ii(ispin))
            CALL dbcsr_multiply("T", "N", 1.0_dp, &
                                m_t(ispin), &
                                m_x_out(ispin), &
                                0.0_dp, m_tmp_oo_1, &
                                retain_sparsity=.TRUE.)
            CALL dbcsr_get_info(m_sig_sqrti_ii(ispin), nfullrows_total=dim0)
            ALLOCATE (tg_diagonal(dim0))
            CALL dbcsr_get_diag(m_tmp_oo_1, tg_diagonal)
            CALL dbcsr_set(m_tmp_oo_1, 0.0_dp)
            CALL dbcsr_set_diag(m_tmp_oo_1, tg_diagonal)
            DEALLOCATE (tg_diagonal)

            CALL dbcsr_multiply("N", "N", -1.0_dp, &
                                m_ST(ispin), &
                                m_tmp_oo_1, &
                                1.0_dp, m_x_out(ispin), &
                                retain_sparsity=.TRUE.)
            CALL dbcsr_copy(m_tmp_no_1, m_x_out(ispin))
            CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                m_tmp_no_1, &
                                m_sig_sqrti_ii(ispin), &
                                0.0_dp, m_x_out(ispin), &
                                retain_sparsity=.TRUE.)

         ENDIF ! normalize_orbitals

         CALL dbcsr_scale(m_x_out(ispin), &
                          2.0_dp*spin_factor)

         CALL dbcsr_release(m_tmp_oo_1)
         CALL dbcsr_release(m_tmp_no_1)
         CALL dbcsr_release(m_tmp_no_2)
         CALL dbcsr_release(m_tmp_x_in)

      ENDDO !ispin

      ! there is one more part of the hessian that comes
      ! from T-dependence of the KS matrix
      ! it is neglected here

      CALL timestop(handle)

   END SUBROUTINE apply_hessian

! *****************************************************************************
!> \brief Serial code that constructs an approximate Hessian
!> \param matrix_grad ...
!> \param matrix_step ...
!> \param matrix_S_ao ...
!> \param matrix_F_ao ...
!> \param matrix_S_mo ...
!> \param matrix_F_mo ...
!> \param matrix_S_vo ...
!> \param matrix_F_vo ...
!> \param quench_t ...
!> \param penalty_occ_vol ...
!> \param penalty_occ_vol_prefactor ...
!> \param penalty_occ_vol_pf2 ...
!> \param spin_factor ...
!> \param eps_zero ...
!> \param m_s ...
!> \param para_env ...
!> \param blacs_env ...
!> \par History
!>       2012.02 created [Rustam Z. Khaliullin]
!> \author Rustam Z. Khaliullin
! **************************************************************************************************
   SUBROUTINE hessian_diag_apply(matrix_grad, matrix_step, matrix_S_ao, &
                                 matrix_F_ao, matrix_S_mo, matrix_F_mo, matrix_S_vo, matrix_F_vo, quench_t, &
                                 penalty_occ_vol, penalty_occ_vol_prefactor, penalty_occ_vol_pf2, &
                                 spin_factor, eps_zero, m_s, para_env, blacs_env)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_grad, matrix_step, matrix_S_ao, &
                                                            matrix_F_ao, matrix_S_mo
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_F_mo
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_S_vo, matrix_F_vo, quench_t
      LOGICAL, INTENT(IN)                                :: penalty_occ_vol
      REAL(KIND=dp), INTENT(IN)                          :: penalty_occ_vol_prefactor, &
                                                            penalty_occ_vol_pf2, spin_factor, &
                                                            eps_zero
      TYPE(dbcsr_type), INTENT(IN)                       :: m_s
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'hessian_diag_apply', &
         routineP = moduleN//':'//routineN

      INTEGER :: ao_hori_offset, ao_vert_offset, block_col, block_row, col, H_size, handle, ii, &
         INFO, jj, lev1_hori_offset, lev1_vert_offset, lev2_hori_offset, lev2_vert_offset, LWORK, &
         nblkcols_tot, nblkrows_tot, ncores, orb_i, orb_j, row, unit_nr, zero_neg_eiv
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ao_block_sizes, ao_domain_sizes, &
                                                            mo_block_sizes
      INTEGER, DIMENSION(:), POINTER                     :: ao_blk_sizes, mo_blk_sizes
      LOGICAL                                            :: found, found_col, found_row
      REAL(KIND=dp)                                      :: penalty_prefactor_local, test_error
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues, Grad_vec, Step_vec, tmp, &
                                                            tmpr, work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: F_ao_block, F_mo_block, H, Hinv, &
                                                            S_ao_block, S_mo_block, test, test2
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block_p, p_new_block
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_distribution_type)                      :: main_dist
      TYPE(dbcsr_type)                                   :: matrix_F_ao_sym, matrix_F_mo_sym, &
                                                            matrix_S_ao_sym, matrix_S_mo_sym

      CALL timeset(routineN, handle)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      !JHU use and test for unused debug variables
      CPASSERT(ASSOCIATED(blacs_env))
      CPASSERT(ASSOCIATED(para_env))
      CALL dbcsr_get_info(m_s, row_blk_size=ao_blk_sizes)
      CALL dbcsr_get_info(matrix_S_vo, row_blk_size=ao_blk_sizes)
      CALL dbcsr_get_info(matrix_F_vo, row_blk_size=ao_blk_sizes)

      ! serial code only
      CALL dbcsr_get_info(matrix=matrix_S_ao, distribution=main_dist)
      CALL dbcsr_distribution_get(main_dist, numnodes=ncores)
      IF (ncores .GT. 1) THEN
         CPABORT("serial code only")
      ENDIF

      nblkrows_tot = dbcsr_nblkrows_total(quench_t)
      nblkcols_tot = dbcsr_nblkcols_total(quench_t)
      CPASSERT(nblkrows_tot == nblkcols_tot)
      CALL dbcsr_get_info(quench_t, row_blk_size=ao_blk_sizes)
      CALL dbcsr_get_info(quench_t, col_blk_size=mo_blk_sizes)
      ALLOCATE (mo_block_sizes(nblkcols_tot), ao_block_sizes(nblkcols_tot))
      ALLOCATE (ao_domain_sizes(nblkcols_tot))
      mo_block_sizes(:) = mo_blk_sizes(:)
      ao_block_sizes(:) = ao_blk_sizes(:)
      ao_domain_sizes(:) = 0

      CALL dbcsr_create(matrix_S_ao_sym, &
                        template=matrix_S_ao, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_desymmetrize(matrix_S_ao, matrix_S_ao_sym)
      CALL dbcsr_scale(matrix_S_ao_sym, 2.0_dp*spin_factor)

      CALL dbcsr_create(matrix_F_ao_sym, &
                        template=matrix_F_ao, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_desymmetrize(matrix_F_ao, matrix_F_ao_sym)
      CALL dbcsr_scale(matrix_F_ao_sym, 2.0_dp*spin_factor)

      CALL dbcsr_create(matrix_S_mo_sym, &
                        template=matrix_S_mo, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_desymmetrize(matrix_S_mo, matrix_S_mo_sym)

      CALL dbcsr_create(matrix_F_mo_sym, &
                        template=matrix_F_mo, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_desymmetrize(matrix_F_mo, matrix_F_mo_sym)

      IF (penalty_occ_vol) THEN
         penalty_prefactor_local = penalty_occ_vol_prefactor/(2.0_dp*spin_factor)
      ELSE
         penalty_prefactor_local = 0.0_dp
      ENDIF

      WRITE (unit_nr, *) "penalty_prefactor_local: ", penalty_prefactor_local
      WRITE (unit_nr, *) "penalty_prefactor_2: ", penalty_occ_vol_pf2

      !CALL dbcsr_print(matrix_grad)
      !CALL dbcsr_print(matrix_F_ao_sym)
      !CALL dbcsr_print(matrix_S_ao_sym)
      !CALL dbcsr_print(matrix_F_mo_sym)
      !CALL dbcsr_print(matrix_S_mo_sym)

      ! loop over domains to find the size of the Hessian
      H_size = 0
      DO col = 1, nblkcols_tot

         ! find sizes of AO submatrices
         DO row = 1, nblkrows_tot

            CALL dbcsr_get_block_p(quench_t, &
                                   row, col, block_p, found)
            IF (found) THEN
               ao_domain_sizes(col) = ao_domain_sizes(col) + ao_blk_sizes(row)
            ENDIF

         ENDDO

         H_size = H_size + ao_domain_sizes(col)*mo_block_sizes(col)

      ENDDO

      ALLOCATE (H(H_size, H_size))
      H(:, :) = 0.0_dp

      ! fill the Hessian matrix
      lev1_vert_offset = 0
      ! loop over all pairs of fragments
      DO row = 1, nblkcols_tot

         lev1_hori_offset = 0
         DO col = 1, nblkcols_tot

            ! prepare blocks for the current row-column fragment pair
            ALLOCATE (F_ao_block(ao_domain_sizes(row), ao_domain_sizes(col)))
            ALLOCATE (S_ao_block(ao_domain_sizes(row), ao_domain_sizes(col)))
            ALLOCATE (F_mo_block(mo_block_sizes(row), mo_block_sizes(col)))
            ALLOCATE (S_mo_block(mo_block_sizes(row), mo_block_sizes(col)))

            F_ao_block(:, :) = 0.0_dp
            S_ao_block(:, :) = 0.0_dp
            F_mo_block(:, :) = 0.0_dp
            S_mo_block(:, :) = 0.0_dp

            ! fill AO submatrices
            ! loop over all blocks of the AO dbcsr matrix
            ao_vert_offset = 0
            DO block_row = 1, nblkcols_tot

               CALL dbcsr_get_block_p(quench_t, &
                                      block_row, row, block_p, found_row)
               IF (found_row) THEN

                  ao_hori_offset = 0
                  DO block_col = 1, nblkcols_tot

                     CALL dbcsr_get_block_p(quench_t, &
                                            block_col, col, block_p, found_col)
                     IF (found_col) THEN

                        CALL dbcsr_get_block_p(matrix_F_ao_sym, &
                                               block_row, block_col, block_p, found)
                        IF (found) THEN
                           ! copy the block into the submatrix
                           F_ao_block(ao_vert_offset + 1:ao_vert_offset + ao_block_sizes(block_row), &
                                      ao_hori_offset + 1:ao_hori_offset + ao_block_sizes(block_col)) &
                              = block_p(:, :)
                        ENDIF

                        CALL dbcsr_get_block_p(matrix_S_ao_sym, &
                                               block_row, block_col, block_p, found)
                        IF (found) THEN
                           ! copy the block into the submatrix
                           S_ao_block(ao_vert_offset + 1:ao_vert_offset + ao_block_sizes(block_row), &
                                      ao_hori_offset + 1:ao_hori_offset + ao_block_sizes(block_col)) &
                              = block_p(:, :)
                        ENDIF

                        ao_hori_offset = ao_hori_offset + ao_block_sizes(block_col)

                     ENDIF

                  ENDDO

                  ao_vert_offset = ao_vert_offset + ao_block_sizes(block_row)

               ENDIF

            ENDDO

            ! fill MO submatrices
            CALL dbcsr_get_block_p(matrix_F_mo_sym, row, col, block_p, found)
            IF (found) THEN
               ! copy the block into the submatrix
               F_mo_block(1:mo_block_sizes(row), 1:mo_block_sizes(col)) = block_p(:, :)
            ENDIF
            CALL dbcsr_get_block_p(matrix_S_mo_sym, row, col, block_p, found)
            IF (found) THEN
               ! copy the block into the submatrix
               S_mo_block(1:mo_block_sizes(row), 1:mo_block_sizes(col)) = block_p(:, :)
            ENDIF

            !WRITE(*,*) "F_AO_BLOCK", row, col, ao_domain_sizes(row), ao_domain_sizes(col)
            !DO ii=1,ao_domain_sizes(row)
            !  WRITE(*,'(100F13.9)') F_ao_block(ii,:)
            !ENDDO
            !WRITE(*,*) "S_AO_BLOCK", row, col
            !DO ii=1,ao_domain_sizes(row)
            !  WRITE(*,'(100F13.9)') S_ao_block(ii,:)
            !ENDDO
            !WRITE(*,*) "F_MO_BLOCK", row, col
            !DO ii=1,mo_block_sizes(row)
            !  WRITE(*,'(100F13.9)') F_mo_block(ii,:)
            !ENDDO
            !WRITE(*,*) "S_MO_BLOCK", row, col, mo_block_sizes(row), mo_block_sizes(col)
            !DO ii=1,mo_block_sizes(row)
            !  WRITE(*,'(100F13.9)') S_mo_block(ii,:)
            !ENDDO

            ! construct tensor products for the current row-column fragment pair
            lev2_vert_offset = 0
            DO orb_j = 1, mo_block_sizes(row)

               lev2_hori_offset = 0
               DO orb_i = 1, mo_block_sizes(col)
                  IF (orb_i .EQ. orb_j .AND. row .EQ. col) THEN
                     H(lev1_vert_offset + lev2_vert_offset + 1:lev1_vert_offset + lev2_vert_offset + ao_domain_sizes(row), &
                       lev1_hori_offset + lev2_hori_offset + 1:lev1_hori_offset + lev2_hori_offset + ao_domain_sizes(col)) &
                        != -penalty_prefactor_local*S_ao_block(:,:)
                        = F_ao_block(:, :) + S_ao_block(:, :)
!=S_ao_block(:,:)
!RZK-warning               =F_ao_block(:,:)+( 1.0_dp + penalty_prefactor_local )*S_ao_block(:,:)
!               =S_mo_block(orb_j,orb_i)*F_ao_block(:,:)&
!               -F_mo_block(orb_j,orb_i)*S_ao_block(:,:)&
!               +penalty_prefactor_local*S_mo_block(orb_j,orb_i)*S_ao_block(:,:)
                  ENDIF
                  !WRITE(*,*) row, col, orb_j, orb_i, lev1_vert_offset+lev2_vert_offset+1, ao_domain_sizes(row),&
                  !   lev1_hori_offset+lev2_hori_offset+1, ao_domain_sizes(col), S_mo_block(orb_j,orb_i)

                  lev2_hori_offset = lev2_hori_offset + ao_domain_sizes(col)

               ENDDO

               lev2_vert_offset = lev2_vert_offset + ao_domain_sizes(row)

            ENDDO

            lev1_hori_offset = lev1_hori_offset + ao_domain_sizes(col)*mo_block_sizes(col)

            DEALLOCATE (F_ao_block)
            DEALLOCATE (S_ao_block)
            DEALLOCATE (F_mo_block)
            DEALLOCATE (S_mo_block)

         ENDDO ! col fragment

         lev1_vert_offset = lev1_vert_offset + ao_domain_sizes(row)*mo_block_sizes(row)

      ENDDO ! row fragment

      CALL dbcsr_release(matrix_S_ao_sym)
      CALL dbcsr_release(matrix_F_ao_sym)
      CALL dbcsr_release(matrix_S_mo_sym)
      CALL dbcsr_release(matrix_F_mo_sym)

!!    ! Two more terms of the Hessian: S_vo.D.F_vo and F_vo.D.S_vo
!!    ! It seems that these terms break positive definite property of the Hessian
!!    ALLOCATE(H1(H_size,H_size))
!!    ALLOCATE(H2(H_size,H_size))
!!    H1=0.0_dp
!!    H2=0.0_dp
!!    DO row = 1, nblkcols_tot
!!
!!       lev1_hori_offset=0
!!       DO col = 1, nblkcols_tot
!!
!!          CALL dbcsr_get_block_p(matrix_F_vo,&
!!                  row, col, block_p, found)
!!          CALL dbcsr_get_block_p(matrix_S_vo,&
!!                  row, col, block_p2, found2)
!!
!!          lev1_vert_offset=0
!!          DO block_col = 1, nblkcols_tot
!!
!!             CALL dbcsr_get_block_p(quench_t,&
!!                     row, block_col, p_new_block, found_row)
!!
!!             IF (found_row) THEN
!!
!!                ! determine offset in this short loop
!!                lev2_vert_offset=0
!!                DO block_row=1,row-1
!!                   CALL dbcsr_get_block_p(quench_t,&
!!                           block_row, block_col, p_new_block, found_col)
!!                   IF (found_col) lev2_vert_offset=lev2_vert_offset+ao_block_sizes(block_row)
!!                ENDDO
!!                !!!!!!!! short loop
!!
!!                ! over all electrons of the block
!!                DO orb_i=1, mo_block_sizes(col)
!!
!!                   ! into all possible locations
!!                   DO orb_j=1, mo_block_sizes(block_col)
!!
!!                      ! column is copied several times
!!                      DO copy=1, ao_domain_sizes(col)
!!
!!                         IF (found) THEN
!!
!!                            !WRITE(*,*) row, col, block_col, orb_i, orb_j, copy,&
!!                            ! lev1_vert_offset+(orb_j-1)*ao_domain_sizes(block_col)+lev2_vert_offset+1,&
!!                            ! lev1_hori_offset+(orb_i-1)*ao_domain_sizes(col)+copy
!!
!!                            H1( lev1_vert_offset+(orb_j-1)*ao_domain_sizes(block_col)+lev2_vert_offset+1:&
!!                                lev1_vert_offset+(orb_j-1)*ao_domain_sizes(block_col)+lev2_vert_offset+ao_block_sizes(row),&
!!                                lev1_hori_offset+(orb_i-1)*ao_domain_sizes(col)+copy )&
!!                              =block_p(:,orb_i)
!!
!!                         ENDIF ! found block in the data matrix
!!
!!                         IF (found2) THEN
!!
!!                            H2( lev1_vert_offset+(orb_j-1)*ao_domain_sizes(block_col)+lev2_vert_offset+1:&
!!                                lev1_vert_offset+(orb_j-1)*ao_domain_sizes(block_col)+lev2_vert_offset+ao_block_sizes(row),&
!!                                lev1_hori_offset+(orb_i-1)*ao_domain_sizes(col)+copy )&
!!                              =block_p2(:,orb_i)
!!
!!                         ENDIF ! found block in the data matrix
!!
!!                      ENDDO
!!
!!                   ENDDO
!!
!!                ENDDO
!!
!!                !lev2_vert_offset=lev2_vert_offset+ao_block_sizes(row)
!!
!!             ENDIF ! found block in the quench matrix
!!
!!             lev1_vert_offset=lev1_vert_offset+&
!!                ao_domain_sizes(block_col)*mo_block_sizes(block_col)
!!
!!          ENDDO
!!
!!          lev1_hori_offset=lev1_hori_offset+&
!!             ao_domain_sizes(col)*mo_block_sizes(col)
!!
!!       ENDDO
!!
!!       !lev2_vert_offset=lev2_vert_offset+ao_block_sizes(row)
!!
!!    ENDDO
!!    H1(:,:)=H1(:,:)*2.0_dp*spin_factor
!!    !!!WRITE(*,*) "F_vo"
!!    !!!DO ii=1,H_size
!!    !!! WRITE(*,'(100F13.9)') H1(ii,:)
!!    !!!ENDDO
!!    !!!WRITE(*,*) "S_vo"
!!    !!!DO ii=1,H_size
!!    !!! WRITE(*,'(100F13.9)') H2(ii,:)
!!    !!!ENDDO
!!    !!!!! add terms to the hessian
!!    DO ii=1,H_size
!!       DO jj=1,H_size
!!! add penalty_occ_vol term
!!          H(ii,jj)=H(ii,jj)-H1(ii,jj)*H2(jj,ii)-H1(jj,ii)*H2(ii,jj)
!!       ENDDO
!!    ENDDO
!!    DEALLOCATE(H1)
!!    DEALLOCATE(H2)

!!    ! S_vo.S_vo diagonal component due to determiant constraint
!!    ! use grad vector temporarily
!!    IF (penalty_occ_vol) THEN
!!       ALLOCATE(Grad_vec(H_size))
!!       Grad_vec(:)=0.0_dp
!!       lev1_vert_offset=0
!!       ! loop over all electron blocks
!!       DO col = 1, nblkcols_tot
!!
!!          ! loop over AO-rows of the dbcsr matrix
!!          lev2_vert_offset=0
!!          DO row = 1, nblkrows_tot
!!
!!             CALL dbcsr_get_block_p(quench_t,&
!!                     row, col, block_p, found_row)
!!             IF (found_row) THEN
!!
!!                CALL dbcsr_get_block_p(matrix_S_vo,&
!!                        row, col, block_p, found)
!!                IF (found) THEN
!!                   ! copy the data into the vector, column by column
!!                   DO orb_i=1, mo_block_sizes(col)
!!                      Grad_vec(lev1_vert_offset+ao_domain_sizes(col)*(orb_i-1)+lev2_vert_offset+1:&
!!                               lev1_vert_offset+ao_domain_sizes(col)*(orb_i-1)+lev2_vert_offset+ao_block_sizes(row))&
!!                               =block_p(:,orb_i)
!!                   ENDDO
!!
!!                ENDIF
!!
!!                lev2_vert_offset=lev2_vert_offset+ao_block_sizes(row)
!!
!!             ENDIF
!!
!!          ENDDO
!!
!!          lev1_vert_offset=lev1_vert_offset+ao_domain_sizes(col)*mo_block_sizes(col)
!!
!!       ENDDO ! loop over electron blocks
!!       ! update H now
!!       DO ii=1,H_size
!!          DO jj=1,H_size
!!             H(ii,jj)=H(ii,jj)+penalty_occ_vol_prefactor*&
!!                      penalty_occ_vol_pf2*Grad_vec(ii)*Grad_vec(jj)
!!          ENDDO
!!       ENDDO
!!       DEALLOCATE(Grad_vec)
!!    ENDIF ! penalty_occ_vol

!S-1.G ! invert S using cholesky
!S-1.G CALL dbcsr_create(m_prec_out,&
!S-1.G         template=m_s,&
!S-1.G         matrix_type=dbcsr_type_no_symmetry)
!S-1.G CALL dbcsr_copy(m_prec_out,m_s)
!S-1.G CALL dbcsr_cholesky_decompose(m_prec_out,&
!S-1.G         para_env=para_env,&
!S-1.G         blacs_env=blacs_env)
!S-1.G CALL dbcsr_cholesky_invert(m_prec_out,&
!S-1.G         para_env=para_env,&
!S-1.G         blacs_env=blacs_env,&
!S-1.G         upper_to_full=.TRUE.)
!S-1.G CALL dbcsr_multiply("N","N",1.0_dp,&
!S-1.G         m_prec_out,&
!S-1.G         matrix_grad,&
!S-1.G         0.0_dp,matrix_step,&
!S-1.G         filter_eps=1.0E-10_dp)
!S-1.G !CALL dbcsr_release(m_prec_out)
!S-1.G ALLOCATE(test3(H_size))

      ! convert gradient from the dbcsr matrix to the vector form
      ALLOCATE (Grad_vec(H_size))
      Grad_vec(:) = 0.0_dp
      lev1_vert_offset = 0
      ! loop over all electron blocks
      DO col = 1, nblkcols_tot

         ! loop over AO-rows of the dbcsr matrix
         lev2_vert_offset = 0
         DO row = 1, nblkrows_tot

            CALL dbcsr_get_block_p(quench_t, &
                                   row, col, block_p, found_row)
            IF (found_row) THEN

               CALL dbcsr_get_block_p(matrix_grad, &
                                      row, col, block_p, found)
               IF (found) THEN
                  ! copy the data into the vector, column by column
                  DO orb_i = 1, mo_block_sizes(col)
                     Grad_vec(lev1_vert_offset + ao_domain_sizes(col)*(orb_i - 1) + lev2_vert_offset + 1: &
                              lev1_vert_offset + ao_domain_sizes(col)*(orb_i - 1) + lev2_vert_offset + ao_block_sizes(row)) &
                        = block_p(:, orb_i)
!WRITE(*,*) "GRAD: ", row, col, orb_i, lev1_vert_offset+ao_domain_sizes(col)*(orb_i-1)+lev2_vert_offset+1, ao_block_sizes(row)
                  ENDDO

               ENDIF

!S-1.G CALL dbcsr_get_block_p(matrix_step,&
!S-1.G         row, col, block_p, found)
!S-1.G IF (found) THEN
!S-1.G    ! copy the data into the vector, column by column
!S-1.G    DO orb_i=1, mo_block_sizes(col)
!S-1.G       test3(lev1_vert_offset+ao_domain_sizes(col)*(orb_i-1)+lev2_vert_offset+1:&
!S-1.G                lev1_vert_offset+ao_domain_sizes(col)*(orb_i-1)+lev2_vert_offset+ao_block_sizes(row))&
!S-1.G                =block_p(:,orb_i)
!S-1.G    ENDDO
!S-1.G ENDIF

               lev2_vert_offset = lev2_vert_offset + ao_block_sizes(row)

            ENDIF

         ENDDO

         lev1_vert_offset = lev1_vert_offset + ao_domain_sizes(col)*mo_block_sizes(col)

      ENDDO ! loop over electron blocks

      !WRITE(*,*) "HESSIAN"
      !DO ii=1,H_size
      ! WRITE(*,*) ii
      ! WRITE(*,'(20F14.10)') H(ii,:)
      !ENDDO

      ! invert the Hessian
      INFO = 0
      ALLOCATE (Hinv(H_size, H_size))
      Hinv(:, :) = H(:, :)

      ! before inverting diagonalize
      ALLOCATE (eigenvalues(H_size))
      ! Query the optimal workspace for dsyev
      LWORK = -1
      ALLOCATE (WORK(MAX(1, LWORK)))
      CALL DSYEV('V', 'L', H_size, Hinv, H_size, eigenvalues, WORK, LWORK, INFO)
      LWORK = INT(WORK(1))
      DEALLOCATE (WORK)
      ! Allocate the workspace and solve the eigenproblem
      ALLOCATE (WORK(MAX(1, LWORK)))
      CALL DSYEV('V', 'L', H_size, Hinv, H_size, eigenvalues, WORK, LWORK, INFO)
      IF (INFO .NE. 0) THEN
         WRITE (unit_nr, *) 'DSYEV ERROR MESSAGE: ', INFO
         CPABORT("DSYEV failed")
      END IF
      DEALLOCATE (WORK)

      ! compute grad vector in the basis of Hessian eigenvectors
      ALLOCATE (Step_vec(H_size))
      ! Step_vec contains Grad_vec here
      Step_vec(:) = MATMUL(TRANSPOSE(Hinv), Grad_vec)

      ! compute U.tr(U)-1 = error
      !ALLOCATE(test(H_size,H_size))
      !test(:,:)=MATMUL(TRANSPOSE(Hinv),Hinv)
      !DO ii=1,H_size
      !   test(ii,ii)=test(ii,ii)-1.0_dp
      !ENDDO
      !test_error=0.0_dp
      !DO ii=1,H_size
      !   DO jj=1,H_size
      !      test_error=test_error+test(jj,ii)*test(jj,ii)
      !   ENDDO
      !ENDDO
      !WRITE(*,*) "U.tr(U)-1 error: ", SQRT(test_error)
      !DEALLOCATE(test)

      ! invert eigenvalues and use eigenvectors to compute the Hessian inverse
      ! project out zero-eigenvalue directions
      ALLOCATE (test(H_size, H_size))
      zero_neg_eiv = 0
      DO jj = 1, H_size
         WRITE (unit_nr, "(I10,F20.10,F20.10)") jj, eigenvalues(jj), Step_vec(jj)
         IF (eigenvalues(jj) .GT. eps_zero) THEN
            test(jj, :) = Hinv(:, jj)/eigenvalues(jj)
         ELSE
            test(jj, :) = Hinv(:, jj)*0.0_dp
            zero_neg_eiv = zero_neg_eiv + 1
         ENDIF
      ENDDO
      WRITE (unit_nr, *) 'ZERO OR NEGATIVE EIGENVALUES: ', zero_neg_eiv
      DEALLOCATE (Step_vec)

      ALLOCATE (test2(H_size, H_size))
      test2(:, :) = MATMUL(Hinv, test)
      Hinv(:, :) = test2(:, :)
      DEALLOCATE (test, test2)

      !! shift to kill singularity
      !shift=0.0_dp
      !IF (eigenvalues(1).lt.0.0_dp) THEN
      !   CPABORT("Negative eigenvalue(s)")
      !   shift=abs(eigenvalues(1))
      !   WRITE(*,*) "Lowest eigenvalue: ", eigenvalues(1)
      !ENDIF
      !DO ii=1, H_size
      !   IF (eigenvalues(ii).gt.eps_zero) THEN
      !      shift=shift+min(1.0_dp,eigenvalues(ii))*1.0E-4_dp
      !      EXIT
      !   ENDIF
      !ENDDO
      !WRITE(*,*) "Hessian shift: ", shift
      !DO ii=1, H_size
      !   H(ii,ii)=H(ii,ii)+shift
      !ENDDO
      !! end shift

      DEALLOCATE (eigenvalues)

!!!!    Hinv=H
!!!!    INFO=0
!!!!    CALL DPOTRF('L', H_size, Hinv, H_size, INFO )
!!!!    IF( INFO.NE.0 ) THEN
!!!!       WRITE(*,*) 'DPOTRF ERROR MESSAGE: ', INFO
!!!!       CPABORT("DPOTRF failed")
!!!!    END IF
!!!!    CALL DPOTRI('L', H_size, Hinv, H_size, INFO )
!!!!    IF( INFO.NE.0 ) THEN
!!!!       WRITE(*,*) 'DPOTRI ERROR MESSAGE: ', INFO
!!!!       CPABORT("DPOTRI failed")
!!!!    END IF
!!!!    ! complete the matrix
!!!!    DO ii=1,H_size
!!!!       DO jj=ii+1,H_size
!!!!          Hinv(ii,jj)=Hinv(jj,ii)
!!!!       ENDDO
!!!!    ENDDO

      ! compute the inversion error
      ALLOCATE (test(H_size, H_size))
      test(:, :) = MATMUL(Hinv, H)
      DO ii = 1, H_size
         test(ii, ii) = test(ii, ii) - 1.0_dp
      ENDDO
      test_error = 0.0_dp
      DO ii = 1, H_size
         DO jj = 1, H_size
            test_error = test_error + test(jj, ii)*test(jj, ii)
         ENDDO
      ENDDO
      WRITE (unit_nr, *) "Hessian inversion error: ", SQRT(test_error)
      DEALLOCATE (test)

      ! prepare the output vector
      ALLOCATE (Step_vec(H_size))
      ALLOCATE (tmp(H_size))
      tmp(:) = MATMUL(Hinv, Grad_vec)
      !tmp(:)=MATMUL(Hinv,test3)
      Step_vec(:) = -1.0_dp*tmp(:)

      ALLOCATE (tmpr(H_size))
      tmpr(:) = MATMUL(H, Step_vec)
      tmp(:) = tmpr(:) + Grad_vec(:)
      DEALLOCATE (tmpr)
      WRITE (unit_nr, *) "NEWTOV step error: ", MAXVAL(ABS(tmp))

      DEALLOCATE (tmp)

      DEALLOCATE (H)
      DEALLOCATE (Hinv)
      DEALLOCATE (Grad_vec)

!S-1.G DEALLOCATE(test3)

      ! copy the step from the vector into the dbcsr matrix

      ! re-create the step matrix to remove all blocks
      CALL dbcsr_create(matrix_step, &
                        template=matrix_grad, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_work_create(matrix_step, work_mutable=.TRUE.)

      lev1_vert_offset = 0
      ! loop over all electron blocks
      DO col = 1, nblkcols_tot

         ! loop over AO-rows of the dbcsr matrix
         lev2_vert_offset = 0
         DO row = 1, nblkrows_tot

            CALL dbcsr_get_block_p(quench_t, &
                                   row, col, block_p, found_row)
            IF (found_row) THEN

               NULLIFY (p_new_block)
               CALL dbcsr_reserve_block2d(matrix_step, row, col, p_new_block)
               CPASSERT(ASSOCIATED(p_new_block))
               ! copy the data column by column
               DO orb_i = 1, mo_block_sizes(col)
                  p_new_block(:, orb_i) = &
                     Step_vec(lev1_vert_offset + ao_domain_sizes(col)*(orb_i - 1) + lev2_vert_offset + 1: &
                              lev1_vert_offset + ao_domain_sizes(col)*(orb_i - 1) + lev2_vert_offset + ao_block_sizes(row))
               ENDDO

               lev2_vert_offset = lev2_vert_offset + ao_block_sizes(row)

            ENDIF

         ENDDO

         lev1_vert_offset = lev1_vert_offset + ao_domain_sizes(col)*mo_block_sizes(col)

      ENDDO ! loop over electron blocks

      DEALLOCATE (Step_vec)

      CALL dbcsr_finalize(matrix_step)

!S-1.G CALL dbcsr_create(m_tmp_no_1,&
!S-1.G         template=matrix_step,&
!S-1.G         matrix_type=dbcsr_type_no_symmetry)
!S-1.G CALL dbcsr_multiply("N","N",1.0_dp,&
!S-1.G         m_prec_out,&
!S-1.G         matrix_step,&
!S-1.G         0.0_dp,m_tmp_no_1,&
!S-1.G         filter_eps=1.0E-10_dp,&
!S-1.G         )
!S-1.G CALL dbcsr_copy(matrix_step,m_tmp_no_1)
!S-1.G CALL dbcsr_release(m_tmp_no_1)
!S-1.G CALL dbcsr_release(m_prec_out)

      DEALLOCATE (mo_block_sizes, ao_block_sizes)
      DEALLOCATE (ao_domain_sizes)

      CALL dbcsr_create(matrix_S_ao_sym, &
                        template=quench_t, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_copy(matrix_S_ao_sym, quench_t)
      CALL dbcsr_multiply("N", "N", 1.0_dp, &
                          matrix_F_ao, &
                          matrix_step, &
                          0.0_dp, matrix_S_ao_sym, &
                          retain_sparsity=.TRUE.)
      CALL dbcsr_create(matrix_F_ao_sym, &
                        template=quench_t, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_copy(matrix_F_ao_sym, quench_t)
      CALL dbcsr_multiply("N", "N", 1.0_dp, &
                          matrix_S_ao, &
                          matrix_step, &
                          0.0_dp, matrix_F_ao_sym, &
                          retain_sparsity=.TRUE.)
      CALL dbcsr_add(matrix_S_ao_sym, matrix_F_ao_sym, &
                     1.0_dp, 1.0_dp)
      CALL dbcsr_scale(matrix_S_ao_sym, 2.0_dp*spin_factor)
      CALL dbcsr_add(matrix_S_ao_sym, matrix_grad, &
                     1.0_dp, 1.0_dp)
      CALL dbcsr_norm(matrix_S_ao_sym, dbcsr_norm_maxabsnorm, &
                      norm_scalar=test_error)
      WRITE (unit_nr, *) "NEWTOL step error: ", test_error
      CALL dbcsr_release(matrix_S_ao_sym)
      CALL dbcsr_release(matrix_F_ao_sym)

      CALL timestop(handle)

   END SUBROUTINE hessian_diag_apply

END MODULE almo_scf_optimizer

