!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!
MODULE qs_tddfpt2_methods
   USE bibliography,                    ONLY: Iannuzzi2005,&
                                              cite_reference
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type,&
                                              tddfpt2_control_type
   USE cp_dbcsr_interface,              ONLY: cp_dbcsr_allocate_matrix_set,&
                                              cp_dbcsr_copy,&
                                              cp_dbcsr_deallocate_matrix,&
                                              cp_dbcsr_deallocate_matrix_set,&
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
        cp_fm_create, cp_fm_get_info, cp_fm_get_submatrix, cp_fm_maxabsval, cp_fm_p_type, &
        cp_fm_release, cp_fm_set_all, cp_fm_set_submatrix, cp_fm_to_fm, cp_fm_type
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
   USE hfx_admm_utils,                  ONLY: tddft_hfx_matrix
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
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integrate_potential,          ONLY: integrate_v_rspace
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

! **************************************************************************************************
!> \brief Ground state molecular orbitals
!> \par History
!>   * 06.2016 created [Sergey Chulkov]
! **************************************************************************************************
   TYPE tddfpt_ground_state_mos
      !> occupied MOs stored in a matrix form [nao x nmo_occ]
      TYPE(cp_fm_type), POINTER :: matrix_mos_occ
      !> virtual MOs stored in a matrix form [nao x nmo_virt]
      TYPE(cp_fm_type), POINTER :: matrix_mos_virt
      !> occupied MOs stored as individual vectors [nao x 1]
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:) :: vectors_mo_occ
      !> virtual MOs stored as individual vectors [nao x 1]
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:) :: vectors_mo_virt
      !> occupied orbital energies
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: evals_occ
      !> virtual orbital energies
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: evals_virt
   END TYPE tddfpt_ground_state_mos

! **************************************************************************************************
!> \brief Time-Dependent DFT guess vectors
!> \par History
!>   * 06.2016 created [Sergey Chulkov]
! **************************************************************************************************
   TYPE tddfpt_guess_vectors
      !> indices of occupied MOs to deoccupy
      INTEGER, ALLOCATABLE, DIMENSION(:) :: imos_occ
      !> indices of virtual MOs to occupy
      INTEGER, ALLOCATABLE, DIMENSION(:) :: imos_virt
      !> spin of excited electron (in spin-polarized case)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ispins
      !> transition energies
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: evals
   END TYPE tddfpt_guess_vectors

! **************************************************************************************************
!> \brief Allocatable list of full matrices
!> \par History
!>   * 06.2016 created [Sergey Chulkov]
! **************************************************************************************************
   TYPE tddfpt_fm_vectors
      !> allocatable list of full matrices
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:) :: vects
   END TYPE tddfpt_fm_vectors

! **************************************************************************************************
!> \brief Allocatable 2-D matrix
!> \par History
!>   * 06.2016 created [Sergey Chulkov]
! **************************************************************************************************
   TYPE tddfpt_allocatable_matrix
      !> allocatable 2-D matrix
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :) :: matrix
   END TYPE tddfpt_allocatable_matrix

CONTAINS

! **************************************************************************************************
!> \brief Perform TDDFPT calculation.
!> \param qs_env  Quickstep environment
!> \par History
!>    * 05.2016 created [Sergey Chulkov]
!>    * 06.2016 refactored to be used with Davidson eigensolver [Sergey Chulkov]
!> \note Based on the subroutines apply_op(), iterative_solver(), tddfpt_env_init(), and
!>       tddfpt_env_deallocate().
! **************************************************************************************************
   SUBROUTINE tddfpt(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt', routineP = moduleN//':'//routineN
      INTEGER, PARAMETER                                 :: max_reortho = 16384

      CHARACTER(len=20)                                  :: navail_str, nstates_str
      INTEGER :: energy_unit, handle, ispin, istate, iter, log_unit, max_krylov_vects, mult, nao, &
         niters, nspins, nstates, nstates_conv, nstates_total, nvects_exists, nvects_new
      INTEGER(kind=int_8)                                :: nstates_total_i8
      LOGICAL                                            :: do_hfx, do_hfx_residual, is_nonortho, &
                                                            is_restarted
      REAL(kind=dp)                                      :: C_hf, conv, t1, t2
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: evals_last, evals_prev
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_ks, matrix_s, rho_ij_ao
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: S_C0, S_C0_C0T
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:, :)   :: Aop_evects, Aop_ritz, evects, K_ia_munu, &
                                                            K_ia_ritz, ritz_vects
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: rho_ij_fm
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_ij_r, rho_ij_r2
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_rho_type), POINTER                         :: rho_ij_struct
      TYPE(section_vals_type), POINTER                   :: hfx_section, input, tddfpt_section, &
                                                            xc_section
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control
      TYPE(tddfpt_ground_state_mos), ALLOCATABLE, &
         DIMENSION(:)                                    :: gs_mos
      TYPE(tddfpt_guess_vectors)                         :: guess_vectors
      TYPE(xc_derivative_set_type), POINTER              :: xc_deriv_set
      TYPE(xc_rho_set_type), POINTER                     :: xc_rho_set

      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()

      CALL cite_reference(Iannuzzi2005)

      NULLIFY (blacs_env, dft_control, input, matrix_ks, matrix_s, mos, para_env, pw_env)
      CALL get_qs_env(qs_env, blacs_env=blacs_env, dft_control=dft_control, &
                      input=input, matrix_ks=matrix_ks, matrix_s=matrix_s, mos=mos, &
                      para_env=para_env, pw_env=pw_env)
      tddfpt_control => dft_control%tddfpt2_control

      CPASSERT(dft_control%nimages <= 1)

      IF (tddfpt_control%nstates <= 0) THEN
         CALL integer_to_string(tddfpt_control%nstates, nstates_str)
         CALL cp_warn(__LOCATION__, "TDDFPT calculation was requested for "// &
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
      do_hfx_residual = do_hfx .AND. tddfpt_control%hfx_residual

      nspins = dft_control%nspins

      ALLOCATE (gs_mos(nspins))
      DO ispin = 1, nspins
         CALL tddfpt_init_ground_state_mos(gs_mos=gs_mos(ispin), mo_set=mos(ispin)%mo_set, &
                                           matrix_ks=matrix_ks(ispin)%matrix, matrix_s=matrix_s(1)%matrix, &
                                           cholesky_method=cholesky_restore, blacs_env=blacs_env, para_env=para_env)
      END DO

      ! check we have enough molecular orbitals to construct the requested number of excited states
      nstates_total_i8 = 0
      DO ispin = 1, nspins
         nstates_total_i8 = nstates_total_i8+SIZE(gs_mos(ispin)%evals_occ, kind=int_8)* &
                            SIZE(gs_mos(ispin)%evals_virt, kind=int_8)
      END DO

      IF (INT(tddfpt_control%nstates, kind=int_8) > nstates_total_i8) THEN
         CALL integer_to_string(tddfpt_control%nstates, nstates_str)
         CALL integer_to_string(INT(nstates_total_i8), navail_str)
         CALL cp_warn(__LOCATION__, TRIM(nstates_str)//" excited states were requested, however only "// &
                      TRIM(navail_str)//" excited states can be constructed using the atomic basis set given.")
         tddfpt_control%nstates = INT(nstates_total_i8)
      END IF

      ! adjust the number of trial vectors
      IF (tddfpt_control%added_states < 0 .OR. &
          INT(tddfpt_control%added_states, kind=int_8)+INT(tddfpt_control%nstates, kind=int_8) > nstates_total_i8) THEN

         tddfpt_control%added_states = INT(nstates_total_i8-INT(tddfpt_control%nstates, kind=int_8))
      END IF
      nstates = tddfpt_control%nstates

      ! reset rks_triplets if UKS is in use
      IF (tddfpt_control%rks_triplets .AND. nspins > 1) THEN
         tddfpt_control%rks_triplets = .FALSE.
         CALL cp_warn(__LOCATION__, "Keyword RKS_TRIPLETS has been ignored for spin-polarized calculations")
      END IF

      ! multiplicity of molecular system
      IF (nspins > 1) THEN
         mult = ABS(SIZE(gs_mos(1)%evals_occ)-SIZE(gs_mos(2)%evals_occ))+1
         IF (mult > 2) &
            CALL cp_warn(__LOCATION__, "There is a convergence issue for multiplicity >= 3")
      ELSE
         IF (tddfpt_control%rks_triplets) THEN
            mult = 3
         ELSE
            mult = 1
         END IF
      END IF

      log_unit = cp_print_key_unit_nr(logger, tddfpt_section, "PRINT%GUESS_VECTORS", extension=".tddfptLog")
      CALL tddfpt_init_guess_vectors(guess_vectors=guess_vectors, nstates=nstates, &
                                     nstates_active=nstates+tddfpt_control%added_states, &
                                     gs_mos=gs_mos, is_add_degenerate=tddfpt_control%add_degenerate, &
                                     degenerate_eps=tddfpt_control%degenerate_eps, &
                                     is_rks_triplets=tddfpt_control%rks_triplets, log_unit=log_unit)
      CALL cp_print_key_finished_output(log_unit, logger, tddfpt_section, "PRINT%GUESS_VECTORS")

      nstates_total = SIZE(guess_vectors%evals)
      ! adjust the number of Krylov vectors
      max_krylov_vects = tddfpt_control%nkvs
      IF (max_krylov_vects < nstates) max_krylov_vects = nstates
      IF (max_krylov_vects > nstates_total) max_krylov_vects = nstates_total

      ! ++ S * C_{0, sigma}
      ALLOCATE (S_C0(nspins))
      DO ispin = 1, nspins
         NULLIFY (S_C0(ispin)%matrix)
         CALL cp_fm_create(S_C0(ispin)%matrix, gs_mos(ispin)%matrix_mos_occ%matrix_struct)
         CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, gs_mos(ispin)%matrix_mos_occ, S_C0(ispin)%matrix, &
                                      ncol=SIZE(gs_mos(ispin)%evals_occ), alpha=1.0_dp, beta=0.0_dp)
      END DO

      ! *** allocate storage for trial vectors and matrix-vector products ***
      ALLOCATE (evects(nspins, max_krylov_vects))
      ALLOCATE (Aop_evects(nspins, max_krylov_vects))
      DO istate = 1, max_krylov_vects
         DO ispin = 1, nspins
            NULLIFY (evects(ispin, istate)%matrix, Aop_evects(ispin, istate)%matrix)
         END DO
      END DO

      ALLOCATE (ritz_vects(nspins, nstates))
      ALLOCATE (Aop_ritz(nspins, nstates))
      DO istate = 1, nstates
         DO ispin = 1, nspins
            ! allocate few initial trial vectors, other vectors will be allocated during Krylov's space expansion
            CALL cp_fm_create(evects(ispin, istate)%matrix, gs_mos(ispin)%matrix_mos_occ%matrix_struct)

            NULLIFY (ritz_vects(ispin, istate)%matrix, Aop_ritz(ispin, istate)%matrix)
            CALL cp_fm_create(ritz_vects(ispin, istate)%matrix, gs_mos(ispin)%matrix_mos_occ%matrix_struct)
            CALL cp_fm_create(Aop_ritz(ispin, istate)%matrix, gs_mos(ispin)%matrix_mos_occ%matrix_struct)
         END DO
      END DO

      ! *** initialise xc_rho_set and xc_deriv_set (needed for kernel in adiabatic approximation) ***
      CALL cp_fm_get_info(gs_mos(1)%matrix_mos_occ, nrow_global=nao)

      NULLIFY (auxbas_pw_pool)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, para_env=para_env, context=blacs_env)

      ! *** allocate storage for half-transformed 4-centre integrals ***
      ALLOCATE (K_ia_munu(nspins, max_krylov_vects))
      DO istate = 1, max_krylov_vects
         DO ispin = 1, nspins
            NULLIFY (K_ia_munu(ispin, istate)%matrix)
         END DO
      END DO
      IF (do_hfx) &
         CALL cp_fm_create(K_ia_munu(1, 1)%matrix, fm_struct)

      ALLOCATE (K_ia_ritz(nspins, nstates))
      DO istate = 1, nstates
         DO ispin = 1, nspins
            NULLIFY (K_ia_ritz(ispin, istate)%matrix)
            IF (do_hfx_residual) &
               CALL cp_fm_create(K_ia_ritz(ispin, istate)%matrix, fm_struct)
         END DO
      END DO

      NULLIFY (rho_ij_fm)
      CALL cp_fm_create(rho_ij_fm, fm_struct)

      ! ++ S * P_{0, sigma}
      ALLOCATE (S_C0_C0T(nspins))
      DO ispin = 1, nspins
         NULLIFY (S_C0_C0T(ispin)%matrix)
         CALL cp_fm_create(S_C0_C0T(ispin)%matrix, fm_struct)
      END DO

      ! rho_ij_struct
      NULLIFY (rho_ij_ao)
      CALL cp_dbcsr_allocate_matrix_set(rho_ij_ao, nspins)
      DO ispin = 1, nspins
         CALL cp_dbcsr_init_p(rho_ij_ao(ispin)%matrix)
         CALL cp_dbcsr_copy(rho_ij_ao(ispin)%matrix, matrix_s(1)%matrix)
      END DO

      NULLIFY (rho_ij_struct, rho_ij_r)
      CALL qs_rho_create(rho_ij_struct)
      CALL qs_rho_set(rho_ij_struct, rho_ao=rho_ij_ao)
      CALL qs_rho_rebuild(rho_ij_struct, qs_env, rebuild_ao=.FALSE., rebuild_grids=.TRUE.)
      CALL qs_rho_get(rho_ij_struct, rho_r=rho_ij_r)

      DO ispin = 1, nspins
         CALL cp_gemm('N', 'T', nao, nao, SIZE(gs_mos(ispin)%evals_occ), 1.0_dp, &
                      gs_mos(ispin)%matrix_mos_occ, gs_mos(ispin)%matrix_mos_occ, 0.0_dp, rho_ij_fm)

         CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, rho_ij_fm, S_C0_C0T(ispin)%matrix, &
                                      ncol=nao, alpha=1.0_dp, beta=0.0_dp)

         ! take into account that all MOs are doubly occupied in spin-restricted case
         IF (nspins == 1 .AND. (.NOT. tddfpt_control%rks_triplets)) &
            CALL cp_fm_scale(2.0_dp, rho_ij_fm)

         CALL copy_fm_to_dbcsr(rho_ij_fm, rho_ij_ao(ispin)%matrix, keep_sparsity=.TRUE.)
      END DO
      CALL qs_rho_update_rho(rho_ij_struct, qs_env)

      NULLIFY (xc_rho_set, xc_deriv_set)

      IF (tddfpt_control%rks_triplets) THEN
         ! going to compute triplet states using spin-restricted MOs;
         ! still need spin-beta density in order to compute TDDFT kernel
         ALLOCATE (rho_ij_r2(2))
         rho_ij_r2(1)%pw => rho_ij_r(1)%pw
         rho_ij_r2(2)%pw => rho_ij_r(1)%pw

         CALL xc_prep_2nd_deriv(xc_deriv_set, xc_rho_set, rho_ij_r2, auxbas_pw_pool, xc_section=xc_section)

         DEALLOCATE (rho_ij_r2)
      ELSE
         CALL xc_prep_2nd_deriv(xc_deriv_set, xc_rho_set, rho_ij_r, auxbas_pw_pool, xc_section=xc_section)
      END IF

      CALL qs_rho_release(rho_ij_struct)
      CALL cp_fm_release(rho_ij_fm)
      CALL cp_fm_struct_release(fm_struct)

      ! *** initial Krylov vectors ***
      nvects_exists = 0
      nvects_new = nstates
      DO istate = 1, nvects_new
         ! initial guess
         DO ispin = 1, nspins
            CALL cp_fm_set_all(evects(ispin, istate)%matrix, 0.0_dp)
         END DO

         ispin = guess_vectors%ispins(istate)
         CALL cp_fm_to_fm(gs_mos(ispin)%matrix_mos_virt, evects(ispin, istate)%matrix, &
                          ncol=1, source_start=guess_vectors%imos_virt(istate), &
                          target_start=guess_vectors%imos_occ(istate))
      END DO

      ALLOCATE (evals_last(max_krylov_vects))
      ALLOCATE (evals_prev(nstates))
      evals_prev(1:nstates) = guess_vectors%evals(1:nstates)

      ! *** perform Davidson iterations ***
      log_unit = cp_print_key_unit_nr(logger, tddfpt_section, "PRINT%ITERATION_INFO", extension=".tddfptLog")
      energy_unit = cp_print_key_unit_nr(logger, tddfpt_section, "PRINT%DETAILED_ENERGY", extension=".tddfptLog")

      IF (log_unit > 0) THEN
         WRITE (log_unit, '(/,1X,A,/)') 'TDDFT WAVEFUNCTION OPTIMIZATION'
         WRITE (log_unit, '(5X,A,T15,A,T24,A,T40,A)') "Step", "Time", "Convergence", "Conv. states"
         WRITE (log_unit, '(1X,50("-"))')
      END IF

      CALL cp_add_iter_level(logger%iter_info, "TDDFT_SCF")

      niters = tddfpt_control%niters
      is_restarted = .FALSE.
      t1 = m_walltime()

      DO iter = 1, niters
         CALL cp_iterate(logger%iter_info, last=.FALSE., iter_nr=iter)

         CALL tddfpt_compute_Aop_evects(Aop_evects=Aop_evects(:, nvects_exists+1:nvects_exists+nvects_new), &
                                        K_ia_munu=K_ia_munu(:, nvects_exists+1:nvects_exists+nvects_new), &
                                        evects=evects(:, nvects_exists+1:nvects_exists+nvects_new), &
                                        is_rks_triplets=tddfpt_control%rks_triplets, K_ia_template=K_ia_munu(1, 1)%matrix, &
                                        gs_mos=gs_mos, qs_env=qs_env, matrix_ks=matrix_ks, matrix_s=matrix_s(1)%matrix, &
                                        xc_rho_set=xc_rho_set, xc_deriv_set=xc_deriv_set, xc_section=xc_section)

         CALL tddfpt_compute_ritz_vects(ritz_vects=ritz_vects(:, 1:nstates), Aop_ritz=Aop_ritz(:, 1:nstates), &
                                        evals=evals_last(1:nvects_exists+nvects_new), &
                                        Aop_evects=Aop_evects(:, 1:nvects_exists+nvects_new), &
                                        K_ia_munu=K_ia_munu(:, 1:nvects_exists+nvects_new), &
                                        evects=evects(:, 1:nvects_exists+nvects_new), gs_mos=gs_mos)

         conv = MAXVAL(ABS(evals_last(1:nstates)-evals_prev(1:nstates)))
         IF (conv > tddfpt_control%conv .OR. .NOT. is_restarted) THEN
            ! expand the active space
            nvects_exists = nvects_exists+nvects_new
            IF (nvects_exists < max_krylov_vects .AND. conv > tddfpt_control%conv) THEN

               IF (nvects_exists+nvects_new > max_krylov_vects) &
                  nvects_new = max_krylov_vects-nvects_exists

               IF (do_hfx_residual) &
                  CALL tddfpt_apply_hfx(K_ia_munu=K_ia_ritz(:, 1:nvects_new), evects=ritz_vects(:, 1:nvects_new), &
                                        gs_mos=gs_mos, qs_env=qs_env, matrix_s=matrix_s(1)%matrix)

               CALL tddfpt_compute_residual_vects(residual_vects=evects(:, nvects_exists+1:nvects_exists+nvects_new), &
                                                  evals=evals_last(1:nvects_new), Aop_evects=Aop_ritz(:, 1:nvects_new), &
                                                  K_ia_munu=K_ia_ritz(:, 1:nvects_new), &
                                                  evects=ritz_vects(:, 1:nvects_new), gs_mos=gs_mos, guess_vectors=guess_vectors, &
                                                  matrix_s=matrix_s(1)%matrix)

               DO istate = tddfpt_control%nreortho, 1, -1
                  CALL tddfpt_orthogonalize_psi1_psi0(evects=evects(:, nvects_exists+1:nvects_exists+nvects_new), &
                                                      S_C0_C0T=S_C0_C0T)
                  CALL tddfpt_orthonormalize_psi1_psi1(evects=evects(:, 1:nvects_exists+nvects_new), nvects_new=nvects_new, &
                                                       matrix_s=matrix_s(1)%matrix)
                  is_nonortho = tddfpt_is_nonorthogonal_psi1_psi0(evects=evects(:, nvects_exists+1:nvects_exists+nvects_new), &
                                                                  S_C0=S_C0, max_norm=tddfpt_control%orthogonal_eps)
                  IF (.NOT. is_nonortho) EXIT
               END DO
            ELSE
               is_nonortho = .FALSE.
            END IF

            ! restart Davidson iterations
            IF (nvects_exists >= max_krylov_vects .OR. is_nonortho .OR. &
                (conv <= tddfpt_control%conv .AND. .NOT. is_restarted)) THEN
               DO istate = 1, nstates
                  DO ispin = 1, nspins
                     CALL cp_fm_to_fm(ritz_vects(ispin, istate)%matrix, evects(ispin, istate)%matrix)
                  END DO
               END DO

               DO istate = 1, max_reortho
                  CALL tddfpt_orthogonalize_psi1_psi0(evects=evects(:, 1:nstates), S_C0_C0T=S_C0_C0T)
                  CALL tddfpt_orthonormalize_psi1_psi1(evects=evects(:, 1:nstates), nvects_new=nstates, &
                                                       matrix_s=matrix_s(1)%matrix)
                  is_nonortho = tddfpt_is_nonorthogonal_psi1_psi0(evects=evects(:, 1:nstates), &
                                                                  S_C0=S_C0, max_norm=tddfpt_control%orthogonal_eps)
                  IF (.NOT. is_nonortho) EXIT
               END DO

               nvects_exists = 0
               nvects_new = nstates
               IF (conv > tddfpt_control%conv) &
                  is_restarted = .TRUE.
            END IF
         END IF

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

         IF (log_unit > 0) THEN
            nstates_conv = 0
            DO istate = 1, nstates
               IF (ABS(evals_last(istate)-evals_prev(istate)) <= tddfpt_control%conv) &
                  nstates_conv = nstates_conv+1
            END DO

            WRITE (log_unit, '(1X,I8,T12,F7.1,T24,ES11.4,T42,I8)') iter, t2-t1, conv, nstates_conv

            IF (nvects_exists == 0) &
               WRITE (log_unit, '(1X,10("-"),1X,A,1X,11("-"))') "Restart Davidson iterations"

            CALL m_flush(log_unit)
         END IF

         t1 = t2
         evals_prev(1:nstates) = evals_last(1:nstates)

         IF (conv <= tddfpt_control%conv) THEN
            IF (is_restarted) THEN
               CALL cp_iterate(logger%iter_info, last=.TRUE., iter_nr=iter)
               EXIT
            ELSE
               is_restarted = .TRUE.
            END IF
         END IF
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

      ! *** print summary information ***
      log_unit = cp_logger_get_default_io_unit()
      CALL tddfpt_orthogonalize_psi1_psi0(evects=ritz_vects(:, 1:nstates), S_C0_C0T=S_C0_C0T)
      CALL tddfpt_orthonormalize_psi1_psi1(evects=ritz_vects(:, 1:nstates), nvects_new=nstates, matrix_s=matrix_s(1)%matrix)

      CALL tddfpt_print_summary(log_unit=log_unit, evects=ritz_vects(:, 1:nstates), evals=evals_last(1:nstates), &
                                mult=mult, gs_mos=gs_mos, matrix_s=matrix_s, &
                                min_amplitude=tddfpt_control%min_excitation_amplitude)

      ! -- clean up all useless stuff
      DEALLOCATE (evals_prev)
      DEALLOCATE (evals_last)

      CALL xc_dset_release(xc_deriv_set)
      CALL xc_rho_set_release(xc_rho_set)

      DO ispin = nspins, 1, -1
         CALL cp_fm_release(S_C0_C0T(ispin)%matrix)
         CALL cp_fm_release(S_C0(ispin)%matrix)
      END DO
      DEALLOCATE (S_C0_C0T, S_C0)

      DO istate = SIZE(K_ia_ritz, 2), 1, -1
         DO ispin = nspins, 1, -1
            IF (ASSOCIATED(K_ia_ritz(ispin, istate)%matrix)) &
               CALL cp_fm_release(K_ia_ritz(ispin, istate)%matrix)
         END DO
      END DO
      DEALLOCATE (K_ia_ritz)

      DO istate = SIZE(K_ia_munu, 2), 1, -1
         DO ispin = nspins, 1, -1
            IF (ASSOCIATED(K_ia_munu(ispin, istate)%matrix)) &
               CALL cp_fm_release(K_ia_munu(ispin, istate)%matrix)
         END DO
      END DO
      DEALLOCATE (K_ia_munu)

      DO istate = SIZE(ritz_vects, 2), 1, -1
         DO ispin = nspins, 1, -1
            CALL cp_fm_release(Aop_ritz(ispin, istate)%matrix)
            CALL cp_fm_release(ritz_vects(ispin, istate)%matrix)
         END DO
      END DO
      DEALLOCATE (Aop_ritz)
      DEALLOCATE (ritz_vects)

      DO istate = SIZE(evects, 2), 1, -1
         DO ispin = nspins, 1, -1
            IF (ASSOCIATED(Aop_evects(ispin, istate)%matrix)) &
               CALL cp_fm_release(Aop_evects(ispin, istate)%matrix)
            IF (ASSOCIATED(evects(ispin, istate)%matrix)) &
               CALL cp_fm_release(evects(ispin, istate)%matrix)
         END DO
      END DO
      DEALLOCATE (Aop_evects)
      DEALLOCATE (evects)

      CALL tddfpt_release_guess_vectors(guess_vectors=guess_vectors)

      DO ispin = nspins, 1, -1
         CALL tddfpt_release_ground_state_mos(gs_mos=gs_mos(ispin))
      END DO
      DEALLOCATE (gs_mos)

      CALL timestop(handle)
   END SUBROUTINE tddfpt

! **************************************************************************************************
!> \brief Generate all virtual molecular orbitals for a given spin by diagonalising
!>        the corresponding Kohn-Sham matrix.
!> \param gs_mos          structure to store occupied and virtual molecular orbitals
!>                        (allocated and initialised on exit)
!> \param mo_set          ground state molecular orbitals for a given spin
!> \param matrix_ks       Kohn-Sham matrix for a given spin
!> \param matrix_s        overlap matrix
!> \param cholesky_method Cholesky method to compute the inverse overlap matrix
!> \param blacs_env       BLACS environment for newly allocated distributed matrices
!> \param para_env        parallel environment for newly allocated distributed matrices
!> \par History
!>    * 05.2016 created as tddfpt_lumos() [Sergey Chulkov]
!>    * 06.2016 renamed, altered prototype [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_init_ground_state_mos(gs_mos, mo_set, matrix_ks, matrix_s, cholesky_method, blacs_env, para_env)
      TYPE(tddfpt_ground_state_mos), INTENT(out)         :: gs_mos
      TYPE(mo_set_type), POINTER                         :: mo_set
      TYPE(cp_dbcsr_type), POINTER                       :: matrix_ks, matrix_s
      INTEGER, INTENT(in)                                :: cholesky_method
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_init_ground_state_mos', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: cholesky_method_inout, handle, imo, nao, &
                                                            nelectrons, nmo_occ, nmo_scf, nmo_virt
      REAL(kind=dp)                                      :: maxocc
      REAL(kind=dp), DIMENSION(:), POINTER               :: mo_evals_extended, mo_occ_extended, &
                                                            mo_occ_scf
      TYPE(cp_fm_pool_type), POINTER                     :: ao_ao_fm_pool
      TYPE(cp_fm_struct_type), POINTER                   :: ao_ao_fm_struct, ao_mo_occ_fm_struct, &
                                                            ao_mo_virt_fm_struct, ao_vect_fm_struct
      TYPE(cp_fm_type), POINTER                          :: matrix_ks_fm, mo_coeff_extended, &
                                                            ortho_fm, work_fm
      TYPE(mo_set_type), POINTER                         :: mos_extended

      CALL timeset(routineN, handle)

      NULLIFY (mo_occ_scf)
      CALL get_mo_set(mo_set, nao=nao, nmo=nmo_scf, homo=nmo_occ, maxocc=maxocc, &
                      nelectron=nelectrons, occupation_numbers=mo_occ_scf)

      nmo_virt = nao-nmo_occ

      IF (nmo_virt <= 0) &
         CALL cp_abort(__LOCATION__, &
                       'Unable to generate virtual molecular orbitals. Please use a larger atomic basis set.')

      ! ++ allocate storage space for gs_mos
      NULLIFY (ao_mo_occ_fm_struct, ao_mo_virt_fm_struct, ao_vect_fm_struct)
      CALL cp_fm_struct_create(ao_mo_occ_fm_struct, nrow_global=nao, ncol_global=nmo_occ, &
                               para_env=para_env, context=blacs_env)
      CALL cp_fm_struct_create(ao_mo_virt_fm_struct, nrow_global=nao, ncol_global=nmo_virt, &
                               para_env=para_env, context=blacs_env)
      CALL cp_fm_struct_create(ao_vect_fm_struct, nrow_global=nao, ncol_global=1, &
                               para_env=para_env, context=blacs_env)

      NULLIFY (gs_mos%matrix_mos_occ, gs_mos%matrix_mos_virt)
      CALL cp_fm_create(gs_mos%matrix_mos_occ, ao_mo_occ_fm_struct)
      CALL cp_fm_create(gs_mos%matrix_mos_virt, ao_mo_virt_fm_struct)

      ALLOCATE (gs_mos%vectors_mo_occ(nmo_occ))
      DO imo = 1, nmo_occ
         NULLIFY (gs_mos%vectors_mo_occ(imo)%matrix)
         CALL cp_fm_create(gs_mos%vectors_mo_occ(imo)%matrix, ao_vect_fm_struct)
      END DO

      ALLOCATE (gs_mos%vectors_mo_virt(nmo_virt))
      DO imo = 1, nmo_virt
         NULLIFY (gs_mos%vectors_mo_virt(imo)%matrix)
         CALL cp_fm_create(gs_mos%vectors_mo_virt(imo)%matrix, ao_vect_fm_struct)
      END DO

      ALLOCATE (gs_mos%evals_occ(nmo_occ))
      ALLOCATE (gs_mos%evals_virt(nmo_virt))

      CALL cp_fm_struct_release(ao_vect_fm_struct)
      CALL cp_fm_struct_release(ao_mo_virt_fm_struct)
      CALL cp_fm_struct_release(ao_mo_occ_fm_struct)

      ! ++ set of molecular orbitals
      NULLIFY (ao_ao_fm_struct, ao_ao_fm_pool)
      CALL cp_fm_struct_create(ao_ao_fm_struct, nrow_global=nao, ncol_global=nao, &
                               para_env=para_env, context=blacs_env)
      CALL fm_pool_create(ao_ao_fm_pool, ao_ao_fm_struct)

      NULLIFY (mos_extended, mo_coeff_extended, mo_evals_extended, mo_occ_extended)
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
      DO imo = 1, nmo_scf
         mo_occ_extended(imo) = mo_occ_scf(imo)
      END DO
      mo_occ_extended(nmo_scf+1:) = 0.0_dp

      ! ++ allocate temporary matrices
      NULLIFY (matrix_ks_fm, ortho_fm, work_fm)
      CALL cp_fm_create(matrix_ks_fm, ao_ao_fm_struct)
      CALL cp_fm_create(ortho_fm, ao_ao_fm_struct)
      CALL cp_fm_create(work_fm, ao_ao_fm_struct)

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
      CALL cp_fm_to_fm(mo_coeff_extended, gs_mos%matrix_mos_occ, ncol=nmo_occ, source_start=1, target_start=1)
      DO imo = 1, nmo_occ
         CALL cp_fm_to_fm(mo_coeff_extended, gs_mos%vectors_mo_occ(imo)%matrix, &
                          ncol=1, source_start=imo, target_start=1)
      END DO
      gs_mos%evals_occ(1:nmo_occ) = mo_evals_extended(1:nmo_occ)

      CALL cp_fm_to_fm(mo_coeff_extended, gs_mos%matrix_mos_virt, ncol=nmo_virt, source_start=nmo_occ+1, target_start=1)
      DO imo = 1, nmo_virt
         CALL cp_fm_to_fm(mo_coeff_extended, gs_mos%vectors_mo_virt(imo)%matrix, &
                          ncol=1, source_start=nmo_occ+imo, target_start=1)
      END DO
      gs_mos%evals_virt(1:nmo_virt) = mo_evals_extended(nmo_occ+1:nmo_occ+nmo_virt)

      CALL deallocate_mo_set(mos_extended)
      CALL timestop(handle)
   END SUBROUTINE tddfpt_init_ground_state_mos

! **************************************************************************************************
!> \brief Release molecular orbitals.
!> \param gs_mos          structure that holds occupied and virtual molecular orbitals
!> \par History
!>    * 06.2016 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_release_ground_state_mos(gs_mos)
      TYPE(tddfpt_ground_state_mos), INTENT(inout)       :: gs_mos

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_release_ground_state_mos', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, imo

      CALL timeset(routineN, handle)

      IF (ALLOCATED(gs_mos%evals_virt)) &
         DEALLOCATE (gs_mos%evals_virt)

      IF (ALLOCATED(gs_mos%evals_occ)) &
         DEALLOCATE (gs_mos%evals_occ)

      IF (ALLOCATED(gs_mos%vectors_mo_virt)) THEN
         DO imo = SIZE(gs_mos%vectors_mo_virt), 1, -1
            IF (ASSOCIATED(gs_mos%vectors_mo_virt(imo)%matrix)) &
               CALL cp_fm_release(gs_mos%vectors_mo_virt(imo)%matrix)
         END DO
         DEALLOCATE (gs_mos%vectors_mo_virt)
      END IF

      IF (ALLOCATED(gs_mos%vectors_mo_occ)) THEN
         DO imo = SIZE(gs_mos%vectors_mo_occ), 1, -1
            IF (ASSOCIATED(gs_mos%vectors_mo_occ(imo)%matrix)) &
               CALL cp_fm_release(gs_mos%vectors_mo_occ(imo)%matrix)
         END DO
         DEALLOCATE (gs_mos%vectors_mo_occ)
      END IF

      IF (ASSOCIATED(gs_mos%matrix_mos_virt)) &
         CALL cp_fm_release(gs_mos%matrix_mos_virt)

      IF (ASSOCIATED(gs_mos%matrix_mos_occ)) &
         CALL cp_fm_release(gs_mos%matrix_mos_occ)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_release_ground_state_mos

! **************************************************************************************************
!> \brief Generate the list of guess vectors.
!> \param guess_vectors     guess vectors (initialised on exit)
!> \param nstates           number of excited states to converge
!> \param nstates_active    the number of excited states in the active space
!> \param gs_mos            occupied and virtual molecular orbitals optimised for the ground state
!> \param is_add_degenerate indicates that all degenerate excited states should be automatically
!>                          added to the active space
!> \param degenerate_eps    energy threshold which controls when excited states
!>                          are considered to be degenerate
!> \param is_rks_triplets   indicates that the triplet excited states calculation using
!>                          spin-unpolarised molecular orbitals has been requested
!> \param log_unit          output unit
!> \par History
!>    * 05.2016 created as tddfpt_guess() [Sergey Chulkov]
!>    * 06.2016 renamed, altered prototype, supports spin-polarised density [Sergey Chulkov]
!> \note Based on the subroutine co_initial_guess() which was originally created by
!>       Thomas Chassaing on 06.2003.
! **************************************************************************************************
   SUBROUTINE tddfpt_init_guess_vectors(guess_vectors, nstates, nstates_active, gs_mos, &
                                        is_add_degenerate, degenerate_eps, is_rks_triplets, log_unit)
      TYPE(tddfpt_guess_vectors), INTENT(out)            :: guess_vectors
      INTEGER, INTENT(in)                                :: nstates, nstates_active
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      LOGICAL, INTENT(in)                                :: is_add_degenerate
      REAL(kind=dp), INTENT(in)                          :: degenerate_eps
      LOGICAL, INTENT(in)                                :: is_rks_triplets
      INTEGER, INTENT(in)                                :: log_unit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_init_guess_vectors', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=1)                                   :: spin_occ, spin_virt
      INTEGER :: handle, imo_occ, imo_virt, ind, ispin, istate, jstate, nmo_occ_avail_spin, &
         nmo_occ_selected_spin, nmo_virt_avail_spin, nmo_virt_selected_spin, nspins, &
         nstates_occ_virt_alpha, nstates_selected, nstates_total
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: inds
      INTEGER, DIMENSION(2)                              :: nmo_occ_avail, nmo_occ_selected, &
                                                            nmo_virt_selected
      REAL(kind=dp)                                      :: e_occ
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: e_virt_minus_occ

      CALL timeset(routineN, handle)

      nspins = SIZE(gs_mos)
      CPASSERT(nspins == 1 .OR. nspins == 2)

      DO ispin = 1, nspins
         nmo_occ_avail(ispin) = SIZE(gs_mos(ispin)%evals_occ)
      END DO

      DO ispin = 1, nspins
         nmo_occ_avail_spin = nmo_occ_avail(ispin)
         nmo_virt_avail_spin = SIZE(gs_mos(ispin)%evals_virt)

         DO nmo_occ_selected_spin = nmo_occ_avail_spin-MIN(nstates_active, nmo_occ_avail_spin), 1, -1
            IF (gs_mos(ispin)%evals_occ(nmo_occ_selected_spin+1)- &
                gs_mos(ispin)%evals_occ(nmo_occ_selected_spin) > degenerate_eps) EXIT
         END DO
         nmo_occ_selected_spin = nmo_occ_avail_spin-nmo_occ_selected_spin

         DO nmo_virt_selected_spin = MIN(nstates_active, nmo_virt_avail_spin), nmo_virt_avail_spin-1
            IF (gs_mos(ispin)%evals_virt(nmo_virt_selected_spin+1)- &
                gs_mos(ispin)%evals_virt(nmo_virt_selected_spin) > degenerate_eps) EXIT
         END DO

         nmo_occ_selected(ispin) = nmo_occ_selected_spin
         nmo_virt_selected(ispin) = nmo_virt_selected_spin
      END DO

      ! TO DO: the variable 'nstates_selected' should probably be declared as INTEGER(kind=int_8),
      !        however we need a special version of the subroutine sort() in order to do so
      nstates_selected = DOT_PRODUCT(nmo_occ_selected(1:nspins), nmo_virt_selected(1:nspins))

      ALLOCATE (inds(nstates_selected))
      ALLOCATE (e_virt_minus_occ(nstates_selected))

      istate = 0
      DO ispin = 1, nspins
         nmo_occ_avail_spin = nmo_occ_avail(ispin)
         nmo_occ_selected_spin = nmo_occ_selected(ispin)
         nmo_virt_selected_spin = nmo_virt_selected(ispin)

         DO imo_occ = 1, nmo_occ_selected_spin
            e_occ = gs_mos(ispin)%evals_occ(nmo_occ_avail_spin-imo_occ+1)

            DO imo_virt = 1, nmo_virt_selected_spin
               istate = istate+1
               e_virt_minus_occ(istate) = gs_mos(ispin)%evals_virt(imo_virt)-e_occ
            END DO
         END DO
      END DO

      CPASSERT(istate == nstates_selected)

      CALL sort(e_virt_minus_occ, nstates_selected, inds)

      IF (nspins == 1) THEN
         ispin = 1

         spin_virt = 'a'
         IF (is_rks_triplets) THEN
            spin_occ = 'b'
         ELSE
            spin_occ = 'a'
         END IF
      END IF

      nstates_occ_virt_alpha = nmo_occ_selected(1)*nmo_virt_selected(1)
      IF (log_unit > 0) THEN
         WRITE (log_unit, '(/,22X,A)') "TDDFT initial guess"
         WRITE (log_unit, '(1X,60("-"))')
         WRITE (log_unit, '(5X,A)') "State       Occupied    ->    Virtual        Excitation"
         WRITE (log_unit, '(5X,A)') "number       orbital          orbital        energy (eV)"
         WRITE (log_unit, '(1X,60("-"))')

         DO istate = 1, nstates
            ind = inds(istate)-1
            IF (nspins > 1) THEN
               IF (ind < nstates_occ_virt_alpha) THEN
                  ispin = 1
                  spin_occ = 'a'
               ELSE
                  ispin = 2
                  ind = ind-nstates_occ_virt_alpha
                  spin_occ = 'b'
               END IF
               spin_virt = spin_occ
            END IF

            imo_occ = nmo_occ_avail(ispin)-ind/nmo_virt_selected(ispin)
            imo_virt = nmo_occ_avail(ispin)+MOD(ind, nmo_virt_selected(ispin))+1

            WRITE (log_unit, '(1X,I8,4X,I8,1X,A1,A1,A1,5X,I8,1X,A1,A1,A1,5X,F14.5)') &
               istate, imo_occ, '(', spin_occ, ')', imo_virt, '(', spin_virt, ')', e_virt_minus_occ(istate)*evolt
         END DO
      END IF

      ! detect degenerate states outside the active space
      nstates_total = nstates_active

      DO istate = nstates_total+1, nstates_selected
         IF (e_virt_minus_occ(istate)-e_virt_minus_occ(istate-1) > degenerate_eps) EXIT
      END DO

      IF (istate > nstates_total+1) THEN
         istate = istate-1

         IF (is_add_degenerate) THEN
            nstates_total = istate
         ELSE
            IF (log_unit > 0) THEN
               DO jstate = nstates_total, 2, -1
                  IF (e_virt_minus_occ(jstate)-e_virt_minus_occ(jstate-1) > degenerate_eps) EXIT
               END DO

               IF (jstate <= nstates) &
                  jstate = nstates+1

               IF (jstate <= nstates_total) THEN
                  WRITE (log_unit, '(1X,19("-"),1X,A,1X,20("-"))') "Other active states"

                  DO jstate = jstate, nstates_total
                     ind = inds(jstate)-1
                     IF (nspins > 1) THEN
                        IF (ind < nstates_occ_virt_alpha) THEN
                           ispin = 1
                           spin_occ = 'a'
                        ELSE
                           ispin = 2
                           ind = ind-nstates_occ_virt_alpha
                           spin_occ = 'b'
                        END IF
                        spin_virt = spin_occ
                     END IF

                     imo_occ = nmo_occ_avail(ispin)-ind/nmo_virt_selected(ispin)
                     imo_virt = nmo_occ_avail(ispin)+MOD(ind, nmo_virt_selected(ispin))+1

                     WRITE (log_unit, '(1X,I8,4X,I8,1X,A1,A1,A1,5X,I8,1X,A1,A1,A1,5X,F14.5)') &
                        istate, imo_occ, '(', spin_occ, ')', imo_virt, '(', spin_virt, ')', e_virt_minus_occ(jstate)*evolt
                  END DO
               END IF

               WRITE (log_unit, '(1X,21("-"),1X,A,1X,22("-"))') "Inactive states"

               DO jstate = nstates_total+1, istate
                  ind = inds(jstate)-1
                  IF (nspins > 1) THEN
                     IF (ind < nstates_occ_virt_alpha) THEN
                        ispin = 1
                        spin_occ = 'a'
                     ELSE
                        ispin = 2
                        ind = ind-nstates_occ_virt_alpha
                        spin_occ = 'b'
                     END IF
                     spin_virt = spin_occ
                  END IF

                  imo_occ = nmo_occ_avail(ispin)-ind/nmo_virt_selected(ispin)
                  imo_virt = nmo_occ_avail(ispin)+MOD(ind, nmo_virt_selected(ispin))+1

                  WRITE (log_unit, '(1X,I8,4X,I8,1X,A1,A1,A1,5X,I8,1X,A1,A1,A1,5X,F14.5)') &
                     istate, imo_occ, '(', spin_occ, ')', imo_virt, '(', spin_virt, ')', e_virt_minus_occ(jstate)*evolt
               END DO
            END IF

            CALL cp_warn(__LOCATION__, "It is advised to include all proposed"// &
                         " degenerate excited states into TDDFT active space")
         END IF
      END IF

      IF (log_unit > 0) &
         WRITE (log_unit, '(/,1X,A,T46,I8)') 'Number of active states:', nstates_total

      ALLOCATE (guess_vectors%imos_occ(nstates_total))
      ALLOCATE (guess_vectors%imos_virt(nstates_total))
      ALLOCATE (guess_vectors%ispins(nstates_total))
      ALLOCATE (guess_vectors%evals(nstates_total))

      DO istate = 1, nstates_total
         ind = inds(istate)-1
         IF (ind < nstates_occ_virt_alpha) THEN
            ispin = 1
         ELSE
            ispin = 2
            ind = ind-nstates_occ_virt_alpha
         END IF

         guess_vectors%ispins(istate) = ispin
         guess_vectors%imos_occ(istate) = nmo_occ_avail(ispin)-ind/nmo_virt_selected(ispin)
         guess_vectors%imos_virt(istate) = MOD(ind, nmo_virt_selected(ispin))+1
      END DO
      guess_vectors%evals(1:nstates_total) = e_virt_minus_occ(1:nstates_total)

      DEALLOCATE (e_virt_minus_occ)
      DEALLOCATE (inds)
      CALL timestop(handle)
   END SUBROUTINE tddfpt_init_guess_vectors

! **************************************************************************************************
!> \brief Release initial guess vectors.
!> \param guess_vectors   structure that holds guess vectors
!> \par History
!>    * 06.2016 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_release_guess_vectors(guess_vectors)
      TYPE(tddfpt_guess_vectors), INTENT(inout)          :: guess_vectors

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_release_guess_vectors', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (ALLOCATED(guess_vectors%evals)) &
         DEALLOCATE (guess_vectors%evals)

      IF (ALLOCATED(guess_vectors%ispins)) &
         DEALLOCATE (guess_vectors%ispins)

      IF (ALLOCATED(guess_vectors%imos_virt)) &
         DEALLOCATE (guess_vectors%imos_virt)

      IF (ALLOCATED(guess_vectors%imos_occ)) &
         DEALLOCATE (guess_vectors%imos_occ)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_release_guess_vectors

! **************************************************************************************************
!> \brief Make TDDFPT trial vectors orthogonal to all occupied molecular orbitals.
!> \param evects   trial vectors (modified on exit)
!> \param S_C0_C0T matrix product S * C_0 * C_0^T, where C_0 is the ground state wave-function
!>            for each spin expressed in atomic basis set, and S is the corresponding overlap matrix
!> \par History
!>    * 05.2016 created [Sergey Chulkov]
!> \note  Based on the subroutine p_preortho() which was created by Thomas Chassaing on 09.2002.
!>        Should be useless when the ground states MOs are computed with extremely high accuracy,
!>        as all virtual orbitals are already orthogonal to the occupied ones by design.
!>        It is vitally important, however, when the norm of residual vectors is relatively small
!>        (e.g. less then SCF_EPS). In this case, new krylov vectors seem to be random and should be
!>        orthogonalised even with respect to the occupied MOs.
! **************************************************************************************************
   SUBROUTINE tddfpt_orthogonalize_psi1_psi0(evects, S_C0_C0T)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: evects
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(in)       :: S_C0_C0T

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_orthogonalize_psi1_psi0', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, ivect, nao, nspins, nvects
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nmo_occ
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: C0_C0T_S_C1
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nvects = SIZE(evects, 2)

      CPASSERT(nspins > 0)
      CPASSERT(SIZE(S_C0_C0T) == nspins)

      IF (nvects > 0) THEN
         ALLOCATE (nmo_occ(nspins))
         ALLOCATE (C0_C0T_S_C1(nspins))

         DO ispin = 1, nspins
            NULLIFY (fm_struct, C0_C0T_S_C1(ispin)%matrix)
            CALL cp_fm_get_info(matrix=evects(ispin, 1)%matrix, nrow_global=nao, &
                                ncol_global=nmo_occ(ispin), matrix_struct=fm_struct)
            CALL cp_fm_create(C0_C0T_S_C1(ispin)%matrix, fm_struct)
         END DO

         DO ivect = 1, nvects
            DO ispin = 1, nspins
               ! C0 * C0^T * S * C1 == (S * C0 * C0^T)^T * C1
               CALL cp_gemm('T', 'N', nao, nmo_occ(ispin), nao, 1.0_dp, S_C0_C0T(ispin)%matrix, &
                            evects(ispin, ivect)%matrix, 0.0_dp, C0_C0T_S_C1(ispin)%matrix)

               CALL cp_fm_scale_and_add(1.0_dp, evects(ispin, ivect)%matrix, -1.0_dp, C0_C0T_S_C1(ispin)%matrix)
            END DO
         END DO

         DO ispin = nspins, 1, -1
            CALL cp_fm_release(C0_C0T_S_C1(ispin)%matrix)
         END DO
         DEALLOCATE (C0_C0T_S_C1)
         DEALLOCATE (nmo_occ)
      END IF

      CALL timestop(handle)
   END SUBROUTINE tddfpt_orthogonalize_psi1_psi0

! **************************************************************************************************
!> \brief Check that orthogonalised TDDFPT trial vectors remain orthogonal to
!>        occupied molecular orbitals.
!> \param evects   trial vectors (modified on exit)
!> \param S_C0     matrix product S * C_0, where C_0 is the ground state wave-function
!>                 for each spin in atomic basis set, and S is the corresponding overlap matrix
!> \param max_norm the largest possible overlap between the ground state and
!>                 excited state wave-functions
!> \retval is_nonortho true if trial vectors are non-orthogonal to occupied molecular orbitals
!> \par History
!>    * 07.2016 created [Sergey Chulkov]
! **************************************************************************************************
   FUNCTION tddfpt_is_nonorthogonal_psi1_psi0(evects, S_C0, max_norm) RESULT(is_nonortho)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: evects
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(in)       :: S_C0
      REAL(kind=dp), INTENT(in)                          :: max_norm
      LOGICAL                                            :: is_nonortho

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_is_nonorthogonal_psi1_psi0', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, ivect, nao, nspins, nvects
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nmo_occ
      REAL(kind=dp)                                      :: maxabs_val
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: C0T_S_C1
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nvects = SIZE(evects, 2)

      ALLOCATE (nmo_occ(nspins))
      CALL cp_fm_get_info(matrix=S_C0(1)%matrix, nrow_global=nao, ncol_global=nmo_occ(1), &
                          context=blacs_env, para_env=para_env)

      DO ispin = 2, nspins
         CALL cp_fm_get_info(matrix=S_C0(ispin)%matrix, ncol_global=nmo_occ(ispin))
      END DO

      NULLIFY (fm_struct)
      ALLOCATE (C0T_S_C1(nspins))
      DO ispin = 1, nspins
         NULLIFY (C0T_S_C1(ispin)%matrix)
         CALL cp_fm_struct_create(fm_struct, nrow_global=nmo_occ(ispin), ncol_global=nmo_occ(ispin), &
                                  para_env=para_env, context=blacs_env)
         CALL cp_fm_create(C0T_S_C1(ispin)%matrix, fm_struct)
         CALL cp_fm_struct_release(fm_struct)
      END DO

      is_nonortho = .FALSE.

      loop: DO ivect = 1, nvects
         DO ispin = 1, nspins
            CALL cp_gemm('T', 'N', nmo_occ(ispin), nmo_occ(ispin), nao, 1.0_dp, S_C0(ispin)%matrix, &
                         evects(ispin, ivect)%matrix, 0.0_dp, C0T_S_C1(ispin)%matrix)

            CALL cp_fm_maxabsval(C0T_S_C1(ispin)%matrix, maxabs_val)
            is_nonortho = maxabs_val > max_norm
            IF (is_nonortho) EXIT loop
         END DO
      END DO loop

      DO ispin = SIZE(C0T_S_C1), 1, -1
         CALL cp_fm_release(C0T_S_C1(ispin)%matrix)
      END DO
      DEALLOCATE (C0T_S_C1)

      DEALLOCATE (nmo_occ)
      CALL timestop(handle)
   END FUNCTION tddfpt_is_nonorthogonal_psi1_psi0

! **************************************************************************************************
!> \brief Make new TDDFPT trial vectors orthonormal to all previous trial vectors.
!> \param evects      trial vectors (modified on exit)
!> \param nvects_new  number of new trial vectors to orthogonalise
!> \param matrix_s    overlap matrix
!> \par History
!>    * 05.2016 created [Sergey Chulkov]
!> \note Based on the subroutines reorthogonalize() and normalize() which were originally created
!>             by Thomas Chassaing on 03.2003.
! **************************************************************************************************
   SUBROUTINE tddfpt_orthonormalize_psi1_psi1(evects, nvects_new, matrix_s)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: evects
      INTEGER, INTENT(in)                                :: nvects_new
      TYPE(cp_dbcsr_type), POINTER                       :: matrix_s

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_orthonormalize_psi1_psi1', &
         routineP = moduleN//':'//routineN
      REAL(kind=dp), PARAMETER                           :: threshold = 16*EPSILON(1.0_dp)

      INTEGER                                            :: handle, ispin, ivect, jvect, nspins, &
                                                            nvects_old, nvects_total
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nmo_occ
      REAL(kind=dp)                                      :: norm, norm_spin, weight
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: S_C1
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nvects_total = SIZE(evects, 2)
      nvects_old = nvects_total-nvects_new
      norm_spin = 1.0_dp

      ALLOCATE (nmo_occ(nspins))
      ALLOCATE (S_C1(nspins))
      DO ispin = 1, nspins
         NULLIFY (fm_struct, S_C1(ispin)%matrix)
         CALL cp_fm_get_info(evects(ispin, 1)%matrix, ncol_global=nmo_occ(ispin), matrix_struct=fm_struct)
         CALL cp_fm_create(S_C1(ispin)%matrix, fm_struct)
      END DO

      ! ensure <psi_{i,alpha} | psi_{j,alpha}> == <psi_{i,beta} | psi_{j,beta}> == 0, for i /= j ,
      ! which leads to <psi_i | psi_j> = 0 for i /= j , where psi_i = psi_{i,alpha} + psi_{i,beta}.
      ! The orthogonality condition <psi_i | psi_j> = delta_{i,j} does not imply that
      ! the spin-components psi_{i,alpha} and psi_{i,beta} are normalised by themselves
      ! (<psi_{i,alpha} | psi_{i,alpha}> in [0..1]; <psi_{i,beta} | psi_{i,beta}> in [0..1]),
      ! only their sum does (<psi_{i,alpha} | psi_{i,alpha}> + <psi_{i,beta} | psi_{i,beta}> == 1)
      DO jvect = nvects_old+1, nvects_total
         DO ivect = 1, jvect-1
            DO ispin = 1, nspins
               ! TO DO: do not recompute S * C1 for old trial wectors
               CALL cp_dbcsr_sm_fm_multiply(matrix_s, evects(ispin, ivect)%matrix, S_C1(ispin)%matrix, &
                                            ncol=nmo_occ(ispin), alpha=1.0_dp, beta=0.0_dp)
               IF (nspins > 1) &
                  CALL cp_fm_trace(evects(ispin, ivect)%matrix, S_C1(ispin)%matrix, norm_spin)

               IF (norm_spin > threshold) THEN
                  CALL cp_fm_trace(evects(ispin, jvect)%matrix, S_C1(ispin)%matrix, weight)

                  CALL cp_fm_scale_and_add(1.0_dp, evects(ispin, jvect)%matrix, -weight/norm_spin, evects(ispin, ivect)%matrix)
               END IF
            END DO
         END DO

         norm = 0.0_dp
         DO ispin = 1, nspins
            CALL cp_dbcsr_sm_fm_multiply(matrix_s, evects(ispin, jvect)%matrix, S_C1(ispin)%matrix, &
                                         ncol=nmo_occ(ispin), alpha=1.0_dp, beta=0.0_dp)
            CALL cp_fm_trace(evects(ispin, jvect)%matrix, S_C1(ispin)%matrix, weight)
            norm = norm+weight
         END DO

         norm = 1.0_dp/SQRT(norm)
         DO ispin = 1, nspins
            CALL cp_fm_scale(norm, evects(ispin, ivect)%matrix)
         END DO
      END DO

      DO ispin = nspins, 1, -1
         CALL cp_fm_release(S_C1(ispin)%matrix)
      END DO
      DEALLOCATE (S_C1)
      DEALLOCATE (nmo_occ)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_orthonormalize_psi1_psi1

! **************************************************************************************************
!> \brief Initialise hat{A} C_{1,i} by energy difference using the following expression:
!>        KS_{sigma} * C_{1,i,sigma} - S * C_{1,i,sigma} * diag(evals_occ) .
!> \param Aop_evects  hat{A} C_{1,i} (initialised on exit)
!> \param evects      trial vectors C_{1,i}
!> \param gs_mos      molecular orbitals optimised for the ground state (only occupied orbital
!>                    energies [field %evals_occ] are needed)
!> \param matrix_ks   Kohn-Sham matrix
!> \param matrix_s    overlap matrix
!> \par History
!>    * 05.2016 initialise all matrix elements in one go [Sergey Chulkov]
!> \note Based on the subroutine p_op_l1() which was originally created by
!>       Thomas Chassaing on 08.2002.
! **************************************************************************************************
   SUBROUTINE tddfpt_init_by_energy_diff(Aop_evects, evects, gs_mos, matrix_ks, matrix_s)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: Aop_evects, evects
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      TYPE(cp_dbcsr_p_type), DIMENSION(:), INTENT(in)    :: matrix_ks
      TYPE(cp_dbcsr_type), POINTER                       :: matrix_s

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_init_by_energy_diff', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, ivect, nspins, nvects
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nmo_occ
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: S_C1
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nvects = SIZE(evects, 2)

      ALLOCATE (nmo_occ(nspins))
      ALLOCATE (S_C1(nspins))
      DO ispin = 1, nspins
         NULLIFY (fm_struct, S_C1(ispin)%matrix)
         CALL cp_fm_get_info(evects(ispin, 1)%matrix, ncol_global=nmo_occ(ispin), matrix_struct=fm_struct)
         CALL cp_fm_create(S_C1(ispin)%matrix, fm_struct)
      END DO

      DO ivect = 1, nvects
         DO ispin = 1, nspins
            CALL cp_dbcsr_sm_fm_multiply(matrix_ks(ispin)%matrix, evects(ispin, ivect)%matrix, &
                                         Aop_evects(ispin, ivect)%matrix, ncol=nmo_occ(ispin), &
                                         alpha=1.0_dp, beta=0.0_dp)
            CALL cp_dbcsr_sm_fm_multiply(matrix_s, evects(ispin, ivect)%matrix, S_C1(ispin)%matrix, &
                                         ncol=nmo_occ(ispin), alpha=1.0_dp, beta=0.0_dp)

            CALL cp_fm_column_scale(S_C1(ispin)%matrix, gs_mos(ispin)%evals_occ)
            ! KS * C1 - S * C1 * occupied_orbital_energies
            CALL cp_fm_scale_and_add(1.0_dp, Aop_evects(ispin, ivect)%matrix, -1.0_dp, S_C1(ispin)%matrix)
         END DO
      END DO

      DO ispin = nspins, 1, -1
         CALL cp_fm_release(S_C1(ispin)%matrix)
      END DO
      DEALLOCATE (S_C1)
      DEALLOCATE (nmo_occ)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_init_by_energy_diff

! **************************************************************************************************
!> \brief  Update the hat{A} C_{1,i} vectors by adding coulomb terms
!> \param Aop_evects  hat{A} C_{1,i} (updated on exit)
!> \param evects      trial vectors C_{1,i}
!> \param gs_mos      molecular orbitals optimised for the ground state (only occupied molecular
!>                    orbitals [field %matrix_mos_occ] are needed)
!> \param qs_env      Quickstep environment
!> \param matrix_s    overlap matrix
!> \par History
!>    * 05.2016 compute all coulomb terms in one go [Sergey Chulkov]
!> \note Based on the subroutine kpp1_calc_k_p_p1() which was originally created by
!>       Mohamed Fawzi on 10.2002.
! **************************************************************************************************
   SUBROUTINE tddfpt_apply_coulomb(Aop_evects, evects, gs_mos, qs_env, matrix_s)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: Aop_evects, evects
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_dbcsr_type), POINTER                       :: matrix_s

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_apply_coulomb', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, ivect, nao, nspins, nvects
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nmo_occ
      REAL(kind=dp)                                      :: alpha, pair_energy
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_dbcsr_p_type)                              :: J_ia_munu
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: rho_ia_ao
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: rho_ia_fm
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: v_gspace, v_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_ia_g
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_rho_type), POINTER                         :: rho_ia_struct

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nvects = SIZE(evects, 2)

      ! *** initialize variables ***
      NULLIFY (blacs_env, dft_control, para_env, pw_env)
      CALL get_qs_env(qs_env, dft_control=dft_control, blacs_env=blacs_env, para_env=para_env, pw_env=pw_env)

      CPASSERT(dft_control%nspins == nspins)
      ! a sum J_i{alpha}a{alpha}_munu + J_i{beta}a{beta}_munu can be computed by solving
      ! the Poisson equation for combined density (rho_{ia,alpha} + rho_{ia,beta}) .
      ! Temporary activate spin-restricted mode for acceleration purpose, as we actually
      ! do not need a separate beta-spin component
      dft_control%nspins = 1

      NULLIFY (auxbas_pw_pool, poisson_env)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, poisson_env=poisson_env)

      ALLOCATE (nmo_occ(nspins))
      DO ispin = 1, nspins
         CALL cp_fm_get_info(gs_mos(ispin)%matrix_mos_occ, nrow_global=nao, ncol_global=nmo_occ(ispin))
      END DO

      ! alpha-spin electron density
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
      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, context=blacs_env, para_env=para_env)
      CALL cp_fm_create(rho_ia_fm, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      NULLIFY (v_gspace%pw, v_rspace%pw)
      CALL pw_pool_create_pw(auxbas_pw_pool, v_gspace%pw, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, v_rspace%pw, use_data=REALDATA3D, in_space=REALSPACE)

      IF (nspins > 1) THEN
         alpha = 1.0_dp
      ELSE
         ! spin-restricted case;
         ! alpha == 2 due to singlet state, for triplet states we should not call this subroutine at all
         alpha = 2.0_dp
      END IF

      ! *** compute coulomb terms ***
      DO ivect = 1, nvects
         CALL cp_gemm('N', 'T', nao, nao, nmo_occ(1), 0.5_dp, &
                      gs_mos(1)%matrix_mos_occ, evects(1, ivect)%matrix, 0.0_dp, rho_ia_fm)
         CALL cp_gemm('N', 'T', nao, nao, nmo_occ(1), 0.5_dp, &
                      evects(1, ivect)%matrix, gs_mos(1)%matrix_mos_occ, 1.0_dp, rho_ia_fm)

         DO ispin = 2, nspins
            CALL cp_gemm('N', 'T', nao, nao, nmo_occ(ispin), 0.5_dp, &
                         gs_mos(ispin)%matrix_mos_occ, evects(ispin, ivect)%matrix, 1.0_dp, rho_ia_fm)
            CALL cp_gemm('N', 'T', nao, nao, nmo_occ(ispin), 0.5_dp, &
                         evects(ispin, ivect)%matrix, gs_mos(ispin)%matrix_mos_occ, 1.0_dp, rho_ia_fm)
         END DO

         CALL copy_fm_to_dbcsr(rho_ia_fm, rho_ia_ao(1)%matrix, keep_sparsity=.TRUE.)
         CALL qs_rho_update_rho(rho_ia_struct, qs_env)

         CALL pw_poisson_solve(poisson_env, rho_ia_g(1)%pw, pair_energy, v_gspace%pw)
         CALL pw_transfer(v_gspace%pw, v_rspace%pw)
         CALL pw_scale(v_rspace%pw, v_rspace%pw%pw_grid%dvol)

         CALL cp_dbcsr_set(J_ia_munu%matrix, 0.0_dp)
         CALL integrate_v_rspace(v_rspace=v_rspace, hmat=J_ia_munu, qs_env=qs_env, calculate_forces=.FALSE., gapw=.FALSE.)

         ! (i a || j b) = ( i_alpha a_alpha + i_beta a_beta || j_alpha b_alpha + j_beta b_beta) =
         !                tr (Cj_alpha^T * [J_i{alpha}a{alpha}_munu + J_i{beta}a{beta}_munu] * Cb_alpha) +
         !                tr (Cj_beta^T * [J_i{alpha}a{alpha}_munu + J_i{beta}a{beta}_munu] * Cb_beta)
         DO ispin = 1, nspins
            CALL cp_dbcsr_sm_fm_multiply(J_ia_munu%matrix, gs_mos(ispin)%matrix_mos_occ, Aop_evects(ispin, ivect)%matrix, &
                                         ncol=nmo_occ(ispin), alpha=alpha, beta=1.0_dp)
         END DO
      END DO

      dft_control%nspins = nspins

      ! *** clean up ***
      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_rspace%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_gspace%pw)

      CALL cp_fm_release(rho_ia_fm)

      CALL cp_dbcsr_deallocate_matrix(J_ia_munu%matrix)
      CALL qs_rho_release(rho_ia_struct)

      DEALLOCATE (nmo_occ)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_apply_coulomb

! **************************************************************************************************
!> \brief Update the action matrix by adding kernel terms.
!> \param Aop_evects      hat{A} C_{1,i} (updated on exit)
!> \param evects          TDDFPT trial vectors
!> \param gs_mos          molecular orbitals optimised for the ground state (only occupied molecular
!>                        orbitals [field %matrix_mos_occ] are needed)
!> \param is_rks_triplets indicates that the triplet excited states calculation using
!>                        spin-unpolarised molecular orbitals has been requested
!> \param xc_rho_set      a variable to compute the 2nd derivative of XC-functional
!> \param xc_deriv_set    a variable to compute the 2nd derivative of XC-functional
!> \param qs_env          Quickstep environment
!> \param matrix_s        overlap matrix
!> \param xc_section      XC input section
!> \par History
!>    * 05.2016 compute all kernel terms in one go [Sergey Chulkov]
!> \note Based on the subroutine kpp1_calc_k_p_p1() which was originally created by
!>       Mohamed Fawzi on 10.2002.
! **************************************************************************************************
   SUBROUTINE tddfpt_apply_kernel(Aop_evects, evects, gs_mos, is_rks_triplets, &
                                  xc_rho_set, xc_deriv_set, qs_env, matrix_s, xc_section)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: Aop_evects, evects
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      LOGICAL, INTENT(in)                                :: is_rks_triplets
      TYPE(xc_rho_set_type), POINTER                     :: xc_rho_set
      TYPE(xc_derivative_set_type), POINTER              :: xc_deriv_set
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_dbcsr_type), POINTER                       :: matrix_s
      TYPE(section_vals_type), POINTER                   :: xc_section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_apply_kernel', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, ivect, nao, nspins, &
                                                            nvects, xc_deriv_method_id, &
                                                            xc_rho_smooth_id
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nmo_occ
      LOGICAL                                            :: lsd
      REAL(kind=dp)                                      :: alpha, tddfpt_fac
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_dbcsr_p_type)                              :: F_ia_munu
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: rho_ia_ao
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: rho_ia_fm
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_ia_g, rho_ia_g2, rho_ia_r, &
                                                            rho_ia_r2, v_xc
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_rho_type), POINTER                         :: rho_ia_struct
      TYPE(section_vals_type), POINTER                   :: xc_fun_section
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type), POINTER                     :: xc_rho1_set

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nvects = SIZE(evects, 2)
      lsd = (nspins > 1) .OR. is_rks_triplets

      ! *** initialize variables ***
      NULLIFY (blacs_env, para_env, pw_env)
      CALL get_qs_env(qs_env, blacs_env=blacs_env, para_env=para_env, pw_env=pw_env)

      NULLIFY (auxbas_pw_pool)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      ALLOCATE (nmo_occ(nspins))
      DO ispin = 1, nspins
         CALL cp_fm_get_info(gs_mos(ispin)%matrix_mos_occ, nrow_global=nao, ncol_global=nmo_occ(ispin))
      END DO

      NULLIFY (rho_ia_ao)
      CALL cp_dbcsr_allocate_matrix_set(rho_ia_ao, nspins)
      DO ispin = 1, nspins
         CALL cp_dbcsr_init_p(rho_ia_ao(ispin)%matrix)
         CALL cp_dbcsr_copy(rho_ia_ao(ispin)%matrix, matrix_s)
      END DO

      NULLIFY (rho_ia_struct, rho_ia_r, rho_ia_g)
      CALL qs_rho_create(rho_ia_struct)
      CALL qs_rho_set(rho_ia_struct, rho_ao=rho_ia_ao)

      CALL qs_rho_rebuild(rho_ia_struct, qs_env, rebuild_ao=.FALSE., rebuild_grids=.TRUE.)
      CALL qs_rho_get(rho_ia_struct, rho_r=rho_ia_r, rho_g=rho_ia_g)

      NULLIFY (F_ia_munu%matrix)
      CALL cp_dbcsr_init_p(F_ia_munu%matrix)
      CALL cp_dbcsr_copy(F_ia_munu%matrix, matrix_s)

      NULLIFY (fm_struct, rho_ia_fm)
      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, context=blacs_env, para_env=para_env)
      CALL cp_fm_create(rho_ia_fm, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      IF (is_rks_triplets) THEN
         ALLOCATE (rho_ia_r2(2))
         ALLOCATE (rho_ia_g2(2))
         rho_ia_r2(1)%pw => rho_ia_r(1)%pw
         rho_ia_r2(2)%pw => rho_ia_r(1)%pw
         rho_ia_g2(1)%pw => rho_ia_g(1)%pw
         rho_ia_g2(2)%pw => rho_ia_g(1)%pw
      ELSE
         ALLOCATE (rho_ia_r2(nspins))
         ALLOCATE (rho_ia_g2(nspins))
         DO ispin = 1, nspins
            rho_ia_r2(ispin)%pw => rho_ia_r(ispin)%pw
            rho_ia_g2(ispin)%pw => rho_ia_g(ispin)%pw
         END DO
      END IF

      ALLOCATE (v_xc(nspins))
      DO ispin = 1, SIZE(v_xc)
         NULLIFY (v_xc(ispin)%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, v_xc(ispin)%pw, use_data=REALDATA3D, in_space=REALSPACE)
      END DO

      NULLIFY (xc_rho1_set)
      CALL xc_rho_set_create(xc_rho1_set, rho_ia_r(1)%pw%pw_grid%bounds_local, &
                             rho_cutoff=section_get_rval(xc_section, "DENSITY_CUTOFF"), &
                             drho_cutoff=section_get_rval(xc_section, "GRADIENT_CUTOFF"), &
                             tau_cutoff=section_get_rval(xc_section, "TAU_CUTOFF"))

      xc_fun_section => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      xc_deriv_method_id = section_get_ival(xc_section, "XC_GRID%XC_DERIV")
      xc_rho_smooth_id = section_get_ival(xc_section, "XC_GRID%XC_SMOOTH_RHO")

      ! get arguments needed by the given functional
      needs = xc_functionals_get_needs(functionals=xc_fun_section, lsd=lsd, add_basic_components=.TRUE.)

      alpha = 1.0_dp
      tddfpt_fac = 0.0_dp
      ! tddfpt_fac is taken into account in spin-restricted case only
      IF (nspins == 1) THEN
         IF (is_rks_triplets) THEN
            ! K_{triplets} = K_{alpha,alpha} - K_{alpha,beta}
            tddfpt_fac = -1.0_dp
         ELSE
            !                                                 alpha              tddfpt_fac
            ! K_{singlets} = K_{alpha,alpha} + K_{alpha,beta} = 2 * K_{alpha,alpha} + 0 * K_{alpha,beta},
            ! due to the following relation : K_{alpha,alpha,singlets} == K_{alpha,beta,singlets}
            alpha = 2.0_dp
         END IF
      END IF

      ! *** compute kernel terms ***
      DO ivect = 1, nvects
         DO ispin = 1, nspins
            CALL cp_gemm('N', 'T', nao, nao, nmo_occ(ispin), 0.5_dp, &
                         gs_mos(ispin)%matrix_mos_occ, evects(ispin, ivect)%matrix, 0.0_dp, rho_ia_fm)
            CALL cp_gemm('N', 'T', nao, nao, nmo_occ(ispin), 0.5_dp, &
                         evects(ispin, ivect)%matrix, gs_mos(ispin)%matrix_mos_occ, 1.0_dp, rho_ia_fm)

            CALL copy_fm_to_dbcsr(rho_ia_fm, rho_ia_ao(ispin)%matrix, keep_sparsity=.TRUE.)
         END DO
         CALL qs_rho_update_rho(rho_ia_struct, qs_env)

         CALL xc_rho_set_update(rho_set=xc_rho1_set, rho_r=rho_ia_r2, rho_g=rho_ia_g2, &
                                tau=null(), needs=needs, xc_deriv_method_id=xc_deriv_method_id, &
                                xc_rho_smooth_id=xc_rho_smooth_id, pw_pool=auxbas_pw_pool)

         DO ispin = 1, SIZE(v_xc)
            CALL pw_zero(v_xc(ispin)%pw)
         END DO

         CALL xc_calc_2nd_deriv(v_xc=v_xc, deriv_set=xc_deriv_set, rho_set=xc_rho_set, &
                                rho1_set=xc_rho1_set, pw_pool=auxbas_pw_pool, &
                                xc_section=xc_section, gapw=.FALSE., tddfpt_fac=tddfpt_fac)

         DO ispin = 1, nspins
            CALL pw_scale(v_xc(ispin)%pw, v_xc(ispin)%pw%pw_grid%dvol)

            CALL cp_dbcsr_set(F_ia_munu%matrix, 0.0_dp)
            CALL integrate_v_rspace(v_rspace=v_xc(ispin), hmat=F_ia_munu, qs_env=qs_env, calculate_forces=.FALSE., gapw=.FALSE.)

            CALL cp_dbcsr_sm_fm_multiply(F_ia_munu%matrix, gs_mos(ispin)%matrix_mos_occ, Aop_evects(ispin, ivect)%matrix, &
                                         ncol=nmo_occ(ispin), alpha=alpha, beta=1.0_dp)
         END DO
      END DO

      ! *** clean up ***
      CALL xc_rho_set_release(xc_rho1_set)

      DO ispin = SIZE(v_xc), 1, -1
         CALL pw_pool_give_back_pw(auxbas_pw_pool, v_xc(ispin)%pw)
      END DO
      DEALLOCATE (v_xc)
      DEALLOCATE (rho_ia_g2)
      DEALLOCATE (rho_ia_r2)

      CALL cp_fm_release(rho_ia_fm)

      CALL cp_dbcsr_deallocate_matrix(F_ia_munu%matrix)
      CALL qs_rho_release(rho_ia_struct)

      DEALLOCATE (nmo_occ)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_apply_kernel

! **************************************************************************************************
!> \brief  Compute linear combination of half-transformed four-centre integrals (mu j || nu b) .
!> \param K_ia_munu  half-transformed integrals (initialised on exit)
!> \param evects     trial vectors C_{1,i}
!> \param gs_mos     molecular orbitals optimised for the ground state (only occupied molecular
!>                   orbitals [field %matrix_mos_occ] are needed)
!> \param qs_env     Quickstep environment
!> \param matrix_s   overlap matrix
!> \par History
!>    * 05.2016 compute all coulomb terms in one go [Sergey Chulkov]
!> \note Based on the subroutine kpp1_calc_k_p_p1() which was originally created by
!>       Mohamed Fawzi on 10.2002.
! **************************************************************************************************
   SUBROUTINE tddfpt_apply_hfx(K_ia_munu, evects, gs_mos, qs_env, matrix_s)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: K_ia_munu, evects
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_dbcsr_type), POINTER                       :: matrix_s

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_apply_hfx', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, ivect, nao, nspins, nvects
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nmo_occ
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: K_ia_munu_dbcsr, rho_ia_ao
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: rho_ia_fm
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_rho_type), POINTER                         :: rho_ia_struct

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nvects = SIZE(evects, 2)

      ! *** initialize variables ***
      NULLIFY (blacs_env, para_env, pw_env)
      CALL get_qs_env(qs_env, blacs_env=blacs_env, para_env=para_env, pw_env=pw_env)

      NULLIFY (auxbas_pw_pool, poisson_env)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, poisson_env=poisson_env)

      ALLOCATE (nmo_occ(nspins))
      DO ispin = 1, nspins
         CALL cp_fm_get_info(gs_mos(ispin)%matrix_mos_occ, nrow_global=nao, ncol_global=nmo_occ(ispin))
      END DO

      NULLIFY (rho_ia_ao)
      CALL cp_dbcsr_allocate_matrix_set(rho_ia_ao, nspins)
      DO ispin = 1, nspins
         CALL cp_dbcsr_init_p(rho_ia_ao(ispin)%matrix)
         CALL cp_dbcsr_copy(rho_ia_ao(ispin)%matrix, matrix_s)
      END DO

      NULLIFY (rho_ia_struct)
      CALL qs_rho_create(rho_ia_struct)
      CALL qs_rho_set(rho_ia_struct, rho_ao=rho_ia_ao)

      NULLIFY (K_ia_munu_dbcsr)
      CALL cp_dbcsr_allocate_matrix_set(K_ia_munu_dbcsr, nspins)
      DO ispin = 1, nspins
         CALL cp_dbcsr_init_p(K_ia_munu_dbcsr(ispin)%matrix)
         CALL cp_dbcsr_copy(K_ia_munu_dbcsr(ispin)%matrix, matrix_s)
      END DO

      NULLIFY (fm_struct, rho_ia_fm)
      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, context=blacs_env, para_env=para_env)
      CALL cp_fm_create(rho_ia_fm, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      ! *** compute hfx terms ***
      DO ivect = 1, nvects
         DO ispin = 1, nspins
            CALL cp_gemm('N', 'T', nao, nao, nmo_occ(ispin), 0.5_dp, &
                         gs_mos(ispin)%matrix_mos_occ, evects(ispin, ivect)%matrix, 0.0_dp, rho_ia_fm)
            CALL cp_gemm('N', 'T', nao, nao, nmo_occ(ispin), 0.5_dp, &
                         evects(ispin, ivect)%matrix, gs_mos(ispin)%matrix_mos_occ, 1.0_dp, rho_ia_fm)

            CALL copy_fm_to_dbcsr(rho_ia_fm, rho_ia_ao(ispin)%matrix, keep_sparsity=.TRUE.)
            CALL cp_dbcsr_set(K_ia_munu_dbcsr(ispin)%matrix, 0.0_dp)
         END DO

         CALL tddft_hfx_matrix(qs_env, rho_ia_struct, K_ia_munu_dbcsr)

         DO ispin = 1, nspins
            CALL copy_dbcsr_to_fm(K_ia_munu_dbcsr(ispin)%matrix, K_ia_munu(ispin, ivect)%matrix)
         END DO
      END DO

      ! *** clean up ***
      CALL cp_fm_release(rho_ia_fm)
      CALL cp_dbcsr_deallocate_matrix_set(K_ia_munu_dbcsr)
      CALL qs_rho_release(rho_ia_struct)

      DEALLOCATE (nmo_occ)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_apply_hfx

! **************************************************************************************************
!> \brief Compute the action matrix-vector product and half-transformed four-centre integrals.
!> \param Aop_evects      hat{A} C_{1,i} (allocated and initialised on exit)
!> \param K_ia_munu       half-transformed 4-centre integrals (allocated and initialised on exit)
!> \param evects          TDDFPT trial vectors
!> \param gs_mos          molecular orbitals optimised for the ground state
!> \param is_rks_triplets indicates that the triplet excited states calculation using
!>                        spin-unpolarised molecular orbitals has been requested
!> \param K_ia_template   a dense matrix with the same dimension as K_ia_munu
!>                        (used as a template to allocate unallocated K_ia_munu-s;
!>                        NULL means no exact-exchange terms are needed)
!> \param qs_env          Quickstep environment
!> \param matrix_ks       Kohn-Sham matrix
!> \param matrix_s        overlap matrix
!> \param xc_rho_set      a variable to compute the 2nd derivative of XC-functional
!> \param xc_deriv_set    a variable to compute the 2nd derivative of XC-functional
!> \param xc_section      XC input section
!> \par History
!>    * 06.2016 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_compute_Aop_evects(Aop_evects, K_ia_munu, evects, gs_mos, is_rks_triplets, K_ia_template, &
                                        qs_env, matrix_ks, matrix_s, xc_rho_set, xc_deriv_set, xc_section)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: Aop_evects, K_ia_munu, evects
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      LOGICAL, INTENT(in)                                :: is_rks_triplets
      TYPE(cp_fm_type), POINTER                          :: K_ia_template
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_dbcsr_p_type), DIMENSION(:), INTENT(in)    :: matrix_ks
      TYPE(cp_dbcsr_type), POINTER                       :: matrix_s
      TYPE(xc_rho_set_type), POINTER                     :: xc_rho_set
      TYPE(xc_derivative_set_type), POINTER              :: xc_deriv_set
      TYPE(section_vals_type), POINTER                   :: xc_section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_compute_Aop_evects', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, ivect, nspins, nvects

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nvects = SIZE(evects, 2)

      CPASSERT(nspins > 0)
      CPASSERT(SIZE(Aop_evects, 1) == nspins)
      CPASSERT(SIZE(Aop_evects, 2) == nvects)
      CPASSERT(SIZE(gs_mos) == nspins)

      IF (nvects > 0) THEN
         DO ivect = 1, nvects
            DO ispin = 1, nspins
               IF (.NOT. ASSOCIATED(Aop_evects(ispin, ivect)%matrix)) &
                  CALL cp_fm_create(Aop_evects(ispin, ivect)%matrix, evects(ispin, 1)%matrix%matrix_struct)
            END DO
         END DO

         IF (ASSOCIATED(K_ia_template)) THEN
            DO ivect = 1, nvects
               DO ispin = 1, nspins
                  IF (.NOT. ASSOCIATED(K_ia_munu(ispin, ivect)%matrix)) &
                     CALL cp_fm_create(K_ia_munu(ispin, ivect)%matrix, K_ia_template%matrix_struct)
               END DO
            END DO
         END IF

         ! initialise hat{A} C_{1,i} by orbital energy difference
         CALL tddfpt_init_by_energy_diff(Aop_evects=Aop_evects, evects=evects, gs_mos=gs_mos, &
                                         matrix_ks=matrix_ks, matrix_s=matrix_s)

         ! add contributions from coulomb terms
         IF (.NOT. is_rks_triplets) &
            CALL tddfpt_apply_coulomb(Aop_evects=Aop_evects, evects=evects, &
                                      gs_mos=gs_mos, qs_env=qs_env, matrix_s=matrix_s)

         ! add contributions from the adiabatic TDDFT kernel
         CALL tddfpt_apply_kernel(Aop_evects=Aop_evects, evects=evects, gs_mos=gs_mos, &
                                  is_rks_triplets=is_rks_triplets, &
                                  xc_rho_set=xc_rho_set, xc_deriv_set=xc_deriv_set, &
                                  qs_env=qs_env, matrix_s=matrix_s, xc_section=xc_section)

         ! add contributions from exact-exchange terms
         IF (ASSOCIATED(K_ia_template)) &
            CALL tddfpt_apply_hfx(K_ia_munu=K_ia_munu, evects=evects, gs_mos=gs_mos, qs_env=qs_env, matrix_s=matrix_s)
      END IF

      CALL timestop(handle)
   END SUBROUTINE tddfpt_compute_Aop_evects

! **************************************************************************************************
!> \brief Solve eigenproblem for the reduced action matrix and find new Ritz eigenvectors and
!>        eigenvalues.
!> \param ritz_vects      Ritz eigenvectors (initialised on exit)
!> \param Aop_ritz        action matrix -- Ritz vector product
!> \param evals           Ritz eigenvalues
!> \param Aop_evects      action matrix -- trial vector product
!> \param K_ia_munu       half-transformed 4-centre integrals
!> \param evects          TDDFPT trial vectors
!> \param gs_mos          molecular orbitals optimized for the ground state
!> \par History
!>    * 06.2016 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_compute_ritz_vects(ritz_vects, Aop_ritz, evals, Aop_evects, K_ia_munu, evects, gs_mos)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: ritz_vects, Aop_ritz
      REAL(kind=dp), DIMENSION(:), INTENT(out)           :: evals
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: Aop_evects, K_ia_munu, evects
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_compute_ritz_vects', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, ivect, jvect, nao, nkvs, &
                                                            nrvs, nspins
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nmo_occ
      LOGICAL                                            :: do_hfx
      REAL(kind=dp)                                      :: act, alpha_hfx, trace
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: Atilde
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: Atilde_fm, evects_Atilde_fm, rho_jb_fm
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nkvs = SIZE(evects, 2)
      nrvs = SIZE(ritz_vects, 2)
      do_hfx = ASSOCIATED(K_ia_munu(1, 1)%matrix)

      NULLIFY (blacs_env, para_env)
      CALL cp_fm_get_info(gs_mos(1)%matrix_mos_occ, nrow_global=nao, context=blacs_env, para_env=para_env)

      ! *** compute upper-diagonal reduced action matrix ***
      ALLOCATE (Atilde(nkvs, nkvs))
      DO ivect = 1, nkvs
         DO jvect = 1, ivect
            act = 0.0_dp
            DO ispin = 1, nspins
               CALL cp_fm_trace(Aop_evects(ispin, jvect)%matrix, evects(ispin, ivect)%matrix, trace)
               act = act+trace
            END DO

            Atilde(jvect, ivect) = act
         END DO
      END DO

      ! *** apply exact exchange terms ***
      IF (do_hfx) THEN
         ALLOCATE (nmo_occ(nspins))
         DO ispin = 1, nspins
            nmo_occ(ispin) = SIZE(gs_mos(ispin)%evals_occ)
         END DO

         NULLIFY (fm_struct, rho_jb_fm)
         CALL cp_fm_get_info(K_ia_munu(1, 1)%matrix, matrix_struct=fm_struct)
         CALL cp_fm_create(rho_jb_fm, fm_struct)

         IF (nspins > 1) THEN
            alpha_hfx = 2.0_dp
         ELSE
            alpha_hfx = 4.0_dp
         END IF

         DO ivect = 1, nkvs
            DO ispin = 1, nspins
               CALL cp_gemm('N', 'T', nao, nao, nmo_occ(ispin), 1.0_dp, &
                            gs_mos(ispin)%matrix_mos_occ, evects(ispin, ivect)%matrix, 0.0_dp, rho_jb_fm)

               DO jvect = 1, ivect
                  CALL cp_fm_trace(K_ia_munu(ispin, jvect)%matrix, rho_jb_fm, trace)
                  Atilde(jvect, ivect) = Atilde(jvect, ivect)+alpha_hfx*trace
               END DO
            END DO
         END DO

         CALL cp_fm_release(rho_jb_fm)
         DEALLOCATE (nmo_occ)
      END IF

      ! *** solve eigenproblem for reduced matrices ***
      NULLIFY (fm_struct, Atilde_fm, evects_Atilde_fm)
      CALL cp_fm_struct_create(fm_struct, nrow_global=nkvs, ncol_global=nkvs, context=blacs_env, para_env=para_env)
      CALL cp_fm_create(Atilde_fm, fm_struct)
      CALL cp_fm_create(evects_Atilde_fm, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      CALL cp_fm_set_submatrix(Atilde_fm, Atilde)
      CALL choose_eigv_solver(Atilde_fm, evects_Atilde_fm, evals(1:nkvs))
      CALL cp_fm_get_submatrix(evects_Atilde_fm, Atilde)

      CALL cp_fm_release(evects_Atilde_fm)
      CALL cp_fm_release(Atilde_fm)

      ! *** compute Ritz vectors ***
      DO ivect = 1, nrvs
         DO ispin = 1, nspins
            CALL cp_fm_set_all(ritz_vects(ispin, ivect)%matrix, 0.0_dp)
            CALL cp_fm_set_all(Aop_ritz(ispin, ivect)%matrix, 0.0_dp)
         END DO

         DO jvect = 1, nkvs
            act = Atilde(jvect, ivect)
            DO ispin = 1, nspins
               CALL cp_fm_scale_and_add(1.0_dp, ritz_vects(ispin, ivect)%matrix, act, evects(ispin, jvect)%matrix)
               CALL cp_fm_scale_and_add(1.0_dp, Aop_ritz(ispin, ivect)%matrix, act, Aop_evects(ispin, jvect)%matrix)
            END DO
         END DO
      END DO

      DEALLOCATE (Atilde)
      CALL timestop(handle)
   END SUBROUTINE tddfpt_compute_ritz_vects

! **************************************************************************************************
!> \brief Expand Krylov space by computing residual vectors.
!> \param residual_vects  residual vectors
!> \param evals           Ritz eigenvalues
!> \param Aop_evects      action matrix -- Ritz vector product
!> \param K_ia_munu       half-transformed 4-centre integrals (or list of NULL-s HFX terms
!>                        should not be taken into account for residual vectors)
!> \param evects          Ritz eigenvectors
!> \param gs_mos          molecular orbitals optimised for the ground state
!> \param guess_vectors   initial guess vectors
!> \param matrix_s        overlap matrix
!> \par History
!>    * 06.2016 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_compute_residual_vects(residual_vects, evals, Aop_evects, K_ia_munu, evects, &
                                            gs_mos, guess_vectors, matrix_s)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: residual_vects
      REAL(kind=dp), DIMENSION(:), INTENT(in)            :: evals
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: Aop_evects, K_ia_munu, evects
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      TYPE(tddfpt_guess_vectors), INTENT(in)             :: guess_vectors
      TYPE(cp_dbcsr_type), POINTER                       :: matrix_s

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_compute_residual_vects', &
         routineP = moduleN//':'//routineN
      REAL(kind=dp), PARAMETER :: eref_scale = 0.99_dp, threshold = 16.0_dp*EPSILON(1.0_dp)

      INTEGER                                            :: handle, imo_occ, imo_virt, irv, ispin, &
                                                            istate, nao, nkvs, nrvs, nspins, &
                                                            nstates_total
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nmo_occ
      LOGICAL                                            :: do_hfx
      REAL(kind=dp)                                      :: alpha_hfx, eref, K_ia_jb, minus_lambda, &
                                                            weight
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: Ab_minus_eb
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: rho_jb_fm
      TYPE(tddfpt_fm_vectors), ALLOCATABLE, DIMENSION(:) :: Ab_cols, residual_cols

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nkvs = SIZE(evects, 2)
      nrvs = SIZE(residual_vects, 2)
      nstates_total = SIZE(guess_vectors%evals)

      do_hfx = ASSOCIATED(K_ia_munu(1, 1)%matrix)
      IF (nspins > 1) THEN
         alpha_hfx = 2.0_dp
      ELSE
         alpha_hfx = 4.0_dp
      END IF

      IF (nrvs > 0) THEN
         DO istate = 1, nrvs
            DO ispin = 1, nspins
               IF (.NOT. ASSOCIATED(residual_vects(ispin, istate)%matrix)) &
                  CALL cp_fm_create(residual_vects(ispin, istate)%matrix, evects(ispin, 1)%matrix%matrix_struct)
            END DO
         END DO

         ALLOCATE (nmo_occ(nspins))
         ALLOCATE (Ab_minus_eb(nspins))
         DO ispin = 1, nspins
            NULLIFY (fm_struct, Ab_minus_eb(ispin)%matrix)
            CALL cp_fm_get_info(gs_mos(ispin)%matrix_mos_occ, ncol_global=nmo_occ(ispin), matrix_struct=fm_struct)
            CALL cp_fm_create(Ab_minus_eb(ispin)%matrix, fm_struct)
         END DO

         ALLOCATE (Ab_cols(nspins), residual_cols(nspins))
         DO ispin = 1, nspins
            NULLIFY (fm_struct)
            CALL cp_fm_get_info(gs_mos(ispin)%vectors_mo_occ(1)%matrix, matrix_struct=fm_struct)

            ALLOCATE (Ab_cols(ispin)%vects(nmo_occ(ispin)))
            ALLOCATE (residual_cols(ispin)%vects(nmo_occ(ispin)))
            DO imo_occ = 1, nmo_occ(ispin)
               NULLIFY (Ab_cols(ispin)%vects(imo_occ)%matrix, residual_cols(ispin)%vects(imo_occ)%matrix)
               CALL cp_fm_create(Ab_cols(ispin)%vects(imo_occ)%matrix, fm_struct)
               CALL cp_fm_create(residual_cols(ispin)%vects(imo_occ)%matrix, fm_struct)
            END DO
         END DO

         NULLIFY (rho_jb_fm)
         IF (do_hfx) THEN
            NULLIFY (fm_struct)
            CALL cp_fm_get_info(K_ia_munu(1, 1)%matrix, nrow_global=nao, matrix_struct=fm_struct)
            CALL cp_fm_create(rho_jb_fm, fm_struct)
         END IF

         ! *** actually compute residual vectors ***
         DO irv = 1, nrvs
            minus_lambda = -evals(irv)
            DO ispin = 1, nspins
               CALL cp_dbcsr_sm_fm_multiply(matrix_s, evects(ispin, irv)%matrix, Ab_minus_eb(ispin)%matrix, &
                                            ncol=nmo_occ(ispin), alpha=minus_lambda, beta=0.0_dp)
               CALL cp_fm_scale_and_add(1.0_dp, Ab_minus_eb(ispin)%matrix, 1.0_dp, Aop_evects(ispin, irv)%matrix)
            END DO

            DO ispin = 1, nspins
               DO imo_occ = nmo_occ(ispin), 1, -1
                  CALL cp_fm_to_fm(Ab_minus_eb(ispin)%matrix, Ab_cols(ispin)%vects(imo_occ)%matrix, &
                                   ncol=1, source_start=imo_occ, target_start=1)
                  CALL cp_fm_set_all(residual_cols(ispin)%vects(imo_occ)%matrix, 0.0_dp)
               END DO
            END DO

            DO istate = 1, nstates_total
               ispin = guess_vectors%ispins(istate)
               imo_occ = guess_vectors%imos_occ(istate)
               imo_virt = guess_vectors%imos_virt(istate)

               eref = guess_vectors%evals(istate)+minus_lambda
               IF (ABS(eref) < threshold) &
                  eref = guess_vectors%evals(istate)+eref_scale*minus_lambda

               CALL cp_fm_trace(gs_mos(ispin)%vectors_mo_virt(imo_virt)%matrix, Ab_cols(ispin)%vects(imo_occ)%matrix, weight)

               IF (do_hfx) THEN
                  CALL cp_gemm('N', 'T', nao, nao, 1, alpha_hfx, gs_mos(ispin)%vectors_mo_occ(imo_occ)%matrix, &
                               gs_mos(ispin)%vectors_mo_virt(imo_virt)%matrix, 0.0_dp, rho_jb_fm)

                  CALL cp_fm_trace(rho_jb_fm, K_ia_munu(ispin, irv)%matrix, K_ia_jb)
                  weight = weight+K_ia_jb
               END IF

               weight = weight/eref
               CALL cp_fm_scale_and_add(1.0_dp, residual_cols(ispin)%vects(imo_occ)%matrix, &
                                        weight, gs_mos(ispin)%vectors_mo_virt(imo_virt)%matrix)
            END DO

            DO ispin = 1, nspins
               DO imo_occ = nmo_occ(ispin), 1, -1
                  CALL cp_fm_to_fm(residual_cols(ispin)%vects(imo_occ)%matrix, residual_vects(ispin, irv)%matrix, &
                                   ncol=1, source_start=1, target_start=imo_occ)
               END DO
            END DO
         END DO

         IF (ASSOCIATED(rho_jb_fm)) &
            CALL cp_fm_release(rho_jb_fm)

         DO ispin = SIZE(Ab_cols), 1, -1
            DO imo_occ = SIZE(Ab_cols(ispin)%vects), 1, -1
               CALL cp_fm_release(residual_cols(ispin)%vects(imo_occ)%matrix)
               CALL cp_fm_release(Ab_cols(ispin)%vects(imo_occ)%matrix)
            END DO
            DEALLOCATE (residual_cols(ispin)%vects, Ab_cols(ispin)%vects)
         END DO
         DEALLOCATE (residual_cols, Ab_cols)

         DO ispin = SIZE(Ab_minus_eb), 1, -1
            CALL cp_fm_release(Ab_minus_eb(ispin)%matrix)
         END DO
         DEALLOCATE (Ab_minus_eb)
         DEALLOCATE (nmo_occ)
      END IF

      CALL timestop(handle)
   END SUBROUTINE tddfpt_compute_residual_vects

! **************************************************************************************************
!> \brief Print final TDDFPT excitation energies and analysis.
!> \param log_unit output unit
!> \param evects   TDDFPT trial vectors (SIZE(evects,1) -- number of spins;
!>                 SIZE(evects,2) -- number of excited states to print)
!> \param evals    TDDFPT eigenvalues
!> \param mult     multiplicity
!> \param gs_mos   molecular orbitals optimised for the ground state
!> \param matrix_s overlap matrix and their firts derivatives over x, y, and z directions
!> \param min_amplitude the smallest excitation amplitude to print
!> \par History
!>    * 05.2016 created [Sergey Chulkov]
!>    * 06.2016 transition dipole moments and oscillator strengths [Sergey Chulkov]
!>    * 07.2016 spin-unpolarised electron density [Sergey Chulkov]
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
   SUBROUTINE tddfpt_print_summary(log_unit, evects, evals, mult, gs_mos, matrix_s, min_amplitude)
      INTEGER, INTENT(in)                                :: log_unit
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: evects
      REAL(kind=dp), DIMENSION(:), INTENT(in)            :: evals
      INTEGER, INTENT(in)                                :: mult
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s
      REAL(kind=dp), INTENT(in)                          :: min_amplitude

      CHARACTER(len=*), PARAMETER :: routineN = 'tddfpt_print_summary', &
         routineP = moduleN//':'//routineN
      INTEGER, PARAMETER                                 :: nderivs = 3

      CHARACTER(len=1)                                   :: lsd_str, spin_occ, spin_virt
      CHARACTER(len=20)                                  :: mult_str
      INTEGER :: ideriv, iexc, imo_occ, imo_virt, ind, ispin, istate, nao, nmo_virt_occ, &
         nmo_virt_occ_alpha, nmo_virt_occ_spin, nspins, nstates, state_spin
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: inds, nmo_occ, nmo_virt
      REAL(kind=dp)                                      :: oscillator_strength
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: weights_abs_1d
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: trans_dipole
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: ediff_inv, ediff_inv_weights, &
                                                            S_mos_virt, scaled_evect, weights_fm
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:, :)   :: dS_mos_occ
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(tddfpt_allocatable_matrix), ALLOCATABLE, &
         DIMENSION(:)                                    :: weights_2d

      EXTERNAL :: dcopy

      nspins = SIZE(evects, 1)
      nstates = SIZE(evects, 2)
      CALL cp_fm_get_info(evects(1, 1)%matrix, nrow_global=nao, context=blacs_env, para_env=para_env)

      IF (nspins > 1) THEN
         lsd_str = 'U'
      ELSE
         lsd_str = 'R'
      END IF

      ALLOCATE (nmo_occ(nspins))
      ALLOCATE (nmo_virt(nspins))
      DO ispin = 1, nspins
         nmo_occ(ispin) = SIZE(gs_mos(ispin)%evals_occ)
         nmo_virt(ispin) = SIZE(gs_mos(ispin)%evals_virt)
      END DO

      ALLOCATE (S_mos_virt(nspins))
      DO ispin = 1, nspins
         NULLIFY (fm_struct, S_mos_virt(ispin)%matrix)
         CALL cp_fm_get_info(gs_mos(ispin)%matrix_mos_virt, matrix_struct=fm_struct)
         CALL cp_fm_create(S_mos_virt(ispin)%matrix, fm_struct)
      END DO

      ALLOCATE (weights_2d(nspins))
      DO ispin = 1, nspins
         ALLOCATE (weights_2d(ispin)%matrix(nmo_virt(ispin), nmo_occ(ispin)))
      END DO

      NULLIFY (fm_struct)
      ALLOCATE (weights_fm(nspins))
      DO ispin = 1, nspins
         CALL cp_fm_struct_create(fm_struct, nrow_global=nmo_virt(ispin), ncol_global=nmo_occ(ispin), &
                                  para_env=para_env, context=blacs_env)
         NULLIFY (weights_fm(ispin)%matrix)
         CALL cp_fm_create(weights_fm(ispin)%matrix, fm_struct)
         CALL cp_fm_struct_release(fm_struct)
      END DO

      ALLOCATE (dS_mos_occ(nderivs, nspins), ediff_inv(nspins), ediff_inv_weights(nspins), scaled_evect(nspins))
      DO ispin = 1, nspins
         NULLIFY (fm_struct, ediff_inv(ispin)%matrix, ediff_inv_weights(ispin)%matrix, scaled_evect(ispin)%matrix)
         CALL cp_fm_get_info(gs_mos(ispin)%matrix_mos_occ, matrix_struct=fm_struct)

         DO ideriv = 1, nderivs
            NULLIFY (dS_mos_occ(ideriv, ispin)%matrix)
            CALL cp_fm_create(dS_mos_occ(ideriv, ispin)%matrix, fm_struct)
         END DO
         CALL cp_fm_create(ediff_inv(ispin)%matrix, fm_struct)
         CALL cp_fm_create(ediff_inv_weights(ispin)%matrix, fm_struct)
         CALL cp_fm_create(scaled_evect(ispin)%matrix, fm_struct)
      END DO

      DO ispin = 1, nspins
         CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, gs_mos(ispin)%matrix_mos_virt, S_mos_virt(ispin)%matrix, &
                                      ncol=nmo_virt(ispin), alpha=1.0_dp, beta=0.0_dp)
         DO ideriv = 1, nderivs
            CALL cp_dbcsr_sm_fm_multiply(matrix_s(ideriv+1)%matrix, gs_mos(ispin)%matrix_mos_occ, &
                                         dS_mos_occ(ideriv, ispin)%matrix, ncol=nmo_occ(ispin), alpha=1.0_dp, beta=0.0_dp)
         END DO
      END DO

      DO ispin = 1, nspins
         DO imo_occ = nmo_occ(ispin), 1, -1
            DO imo_virt = nmo_virt(ispin), 1, -1
               weights_2d(ispin)%matrix(imo_virt, imo_occ) = &
                  1.0_dp/(gs_mos(ispin)%evals_virt(imo_virt)-gs_mos(ispin)%evals_occ(imo_occ))
            END DO
         END DO
         CALL cp_fm_set_submatrix(ediff_inv(ispin)%matrix, weights_2d(ispin)%matrix)
      END DO

      ! *** summary information ***
      IF (log_unit > 0) THEN
         CALL integer_to_string(mult, mult_str)
         WRITE (log_unit, '(/,1X,A1,A,1X,A,/)') lsd_str, "-TDDFT states of multiplicity", TRIM(mult_str)

         WRITE (log_unit, '(3X,A,T14,A,T31,A,T63,A)') "State", "Excitation", "Transition dipole (a.u.)", "Oscillator"
         WRITE (log_unit, '(3X,A,T14,A,T33,A,T43,A,T53,A,T61,A)') "number", "energy (eV)", "x", "y", "z", "strength (a.u.)"
         WRITE (log_unit, '(1X,74("-"))')
      END IF

      ALLOCATE (trans_dipole(nderivs, nspins))
      DO istate = 1, nstates
         DO ispin = 1, nspins
            CALL cp_gemm('T', 'N', nmo_virt(ispin), nmo_occ(ispin), nao, 1.0_dp, S_mos_virt(ispin)%matrix, &
                         evects(ispin, istate)%matrix, 0.0_dp, weights_fm(ispin)%matrix)

            CALL cp_fm_schur_product(weights_fm(ispin)%matrix, ediff_inv(ispin)%matrix, ediff_inv_weights(ispin)%matrix)
            CALL cp_gemm('N', 'N', nao, nmo_occ(ispin), nmo_virt(ispin), 1.0_dp, gs_mos(ispin)%matrix_mos_virt, &
                         ediff_inv_weights(ispin)%matrix, 0.0_dp, scaled_evect(ispin)%matrix)

            DO ideriv = 1, nderivs
               CALL cp_fm_trace(dS_mos_occ(ideriv, ispin)%matrix, scaled_evect(ispin)%matrix, trans_dipole(ideriv, ispin))
            END DO
         END DO

         DO ispin = 2, nspins
            DO ideriv = 1, nderivs
               trans_dipole(ideriv, 1) = trans_dipole(ideriv, 1)+trans_dipole(ideriv, ispin)
            END DO
         END DO
         oscillator_strength = 2.0_dp/3.0_dp*evals(istate)*DOT_PRODUCT(trans_dipole(:, 1), trans_dipole(:, 1))

         IF (log_unit > 0) &
            WRITE (log_unit, '(1X,I8,T14,F11.5,T28,3(1X,F9.4),T63,F10.5)') &
            istate, evals(istate)*evolt, trans_dipole(1:nderivs, 1), oscillator_strength
      END DO

      DO ispin = nspins, 1, -1
         CALL cp_fm_release(scaled_evect(ispin)%matrix)
         CALL cp_fm_release(ediff_inv_weights(ispin)%matrix)
         CALL cp_fm_release(ediff_inv(ispin)%matrix)

         DO ideriv = nderivs, 1, -1
            CALL cp_fm_release(dS_mos_occ(ideriv, ispin)%matrix)
         END DO
      END DO
      DEALLOCATE (trans_dipole, scaled_evect, ediff_inv_weights, ediff_inv, dS_mos_occ)

      ! *** excitation analysis ***
      IF (min_amplitude < 1.0_dp) THEN
         IF (nspins == 1) THEN
            state_spin = 1
            spin_virt = 'a'

            IF (mult == 3) THEN
               spin_occ = 'b'
            ELSE
               spin_occ = 'a'
            END IF
         END IF

         nmo_virt_occ_alpha = nmo_virt(1)*nmo_occ(1)
         nmo_virt_occ = DOT_PRODUCT(nmo_virt(:), nmo_occ(:))
         ALLOCATE (weights_abs_1d(nmo_virt_occ))
         ALLOCATE (inds(nmo_virt_occ))

         IF (log_unit > 0) THEN
            WRITE (log_unit, '(/,1X,A,/)') "Excitation analysis"

            WRITE (log_unit, '(3X,A,T17,A,T34,A,T49,A)') "State", "Occupied", "Virtual", "Excitation"
            WRITE (log_unit, '(3X,A,T18,A,T34,A,T49,A)') "number", "orbital", "orbital", "amplitude"
            WRITE (log_unit, '(1X,57("-"))')
         END IF

         DO istate = 1, nstates
            nmo_virt_occ = 0
            DO ispin = 1, nspins
               CALL cp_gemm('T', 'N', nmo_virt(ispin), nmo_occ(ispin), nao, 1.0_dp, S_mos_virt(ispin)%matrix, &
                            evects(ispin, istate)%matrix, 0.0_dp, weights_fm(ispin)%matrix)
               CALL cp_fm_get_submatrix(weights_fm(ispin)%matrix, weights_2d(ispin)%matrix)

               nmo_virt_occ_spin = nmo_virt(ispin)*nmo_occ(ispin)
               CALL dcopy(nmo_virt_occ_spin, weights_2d(ispin)%matrix, 1, &
                          weights_abs_1d(nmo_virt_occ+1:nmo_virt_occ+nmo_virt_occ_spin), 1)
               nmo_virt_occ = nmo_virt_occ+nmo_virt_occ_spin
            END DO

            weights_abs_1d = ABS(weights_abs_1d)
            CALL sort(weights_abs_1d, nmo_virt_occ, inds)

            IF (log_unit > 0) &
               WRITE (log_unit, '(1X,I8)') istate

            DO iexc = nmo_virt_occ, 1, -1
               IF (weights_abs_1d(iexc) < min_amplitude) EXIT

               ind = inds(iexc)-1

               IF (nspins > 1) THEN
                  IF (ind < nmo_virt_occ_alpha) THEN
                     state_spin = 1
                     spin_occ = 'a'
                  ELSE
                     state_spin = 2
                     ind = ind-nmo_virt_occ_alpha
                     spin_occ = 'b'
                  END IF
                  spin_virt = spin_occ
               END IF

               imo_occ = ind/nmo_virt(state_spin)+1
               imo_virt = MOD(ind, nmo_virt(state_spin))+1

               IF (log_unit > 0) &
                  WRITE (log_unit, '(T14,I8,1X,A1,A1,A1,T30,I8,1X,A1,A1,A1,T50,F9.6)') imo_occ, '(', spin_occ, ')', &
                  nmo_occ(state_spin)+imo_virt, '(', spin_virt, ')', weights_2d(state_spin)%matrix(imo_virt, imo_occ)
            END DO
         END DO

         DEALLOCATE (inds)
         DEALLOCATE (weights_abs_1d)
      END IF

      DO ispin = nspins, 1, -1
         CALL cp_fm_release(weights_fm(ispin)%matrix)
         DEALLOCATE (weights_2d(ispin)%matrix)
      END DO
      DEALLOCATE (weights_fm, weights_2d)

      DO ispin = nspins, 1, -1
         CALL cp_fm_release(S_mos_virt(ispin)%matrix)
      END DO
      DEALLOCATE (S_mos_virt)

      DEALLOCATE (nmo_virt, nmo_occ)
   END SUBROUTINE tddfpt_print_summary
END MODULE qs_tddfpt2_methods
