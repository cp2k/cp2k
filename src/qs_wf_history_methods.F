!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Storage of past states of the qs_env.
!>      Methods to interpolate (or actually normally extrapolate) the
!>      new guess for density and wavefunctions.
!> \note
!>      Most of the last snapshot should actually be in qs_env, but taking
!>      advantage of it would make the programming much convoluted
!> \par History
!>      02.2003 created [fawzi]
!>      11.2003 Joost VandeVondele : Implemented Nth order PS extrapolation
!>      02.2005 modified for KG_GPW [MI]
!> \author fawzi
! **************************************************************************************************
MODULE qs_wf_history_methods
   USE bibliography,                    ONLY: Kolafa2004,&
                                              Kuhne2007,&
                                              VandeVondele2005a,&
                                              cite_reference
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale,&
                                              cp_fm_scale_and_add
   USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                              cp_fm_pool_type,&
                                              fm_pools_create_fm_vect,&
                                              fm_pools_give_back_fm_vect
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr,&
                                              low_print_level
   USE dbcsr_api,                       ONLY: dbcsr_add,&
                                              dbcsr_copy,&
                                              dbcsr_deallocate_matrix,&
                                              dbcsr_p_type
   USE input_constants,                 ONLY: &
        wfi_aspc_nr, wfi_frozen_method_nr, wfi_linear_p_method_nr, wfi_linear_ps_method_nr, &
        wfi_linear_wf_method_nr, wfi_ps_method_nr, wfi_use_guess_method_nr, &
        wfi_use_prev_p_method_nr, wfi_use_prev_rho_r_method_nr, wfi_use_prev_wf_method_nr
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: binomial
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_copy
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_p_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_ks_types,                     ONLY: qs_ks_did_change
   USE qs_matrix_pools,                 ONLY: mpools_get,&
                                              qs_matrix_pools_type
   USE qs_mo_methods,                   ONLY: calculate_density_matrix,&
                                              make_basis_cholesky,&
                                              make_basis_lowdin,&
                                              make_basis_simple,&
                                              make_basis_sm
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_rho_methods,                  ONLY: qs_rho_update_rho
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_set,&
                                              qs_rho_type
   USE qs_scf_types,                    ONLY: ot_method_nr,&
                                              qs_scf_env_type
   USE qs_wf_history_types,             ONLY: qs_wf_history_type,&
                                              qs_wf_snapshot_type,&
                                              wfi_get_snapshot,&
                                              wfi_release
   USE scf_control_types,               ONLY: scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_wf_history_methods'
   INTEGER, SAVE, PRIVATE :: last_wfs_id = 0, last_wfi_id = 0

   PUBLIC :: wfi_create, wfi_update, wfi_create_for_kp, &
             wfi_extrapolate, wfi_get_method_label, &
             reorthogonalize_vectors, wfi_purge_history

CONTAINS

! **************************************************************************************************
!> \brief allocates and initialize a wavefunction snapshot
!> \param snapshot the snapshot to create
!> \par History
!>      02.2003 created [fawzi]
!>      02.2005 added wf_mol [MI]
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE wfs_create(snapshot)
      TYPE(qs_wf_snapshot_type), POINTER                 :: snapshot

      CHARACTER(len=*), PARAMETER :: routineN = 'wfs_create', routineP = moduleN//':'//routineN

      ALLOCATE (snapshot)
      last_wfs_id = last_wfs_id + 1
      snapshot%id_nr = last_wfs_id
      NULLIFY (snapshot%wf, snapshot%rho_r, &
               snapshot%rho_g, snapshot%rho_ao, snapshot%rho_ao_kp, &
               snapshot%overlap, snapshot%rho_frozen)
      snapshot%dt = 1.0_dp
      snapshot%ref_count = 1
   END SUBROUTINE wfs_create

! **************************************************************************************************
!> \brief updates the given snapshot
!> \param snapshot the snapshot to be updated
!> \param wf_history the history
!> \param qs_env the qs_env that should be snapshotted
!> \param dt the time of the snapshot (wrt. to the previous snapshot)
!> \par History
!>      02.2003 created [fawzi]
!>      02.2005 added kg_fm_mol_set for KG_GPW [MI]
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE wfs_update(snapshot, wf_history, qs_env, dt)
      TYPE(qs_wf_snapshot_type), POINTER                 :: snapshot
      TYPE(qs_wf_history_type), POINTER                  :: wf_history
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), INTENT(in), OPTIONAL                :: dt

      CHARACTER(len=*), PARAMETER :: routineN = 'wfs_update', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, img, ispin, nimg, nspins
      TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER     :: ao_mo_pools
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, rho_ao
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_g, rho_r
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_rho_type), POINTER                         :: rho

      CALL timeset(routineN, handle)

      NULLIFY (pw_env, auxbas_pw_pool, ao_mo_pools, dft_control, mos, mo_coeff, &
               rho, rho_r, rho_g, rho_ao, matrix_s)
      CALL get_qs_env(qs_env, pw_env=pw_env, &
                      dft_control=dft_control, rho=rho)
      CALL mpools_get(qs_env%mpools, ao_mo_fm_pools=ao_mo_pools)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      CPASSERT(ASSOCIATED(wf_history))
      CPASSERT(ASSOCIATED(dft_control))
      IF (.NOT. ASSOCIATED(snapshot)) THEN
         CALL wfs_create(snapshot)
      END IF
      CPASSERT(wf_history%ref_count > 0)
      CPASSERT(snapshot%ref_count > 0)

      nspins = dft_control%nspins
      snapshot%dt = 1.0_dp
      IF (PRESENT(dt)) snapshot%dt = dt
      IF (wf_history%store_wf) THEN
         CALL get_qs_env(qs_env, mos=mos)
         IF (.NOT. ASSOCIATED(snapshot%wf)) THEN
            CALL fm_pools_create_fm_vect(ao_mo_pools, snapshot%wf, &
                                         name="ws_snap"//TRIM(ADJUSTL(cp_to_string(snapshot%id_nr)))// &
                                         "ws")
            CPASSERT(nspins == SIZE(snapshot%wf))
         END IF
         DO ispin = 1, nspins
            CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff)
            CALL cp_fm_to_fm(mo_coeff, snapshot%wf(ispin)%matrix)
         END DO
      ELSE IF (ASSOCIATED(snapshot%wf)) THEN
         CALL fm_pools_give_back_fm_vect(ao_mo_pools, snapshot%wf)
      END IF

      IF (wf_history%store_rho_r) THEN
         CALL qs_rho_get(rho, rho_r=rho_r)
         CPASSERT(ASSOCIATED(rho_r))
         IF (.NOT. ASSOCIATED(snapshot%rho_r)) THEN
            ALLOCATE (snapshot%rho_r(nspins))
            DO ispin = 1, nspins
               NULLIFY (snapshot%rho_r(ispin)%pw)
               CALL pw_pool_create_pw(auxbas_pw_pool, snapshot%rho_r(ispin)%pw, &
                                      in_space=REALSPACE, use_data=REALDATA3D)
            END DO
         END IF
         DO ispin = 1, nspins
            CALL pw_copy(rho_r(ispin)%pw, snapshot%rho_r(ispin)%pw)
         END DO
      ELSE IF (ASSOCIATED(snapshot%rho_r)) THEN
         DO ispin = 1, SIZE(snapshot%rho_r)
            CALL pw_pool_give_back_pw(auxbas_pw_pool, snapshot%rho_r(ispin)%pw)
         END DO
         DEALLOCATE (snapshot%rho_r)
      END IF

      IF (wf_history%store_rho_g) THEN
         CALL qs_rho_get(rho, rho_g=rho_g)
         CPASSERT(ASSOCIATED(rho_g))
         IF (.NOT. ASSOCIATED(snapshot%rho_g)) THEN
            ALLOCATE (snapshot%rho_g(nspins))
            DO ispin = 1, nspins
               NULLIFY (snapshot%rho_g(ispin)%pw)
               CALL pw_pool_create_pw(auxbas_pw_pool, snapshot%rho_g(ispin)%pw, &
                                      in_space=RECIPROCALSPACE, use_data=COMPLEXDATA1D)
            END DO
         END IF
         DO ispin = 1, nspins
            CALL pw_copy(rho_g(ispin)%pw, snapshot%rho_g(ispin)%pw)
         END DO
      ELSE IF (ASSOCIATED(snapshot%rho_g)) THEN
         DO ispin = 1, SIZE(snapshot%rho_g)
            CALL pw_pool_give_back_pw(auxbas_pw_pool, snapshot%rho_g(ispin)%pw)
         END DO
         DEALLOCATE (snapshot%rho_g)
      END IF

      IF (ASSOCIATED(snapshot%rho_ao)) THEN ! the sparsity might be different
         ! (future struct:check)
         CALL dbcsr_deallocate_matrix_set(snapshot%rho_ao)
      END IF
      IF (wf_history%store_rho_ao) THEN
         CALL qs_rho_get(rho, rho_ao=rho_ao)
         CPASSERT(ASSOCIATED(rho_ao))

         CALL dbcsr_allocate_matrix_set(snapshot%rho_ao, nspins)
         DO ispin = 1, nspins
            ALLOCATE (snapshot%rho_ao(ispin)%matrix)
            CALL dbcsr_copy(snapshot%rho_ao(ispin)%matrix, rho_ao(ispin)%matrix)
         END DO
      END IF

      IF (ASSOCIATED(snapshot%rho_ao_kp)) THEN ! the sparsity might be different
         ! (future struct:check)
         CALL dbcsr_deallocate_matrix_set(snapshot%rho_ao_kp)
      END IF
      IF (wf_history%store_rho_ao_kp) THEN
         CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)
         CPASSERT(ASSOCIATED(rho_ao_kp))

         nimg = dft_control%nimages
         CALL dbcsr_allocate_matrix_set(snapshot%rho_ao_kp, nspins, nimg)
         DO ispin = 1, nspins
            DO img = 1, nimg
               ALLOCATE (snapshot%rho_ao_kp(ispin, img)%matrix)
               CALL dbcsr_copy(snapshot%rho_ao_kp(ispin, img)%matrix, &
                               rho_ao_kp(ispin, img)%matrix)
            END DO
         END DO
      END IF

      IF (ASSOCIATED(snapshot%overlap)) THEN ! the sparsity might be different
         ! (future struct:check)
         CALL dbcsr_deallocate_matrix(snapshot%overlap)
      END IF
      IF (wf_history%store_overlap) THEN
         CALL get_qs_env(qs_env, matrix_s=matrix_s)
         CPASSERT(ASSOCIATED(matrix_s))
         CPASSERT(ASSOCIATED(matrix_s(1)%matrix))
         ALLOCATE (snapshot%overlap)
         CALL dbcsr_copy(snapshot%overlap, matrix_s(1)%matrix)
      END IF

      IF (wf_history%store_frozen_density) THEN
         ! do nothing
         ! CALL deallocate_matrix_set(snapshot%rho_frozen%rho_ao)
      END IF

      CALL timestop(handle)

   END SUBROUTINE wfs_update

! **************************************************************************************************
!> \brief ...
!> \param wf_history ...
!> \param interpolation_method_nr the tag of the method used for
!>        the extrapolation of the initial density for the next md step
!>        (see qs_wf_history_types:wfi_*_method_nr)
!> \param extrapolation_order ...
!> \param has_unit_metric ...
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE wfi_create(wf_history, interpolation_method_nr, extrapolation_order, &
                         has_unit_metric)
      TYPE(qs_wf_history_type), POINTER                  :: wf_history
      INTEGER, INTENT(in)                                :: interpolation_method_nr, &
                                                            extrapolation_order
      LOGICAL, INTENT(IN)                                :: has_unit_metric

      CHARACTER(len=*), PARAMETER :: routineN = 'wfi_create', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i

      CALL timeset(routineN, handle)

      ALLOCATE (wf_history)
      last_wfi_id = last_wfi_id + 1
      wf_history%id_nr = last_wfi_id
      wf_history%ref_count = 1
      wf_history%memory_depth = 0
      wf_history%snapshot_count = 0
      wf_history%last_state_index = 1
      wf_history%store_wf = .FALSE.
      wf_history%store_rho_r = .FALSE.
      wf_history%store_rho_g = .FALSE.
      wf_history%store_rho_ao = .FALSE.
      wf_history%store_rho_ao_kp = .FALSE.
      wf_history%store_overlap = .FALSE.
      wf_history%store_frozen_density = .FALSE.
      NULLIFY (wf_history%past_states)

      wf_history%interpolation_method_nr = interpolation_method_nr

      SELECT CASE (wf_history%interpolation_method_nr)
      CASE (wfi_use_guess_method_nr)
         wf_history%memory_depth = 0
      CASE (wfi_use_prev_wf_method_nr)
         wf_history%memory_depth = 0
      CASE (wfi_use_prev_p_method_nr)
         wf_history%memory_depth = 1
         wf_history%store_rho_ao = .TRUE.
      CASE (wfi_use_prev_rho_r_method_nr)
         wf_history%memory_depth = 1
         wf_history%store_rho_ao = .TRUE.
      CASE (wfi_linear_wf_method_nr)
         wf_history%memory_depth = 2
         wf_history%store_wf = .TRUE.
      CASE (wfi_linear_p_method_nr)
         wf_history%memory_depth = 2
         wf_history%store_rho_ao = .TRUE.
      CASE (wfi_linear_ps_method_nr)
         wf_history%memory_depth = 2
         wf_history%store_wf = .TRUE.
         IF (.NOT. has_unit_metric) wf_history%store_overlap = .TRUE.
      CASE (wfi_ps_method_nr)
         CALL cite_reference(VandeVondele2005a)
         wf_history%memory_depth = extrapolation_order + 1
         wf_history%store_wf = .TRUE.
         IF (.NOT. has_unit_metric) wf_history%store_overlap = .TRUE.
      CASE (wfi_frozen_method_nr)
         wf_history%memory_depth = 1
         wf_history%store_frozen_density = .TRUE.
      CASE (wfi_aspc_nr)
         wf_history%memory_depth = extrapolation_order + 2
         wf_history%store_wf = .TRUE.
         IF (.NOT. has_unit_metric) wf_history%store_overlap = .TRUE.
      CASE default
         CALL cp_abort(__LOCATION__, &
                       "Unknown interpolation method: "// &
                       TRIM(ADJUSTL(cp_to_string(interpolation_method_nr))))
         wf_history%interpolation_method_nr = wfi_use_prev_rho_r_method_nr
      END SELECT
      ALLOCATE (wf_history%past_states(wf_history%memory_depth))

      DO i = 1, SIZE(wf_history%past_states)
         NULLIFY (wf_history%past_states(i)%snapshot)
      END DO

      CALL timestop(handle)
   END SUBROUTINE wfi_create

! **************************************************************************************************
!> \brief ...
!> \param wf_history ...
!> \par History
!>      06.2015 created [jhu]
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE wfi_create_for_kp(wf_history)
      TYPE(qs_wf_history_type), POINTER                  :: wf_history

      CHARACTER(len=*), PARAMETER :: routineN = 'wfi_create_for_kp', &
         routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(wf_history))
      IF (wf_history%store_rho_ao) THEN
         wf_history%store_rho_ao_kp = .TRUE.
         wf_history%store_rho_ao = .FALSE.
      END IF
      ! Check for KP compatible methods
      IF (wf_history%store_wf) THEN
         CPABORT("WFN based interpolation method not possible for kpoints.")
      END IF
      IF (wf_history%store_frozen_density) THEN
         CPABORT("Frozen density initialization method not possible for kpoints.")
      END IF
      IF (wf_history%store_overlap) THEN
         CPABORT("Inconsistent interpolation method for kpoints.")
      END IF

   END SUBROUTINE wfi_create_for_kp

! **************************************************************************************************
!> \brief returns a string describing the interpolation method
!> \param method_nr ...
!> \return ...
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
   FUNCTION wfi_get_method_label(method_nr) RESULT(res)
      INTEGER, INTENT(in)                                :: method_nr
      CHARACTER(len=30)                                  :: res

      CHARACTER(len=*), PARAMETER :: routineN = 'wfi_get_method_label', &
         routineP = moduleN//':'//routineN

      res = "unknown"
      SELECT CASE (method_nr)
      CASE (wfi_use_prev_p_method_nr)
         res = "previous_p"
      CASE (wfi_use_prev_wf_method_nr)
         res = "previous_wf"
      CASE (wfi_use_prev_rho_r_method_nr)
         res = "previous_rho_r"
      CASE (wfi_use_guess_method_nr)
         res = "initial_guess"
      CASE (wfi_linear_wf_method_nr)
         res = "mo linear"
      CASE (wfi_linear_p_method_nr)
         res = "P linear"
      CASE (wfi_linear_ps_method_nr)
         res = "PS linear"
      CASE (wfi_ps_method_nr)
         res = "PS Nth order"
      CASE (wfi_frozen_method_nr)
         res = "frozen density approximation"
      CASE (wfi_aspc_nr)
         res = "ASPC"
      CASE default
         CALL cp_abort(__LOCATION__, &
                       "Unknown interpolation method: "// &
                       TRIM(ADJUSTL(cp_to_string(method_nr))))
      END SELECT
   END FUNCTION wfi_get_method_label

! **************************************************************************************************
!> \brief calculates the new starting state for the scf for the next
!>      wf optimization
!> \param wf_history the previous history needed to extrapolate
!> \param qs_env the qs env with the latest result, and that will contain
!>        the new starting state
!> \param dt the time at which to extrapolate (wrt. to the last snapshot)
!> \param extrapolation_method_nr returns the extrapolation method used
!> \param orthogonal_wf ...
!> \par History
!>      02.2003 created [fawzi]
!>      11.2003 Joost VandeVondele : Implemented Nth order PS extrapolation
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE wfi_extrapolate(wf_history, qs_env, dt, extrapolation_method_nr, &
                              orthogonal_wf)
      TYPE(qs_wf_history_type), POINTER                  :: wf_history
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), INTENT(IN)                          :: dt
      INTEGER, INTENT(OUT), OPTIONAL                     :: extrapolation_method_nr
      LOGICAL, INTENT(OUT), OPTIONAL                     :: orthogonal_wf

      CHARACTER(len=*), PARAMETER :: routineN = 'wfi_extrapolate', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: actual_extrapolation_method_nr, handle, &
                                                            i, img, io_unit, ispin, k, n, nmo, &
                                                            nvec, print_level
      LOGICAL                                            :: do_kpoints, my_orthogonal_wf, use_overlap
      REAL(KIND=dp)                                      :: alpha, t0, t1, t2
      TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER     :: ao_mo_fm_pools
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct, matrix_struct_new
      TYPE(cp_fm_type), POINTER                          :: csc, fm_tmp, mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao, rho_frozen_ao
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_kp
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_wf_snapshot_type), POINTER                 :: t0_state, t1_state

      NULLIFY (mos, ao_mo_fm_pools, t0_state, t1_state, mo_coeff, &
               rho, rho_ao, rho_frozen_ao)

      use_overlap = wf_history%store_overlap

      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()
      print_level = logger%iter_info%print_level
      io_unit = cp_print_key_unit_nr(logger, qs_env%input, "DFT%SCF%PRINT%PROGRAM_RUN_INFO", &
                                     extension=".scfLog")

      CPASSERT(ASSOCIATED(wf_history))
      CPASSERT(wf_history%ref_count > 0)
      CPASSERT(ASSOCIATED(qs_env))
      CPASSERT(qs_env%ref_count > 0)
      CALL get_qs_env(qs_env, mos=mos, rho=rho, do_kpoints=do_kpoints)
      CALL mpools_get(qs_env%mpools, ao_mo_fm_pools=ao_mo_fm_pools)
      ! chooses the method for this extrapolation
      IF (wf_history%snapshot_count < 1) THEN
         actual_extrapolation_method_nr = wfi_use_guess_method_nr
      ELSE
         actual_extrapolation_method_nr = wf_history%interpolation_method_nr
      END IF

      SELECT CASE (actual_extrapolation_method_nr)
      CASE (wfi_linear_wf_method_nr)
         IF (wf_history%snapshot_count < 2) THEN
            actual_extrapolation_method_nr = wfi_use_prev_wf_method_nr
         END IF
      CASE (wfi_linear_p_method_nr)
         IF (wf_history%snapshot_count < 2) THEN
            IF (do_kpoints) THEN
               actual_extrapolation_method_nr = wfi_use_guess_method_nr
            ELSE
               actual_extrapolation_method_nr = wfi_use_prev_wf_method_nr
            ENDIF
         END IF
      CASE (wfi_linear_ps_method_nr)
         IF (wf_history%snapshot_count < 2) THEN
            actual_extrapolation_method_nr = wfi_use_prev_wf_method_nr
         END IF
      END SELECT

      IF (PRESENT(extrapolation_method_nr)) &
         extrapolation_method_nr = actual_extrapolation_method_nr
      my_orthogonal_wf = .FALSE.

      SELECT CASE (actual_extrapolation_method_nr)
      CASE (wfi_frozen_method_nr)
         CPASSERT(.NOT. do_kpoints)
         t0_state => wfi_get_snapshot(wf_history, wf_index=1)
         CPASSERT(ASSOCIATED(t0_state%rho_frozen))

         nvec = MIN(wf_history%memory_depth, wf_history%snapshot_count)
         CALL wfi_set_history_variables(qs_env=qs_env, nvec=nvec)

         CALL qs_rho_get(t0_state%rho_frozen, rho_ao=rho_frozen_ao)
         CALL qs_rho_get(rho, rho_ao=rho_ao)
         DO ispin = 1, SIZE(rho_frozen_ao)
            CALL dbcsr_copy(rho_ao(ispin)%matrix, &
                            rho_frozen_ao(ispin)%matrix, &
                            keep_sparsity=.TRUE.)
         END DO
         !FM updating rho_ao directly with t0_state%rho_ao would have the
         !FM wrong matrix structure
         CALL qs_rho_update_rho(rho, qs_env=qs_env)
         CALL qs_ks_did_change(qs_env%ks_env, rho_changed=.TRUE.)

         my_orthogonal_wf = .FALSE.
      CASE (wfi_use_prev_rho_r_method_nr)
         t0_state => wfi_get_snapshot(wf_history, wf_index=1)
         nvec = MIN(wf_history%memory_depth, wf_history%snapshot_count)
         CALL wfi_set_history_variables(qs_env=qs_env, nvec=nvec)
         IF (do_kpoints) THEN
            CPASSERT(ASSOCIATED(t0_state%rho_ao_kp))
            CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)
            DO ispin = 1, SIZE(t0_state%rho_ao_kp, 1)
               DO img = 1, SIZE(t0_state%rho_ao_kp, 2)
                  IF (img > SIZE(rho_ao_kp, 2)) THEN
                     CPWARN("Change in cell neighborlist: might affect quality of initial guess")
                  ELSE
                     CALL dbcsr_copy(rho_ao_kp(ispin, img)%matrix, &
                                     t0_state%rho_ao_kp(ispin, img)%matrix, &
                                     keep_sparsity=.TRUE.)
                  END IF
               END DO
            END DO
         ELSE
            CPASSERT(ASSOCIATED(t0_state%rho_ao))
            CALL qs_rho_get(rho, rho_ao=rho_ao)
            DO ispin = 1, SIZE(t0_state%rho_ao)
               CALL dbcsr_copy(rho_ao(ispin)%matrix, &
                               t0_state%rho_ao(ispin)%matrix, &
                               keep_sparsity=.TRUE.)
            END DO
         END IF
         ! Why is rho_g valid at this point ?
         CALL qs_rho_set(rho, rho_g_valid=.TRUE.)

         ! does nothing
      CASE (wfi_use_prev_wf_method_nr)
         CPASSERT(.NOT. do_kpoints)
         my_orthogonal_wf = .TRUE.
         nvec = MIN(wf_history%memory_depth, wf_history%snapshot_count)
         CALL wfi_set_history_variables(qs_env=qs_env, nvec=nvec)
         CALL qs_rho_get(rho, rho_ao=rho_ao)
         DO ispin = 1, SIZE(mos)
            CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff, &
                            nmo=nmo)
            CALL reorthogonalize_vectors(qs_env, &
                                         v_matrix=mo_coeff, &
                                         n_col=nmo)
            CALL calculate_density_matrix(mo_set=mos(ispin)%mo_set, &
                                          density_matrix=rho_ao(ispin)%matrix)
         END DO
         CALL qs_rho_update_rho(rho, qs_env=qs_env)
         CALL qs_ks_did_change(qs_env%ks_env, rho_changed=.TRUE.)

      CASE (wfi_use_prev_p_method_nr)
         t0_state => wfi_get_snapshot(wf_history, wf_index=1)
         nvec = MIN(wf_history%memory_depth, wf_history%snapshot_count)
         CALL wfi_set_history_variables(qs_env=qs_env, nvec=nvec)
         IF (do_kpoints) THEN
            CPASSERT(ASSOCIATED(t0_state%rho_ao_kp))
            CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)
            DO ispin = 1, SIZE(t0_state%rho_ao_kp, 1)
               DO img = 1, SIZE(t0_state%rho_ao_kp, 2)
                  IF (img > SIZE(rho_ao_kp, 2)) THEN
                     CPWARN("Change in cell neighborlist: might affect quality of initial guess")
                  ELSE
                     CALL dbcsr_copy(rho_ao_kp(ispin, img)%matrix, &
                                     t0_state%rho_ao_kp(ispin, img)%matrix, &
                                     keep_sparsity=.TRUE.)
                  END IF
               END DO
            END DO
         ELSE
            CPASSERT(ASSOCIATED(t0_state%rho_ao))
            CALL qs_rho_get(rho, rho_ao=rho_ao)
            DO ispin = 1, SIZE(t0_state%rho_ao)
               CALL dbcsr_copy(rho_ao(ispin)%matrix, &
                               t0_state%rho_ao(ispin)%matrix, &
                               keep_sparsity=.TRUE.)
            END DO
         END IF
         !FM updating rho_ao directly with t0_state%rho_ao would have the
         !FM wrong matrix structure
         CALL qs_rho_update_rho(rho, qs_env=qs_env)
         CALL qs_ks_did_change(qs_env%ks_env, rho_changed=.TRUE.)

      CASE (wfi_use_guess_method_nr)
         !FM more clean to do it here, but it
         !FM might need to read a file (restart) and thus globenv
         !FM I do not want globenv here, thus done by the caller
         !FM (btw. it also needs the eigensolver, and unless you relocate it
         !FM gives circular dependencies)
         nvec = MIN(wf_history%memory_depth, wf_history%snapshot_count)
         CALL wfi_set_history_variables(qs_env=qs_env, nvec=nvec)
      CASE (wfi_linear_wf_method_nr)
         CPASSERT(.NOT. do_kpoints)
         t0_state => wfi_get_snapshot(wf_history, wf_index=2)
         t1_state => wfi_get_snapshot(wf_history, wf_index=1)
         CPASSERT(ASSOCIATED(t0_state))
         CPASSERT(ASSOCIATED(t1_state))
         CPASSERT(ASSOCIATED(t0_state%wf))
         CPASSERT(ASSOCIATED(t1_state%wf))
         nvec = MIN(wf_history%memory_depth, wf_history%snapshot_count)
         CALL wfi_set_history_variables(qs_env=qs_env, nvec=nvec)

         my_orthogonal_wf = .TRUE.
         t0 = 0.0_dp
         t1 = t1_state%dt
         t2 = t1 + dt
         CALL qs_rho_get(rho, rho_ao=rho_ao)
         DO ispin = 1, SIZE(mos)
            CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff, &
                            nmo=nmo)
            CALL cp_fm_scale_and_add(alpha=0.0_dp, &
                                     matrix_a=mo_coeff, &
                                     matrix_b=t1_state%wf(ispin)%matrix, &
                                     beta=(t2 - t0)/(t1 - t0))
            ! this copy should be unnecessary
            CALL cp_fm_scale_and_add(alpha=1.0_dp, &
                                     matrix_a=mo_coeff, &
                                     beta=(t1 - t2)/(t1 - t0), matrix_b=t0_state%wf(ispin)%matrix)
            CALL reorthogonalize_vectors(qs_env, &
                                         v_matrix=mo_coeff, &
                                         n_col=nmo)
            CALL calculate_density_matrix(mo_set=mos(ispin)%mo_set, &
                                          density_matrix=rho_ao(ispin)%matrix)
         END DO
         CALL qs_rho_update_rho(rho, qs_env=qs_env)

         CALL qs_ks_did_change(qs_env%ks_env, &
                               rho_changed=.TRUE.)
      CASE (wfi_linear_p_method_nr)
         t0_state => wfi_get_snapshot(wf_history, wf_index=2)
         t1_state => wfi_get_snapshot(wf_history, wf_index=1)
         CPASSERT(ASSOCIATED(t0_state))
         CPASSERT(ASSOCIATED(t1_state))
         IF (do_kpoints) THEN
            CPASSERT(ASSOCIATED(t0_state%rho_ao_kp))
            CPASSERT(ASSOCIATED(t1_state%rho_ao_kp))
         ELSE
            CPASSERT(ASSOCIATED(t0_state%rho_ao))
            CPASSERT(ASSOCIATED(t1_state%rho_ao))
         END IF
         nvec = MIN(wf_history%memory_depth, wf_history%snapshot_count)
         CALL wfi_set_history_variables(qs_env=qs_env, nvec=nvec)

         t0 = 0.0_dp
         t1 = t1_state%dt
         t2 = t1 + dt
         IF (do_kpoints) THEN
            CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)
            DO ispin = 1, SIZE(rho_ao_kp, 1)
               DO img = 1, SIZE(rho_ao_kp, 2)
                  IF (img > SIZE(t0_state%rho_ao_kp, 2) .OR. &
                      img > SIZE(t1_state%rho_ao_kp, 2)) THEN
                     CPWARN("Change in cell neighborlist: might affect quality of initial guess")
                  ELSE
                     CALL dbcsr_add(rho_ao_kp(ispin, img)%matrix, t1_state%rho_ao_kp(ispin, img)%matrix, &
                                    alpha_scalar=0.0_dp, beta_scalar=(t2 - t0)/(t1 - t0)) ! this copy should be unnecessary
                     CALL dbcsr_add(rho_ao_kp(ispin, img)%matrix, t0_state%rho_ao_kp(ispin, img)%matrix, &
                                    alpha_scalar=1.0_dp, beta_scalar=(t1 - t2)/(t1 - t0))
                  END IF
               END DO
            END DO
         ELSE
            CALL qs_rho_get(rho, rho_ao=rho_ao)
            DO ispin = 1, SIZE(rho_ao)
               CALL dbcsr_add(rho_ao(ispin)%matrix, t1_state%rho_ao(ispin)%matrix, &
                              alpha_scalar=0.0_dp, beta_scalar=(t2 - t0)/(t1 - t0)) ! this copy should be unnecessary
               CALL dbcsr_add(rho_ao(ispin)%matrix, t0_state%rho_ao(ispin)%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=(t1 - t2)/(t1 - t0))
            END DO
         END IF
         CALL qs_rho_update_rho(rho, qs_env=qs_env)
         CALL qs_ks_did_change(qs_env%ks_env, rho_changed=.TRUE.)
      CASE (wfi_linear_ps_method_nr)
         ! wf not calculated, extract with PSC renormalized?
         ! use wf_linear?
         CPASSERT(.NOT. do_kpoints)
         t0_state => wfi_get_snapshot(wf_history, wf_index=2)
         t1_state => wfi_get_snapshot(wf_history, wf_index=1)
         CPASSERT(ASSOCIATED(t0_state))
         CPASSERT(ASSOCIATED(t1_state))
         CPASSERT(ASSOCIATED(t0_state%wf))
         CPASSERT(ASSOCIATED(t1_state%wf))
         IF (wf_history%store_overlap) THEN
            CPASSERT(ASSOCIATED(t0_state%overlap))
            CPASSERT(ASSOCIATED(t1_state%overlap))
         END IF
         nvec = MIN(wf_history%memory_depth, wf_history%snapshot_count)
         IF (nvec >= wf_history%memory_depth) THEN
            IF ((qs_env%scf_control%max_scf_hist .NE. 0) .AND. (qs_env%scf_control%eps_scf_hist .NE. 0)) THEN
               qs_env%scf_control%max_scf = qs_env%scf_control%max_scf_hist
               qs_env%scf_control%eps_scf = qs_env%scf_control%eps_scf_hist
               qs_env%scf_control%outer_scf%have_scf = .FALSE.
            ELSE IF (qs_env%scf_control%max_scf_hist .NE. 0) THEN
               qs_env%scf_control%max_scf = qs_env%scf_control%max_scf_hist
               qs_env%scf_control%outer_scf%have_scf = .FALSE.
            ELSE IF (qs_env%scf_control%eps_scf_hist .NE. 0) THEN
               qs_env%scf_control%eps_scf = qs_env%scf_control%eps_scf_hist
            END IF
         END IF

         my_orthogonal_wf = .TRUE.
         ! use PS_2=2 PS_1-PS_0
         ! C_2 comes from using PS_2 as a projector acting on C_1
         CALL qs_rho_get(rho, rho_ao=rho_ao)
         DO ispin = 1, SIZE(mos)
            NULLIFY (mo_coeff, matrix_struct, matrix_struct_new, csc)
            CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff)
            CALL cp_fm_get_info(mo_coeff, nrow_global=n, ncol_global=k, &
                                matrix_struct=matrix_struct)
            CALL cp_fm_struct_create(matrix_struct_new, template_fmstruct=matrix_struct, &
                                     nrow_global=k, ncol_global=k)
            CALL cp_fm_create(csc, matrix_struct_new)
            CALL cp_fm_struct_release(matrix_struct_new)

            IF (use_overlap) THEN
               CALL cp_dbcsr_sm_fm_multiply(t0_state%overlap, t1_state%wf(ispin)%matrix, mo_coeff, k)
               CALL cp_gemm('T', 'N', k, k, n, 1.0_dp, t0_state%wf(ispin)%matrix, mo_coeff, 0.0_dp, csc)
            ELSE
               CALL cp_gemm('T', 'N', k, k, n, 1.0_dp, t0_state%wf(ispin)%matrix, &
                            t1_state%wf(ispin)%matrix, 0.0_dp, csc)
            END IF
            CALL cp_gemm('N', 'N', n, k, k, 1.0_dp, t0_state%wf(ispin)%matrix, csc, 0.0_dp, mo_coeff)
            CALL cp_fm_release(csc)
            CALL cp_fm_scale_and_add(-1.0_dp, mo_coeff, 2.0_dp, t1_state%wf(ispin)%matrix)
            CALL reorthogonalize_vectors(qs_env, &
                                         v_matrix=mo_coeff, &
                                         n_col=k)
            CALL calculate_density_matrix(mo_set=mos(ispin)%mo_set, &
                                          density_matrix=rho_ao(ispin)%matrix)
         END DO
         CALL qs_rho_update_rho(rho, qs_env=qs_env)
         CALL qs_ks_did_change(qs_env%ks_env, rho_changed=.TRUE.)

      CASE (wfi_ps_method_nr)
         CPASSERT(.NOT. do_kpoints)
         ! figure out the actual number of vectors to use in the extrapolation:
         nvec = MIN(wf_history%memory_depth, wf_history%snapshot_count)
         CPASSERT(nvec .GT. 0)
         IF (nvec >= wf_history%memory_depth) THEN
            IF ((qs_env%scf_control%max_scf_hist .NE. 0) .AND. (qs_env%scf_control%eps_scf_hist .NE. 0)) THEN
               qs_env%scf_control%max_scf = qs_env%scf_control%max_scf_hist
               qs_env%scf_control%eps_scf = qs_env%scf_control%eps_scf_hist
               qs_env%scf_control%outer_scf%have_scf = .FALSE.
            ELSE IF (qs_env%scf_control%max_scf_hist .NE. 0) THEN
               qs_env%scf_control%max_scf = qs_env%scf_control%max_scf_hist
               qs_env%scf_control%outer_scf%have_scf = .FALSE.
            ELSE IF (qs_env%scf_control%eps_scf_hist .NE. 0) THEN
               qs_env%scf_control%eps_scf = qs_env%scf_control%eps_scf_hist
            END IF
         END IF

         my_orthogonal_wf = .TRUE.
         DO ispin = 1, SIZE(mos)
            NULLIFY (mo_coeff, matrix_struct, matrix_struct_new, csc, fm_tmp)
            CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff)
            CALL cp_fm_get_info(mo_coeff, nrow_global=n, ncol_global=k, &
                                matrix_struct=matrix_struct)
            CALL cp_fm_create(fm_tmp, matrix_struct)
            CALL cp_fm_struct_create(matrix_struct_new, template_fmstruct=matrix_struct, &
                                     nrow_global=k, ncol_global=k)
            CALL cp_fm_create(csc, matrix_struct_new)
            CALL cp_fm_struct_release(matrix_struct_new)
            ! first the most recent
            t1_state => wfi_get_snapshot(wf_history, wf_index=1)
            CALL cp_fm_to_fm(t1_state%wf(ispin)%matrix, mo_coeff)
            alpha = nvec
            CALL cp_fm_scale(alpha, mo_coeff)
            CALL qs_rho_get(rho, rho_ao=rho_ao)
            DO i = 2, nvec
               t0_state => wfi_get_snapshot(wf_history, wf_index=i)
               IF (use_overlap) THEN
                  CALL cp_dbcsr_sm_fm_multiply(t0_state%overlap, t1_state%wf(ispin)%matrix, fm_tmp, k)
                  CALL cp_gemm('T', 'N', k, k, n, 1.0_dp, t0_state%wf(ispin)%matrix, fm_tmp, 0.0_dp, csc)
               ELSE
                  CALL cp_gemm('T', 'N', k, k, n, 1.0_dp, t0_state%wf(ispin)%matrix, &
                               t1_state%wf(ispin)%matrix, 0.0_dp, csc)
               END IF
               CALL cp_gemm('N', 'N', n, k, k, 1.0_dp, t0_state%wf(ispin)%matrix, csc, 0.0_dp, fm_tmp)
               alpha = -1.0_dp*alpha*REAL(nvec - i + 1, dp)/REAL(i, dp)
               CALL cp_fm_scale_and_add(1.0_dp, mo_coeff, alpha, fm_tmp)
            ENDDO

            CALL cp_fm_release(csc)
            CALL cp_fm_release(fm_tmp)
            CALL reorthogonalize_vectors(qs_env, &
                                         v_matrix=mo_coeff, &
                                         n_col=k)
            CALL calculate_density_matrix(mo_set=mos(ispin)%mo_set, &
                                          density_matrix=rho_ao(ispin)%matrix)
         END DO
         CALL qs_rho_update_rho(rho, qs_env=qs_env)
         CALL qs_ks_did_change(qs_env%ks_env, rho_changed=.TRUE.)

      CASE (wfi_aspc_nr)
         CPASSERT(.NOT. do_kpoints)
         CALL cite_reference(Kolafa2004)
         CALL cite_reference(Kuhne2007)
         ! figure out the actual number of vectors to use in the extrapolation:
         nvec = MIN(wf_history%memory_depth, wf_history%snapshot_count)
         CPASSERT(nvec .GT. 0)
         IF (nvec >= wf_history%memory_depth) THEN
            IF ((qs_env%scf_control%max_scf_hist .NE. 0) .AND. &
                (qs_env%scf_control%eps_scf_hist .NE. 0)) THEN
               qs_env%scf_control%max_scf = qs_env%scf_control%max_scf_hist
               qs_env%scf_control%eps_scf = qs_env%scf_control%eps_scf_hist
               qs_env%scf_control%outer_scf%have_scf = .FALSE.
            ELSE IF (qs_env%scf_control%max_scf_hist .NE. 0) THEN
               qs_env%scf_control%max_scf = qs_env%scf_control%max_scf_hist
               qs_env%scf_control%outer_scf%have_scf = .FALSE.
            ELSE IF (qs_env%scf_control%eps_scf_hist .NE. 0) THEN
               qs_env%scf_control%eps_scf = qs_env%scf_control%eps_scf_hist
            END IF
         END IF

         my_orthogonal_wf = .TRUE.
         CALL qs_rho_get(rho, rho_ao=rho_ao)
         DO ispin = 1, SIZE(mos)
            NULLIFY (mo_coeff, matrix_struct, matrix_struct_new, csc, fm_tmp)
            CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff)
            CALL cp_fm_get_info(mo_coeff, &
                                nrow_global=n, &
                                ncol_global=k, &
                                matrix_struct=matrix_struct)
            CALL cp_fm_create(fm_tmp, matrix_struct)
            CALL cp_fm_struct_create(matrix_struct_new, &
                                     template_fmstruct=matrix_struct, &
                                     nrow_global=k, &
                                     ncol_global=k)
            CALL cp_fm_create(csc, matrix_struct_new)
            CALL cp_fm_struct_release(matrix_struct_new)
            ! first the most recent
            t1_state => wfi_get_snapshot(wf_history, &
                                         wf_index=1)
            CALL cp_fm_to_fm(t1_state%wf(ispin)%matrix, mo_coeff)
            alpha = REAL(4*nvec - 2, KIND=dp)/REAL(nvec + 1, KIND=dp)
            IF ((io_unit > 0) .AND. (print_level > low_print_level)) THEN
               WRITE (UNIT=io_unit, FMT="(/,T2,A,/,/,T3,A,I0,/,/,T3,A2,I0,A4,F10.6)") &
                  "Parameters for the always stable predictor-corrector (ASPC) method:", &
                  "ASPC order: ", MAX(nvec - 2, 0), &
                  "B(", 1, ") = ", alpha
            END IF
            CALL cp_fm_scale(alpha, mo_coeff)

            DO i = 2, nvec
               t0_state => wfi_get_snapshot(wf_history, wf_index=i)
               IF (use_overlap) THEN
                  CALL cp_dbcsr_sm_fm_multiply(t0_state%overlap, t1_state%wf(ispin)%matrix, fm_tmp, k)
                  CALL cp_gemm('T', 'N', k, k, n, 1.0_dp, t0_state%wf(ispin)%matrix, fm_tmp, 0.0_dp, csc)
               ELSE
                  CALL cp_gemm('T', 'N', k, k, n, 1.0_dp, t0_state%wf(ispin)%matrix, &
                               t1_state%wf(ispin)%matrix, 0.0_dp, csc)
               END IF
               CALL cp_gemm('N', 'N', n, k, k, 1.0_dp, t0_state%wf(ispin)%matrix, csc, 0.0_dp, fm_tmp)
               alpha = (-1.0_dp)**(i + 1)*REAL(i, KIND=dp)* &
                       binomial(2*nvec, nvec - i)/binomial(2*nvec - 2, nvec - 1)
               IF ((io_unit > 0) .AND. (print_level > low_print_level)) THEN
                  WRITE (UNIT=io_unit, FMT="(T3,A2,I0,A4,F10.6)") &
                     "B(", i, ") = ", alpha
               END IF
               CALL cp_fm_scale_and_add(1.0_dp, mo_coeff, alpha, fm_tmp)
            END DO
            CALL cp_fm_release(csc)
            CALL cp_fm_release(fm_tmp)
            CALL reorthogonalize_vectors(qs_env, &
                                         v_matrix=mo_coeff, &
                                         n_col=k)
            CALL calculate_density_matrix(mo_set=mos(ispin)%mo_set, &
                                          density_matrix=rho_ao(ispin)%matrix)
         END DO
         CALL qs_rho_update_rho(rho, qs_env=qs_env)
         CALL qs_ks_did_change(qs_env%ks_env, rho_changed=.TRUE.)

      CASE default
         CALL cp_abort(__LOCATION__, &
                       "Unknown interpolation method: "// &
                       TRIM(ADJUSTL(cp_to_string(wf_history%interpolation_method_nr))))
      END SELECT
      IF (PRESENT(orthogonal_wf)) orthogonal_wf = my_orthogonal_wf
      CALL cp_print_key_finished_output(io_unit, logger, qs_env%input, &
                                        "DFT%SCF%PRINT%PROGRAM_RUN_INFO")
      CALL timestop(handle)
   END SUBROUTINE wfi_extrapolate

! **************************************************************************************************
!> \brief Decides if scf control variables has to changed due
!>      to using a WF extrapolation.
!> \param qs_env The QS environment
!> \param nvec ...
!> \par History
!>      11.2006 created [TdK]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! **************************************************************************************************
   SUBROUTINE wfi_set_history_variables(qs_env, nvec)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: nvec

      CHARACTER(len=*), PARAMETER :: routineN = 'wfi_set_history_variables', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(qs_env))
      CPASSERT(qs_env%ref_count > 0)

      IF (nvec >= qs_env%wf_history%memory_depth) THEN
         IF ((qs_env%scf_control%max_scf_hist .NE. 0) .AND. (qs_env%scf_control%eps_scf_hist .NE. 0)) THEN
            qs_env%scf_control%max_scf = qs_env%scf_control%max_scf_hist
            qs_env%scf_control%eps_scf = qs_env%scf_control%eps_scf_hist
            qs_env%scf_control%outer_scf%have_scf = .FALSE.
         ELSE IF (qs_env%scf_control%max_scf_hist .NE. 0) THEN
            qs_env%scf_control%max_scf = qs_env%scf_control%max_scf_hist
            qs_env%scf_control%outer_scf%have_scf = .FALSE.
         ELSE IF (qs_env%scf_control%eps_scf_hist .NE. 0) THEN
            qs_env%scf_control%eps_scf = qs_env%scf_control%eps_scf_hist
            qs_env%scf_control%outer_scf%eps_scf = qs_env%scf_control%eps_scf_hist
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE wfi_set_history_variables

! **************************************************************************************************
!> \brief updates the snapshot buffer, taking a new snapshot
!> \param wf_history the history buffer to update
!> \param qs_env the qs_env we get the info from
!> \param dt ...
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE wfi_update(wf_history, qs_env, dt)
      TYPE(qs_wf_history_type), POINTER                  :: wf_history
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), INTENT(in)                          :: dt

      CHARACTER(len=*), PARAMETER :: routineN = 'wfi_update', routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(wf_history))
      CPASSERT(wf_history%ref_count > 0)
      CPASSERT(ASSOCIATED(qs_env))
      CPASSERT(qs_env%ref_count > 0)

      wf_history%snapshot_count = wf_history%snapshot_count + 1
      IF (wf_history%memory_depth > 0) THEN
         wf_history%last_state_index = MODULO(wf_history%snapshot_count, &
                                              wf_history%memory_depth) + 1
         CALL wfs_update(snapshot=wf_history%past_states &
                         (wf_history%last_state_index)%snapshot, wf_history=wf_history, &
                         qs_env=qs_env, dt=dt)
      END IF
   END SUBROUTINE wfi_update

! **************************************************************************************************
!> \brief reorthogonalizes the mos
!> \param qs_env the qs_env in which to orthogonalize
!> \param v_matrix the vectors to orthogonalize
!> \param n_col number of column of v to orthogonalize
!> \par History
!>      04.2003 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE reorthogonalize_vectors(qs_env, v_matrix, n_col)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), POINTER                          :: v_matrix
      INTEGER, INTENT(in), OPTIONAL                      :: n_col

      CHARACTER(len=*), PARAMETER :: routineN = 'reorthogonalize_vectors', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, my_n_col
      LOGICAL                                            :: has_unit_metric, &
                                                            ortho_contains_cholesky, &
                                                            smearing_is_used
      TYPE(cp_fm_pool_type), POINTER                     :: maxao_maxmo_fm_pool
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_matrix_pools_type), POINTER                :: mpools
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control

      NULLIFY (scf_env, scf_control, maxao_maxmo_fm_pool, matrix_s, mpools, dft_control)
      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(qs_env))
      CPASSERT(ASSOCIATED(v_matrix))

      CALL cp_fm_get_info(v_matrix, ncol_global=my_n_col)
      IF (PRESENT(n_col)) my_n_col = n_col
      CALL get_qs_env(qs_env, mpools=mpools, &
                      scf_env=scf_env, &
                      scf_control=scf_control, &
                      matrix_s=matrix_s, &
                      dft_control=dft_control)
      CALL mpools_get(mpools, maxao_maxmo_fm_pool=maxao_maxmo_fm_pool)
      IF (ASSOCIATED(scf_env)) THEN
         ortho_contains_cholesky = (scf_env%method /= ot_method_nr) .AND. &
                                   (scf_env%cholesky_method > 0) .AND. &
                                   ASSOCIATED(scf_env%ortho)
      ELSE
         ortho_contains_cholesky = .FALSE.
      END IF

      CALL get_qs_env(qs_env, has_unit_metric=has_unit_metric)
      smearing_is_used = .FALSE.
      IF (dft_control%smear) THEN
         smearing_is_used = .TRUE.
      END IF

      IF (has_unit_metric) THEN
         CALL make_basis_simple(v_matrix, my_n_col)
      ELSE IF (smearing_is_used) THEN
         CALL make_basis_lowdin(vmatrix=v_matrix, ncol=my_n_col, &
                                matrix_s=matrix_s(1)%matrix)
      ELSE IF (ortho_contains_cholesky) THEN
         CALL make_basis_cholesky(vmatrix=v_matrix, ncol=my_n_col, &
                                  ortho=scf_env%ortho)
      ELSE
         CALL make_basis_sm(v_matrix, my_n_col, matrix_s(1)%matrix)
      END IF
      CALL timestop(handle)
   END SUBROUTINE reorthogonalize_vectors

! **************************************************************************************************
!> \brief purges wf_history retaining only the latest snapshot
!> \param qs_env the qs env with the latest result, and that will contain
!>        the purged wf_history
!> \par History
!>      05.2016 created [Nico Holmberg]
!> \author Nico Holmberg
! **************************************************************************************************
   SUBROUTINE wfi_purge_history(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'wfi_purge_history', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, io_unit, print_level
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_wf_history_type), POINTER                  :: wf_history

      NULLIFY (dft_control, wf_history)

      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()
      print_level = logger%iter_info%print_level
      io_unit = cp_print_key_unit_nr(logger, qs_env%input, "DFT%SCF%PRINT%PROGRAM_RUN_INFO", &
                                     extension=".scfLog")

      CPASSERT(ASSOCIATED(qs_env))
      CPASSERT(qs_env%ref_count > 0)
      CPASSERT(ASSOCIATED(qs_env%wf_history))
      CPASSERT(qs_env%wf_history%ref_count > 0)
      CALL get_qs_env(qs_env, dft_control=dft_control)

      SELECT CASE (qs_env%wf_history%interpolation_method_nr)
      CASE (wfi_use_guess_method_nr, wfi_use_prev_wf_method_nr, &
            wfi_use_prev_p_method_nr, wfi_use_prev_rho_r_method_nr, &
            wfi_frozen_method_nr)
         ! do nothing
      CASE (wfi_linear_wf_method_nr, wfi_linear_p_method_nr, &
            wfi_linear_ps_method_nr, wfi_ps_method_nr, &
            wfi_aspc_nr)
         IF (qs_env%wf_history%snapshot_count .GE. 2) THEN
            IF (debug_this_module .AND. io_unit > 0) &
               WRITE (io_unit, FMT="(T2,A)") "QS| Purging WFN history"
            CALL wfi_create(wf_history, interpolation_method_nr= &
                            dft_control%qs_control%wf_interpolation_method_nr, &
                            extrapolation_order=dft_control%qs_control%wf_extrapolation_order, &
                            has_unit_metric=qs_env%has_unit_metric)
            CALL set_qs_env(qs_env=qs_env, &
                            wf_history=wf_history)
            CALL wfi_release(wf_history)
            CALL wfi_update(qs_env%wf_history, qs_env=qs_env, dt=1.0_dp)
         END IF
      CASE DEFAULT
         CPABORT("Unknown extrapolation method.")
      END SELECT
      CALL timestop(handle)

   END SUBROUTINE wfi_purge_history

END MODULE qs_wf_history_methods
