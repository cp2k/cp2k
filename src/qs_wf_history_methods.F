!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
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
! *****************************************************************************
MODULE qs_wf_history_methods

  USE bibliography,                    ONLY: Kolafa2004,&
                                             VandeVondele2005a,&
                                             cite_reference
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_add,&
                                             cp_dbcsr_allocate_matrix_set,&
                                             cp_dbcsr_copy,&
                                             cp_dbcsr_deallocate_matrix,&
                                             cp_dbcsr_deallocate_matrix_set,&
                                             cp_dbcsr_init,&
                                             cp_dbcsr_p_type,&
                                             cp_dbcsr_set
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply
  USE cp_fm_basic_linalg,              ONLY: cp_fm_scale,&
                                             cp_fm_scale_and_add
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
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
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE f77_blas
  USE input_constants,                 ONLY: &
       wfi_aspc_nr, wfi_frozen_method_nr, wfi_linear_p_method_nr, &
       wfi_linear_ps_method_nr, wfi_linear_wf_method_nr, wfi_ps_method_nr, &
       wfi_use_guess_method_nr, wfi_use_prev_p_method_nr, &
       wfi_use_prev_rho_r_method_nr, wfi_use_prev_wf_method_nr
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
                                             qs_env_reorthogonalize_vectors,&
                                             qs_environment_type
  USE qs_ks_methods,                   ONLY: qs_ks_did_change
  USE qs_matrix_pools,                 ONLY: mpools_get
  USE qs_mo_methods,                   ONLY: calculate_density_matrix
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_rho_methods,                  ONLY: duplicate_rho_type,&
                                             qs_rho_update_rho
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_type
  USE qs_wf_history_types,             ONLY: qs_wf_history_type,&
                                             qs_wf_snapshot_p_type,&
                                             qs_wf_snapshot_type,&
                                             wfi_get_snapshot,&
                                             wfs_release
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_wf_history_methods'
  INTEGER, SAVE, PRIVATE :: last_wfs_id=0, last_wfi_id=0

  PUBLIC :: wfs_create, wfi_create, wfi_update, &
            wfi_extrapolate, wfi_get_method_label, &
            wfs_duplicate_snapshot, &
            wfi_change_memory_depth

CONTAINS

! *****************************************************************************
!> \brief allocates and initialize a wavefunction snapshot
!> \param snapshot the snapshot to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.2003 created [fawzi]
!>      02.2005 added wf_mol [MI]
!> \author fawzi
! *****************************************************************************
SUBROUTINE wfs_create(snapshot, error)
    TYPE(qs_wf_snapshot_type), POINTER       :: snapshot
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'wfs_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.

  ALLOCATE(snapshot, stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     last_wfs_id=last_wfs_id+1
     snapshot%id_nr=last_wfs_id
     NULLIFY(snapshot%wf, snapshot%rho_r, &
             snapshot%rho_g, snapshot%rho_ao,&
             snapshot%overlap, snapshot%rho_frozen)
     snapshot%dt=1.0_dp
     snapshot%ref_count=1
  END IF
END SUBROUTINE wfs_create

! *****************************************************************************
!> \brief duplicates a single qs_wf_snapshot_type
!> \param input_snapshot the snapshot to be duplicated
!> \param output_snapshot the out coming duplicate
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10.2005 created [TdK]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  SUBROUTINE wfs_duplicate_snapshot(input_snapshot, output_snapshot, qs_env, error)
    TYPE(qs_wf_snapshot_type), POINTER       :: input_snapshot, &
                                                output_snapshot
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'wfs_duplicate_snapshot', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, nspins, stat
    LOGICAL                                  :: failure
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool

    CALL timeset(routineN,handle)
    failure = .FALSE.

    NULLIFY(auxbas_pw_pool, pw_env, dft_control)

    CPPrecondition(ASSOCIATED(input_snapshot), cp_failure_level, routineP, error, failure)
    CPPrecondition(input_snapshot%ref_count>0, cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(qs_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(qs_env%ref_count>0, cp_failure_level, routineP, error, failure)

    IF (.NOT. failure) THEN
      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env, dft_control=dft_control, error=error)
      CALL pw_env_get(pw_env=pw_env, auxbas_pw_pool=auxbas_pw_pool, error=error)

      nspins = dft_control%nspins

      IF (ASSOCIATED(output_snapshot)) THEN
        CALL wfs_release(input_snapshot, error=error)
      ELSE
        CALL wfs_create(output_snapshot, error=error)
      END IF

      ! wf
      IF (ASSOCIATED(input_snapshot%wf)) THEN
        DO i = 1,SIZE(input_snapshot%wf)
          ALLOCATE(output_snapshot%wf(i), stat=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          CALL cp_fm_create(matrix=output_snapshot%wf(i)%matrix, &
                            matrix_struct=input_snapshot%wf(i)%matrix%matrix_struct, &
                            error=error)
          CALL cp_fm_to_fm(source=input_snapshot%wf(i)%matrix, &
                           destination=output_snapshot%wf(i)%matrix, &
                           error=error)
        END DO
      END IF

      ! rho_r
      IF (ASSOCIATED(input_snapshot%rho_r)) THEN
        DO i = 1,SIZE(input_snapshot%rho_r)
          ALLOCATE(output_snapshot%rho_r(nspins), stat=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          IF (.NOT. failure) THEN
            CALL pw_pool_create_pw(auxbas_pw_pool, output_snapshot%rho_r(i)%pw, &
                                   use_data=REALDATA3D, in_space=REALSPACE, &
                                   error=error)
            output_snapshot%rho_r(i)%pw%cr3d(:,:,:) = input_snapshot%rho_r(i)%pw%cr3d(:,:,:)
          END IF
        END DO
      END IF

      ! rho_g
      IF (ASSOCIATED(input_snapshot%rho_g)) THEN
        DO i = 1,SIZE(input_snapshot%rho_g)
          ALLOCATE(output_snapshot%rho_g(nspins), stat=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          IF (.NOT. failure) THEN
            CALL pw_pool_create_pw(auxbas_pw_pool, output_snapshot%rho_g(i)%pw, &
                                   use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE, &
                                   error=error)
            output_snapshot%rho_g(i)%pw%cc(:) = input_snapshot%rho_g(i)%pw%cc(:)
          END IF
        END DO
      END IF

      ! rho_ao
      IF (ASSOCIATED(input_snapshot%rho_ao)) THEN
        CALL cp_dbcsr_allocate_matrix_set(output_snapshot%rho_ao, SIZE(input_snapshot%rho_ao),error=error)
        DO i = 1,SIZE(input_snapshot%rho_ao)
           CALL cp_dbcsr_copy(output_snapshot%rho_ao(i)%matrix,input_snapshot%rho_ao(i)%matrix, &
                name="myDensityMatrix_for_Spin_"//TRIM(ADJUSTL(cp_to_string(i))),error=error)
           CALL cp_dbcsr_set(output_snapshot%rho_ao(i)%matrix,0.0_dp,error=error)
        END DO
      END IF

      ! overlap
      IF (ASSOCIATED(input_snapshot%overlap)) THEN
         CALL cp_dbcsr_copy(output_snapshot%overlap,input_snapshot%overlap,&
              name="myOverlapMatrix",error=error)
         CALL cp_dbcsr_set(output_snapshot%overlap,0.0_dp,error=error)
         !CALL replicate_matrix(source=input_snapshot%overlap, &
         !                     TARGET=output_snapshot%overlap, &
         !                     target_name="myOverlapMatrix", &
         !                     allocate_blocks=.TRUE.,error=error)
      END IF

      ! rho_frozen
      IF (ASSOCIATED(input_snapshot%rho_frozen)) THEN
        CALL duplicate_rho_type(rho_input=input_snapshot%rho_frozen, &
                                rho_output=output_snapshot%rho_frozen, &
                                qs_env=qs_env, error=error)
      END IF

    END IF

    CALL timestop(handle)

  END SUBROUTINE wfs_duplicate_snapshot

! *****************************************************************************
!> \brief updates the given snapshot
!> \param snapshot the snapshot to be updated
!> \param wf_history the history
!> \param qs_env the qs_env that should be snapshotted
!> \param dt the time of the snapshot (wrt. to the previous snapshot)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.2003 created [fawzi]
!>      02.2005 added kg_fm_mol_set for KG_GPW [MI]
!> \author fawzi
! *****************************************************************************
SUBROUTINE wfs_update(snapshot,wf_history,qs_env,dt,error)
    TYPE(qs_wf_snapshot_type), POINTER       :: snapshot
    TYPE(qs_wf_history_type), POINTER        :: wf_history
    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(KIND=dp), INTENT(in), OPTIONAL      :: dt
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'wfs_update', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nspins, stat
    LOGICAL                                  :: failure
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: rho_ao, s_array
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: ao_mo_pools
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_g, rho_r
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_rho_type), POINTER               :: rho

  CALL timeset(routineN,handle)

  failure=.FALSE.
  NULLIFY(pw_env, auxbas_pw_pool, ao_mo_pools, dft_control, mos, mo_coeff,&
       rho, rho_r,rho_g,rho_ao, s_array)
  CALL get_qs_env(qs_env, pw_env=pw_env,&
       dft_control=dft_control,rho=rho,matrix_s=s_array,&
       error=error)
  CALL mpools_get(qs_env%mpools, ao_mo_fm_pools=ao_mo_pools, &
       error=error)
  CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, error=error)

  CPPrecondition(ASSOCIATED(wf_history),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(dft_control),cp_failure_level,routineP,error,failure)
  IF (.not.ASSOCIATED(snapshot)) THEN
     CALL wfs_create(snapshot,error=error)
     CALL cp_error_check(error,failure)
  END IF
  IF (.not.failure) THEN
     CPPrecondition(wf_history%ref_count>0,cp_failure_level,routineP,error,failure)
     CPPrecondition(snapshot%ref_count>0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT. failure) THEN
     nspins=dft_control%nspins
     snapshot%dt=1.0_dp
     IF (PRESENT(dt)) snapshot%dt=dt
     IF (wf_history%store_wf) THEN
        CALL get_qs_env(qs_env,mos=mos,error=error)
        IF (.NOT.ASSOCIATED(snapshot%wf)) THEN
           CALL fm_pools_create_fm_vect(ao_mo_pools,snapshot%wf,&
                name="ws_snap"//TRIM(ADJUSTL(cp_to_string(snapshot%id_nr)))//&
                "ws",error=error)
           CPPostcondition(nspins==SIZE(snapshot%wf),cp_failure_level,routineP,error,failure)
        END IF
        DO ispin=1,nspins
           CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff)
           CALL cp_fm_to_fm(mo_coeff,snapshot%wf(ispin)%matrix,error=error)
        END DO
     ELSE IF (ASSOCIATED(snapshot%wf)) THEN
        CALL fm_pools_give_back_fm_vect(ao_mo_pools,snapshot%wf,&
             error=error)
     END IF

     IF (wf_history%store_rho_r) THEN
        CALL qs_rho_get(rho, rho_r=rho_r,error=error)
        CPPrecondition(ASSOCIATED(rho_r),cp_failure_level,routineP,error,failure)
        IF (.NOT.ASSOCIATED(snapshot%rho_r)) THEN
           ALLOCATE(snapshot%rho_r(nspins),stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           DO ispin=1,nspins
              NULLIFY(snapshot%rho_r(ispin)%pw)
              CALL pw_pool_create_pw(auxbas_pw_pool,snapshot%rho_r(ispin)%pw,&
                   in_space=REALSPACE, use_data=REALDATA3D,error=error)
           END DO
        END IF
        DO ispin=1,nspins
           CALL pw_copy(rho_r(ispin)%pw,snapshot%rho_r(ispin)%pw,error=error)
        END DO
     ELSE IF (ASSOCIATED(snapshot%rho_r)) THEN
        DO ispin=1,SIZE(snapshot%rho_r)
           CALL pw_pool_give_back_pw(auxbas_pw_pool,snapshot%rho_r(ispin)%pw,&
                error=error)
        END DO
        DEALLOCATE(snapshot%rho_r,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF

     IF (wf_history%store_rho_g) THEN
        CALL qs_rho_get(rho, rho_g=rho_g,error=error)
        CPPrecondition(ASSOCIATED(rho_g),cp_failure_level,routineP,error,failure)
        IF (.NOT.ASSOCIATED(snapshot%rho_g)) THEN
           ALLOCATE(snapshot%rho_g(nspins),stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           DO ispin=1,nspins
              NULLIFY(snapshot%rho_g(ispin)%pw)
              CALL pw_pool_create_pw(auxbas_pw_pool,snapshot%rho_g(ispin)%pw,&
                   in_space=RECIPROCALSPACE,use_data=COMPLEXDATA1D,&
                   error=error)
           END DO
        END IF
        DO ispin=1,nspins
           CALL pw_copy(rho_g(ispin)%pw,snapshot%rho_g(ispin)%pw,error=error)
        END DO
     ELSE IF (ASSOCIATED(snapshot%rho_g)) THEN
        DO ispin=1,SIZE(snapshot%rho_g)
           CALL pw_pool_give_back_pw(auxbas_pw_pool,snapshot%rho_g(ispin)%pw,&
                error=error)
        END DO
        DEALLOCATE(snapshot%rho_g,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF

     IF (ASSOCIATED(snapshot%rho_ao)) THEN ! the sparsity might be different
        ! (future struct:check)
        CALL cp_dbcsr_deallocate_matrix_set(snapshot%rho_ao,error=error)
     END IF
     IF (wf_history%store_rho_ao) THEN
        CALL qs_rho_get(rho, rho_ao=rho_ao,error=error)
        CPPrecondition(ASSOCIATED(rho_ao),cp_failure_level,routineP,error,failure)

        CALL cp_dbcsr_allocate_matrix_set(snapshot%rho_ao,nspins,error=error)
        DO ispin=1,nspins
           ALLOCATE(snapshot%rho_ao(ispin)%matrix)
           CALL cp_dbcsr_init(snapshot%rho_ao(ispin)%matrix,error=error)
           CALL cp_dbcsr_copy(snapshot%rho_ao(ispin)%matrix,rho_ao(ispin)%matrix,&
                error=error)
        END DO
     END IF

     IF (ASSOCIATED(snapshot%overlap)) THEN ! the sparsity might be different
        ! (future struct:check)
        CALL cp_dbcsr_deallocate_matrix(snapshot%overlap,error=error)
     END IF
     IF (wf_history%store_overlap) THEN
        CPPrecondition(ASSOCIATED(s_array),cp_failure_level,routineP,error,failure)
        CPPrecondition(ASSOCIATED(s_array(1)%matrix),cp_failure_level,routineP,error,failure)
        ALLOCATE(snapshot%overlap)
        CALL cp_dbcsr_init(snapshot%overlap,error=error)
        CALL cp_dbcsr_copy(snapshot%overlap,s_array(1)%matrix,error=error)
     END IF

     IF (wf_history%store_frozen_density) THEN
       ! do nothing
       ! CALL deallocate_matrix_set(snapshot%rho_frozen%rho_ao)
     END IF

  END IF
  CALL timestop(handle)

END SUBROUTINE wfs_update

! *****************************************************************************
!> \param interpolation_method_nr the tag of the method used for
!>        the extrapolation of the intial density for the next md step
!>        (see qs_wf_history_types:wfi_*_method_nr)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
SUBROUTINE wfi_create(wf_history, interpolation_method_nr, extrapolation_order, has_unit_metric, error)
    TYPE(qs_wf_history_type), POINTER        :: wf_history
    INTEGER, INTENT(in)                      :: interpolation_method_nr
    INTEGER, INTENT(in), OPTIONAL            :: extrapolation_order
    LOGICAL, INTENT(IN)                      :: has_unit_metric
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'wfi_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, &
                                                my_extrapolation_order, stat
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)

  my_extrapolation_order = 1
  IF (PRESENT(extrapolation_order)) my_extrapolation_order=extrapolation_order

  failure=.FALSE.

  ALLOCATE(wf_history, stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     last_wfi_id=last_wfi_id+1
     wf_history%id_nr=last_wfi_id
     wf_history%ref_count=1
     wf_history%memory_depth=0
     wf_history%snapshot_count=0
     wf_history%last_state_index=1
     wf_history%store_wf=.FALSE.
     wf_history%store_rho_r=.FALSE.
     wf_history%store_rho_g=.FALSE.
     wf_history%store_rho_ao=.FALSE.
     wf_history%store_overlap=.FALSE.
     wf_history%store_frozen_density=.FALSE.
     NULLIFY(wf_history%past_states)

     wf_history%interpolation_method_nr=interpolation_method_nr

     SELECT CASE(wf_history%interpolation_method_nr)
     CASE(wfi_use_guess_method_nr)
        wf_history%memory_depth=0
     CASE(wfi_use_prev_wf_method_nr)
        wf_history%memory_depth=0
     CASE(wfi_use_prev_p_method_nr)
        wf_history%memory_depth=1
        wf_history%store_rho_ao=.TRUE.
     CASE(wfi_use_prev_rho_r_method_nr)
        wf_history%memory_depth=1
        wf_history%store_rho_ao=.TRUE.
     CASE(wfi_linear_wf_method_nr)
        wf_history%memory_depth=2
        wf_history%store_wf=.TRUE.
     CASE(wfi_linear_p_method_nr)
        wf_history%memory_depth=2
        wf_history%store_rho_ao=.TRUE.
     CASE(wfi_linear_ps_method_nr)
        wf_history%memory_depth=2
        wf_history%store_wf=.TRUE.
        IF(.NOT.has_unit_metric) wf_history%store_overlap=.TRUE.
     CASE(wfi_ps_method_nr)
        CALL cite_reference(VandeVondele2005a)
        wf_history%memory_depth=my_extrapolation_order+1
        wf_history%store_wf=.TRUE.
        IF(.NOT.has_unit_metric) wf_history%store_overlap=.TRUE.
     CASE(wfi_frozen_method_nr)
        wf_history%memory_depth=1
        wf_history%store_frozen_density=.TRUE.
     CASE (wfi_aspc_nr)
       wf_history%memory_depth = my_extrapolation_order + 2
       wf_history%store_wf = .TRUE.
       IF(.NOT.has_unit_metric) wf_history%store_overlap = .TRUE.
     CASE default
        CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
             routineP,"Unknown interpolation method: "//&
             TRIM(ADJUSTL(cp_to_string(interpolation_method_nr)))//" in "//&
CPSourceFileRef,&
             error, failure)
        wf_history%interpolation_method_nr=wfi_use_prev_rho_r_method_nr
     END SELECT
     ALLOCATE(wf_history%past_states(wf_history%memory_depth),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     IF (.NOT. failure) THEN
        DO i=1,SIZE(wf_history%past_states)
           NULLIFY(wf_history%past_states(i)%snapshot)
        END DO
     END IF
  END IF

  CALL timestop(handle)
END SUBROUTINE wfi_create

! *****************************************************************************
!> \brief returns a string describing the interpolation method
!> \param wf_history the wf history object you want information about
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
FUNCTION wfi_get_method_label(method_nr,error) RESULT(res)
    INTEGER, INTENT(in)                      :: method_nr
    TYPE(cp_error_type), INTENT(inout)       :: error
    CHARACTER(len=30)                        :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'wfi_get_method_label', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  res="unknown"
  IF (.NOT. failure) THEN
     SELECT CASE(method_nr)
     CASE(wfi_use_prev_p_method_nr)
        res="previous_p"
     CASE(wfi_use_prev_wf_method_nr)
        res="previous_wf"
     CASE(wfi_use_prev_rho_r_method_nr)
        res="previous_rho_r"
     CASE(wfi_use_guess_method_nr)
        res="initial_guess"
     CASE(wfi_linear_wf_method_nr)
        res="mo linear"
     CASE(wfi_linear_p_method_nr)
        res="P linear"
     CASE(wfi_linear_ps_method_nr)
        res="PS linear"
     CASE(wfi_ps_method_nr)
        res="PS Nth order"
     CASE(wfi_frozen_method_nr)
        res="frozen density approximation"
     CASE(wfi_aspc_nr)
       res = "ASPC"
     CASE default
        CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
             routineP,"Unknown interpolation method: "//&
             TRIM(ADJUSTL(cp_to_string(method_nr)))//&
             " in "//&
CPSourceFileRef,&
             error, failure)
     END SELECT
  END IF
END FUNCTION wfi_get_method_label

! *****************************************************************************
!> \brief calculates the new starting state for the scf for the next
!>      wf optimization
!> \param wf_history the previous history needed to extrapolate
!> \param qs_env the qs env with the latest result, and that will contain
!>        the new starting state
!> \param dt the time at which to extrapolate (wrt. to the last snapshot)
!> \param extrapolation_method_nr returns the extrapolation method used
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.2003 created [fawzi]
!>      11.2003 Joost VandeVondele : Implemented Nth order PS extrapolation
!> \author fawzi
! *****************************************************************************
  SUBROUTINE wfi_extrapolate(wf_history, qs_env, dt, extrapolation_method_nr, &
       orthogonal_wf, error)
    TYPE(qs_wf_history_type), POINTER        :: wf_history
    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(KIND=dp), INTENT(IN)                :: dt
    INTEGER, INTENT(OUT), OPTIONAL           :: extrapolation_method_nr
    LOGICAL, INTENT(OUT), OPTIONAL           :: orthogonal_wf
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'wfi_extrapolate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: actual_extrapolation_method_nr&
                                                , handle, i, ispin, k, n, &
                                                nmo, nvec, output_unit
    LOGICAL                                  :: failure, my_orthogonal_wf, &
                                                use_overlap
    REAL(KIND=dp)                            :: alpha, t0, t1, t2
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: ao_mo_fm_pools
    TYPE(cp_fm_struct_type), POINTER         :: matrix_struct, &
                                                matrix_struct_new
    TYPE(cp_fm_type), POINTER                :: csc, fm_tmp, mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_rho_type), POINTER               :: rho, rho_xc
    TYPE(qs_wf_snapshot_type), POINTER       :: t0_state, t1_state

    NULLIFY(mos, ao_mo_fm_pools, t0_state, t1_state, matrix_s, mo_coeff, &
         rho, rho_xc)
    failure=.FALSE.

    use_overlap = wf_history%store_overlap

    CALL timeset(routineN,handle)
    logger => cp_error_get_logger(error)
    output_unit = cp_print_key_unit_nr(logger,qs_env%input,"DFT%SCF%PRINT%PROGRAM_RUN_INFO",&
                  extension=".scfLog",error=error)

    CPPrecondition(ASSOCIATED(wf_history),cp_failure_level,routineP,error,failure)
    CPPrecondition(wf_history%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(qs_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL get_qs_env(qs_env, mos=mos,&
            rho=rho,rho_xc=rho_xc,matrix_s=matrix_s,&
            error=error)
       CALL mpools_get(qs_env%mpools, ao_mo_fm_pools=ao_mo_fm_pools,&
            error=error)

       ! chooses the method for this extrapolation
       IF (wf_history%snapshot_count<1) THEN
          actual_extrapolation_method_nr=wfi_use_guess_method_nr
       ELSE
          actual_extrapolation_method_nr=wf_history%interpolation_method_nr
       END IF

       SELECT CASE(actual_extrapolation_method_nr)
       CASE(wfi_linear_wf_method_nr)
           IF (wf_history%snapshot_count<2) THEN
             actual_extrapolation_method_nr=wfi_use_prev_wf_method_nr
          END IF
       CASE(wfi_linear_p_method_nr)
          IF (wf_history%snapshot_count<2) THEN
             actual_extrapolation_method_nr=wfi_use_prev_wf_method_nr
          END IF
       CASE(wfi_linear_ps_method_nr)
          IF (wf_history%snapshot_count<2) THEN
             actual_extrapolation_method_nr=wfi_use_prev_wf_method_nr
          END IF
       END SELECT

       IF (PRESENT(extrapolation_method_nr)) &
            extrapolation_method_nr=actual_extrapolation_method_nr
       my_orthogonal_wf=.FALSE.

       SELECT CASE (actual_extrapolation_method_nr)
       CASE(wfi_frozen_method_nr)
          t0_state => wfi_get_snapshot(wf_history, index=1, error=error)
          CPPrecondition(ASSOCIATED(t0_state%rho_frozen),cp_failure_level,routineP,error,failure)

          nvec = MIN(wf_history%memory_depth, wf_history%snapshot_count)
          CALL wfi_set_history_variables(qs_env=qs_env, nvec=nvec, error=error)

          DO ispin=1,SIZE(t0_state%rho_frozen%rho_ao)
             CALL cp_dbcsr_copy(rho%rho_ao(ispin)%matrix,&
                             t0_state%rho_frozen%rho_ao(ispin)%matrix,&
                             keep_sparsity=.TRUE.,&
                             error=error)
          END DO
          !FM updating rho_ao directly with t0_state%rho_ao would have the
          !FM wrong matrix structure
          CALL qs_rho_update_rho(rho, qs_env=qs_env, error=error)
          CALL qs_ks_did_change(qs_env%ks_env, rho_changed=.TRUE., error=error)

          my_orthogonal_wf=.FALSE.
       CASE(wfi_use_prev_rho_r_method_nr)
          t0_state => wfi_get_snapshot(wf_history, index=1, error=error)
          CPPrecondition(ASSOCIATED(t0_state%rho_ao),cp_failure_level,routineP,error,failure)
          nvec = MIN(wf_history%memory_depth, wf_history%snapshot_count)
          CALL wfi_set_history_variables(qs_env=qs_env, nvec=nvec, error=error)

          DO ispin=1,SIZE(t0_state%rho_ao)
             CALL cp_dbcsr_copy(rho%rho_ao(ispin)%matrix,&
                             t0_state%rho_ao(ispin)%matrix,&
                             keep_sparsity=.TRUE.,&
                             error=error)
          END DO
          rho%rho_g_valid=.TRUE.

          ! does nothing
       CASE(wfi_use_prev_wf_method_nr)
          my_orthogonal_wf=.TRUE.
          nvec = MIN(wf_history%memory_depth, wf_history%snapshot_count)
          CALL wfi_set_history_variables(qs_env=qs_env, nvec=nvec, error=error)

          DO ispin=1,SIZE(mos)
             CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff,&
                  nmo=nmo)
             CALL qs_env_reorthogonalize_vectors(qs_env,&
                  v_matrix=mo_coeff,&
                  n_col=nmo, error=error)
             CALL calculate_density_matrix(mo_set=mos(ispin)%mo_set,&
                  density_matrix=rho%rho_ao(ispin)%matrix,error=error)
          END DO
          CALL qs_rho_update_rho(rho, qs_env=qs_env, error=error)
          CALL qs_ks_did_change(qs_env%ks_env,rho_changed=.TRUE., error=error)

       CASE(wfi_use_prev_p_method_nr)
          t0_state => wfi_get_snapshot(wf_history, index=1, error=error)
          CPPrecondition(ASSOCIATED(t0_state%rho_ao),cp_failure_level,routineP,error,failure)
          nvec = MIN(wf_history%memory_depth, wf_history%snapshot_count)
          CALL wfi_set_history_variables(qs_env=qs_env, nvec=nvec, error=error)

          DO ispin=1,SIZE(t0_state%rho_ao)
             CALL cp_dbcsr_copy(rho%rho_ao(ispin)%matrix,&
                             t0_state%rho_ao(ispin)%matrix,&
                             keep_sparsity=.TRUE.,&
                             error=error)
          END DO
          !FM updating rho_ao directly with t0_state%rho_ao would have the
          !FM wrong matrix structure
          CALL qs_rho_update_rho(rho, qs_env=qs_env, error=error)
          CALL qs_ks_did_change(qs_env%ks_env, rho_changed=.TRUE., error=error)

       CASE(wfi_use_guess_method_nr)
          !FM more clean to do it here, but it
          !FM might need to read a file (restart) and thus globenv
          !FM I do not want globenv here, thus done by the caller
          !FM (btw. it also needs the eigensolver, and unless you relocate it
          !FM gives circular dependencies)
          nvec = MIN(wf_history%memory_depth, wf_history%snapshot_count)
          CALL wfi_set_history_variables(qs_env=qs_env, nvec=nvec, error=error)
       CASE(wfi_linear_wf_method_nr)
          t0_state => wfi_get_snapshot(wf_history, index=2, error=error)
          t1_state => wfi_get_snapshot(wf_history, index=1, error=error)
          CPPrecondition(ASSOCIATED(t0_state),cp_warning_level,routineP,error,failure)
          CPPrecondition(ASSOCIATED(t1_state),cp_warning_level,routineP,error,failure)
          IF (.not.failure) THEN
             CPPrecondition(ASSOCIATED(t0_state%wf),cp_warning_level,routineP,error,failure)
             CPPrecondition(ASSOCIATED(t1_state%wf),cp_warning_level,routineP,error,failure)
          END IF
          nvec = MIN(wf_history%memory_depth, wf_history%snapshot_count)
          CALL wfi_set_history_variables(qs_env=qs_env, nvec=nvec, error=error)

          IF (.NOT.failure) THEN
             my_orthogonal_wf=.TRUE.
             t0=0.0_dp
             t1=t1_state%dt
             t2=t1+dt
             DO ispin=1,SIZE(mos)
                CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff,&
                     nmo=nmo)
                CALL cp_fm_scale_and_add(alpha=0.0_dp,&
                     matrix_a=mo_coeff,&
                     matrix_b=t1_state%wf(ispin)%matrix,&
                     beta=(t2-t0)/(t1-t0),&
                     error=error) ! this copy should be unnecessary
                CALL cp_fm_scale_and_add(alpha=1.0_dp,&
                     matrix_a=mo_coeff,&
                     beta=(t1-t2)/(t1-t0), matrix_b=t0_state%wf(ispin)%matrix,&
                     error=error)
                CALL qs_env_reorthogonalize_vectors(qs_env,&
                     v_matrix=mo_coeff,&
                     n_col=nmo, error=error)
                CALL calculate_density_matrix(mo_set=mos(ispin)%mo_set,&
                     density_matrix=rho%rho_ao(ispin)%matrix,error=error)
             END DO
             CALL qs_rho_update_rho(rho, qs_env=qs_env,  error=error)

             CALL qs_ks_did_change(qs_env%ks_env,&
                  rho_changed=.TRUE., error=error)
          END IF
       CASE(wfi_linear_p_method_nr)
          t0_state => wfi_get_snapshot(wf_history, index=2, error=error)
          t1_state => wfi_get_snapshot(wf_history, index=1, error=error)
          CPPrecondition(ASSOCIATED(t0_state),cp_warning_level,routineP,error,failure)
          CPPrecondition(ASSOCIATED(t1_state),cp_warning_level,routineP,error,failure)
          IF (.not.failure) THEN
             CPPrecondition(ASSOCIATED(t0_state%rho_ao),cp_warning_level,routineP,error,failure)
             CPPrecondition(ASSOCIATED(t1_state%rho_ao),cp_warning_level,routineP,error,failure)
          END IF
          nvec = MIN(wf_history%memory_depth, wf_history%snapshot_count)
          CALL wfi_set_history_variables(qs_env=qs_env, nvec=nvec, error=error)

          IF (.not.failure) THEN
            t0=0.0_dp
            t1=t1_state%dt
            t2=t1+dt
            DO ispin=1,SIZE(rho%rho_ao)
               CALL cp_dbcsr_add(rho%rho_ao(ispin)%matrix,t1_state%rho_ao(ispin)%matrix,&
                    alpha_scalar=0.0_dp,beta_scalar=(t2-t0)/(t1-t0),error=error) ! this copy should be unnecessary
               CALL cp_dbcsr_add(rho%rho_ao(ispin)%matrix,t0_state%rho_ao(ispin)%matrix,&
                    alpha_scalar=1.0_dp,beta_scalar=(t1-t2)/(t1-t0),error=error)
            END DO
            CALL qs_rho_update_rho(rho, qs_env=qs_env, error=error)
            CALL qs_ks_did_change(qs_env%ks_env, rho_changed=.TRUE., error=error)

            ! wf not calculated, extract with PSC renormalized?
            ! use wf_linear?
          END IF
       CASE(wfi_linear_ps_method_nr)
          t0_state => wfi_get_snapshot(wf_history, index=2, error=error)
          t1_state => wfi_get_snapshot(wf_history, index=1, error=error)
          CPPrecondition(ASSOCIATED(t0_state),cp_warning_level,routineP,error,failure)
          CPPrecondition(ASSOCIATED(t1_state),cp_warning_level,routineP,error,failure)
          IF (.not.failure) THEN
             CPPrecondition(ASSOCIATED(t0_state%wf),cp_warning_level,routineP,error,failure)
             CPPrecondition(ASSOCIATED(t1_state%wf),cp_warning_level,routineP,error,failure)
             IF (wf_history%store_overlap) THEN
               CPPrecondition(ASSOCIATED(t0_state%overlap),cp_warning_level,routineP,error,failure)
               CPPrecondition(ASSOCIATED(t1_state%overlap),cp_warning_level,routineP,error,failure)
             END IF
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

          IF (.not.failure) THEN
             my_orthogonal_wf=.TRUE.
             ! use PS_2=2 PS_1-PS_0
             ! C_2 comes from using PS_2 as a projector acting on C_1
             DO ispin=1,SIZE(mos)
                NULLIFY(mo_coeff,matrix_struct,matrix_struct_new,csc)
                CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff)
                CALL cp_fm_get_info(mo_coeff,nrow_global=n,ncol_global=k, &
                     matrix_struct=matrix_struct,error=error)
                CALL cp_fm_struct_create(matrix_struct_new,template_fmstruct=matrix_struct, &
                                         nrow_global=k,ncol_global=k,error=error)
                CALL cp_fm_create(csc,matrix_struct_new,error=error)
                CALL cp_fm_struct_release(matrix_struct_new,error=error)

                IF ( use_overlap ) THEN
                  CALL cp_dbcsr_sm_fm_multiply(t0_state%overlap,t1_state%wf(ispin)%matrix,mo_coeff, k,error=error)
                  CALL cp_gemm('T','N',k,k,n,1.0_dp,t0_state%wf(ispin)%matrix,mo_coeff,0.0_dp,csc,error=error)
                ELSE
                  CALL cp_gemm('T','N',k,k,n,1.0_dp,t0_state%wf(ispin)%matrix,&
                                  t1_state%wf(ispin)%matrix,0.0_dp,csc,error=error)
                END IF
                CALL cp_gemm('N','N',n,k,k,1.0_dp,t0_state%wf(ispin)%matrix,csc,0.0_dp,mo_coeff,error=error)
                CALL cp_fm_release(csc,error=error)
                CALL cp_fm_scale_and_add(-1.0_dp,mo_coeff,2.0_dp,t1_state%wf(ispin)%matrix,error=error)
                CALL qs_env_reorthogonalize_vectors(qs_env,&
                     v_matrix=mo_coeff,&
                     n_col=k, error=error)
                CALL calculate_density_matrix(mo_set=mos(ispin)%mo_set,&
                     density_matrix=rho%rho_ao(ispin)%matrix,error=error)
             END DO
             CALL qs_rho_update_rho(rho, qs_env=qs_env, error=error)
             CALL qs_ks_did_change(qs_env%ks_env,rho_changed=.TRUE., error=error)

          END IF
      CASE(wfi_ps_method_nr)
          ! figure out the actual number of vectors to use in the extrapolation:
          nvec = MIN(wf_history%memory_depth, wf_history%snapshot_count)
          CPPrecondition(nvec .GT. 0,cp_failure_level,routineP,error,failure)
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

          my_orthogonal_wf=.TRUE.
          DO ispin=1,SIZE(mos)
             NULLIFY(mo_coeff,matrix_struct,matrix_struct_new,csc,fm_tmp)
             CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff)
             CALL cp_fm_get_info(mo_coeff,nrow_global=n,ncol_global=k, &
                  matrix_struct=matrix_struct,error=error)
             CALL cp_fm_create(fm_tmp,matrix_struct,error=error)
             CALL cp_fm_struct_create(matrix_struct_new,template_fmstruct=matrix_struct, &
                                      nrow_global=k,ncol_global=k,error=error)
             CALL cp_fm_create(csc,matrix_struct_new,error=error)
             CALL cp_fm_struct_release(matrix_struct_new,error=error)
             ! first the most recent
             t1_state => wfi_get_snapshot(wf_history, index=1, error=error)
             CALL cp_fm_to_fm(t1_state%wf(ispin)%matrix,mo_coeff,error=error)
             alpha = nvec
             CALL cp_fm_scale(alpha,mo_coeff,error=error)
             DO i=2,nvec
                t0_state => wfi_get_snapshot(wf_history, index=i, error=error)
                IF ( use_overlap ) THEN
                  CALL cp_dbcsr_sm_fm_multiply(t0_state%overlap,t1_state%wf(ispin)%matrix,fm_tmp, k,error=error)
                  CALL cp_gemm('T','N',k,k,n,1.0_dp,t0_state%wf(ispin)%matrix,fm_tmp,0.0_dp,csc,error=error)
                ELSE
                  CALL cp_gemm('T','N',k,k,n,1.0_dp,t0_state%wf(ispin)%matrix,&
                                  t1_state%wf(ispin)%matrix,0.0_dp,csc,error=error)
                END IF
                CALL cp_gemm('N','N',n,k,k,1.0_dp,t0_state%wf(ispin)%matrix,csc,0.0_dp,fm_tmp,error=error)
                alpha = -1.0_dp * alpha * REAL(nvec - i + 1 , dp ) / REAL(i , dp )
                CALL cp_fm_scale_and_add(1.0_dp,mo_coeff,alpha,fm_tmp,error=error)
             ENDDO

             CALL cp_fm_release(csc,error=error)
             CALL cp_fm_release(fm_tmp,error=error)
             CALL qs_env_reorthogonalize_vectors(qs_env,&
                  v_matrix=mo_coeff,&
                  n_col=k, error=error)
             CALL calculate_density_matrix(mo_set=mos(ispin)%mo_set,&
                  density_matrix=rho%rho_ao(ispin)%matrix,error=error)
          END DO
          CALL qs_rho_update_rho(rho, qs_env=qs_env, error=error)
          CALL qs_ks_did_change(qs_env%ks_env,rho_changed=.TRUE., error=error)

      CASE (wfi_aspc_nr)
        CALL cite_reference(Kolafa2004)
        ! figure out the actual number of vectors to use in the extrapolation:
        nvec = MIN(wf_history%memory_depth,wf_history%snapshot_count)
        CPPrecondition(nvec.GT.0,cp_failure_level,routineP,error,failure)
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
        DO ispin=1,SIZE(mos)
          NULLIFY (mo_coeff,matrix_struct,matrix_struct_new,csc,fm_tmp)
          CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff)
          CALL cp_fm_get_info(mo_coeff,&
                              nrow_global=n,&
                              ncol_global=k,&
                              matrix_struct=matrix_struct,error=error)
          CALL cp_fm_create(fm_tmp,matrix_struct,error=error)
          CALL cp_fm_struct_create(matrix_struct_new,&
                                   template_fmstruct=matrix_struct,&
                                   nrow_global=k,&
                                   ncol_global=k,error=error)
          CALL cp_fm_create(csc,matrix_struct_new,error=error)
          CALL cp_fm_struct_release(matrix_struct_new,error=error)
          ! first the most recent
          t1_state => wfi_get_snapshot(wf_history,&
                                       index=1,&
                                       error=error)
          CALL cp_fm_to_fm(t1_state%wf(ispin)%matrix,mo_coeff,error=error)
          alpha = REAL(4*nvec - 2,KIND=dp)/REAL(nvec + 1,KIND=dp)
          IF (output_unit>0) THEN
            WRITE (UNIT=output_unit,FMT="(/,T3,A,/,/,T3,A,I0,/,/,T3,A2,I0,A4,F10.6)")&
              "Parameters for the always stable predictor-corrector (ASPC) method:",&
              "ASPC order: ",MAX(nvec - 2,0),&
              "B(",1,") = ",alpha
          END IF
          CALL cp_fm_scale(alpha,mo_coeff,error=error)

          DO i=2,nvec
            t0_state => wfi_get_snapshot(wf_history,index=i,error=error)
            IF ( use_overlap ) THEN
              CALL cp_dbcsr_sm_fm_multiply(t0_state%overlap,t1_state%wf(ispin)%matrix,fm_tmp,k,error=error)
              CALL cp_gemm('T','N',k,k,n,1.0_dp,t0_state%wf(ispin)%matrix,fm_tmp,0.0_dp,csc,error=error)
            ELSE
              CALL cp_gemm('T','N',k,k,n,1.0_dp,t0_state%wf(ispin)%matrix,&
                              t1_state%wf(ispin)%matrix,0.0_dp,csc,error=error)
            END IF
            CALL cp_gemm('N','N',n,k,k,1.0_dp,t0_state%wf(ispin)%matrix,csc,0.0_dp,fm_tmp,error=error)
            alpha = (-1.0_dp)**(i + 1)*REAL(i,KIND=dp)*&
                    binomial(2*nvec,nvec - i)/binomial(2*nvec - 2,nvec -1)
            IF (output_unit>0) THEN
              WRITE (UNIT=output_unit,FMT="(T3,A2,I0,A4,F10.6)")&
                "B(",i,") = ",alpha
            END IF
            CALL cp_fm_scale_and_add(1.0_dp,mo_coeff,alpha,fm_tmp,error=error)
          END DO
          CALL cp_fm_release(csc,error=error)
          CALL cp_fm_release(fm_tmp,error=error)
          CALL qs_env_reorthogonalize_vectors(qs_env,&
                                              v_matrix=mo_coeff,&
                                              n_col=k,&
                                              error=error)
          CALL calculate_density_matrix(mo_set=mos(ispin)%mo_set,&
                                        density_matrix=rho%rho_ao(ispin)%matrix,&
                                        error=error)
        END DO
        CALL qs_rho_update_rho(rho, qs_env=qs_env, error=error)
        CALL qs_ks_did_change(qs_env%ks_env,rho_changed=.TRUE.,error=error)

      CASE default
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"Unknown interpolation method: "//&
               TRIM(ADJUSTL(cp_to_string(wf_history%interpolation_method_nr)))//&
               " in "//&
CPSourceFileRef,&
               error, failure)
       END SELECT
       IF (PRESENT(orthogonal_wf)) orthogonal_wf=my_orthogonal_wf
    END IF
    CALL cp_print_key_finished_output(output_unit,logger,qs_env%input,&
         "DFT%SCF%PRINT%PROGRAM_RUN_INFO",error=error)
    CALL timestop(handle)
  END SUBROUTINE wfi_extrapolate

! *****************************************************************************
!> \brief Decides if scf control variables has to changed due
!>      to using a WF extrapolation.
!> \param qs_env The QS environment
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      11.2006 created [TdK]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  SUBROUTINE wfi_set_history_variables(qs_env, nvec, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: nvec
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'wfi_set_history_variables', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)
    failure = .FALSE.

    CPPrecondition(ASSOCIATED(qs_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(qs_env%ref_count>0, cp_failure_level, routineP, error, failure)

    IF (.NOT. failure) THEN
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
    END IF

    CALL timestop(handle)

  END SUBROUTINE wfi_set_history_variables

! *****************************************************************************
!> \brief updates the snapshot buffer, taking a new snapshot
!> \param wf_history the history buffer to update
!> \param qs_env the qs_env we get the info from
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
SUBROUTINE wfi_update(wf_history, qs_env, dt, error)
    TYPE(qs_wf_history_type), POINTER        :: wf_history
    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(KIND=dp), INTENT(in)                :: dt
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'wfi_update', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(wf_history),cp_failure_level,routineP,error,failure)
  CPPrecondition(wf_history%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(qs_env%ref_count>0,cp_failure_level,routineP,error,failure)

  IF (.NOT. failure) THEN
     wf_history%snapshot_count=wf_history%snapshot_count+1
     IF (wf_history%memory_depth>0) THEN
        wf_history%last_state_index=MODULO(wf_history%snapshot_count,&
             wf_history%memory_depth)+1
        CALL wfs_update(snapshot=wf_history%past_states &
             (wf_history%last_state_index)%snapshot,wf_history=wf_history,&
             qs_env=qs_env,dt=dt,error=error)
     END IF
  END IF
END SUBROUTINE wfi_update

! *****************************************************************************
!> \brief chenge the memory depth stored in the history
!> \param wf_history the wf_history to modify
!> \param memory_depth the new memory_depth
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
SUBROUTINE wfi_change_memory_depth(wf_history,memory_depth,error)
    TYPE(qs_wf_history_type), POINTER        :: wf_history
    INTEGER, INTENT(in)                      :: memory_depth
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'wfi_change_memory_depth', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, stat
    LOGICAL                                  :: failure
    TYPE(qs_wf_snapshot_p_type), &
      DIMENSION(:), POINTER                  :: new_past_states

  failure=.FALSE.

  CALL timeset(routineN,handle)
  CPPrecondition(ASSOCIATED(wf_history),cp_failure_level,routineP,error,failure)
  CPPrecondition(wf_history%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     IF (wf_history%memory_depth/=memory_depth)THEN
        ALLOCATE(new_past_states(memory_depth),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        IF (.NOT. failure) THEN
           DO i=1,memory_depth
              NULLIFY(new_past_states(i)%snapshot)
              IF (i<wf_history%memory_depth) THEN
                 new_past_states(MODULO(wf_history%snapshot_count+1-i,&
                      wf_history%memory_depth)+1)%snapshot => &
                      wf_history%past_states(&
                      MODULO(wf_history%snapshot_count+1-i,&
                      wf_history%memory_depth)+1)%snapshot
              END IF
           END DO
           IF (ASSOCIATED(wf_history%past_states)) THEN
              DEALLOCATE(wf_history%past_states,stat=stat)
              CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
           END IF
           wf_history%past_states => new_past_states
           wf_history%memory_depth = memory_depth
        END IF
     END IF
  END IF
  CALL timestop(handle)
END SUBROUTINE wfi_change_memory_depth

END MODULE qs_wf_history_methods
