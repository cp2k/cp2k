!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief module that contains the definitions of the scf types
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
MODULE qs_scf_types
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_deallocate_matrix_set
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_p_type
  USE cp_fm_types,                     ONLY: cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE cp_fm_vect,                      ONLY: cp_fm_vect_dealloc
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE preconditioner_types,            ONLY: destroy_preconditioner,&
                                             preconditioner_p_type
  USE qs_block_davidson_types,         ONLY: block_davidson_release,&
                                             davidson_type
  USE qs_density_mixing_types,         ONLY: direct_mixing_nr,&
                                             mixing_storage_create,&
                                             mixing_storage_release,&
                                             mixing_storage_type
  USE qs_diis_types,                   ONLY: qs_diis_b_release,&
                                             qs_diis_buffer_type
  USE qs_ot_types,                     ONLY: qs_ot_destroy,&
                                             qs_ot_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf_types'

  INTEGER, PARAMETER, PUBLIC :: &
                                general_diag_method_nr=1,&
                                special_diag_method_nr=2,&
                                ot_diag_method_nr=3,&
                                block_krylov_diag_method_nr=4,&
                                block_davidson_diag_method_nr=5,&
                                ot_method_nr=10

  INTEGER, SAVE, PRIVATE :: last_scf_env_id=0

  PUBLIC :: qs_scf_env_type,qs_outer_scf_type
  PUBLIC :: scf_env_retain, scf_env_release, scf_env_create, scf_env_did_change
  PUBLIC :: krylov_space_type, krylov_space_create, subspace_env_type
  PUBLIC :: diag_subspace_env_create

! *****************************************************************************
!> \brief wrapper for temporary and cached objects used in the scf iteration
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
  TYPE qs_outer_scf_type
    INTEGER :: iter_count

    ! these are the variable of outer loop.
    ! right now, we assume that they can be easily written as
    ! small arrays, but we might want to go the cp_fm_types
    ! at a later stage
    ! also, we just store the full iteration history
    REAL(KIND=dp), DIMENSION(:), POINTER :: energy
    REAL(KIND=dp), DIMENSION(:,:), POINTER :: variables
    REAL(KIND=dp), DIMENSION(:,:), POINTER :: gradient
    INTEGER, DIMENSION(:), POINTER :: count

  END TYPE qs_outer_scf_type

  TYPE krylov_space_type
     INTEGER :: nkrylov, nblock, nmo_conv, nmo_nc, max_iter
     LOGICAL :: always_check_conv
     REAL(dp) :: eps_std_diag, eps_conv, eps_adapt, max_res_norm, min_res_norm
     REAL(dp), DIMENSION(:), POINTER :: c_eval, t_eval
     TYPE(cp_fm_p_type), DIMENSION(:), POINTER :: v_mat, mo_conv, mo_refine
     TYPE(cp_fm_type),  POINTER ::  tmp_mat
!NOTE: the following matrices are small and could be used as standard array rather than distributed fm
     TYPE(cp_fm_type),  POINTER :: block1_mat, block2_mat, block3_mat, block4_mat, block5_mat
     TYPE(cp_fm_p_type), DIMENSION(:), POINTER ::  c_vec, chc_mat
  END TYPE

  TYPE subspace_env_type
      INTEGER :: max_iter, mixing_method
      REAL(dp) :: eps_diag_sub, eps_ene, eps_adapt
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER :: p_matrix_store,p_matrix_mix
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER :: chc_mat, c_vec, c0
      TYPE(mixing_storage_type), POINTER :: mixing_store
  END TYPE



! *****************************************************************************
TYPE qs_scf_env_type
   TYPE(qs_outer_scf_type) :: outer_scf
   INTEGER :: ref_count, id_nr, print_count, iter_count
   INTEGER :: cholesky_method, nelectron, method, mixing_method, nskip_mixing
   REAL(KIND = dp) :: iter_param, iter_delta, p_mix_alpha
   CHARACTER (len=15) :: iter_method
   COMPLEX(KIND=dp), DIMENSION(:,:,:), POINTER :: cc_buffer
   LOGICAL :: do_diag_sub, print_iter_line, skip_mixing, skip_diis, needs_ortho
   TYPE(mixing_storage_type), POINTER :: mixing_store
   TYPE(cp_fm_p_type), DIMENSION(:), POINTER :: scf_work1
   TYPE(cp_fm_type), POINTER  :: scf_work2, ortho, ortho_m1, s_half, s_minus_one
   TYPE(krylov_space_type), POINTER :: krylov_space
   TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER :: p_delta, p_mix_new
   TYPE(preconditioner_p_type), DIMENSION(:), POINTER :: ot_preconditioner
   TYPE(qs_ot_type), POINTER, DIMENSION(:)  :: qs_ot_env
   TYPE(qs_diis_buffer_type), POINTER :: scf_diis_buffer
   TYPE(subspace_env_type), POINTER :: subspace_env
   TYPE(davidson_type), POINTER, DIMENSION(:) :: block_davidson_env
END TYPE qs_scf_env_type

CONTAINS

! *****************************************************************************
!> \brief allocates and initialize an scf_env
!> \param scf_env the scf env to initialize
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE scf_env_create(scf_env,  error)
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scf_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)

    failure=.FALSE.

    ALLOCATE(scf_env, stat=stat)
    CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    IF (.NOT.failure) THEN

       scf_env%ref_count=1
       scf_env%print_count=0
       last_scf_env_id=last_scf_env_id+1
       scf_env%id_nr=last_scf_env_id
       scf_env%print_count=0
       scf_env%iter_count=0
       scf_env%nelectron=0
       scf_env%iter_param=0.0_dp
       scf_env%iter_delta=0.0_dp
       scf_env%iter_method=""
       scf_env%print_iter_line=.TRUE.
       scf_env%skip_mixing=.FALSE.
       scf_env%skip_diis=.FALSE.
       scf_env%do_diag_sub=.FALSE.
       scf_env%needs_ortho=.FALSE.

       scf_env%method= -1
       scf_env%mixing_method = -1
       scf_env%p_mix_alpha=1.0_dp
       scf_env%cholesky_method = -1

       scf_env%outer_scf%iter_count=0
       NULLIFY(scf_env%outer_scf%variables, &
            scf_env%outer_scf%gradient,scf_env%outer_scf%energy,scf_env%outer_scf%count)

       NULLIFY(scf_env%scf_work1, scf_env%scf_work2, scf_env%ortho, scf_env%ortho_m1,&
            scf_env%p_mix_new, scf_env%ot_preconditioner, scf_env%qs_ot_env,&
            scf_env%scf_diis_buffer)
       NULLIFY (scf_env%s_half, scf_env%p_delta,  scf_env%s_minus_one)
       NULLIFY(scf_env%cc_buffer)
       NULLIFY(scf_env%mixing_store)
       NULLIFY(scf_env%krylov_space)
       NULLIFY(scf_env%subspace_env)
       NULLIFY(scf_env%block_davidson_env)

    END IF

    CALL timestop(handle)

  END SUBROUTINE scf_env_create

! *****************************************************************************
!> \brief retains an scf_env (see doc/ReferenceCounting.html)
!> \param scf_env the environment to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE scf_env_retain(scf_env,error)
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scf_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
     scf_env%ref_count=scf_env%ref_count+1
  END IF
  END SUBROUTINE scf_env_retain

! *****************************************************************************
!> \brief function to be called to inform the scf_env about changes
!> \param scf_env the scf env to inform
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      03.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE scf_env_did_change(scf_env,error)
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scf_env_did_change', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)

    failure=.FALSE.

  CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     IF (ASSOCIATED(scf_env%p_mix_new)) THEN
        CALL cp_dbcsr_deallocate_matrix_set(scf_env%p_mix_new,error=error)
     END IF
     IF (ASSOCIATED(scf_env%p_delta)) THEN
        CALL cp_dbcsr_deallocate_matrix_set(scf_env%p_delta,error=error)
     END IF
  END IF

    CALL timestop(handle)

  END SUBROUTINE scf_env_did_change

! *****************************************************************************
!> \brief releases an scf_env (see doc/ReferenceCounting.html)
!> \param scf_env the environment to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE scf_env_release(scf_env,error)
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scf_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(scf_env)) THEN
     CPPreconditionNoFail(scf_env%ref_count>0,cp_failure_level,routineP,error)
     scf_env%ref_count=scf_env%ref_count-1
     IF (scf_env%ref_count==0) THEN
        CALL cp_fm_vect_dealloc(scf_env%scf_work1,error=error)
        CALL cp_fm_release(scf_env%scf_work2,error=error)
        CALL cp_fm_release(scf_env%ortho,error=error)
        CALL cp_fm_release(scf_env%ortho_m1,error=error)
        CALL cp_fm_release(scf_env%s_half,error=error)
        CALL cp_fm_release(scf_env%s_minus_one,error=error)
        IF (ASSOCIATED(scf_env%p_mix_new)) THEN
           CPPreconditionNoFail(.TRUE.,cp_failure_level,routineP,error) ! we should not end up here, and give back using the pools
           CALL cp_dbcsr_deallocate_matrix_set(scf_env%p_mix_new,error=error)
        ENDIF
        IF (ASSOCIATED(scf_env%p_delta)) THEN
           CPPreconditionNoFail(.TRUE.,cp_failure_level,routineP,error) ! we should not end up here, and give back using the pools
           CALL cp_dbcsr_deallocate_matrix_set(scf_env%p_delta,error=error)
        ENDIF
        IF (ASSOCIATED(scf_env%ot_preconditioner)) THEN
           DO i=1,SIZE(scf_env%ot_preconditioner)
              CALL destroy_preconditioner(scf_env%ot_preconditioner(i)%preconditioner,error=error)
              DEALLOCATE(scf_env%ot_preconditioner(i)%preconditioner)
           ENDDO
           DEALLOCATE(scf_env%ot_preconditioner,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(scf_env%qs_ot_env)) THEN
           DO i=1,SIZE(scf_env%qs_ot_env)
              CALL qs_ot_destroy(scf_env%qs_ot_env(i),error=error)
           END DO
           DEALLOCATE(scf_env%qs_ot_env, stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        CALL qs_diis_b_release(scf_env%scf_diis_buffer,error=error)

        IF (ASSOCIATED(scf_env%outer_scf%variables)) THEN
           DEALLOCATE(scf_env%outer_scf%variables,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        ENDIF
        IF (ASSOCIATED(scf_env%outer_scf%count)) THEN
           DEALLOCATE(scf_env%outer_scf%count,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        ENDIF
        IF (ASSOCIATED(scf_env%outer_scf%gradient)) THEN
           DEALLOCATE(scf_env%outer_scf%gradient,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        ENDIF
        IF (ASSOCIATED(scf_env%outer_scf%energy)) THEN
           DEALLOCATE(scf_env%outer_scf%energy,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        ENDIF

        IF (ASSOCIATED(scf_env%cc_buffer)) THEN
           DEALLOCATE(scf_env%cc_buffer,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        ENDIF

        IF (ASSOCIATED(scf_env%mixing_store)) THEN
           CALL mixing_storage_release(scf_env%mixing_store,error=error)
        END IF

        IF (ASSOCIATED(scf_env%krylov_space)) THEN
           CALL krylov_space_release(scf_env%krylov_space,error=error)
        END IF

        IF (ASSOCIATED(scf_env%subspace_env)) THEN
           CALL diag_subspace_env_release(scf_env%subspace_env,error=error)
        END IF

        IF (ASSOCIATED(scf_env%block_davidson_env)) THEN
           CALL block_davidson_release(scf_env%block_davidson_env,error=error)
        END IF

        DEALLOCATE(scf_env, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
    END IF
    NULLIFY(scf_env)

  END SUBROUTINE scf_env_release

! *****************************************************************************
!> \brief  creates krylov space
!> \par History
!>      05.2009 created [MI]
!> \author [MI]
! *****************************************************************************
 SUBROUTINE krylov_space_create(krylov_space,scf_section, error)

    TYPE(krylov_space_type), POINTER         :: krylov_space
    TYPE(section_vals_type), POINTER         :: scf_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'krylov_space_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(krylov_space),cp_failure_level,routineP,error,failure)
    IF(.NOT. failure) THEN
      ALLOCATE(krylov_space, stat=stat)
      CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      IF (.NOT.failure) THEN


        NULLIFY(krylov_space%c_eval,krylov_space%t_eval)
        NULLIFY(krylov_space%v_mat)
        NULLIFY(krylov_space%mo_conv, krylov_space%mo_refine)
        NULLIFY(krylov_space%chc_mat, krylov_space%c_vec)
        NULLIFY(krylov_space%tmp_mat)
        NULLIFY(krylov_space%block1_mat, krylov_space%block2_mat)
        NULLIFY(krylov_space%block3_mat, krylov_space%block4_mat, krylov_space%block5_mat)

        CALL section_vals_val_get(scf_section,"DIAGONALIZATION%MAX_ITER",&
             i_val=krylov_space%max_iter,error=error)
        CALL section_vals_val_get(scf_section,"DIAGONALIZATION%KRYLOV%NKRYLOV",&
             i_val=krylov_space%nkrylov,error=error)
        CALL section_vals_val_get(scf_section,"DIAGONALIZATION%KRYLOV%NBLOCK",&
             i_val=krylov_space%nblock,error=error)
        CALL section_vals_val_get(scf_section,"DIAGONALIZATION%KRYLOV%EPS_KRYLOV",&
             r_val=krylov_space%eps_conv,error=error)
        CALL section_vals_val_get(scf_section,"DIAGONALIZATION%KRYLOV%EPS_STD_DIAG",&
             r_val=krylov_space%eps_std_diag,error=error)
        CALL section_vals_val_get(scf_section,"DIAGONALIZATION%EPS_ADAPT",&
             r_val=krylov_space%eps_adapt,error=error)
        CALL section_vals_val_get(scf_section,"DIAGONALIZATION%KRYLOV%CHECK_MOS_CONV",&
             l_val=krylov_space%always_check_conv,error=error)
      END IF
    END IF

 END SUBROUTINE krylov_space_create

! *****************************************************************************
!> \brief releases krylov space
!> \par History
!>      05.2009 created [MI]
!> \author [MI]
! *****************************************************************************
 SUBROUTINE krylov_space_release(krylov_space,error)
    TYPE(krylov_space_type), POINTER         :: krylov_space
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'krylov_space_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure = .FALSE.

    IF(ASSOCIATED(krylov_space)) THEN

        DEALLOCATE(krylov_space%c_eval,STAT=stat)
      CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DEALLOCATE(krylov_space%t_eval,STAT=stat)
      CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

        DO i = 1,SIZE(krylov_space%v_mat)
          CALL cp_fm_release(krylov_space%v_mat(i)%matrix,error=error)
        END DO
        DEALLOCATE(krylov_space%v_mat,STAT=stat)
      CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

        DO i = 1,SIZE(krylov_space%mo_conv)
          CALL cp_fm_release(krylov_space%mo_conv(i)%matrix,error=error)
        END DO
        DEALLOCATE(krylov_space%mo_conv,STAT=stat)
      CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

        DO i = 1,SIZE(krylov_space%mo_refine)
          CALL cp_fm_release(krylov_space%mo_refine(i)%matrix,error=error)
        END DO
        DEALLOCATE(krylov_space%mo_refine,STAT=stat)
      CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

        DO i = 1,SIZE(krylov_space%chc_mat)
          CALL cp_fm_release(krylov_space%chc_mat(i)%matrix,error=error)
        END DO
        DEALLOCATE(krylov_space%chc_mat,STAT=stat)
        DO i = 1,SIZE(krylov_space%c_vec)
          CALL cp_fm_release(krylov_space%c_vec(i)%matrix,error=error)
        END DO
        DEALLOCATE(krylov_space%c_vec,STAT=stat)
        CALL cp_fm_release(krylov_space%tmp_mat,error=error)
        CALL cp_fm_release(krylov_space%block1_mat,error=error)
        CALL cp_fm_release(krylov_space%block2_mat,error=error)
        CALL cp_fm_release(krylov_space%block3_mat,error=error)
        CALL cp_fm_release(krylov_space%block4_mat,error=error)
        CALL cp_fm_release(krylov_space%block5_mat,error=error)

        DEALLOCATE(krylov_space,STAT=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

      NULLIFY(krylov_space)
    END IF

 END SUBROUTINE krylov_space_release

! *****************************************************************************
!> \brief creates subspace-rotation environment
!> \par History
!>      09.2009 created [MI]
!> \author [MI]
! *****************************************************************************
 SUBROUTINE diag_subspace_env_create(subspace_env,scf_section,ecut,error)

    TYPE(subspace_env_type), POINTER         :: subspace_env
    TYPE(section_vals_type), POINTER         :: scf_section
    REAL(dp), INTENT(IN)                     :: ecut
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'diag_subspace_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: do_mixing, failure
    TYPE(section_vals_type), POINTER         :: mixing_section

    failure = .FALSE.

    CPPrecondition(.NOT.ASSOCIATED(subspace_env),cp_failure_level,routineP,error,failure)
    IF(.NOT. failure) THEN
      ALLOCATE(subspace_env, stat=stat)
      CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      IF (.NOT.failure) THEN


        NULLIFY(subspace_env%p_matrix_store)
        NULLIFY(subspace_env%p_matrix_mix)
        NULLIFY(subspace_env%chc_mat)
        NULLIFY(subspace_env%c_vec)
        NULLIFY(subspace_env%c0)
        NULLIFY(subspace_env%mixing_store)
        NULLIFY(mixing_section)

        CALL section_vals_val_get(scf_section,"DIAGONALIZATION%DIAG_SUB_SCF%MAX_ITER",&
             i_val=subspace_env%max_iter,error=error)
        CALL section_vals_val_get(scf_section,"DIAGONALIZATION%DIAG_SUB_SCF%EPS_ENE",&
             r_val=subspace_env%eps_ene,error=error)
        CALL section_vals_val_get(scf_section,"DIAGONALIZATION%DIAG_SUB_SCF%EPS_SKIP_SUB_DIAG",&
             r_val=subspace_env%eps_diag_sub,error=error)
        CALL section_vals_val_get(scf_section,"DIAGONALIZATION%DIAG_SUB_SCF%EPS_ADAPT_SCF",&
             r_val=subspace_env%eps_adapt,error=error)
        subspace_env%mixing_method = 0
        do_mixing = .FALSE.
        mixing_section => section_vals_get_subs_vals(scf_section,"DIAGONALIZATION%DIAG_SUB_SCF%MIXING",error=error)
        CALL section_vals_val_get(mixing_section,"_SECTION_PARAMETERS_",&
                                 l_val=do_mixing,&
                                 error=error)
        IF (do_mixing) THEN
          CALL section_vals_val_get(mixing_section,"METHOD",&
                                    i_val=subspace_env%mixing_method,&
                                    error=error)
          IF(subspace_env%mixing_method>=direct_mixing_nr)&
             CALL mixing_storage_create(subspace_env%mixing_store, mixing_section, &
                  subspace_env%mixing_method, ecut=ecut, error=error)
        END IF

      END IF
    END IF

 END SUBROUTINE diag_subspace_env_create

! *****************************************************************************
!> \brief releases subspace-rotation environment
!> \par History
!>      09.2009 created [MI]
!> \author [MI]
! *****************************************************************************
 SUBROUTINE diag_subspace_env_release(subspace_env,error)
    TYPE(subspace_env_type), POINTER         :: subspace_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'diag_subspace_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure = .FALSE.

    IF(ASSOCIATED(subspace_env)) THEN

       IF (ASSOCIATED(subspace_env%p_matrix_store)) THEN

           CPPreconditionNoFail(.TRUE.,cp_failure_level,routineP,error) ! should be given back using the pools
           CALL cp_dbcsr_deallocate_matrix_set(subspace_env%p_matrix_store,error=error)
       ENDIF
       DO i = 1,SIZE(subspace_env%chc_mat)
         CALL cp_fm_release(subspace_env%chc_mat(i)%matrix,error=error)
       END DO
       DEALLOCATE(subspace_env%chc_mat,STAT=stat)
       DO i = 1,SIZE(subspace_env%c_vec)
         CALL cp_fm_release(subspace_env%c_vec(i)%matrix,error=error)
       END DO
       DEALLOCATE(subspace_env%c_vec,STAT=stat)
       DO i = 1,SIZE(subspace_env%c0)
         CALL cp_fm_release(subspace_env%c0(i)%matrix,error=error)
       END DO
       DEALLOCATE(subspace_env%c0,STAT=stat)

        IF (ASSOCIATED(subspace_env%mixing_store)) THEN
           CALL mixing_storage_release(subspace_env%mixing_store,error=error)
        END IF

       DEALLOCATE(subspace_env, stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF

 END SUBROUTINE diag_subspace_env_release

END MODULE qs_scf_types
