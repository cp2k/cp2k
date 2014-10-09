!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      11.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE qs_ks_types
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_deallocate_matrix_set,&
                                             cp_dbcsr_p_type
  USE f77_blas
  USE pw_types,                        ONLY: pw_release,&
                                             pw_type
  USE qs_rho_types,                    ONLY: qs_rho_create,&
                                             qs_rho_release,&
                                             qs_rho_retain,&
                                             qs_rho_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ks_types'

  PUBLIC :: qs_ks_env_type, qs_ks_env_create
  PUBLIC :: qs_ks_release, qs_ks_retain, get_ks_env, set_ks_env

! *****************************************************************************
!> \brief calculation environement to calculate the ks matrix,
!>      holds all the needed vars.
!>      assumes that the core hamiltonian and energy are up to date.
!> \param s_mstruct_changed if the sparsity structure of s in qs_env where this
!>        ks_env lives has changed
!> \param rho_changed if the rho stored in the qs_env has changed
!> \param forces_up_to_date if the forces are up to date
!> \param potential_changed if some external potentian changed (such as due to the constraint)
!> \par History
!>      04.2002 created [fawzi]
!>      10.2002 introduced pools, uses updated rho as input,
!>              removed most temporary variables, renamed may vars,
!>              began conversion to LSD [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE qs_ks_env_type
     PRIVATE
     LOGICAL                                               :: s_mstruct_changed = .TRUE.,&
                                                              rho_changed = .TRUE.,&
                                                              potential_changed = .TRUE.,&
                                                              forces_up_to_date = .FALSE.
     INTEGER                                               :: n_evals = 0,&
                                                              id_nr = -1,&
                                                              ref_count = 1

     ! hartree is supposed to contain the hartree potential (for cube ouput)
     ! ugly to keep it always around only for a cube output...
     TYPE(pw_type), POINTER                                :: v_hartree_rspace => Null()
     TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_h => Null(),&
                                                              matrix_ks => Null(),&
                                                              matrix_vxc => Null(),&
                                                              matrix_s => Null(),&
                                                              kinetic => Null(),&
                                                              matrix_w => Null(),&
                                                              matrix_w_mp2 => Null(),&
                                                              matrix_p_mp2 => Null(),&
                                                              gamma_matrix => Null(),&
                                                              matrix_s_aux_fit => Null(),&
                                                              matrix_s_aux_fit_vs_orb => Null(),&
                                                              matrix_ks_aux_fit => Null()


    TYPE(qs_rho_type), POINTER                            :: rho => Null(),&
                                                             rho_buffer => Null(),&
                                                             rho_xc => Null(),&
                                                             rho_aux_fit => Null()

  END TYPE qs_ks_env_type

! *****************************************************************************
!> \brief type to build arrays of pointers
!> \param ks_env the ks_env pointer
!> \par History
!>      11.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE qs_ks_env_p_type
     TYPE(qs_ks_env_type), POINTER :: ks_env
  END TYPE qs_ks_env_p_type

CONTAINS

! *****************************************************************************
!> \brief Allocates a new instance of ks_env.
!> \param ks_env ...
!> \param error ...
!> \author Ole Schuett
! *****************************************************************************
SUBROUTINE qs_ks_env_create(ks_env, error)
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_ks_env_create', &
      routineP = moduleN//':'//routineN

    IF(ASSOCIATED(ks_env)) STOP routineP//"ks_env already associated"
    ALLOCATE(ks_env)
    CALL qs_rho_create(ks_env%rho, error)
    CALL qs_rho_create(ks_env%rho_aux_fit, error)
    CALL qs_rho_create(ks_env%rho_xc, error)
    CALL qs_rho_create(ks_env%rho_buffer, error)
END SUBROUTINE qs_ks_env_create


! *****************************************************************************
!> \brief ...
!> \param ks_env ...
!> \param v_hartree_rspace ...
!> \param s_mstruct_changed ...
!> \param rho_changed ...
!> \param potential_changed ...
!> \param forces_up_to_date ...
!> \param matrix_h ...
!> \param matrix_ks ...
!> \param matrix_vxc ...
!> \param matrix_ks_aux_fit ...
!> \param kinetic ...
!> \param matrix_s ...
!> \param matrix_s_aux_fit ...
!> \param matrix_s_aux_fit_vs_orb ...
!> \param matrix_w ...
!> \param matrix_w_mp2 ...
!> \param matrix_p_mp2 ...
!> \param gamma_matrix ...
!> \param rho ...
!> \param rho_aux_fit ...
!> \param rho_buffer ...
!> \param rho_xc ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE get_ks_env(ks_env,v_hartree_rspace,&
                        s_mstruct_changed,rho_changed,&
                        potential_changed,forces_up_to_date,&
                        matrix_h,matrix_ks,matrix_vxc,matrix_ks_aux_fit,&
                        kinetic,matrix_s,matrix_s_aux_fit,matrix_s_aux_fit_vs_orb,&
                        matrix_w,matrix_w_mp2,matrix_p_mp2,gamma_matrix,&
                        rho, rho_aux_fit, rho_buffer, rho_xc,&
                        error)

    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(pw_type), OPTIONAL, POINTER         :: v_hartree_rspace
    LOGICAL, OPTIONAL                        :: s_mstruct_changed, &
                                                rho_changed, &
                                                potential_changed, &
                                                forces_up_to_date
    TYPE(cp_dbcsr_p_type), DIMENSION(:), OPTIONAL, POINTER :: matrix_h, &
      matrix_ks, matrix_vxc, matrix_ks_aux_fit, kinetic, matrix_s, &
      matrix_s_aux_fit, matrix_s_aux_fit_vs_orb, matrix_w, matrix_w_mp2, &
      matrix_p_mp2, gamma_matrix
    TYPE(qs_rho_type), OPTIONAL, POINTER     :: rho, rho_aux_fit, rho_buffer, &
                                                rho_xc
    TYPE(cp_error_type), INTENT(INOUT)       :: error

     IF(.NOT. ASSOCIATED(ks_env)) STOP "get_ks_env: not associated"
     IF(ks_env%ref_count<1) STOP "get_ks_env: ks_env%ref_count<1"

     IF (PRESENT(s_mstruct_changed)) s_mstruct_changed = ks_env%s_mstruct_changed
     IF (PRESENT(rho_changed)) rho_changed = ks_env%rho_changed
     IF (PRESENT(potential_changed)) potential_changed = ks_env%potential_changed
     IF (PRESENT(forces_up_to_date)) forces_up_to_date = ks_env%forces_up_to_date
     IF (PRESENT(v_hartree_rspace)) v_hartree_rspace => ks_env%v_hartree_rspace
     IF (PRESENT(matrix_h)) matrix_h => ks_env%matrix_h
     IF (PRESENT(matrix_ks)) matrix_ks => ks_env%matrix_ks
     IF (PRESENT(matrix_vxc)) matrix_vxc => ks_env%matrix_vxc
     IF (PRESENT(matrix_ks_aux_fit)) matrix_ks_aux_fit => ks_env%matrix_ks_aux_fit
     IF (PRESENT(kinetic)) kinetic => ks_env%kinetic
     IF (PRESENT(matrix_s)) matrix_s => ks_env%matrix_s
     IF (PRESENT(matrix_s_aux_fit)) matrix_s_aux_fit => ks_env%matrix_s_aux_fit
     IF (PRESENT(matrix_s_aux_fit_vs_orb)) matrix_s_aux_fit_vs_orb => ks_env%matrix_s_aux_fit_vs_orb
     IF (PRESENT(matrix_w)) matrix_w => ks_env%matrix_w
     IF (PRESENT(matrix_w_mp2)) matrix_w_mp2 => ks_env%matrix_w_mp2
     IF (PRESENT(matrix_p_mp2)) matrix_p_mp2 => ks_env%matrix_p_mp2
     IF (PRESENT(gamma_matrix)) gamma_matrix => ks_env%gamma_matrix
     IF (PRESENT(rho)) rho => ks_env%rho
     IF (PRESENT(rho_aux_fit)) rho_aux_fit => ks_env%rho_aux_fit
     IF (PRESENT(rho_buffer)) rho_buffer => ks_env%rho_buffer
     IF (PRESENT(rho_xc)) rho_xc => ks_env%rho_xc

  END SUBROUTINE get_ks_env


! *****************************************************************************
!> \brief ...
!> \param ks_env ...
!> \param v_hartree_rspace ...
!> \param s_mstruct_changed ...
!> \param rho_changed ...
!> \param potential_changed ...
!> \param forces_up_to_date ...
!> \param matrix_h ...
!> \param matrix_ks ...
!> \param matrix_vxc ...
!> \param matrix_ks_aux_fit ...
!> \param kinetic ...
!> \param matrix_s ...
!> \param matrix_s_aux_fit ...
!> \param matrix_s_aux_fit_vs_orb ...
!> \param matrix_w ...
!> \param matrix_w_mp2 ...
!> \param matrix_p_mp2 ...
!> \param gamma_matrix ...
!> \param rho ...
!> \param rho_aux_fit ...
!> \param rho_buffer ...
!> \param rho_xc ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE set_ks_env(ks_env,v_hartree_rspace,&
                        s_mstruct_changed,rho_changed,&
                        potential_changed,forces_up_to_date,&
                        matrix_h,matrix_ks,matrix_vxc,matrix_ks_aux_fit,&
                        kinetic,matrix_s,matrix_s_aux_fit,matrix_s_aux_fit_vs_orb,&
                        matrix_w,matrix_w_mp2,matrix_p_mp2,gamma_matrix,&
                        rho, rho_aux_fit, rho_buffer, rho_xc,&
                        error)

    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(pw_type), OPTIONAL, POINTER         :: v_hartree_rspace
    LOGICAL, OPTIONAL                        :: s_mstruct_changed, &
                                                rho_changed, &
                                                potential_changed, &
                                                forces_up_to_date
    TYPE(cp_dbcsr_p_type), DIMENSION(:), OPTIONAL, POINTER :: matrix_h, &
      matrix_ks, matrix_vxc, matrix_ks_aux_fit, kinetic, matrix_s, &
      matrix_s_aux_fit, matrix_s_aux_fit_vs_orb, matrix_w, matrix_w_mp2, &
      matrix_p_mp2, gamma_matrix
    TYPE(qs_rho_type), OPTIONAL, POINTER     :: rho, rho_aux_fit, rho_buffer, &
                                                rho_xc
    TYPE(cp_error_type), INTENT(INOUT)       :: error

     IF(.NOT. ASSOCIATED(ks_env)) STOP "set_ks_env: not associated"
     IF(ks_env%ref_count<1) STOP "set_ks_env: ks_env%ref_count<1"

     IF (PRESENT(s_mstruct_changed)) ks_env%s_mstruct_changed=s_mstruct_changed
     IF (PRESENT(rho_changed)) ks_env%rho_changed=rho_changed
     IF (PRESENT(potential_changed)) ks_env%potential_changed=potential_changed
     IF (PRESENT(forces_up_to_date)) ks_env%forces_up_to_date=forces_up_to_date
     IF (ks_env%s_mstruct_changed .OR. ks_env%potential_changed .OR. ks_env%rho_changed)&
          ks_env%forces_up_to_date=.FALSE.
     IF (PRESENT(v_hartree_rspace)) ks_env%v_hartree_rspace => v_hartree_rspace
     IF (PRESENT(matrix_h)) ks_env%matrix_h => matrix_h
     IF (PRESENT(matrix_ks)) ks_env%matrix_ks => matrix_ks
     IF (PRESENT(matrix_vxc)) ks_env%matrix_vxc => matrix_vxc
     IF (PRESENT(matrix_ks_aux_fit)) ks_env%matrix_ks_aux_fit => matrix_ks_aux_fit
     IF (PRESENT(matrix_s)) ks_env%matrix_s => matrix_s
     IF (PRESENT(matrix_s_aux_fit)) ks_env%matrix_s_aux_fit => matrix_s_aux_fit
     IF (PRESENT(matrix_s_aux_fit_vs_orb)) ks_env%matrix_s_aux_fit_vs_orb => matrix_s_aux_fit_vs_orb
     IF (PRESENT(matrix_w)) ks_env%matrix_w => matrix_w
     IF (PRESENT(matrix_w_mp2)) ks_env%matrix_w_mp2 => matrix_w_mp2
     IF (PRESENT(matrix_p_mp2)) ks_env%matrix_p_mp2 => matrix_p_mp2
     IF (PRESENT(gamma_matrix)) ks_env%gamma_matrix => gamma_matrix
     IF (PRESENT(kinetic)) ks_env%kinetic => kinetic

     IF (PRESENT(rho)) THEN ! accepts also null pointers !
      IF (ASSOCIATED(rho)) CALL qs_rho_retain(rho,error=error)
      CALL qs_rho_release(ks_env%rho,error=error)
      ks_env%rho => rho
     END IF
     IF (PRESENT(rho_aux_fit)) THEN ! accepts also null pointers !
       IF (ASSOCIATED(rho_aux_fit)) CALL qs_rho_retain(rho_aux_fit,error=error)
       CALL qs_rho_release(ks_env%rho_aux_fit,error=error)
       ks_env%rho_aux_fit => rho_aux_fit
     END IF
     IF (PRESENT(rho_buffer)) THEN ! accepts also null pointers !
       IF (ASSOCIATED(rho_buffer)) CALL qs_rho_retain(rho_buffer,error=error)
       CALL qs_rho_release(ks_env%rho_buffer,error=error)
       ks_env%rho_buffer => rho_buffer
     END IF
     IF (PRESENT(rho_xc)) THEN ! accepts also null pointers !
       IF (ASSOCIATED(rho_xc)) CALL qs_rho_retain(rho_xc,error=error)
       CALL qs_rho_release(ks_env%rho_xc,error=error)
       ks_env%rho_xc => rho_xc
     END IF

  END SUBROUTINE set_ks_env

! *****************************************************************************
!> \brief releases the ks_env (see doc/ReferenceCounting.html)
!> \param ks_env the ks_env to be released
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      05.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE qs_ks_release(ks_env,error)
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ks_release', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(ks_env)) THEN
       CPPrecondition(ks_env%ref_count>0,cp_failure_level,routineP,error,failure)
       ks_env%ref_count=ks_env%ref_count-1

       IF (ks_env%ref_count<1) THEN
          IF (ASSOCIATED(ks_env%v_hartree_rspace)) &
            CALL pw_release(ks_env%v_hartree_rspace,error=error)
          IF (ASSOCIATED(ks_env%matrix_h)) &
            CALL cp_dbcsr_deallocate_matrix_set(ks_env%matrix_h,error)
          IF (ASSOCIATED(ks_env%matrix_ks)) &
            CALL cp_dbcsr_deallocate_matrix_set(ks_env%matrix_ks,error)
          IF (ASSOCIATED(ks_env%matrix_vxc)) &
            CALL cp_dbcsr_deallocate_matrix_set(ks_env%matrix_vxc,error)
          IF (ASSOCIATED(ks_env%matrix_ks_aux_fit)) &
            CALL cp_dbcsr_deallocate_matrix_set(ks_env%matrix_ks_aux_fit,error)
          IF (ASSOCIATED(ks_env%matrix_s)) &
            CALL cp_dbcsr_deallocate_matrix_set(ks_env%matrix_s,error)
          IF (ASSOCIATED(ks_env%matrix_s_aux_fit)) &
            CALL cp_dbcsr_deallocate_matrix_set(ks_env%matrix_s_aux_fit,error)
          IF (ASSOCIATED(ks_env%matrix_s_aux_fit_vs_orb)) &
            CALL cp_dbcsr_deallocate_matrix_set(ks_env%matrix_s_aux_fit_vs_orb,error)
          IF (ASSOCIATED(ks_env%matrix_w)) &
            CALL cp_dbcsr_deallocate_matrix_set(ks_env%matrix_w,error)
          IF (ASSOCIATED(ks_env%matrix_w_mp2)) &
            CALL cp_dbcsr_deallocate_matrix_set(ks_env%matrix_w_mp2,error)
          IF (ASSOCIATED(ks_env%matrix_p_mp2)) &
            CALL cp_dbcsr_deallocate_matrix_set(ks_env%matrix_p_mp2,error)
          IF (ASSOCIATED(ks_env%gamma_matrix)) &
            CALL cp_dbcsr_deallocate_matrix_set(ks_env%gamma_matrix,error)
          IF (ASSOCIATED(ks_env%kinetic))&
            CALL cp_dbcsr_deallocate_matrix_set(ks_env%kinetic,error)
          IF(ASSOCIATED(ks_env%rho))&
            CALL qs_rho_release(ks_env%rho, error=error)
          IF(ASSOCIATED(ks_env%rho_buffer))&
            CALL qs_rho_release(ks_env%rho_buffer,error=error)
          IF(ASSOCIATED(ks_env%rho_xc))&
            CALL qs_rho_release(ks_env%rho_xc, error=error)
          IF(ASSOCIATED(ks_env%rho_aux_fit))&
            CALL qs_rho_release(ks_env%rho_aux_fit, error=error)

          DEALLOCATE(ks_env)
       END IF
    END IF
    NULLIFY(ks_env)
  END SUBROUTINE qs_ks_release

! *****************************************************************************
!> \brief retains the given ks_environment
!> \param ks_env the KohnSham environment to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      11.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE qs_ks_retain(ks_env, error)
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ks_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(ks_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(ks_env%ref_count>0,cp_failure_level,routineP,error,failure)
     ks_env%ref_count=ks_env%ref_count+1
  END IF
END SUBROUTINE qs_ks_retain

END MODULE qs_ks_types
