!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Contains types used for a Dimer Method calculations
!> \par History
!>      none
!> \author Luca Bellucci and Teodoro Laino - created [tlaino] - 01.2008
! *****************************************************************************
MODULE dimer_types
  
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: do_first_rotation_step
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE parallel_rng_types,              ONLY: random_numbers
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dimer_types'
  INTEGER, PRIVATE, SAVE :: last_dimer_id=0

  PUBLIC :: dimer_env_type,&
            dimer_env_create,&
            dimer_env_retain,&
            dimer_env_release

! *****************************************************************************
!> \brief Type containing all informations abour the rotation of the Dimer
!> \par History
!>      none
!> \author Luca Bellucci and Teodoro Laino - created [tlaino] - 01.2008
! *****************************************************************************
  TYPE dimer_rotational_type
     ! Rotational parameters
     INTEGER                                    :: rotation_step
     LOGICAL                                    :: interpolate_gradient
     REAL(KIND=dp)                              :: angle_tol, angle1, angle2, dCdp, curvature
     REAL(KIND=dp), POINTER, DIMENSION(:)       :: g0, g1, g1p
  END TYPE dimer_rotational_type

! *****************************************************************************
!> \brief Type containing all informations abour the translation of the Dimer
!> \par History
!>      none
!> \author Luca Bellucci and Teodoro Laino - created [tlaino] - 01.2008
! *****************************************************************************
  TYPE dimer_translational_type
     ! Translational parameters
     REAL(KIND=dp), POINTER, DIMENSION(:)       :: tls_vec
  END TYPE dimer_translational_type

! *****************************************************************************
!> \brief Conjugate Directions type
!> \par History
!>      none
!> \author Luca Bellucci and Teodoro Laino - created [tlaino] - 01.2008
! *****************************************************************************
  TYPE dimer_cg_rot_type
     REAL(KIND=dp)                              :: norm_theta, norm_theta_old, norm_h
     REAL(KIND=dp), POINTER, DIMENSION(:)       :: nvec_old
  END TYPE dimer_cg_rot_type

! *****************************************************************************
!> \brief Defines the environment for a Dimer Method calculation
!> \par History
!>      none
!> \author Luca Bellucci and Teodoro Laino - created [tlaino] - 01.2008
! *****************************************************************************
  TYPE dimer_env_type
     INTEGER                                    :: ref_count, id_nr
     REAL(KIND=dp)                              :: dr
     REAL(KIND=dp), POINTER, DIMENSION(:)       :: nvec
     TYPE(dimer_rotational_type)                :: rot
     TYPE(dimer_translational_type)             :: tsl
     TYPE(dimer_cg_rot_type)                    :: cg_rot
  END TYPE dimer_env_type

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param dimer_env ...
!> \param natom ...
!> \param globenv ...
!> \param dimer_section ...
!> \param error ...
!> \par History
!>      none
!> \author Luca Bellucci and Teodoro Laino - created [tlaino] - 01.2008
! *****************************************************************************
  SUBROUTINE dimer_env_create(dimer_env, natom, globenv, dimer_section, error)
    TYPE(dimer_env_type), POINTER            :: dimer_env
    INTEGER, INTENT(IN)                      :: natom
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: dimer_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dimer_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, isize, j, k, n_rep_val, &
                                                stat
    LOGICAL                                  :: explicit, failure
    REAL(KIND=dp)                            :: norm, xval(3)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: array
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: nvec_section

    failure=.FALSE.
    CPPostcondition(.NOT.ASSOCIATED(dimer_env),cp_failure_level,routineP,error,failure)
    ALLOCATE(dimer_env, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    logger => cp_error_get_logger(error)
    dimer_env%ref_count=1
    last_dimer_id=last_dimer_id+1
    dimer_env%id_nr=last_dimer_id
    ! Setup NVEC
    NULLIFY(dimer_env%nvec, dimer_env%rot%g0, dimer_env%rot%g1, dimer_env%rot%g1p,&
         dimer_env%tsl%tls_vec)
    ! Allocate the working arrays
    ALLOCATE(dimer_env%nvec(natom*3),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(dimer_env%rot%g0(natom*3),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(dimer_env%rot%g1(natom*3),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(dimer_env%rot%g1p(natom*3),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ! Check if the dimer vector is available in the input or not..
    nvec_section => section_vals_get_subs_vals(dimer_section,"DIMER_VECTOR",error=error)
    CALL section_vals_get(nvec_section, explicit=explicit, error=error)
    IF (explicit) THEN
       CALL cp_log(logger=logger, level=cp_note_level, fromWhere=routineP,&
            message="Reading Dimer Vector from file!",local=.FALSE.)
       NULLIFY(array)
       CALL section_vals_val_get(nvec_section,"_DEFAULT_KEYWORD_",n_rep_val=n_rep_val,error=error)
       isize = 0
       DO i = 1, n_rep_val
          CALL section_vals_val_get(nvec_section,"_DEFAULT_KEYWORD_",r_vals=array,i_rep_val=i,error=error)
          DO j = 1, SIZE(array)
             isize = isize + 1
             dimer_env%nvec(isize) = array(j)
          END DO
       END DO
       CPPostcondition(isize==SIZE(dimer_env%nvec),cp_failure_level,routineP,error,failure)
    ELSE
       CALL random_numbers(dimer_env%nvec, globenv%gaussian_rng_stream, error)
    END IF
    ! Check for translation in the dimer vector and remove them
    IF (natom>1) THEN
       xval = 0.0_dp
       DO j = 1, natom
          DO k = 1, 3
             i = (j-1)*3+k
             xval(k) = xval(k)+dimer_env%nvec(i)
          END DO
       END DO
       ! Subtract net translations
       xval = xval / REAL(natom*3,KIND=dp)
       DO j = 1, natom
          DO k = 1, 3
             i = (j-1)*3+k
             dimer_env%nvec(i)=dimer_env%nvec(i)-xval(k)
          END DO
       END DO
    END IF
    norm = SQRT(SUM(dimer_env%nvec**2))
    CALL cp_assert(norm>EPSILON(0.0_dp),cp_failure_level,cp_assertion_failed,routineP,&
         "The norm of the dimer vector is 0! Calculation cannot proceed further. "//&
CPSourceFileRef,&
         only_ionode=.TRUE.)
    dimer_env%nvec = dimer_env%nvec / norm
    dimer_env%rot%rotation_step= do_first_rotation_step
    CALL section_vals_val_get(dimer_section,"DR",r_val=dimer_env%dr,error=error)
    CALL section_vals_val_get(dimer_section,"INTERPOLATE_GRADIENT",&
         l_val=dimer_env%rot%interpolate_gradient,error=error)
    CALL section_vals_val_get(dimer_section,"ANGLE_TOLERANCE",&
         r_val=dimer_env%rot%angle_tol,error=error)
    dimer_env%cg_rot%norm_h         = 1.0_dp
    dimer_env%cg_rot%norm_theta     = 0.0_dp
    dimer_env%cg_rot%norm_theta_old = 0.0_dp
    ALLOCATE(dimer_env%cg_rot%nvec_old(natom*3),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE dimer_env_create

! *****************************************************************************
!> \brief ...
!> \param dimer_env ...
!> \param error ...
!> \par History
!>      none
!> \author Luca Bellucci and Teodoro Laino - created [tlaino] - 01.2008
! *****************************************************************************
  SUBROUTINE dimer_env_retain(dimer_env, error)
    TYPE(dimer_env_type), POINTER            :: dimer_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dimer_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(dimer_env),cp_failure_level,routineP,error,failure)
    CPPreconditionNoFail(dimer_env%ref_count>0,cp_failure_level,routineP,error)
    dimer_env%ref_count=dimer_env%ref_count+1
  END SUBROUTINE dimer_env_retain

! *****************************************************************************
!> \brief ...
!> \param dimer_env ...
!> \param error ...
!> \par History
!>      none
!> \author Luca Bellucci and Teodoro Laino - created [tlaino] - 01.2008
! *****************************************************************************
  SUBROUTINE dimer_env_release(dimer_env, error)
    TYPE(dimer_env_type), POINTER            :: dimer_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dimer_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(dimer_env)) THEN
       CPPreconditionNoFail(dimer_env%ref_count>0,cp_failure_level,routineP,error)
       dimer_env%ref_count=dimer_env%ref_count-1
       IF (dimer_env%ref_count==0) THEN
          IF (ASSOCIATED(dimer_env%nvec)) THEN
             DEALLOCATE(dimer_env%nvec, stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(dimer_env%rot%g0)) THEN
             DEALLOCATE(dimer_env%rot%g0, stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(dimer_env%rot%g1)) THEN
             DEALLOCATE(dimer_env%rot%g1, stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(dimer_env%rot%g1p)) THEN
             DEALLOCATE(dimer_env%rot%g1p, stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(dimer_env%cg_rot%nvec_old)) THEN
             DEALLOCATE(dimer_env%cg_rot%nvec_old, stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          ! No need to deallocate tls_vec (just a pointer to aother local array)
          NULLIFY(dimer_env%tsl%tls_vec)
          DEALLOCATE(dimer_env, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
    END IF
  END SUBROUTINE dimer_env_release

END MODULE dimer_types
