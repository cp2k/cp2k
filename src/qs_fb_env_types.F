MODULE qs_fb_env_types

  USE kinds,                           ONLY: dp
  USE qs_fb_atomic_halo_types,         ONLY: fb_atomic_halo_list_associate,&
                                             fb_atomic_halo_list_has_data,&
                                             fb_atomic_halo_list_nullify,&
                                             fb_atomic_halo_list_obj,&
                                             fb_atomic_halo_list_release,&
                                             fb_atomic_halo_list_retain
  USE qs_fb_trial_fns_types,           ONLY: fb_trial_fns_associate,&
                                             fb_trial_fns_has_data,&
                                             fb_trial_fns_nullify,&
                                             fb_trial_fns_obj,&
                                             fb_trial_fns_release,&
                                             fb_trial_fns_retain
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! public types
  PUBLIC :: fb_env_obj

! public methods
  PUBLIC :: fb_env_retain, &
            fb_env_release, &
            fb_env_nullify, &
            fb_env_associate, &
            fb_env_has_data, &
            fb_env_create, &
            fb_env_init, &
            fb_env_get, &
            fb_env_set


  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_fb_env_types'
  INTEGER, PRIVATE, SAVE :: last_fb_env_id = 0


! **********************************************************************
!> \brief wrapper to the simulation parameters used for filtered basis
!>        method
!> \param rcut   : cutoff for included filtered basis set centred at
!>                 each atom. These defines the ranges of the atomic
!>                 halos. rcut(ikind) gives the range for atom of
!>                 global kind ikind
!> \param atomic_halos  : stores information on the neighbors of each
!>                        atom ii, which are defined by rcut
!> \param filter_temperature : parameter controlling the smoothness of
!>                             the filter function during the construction
!>                             of the filter matrix
!> \param auto_cutoff_scale : scale multiplied to max atomic orbital
!>                            radii used for automatic construction of
!>                            rcut
!> \param id_nr : unique id of this object
!> \param ref_count : reference counter of this object
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! **********************************************************************
  TYPE fb_env_data
     INTEGER :: id_nr, ref_count
     REAL(KIND=dp), DIMENSION(:), POINTER :: rcut
     TYPE(fb_atomic_halo_list_obj) :: atomic_halos
     TYPE(fb_trial_fns_obj) :: trial_fns
     REAL(KIND=dp) :: filter_temperature
     REAL(KIND=dp) :: auto_cutoff_scale
  END TYPE fb_env_data


! *****************************************************************************
!> \brief the object container which allows for the creation of an array of
!>        pointers to fb_env
!> \param obj : pointer to a filtered basis environment
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  TYPE fb_env_obj
     TYPE(fb_env_data), POINTER, PRIVATE :: obj
  END TYPE fb_env_obj


CONTAINS


! **********************************************************************
!> \brief retains the given fb_env
!> \param fb_env : the fb_env to retain
!> \param error  : variable to control error logging, stopping, ...
!>                 see module cp_error_handling
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_env_retain(fb_env, error)
    TYPE(fb_env_obj), INTENT(in)             :: fb_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(fb_env%obj), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       CPPreconditionNoFail(fb_env%obj%ref_count>0, cp_failure_level, routineP, error)
       fb_env%obj%ref_count = fb_env%obj%ref_count + 1
    END IF
  END SUBROUTINE fb_env_retain


! **********************************************************************
!> \brief releases a given fb_env
!> \brief ...
!> \param fb_env : the fb_env to release
!> \param error  : variable to control error logging, stopping, ...
!>                 see module cp_error_handling
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_env_release(fb_env, error)
    TYPE(fb_env_obj), INTENT(inout)          :: fb_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(fb_env%obj)) THEN
       CPPreconditionNoFail(fb_env%obj%ref_count>0, cp_failure_level, routineP, error)
       fb_env%obj%ref_count = fb_env%obj%ref_count - 1
       IF (fb_env%obj%ref_count == 0) THEN
          fb_env%obj%ref_count = 1
          IF (ASSOCIATED(fb_env%obj%rcut)) THEN
             DEALLOCATE(fb_env%obj%rcut, STAT=stat)
             CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          END IF
          CALL fb_atomic_halo_list_release(fb_env%obj%atomic_halos, error)
          CALL fb_trial_fns_release(fb_env%obj%trial_fns, error)
          fb_env%obj%ref_count = 0
          DEALLOCATE(fb_env%obj, STAT=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
       END IF
    END IF
    NULLIFY(fb_env%obj)
  END SUBROUTINE fb_env_release


! **********************************************************************
!> \brief nullifies a fb_env object, note that this does not
!>        release the original object. This procedure is used mainly
!>        to nullify the pointer inside the object which is used to
!>        point to the actual data content of the object.
!> \param fb_env : its content must be a NULL fb_env pointer on input,
!>                 and the output returns an empty fb_env object
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_env_nullify(fb_env)
    TYPE(fb_env_obj), INTENT(inout)          :: fb_env

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_env_nullify', &
      routineP = moduleN//':'//routineN

    NULLIFY(fb_env%obj)
  END SUBROUTINE fb_env_nullify


! **********************************************************************
!> \brief Associates one fb_env object to another
!> \param a the fb_env object to be associated  
!> \param b the fb_env object that a is to be associated to  
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_env_associate(a, b)
    TYPE(fb_env_obj), INTENT(out)            :: a
    TYPE(fb_env_obj), INTENT(in)             :: b

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_env_associate', &
      routineP = moduleN//':'//routineN

    a%obj => b%obj
  END SUBROUTINE fb_env_associate


! **********************************************************************
!> \brief Checks if a fb_env object is associated with an actual
!>        data content or not
!> \param fb_env the fb_env object  
!> \retval res ...
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  FUNCTION fb_env_has_data(fb_env) RESULT(res)
    TYPE(fb_env_obj), INTENT(in)             :: fb_env
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_env_has_data', &
      routineP = moduleN//':'//routineN

    res = ASSOCIATED(fb_env%obj)
  END FUNCTION fb_env_has_data


! **********************************************************************
!> \brief creates an empty fb_env object
!> \param fb_env : its content must be a NULL fb_env pointer on input,
!>                 and the output returns an empty fb_env object
!> \param error  : variable to control error logging, stopping, ...
!>                 see module cp_error_handling
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_env_create(fb_env, error)
    TYPE(fb_env_obj), INTENT(inout)          :: fb_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(.NOT.ASSOCIATED(fb_env%obj), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       ALLOCATE(fb_env%obj, STAT=stat)
       CPPostconditionNoFail(stat==0, cp_failure_level, routineP, error)
       NULLIFY(fb_env%obj%rcut)
       CALL fb_atomic_halo_list_nullify(fb_env%obj%atomic_halos)
       CALL fb_trial_fns_nullify(fb_env%obj%trial_fns)
       fb_env%obj%filter_temperature = 0.0_dp
       fb_env%obj%auto_cutoff_scale = 1.0_dp
       fb_env%obj%ref_count = 1
       fb_env%obj%id_nr = last_fb_env_id + 1
       last_fb_env_id = fb_env%obj%id_nr
    END IF
  END SUBROUTINE fb_env_create


! **********************************************************************
!> \brief initialises a fb_env object to become empty
!> \brief ...
!> \param fb_env : the fb_env object, which must not be NULL or
!>                 UNDEFINED upon entry
!> \param error  : variable to control error logging, stopping, ...
!>                 see module cp_error_handling
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_env_init(fb_env, error)
    TYPE(fb_env_obj), INTENT(in)             :: fb_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_env_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(fb_env%obj), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       IF (ASSOCIATED(fb_env%obj%rcut)) THEN
          DEALLOCATE(fb_env%obj%rcut, STAT=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          NULLIFY(fb_env%obj%rcut)
       END IF
       CALL fb_atomic_halo_list_release(fb_env%obj%atomic_halos, error)
       CALL fb_trial_fns_release(fb_env%obj%trial_fns, error)
       fb_env%obj%filter_temperature = 0.0_dp
       fb_env%obj%auto_cutoff_scale = 1.0_dp
    END IF
  END SUBROUTINE fb_env_init


! **********************************************************************
!> \brief method to get attributes from a given fb_env object
!> \brief ...
!> \param fb_env : the fb_env object in question
!> \param rcut   : outputs pointer to rcut attribute of fb_env (optional)
!> \param filter_temperature : outputs pointer to filter_temperature
!>                             attribute of fb_env (optional)
!> \param auto_cutoff_scale  : outputs pointer to auto_cutoff_scale
!>                             attribure of fb_env (optional)
!> \param atomic_halos       : outputs pointer to atomic_halos
!>                             attribute of fb_env (optional)
!> \param trial_fns          : outputs pointer to trial_fns
!>                             attribute of fb_env (optional)
!> \param error  : variable to control error logging, stopping, ...
!>                 see module cp_error_handling
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_env_get(fb_env, &
                        rcut, &
                        filter_temperature, &
                        auto_cutoff_scale, &
                        atomic_halos, &
                        trial_fns, &
                        error)
    TYPE(fb_env_obj), INTENT(in)             :: fb_env
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: rcut
    REAL(KIND=dp), INTENT(out), OPTIONAL     :: filter_temperature, &
                                                auto_cutoff_scale
    TYPE(fb_atomic_halo_list_obj), &
      INTENT(out), OPTIONAL                  :: atomic_halos
    TYPE(fb_trial_fns_obj), INTENT(out), &
      OPTIONAL                               :: trial_fns
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_env_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(fb_env%obj), cp_failure_level, routineP, error, failure)
    CPPrecondition(fb_env%obj%ref_count>0, cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(rcut)) rcut => fb_env%obj%rcut
       IF (PRESENT(filter_temperature)) filter_temperature = fb_env%obj%filter_temperature
       IF (PRESENT(auto_cutoff_scale)) auto_cutoff_scale = fb_env%obj%auto_cutoff_scale
       IF (PRESENT(atomic_halos)) &
            CALL fb_atomic_halo_list_associate(atomic_halos, fb_env%obj%atomic_halos)
       IF (PRESENT(trial_fns)) &
            CALL fb_trial_fns_associate(trial_fns, fb_env%obj%trial_fns)
    END IF
  END SUBROUTINE fb_env_get


! **********************************************************************
!> \brief method to set attributes from a given fb_env object
!> \brief ...
!> \param fb_env : the fb_env object in question
!> \param rcut   : sets rcut attribute of fb_env (optional)
!> \param filter_temperature : sets filter_temperature attribute of fb_env (optional)
!> \param auto_cutoff_scale  : sets auto_cutoff_scale attribute of fb_env (optional)
!> \param atomic_halos       : sets atomic_halos attribute of fb_env (optional)
!> \param trial_fns          : sets trial_fns attribute of fb_env (optional)
!> \param error  : variable to control error logging, stopping, ...
!>                 see module cp_error_handling
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_env_set(fb_env, &
                        rcut, &
                        filter_temperature, &
                        auto_cutoff_scale, &
                        atomic_halos, &
                        trial_fns, &
                        error)
    TYPE(fb_env_obj), INTENT(in)             :: fb_env
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: rcut
    REAL(KIND=dp), INTENT(in), OPTIONAL      :: filter_temperature, &
                                                auto_cutoff_scale
    TYPE(fb_atomic_halo_list_obj), &
      INTENT(in), OPTIONAL                   :: atomic_halos
    TYPE(fb_trial_fns_obj), INTENT(in), &
      OPTIONAL                               :: trial_fns
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_env_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(fb_env%obj), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(rcut)) THEN
          IF (ASSOCIATED(fb_env%obj%rcut)) THEN
             DEALLOCATE(fb_env%obj%rcut, STAT=stat)
             CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          END IF
          fb_env%obj%rcut => rcut
       END IF
       IF (PRESENT(filter_temperature)) &
            fb_env%obj%filter_temperature = filter_temperature
       IF (PRESENT(auto_cutoff_scale)) &
            fb_env%obj%auto_cutoff_scale = auto_cutoff_scale
       IF (PRESENT(atomic_halos)) THEN
          IF (fb_atomic_halo_list_has_data(atomic_halos)) &
               CALL fb_atomic_halo_list_retain(atomic_halos, error)
          CALL fb_atomic_halo_list_release(fb_env%obj%atomic_halos, error)
          CALL fb_atomic_halo_list_associate(fb_env%obj%atomic_halos, atomic_halos)
       END IF
       IF (PRESENT(trial_fns)) THEN
          IF (fb_trial_fns_has_data(trial_fns)) &
               CALL fb_trial_fns_retain(trial_fns, error)
          CALL fb_trial_fns_release(fb_env%obj%trial_fns, error)
          CALL fb_trial_fns_associate(fb_env%obj%trial_fns, trial_fns)
       END IF
    END IF
  END SUBROUTINE fb_env_set

END MODULE qs_fb_env_types
