!----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations        !
!   Copyright (C) 2001  CP2K developers group                                !
!----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/pao_glob_angles_methods [1.0] *
!!
!!   NAME
!!     pao_glob_angles_methods
!!
!!   FUNCTION
!!     metodhs to access and modify the angles that define the unitary
!!     transformation
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     2.2002 created (moved off from pao_types) [fawzi]
!!     5.2002 adapted to new qs_* structure, updated doc [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE pao_glob_angles_methods
  USE cp_array_utils,                  ONLY: cp_1d_logical_guarantee_size,&
                                             cp_1d_r_guarantee_size
  USE cp_b_matrix_structure,           ONLY: cp_b_matrix_struct_type
  USE cp_block_matrix,                 ONLY: cp_block_matrix_p_type,&
                                             cp_block_matrix_type,&
                                             cp_dealloc,&
                                             cp_release,&
                                             cp_retain,&
                                             cp_sp_create
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_level,&
                                             cp_error_get_print_level,&
                                             cp_error_get_stop_level,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_propagate_error,&
                                             cp_error_set,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE pao_types,                       ONLY: cp_get,&
                                             pao_env_type,&
                                             pao_glob_angles_type
  USE particle_list_types,             ONLY: particle_list_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='pao_glob_angles_methods'

! core procedures
  PUBLIC :: cp_init, cp_dealloc_ref, cp_set, cp_get, cp_create, cp_dealloc,&
       cp_did_change

! core functions
  PUBLIC :: cp_validate, cp_valid

! the underlying functions
  PUBLIC :: pao_g_ang_init, pao_g_ang_dealloc_ref, &
       pao_g_ang_create, pao_g_ang_dealloc,&
       pao_g_ang_validate, pao_g_ang_valid,&
       pao_g_ang_set, pao_g_ang_get, pao_g_ang_did_change

! ============== interfaces ===========

! +++ begin of core methods +++

!! initializes the first, argument.
!! Support optional initial values that depend on the 
!! first argument 
  INTERFACE cp_init
     MODULE PROCEDURE pao_g_ang_init
  END INTERFACE

!! Deallocates the memory allocated by the first argument.
!! Supports error=error
  INTERFACE cp_dealloc_ref
     MODULE PROCEDURE pao_g_ang_dealloc_ref
  END INTERFACE

!! Allocates and initialize a new object returning a pointer to it in 
!! the first argument
  INTERFACE cp_create
     MODULE PROCEDURE pao_g_ang_create
  END INTERFACE

!! Deallocates the memory used by the first argument, and the object
!! itself
  INTERFACE cp_dealloc
     MODULE PROCEDURE pao_g_ang_dealloc
  END INTERFACE

!! returns true if the content of the first argumenst is valid.
!! writes the errors ar warnings, full validation (long)
  INTERFACE cp_validate
     MODULE PROCEDURE pao_g_ang_validate
  END INTERFACE

!! returns true if the first argument is valid. Only a fast minimal 
!! validation is performed (no dangling pointers)
  INTERFACE cp_valid
     MODULE PROCEDURE pao_g_ang_valid
  END INTERFACE

!! sets the value of various optonal attributes of the firs argument.
!! the optional attributes depend on the type of the first
  INTERFACE cp_set
     MODULE PROCEDURE pao_g_ang_set
  END INTERFACE

!! gets values from the first argument via optional keyword arguments
!! that depend on the type of the first argument
  INTERFACE cp_get
     MODULE PROCEDURE pao_g_ang_get
  END INTERFACE

!! tells the object that its internal values have changes and its cached
!! values may be invalid
  INTERFACE cp_did_change
     MODULE PROCEDURE pao_g_ang_did_change
  END INTERFACE

! +++ end of the core methods +++

!!***
!****************************************************************************

CONTAINS

! ======= core methods ======

!!****f* pao_glob_angles_methods/pao_g_ang_init [1.0] *
!!
!!   NAME
!!     pao_g_ang_init
!!
!!   SYNOPSIS
!!     Subroutine pao_g_ang_init(glob_angles, qs_env, global_env, angles,&
!!         source, min_overlap_m, min_hamiltonian_m, error)
!!       Type(pao_glob_angles_type), Intent (OUT):: glob_angles
!!       Type(qs_environment_type), Intent (IN), Target:: qs_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Real(KIND = dp), Optional, Dimension(:), Intent (IN):: angles
!!       Type(pao_glob_angles_type), Intent (IN), Optional:: source
!!       Type(cp_block_matrix_p_type), Intent (IN), Optional, Dimension(:)::&
!!         min_overlap_m, min_hamiltonian_m
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_g_ang_init
!!
!!   FUNCTION
!!     initializes the global angles
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     glob_angles: the global angles to be initialized
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE pao_g_ang_init(glob_angles, qs_env, global_env, angles, source,&
       min_overlap_m, min_hamiltonian_m, error)
    TYPE(pao_glob_angles_type), INTENT(out)  :: glob_angles
    TYPE(qs_environment_type), INTENT(in), &
      TARGET                                 :: qs_env
    TYPE(global_environment_type), &
      INTENT(in), TARGET                     :: global_env
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(in), OPTIONAL                   :: angles
    TYPE(pao_glob_angles_type), INTENT(in), &
      OPTIONAL                               :: source
    TYPE(cp_block_matrix_p_type), &
      DIMENSION(:), INTENT(in), OPTIONAL     :: min_overlap_m, &
                                                min_hamiltonian_m
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'pao_glob_angles_methods:pao_g_ang_init'

    LOGICAL                                  :: failure

    failure=.FALSE.

    NULLIFY(glob_angles%angles, glob_angles%unitary_t,&
         glob_angles%ut_calculated, glob_angles%min_overlap_m, &
         glob_angles%min_hamiltonian_m, glob_angles%NUi_injection_calculated,&
         glob_angles%NUi_injection)
    glob_angles%m_s_m_valid=.FALSE.
    glob_angles%m_h_m_valid=.FALSE.
    glob_angles%initialized = .TRUE.
    IF (.not.failure) THEN
       CALL pao_g_ang_set(glob_angles,qs_env=qs_env, global_env=global_env,&
            angles=angles,&
            source=source, min_hamiltonian_m=min_hamiltonian_m,&
            min_overlap_m=min_overlap_m, error=error)
    END IF
  END SUBROUTINE pao_g_ang_init
!***************************************************************************

!!****f* pao_glob_angles_methods/pao_g_ang_dealloc_ref [1.0] *
!!
!!   NAME
!!     pao_g_ang_dealloc_ref
!!
!!   SYNOPSIS
!!     Subroutine pao_g_ang_dealloc_ref(glob_angles, global_env, error)
!!       Type(pao_glob_angles_type), Intent (INOUT):: glob_angles
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_g_ang_dealloc_ref
!!
!!   FUNCTION
!!     deallocates the memory allocated by a glob_angles object.
!!     does NOT deallocate the object itself
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     pao_glob_angles: the object to release
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE pao_g_ang_dealloc_ref(glob_angles,global_env, error)
    TYPE(pao_glob_angles_type), &
      INTENT(inout)                          :: glob_angles
    TYPE(global_environment_type), &
      INTENT(in), TARGET                     :: global_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'pao_glob_angles_methods:pao_g_ang_dealloc_ref'

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(glob_angles%initialized, cp_failure_level,routineP,error,failure)
    failureIf: IF (.NOT. failure) THEN
       IF (ASSOCIATED(glob_angles%angles)) THEN
          DEALLOCATE(glob_angles%angles, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(glob_angles%ut_calculated)) THEN
          DEALLOCATE(glob_angles%ut_calculated,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(glob_angles%NUi_injection_calculated)) THEN
          DEALLOCATE(glob_angles%NUi_injection_calculated)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
       CALL cp_release(glob_angles%NUi_injection,error=error)
       CALL cp_release(glob_angles%unitary_t,error=error)
       CALL cp_dealloc(glob_angles%min_hamiltonian_m,error=error)
       CALL cp_dealloc(glob_angles%min_overlap_m,error=error)
    END IF failureIf
  END SUBROUTINE pao_g_ang_dealloc_ref
!***************************************************************************

!!****f* pao_glob_angles_methods/pao_g_ang_dealloc [1.0] *
!!
!!   NAME
!!     pao_g_ang_dealloc
!!
!!   SYNOPSIS
!!     Subroutine pao_g_ang_dealloc(glob_angles, global_env, error)
!!       Type(pao_glob_angles_type), Pointer:: glob_angles
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_g_ang_dealloc
!!
!!   FUNCTION
!!     deallocated the memory allocated by this object and the oject itself
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     glob_angles: the object to be deallocated
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE pao_g_ang_dealloc(glob_angles,global_env,error)
    TYPE(pao_glob_angles_type), POINTER      :: glob_angles
    TYPE(global_environment_type), &
      INTENT(in), TARGET                     :: global_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'pao_glob_angles_methods:pao_g_ang_dealloc'

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(glob_angles)) THEN
       CPPrecondition(glob_angles%initialized,cp_failure_level,routineP,error,failure)
       IF (.NOT. failure) THEN
          CALL pao_g_ang_dealloc_ref(glob_angles, global_env=global_env,&
               error=error)
       END IF
       DEALLOCATE(glob_angles,stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    END IF
  END SUBROUTINE pao_g_ang_dealloc
!***************************************************************************

!!****f* pao_glob_angles_methods/pao_g_ang_create [1.0] *
!!
!!   NAME
!!     pao_g_ang_create
!!
!!   SYNOPSIS
!!     Subroutine pao_g_ang_create(glob_angles, qs_env, global_env, angles,&
!!         source, min_overlap_m, min_hamiltonian_m, error)
!!       Type(pao_glob_angles_type), Pointer:: glob_angles
!!       Type(qs_environment_type), Intent (IN), Target:: qs_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Real(KIND = dp), Optional, Dimension(:), Intent (IN):: angles
!!       Type(pao_glob_angles_type), Intent (IN), Optional:: source
!!       Type(cp_block_matrix_p_type), Intent (IN), Dimension(:), Optional::&
!!         min_hamiltonian_m, min_overlap_m
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_g_ang_create
!!
!!   FUNCTION
!!     creates a new glob_angles object
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     glob_angles: will be a pointer to the newly created object
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!     for the others inputs see pao_g_ang_init
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE pao_g_ang_create(glob_angles, qs_env, global_env, angles, source,&
       min_overlap_m, min_hamiltonian_m, error)
    TYPE(pao_glob_angles_type), POINTER      :: glob_angles
    TYPE(qs_environment_type), INTENT(in), &
      TARGET                                 :: qs_env
    TYPE(global_environment_type), &
      INTENT(in), TARGET                     :: global_env
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(in), OPTIONAL                   :: angles
    TYPE(pao_glob_angles_type), INTENT(in), &
      OPTIONAL                               :: source
    TYPE(cp_block_matrix_p_type), &
      DIMENSION(:), INTENT(in), OPTIONAL     :: min_overlap_m, &
                                                min_hamiltonian_m
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'pao_glob_angles_methods:pao_g_ang_create'

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    ALLOCATE(glob_angles, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL pao_g_ang_init(glob_angles, qs_env=qs_env, global_env=global_env,&
            angles=angles, source=source,&
            min_overlap_m=min_overlap_m, min_hamiltonian_m=min_hamiltonian_m,&
            error=error)
    END IF
  END SUBROUTINE pao_g_ang_create
!***************************************************************************

!!****f* pao_glob_angles_methods/pao_g_ang_validate [1.0] *
!!
!!   NAME
!!     pao_g_ang_validate
!!
!!   SYNOPSIS
!!     Function pao_g_ang_validate(glob_angles, qs_env, global_env, error)&
!!         Result(res)
!!       Logical:: res
!!       Type(pao_glob_angles_type), Intent (IN):: glob_angles
!!       Type(qs_environment_type), Intent (IN), Target:: qs_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function pao_g_ang_validate
!!
!!   FUNCTION
!!     returns true if the angles are valid
!!
!!   NOTES
!!     Writes out all the failures as warning, use error to change
!!     this behaviour
!!
!!   INPUTS
!!     glob_angles: the gangles to be tested
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  FUNCTION pao_g_ang_validate(glob_angles,qs_env,global_env, error) RESULT(res)
    TYPE(pao_glob_angles_type), INTENT(in)   :: glob_angles
    TYPE(qs_environment_type), INTENT(in), &
      TARGET                                 :: qs_env
    TYPE(global_environment_type), &
      INTENT(in), TARGET                     :: global_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'pao_glob_angles_methods:pao_g_ang_validate'

    LOGICAL                                  :: failure
    TYPE(pao_env_type), POINTER              :: pao_env

    failure=.FALSE.
    CPPrecondition(glob_angles%initialized,cp_failure_level,routineP,error,failure)
    failureIf: IF (.NOT. failure) THEN
       CPAssert(ASSOCIATED(glob_angles%angles),cp_warning_level,routineP,error,failure)
       CALL get_qs_env(qs_env,pao_env=pao_env)
       CPAssert(ASSOCIATED(pao_env),cp_warning_level,routineP,error,failure)
       IF (.not.failure) THEN
          CALL cp_assert(pao_env%angles_begin_at &
               (SIZE(pao_env%angles_begin_at))-1== &
               SIZE(glob_angles%angles),cp_warning_level,&
               cp_assertion_failed, routineP,&
               "number of angles different from last begin_at in "//&
               CPSourceFileRef,&
               error,failure)
       END IF
!!FM       if (associated(glob_angles%ut_calculated)) then
!!FM          call cp_assert(size(glob_angles%ut_calculated)==&
!!FM               glob_angles%pao_env%n_atoms, cp_warning_level,&
!!FM               cp_assertion_failed, routineP,&
!!FM               "cached invalid ut_calculated in "//&
!!FM               CPSourceFileRef,&
!!FM               error,failure)
!!FM       end if
!!FM       if (associated(glob_angles%unitary_t)) then
!!FM          call cp_assert(cp_validate(glob_angles%unitary_t,&
!!FM               error=error), cp_warning_level,&
!!FM               cp_assertion_failed, routineP,&
!!FM               "cached invalid unitary_t in "//&
!!FM               CPSourceFileRef,&
!!FM               error,failure)
!!FM          CPAssert(associated(glob_angles%ut_calculated),cp_warning_level,routineP,error,failure)
!!FM       end if
       IF (ASSOCIATED(glob_angles%min_hamiltonian_m)) THEN
!!FM          call cp_assert(cp_validate(glob_angles%min_hamiltonian_m,&
!!FM               error=error), cp_warning_level,&
!!FM               cp_assertion_failed, routineP,&
!!FM               "cached invalid min_hamiltonian_m in "//&
!!FM               CPSourceFileRef,&
!!FM               error,failure)
       ELSE
          CPAssert(.not.glob_angles%m_h_m_valid,cp_warning_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(glob_angles%min_overlap_m)) THEN
!!FM          call cp_assert(cp_validate(glob_angles%min_overlap_m,&
!!FM               error=error), cp_warning_level,&
!!FM               cp_assertion_failed, routineP,&
!!FM               "cached invalid min_overlap_m in "//&
!!FM               CPSourceFileRef,&
!!FM               error,failure)
       ELSE
          CPAssert(.not.glob_angles%m_s_m_valid,cp_warning_level,routineP,error,failure)
       END IF
    END IF failureIf
    res=.not.failure
  END FUNCTION pao_g_ang_validate
!***************************************************************************

!!****f* pao_glob_angles_methods/pao_g_ang_set [1.0] *
!!
!!   NAME
!!     pao_g_ang_set
!!
!!   SYNOPSIS
!!     Subroutine pao_g_ang_set(glob_angles, qs_env, global_env, angles,&
!!         source, min_overlap_m, min_hamiltonian_m, m_s_m_valid,&
!!         m_h_m_valid, angles_val, error)
!!       Type(pao_glob_angles_type), Intent (INOUT):: glob_angles
!!       Type(qs_environment_type), Intent (IN), Target:: qs_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Real(KIND = dp), Optional, Dimension(:), Intent (IN):: angles
!!       Type(pao_glob_angles_type), Optional, Intent (IN):: source
!!       Type(cp_block_matrix_p_type), Optional, Intent (IN), Dimension(:)::&
!!         min_overlap_m
!!       Type(cp_block_matrix_p_type), Dimension(:), Intent (IN), Optional::&
!!         min_hamiltonian_m
!!       Logical, Intent (IN), Optional:: m_h_m_valid, m_s_m_valid
!!       Real(KIND = dp), Intent (IN), Optional:: angles_val
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_g_ang_set
!!
!!   FUNCTION
!!     can change the value of the various components of pao_glob_angles
!!
!!   NOTES
!!     you should almost never need this function
!!
!!   INPUTS
!!     glob_angles: the structure wose values are set
!!     pao_env: the pao env that this glob_angles refers to (SHARED)
!!     angles: the angles themselves as a long 1d array (copied)
!!     source: a glob angles that is copied into this (other arguments
!!             have the precedence for setting the values, the values
!!             of source are taken as default, if the value is not directly
!!             set)
!!     ut_calculated: flags to see which blocka have been calculated
!!     unitary_t: the cached unitary matrix. Be careful: only the 
!!             blocks that have been calculated are valid.
!!             Every processor calculates oly the blocks he needs.
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE pao_g_ang_set(glob_angles, qs_env, global_env, angles,source,&
       min_overlap_m, min_hamiltonian_m,&
       m_s_m_valid, m_h_m_valid, angles_val, error)
    TYPE(pao_glob_angles_type), &
      INTENT(inout)                          :: glob_angles
    TYPE(qs_environment_type), INTENT(in), &
      TARGET                                 :: qs_env
    TYPE(global_environment_type), &
      INTENT(in), TARGET                     :: global_env
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(in), OPTIONAL                   :: angles
    TYPE(pao_glob_angles_type), INTENT(in), &
      OPTIONAL                               :: source
    TYPE(cp_block_matrix_p_type), &
      DIMENSION(:), INTENT(in), OPTIONAL     :: min_overlap_m, &
                                                min_hamiltonian_m
    LOGICAL, INTENT(in), OPTIONAL            :: m_s_m_valid, m_h_m_valid
    REAL(KIND=dp), INTENT(in), OPTIONAL      :: angles_val
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'pao_glob_angles_methods:pao_g_ang_set'

    INTEGER                                  :: stat
    LOGICAL                                  :: failure
    TYPE(cp_error_type)                      :: iError
    TYPE(pao_env_type), POINTER              :: pao_env

    CALL cp_error_init(iError,template_error=error)
    failure=.FALSE.

    CALL get_qs_env(qs_env, pao_env=pao_env)
    CPPrecondition(glob_angles%initialized,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(pao_env),cp_failure_level,routineP,error,failure)
    failureIf: IF (.NOT. failure) THEN
       IF (PRESENT(source)) THEN
          CPPrecondition(source%initialized,cp_failure_level,routineP,error,failure)
          IF (.not.failure) THEN
             CALL cp_error_set(iError, print_level=MAX(cp_failure_level,&
                  cp_error_get_print_level(iError)),&
                  stop_level=MAX(cp_failure_level,&
                  cp_error_get_stop_level(iError)))
             IF (ASSOCIATED(source%angles)) THEN
                CALL cp_1d_r_guarantee_size(glob_angles%angles,&
                     SIZE(source%angles),init=0.0_dp,&
                     error=iError)
                CALL cp_assert(cp_error_get_level(iError)<cp_failure_level,&
                     cp_error_get_level(iError),cp_assertion_failed,routineP,&
                     "sub routine cp_guarantee_size error in"//&
                     CPSourceFileRef,&
                     error,failure)
                glob_angles%angles=source%angles
             ELSE IF (ASSOCIATED(glob_angles%angles)) THEN
                DEALLOCATE(glob_angles%angles, stat=stat)
                CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
             END IF
             CALL pao_g_angles_set_unitary_t(glob_angles,&
                  ut_calculated=source%ut_calculated,&
                  new_value=source%unitary_t,&
                  error=error)
             CALL pao_g_angles_set_min_s_m(glob_angles,&
                  new_value=source%min_overlap_m,&
                  error=error)
             CALL pao_g_angles_set_min_h_m(glob_angles,&
                  new_value=source%min_hamiltonian_m,&
                  error=error)
             glob_angles%m_h_m_valid=source%m_h_m_valid
             glob_angles%m_s_m_valid=source%m_s_m_valid
             CALL cp_error_dealloc_ref(iError,error)
             CALL cp_error_init(iError, template_error=error)
          END IF
       END IF
       IF (PRESENT(angles)) THEN
          CALL cp_1d_r_guarantee_size(glob_angles%angles,SIZE(angles),&
               init=0.0_dp,error=iError)
          CALL cp_assert(cp_error_get_level(iError)<cp_failure_level,&
               cp_failure_level,cp_assertion_failed,routineP,&
               "error in cp_guarantee_size called in "//&
               CPSourceFileRef,&
               error,failure)
          IF (.NOT. failure) THEN
             glob_angles%angles = angles
          END IF
       END IF
       IF (PRESENT(angles_val)) THEN
          CALL cp_1d_r_guarantee_size(glob_angles%angles,&
               pao_env%angles_begin_at(SIZE(pao_env%angles_begin_at)),&
               init=0.0_dp,error=iError)
          CALL cp_error_propagate_error(iError,routineP,&
               message="error in cp_guarantee_size called in "//&
               CPSourceFileRef,&
               error=error,failure_level=cp_failure_level,failure=failure)
          glob_angles%angles(:)=angles_val
       END IF
       IF (PRESENT(min_overlap_m)) THEN
          CALL pao_g_angles_set_min_s_m(glob_angles,&
               new_value=min_overlap_m,&
               error=error)
       END IF
       IF (PRESENT(min_hamiltonian_m)) THEN
          CALL pao_g_angles_set_min_h_m(glob_angles,&
               new_value=min_hamiltonian_m,&
               error=error)
       END IF
       IF (.not.failure) THEN
          CALL cp_1d_r_guarantee_size(glob_angles%angles,&
               pao_env%angles_begin_at(SIZE(pao_env%angles_begin_at)),&
               init=0.0_dp,error=error)
       END IF
    END IF failureIf
    CALL cp_error_dealloc_ref(iError,error=error)
  END SUBROUTINE pao_g_ang_set
!***************************************************************************

!!****f* pao_glob_angles_methods/pao_g_ang_get [1.0] *
!!
!!   NAME
!!     pao_g_ang_get
!!
!!   SYNOPSIS
!!     Subroutine pao_g_ang_get(glob_angles, qs_env, global_env, angles,&
!!         m_s_m_valid, m_h_m_valid, ut_calculated, unitary_t,&
!!         nui_injection_calculated, nui_injection, error)
!!       Type(pao_glob_angles_type), Intent (INOUT):: glob_angles
!!       Type(qs_environment_type), Intent (IN), Target:: qs_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Real(KIND = dp), Optional, Dimension(:), Pointer:: angles
!!       Logical, Intent (OUT), Optional:: m_h_m_valid, m_s_m_valid
!!       Logical, Dimension(:), Pointer, Optional:: ut_calculated,&
!!         nui_injection_calculated
!!       Type(cp_block_matrix_type), Pointer, Optional:: unitary_t,&
!!         nui_injection
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_g_ang_get
!!
!!   FUNCTION
!!     retrive the various attributes of the glob angles object
!!
!!   NOTES
!!     try not to use the angles attribute
!!
!!   INPUTS
!!     - glob_angles: the object you get the info from
!!     - ut_calculated: which blocks of unitary_t are valid . Allocated
!!       if necessary.
!!     - unitary_t: the unitary transformation (be careful, only the
!!       blocks flagged in ut_calculated are valid). Allocated if 
!!       necessary.
!!     - NUi_injection_calculated: which blocks of the NUi_injection
!!       are valid. Allocated if necessary.
!!     - NUi_injection: the projection that goes from the minimal
!!       basis to the full basis. N is an atomic orthogonalisation term
!!       U is the unitary transformation, and "i" is the basis injection.
!!       Only the blocks that have true in NUi_injection_calculated are
!!       valid.. Allocated if necessary.
!!
!!     for the other variables see the pao_glob_angles_type in the
!!     pao_types module
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE pao_g_ang_get(glob_angles, qs_env, global_env, angles,&
       m_s_m_valid, m_h_m_valid, ut_calculated, unitary_t, &
       NUi_injection_calculated, NUi_injection, error)
    TYPE(pao_glob_angles_type), &
      INTENT(inout)                          :: glob_angles
    TYPE(qs_environment_type), INTENT(in), &
      TARGET                                 :: qs_env
    TYPE(global_environment_type), &
      INTENT(in), TARGET                     :: global_env
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: angles
    LOGICAL, INTENT(out), OPTIONAL           :: m_s_m_valid, m_h_m_valid
    LOGICAL, DIMENSION(:), OPTIONAL, POINTER :: ut_calculated
    TYPE(cp_block_matrix_type), OPTIONAL, &
      POINTER                                :: unitary_t
    LOGICAL, DIMENSION(:), OPTIONAL, POINTER :: NUi_injection_calculated
    TYPE(cp_block_matrix_type), OPTIONAL, &
      POINTER                                :: NUi_injection
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'pao_glob_angles_methods:pao_g_ang_get'

    INTEGER                                  :: n_atoms, stat
    LOGICAL                                  :: failure
    TYPE(cp_b_matrix_struct_type), POINTER   :: m_struct
    TYPE(cp_error_type)                      :: iError
    TYPE(cp_subsystem_type), POINTER         :: subsys
    TYPE(pao_env_type), POINTER              :: pao_env
    TYPE(particle_list_type), POINTER        :: particles

    NULLIFY(subsys, particles, pao_env, m_struct)
    CALL cp_error_init(iError,template_error=error)
    failure=.FALSE.

    CPPrecondition(glob_angles%initialized,cp_failure_level,routineP,error,failure)
    CALL get_qs_env(qs_env,subsys=subsys,error=error)
    CALL cp_subsys_get(subsys, particles=particles, error=error)
    n_atoms=particles%n_els
    failureIf: IF (.NOT. failure) THEN
       IF (PRESENT(angles)) angles => glob_angles%angles
       IF (PRESENT(m_s_m_valid)) m_s_m_valid = glob_angles%m_s_m_valid
       IF (PRESENT(m_h_m_valid)) m_h_m_valid=glob_angles%m_h_m_valid
       IF (PRESENT(unitary_t)) THEN
          CALL get_qs_env(qs_env, pao_env=pao_env)
          CALL cp_get(pao_env, full_bas_diag_m_struct=m_struct)
          IF (.not.ASSOCIATED(glob_angles%unitary_t)) THEN
             CALL cp_sp_create(glob_angles%unitary_t,&
                  matrix_struct=m_struct,&
                  error=error)
          END IF
          IF (.not.ASSOCIATED(glob_angles%ut_calculated)) THEN
             ALLOCATE(glob_angles%ut_calculated(n_atoms),&
                  stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             IF (stat==0) glob_angles%ut_calculated=.FALSE.
          END IF
          unitary_t => glob_angles%unitary_t
       END IF
       IF (PRESENT(ut_calculated).AND.&
            .not.ASSOCIATED(glob_angles%ut_calculated)) THEN
          ALLOCATE(glob_angles%ut_calculated(n_atoms),&
               stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          IF (stat==0) glob_angles%ut_calculated=.FALSE.
          ut_calculated => glob_angles%ut_calculated
       END IF
       IF (PRESENT(NUi_injection)) THEN
          CALL get_qs_env(qs_env, pao_env=pao_env)
          CALL cp_get(pao_env, full_min_bas_m_struct=m_struct)
          IF (.not.ASSOCIATED(glob_angles%NUi_injection)) THEN
             CALL cp_sp_create(glob_angles%NUi_injection,&
                  matrix_struct=m_struct,&
                  error=error)
          END IF
          IF (.not.ASSOCIATED(glob_angles%NUi_injection_calculated)) THEN
             ALLOCATE(glob_angles%NUi_injection_calculated(n_atoms),&
                  stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             IF (stat==0) glob_angles%NUi_injection_calculated=.FALSE.
          END IF
          NUi_injection => glob_angles%NUi_injection
       END IF
       IF (PRESENT(NUi_injection_calculated).AND.&
            .not.ASSOCIATED(glob_angles%NUi_injection_calculated)) THEN
          ALLOCATE(glob_angles%NUi_injection_calculated(n_atoms),&
               stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          IF (stat==0) glob_angles%NUi_injection_calculated=.FALSE.
          NUi_injection_calculated => glob_angles%NUi_injection_calculated
       END IF
    END IF failureIf
  END SUBROUTINE pao_g_ang_get
!***************************************************************************

!!****f* pao_glob_angles_methods/pao_g_ang_valid [1.0] *
!!
!!   NAME
!!     pao_g_ang_valid
!!
!!   SYNOPSIS
!!     Function pao_g_ang_valid(glob_angles, qs_env, global_env, error)&
!!         Result(res)
!!       Logical:: res
!!       Type(pao_glob_angles_type), Intent (IN):: glob_angles
!!       Type(qs_environment_type), Intent (IN), Target:: qs_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function pao_g_ang_valid
!!
!!   FUNCTION
!!     minimal validation on the glob angles object (no dangling pointers)
!!
!!   NOTES
!!     do not check pao_env?
!!
!!   INPUTS
!!     - glob_angles: the globa angles object to test
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  FUNCTION pao_g_ang_valid(glob_angles, qs_env, global_env, error) RESULT(res)
    TYPE(pao_glob_angles_type), INTENT(in)   :: glob_angles
    TYPE(qs_environment_type), INTENT(in), &
      TARGET                                 :: qs_env
    TYPE(global_environment_type), &
      INTENT(in), TARGET                     :: global_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'pao_glob_angles_methods:pao_g_ang_valid'

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(glob_angles%initialized,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (.NOT.(ASSOCIATED(glob_angles%angles).AND.&
            (.not.ASSOCIATED(glob_angles%unitary_t).OR.&
            ASSOCIATED(glob_angles%ut_calculated)) .AND.&
            (.not.glob_angles%m_s_m_valid .OR.&
            ASSOCIATED(glob_angles%min_overlap_m)).AND.&
            (.not.glob_angles%m_h_m_valid .OR.&
            ASSOCIATED(glob_angles%min_hamiltonian_m))&
            )) THEN
          CPAssert(ASSOCIATED(glob_angles%angles),cp_warning_level,routineP,error,failure)
          IF (ASSOCIATED(glob_angles%unitary_t)) THEN
             CPAssert(ASSOCIATED(glob_angles%ut_calculated),cp_warning_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(glob_angles%NUi_injection_calculated)) THEN
             CPAssert(ASSOCIATED(glob_angles%NUi_injection),cp_warning_level,routineP,error,failure)
          END IF
          IF (glob_angles%m_s_m_valid) THEN
             CPAssert(ASSOCIATED(glob_angles%min_overlap_m),cp_warning_level,routineP,error,failure)
          END IF
          IF (glob_angles%m_h_m_valid) THEN
             CPAssert(ASSOCIATED(glob_angles%min_hamiltonian_m),cp_warning_level,routineP,error,failure)
          END IF
       END IF
    END IF
    res=.not.failure
  END FUNCTION pao_g_ang_valid
!***************************************************************************

! ===== special methods ====

!!****f* pao_glob_angles_methods/pao_g_angles_set_unitary_t [1.0] *
!!
!!   NAME
!!     pao_g_angles_set_unitary_t
!!
!!   FUNCTION
!!     sets the unitary transformation
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     glob_angles: the structure wose values are set
!!     new_value: the new cached unitary matrix.
!!     ut_calculated: which blocks are valid
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE pao_g_angles_set_unitary_t(glob_angles,&
       new_value, ut_calculated,error)
    TYPE(pao_glob_angles_type), &
      INTENT(inout)                          :: glob_angles
    TYPE(cp_block_matrix_type), INTENT(in), &
      OPTIONAL, TARGET                       :: new_value
    LOGICAL, DIMENSION(:), INTENT(in), &
      OPTIONAL, TARGET                       :: ut_calculated
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'pao_glob_angles_methods:pao_g_angles_set_unitary_t'

    LOGICAL                                  :: failure
    LOGICAL, DIMENSION(:), POINTER           :: ut_c
    TYPE(cp_block_matrix_type), POINTER      :: n_value
    TYPE(cp_error_type)                      :: iError

    CALL cp_error_init(iError, template_error=error, &
         print_level=MAX(cp_failure_level, cp_error_get_print_level(error)))
    failure=.FALSE.

    IF (PRESENT(ut_calculated)) THEN
       ut_c => ut_calculated
       IF (ASSOCIATED(ut_c)) THEN
          CALL cp_1d_logical_guarantee_size(glob_angles%ut_calculated,&
               SIZE(ut_calculated),error=iError)
          CALL cp_error_propagate_error(iError,routineP,&
               CPSourceFileRef,&
               error,failure,cp_failure_level)
          glob_angles%ut_calculated=ut_calculated
       ELSE IF (ASSOCIATED(glob_angles%ut_calculated)) THEN
          DEALLOCATE(glob_angles%ut_calculated)
       END IF
    END IF
    IF (PRESENT(new_value)) THEN
       n_value => new_value
       IF (ASSOCIATED(n_value)) THEN
          CALL cp_retain(n_value,error=error)
          CALL cp_release(glob_angles%unitary_t, error=error)
          glob_angles%unitary_t => new_value
       ELSE IF (ASSOCIATED(glob_angles%unitary_t)) THEN
! do not dealloc ??
          CALL cp_release(glob_angles%unitary_t,error=error)
       END IF
    END IF
  END SUBROUTINE pao_g_angles_set_unitary_t
!***************************************************************************

!!****f* pao_glob_angles_methods/pao_g_angles_set_min_s_m [1.0] *
!!
!!   NAME
!!     pao_g_angles_set_min_s_m
!!
!!   FUNCTION
!!     sets the unitary transformation
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     glob_angles: the structure wose values are set
!!     min_overlap_m: the overlap matrix in the minimal basis
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE pao_g_angles_set_min_s_m(glob_angles, new_value,error)
    TYPE(pao_glob_angles_type), &
      INTENT(inout)                          :: glob_angles
    TYPE(cp_block_matrix_p_type), &
      DIMENSION(:), INTENT(in), OPTIONAL     :: new_value
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'pao_glob_angles_methods:pao_g_angles_set_min_s_m'

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure
    TYPE(cp_block_matrix_type), POINTER      :: n_val

    failure=.FALSE.

    IF (PRESENT(new_value)) THEN
       glob_angles%m_s_m_valid=.TRUE.
       DO i=1,SIZE(new_value)
          n_val => new_value(i)%matrix
          CALL cp_retain(n_val, error=error)
       END DO
       IF (ASSOCIATED(glob_angles%min_overlap_m)) THEN
          CALL cp_dealloc(glob_angles%min_overlap_m,error=error)
       END IF
       ALLOCATE(glob_angles%min_overlap_m(SIZE(new_value)),stat=stat)
       DO i=1,SIZE(new_value)
          glob_angles%min_overlap_m(i)%matrix => new_value(i)%matrix
       END DO
    ELSE IF (ASSOCIATED(glob_angles%min_overlap_m)) THEN
! do not dealloc ??
       CALL cp_dealloc(glob_angles%min_overlap_m,error=error)
    END IF
  END SUBROUTINE pao_g_angles_set_min_s_m
!***************************************************************************

!!****f* pao_glob_angles_methods/pao_g_angles_set_min_h_m [1.0] *
!!
!!   NAME
!!     pao_g_angles_set_min_h_m
!!
!!   FUNCTION
!!     sets the unitary transformation
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     glob_angles: the structure wose values are set
!!     min_hamiltonian_m: the hamiltonian matrix in the minimal basis
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE pao_g_angles_set_min_h_m(glob_angles, new_value,error)
    TYPE(pao_glob_angles_type), &
      INTENT(inout)                          :: glob_angles
    TYPE(cp_block_matrix_p_type), &
      DIMENSION(:), INTENT(in), OPTIONAL     :: new_value
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_g_angles_set_min_h_m', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure
    TYPE(cp_block_matrix_type), POINTER      :: n_val

    failure=.FALSE.

    IF (PRESENT(new_value)) THEN
       glob_angles%m_s_m_valid=.TRUE.
       DO i=1,SIZE(new_value)
          n_val => new_value(i)%matrix
          CALL cp_retain(n_val, error=error)
       END DO
       IF (ASSOCIATED(glob_angles%min_overlap_m)) THEN
          CALL cp_dealloc(glob_angles%min_overlap_m,error=error)
       END IF
       ALLOCATE(glob_angles%min_overlap_m(SIZE(new_value)),stat=stat)
       DO i=1,SIZE(new_value)
          glob_angles%min_overlap_m(i)%matrix => new_value(i)%matrix
       END DO
    ELSE IF (ASSOCIATED(glob_angles%min_overlap_m)) THEN
! do not dealloc ??
       CALL cp_dealloc(glob_angles%min_overlap_m,error=error)
    END IF
  END SUBROUTINE pao_g_angles_set_min_h_m
!***************************************************************************

!!****f* pao_glob_angles_methods/pao_g_ang_did_change [1.0] *
!!
!!   NAME
!!     pao_g_ang_did_change
!!
!!   SYNOPSIS
!!     Subroutine pao_g_ang_did_change(glob_angles, qs_env, global_env,&
!!         h_changed, s_changed, angles_changed, full_reset, error)
!!       Type(pao_glob_angles_type), Intent (INOUT), Target:: glob_angles
!!       Type(qs_environment_type), Intent (IN), Target:: qs_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Logical, Intent (IN), Optional:: h_changed, s_changed,&
!!         angles_changed, full_reset
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_g_ang_did_change
!!
!!   FUNCTION
!!     Invalidate all the cached values.
!!     Should be called when the value of the angles change.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     glob_angles: the globa angles object that should reset his cache
!!     full_reset: if true deallocates all the memory of the cache
!!            nedeed when the matrix structure changes (i.e. when the 
!!            atoms move). defaults to false
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE pao_g_ang_did_change(glob_angles,qs_env,global_env,&
       h_changed, s_changed,&
       angles_changed, full_reset, error)
    TYPE(pao_glob_angles_type), &
      INTENT(inout), TARGET                  :: glob_angles
    TYPE(qs_environment_type), INTENT(in), &
      TARGET                                 :: qs_env
    TYPE(global_environment_type), &
      INTENT(in), TARGET                     :: global_env
    LOGICAL, INTENT(in), OPTIONAL            :: h_changed, s_changed, &
                                                angles_changed, full_reset
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_g_ang_did_change', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: a_c, failure, full_r, h_c, s_c

!call timeset(routineN,'I',moduleN,handle)

    failure=.FALSE.
    h_c=.TRUE. ; s_c=.FALSE. ; a_c=.TRUE. ; full_r=.FALSE.

    CPPrecondition(glob_angles%initialized,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(h_changed)) h_c=h_changed
       IF (PRESENT(s_changed)) s_c=s_changed
       IF (PRESENT(angles_changed)) a_c=angles_changed
       IF (PRESENT(full_reset)) full_r=full_reset
       glob_angles%m_h_m_valid=glob_angles%m_h_m_valid.and.&
            .NOT.(h_c.or.s_c.or.a_c.or.full_r)
       glob_angles%m_s_m_valid=glob_angles%m_s_m_valid.and.&
            .NOT.(s_c.or.a_c.or.full_r)
       IF (ASSOCIATED(glob_angles%ut_calculated).and.a_c) THEN
          glob_angles%ut_calculated=.FALSE.
       END IF
       IF (ASSOCIATED(glob_angles%NUi_injection_calculated).and.a_c) THEN
          glob_angles%NUi_injection_calculated=.FALSE.
       END IF
       IF (full_r.or.s_c) THEN
          CALL cp_release(glob_angles%unitary_t,error=error)
          CALL cp_release(glob_angles%NUi_injection,error=error)
          CALL cp_dealloc(glob_angles%min_hamiltonian_m,error=error)
          CALL cp_dealloc(glob_angles%min_overlap_m)
          IF (ASSOCIATED(glob_angles%ut_calculated)) THEN
             IF (full_r) THEN
                DEALLOCATE(glob_angles%ut_calculated, stat=stat)
                CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
             ELSE
                glob_angles%ut_calculated=.FALSE.
             END IF
          END IF
          IF (ASSOCIATED(glob_angles%NUi_injection_calculated)) THEN
             IF (full_r) THEN
                DEALLOCATE(glob_angles%NUi_injection_calculated,stat=stat)
                CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
             ELSE
                glob_angles%NUi_injection_calculated=.FALSE.
             END IF
          END IF
       END IF
    END IF
    !call timestop(0.0_dp,handle)
  END SUBROUTINE pao_g_ang_did_change
!***************************************************************************

END MODULE pao_glob_angles_methods
