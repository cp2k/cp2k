!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/pao_l_angles_methods [1.0] *
!!
!!   NAME
!!     pao_l_angles_methods
!!
!!   FUNCTION
!!     methods to access and modify the handle to the local angles
!!
!!   NOTES
!!     get methods with inout could be changed to in by using a ptr var 
!!     to modify the subcomponents
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!*****************************************************************************
MODULE pao_l_angles_methods
  USE atomic_kind_types,               ONLY: get_atomic_kind
  USE cp_array_utils,                  ONLY: cp_2d_r_output
  USE cp_block_matrix,                 ONLY: cp_block_matrix_type,&
                                             cp_get_matrix,&
                                             cp_sp_get
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error,&
                                             cp_precondition_failed
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_logger_type,&
                                             cp_to_string,&
                                             cp_warning_level
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE pao_bas_projection_types,        ONLY: cp_proj_get,&
                                             pao_bas_projection_type
  USE pao_glob_angles_methods,         ONLY: cp_validate,&
                                             pao_g_ang_get
  USE pao_plain_rot,                   ONLY: cp_left_compose_p_rot
  USE pao_qs_env_methods,              ONLY: pao_qs_env_get,&
                                             pao_qs_env_get_atomic_ortho
  USE pao_types,                       ONLY: pao_env_type,&
                                             pao_glob_angles_type
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE sparse_matrix_types,             ONLY: add_block_node,&
                                             get_block_node,&
                                             put_block_node
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='pao_l_angles_methods'
! types
  PUBLIC :: pao_local_angles_type

! core procedures
  PUBLIC :: cp_init, cp_dealloc_ref, cp_set, cp_get, &
       cp_did_change

! core functions
  PUBLIC :: cp_validate, cp_valid

! special procedures
  PUBLIC :: cp_set_angle

! special functions
  PUBLIC :: cp_get_angle, &
       cp_get_size_min_bas, cp_get_size_full_bas, cp_get_size_excl_bas, &
       cp_get_p_indexes, cp_get_p_index, cp_get_non_p_index,&
       cp_get_non_p_indexes, cp_get_unitary_t,&
       cp_get_g_index

  PUBLIC :: pao_l_compose_u ! to do remove from here and plain rot tests

! underlying functions
  PUBLIC :: pao_l_ang_init, pao_l_ang_dealloc_ref,&
       pao_l_ang_validate, pao_l_ang_set,&
       pao_l_ang_get, pao_l_ang_did_change,&
       pao_l_ang_get_angle,&
       pao_l_ang_set_angle, pao_l_ang_get_size_min_bas,&
       pao_l_ang_get_size_full_bas, pao_l_ang_get_size_excl_bas,&
       pao_l_ang_get_p_indexes, pao_l_ang_get_p_index, &
       pao_l_ang_get_non_p_indexes, pao_l_ang_get_non_p_index,&
       pao_l_ang_get_unitary_t, pao_l_ang_set_unitary_t, &
       pao_l_ang_get_g_index, pao_guarantee_owned_NUi

! ============== interfaces ===========

! +++ begin of core methods +++

!! initializes the first, argument.
!! Support optional initial values that depend on the 
!! first argument 
  INTERFACE cp_init
     MODULE PROCEDURE pao_l_ang_init
  END INTERFACE

!! Deallocates the memory allocated by the first argument.
!! Supports error=error
  INTERFACE cp_dealloc_ref
     MODULE PROCEDURE pao_l_ang_dealloc_ref
  END INTERFACE

!! returns true if the content of the first argumenst is valid.
!! writes the errors ar warnings, full validation (long)
  INTERFACE cp_validate
     MODULE PROCEDURE pao_l_ang_validate
  END INTERFACE

!! returns true if the first argument is valid. Only a fast minimal 
!! validation is performed (no dangling pointers)
  INTERFACE cp_valid
     MODULE PROCEDURE pao_l_angles_valid
  END INTERFACE

!! sets the value of various optonal attributes of the firs argument.
!! the optional attributes depend on the type of the first
  INTERFACE cp_set
     MODULE PROCEDURE pao_l_ang_set
  END INTERFACE

!! gets values from the first argument via optional keyword arguments
!! that depend on the type of the first argument
  INTERFACE cp_get
     MODULE PROCEDURE pao_l_ang_get
  END INTERFACE

!! tells the object that its internal values have changes and its cached
!! values may be invalid
  INTERFACE cp_did_change
     MODULE PROCEDURE pao_l_ang_did_change
  END INTERFACE

! +++ end of the core methods +++

!! returns the angle in the first arg identified by the second (min_bas)
!! and the third arg (full_bas) 
  INTERFACE cp_get_angle
     MODULE PROCEDURE pao_l_ang_get_angle
  END INTERFACE

!! sets the angle in the first arg identified by the second (min_bas)
!! and the third arg (full_bas) to the value in the fourth (value)
  INTERFACE cp_set_angle
     MODULE PROCEDURE pao_l_ang_set_angle
  END INTERFACE

!! returns the size of the minimal basis
  INTERFACE cp_get_size_min_bas
     MODULE PROCEDURE pao_l_ang_get_size_min_bas
  END INTERFACE

!! returns the size of the full basis
  INTERFACE cp_get_size_full_bas
     MODULE PROCEDURE pao_l_ang_get_size_full_bas
  END INTERFACE

!! returns the size of the excluded basis
  INTERFACE cp_get_size_excl_bas
     MODULE PROCEDURE pao_l_ang_get_size_excl_bas
  END INTERFACE

!! returns the projection indexes
  INTERFACE cp_get_p_indexes
     MODULE PROCEDURE pao_l_ang_get_p_indexes
  END INTERFACE

!! returns the indexes of the excluded basis
  INTERFACE cp_get_non_p_indexes
     MODULE PROCEDURE pao_l_ang_get_non_p_indexes
  END INTERFACE

!! returns an index of the projection
  INTERFACE cp_get_p_index
     MODULE PROCEDURE pao_l_ang_get_p_index
  END INTERFACE

!! returns an index of the excluded basis
  INTERFACE cp_get_non_p_index
     MODULE PROCEDURE pao_l_ang_get_non_p_index
  END INTERFACE

!! gets the unitary transformation
  INTERFACE cp_get_unitary_t
     MODULE PROCEDURE pao_l_ang_get_unitary_t
  END INTERFACE

!! sets the unitary transformation
  INTERFACE cp_set_unitary_t
     MODULE PROCEDURE pao_l_ang_set_unitary_t
  END INTERFACE

!! returns the global indexes
INTERFACE cp_get_g_index
   MODULE PROCEDURE pao_l_ang_get_g_index
END INTERFACE

!!***
!****************************************************************************

!!****s* pao_types/pao_local_angles_type [1.0] *
!!
!!   NAME
!!     pao_local_angles_type
!!
!!   FUNCTION
!!     This structure hold the angles that regard one atom.
!!     it should be used to access the angles when you want an "atom view"
!!     (for example to build the unitary transformation).
!!
!!   NOTES
!!      Do not use the global structure, the idea is that parallelization
!!      (distributing various atoms to various processors) should be easier
!!      using this structure.
!!
!!   ATTRIBUTES
!!     initialized: internal flag to (maybe) find errors that come
!!          from not having initialized the structure. is true if 
!!          this type has been initialized
!!     atom_nr : the atom this object focalizes on
!!     glob_arg: a pointer to the global angles
!!     tmp_m: tmeporary array, used for internal calculations
!!      
!!     readonly (virtual) attributes:
!!
!!     n_min_bas: the size of the minimal basis on this atom.
!!     n_full_bas: the size of the full basis on this atom.
!!     pao_env: the pao environment this local angles lives in.
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  TYPE pao_local_angles_type
     PRIVATE
     LOGICAL :: initialized
     INTEGER :: atom_nr, begins_at, size_min_bas
     REAL(KIND = dp), DIMENSION(:,:), POINTER :: tmp_m
     TYPE(qs_environment_type), POINTER :: qs_env
     TYPE(global_environment_type), POINTER :: global_env
     TYPE(pao_bas_projection_type), POINTER :: projection
     TYPE(pao_glob_angles_type), POINTER :: glob_angles
  END TYPE pao_local_angles_type
!!***
!****************************************************************************

CONTAINS

! ====== core methods ======

!!****f* pao_l_angles_methods/pao_l_ang_init [1.0] *
!!
!!   NAME
!!     pao_l_ang_init
!!
!!   SYNOPSIS
!!     Subroutine pao_l_ang_init(local_angles, qs_env, global_env, atom_nr,&
!!         glob_angles, error)
!!       Type(pao_local_angles_type), Intent (OUT):: local_angles
!!       Type(qs_environment_type), Intent (IN), Target:: qs_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Integer, Intent (IN):: atom_nr
!!       Type(pao_glob_angles_type), Target:: glob_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_l_ang_init
!!
!!   FUNCTION
!!     Initializes the pao angles structure
!!
!!   NOTES
!!     add source arg (i.e. copy)?
!!
!!   INPUTS
!!     local_angles: the pao_local angles to be initialized
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling
!!
!!     For the other arguments see pao_l_ang_set
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE pao_l_ang_init(local_angles, qs_env, global_env, atom_nr,&
       glob_angles, error)
    TYPE(pao_local_angles_type), INTENT(out) :: local_angles
    TYPE(qs_environment_type), INTENT(in), &
      TARGET                                 :: qs_env
    TYPE(global_environment_type), &
      INTENT(in), TARGET                     :: global_env
    INTEGER, INTENT(in)                      :: atom_nr
    TYPE(pao_glob_angles_type), TARGET       :: glob_angles
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'pao_l_angles_methods:pao_l_ang_init'

    LOGICAL                                  :: failure

    failure=.FALSE.

    local_angles%initialized=.TRUE.
    local_angles%atom_nr=-1
    local_angles%global_env => global_env
    local_angles%qs_env => qs_env
    local_angles%glob_angles => glob_angles
    NULLIFY(local_angles%projection, local_angles%tmp_m)
    CALL pao_l_ang_set(local_angles, atom_nr, error=error)
  END SUBROUTINE pao_l_ang_init
!***************************************************************************

!!****f* pao_l_angles_methods/pao_l_ang_dealloc_ref [1.0] *
!!
!!   NAME
!!     pao_l_ang_dealloc_ref
!!
!!   SYNOPSIS
!!     Subroutine pao_l_ang_dealloc_ref(local_angles, error)
!!       Type(pao_local_angles_type), Intent (INOUT):: local_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_l_ang_dealloc_ref
!!
!!   FUNCTION
!!     deallocates the memory that the given local angles structure
!!     has allocated. Does not deallocate the local angles structure itself
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the local angles structure to be deallocated
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
  SUBROUTINE pao_l_ang_dealloc_ref(local_angles, error)
    TYPE(pao_local_angles_type), &
      INTENT(inout)                          :: local_angles
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'pao_l_angles_methods:pao_l_ang_dealloc_ref'

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(local_angles%initialized,cp_failure_level,routineP,error,failure)
    failureIf: IF (.NOT. failure) THEN
       IF (ASSOCIATED(local_angles%tmp_m)) THEN
          DEALLOCATE(local_angles%tmp_m,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
       NULLIFY(local_angles%glob_angles, local_angles%projection, &
            local_angles%qs_env, local_angles%global_env, local_angles%tmp_m)
       local_angles%initialized=.FALSE.
    END IF failureIf
  END SUBROUTINE pao_l_ang_dealloc_ref
!***************************************************************************

!!****f* pao_l_angles_methods/pao_l_ang_validate [1.0] *
!!
!!   NAME
!!     pao_l_ang_validate
!!
!!   SYNOPSIS
!!     Function pao_l_ang_validate(local_angles, error) Result(res)
!!       Logical:: res
!!       Type(pao_local_angles_type), Intent (IN):: local_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function pao_l_ang_validate
!!
!!   FUNCTION
!!     checks if the local angles object is valid (returns .true. if it is)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the local angles object to check
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
  FUNCTION pao_l_ang_validate(local_angles, error) RESULT(res)
    TYPE(pao_local_angles_type), INTENT(in)  :: local_angles
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'pao_l_angles_methods:pao_l_ang_validate'

    INTEGER                                  :: beginAt
    LOGICAL                                  :: failure
    TYPE(pao_env_type), POINTER              :: pao_env

    failure=.FALSE.
    NULLIFY(pao_env)

    CPPrecondition(local_angles%initialized,cp_failure_level,routineP,error,failure)
    CPAssert(ASSOCIATED(local_angles%glob_angles),cp_warning_level,routineP,error,failure)
    CPAssert(ASSOCIATED(local_angles%qs_env),cp_warning_level,routineP,error,failure)
    CPAssert(ASSOCIATED(local_angles%global_env),cp_warning_level,routineP,error,failure)
    IF (cp_debug .AND. .NOT. failure) THEN
       CALL cp_assert(cp_validate(local_angles%glob_angles,&
            qs_env=local_angles%qs_env,global_env=local_angles%global_env,&
            error=error),cp_warning_level,cp_assertion_failed,routineP,&
            "referenced pao_glob_angles is invalid, in file "//&
            CPSourceFileRef,&
            error,failure)
    END IF
    IF (.NOT.failure) THEN
       CPAssert(local_angles%atom_nr>0,cp_warning_level,routineP,error,failure)
!       CPAssert(local_angles%atom_nr<=cp_get_n_atom(pao_env),cp_warning_level,routineP,error,failure)
       IF (.NOT.failure) THEN
          CPAssert(local_angles%begins_at>=0,cp_warning_level,routineP,error,failure)
          CALL cp_assert(local_angles%begins_at < &
               SIZE(local_angles%glob_angles%angles),&
               cp_warning_level,cp_assertion_failed,routineP,&
               "begins_at out of range in "//&
               CPSourceFileRef,&
               error,failure)
          CALL get_qs_env(local_angles%qs_env, pao_env=pao_env)
          beginAt=pao_env%angles_begin_at(local_angles%atom_nr)
          CPAssert(local_angles%begins_at==beginAt,cp_warning_level,routineP,error,failure)
       END IF
    END IF
    res=.NOT.failure
  END FUNCTION pao_l_ang_validate
!***************************************************************************


!!****f* pao_l_angles_methods/pao_l_ang_set [1.0] *
!!
!!   NAME
!!     pao_l_ang_set
!!
!!   SYNOPSIS
!!     Subroutine pao_l_ang_set(local_angles, atom_nr, error)
!!       Type(pao_local_angles_type), Intent (INOUT):: local_angles
!!       Integer, Intent (IN):: atom_nr
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_l_ang_set
!!
!!   FUNCTION
!!     sets the values of the various attributes of the local angles
!!     structure.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the object to change
!!     atom_nr: the numbre of the atom this object focalizes on
!!     glob_angles: the global angles object that this object references
!!          SHARED.
!!     source: a local angles object that is copied into this(other arguments
!!             have the precedence for setting the values, the values
!!             of source are taken as default, if the value is not directly
!!             set)
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
  SUBROUTINE pao_l_ang_set(local_angles, atom_nr,error)
    TYPE(pao_local_angles_type), &
      INTENT(inout)                          :: local_angles
    INTEGER, INTENT(in)                      :: atom_nr
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'pao_l_angles_methods:pao_l_ang_set'

    LOGICAL                                  :: failure
    TYPE(cp_error_type)                      :: iError
    TYPE(pao_env_type), POINTER              :: pao_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    failure=.FALSE.
    NULLIFY(pao_env)
    CALL cp_error_init(iError, template_error=error)

    CPPrecondition(local_angles%initialized,cp_failure_level,routineP,error,failure)
    failureIf: IF (.NOT. failure) THEN
       local_angles%atom_nr=atom_nr
       CPAssert(atom_nr>0,cp_warning_level,routineP,error,failure)
       IF (local_angles%atom_nr>0) THEN
          CALL cp_assert(ASSOCIATED(local_angles%glob_angles),&
               cp_failure_level,cp_precondition_failed,routineP,&
               "global_angles not associated in "//&
               CPSourceFileRef,&
               error,failure)
          CALL get_qs_env(local_angles%qs_env,pao_env=pao_env,&
               particle_set=particle_set)
          CPPrecondition(ASSOCIATED(pao_env),cp_failure_level,routineP,error,failure)
          CPPrecondition(ASSOCIATED(particle_set),cp_failure_level,routineP,error,failure)
          IF (.NOT.failure) THEN
             local_angles%begins_at= pao_env%angles_begin_at(atom_nr)
             CALL get_atomic_kind(particle_set(atom_nr)%atomic_kind,&
                  pao_bas_proj=local_angles%projection)
             CPPostcondition(ASSOCIATED(local_angles%projection),cp_failure_level,routineP,error,failure)
          END IF
          IF (.NOT.failure) THEN
             CALL cp_proj_get(local_angles%projection,&
                  size_min_bas=local_angles%size_min_bas, error=error)
          END IF
       ELSE
          NULLIFY(local_angles%projection)
       END IF
    END IF failureIf
    IF (failure) THEN
       NULLIFY(local_angles%projection)
    END IF
    CALL cp_error_dealloc_ref(iError,error=error)
  END SUBROUTINE pao_l_ang_set
!***************************************************************************

!!****f* pao_l_angles_methods/pao_l_ang_get [1.0] *
!!
!!   NAME
!!     pao_l_ang_get
!!
!!   SYNOPSIS
!!     Subroutine pao_l_ang_get(local_angles, atom_nr, pao_env, unitary_t,&
!!         cache_valid, size_min_bas, size_full_bas, size_excl_bas,&
!!         proj_indexes, non_proj_indexes, projection, angles,&
!!         nui_injection, atomic_ortho, error)
!!       Type(pao_local_angles_type), Intent (INOUT):: local_angles
!!       Integer, Optional, Intent (OUT):: atom_nr
!!       Integer, Optional, Intent (OUT):: size_excl_bas, size_full_bas,&
!!         size_min_bas
!!       Integer, Optional, Dimension(:), Pointer:: proj_indexes,&
!!         non_proj_indexes
!!       Type(pao_env_type), Optional, Pointer:: pao_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!       Real(KIND = dp), Dimension(:,:), Pointer, Optional:: unitary_t,&
!!         atomic_ortho, nui_injection
!!       Logical, Intent (OUT), Optional:: cache_valid
!!       Type(pao_bas_projection_type), Pointer, Optional:: projection
!!       Real(KIND = dp), Dimension(:), Pointer, Optional:: angles
!!     End Subroutine pao_l_ang_get
!!
!!   FUNCTION
!!     gets various attributes of the pao local angles structure
!!
!!   NOTES
!!     I try to keep the glob_angles private.
!!     I have not added al local angles matrix
!!
!!   INPUTS
!!     local_angles: the local angles to get the info from.
!!     atom_nr: the atom this local angles is focalized on.
!!     n_min_bas: the size of the minimal basis on this atom.
!!     n_full_bas: the size of the full basis on this atom.
!!     n_excl_bas: the size of the excluded basis on this atom
!!            (n_full_bas-n_min_bas).
!!     pao_env: the pao environment this local angles lives in.
!!     angles: the local angles in a min_bas first order
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
  SUBROUTINE pao_l_ang_get(local_angles, atom_nr, &
       pao_env, unitary_t,cache_valid,&
       size_min_bas, size_full_bas, &
       size_excl_bas, proj_indexes, non_proj_indexes,&
       projection, angles, NUi_injection, atomic_ortho, error)
    TYPE(pao_local_angles_type), &
      INTENT(inout)                          :: local_angles
    INTEGER, INTENT(out), OPTIONAL           :: atom_nr
    TYPE(pao_env_type), OPTIONAL, POINTER    :: pao_env
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: unitary_t
    LOGICAL, INTENT(out), OPTIONAL           :: cache_valid
    INTEGER, INTENT(out), OPTIONAL           :: size_min_bas, size_full_bas, &
                                                size_excl_bas
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: proj_indexes, non_proj_indexes
    TYPE(pao_bas_projection_type), &
      OPTIONAL, POINTER                      :: projection
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: angles
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: NUi_injection, atomic_ortho
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'pao_l_angles_methods:pao_l_ang_get'

    LOGICAL                                  :: failure
    TYPE(pao_env_type), POINTER              :: my_pao_env

    failure=.FALSE.
    NULLIFY(my_pao_env)

    CPAssert(local_angles%initialized,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(atom_nr)) THEN
          atom_nr=local_angles%atom_nr
       END IF
       IF (PRESENT(pao_env)) THEN
          CALL get_qs_env(local_angles%qs_env, pao_env=pao_env)
       END IF
       IF (PRESENT(projection)) projection => local_angles%projection
       CPPrecondition(ASSOCIATED(local_angles%projection),cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT.failure) THEN
       CALL cp_proj_get(local_angles%projection, size_min_bas=size_min_bas,&
            size_full_bas=size_full_bas, size_excl_bas=size_excl_bas,&
            proj_indexes=proj_indexes, non_proj_indexes=non_proj_indexes,&
            error=error)
       IF (PRESENT(unitary_t)) THEN
          unitary_t => pao_l_ang_get_unitary_t(local_angles,error=error)
       END IF
       IF (PRESENT(cache_valid)) THEN
          IF (ASSOCIATED(local_angles%glob_angles%ut_calculated).AND.&
               local_angles%atom_nr>0) THEN
             cache_valid=local_angles%glob_angles%ut_calculated(&
                  local_angles%atom_nr)
          ELSE
             cache_valid=.FALSE.
          END IF
       END IF
       IF (PRESENT(angles)) THEN
          NULLIFY(angles)
          ! seems to give problems on some platforms
          CPAssert(.FALSE.,cp_warning_level,routineP,error,failure)
          CALL get_qs_env(local_angles%qs_env, pao_env=my_pao_env)
          angles => local_angles%glob_angles%angles(local_angles%begins_at:&
               (my_pao_env%angles_begin_at(local_angles%atom_nr+1)-1))
       END IF
       IF (PRESENT(NUi_injection)) THEN
          NUi_injection => pao_l_ang_get_NUi_injection(local_angles,error=error)
       END IF
       IF (PRESENT(atomic_ortho)) THEN
          atomic_ortho => pao_l_ang_get_atomic_ortho(local_angles,error=error)
       END IF
    END IF
  END SUBROUTINE pao_l_ang_get
!***************************************************************************

!!****f* pao_l_angles_methods/pao_l_angles_valid [1.0] *
!!
!!   NAME
!!     pao_l_angles_valid
!!
!!   FUNCTION
!!     performs minimal validation on the object
!!
!!   NOTES
!!     do not check pao_env?
!!
!!   INPUTS
!!     local_angles: the object to validate
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
  FUNCTION pao_l_angles_valid(local_angles, error) RESULT(res)
    TYPE(pao_local_angles_type), INTENT(in)  :: local_angles
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'pao_l_angles_methods:pao_l_angles_valid'

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(local_angles%initialized,cp_failure_level,routineP,error,failure)
    CPAssert(ASSOCIATED(local_angles%glob_angles),cp_warning_level,routineP,error,failure)
    res=.NOT.failure
  END FUNCTION pao_l_angles_valid
!***************************************************************************

! ====== special methods =====

!!****f* pao_l_angles_methods/pao_l_ang_get_p_indexes [1.0] *
!!
!!   NAME
!!     pao_l_ang_get_p_indexes
!!
!!   SYNOPSIS
!!     Function pao_l_ang_get_p_indexes(local_angles, error) Result(res)
!!       Integer, Dimension(:), Pointer:: res
!!       Type(pao_local_angles_type), Intent (IN):: local_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function pao_l_ang_get_p_indexes
!!
!!   FUNCTION
!!     returns the indexes of the projection
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the local angles you get values from
!!     i: number of the index you want to get
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     1.2002 created [fawzi]
!!     4.2002 use bas_projection [fawzi]
!!
!!*** **********************************************************************
  FUNCTION pao_l_ang_get_p_indexes(local_angles,error) RESULT(res)
    TYPE(pao_local_angles_type), INTENT(in)  :: local_angles
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    INTEGER, DIMENSION(:), POINTER           :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_l_ang_get_p_indexes', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    IF (cp_debug .OR. debug_this_module) THEN
       CPPrecondition(local_angles%initialized,cp_failure_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(local_angles%projection),cp_failure_level,routineP,error,failure)
    END IF

    IF (.NOT. failure) THEN
       CALL cp_proj_get(local_angles%projection, proj_indexes=res,error=error)
    ELSE
       NULLIFY(res)
    END IF
  END FUNCTION pao_l_ang_get_p_indexes
!***************************************************************************

!!****f* pao_l_angles_methods/pao_l_ang_get_non_p_indexes [1.0] *
!!
!!   NAME
!!     pao_l_ang_get_non_p_indexes
!!
!!   SYNOPSIS
!!     Function pao_l_ang_get_non_p_indexes(local_angles, error)&
!!         Result(res)
!!       Integer, Dimension(:), Pointer:: res
!!       Type(pao_local_angles_type), Intent (IN):: local_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function pao_l_ang_get_non_p_indexes
!!
!!   FUNCTION
!!     returns the indexes of the excluded basis
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the local angles you get values from
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
  FUNCTION pao_l_ang_get_non_p_indexes(local_angles,error) RESULT(res)
    TYPE(pao_local_angles_type), INTENT(in)  :: local_angles
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    INTEGER, DIMENSION(:), POINTER           :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_l_ang_get_p_indexes', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    IF (cp_debug .OR. debug_this_module) THEN
       CPPrecondition(local_angles%initialized,cp_failure_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(local_angles%projection),cp_failure_level,routineP,error,failure)
    END IF

    IF (.NOT. failure) THEN
       CALL cp_proj_get(local_angles%projection, non_proj_indexes=res,error=error)
    ELSE
       NULLIFY(res)
    END IF
  END FUNCTION pao_l_ang_get_non_p_indexes
!***************************************************************************

!!****f* pao_l_angles_methods/pao_l_ang_get_p_index [1.0] *
!!
!!   NAME
!!     pao_l_ang_get_p_index
!!
!!   SYNOPSIS
!!     Function pao_l_ang_get_p_index(local_angles, min_bas, error)&
!!         Result(res)
!!       Integer:: res
!!       Type(pao_local_angles_type), Intent (IN):: local_angles
!!       Integer, Intent (IN):: min_bas
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function pao_l_ang_get_p_index
!!
!!   FUNCTION
!!     returns the index min_bas of the projection
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the local angles you get values from
!!     min_bas: the index of the minimal basis (1:size_min_bas)
!!     i: number of the index you want to get
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     1.2002 created [fawzi]
!!     4.2002 use bas_projection [fawzi]
!!
!!*** **********************************************************************
  FUNCTION pao_l_ang_get_p_index(local_angles, min_bas, error) RESULT(res)
    TYPE(pao_local_angles_type), INTENT(in)  :: local_angles
    INTEGER, INTENT(in)                      :: min_bas
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    INTEGER                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_l_ang_get_p_index', &
      routineP = moduleN//':'//routineN

    INTEGER, DIMENSION(:), POINTER           :: p
    LOGICAL                                  :: failure

    failure=.FALSE.
    IF (cp_debug .OR. debug_this_module) THEN
       CPPrecondition(local_angles%initialized,cp_failure_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(local_angles%projection),cp_failure_level,routineP,error,failure)
    END IF

    IF (.NOT. failure) THEN
       CALL cp_proj_get(local_angles%projection, proj_indexes=p,error=error)
       res=p(min_bas)
    ELSE
       res=-1
    END IF
  END FUNCTION pao_l_ang_get_p_index
!***************************************************************************

!!****f* pao_l_angles_methods/pao_l_ang_get_non_p_index [1.0] *
!!
!!   NAME
!!     pao_l_ang_get_non_p_index
!!
!!   SYNOPSIS
!!     Function pao_l_ang_get_non_p_index(local_angles, excl_bas, error)&
!!         Result(res)
!!       Integer:: res
!!       Integer, Intent (IN):: excl_bas
!!       Type(pao_local_angles_type), Intent (IN):: local_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function pao_l_ang_get_non_p_index
!!
!!   FUNCTION
!!     returns the index excl_bas of the excluded basis
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - local_angles: the local angles you get values from
!!     - excl_bas: the local index on the excluded basis (1:size_excl_bas)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     1.2002 created [fawzi]
!!     4.2002 use bas projection [fawzi]
!!
!!*** **********************************************************************
  FUNCTION pao_l_ang_get_non_p_index(local_angles,excl_bas,error) RESULT(res)
    TYPE(pao_local_angles_type), INTENT(in)  :: local_angles
    INTEGER, INTENT(in)                      :: excl_bas
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    INTEGER                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_l_ang_get_p_index', &
      routineP = moduleN//':'//routineN

    INTEGER, DIMENSION(:), POINTER           :: np
    LOGICAL                                  :: failure

    failure=.FALSE.
    IF (cp_debug .OR. debug_this_module) THEN
       CPPrecondition(local_angles%initialized,cp_failure_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(local_angles%projection),cp_failure_level,routineP,error,failure)
    END IF

    IF (.NOT. failure) THEN
       CALL cp_proj_get(local_angles%projection, non_proj_indexes=np,error=error)
       res=np(excl_bas)
    ELSE
       res=-1
    END IF
  END FUNCTION pao_l_ang_get_non_p_index
!***************************************************************************

!!****f* pao_l_angles_methods/pao_l_ang_get_angle [1.0] *
!!
!!   NAME
!!     pao_l_ang_get_angle
!!
!!   SYNOPSIS
!!     Function pao_l_ang_get_angle(local_angles, min_bas, excl_bas, error)&
!!         Result(res)
!!       Real(KIND = dp):: res
!!       Type(pao_local_angles_type), Intent (IN):: local_angles
!!       Integer, Intent (IN):: min_bas, excl_bas
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function pao_l_ang_get_angle
!!
!!   FUNCTION
!!     returns the requested angle
!!
!!   NOTES
!!     should be fast, and hopefully inlined
!!
!!   INPUTS
!!     local_angles: the local angles object you want to get the angles from
!!     min_bas: the local index of the minimal basis of which you want 
!!            the angle
!!     excl_bas: the local index of the excluded basis of which you want the
!!            angle
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
  FUNCTION pao_l_ang_get_angle(local_angles, min_bas, excl_bas, &
       error) RESULT(res)
    TYPE(pao_local_angles_type), INTENT(in)  :: local_angles
    INTEGER, INTENT(in)                      :: min_bas, excl_bas
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    REAL(KIND=dp)                            :: res

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'pao_l_angles_methods:pao_l_ang_get_angle'

    INTEGER                                  :: size_excl_bas
    LOGICAL                                  :: failure
    TYPE(pao_glob_angles_type), POINTER      :: glob_angles

    failure=.FALSE.
    glob_angles => local_angles%glob_angles

    debugIf: IF (debug_this_module.AND.cp_debug) THEN
       CPPrecondition(local_angles%initialized,cp_failure_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(glob_angles),cp_failure_level,routineP,error,failure)
       IF (.NOT.failure) THEN
          CPPrecondition(ASSOCIATED(local_angles%projection),cp_failure_level,routineP,error,failure)
       END IF
       IF (.NOT.failure) THEN
          CALL cp_proj_get(local_angles%projection,size_excl_bas=size_excl_bas,&
               error=error)
          CALL cp_assert(min_bas>0.AND.min_bas <= &
               local_angles%size_min_bas,&
               cp_failure_level,cp_assertion_failed,routineP,&
               "min_bas out of range in "//&
               CPSourceFileRef,&
               error,failure)
          CALL cp_assert(excl_bas>0.AND.excl_bas <= &
               size_excl_bas,&
               cp_failure_level,cp_assertion_failed,routineP,&
               "excl_bas out of range in "//&
               CPSourceFileRef,&
               error,failure)
       END IF
    END IF debugIf
    failureIf: IF (.NOT. failure) THEN
       res=glob_angles%angles(local_angles%begins_at + &
            local_angles%size_min_bas*(excl_bas-1) + &
            min_bas)
    ELSE failureIf
       res=-1
    END IF failureIf
  END FUNCTION pao_l_ang_get_angle
!***************************************************************************

!!****f* pao_l_angles_methods/pao_l_ang_set_angle [1.0] *
!!
!!   NAME
!!     pao_l_ang_set_angle
!!
!!   SYNOPSIS
!!     Subroutine pao_l_ang_set_angle(local_angles, min_bas, excl_bas,&
!!         value, error)
!!       Type(pao_local_angles_type), Intent (IN):: local_angles
!!       Integer, Intent (IN):: min_bas, excl_bas
!!       Real(KIND = dp), Intent (IN):: value
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_l_ang_set_angle
!!
!!   FUNCTION
!!     returns the requested angle
!!
!!   NOTES
!!     should be fast, and hopefully inlined
!!
!!   INPUTS
!!     local_angles: the local angles object you want to set the angles from
!!     min_bas: the local index of the minimal basis of which you want 
!!            the angle
!!     full_bas: the local index of the full basis of which you want the
!!            angle
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
  SUBROUTINE pao_l_ang_set_angle(local_angles, min_bas, excl_bas, &
       value,error)
    TYPE(pao_local_angles_type), INTENT(in)  :: local_angles
    INTEGER, INTENT(in)                      :: min_bas, excl_bas
    REAL(KIND=dp), INTENT(in)                :: value
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_l_ang_set_angle', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: size_excl_bas
    LOGICAL                                  :: failure
    TYPE(pao_glob_angles_type), POINTER      :: glob_angles

!actually the substructure globalangles is modifed

    failure=.FALSE.
    glob_angles => local_angles%glob_angles

    debugIf: IF (debug_this_module.AND.cp_debug) THEN
       CPPrecondition(local_angles%initialized,cp_failure_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(glob_angles),cp_failure_level,routineP,error,failure)
       IF (.NOT.failure) THEN
          CALL cp_assert(min_bas>0.AND.min_bas <= &
               local_angles%size_min_bas,&
               cp_failure_level,cp_assertion_failed,routineP,&
               "min_bas out of range in "//&
               CPSourceFileRef,&
               error,failure)
          CALL cp_proj_get(local_angles%projection,size_excl_bas=size_excl_bas,&
               error=error)
          CALL cp_assert(excl_bas>0.AND.excl_bas <= &
               size_excl_bas,&
               cp_failure_level,cp_assertion_failed,routineP,&
               "excl_bas out of range in "//&
               CPSourceFileRef,&
               error,failure)
       END IF
    END IF debugIf
    failureIf: IF (.NOT. failure) THEN
       glob_angles%angles(local_angles%begins_at + &
            local_angles%size_min_bas*(excl_bas-1) + &
            min_bas) = value
    END IF failureIf
    NULLIFY(glob_angles)
  END SUBROUTINE pao_l_ang_set_angle
!***************************************************************************

!!****f* pao_l_angles_methods/pao_l_ang_get_size_excl_bas [1.0] *
!!
!!   NAME
!!     pao_l_ang_get_size_excl_bas
!!
!!   SYNOPSIS
!!     Function pao_l_ang_get_size_excl_bas(local_angles, error)&
!!         Result(res)
!!       Integer:: res
!!       Type(pao_local_angles_type), Intent (IN):: local_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function pao_l_ang_get_size_excl_bas
!!
!!   FUNCTION
!!     returns the size of the excluded basis
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the local angles object you get the information from
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
  FUNCTION pao_l_ang_get_size_excl_bas(local_angles, error) RESULT(res)
    TYPE(pao_local_angles_type), INTENT(in)  :: local_angles
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    INTEGER                                  :: res

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'pao_l_angles_methods:pao_l_ang_get_size_excl_bas'

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(local_angles%initialized, cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(local_angles%projection),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CALL cp_proj_get(local_angles%projection,size_excl_bas=res,error=error)
    ELSE
       res=-1
    END IF
  END FUNCTION pao_l_ang_get_size_excl_bas
!***************************************************************************

!!****f* pao_l_angles_methods/pao_l_ang_get_size_min_bas [1.0] *
!!
!!   NAME
!!     pao_l_ang_get_size_min_bas
!!
!!   SYNOPSIS
!!     Function pao_l_ang_get_size_min_bas(local_angles, error) Result(res)
!!       Integer:: res
!!       Type(pao_local_angles_type), Intent (IN):: local_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function pao_l_ang_get_size_min_bas
!!
!!   FUNCTION
!!     returns the size of the minimal basis
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the local angles object you get the information from
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
  FUNCTION pao_l_ang_get_size_min_bas(local_angles, error) RESULT(res)
    TYPE(pao_local_angles_type), INTENT(in)  :: local_angles
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    INTEGER                                  :: res

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'pao_l_angles_methods:pao_l_ang_get_size_min_bas'

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(local_angles%initialized, cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       res=local_angles%size_min_bas
    ELSE
       res=-1
    END IF
  END FUNCTION pao_l_ang_get_size_min_bas
!***************************************************************************

!!****f* pao_l_angles_methods/pao_l_ang_get_size_full_bas [1.0] *
!!
!!   NAME
!!     pao_l_ang_get_size_full_bas
!!
!!   SYNOPSIS
!!     Function pao_l_ang_get_size_full_bas(local_angles, error)&
!!         Result(res)
!!       Integer:: res
!!       Type(pao_local_angles_type), Intent (IN):: local_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function pao_l_ang_get_size_full_bas
!!
!!   FUNCTION
!!     returns the size of the full basis
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the local angles object you get the information from
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
  FUNCTION pao_l_ang_get_size_full_bas(local_angles, error) RESULT(res)
    TYPE(pao_local_angles_type), INTENT(in)  :: local_angles
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    INTEGER                                  :: res

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'pao_l_angles_methods:pao_l_ang_get_size_full_bas'

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(local_angles%initialized, cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(local_angles%projection),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CALL cp_proj_get(local_angles%projection,size_full_bas=res,error=error)
    ELSE
       res=-1
    END IF
  END FUNCTION pao_l_ang_get_size_full_bas
!***************************************************************************

!!****f* pao_l_angles_methods/pao_l_ang_did_change [1.0] *
!!
!!   NAME
!!     pao_l_ang_did_change
!!
!!   SYNOPSIS
!!     Subroutine pao_l_ang_did_change(local_angles, error)
!!       Type(pao_local_angles_type), Intent (INOUT):: local_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_l_ang_did_change
!!
!!   FUNCTION
!!     invalidate the cached unitary transformation of this atom.
!!     to be called when the angles change
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the local angles object you get the info from
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
  SUBROUTINE pao_l_ang_did_change(local_angles,error)
! changes in local_angles%glob_angles
    TYPE(pao_local_angles_type), &
      INTENT(inout)                          :: local_angles
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'pao_l_angles_methods:pao_l_ang_did_change'

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(local_angles%initialized,cp_failure_level,routineP,error,failure)
    CPPrecondition(local_angles%atom_nr>0,cp_warning_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (ASSOCIATED(local_angles%glob_angles%ut_calculated)) THEN
          local_angles%glob_angles%ut_calculated(local_angles%atom_nr)=.FALSE.
       END IF
       IF (ASSOCIATED(local_angles%glob_angles%NUi_injection_calculated)) THEN
          local_angles%glob_angles%NUi_injection_calculated&
               (local_angles%atom_nr)=.FALSE.
       END IF
    END IF
  END SUBROUTINE pao_l_ang_did_change
!***************************************************************************

!!****f* pao_l_angles_methods/pao_l_ang_set_unitary_t [1.0] *
!!
!!   NAME
!!     pao_l_ang_set_unitary_t
!!
!!   SYNOPSIS
!!     Subroutine pao_l_ang_set_unitary_t(local_angles, unitary_t, error)
!!       Type(pao_local_angles_type), Intent (INOUT):: local_angles
!!       Real(KIND = dp), Dimension(:,:), Intent (IN), Target:: unitary_t
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_l_ang_set_unitary_t
!!
!!   FUNCTION
!!     changes the cached unitary transformation
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the local angles you want to change
!!     unitary_t: the new unitary trasformation
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
  SUBROUTINE pao_l_ang_set_unitary_t(local_angles,unitary_t,error)
! changes in glob_angles
    TYPE(pao_local_angles_type), &
      INTENT(inout)                          :: local_angles
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(in), TARGET                     :: unitary_t
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'pao_l_angles_methods:pao_l_ang_set_unitary_t'

    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: unitary_t_ptr
    TYPE(cp_block_matrix_type), POINTER      :: full_u_t

    failure=.FALSE.

    CPPrecondition(local_angles%initialized,cp_failure_level,routineP,error,failure)
    CPPrecondition(local_angles%atom_nr>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL pao_g_ang_get(local_angles%glob_angles,qs_env=local_angles%qs_env,&
            global_env=local_angles%global_env, unitary_t=full_u_t,&
            error=error)
       CPPrecondition(ASSOCIATED(full_u_t),cp_warning_level,routineP,error,failure)
       CALL cp_assert(ASSOCIATED(local_angles%glob_angles%ut_calculated),&
            cp_failure_level,cp_assertion_failed,routineP,&
            "ut_calculated not allocated in "//&
            CPSourceFileRef,&
            error,failure)
       IF (.NOT.failure) THEN
          local_angles%glob_angles%ut_calculated(local_angles%atom_nr)=.TRUE.
          unitary_t_ptr => unitary_t
          CALL put_block_node(cp_get_matrix(full_u_t),&
               block_row=local_angles%atom_nr,&
               block_col=local_angles%atom_nr, BLOCK=unitary_t_ptr)
       END IF
    END IF
  END SUBROUTINE pao_l_ang_set_unitary_t
!***************************************************************************

!!****f* pao_plain_rot/pao_l_compose_u [1.0] *
!!
!!   NAME
!!     pao_l_compose_u
!!
!!   SYNOPSIS
!!     Subroutine pao_l_compose_u(matrix, local_angles, error)
!!       Real(KIND = dp), Dimension(:,:), Intent (INOUT):: matrix
!!       Type(pao_local_angles_type), Intent (INOUT):: local_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_l_compose_u
!!
!!   FUNCTION
!!     calculates U*matrix where U is the unitary transformation defined by
!!     the given angles. the sequence is: smallbasis then full basis in 
!!     (local) growing order: u(1,1)*u(2,1)*...*u(size_small_bas,1)*...*
!!      u(size_small_bas,size_full_bas)*matrix
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     matrix: the matrix that will be modified
!!     local_angles: the angles that define U
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
  SUBROUTINE pao_l_compose_u(matrix,local_angles,error)
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(inout)                          :: matrix
    TYPE(pao_local_angles_type), &
      INTENT(inout)                          :: local_angles
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineP = 'pao_plain_rot:pao_l_compose_u'

    INTEGER                                  :: i_excl_bas, i_min_bas, &
                                                n_full_bas
    INTEGER, DIMENSION(:), POINTER           :: non_proj_indexes, proj_indexes
    LOGICAL                                  :: failure

    failure=.FALSE.

    CALL pao_l_ang_get(local_angles,&
         size_full_bas=n_full_bas, proj_indexes=proj_indexes,&
         non_proj_indexes=non_proj_indexes)
    CPPrecondition(SIZE(proj_indexes)>0,cp_warning_level,routineP,error,failure)
    CPPrecondition(SIZE(non_proj_indexes)>=0,cp_warning_level,routineP,error,failure)
    CPPrecondition(SIZE(matrix,1)==n_full_bas,cp_warning_level,routineP,error,failure)
    CPPrecondition(SIZE(matrix,2)==n_full_bas,cp_warning_level,routineP,error,failure)
    failureIf: IF (.NOT. failure) THEN
       DO i_excl_bas=SIZE(non_proj_indexes),1,-1
          DO i_min_bas=SIZE(proj_indexes),1,-1
             CALL cp_left_compose_p_rot(matrix=matrix,&
                  angle=cp_get_angle(local_angles,min_bas=i_min_bas,&
                  excl_bas=i_excl_bas,error=error),&
                  i=proj_indexes(i_min_bas),&
                  j=non_proj_indexes(i_excl_bas),&
                  error=error)
          END DO
       END DO
    END IF failureIf
  END SUBROUTINE pao_l_compose_u
!***************************************************************************

!!****f* pao_l_angles_methods/pao_l_ang_get_unitary_t [1.0] *
!!
!!   NAME
!!     pao_l_ang_get_unitary_t
!!
!!   SYNOPSIS
!!     Function pao_l_ang_get_unitary_t(local_angles, error) Result(res)
!!       Type(pao_local_angles_type), Intent (INOUT):: local_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!       Real(KIND = dp), Dimension(:,:), Pointer:: res
!!     End Function pao_l_ang_get_unitary_t
!!
!!   FUNCTION
!!     returns the unitary trasformation defined by these angles
!!
!!   NOTES
!!     the trasformation is cached
!!
!!   INPUTS
!!     local_angles: the local angles you get the information from
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
  FUNCTION pao_l_ang_get_unitary_t(local_angles,error) RESULT(res)
    TYPE(pao_local_angles_type), &
      INTENT(inout)                          :: local_angles
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: res

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'pao_l_angles_methods:pao_l_ang_get_unitary_t'

    INTEGER                                  :: i
    LOGICAL                                  :: failure
    TYPE(cp_block_matrix_type), POINTER      :: unitary_t

    failure=.FALSE.
    NULLIFY(res)

    IF (cp_debug .AND. debug_this_module) THEN
       CPPrecondition(cp_valid(local_angles),cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT. failure) THEN
       CALL pao_g_ang_get(local_angles%glob_angles,qs_env=local_angles%qs_env,&
            global_env=local_angles%global_env, unitary_t=unitary_t,&
            error=error)
       CPPrecondition(ASSOCIATED(unitary_t),cp_failure_level,routineP,error,failure)
       CALL cp_assert(ASSOCIATED(local_angles%glob_angles%ut_calculated),&
            cp_failure_level,cp_assertion_failed,routineP,&
            "ut_calculated not associated in "//&
            CPSourceFileRef,&
            error,failure)
       IF (.NOT.failure) THEN
          IF (local_angles%glob_angles%ut_calculated(local_angles%atom_nr)) THEN
             CALL get_block_node(cp_get_matrix(unitary_t),&
                  BLOCK=res,block_row=local_angles%atom_nr,&
                  block_col=local_angles%atom_nr)
             CPPostcondition(ASSOCIATED(res),cp_failure_level,routineP,error,failure)
          ELSE
             CALL get_block_node(cp_get_matrix(unitary_t),&
                  BLOCK=res,block_row=local_angles%atom_nr,&
                  block_col=local_angles%atom_nr)
             IF (.NOT.ASSOCIATED(res)) THEN
                CALL add_block_node(cp_get_matrix(unitary_t),&
                     local_angles%atom_nr,local_angles%atom_nr,BLOCK=res)
             END IF
             CPPostcondition(ASSOCIATED(res),cp_failure_level,routineP,error,failure)
             IF (.NOT.failure) THEN
                CPPrecondition(SIZE(res,1)==SIZE(res,2),cp_failure_level,routineP,error,failure)
             END IF
             IF (.NOT.failure) THEN
                res=0.0_dp
                FORALL (i=1:SIZE(res,1)) res(i,i)=1.0_dp
                CALL pao_l_compose_u(res,local_angles,error=error)
                local_angles%glob_angles%ut_calculated(local_angles%atom_nr)=.TRUE.
             END IF
          END IF
       END IF
    END IF
    IF (failure) NULLIFY(res)
  END FUNCTION pao_l_ang_get_unitary_t
!***************************************************************************

!!****f* pao_l_angles_methods/pao_l_ang_get_g_index [1.0] *
!!
!!   NAME
!!     pao_l_ang_get_g_index
!!
!!   SYNOPSIS
!!     Function pao_l_ang_get_g_index(local_angles, min_bas, excl_bas,&
!!         error) Result(res)
!!       Integer:: res
!!       Type(pao_local_angles_type), Intent (INOUT):: local_angles
!!       Integer, Intent (IN):: min_bas, excl_bas
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function pao_l_ang_get_g_index
!!
!!   FUNCTION
!!     returns the global angles index for the angle described with
!!     the given local indexes. 
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the actual local angles object
!!     min_bas: the local minimal basis index
!!     excl_bas: the local excluded basis index
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
FUNCTION pao_l_ang_get_g_index(local_angles, min_bas, excl_bas, error)&
     RESULT(res)
    TYPE(pao_local_angles_type), &
      INTENT(inout)                          :: local_angles
    INTEGER, INTENT(in)                      :: min_bas, excl_bas
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    INTEGER                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_l_ang_get_g_index', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (cp_debug .AND. debug_this_module) THEN
     CPPrecondition(local_angles%initialized,cp_failure_level,routineP,error,failure)
     CPPrecondition(local_angles%atom_nr>0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT. failure) THEN
     res = local_angles%begins_at+local_angles%size_min_bas*(excl_bas-1)&
          + min_bas
  ELSE
     res=-1
  END IF
END FUNCTION pao_l_ang_get_g_index
!***************************************************************************

!!****f* pao_l_angles_methods/pao_l_ang_get_atomic_ortho [1.0] *
!!
!!   NAME
!!     pao_l_ang_get_atomic_ortho
!!
!!   FUNCTION
!!     returns the atomic orthogonalization matrix, calculating it
!!     if necessary
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - local_angles: the local angles object that defines the actual atom
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
FUNCTION pao_l_ang_get_atomic_ortho(local_angles,error) RESULT(res)
    TYPE(pao_local_angles_type), INTENT(in)  :: local_angles
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_l_ang_get_atomic_ortho', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (.NOT.cp_valid(local_angles,error=error)) THEN
     CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT. failure) THEN
     res => pao_qs_env_get_atomic_ortho(local_angles%qs_env,&
          global_env=local_angles%global_env,&
          atom_nr=local_angles%atom_nr,&
          error=error)
  ELSE
     NULLIFY(res)
  END IF
END FUNCTION pao_l_ang_get_atomic_ortho
!***************************************************************************

!!****f* pao_l_angles_methods/pao_l_ang_get_NUi_injection [1.0] *
!!
!!   NAME
!!     pao_l_ang_get_NUi_injection
!!
!!   FUNCTION
!!     returns the injection from the small basis to the full basis,
!!     calculating it if necessary (N ortho, U unitary, i basis injection)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - local_angles: the actual local angles object that defines
!!       the actual atom 
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
FUNCTION pao_l_ang_get_NUi_injection(local_angles, error) RESULT(res)
    TYPE(pao_local_angles_type), &
      INTENT(inout)                          :: local_angles
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_l_ang_get_NUi_injection', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, max_full_bas, max_min_bas, &
                                                stat
    INTEGER, DIMENSION(:), POINTER           :: proj_indexes
    INTEGER, SAVE                            :: iter = 0
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: n_m, u_m
    TYPE(cp_block_matrix_type), POINTER      :: NUi_injection
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(pao_env_type), POINTER              :: pao_env

  failure=.FALSE.
  NULLIFY(res, NUi_injection, n_m, u_m, pao_env)

  CPPrecondition(cp_valid(local_angles),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL pao_qs_env_get(local_angles%qs_env,&
          global_env=local_angles%global_env, pao_env=pao_env,&
          pao_max_full_bas=max_full_bas,pao_max_min_bas=max_min_bas,&
          error=error)
     CALL pao_g_ang_get(local_angles%glob_angles,qs_env=local_angles%qs_env,&
          global_env=local_angles%global_env,NUi_injection=NUi_injection,&
          error=error)
     CPPrecondition(ASSOCIATED(NUi_injection),cp_failure_level,routineP,error,failure)
     CPPrecondition(ASSOCIATED(pao_env),cp_failure_level,routineP,error,failure)
     IF (.NOT.failure) THEN
        CALL cp_assert(ASSOCIATED(local_angles%glob_angles%NUi_injection_calculated),&
             cp_failure_level,cp_assertion_failed,routineP,&
             "NUi_injection_calculated not associated in "//&
             CPSourceFileRef,&
             error,failure)
     END IF
  END IF
  IF (.NOT.failure) THEN
     IF (local_angles%glob_angles%NUi_injection_calculated&
          (local_angles%atom_nr)) THEN
        CALL get_block_node(cp_get_matrix(NUi_injection),&
             BLOCK=res,block_row=local_angles%atom_nr,&
             block_col=local_angles%atom_nr)
        CPPostcondition(ASSOCIATED(res),cp_failure_level,routineP,error,failure)
     ELSE
        CALL get_block_node(cp_get_matrix(NUi_injection),&
             BLOCK=res,block_row=local_angles%atom_nr,&
             block_col=local_angles%atom_nr)
        IF (.NOT.ASSOCIATED(res)) THEN
           CALL add_block_node(cp_get_matrix(NUi_injection),&
                local_angles%atom_nr,local_angles%atom_nr,BLOCK=res)
        END IF
        CPPostcondition(ASSOCIATED(res),cp_failure_level,routineP,error,failure)
        IF (.NOT.failure) THEN
           IF (.NOT.ASSOCIATED(local_angles%tmp_m)) THEN
              ALLOCATE(local_angles%tmp_m(max_full_bas,max_min_bas),stat=stat)
              CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           END IF
           CALL cp_assert(SIZE(local_angles%tmp_m,1) >= SIZE(res,1) .AND.&
                SIZE(local_angles%tmp_m,2) >= SIZE(res,2),cp_failure_level,&
                cp_internal_error, routineP, &
                "local_angles%tmp_m not big enough, in"//&
                CPSourceFileRef,&
                error=error,failure=failure)
        END IF
        u_m => pao_l_ang_get_unitary_t(local_angles, error=error)
        CPPrecondition(ASSOCIATED(u_m),cp_failure_level,routineP,error,failure)
        n_m => pao_l_ang_get_atomic_ortho(local_angles, error=error)
        CPPrecondition(ASSOCIATED(n_m),cp_failure_level,routineP,error,failure)
        IF (.NOT.failure) THEN
           CPPrecondition(SIZE(u_m,1)==SIZE(n_m,2),cp_failure_level,routineP,error,failure)
           CALL cp_proj_get(local_angles%projection,proj_indexes=proj_indexes,&
                error=error)
           CPPrecondition(ASSOCIATED(proj_indexes),cp_failure_level,routineP,error,failure)
        END IF
        IF (.NOT.failure) THEN
           DO i=1,SIZE(proj_indexes)
              local_angles%tmp_m(1:SIZE(u_m,1),i)=u_m(:,proj_indexes(i))
           END DO
           CALL dgemm('N','N',SIZE(n_m,1),SIZE(proj_indexes),&
                SIZE(u_m,1),1.0_dp,&
                n_m,SIZE(n_m,1),local_angles%tmp_m,SIZE(local_angles%tmp_m,1),&
                0.0_dp,res,SIZE(res,1))
           IF (cp_debug.AND.debug_this_module) THEN
              iter=iter+1
              logger => cp_error_get_logger(error)
              CALL cp_2d_r_output(logger,&
                   outputName="atomic_NUi",&
                   fromWhere=routineP, iter=iter,array=res,&
                   comment="atom_nr="//cp_to_string(local_angles%atom_nr),&
                   error=error)
              local_angles%glob_angles%NUi_injection_calculated &
                   (local_angles%atom_nr)=.TRUE.
           END IF
        END IF
     END IF
  END IF
  IF (failure) NULLIFY(res)
END FUNCTION pao_l_ang_get_NUi_injection
!***************************************************************************

!!****f* pao_l_angles_methods/pao_guarantee_owned_NUi [1.0] *
!!
!!   NAME
!!     pao_guarantee_owned_NUi
!!
!!   FUNCTION
!!     guarantees that every processor has the blocks of NUi he owns
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - global_angles: the global angles that has the NUi that must be
!!       guaranteed
!!     - qs_env: the qs enironement in which the global_angles lives
!!     - global_env: the global environement
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE pao_guarantee_owned_NUi(glob_angles, qs_env,&
     global_env,error)
    TYPE(pao_glob_angles_type), &
      INTENT(inout), TARGET                  :: glob_angles
    TYPE(qs_environment_type), INTENT(in), &
      TARGET                                 :: qs_env
    TYPE(global_environment_type), &
      INTENT(in), TARGET                     :: global_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_guarantee_owned_NUi', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, n_block_rows
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: NUi_block
    TYPE(cp_block_matrix_type), POINTER      :: NUi_injection
    TYPE(pao_local_angles_type)              :: local_angles

  NULLIFY(NUi_block, NUi_injection)
  failure=.FALSE.
  
  CPPrecondition(glob_angles%initialized,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL pao_g_ang_get(glob_angles,qs_env=qs_env,&
          global_env=global_env,NUi_injection=NUi_injection,&
          error=error)

     CPPrecondition(ASSOCIATED(NUi_injection),cp_failure_level,routineP,error,failure)
     CALL pao_l_ang_init(local_angles,glob_angles=glob_angles, qs_env=qs_env,&
          global_env=global_env, atom_nr=1, error=error)
     CALL cp_sp_get(NUi_injection, n_block_rows=n_block_rows,error=error)
     DO i=1,n_block_rows
        IF (NUi_injection%matrix_struct%distribution%block_owner(i,i) == &
             NUi_injection%matrix_struct%distribution%para_env%mepos) THEN
           CALL pao_l_ang_set(local_angles, atom_nr=i,error=error)
           NUi_block => pao_l_ang_get_NUi_injection(local_angles,error=error)
        END IF
     END DO
  END IF
END SUBROUTINE pao_guarantee_owned_NUi
!***************************************************************************

END MODULE pao_l_angles_methods
