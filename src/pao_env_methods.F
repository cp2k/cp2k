!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/pao_env_methods [1.0] *
!!
!!   NAME
!!     pao_env_methods
!!
!!   FUNCTION
!!     methods to access and modify the pao environement
!!
!!   NOTES
!!     the get methods are defined in pao_types, because of the circular
!!     dependencies in the types
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     02.2002 created [fawzi]
!!
!!   SOURCE
!*****************************************************************************
MODULE pao_env_methods
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE cp_b_matrix_structure,           ONLY: cp_bmstruct_create,&
                                             cp_bmstruct_release
  USE cp_block_matrix,                 ONLY: cp_block_matrix_p_type,&
                                             cp_release,&
                                             cp_sp_create
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE pao_bas_projection_types,        ONLY: cp_proj_get,&
                                             pao_bas_projection_type
  USE pao_glob_angles_methods,         ONLY: cp_dealloc,&
                                             cp_set,&
                                             pao_g_ang_init
  USE pao_types,                       ONLY: pao_env_type,&
                                             pao_env_valid,&
                                             pao_glob_angles_type
  USE particle_list_types,             ONLY: particle_list_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  LOGICAL, PRIVATE, PARAMETER :: standalone=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='pao_env_methods'

! core procedures
  PUBLIC :: cp_init, cp_dealloc_ref, cp_set

! core functions
  PUBLIC :: cp_validate

! underlying routines
  PUBLIC :: pao_env_init, pao_env_dealloc_ref, pao_env_validate,&
       pao_env_set

! ============== interfaces ===========

! +++ begin of core methods +++

!! initializes the first, argument.
!! Support optional initial values that depend on the 
!! first argument 
  INTERFACE cp_init
     MODULE PROCEDURE pao_env_init
  END INTERFACE

!! Deallocates the memory allocated by the first argument.
!! Supports error=error
  INTERFACE cp_dealloc_ref
     MODULE PROCEDURE pao_env_dealloc_ref
  END INTERFACE

!! returns true if the content of the first argumenst is valid.
!! writes the errors ar warnings, full validation (long)
  INTERFACE cp_validate
     MODULE PROCEDURE pao_env_validate
  END INTERFACE

!! sets the value of various optonal attributes of the firs argument.
!! the optional attributes depend on the type of the first
  INTERFACE cp_set
     MODULE PROCEDURE pao_env_set
  END INTERFACE

! +++ end of the core methods +++

!!***
!****************************************************************************

CONTAINS

! ====== core methods =====

!!****f* pao_env_methods/pao_env_init [1.0] *
!!
!!   NAME
!!     pao_env_init
!!
!!   SYNOPSIS
!!     Subroutine pao_env_init(pao_env, qs_env, global_env,&
!!         wanted_subspace_pos_precision, wanted_in_subspace_precision,&
!!         error)
!!       Type(pao_env_type), Intent (OUT):: pao_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Type(qs_environment_type), Intent (INOUT), Target:: qs_env
!!       Real(KIND = dp), Optional, Intent (IN)::&
!!         wanted_in_subspace_precision, wanted_subspace_pos_precision
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_env_init
!!
!!   FUNCTION
!!     initialize the pao environement
!!
!!   NOTES
!!     every environement HAS to be initialized with this routine before 
!!     using.
!!     Unless said all arguments are COPIED if nedeed, not shared.
!!
!!   INPUTS
!!     - pao_env: the environement that has to be initialized
!!     - qs_env: the quickstep environement that contains this pao_env
!!     - global_env: the global environement
!!     - wanted_subspace_pos_precision: the wanted precision to find 
!!       the subspace
!!     - wanted_in_subspace_precision: the wanted precision in the subspace
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     12.2001 created [fawzi]
!!     05.2002 updated doc [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE pao_env_init(pao_env, qs_env, global_env,&
       wanted_subspace_pos_precision,&
       wanted_in_subspace_precision,&
       error)
    TYPE(pao_env_type), INTENT(out)          :: pao_env
    TYPE(qs_environment_type), &
      INTENT(inout), TARGET                  :: qs_env
    TYPE(global_environment_type), &
      INTENT(in), TARGET                     :: global_env
    REAL(KIND=dp), INTENT(in), OPTIONAL      :: wanted_subspace_pos_precision,&
                                                wanted_in_subspace_precision
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineP = 'pao_env_methods:pao_env_init'

    CHARACTER(len=default_string_length)     :: name
    INTEGER                                  :: i, j, n_atoms, size_full_bas, &
                                                size_min_bas, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: full_bas_sizes, min_bas_sizes
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: failure
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_error_type)                      :: iError
    TYPE(cp_subsystem_type), POINTER         :: subsys
    TYPE(pao_bas_projection_type), POINTER   :: bproj
    TYPE(particle_list_type), POINTER        :: particles

    NULLIFY(atom_list, subsys, particles, atomic_kinds)
    failure=.FALSE.
    CALL cp_error_init(iError, template_error=error)

    NULLIFY(pao_env%angles_begin_at, pao_env%atomic_ortho_valid,&
         pao_env%atomic_ortho, pao_env%min_density_m,&
         pao_env%min_bas_m_struct, pao_env%full_bas_m_struct,&
         pao_env%full_min_bas_m_struct, pao_env%min_full_bas_m_struct,&
         pao_env%full_bas_diag_m_struct)
    pao_env%min_d_m_valid=.FALSE.

    ! check projection
    CALL get_qs_env(qs_env, subsys=subsys)
    CALL cp_subsys_get(subsys, particles=particles, atomic_kinds=atomic_kinds,&
         error=error)
    n_atoms=particles%n_els
    CPPrecondition(ASSOCIATED(atomic_kinds),cp_failure_level,routineP,error,failure)
    IF (.not.failure) THEN
       ALLOCATE(min_bas_sizes(n_atoms), full_bas_sizes(n_atoms),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF (.not.failure) THEN
       DO i=1,atomic_kinds%n_els
          atomic_kind => atomic_kinds%els(i)
          CALL get_atomic_kind(atomic_kind,&
                               pao_bas_proj=bproj,&
                               atom_list=atom_list)
          IF (.not.ASSOCIATED(bproj)) THEN
             CALL get_atomic_kind(atomic_kind,name=name)
             CALL cp_assert(.FALSE.,cp_failure_level, cp_internal_error,&
                  routineP,&
                  "projection (PAO_MIN_BAS) not defined in kind "//name &
                  //", "//&
                  CPSourceFileRef,&
                  error=error,failure=failure)
          ELSE
             CALL cp_proj_get(bproj,size_min_bas=size_min_bas,&
                  size_full_bas=size_full_bas,error=error)
             DO j=1,SIZE(atom_list)
                min_bas_sizes(atom_list(j))=size_min_bas
                full_bas_sizes(atom_list(j))=size_full_bas
             END DO
          END IF
       END DO
       
       ! min_bas_m_struct / full_bas_m_struct
       IF (.not.failure) THEN
          CALL cp_bmstruct_create(pao_env%min_bas_m_struct,&
               symmetric=.TRUE.,&
               para_env=global_env%para_env,&
               b_row_sizes=min_bas_sizes,&
               b_col_sizes=min_bas_sizes,&
               error=error)
          CALL cp_bmstruct_create(pao_env%full_bas_m_struct,&
               symmetric=.FALSE.,&
               para_env=global_env%para_env,&
               b_row_sizes=full_bas_sizes,&
               b_col_sizes=full_bas_sizes,&
               error=error)
          CALL cp_bmstruct_create(pao_env%full_min_bas_m_struct,&
               symmetric=.FALSE.,&
               para_env=global_env%para_env,&
               b_row_sizes=full_bas_sizes,&
               b_col_sizes=min_bas_sizes,&
               error=error)
          CALL cp_bmstruct_create(pao_env%min_full_bas_m_struct,&
               symmetric=.FALSE.,&
               para_env=global_env%para_env,&
               b_row_sizes=min_bas_sizes,&
               b_col_sizes=full_bas_sizes,&
               error=error)
          CALL cp_bmstruct_create(pao_env%full_bas_diag_m_struct,&
               symmetric=.FALSE.,&
               para_env=global_env%para_env,&
               b_row_sizes=full_bas_sizes,&
               b_col_sizes=full_bas_sizes,&
               error=error)
          pao_env%full_bas_diag_m_struct%sparsity=.FALSE.
          DO i=1,SIZE(pao_env%full_bas_diag_m_struct%sparsity,1)
             pao_env%full_bas_diag_m_struct%sparsity(i,i)=.TRUE.
          END DO
       END IF
       ! angles_begin_at
       IF (.not.failure) THEN
          ALLOCATE(pao_env%angles_begin_at(n_atoms+1),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF (.not.failure) THEN
          IF (n_atoms>0) pao_env%angles_begin_at(1)=0
          DO i=1,n_atoms
             pao_env%angles_begin_at(i+1)=pao_env%angles_begin_at(i)+&
                  min_bas_sizes(i)*(full_bas_sizes(i)-min_bas_sizes(i))
          END DO
       END IF

       DEALLOCATE(min_bas_sizes, full_bas_sizes, stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    END IF

    ! cache_max vals
    CALL get_atomic_kind_set(atomic_kinds%els,&
         pao_max_min_bas=pao_env%max_min_bas,&
         pao_max_full_bas=pao_env%max_full_bas,&
         pao_max_excl_bas=pao_env%max_excl_bas,&
         pao_tot_min_bas=pao_env%tot_min_bas,&
         pao_tot_full_bas=pao_env%tot_full_bas,&
         pao_tot_excl_bas=pao_env%tot_excl_bas,&
         pao_max_l_angles=pao_env%max_l_angles)
        
    ! others
    n_atoms=particles%n_els
    ALLOCATE(pao_env%min_density_m(1),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.not.failure) THEN
       DO i=1,SIZE(pao_env%min_density_m)
          NULLIFY(pao_env%min_density_m(i)%matrix)
       END DO
    END IF
    ALLOCATE(pao_env%angles_att,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    pao_env%wanted_subspace_pos_precision=1.0e-4_dp
    pao_env%wanted_in_subspace_precision=1.0e-4_dp
    pao_env%initialized=.TRUE.
    ! min_d
    IF (.not.failure) THEN
       DO i=1,SIZE(pao_env%min_density_m)
          CALL cp_sp_create(pao_env%min_density_m(i)%matrix,&
               matrix_struct=pao_env%min_bas_m_struct,&
               error=error)
       END DO
    END IF
    ! angles_att
    IF (.not.failure) THEN
       CALL pao_g_ang_init(pao_env%angles_att, global_env=global_env,&
            qs_env=qs_env, error=error)
    END IF

    CALL pao_env_set(pao_env,qs_env=qs_env, global_env=global_env,&
         wanted_subspace_pos_precision=wanted_subspace_pos_precision,&
         wanted_in_subspace_precision=wanted_in_subspace_precision,&
         error=error)
    IF (cp_debug.and.debug_this_module) THEN
       CALL cp_assert(cp_validate(pao_env,qs_env=qs_env,&
            global_env=global_env,&
            test_pos_att=.FALSE.,error=error),&
            cp_failure_level,cp_assertion_failed,routineP,&
            "pao_env is not valid after init in "//&
            CPSourceFileRef,&
            error,failure)
    END IF
    CALL cp_error_dealloc_ref(iError,error)
  END SUBROUTINE pao_env_init
!***************************************************************************

!!****f* pao_env_methods/pao_env_dealloc_ref [1.0] *
!!
!!   NAME
!!     pao_env_dealloc_ref
!!
!!   SYNOPSIS
!!     Subroutine pao_env_dealloc_ref(pao_env, qs_env, global_env, error)
!!       Type(pao_env_type), Intent (INOUT):: pao_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Type(qs_environment_type), Intent (INOUT), Target:: qs_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_env_dealloc_ref
!!
!!   FUNCTION
!!     deallocate the memory that the pao environement has allocated,
!!     but NOT the pao environement itself.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pao_env: the environement to be deallocated
!!     - qs_env: the quickstep environement that contains this pao_env
!!     - global_env: the global environement
!!     - error: variable to control error logging, stopping,... 
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
  SUBROUTINE pao_env_dealloc_ref(pao_env,qs_env,global_env, error)
    TYPE(pao_env_type), INTENT(inout)        :: pao_env
    TYPE(qs_environment_type), &
      INTENT(inout), TARGET                  :: qs_env
    TYPE(global_environment_type), &
      INTENT(in), TARGET                     :: global_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'pao_env_methods:pao_env_dealloc_ref'

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(pao_env%initialized,cp_failure_level,routineP,error,failure)
    failureIf: IF (.NOT. failure) THEN
       IF (ASSOCIATED(pao_env%angles_att)) THEN
          CALL cp_dealloc(pao_env%angles_att,global_env=global_env,&
               error=error)
       END IF
       IF (ASSOCIATED(pao_env%angles_begin_at)) THEN
          DEALLOCATE(pao_env%angles_begin_at,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(pao_env%min_bas_m_struct)) THEN
          CALL cp_bmstruct_release(pao_env%min_bas_m_struct,error=error)
          NULLIFY(pao_env%min_bas_m_struct)
       END IF
       IF (ASSOCIATED(pao_env%full_bas_m_struct)) THEN
          CALL cp_bmstruct_release(pao_env%full_bas_m_struct, error=error)
          NULLIFY(pao_env%full_bas_m_struct)
       END IF
       IF (ASSOCIATED(pao_env%full_min_bas_m_struct)) THEN
          CALL cp_bmstruct_release(pao_env%full_min_bas_m_struct, error=error)
          NULLIFY(pao_env%full_min_bas_m_struct)
       END IF
       IF (ASSOCIATED(pao_env%min_full_bas_m_struct)) THEN
          CALL cp_bmstruct_release(pao_env%min_full_bas_m_struct, error=error)
          NULLIFY(pao_env%min_full_bas_m_struct)
       END IF
       IF (ASSOCIATED(pao_env%full_bas_diag_m_struct)) THEN
          CALL cp_bmstruct_release(pao_env%full_bas_diag_m_struct, error=error)
          NULLIFY(pao_env%full_bas_diag_m_struct)
       END IF
       IF (ASSOCIATED(pao_env%atomic_ortho)) THEN
          CALL cp_release(pao_env%atomic_ortho,error=error)
       END IF
       IF (ASSOCIATED(pao_env%atomic_ortho_valid)) THEN
          DEALLOCATE(pao_env%atomic_ortho_valid, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
       pao_env%initialized=.FALSE.
    END IF failureIf
  END SUBROUTINE pao_env_dealloc_ref
!***************************************************************************


!!****f* pao_env_methods/pao_env_validate [1.0] *
!!
!!   NAME
!!     pao_env_validate
!!
!!   SYNOPSIS
!!     Function pao_env_validate(pao_env, qs_env, global_env, test_pos_att,&
!!         error) Result(res)
!!       Logical:: res
!!       Type(pao_env_type), Intent (IN), Target:: pao_env
!!       Type(qs_environment_type), Intent (IN), Target:: qs_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Logical, Optional:: test_pos_att
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function pao_env_validate
!!
!!   FUNCTION
!!     returns true if the pao environement is valid
!!
!!   NOTES
!!     Writes out all the failures as warning, use error to change
!!     this behaviour
!!
!!   INPUTS
!!     - pao_env: the environement to be tested
!!     - qs_env: the quickstep environement that contains this pao_env
!!     - global_env: the global environement
!!     - test_pos_att: if it is false does not test the validity of the
!!       actual angles and unitary transformation. defaults to true
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  FUNCTION pao_env_validate(pao_env,qs_env,global_env, test_pos_att, error)&
       RESULT(res)
    TYPE(pao_env_type), INTENT(in), TARGET   :: pao_env
    TYPE(qs_environment_type), INTENT(in), &
      TARGET                                 :: qs_env
    TYPE(global_environment_type), &
      INTENT(in), TARGET                     :: global_env
    LOGICAL, OPTIONAL                        :: test_pos_att
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'pao_env_methods:pao_env_validate'

    LOGICAL                                  :: failure, test_pos

    failure=.FALSE.
    IF (PRESENT(test_pos_att)) THEN
       test_pos=test_pos_att
    ELSE
       test_pos=.TRUE.
    END IF
    CPPrecondition(pao_env%initialized,cp_failure_level,routineP,error,failure)
    failure=failure.and..not.pao_env_valid(pao_env,error=error)
    ! to do: more tests using qs_env
    res=.not.failure
  END FUNCTION pao_env_validate
!***************************************************************************

!!****f* pao_env_methods/pao_env_set [1.0] *
!!
!!   NAME
!!     pao_env_set
!!
!!   SYNOPSIS
!!     Subroutine pao_env_set(pao_env, qs_env, global_env,&
!!         wanted_subspace_pos_precision, wanted_in_subspace_precision,&
!!         angles_att, min_density_m, min_d_valid, error)
!!       Type(pao_env_type), Target, Intent (INOUT):: pao_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Type(qs_environment_type), Intent (IN), Target:: qs_env
!!       Real(KIND = dp), Optional, Intent (IN)::&
!!         wanted_in_subspace_precision, wanted_subspace_pos_precision
!!       Type(pao_glob_angles_type), Optional, Intent (IN):: angles_att
!!       Type(cp_block_matrix_p_type), Optional, Intent (IN), Dimension(:),&
!!         Target:: min_density_m
!!       Logical, Optional:: min_d_valid
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_env_set
!!
!!   FUNCTION
!!     changes some of the parameters of the environement
!!
!!   NOTES
!!     arguments, unless otherwise stated are always copied
!!
!!   INPUTS
!!     - pao_env: the environement that has to be initialized
!!     - qs_env: the quickstep environement where the pao_env lives
!!     - global_env: the global environement
!!     min_density_m: the density in the minimal basis (SHARED)
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!     for the description of the other variables see the attributes 
!!     of pao_env_type
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE pao_env_set(pao_env, qs_env, global_env, &
       wanted_subspace_pos_precision,&
       wanted_in_subspace_precision, angles_att, &
       min_density_m, min_d_valid,error)
    TYPE(pao_env_type), INTENT(inout), &
      TARGET                                 :: pao_env
    TYPE(qs_environment_type), INTENT(in), &
      TARGET                                 :: qs_env
    TYPE(global_environment_type), &
      INTENT(in), TARGET                     :: global_env
    REAL(KIND=dp), INTENT(in), OPTIONAL      :: wanted_subspace_pos_precision,&
                                                wanted_in_subspace_precision
    TYPE(pao_glob_angles_type), INTENT(in), &
      OPTIONAL                               :: angles_att
    TYPE(cp_block_matrix_p_type), &
      DIMENSION(:), INTENT(in), OPTIONAL, &
      TARGET                                 :: min_density_m
    LOGICAL, OPTIONAL                        :: min_d_valid
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineP = 'pao_env_methods:pao_env_set'

    LOGICAL                                  :: failure
    TYPE(cp_error_type)                      :: iError

    failure=.FALSE.
    CALL cp_error_init(iError,template_error=error)

    CPPrecondition(pao_env%initialized,cp_failure_level,routineP,error,failure)
    IF (PRESENT(wanted_subspace_pos_precision)) &
         pao_env%wanted_subspace_pos_precision=wanted_subspace_pos_precision
    IF (PRESENT(wanted_in_subspace_precision)) &
         pao_env%wanted_in_subspace_precision=wanted_in_subspace_precision
    IF (PRESENT(angles_att)) THEN
       CALL cp_set(pao_env%angles_att, qs_env=qs_env, global_env=global_env,&
            source=angles_att, error=error)
    END IF
    IF (PRESENT(min_density_m)) THEN
       pao_env%min_d_m_valid=.TRUE.
       ! to do
       CPAssert(.FALSE.,cp_warning_level,routineP,error,failure)
       pao_env%min_density_m => min_density_m
    END IF
    IF (PRESENT(min_d_valid)) THEN
       pao_env%min_d_m_valid=min_d_valid
    END IF
    CALL cp_error_dealloc_ref(iError,error=error)
  END SUBROUTINE pao_env_set
!***************************************************************************

END MODULE pao_env_methods
