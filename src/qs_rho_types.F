!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/qs_rho_types [1.0] *
!!
!!   NAME
!!     qs_rho_types
!!
!!   FUNCTION
!!     superstucture that hold various representations of the density and 
!!     keeps track of which ones are valid
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE qs_rho_types
  USE coefficient_types,               ONLY: coeff_deallocate,&
                                             coeff_type
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
                                             cp_unimplemented_error
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cube_utils,                      ONLY: cube_info_type
  USE gaussian_gridlevels,             ONLY: gridlevel_info_type
  USE global_types,                    ONLY: DEBUG,&
                                             global_environment_type
  USE kinds,                           ONLY: dp
  USE l_utils,                         ONLY: l_info_type
  USE sparse_matrix_types,             ONLY: add_matrices,&
                                             deallocate_matrix_set,&
                                             first_block_node,&
                                             get_block_node,&
                                             get_matrix_info,&
                                             next_block_node,&
                                             real_block_node_type,&
                                             real_matrix_p_type,&
                                             real_matrix_type,&
                                             replicate_matrix_structure,&
                                             set_matrix,&
                                             symmetrise_diagonal_blocks,&
                                             transfer_matrix
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qs_rho_types'

  PUBLIC :: qs_rho_type
  PUBLIC :: qs_rho_retain, qs_rho_release, qs_rho_did_change,&
       qs_rho_get

!!***
!****************************************************************************

!!****s* qs_rho_types/qs_rho_type [1.0] *
!!
!!   NAME
!!     qs_rho_type
!!
!!   FUNCTION
!!     keeps the density in various representations, keeping track of
!!     which ones are valid.
!!
!!   NOTES
!!     If coeff_type would implement retain/release it would be nice to
!!     store also the core charge and the qs_charges in this structure...
!!
!!   ATTRIBUTES
!!     most attributes are array with either lda or lsd_alpha,lsd_beta.
!!     - rho_ao: the filterd rho in the localized atom basis (to have rho(r)
!!       the filtred matrix is enough, but rho(r,r') is lost).
!!     - rho_r: grids with rho in the real space
!!     - rho_g: grids with rho in the g space
!!     - rho_ao_valid, rho_g_valid, rho_r_valid: if the corresponding
!!       component is valid
!!     - ref_count: the reference count, when it becomes 0 the type
!!       is deallocated.
!!     - updates_since_rebuild: how many updates were done since the
!!       last rebuild (-1 if rebuild is needed)
!!     - rebuild_each: how often a rebuild should be done by default
!!     - tot_rho_r: the total charge in r space (valid only if rho_r is)
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!   SOURCE
!***************************************************************************
  TYPE qs_rho_type
     TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: rho_ao
     TYPE(coeff_type), DIMENSION(:), POINTER :: rho_g, rho_r
     LOGICAL :: rho_ao_valid, &
          rho_g_valid, &
          rho_r_valid
     INTEGER :: ref_count, &
          id_nr, &
          updates_since_rebuild, &
          rebuild_each
     REAL(KIND = dp), DIMENSION(:), POINTER :: tot_rho_r, tot_rho_g
  END TYPE qs_rho_type
!!***
!****************************************************************************

CONTAINS

!!****f* qs_rho_methods/qs_rho_retain [1.0] *
!!
!!   NAME
!!     qs_rho_retain
!!
!!   SYNOPSIS
!!     Subroutine qs_rho_retain(rho_struct, error)
!!       Type(qs_rho_type), Pointer:: rho_struct
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine qs_rho_retain
!!
!!   FUNCTION
!!     retains a rho_struct by increasing the reference count by one
!!     (to be called when you what to keep a shared copy)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - rho_struct: the structure to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE qs_rho_retain(rho_struct,error)
    TYPE(qs_rho_type), POINTER :: rho_struct
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='qs_rho_retain',&
         routineP=moduleN//':'//routineN
    failure=.FALSE.

    CPPrecondition(ASSOCIATED(rho_struct),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CPPrecondition(rho_struct%ref_count>0,cp_failure_level,routineP,error,failure)
       rho_struct%ref_count=rho_struct%ref_count+1
    END IF
  END SUBROUTINE qs_rho_retain
!***************************************************************************

!!****f* qs_rho_methods/qs_rho_release [1.0] *
!!
!!   NAME
!!     qs_rho_release
!!
!!   SYNOPSIS
!!     Subroutine qs_rho_release(rho_struct, error)
!!       Type(qs_rho_type), Pointer:: rho_struct
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine qs_rho_release
!!
!!   FUNCTION
!!     releases a rho_struct by decreasing the reference count by one
!!     and deallocating if it reaces 0 (to be called when you don't whant
!!     anymore a shared copy)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - rho_struct: the structure to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE qs_rho_release(rho_struct,error)
    TYPE(qs_rho_type), POINTER :: rho_struct
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='qs_rho_release',&
         routineP=moduleN//':'//routineN
    INTEGER :: stat, i, j
    failure=.FALSE.

    IF (ASSOCIATED(rho_struct)) THEN
       CPPrecondition(rho_struct%ref_count>0,cp_failure_level,routineP,error,failure)
       rho_struct%ref_count=rho_struct%ref_count-1
       IF (rho_struct%ref_count<1) THEN
          IF (ASSOCIATED(rho_struct%rho_ao)) THEN
             CALL deallocate_matrix_set(rho_struct%rho_ao)
          END IF
          IF (ASSOCIATED(rho_struct%rho_r)) THEN
             DO i=1,SIZE(rho_struct%rho_r)
                CALL coeff_deallocate(rho_struct%rho_r(i))
             END DO
             DEALLOCATE(rho_struct%rho_r,stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(rho_struct%rho_g)) THEN
             DO i=1,SIZE(rho_struct%rho_g)
                CALL coeff_deallocate(rho_struct%rho_g(i))
             END DO
             DEALLOCATE(rho_struct%rho_g,stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          DEALLOCATE(rho_struct,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
    END IF
    NULLIFY(rho_struct)
  END SUBROUTINE qs_rho_release
!***************************************************************************

!!****f* qs_rho_methods/qs_rho_did_change [1.0] *
!!
!!   NAME
!!     qs_rho_did_change
!!
!!   SYNOPSIS
!!     Subroutine qs_rho_did_change(rho_struct, rho_ao_valid, rho_r_valid,&
!!         rho_g_valid, error)
!!       Type(qs_rho_type), Pointer:: rho_struct
!!       Logical, Intent (IN), Optional:: rho_ao_valid, rho_r_valid,&
!!         rho_g_valid
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine qs_rho_did_change
!!
!!   FUNCTION
!!     this function should be called wehn something changes to update the
!!     internal status 
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - rho_struct: the density structure that changed
!!     - rho_ao_valid,rho_r_valid, rho_g_valid: actual valid
!!       representations
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE qs_rho_did_change(rho_struct,rho_ao_valid,rho_r_valid, &
       rho_g_valid, error)
    TYPE(qs_rho_type), POINTER :: rho_struct
    LOGICAL, INTENT(in), OPTIONAL :: rho_ao_valid, rho_r_valid, rho_g_valid
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='qs_rho_did_change',&
         routineP=moduleN//':'//routineN
    failure=.FALSE.

    CPPrecondition(ASSOCIATED(rho_struct),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_struct%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(rho_ao_valid)) rho_struct%rho_ao_valid=rho_ao_valid
       IF (PRESENT(rho_r_valid)) rho_struct%rho_r_valid=rho_r_valid
       IF (PRESENT(rho_g_valid)) rho_struct%rho_g_valid=rho_g_valid
       CALL cp_assert(rho_struct%rho_ao_valid.OR.&
            rho_struct%rho_r_valid.OR.rho_struct%rho_g_valid,&
            cp_failure_level,cp_assertion_failed,fromWhere=routineP,&
            message="at least one representation must be valid "//&
            CPSourceFileRef,&
            error=error,failure=failure)
    END IF
  END SUBROUTINE qs_rho_did_change
!***************************************************************************

!!****f* qs_rho_methods/qs_rho_get [1.0] *
!!
!!   NAME
!!     qs_rho_get
!!
!!   SYNOPSIS
!!     Subroutine qs_rho_get(rho_struct, rho_ao, rho_r, rho_g, rho_ao_valid,&
!!         rho_r_valid, rho_g_valid, rebuild_each, tot_rho_r, error)
!!       Type(qs_rho_type), Pointer:: rho_struct
!!       Type(real_matrix_p_type), Dimension(:), Pointer, Optional:: rho_ao
!!       Type(coeff_type), Dimension(:), Pointer, Optional:: rho_g, rho_r
!!       Logical, Intent (OUT), Optional:: rho_ao_valid, rho_g_valid,&
!!         rho_r_valid
!!       Integer, Intent (OUT), Optional:: rebuild_each
!!       Real(KIND = dp), Dimension(:), Pointer, Optional:: tot_rho_r
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine qs_rho_get
!!
!!   FUNCTION
!!     returns info about the density described by this object.
!!     If some representation is not available an error is issued
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!     For the other inputs the qs_rho_type attributes
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE qs_rho_get(rho_struct, rho_ao, rho_r, rho_g, rho_ao_valid,&
       rho_r_valid, rho_g_valid, rebuild_each, tot_rho_r, error)
    TYPE(qs_rho_type), POINTER :: rho_struct
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER,OPTIONAL :: rho_ao
    TYPE(coeff_type), DIMENSION(:), POINTER,OPTIONAL :: rho_g, rho_r
    LOGICAL, INTENT(out), OPTIONAL :: rho_ao_valid, rho_g_valid, rho_r_valid
    INTEGER, INTENT(out), OPTIONAL :: rebuild_each
    REAL(KIND = dp), DIMENSION(:), POINTER, OPTIONAL :: tot_rho_r
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='qs_rho_get',&
         routineP=moduleN//':'//routineN
    failure=.FALSE.

    CPPrecondition(ASSOCIATED(rho_struct),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_struct%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(rho_ao)) THEN
          IF (.NOT.rho_struct%rho_ao_valid) THEN
             CALL cp_unimplemented_error(routineP,"to do",error=error)
             NULLIFY(rho_ao)
          ELSE
             rho_ao => rho_struct%rho_ao
          END IF
       END IF
       IF (PRESENT(rho_r)) THEN
          IF (PRESENT(rho_r_valid).or.rho_struct%rho_r_valid) THEN
             rho_r => rho_struct%rho_r
          ELSE
             CALL cp_unimplemented_error(routineP,"to do",error=error)
             NULLIFY(rho_r)
          END IF
       END IF
       IF (PRESENT(rho_g)) THEN
          IF (PRESENT(rho_g_valid).or.rho_struct%rho_g_valid) THEN
             rho_g => rho_struct%rho_g
          ELSE
             CALL cp_unimplemented_error(routineP,"to do",error=error)
             NULLIFY(rho_g)
          END IF
       END IF
       IF (PRESENT(rho_ao_valid)) rho_ao_valid=rho_struct%rho_ao_valid
       IF (PRESENT(rho_r_valid)) rho_r_valid=rho_struct%rho_r_valid
       IF (PRESENT(rho_g_valid)) rho_g_valid=rho_struct%rho_g_valid
       IF (PRESENT(rebuild_each)) rebuild_each=rho_struct%rebuild_each
       ! give error if not valid?
       IF (PRESENT(tot_rho_r)) tot_rho_r => rho_struct%tot_rho_r
    END IF
  END SUBROUTINE qs_rho_get
!***************************************************************************

END MODULE qs_rho_types

