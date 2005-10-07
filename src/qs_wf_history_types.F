!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/qs_wf_history_types [1.0] *
!!
!!   NAME
!!     qs_wf_history_types
!!
!!   FUNCTION
!!     interpolate the wavefunctions to speed up the convergence when 
!!     doing MD
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     12.2002 created [fawzi]
!!     02.2005 wf_mol added [MI]
!!
!!   SOURCE
!****************************************************************************
MODULE qs_wf_history_types
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_fm_types,                     ONLY: cp_fm_p_type,&
                                             cp_fm_release
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE kg_gpw_fm_mol_types,             ONLY: kg_fm_mol_set_release,&
                                             kg_fm_p_type
  USE kinds,                           ONLY: dp
  USE pw_types,                        ONLY: pw_p_type
  USE qs_rho_types,                    ONLY: qs_rho_type, &
                                             qs_rho_release
  USE sparse_matrix_types,             ONLY: deallocate_matrix,&
                                             deallocate_matrix_set,&
                                             real_matrix_p_type,&
                                             real_matrix_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE termination,                     ONLY: stop_program
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qs_wf_history_types'

  INTEGER, PARAMETER, PUBLIC :: wfi_use_guess_method_nr=0,&
       wfi_use_prev_p_method_nr=1, wfi_use_prev_rho_r_method_nr=2,&
       wfi_linear_wf_method_nr=3, wfi_linear_p_method_nr=4,&
       wfi_linear_ps_method_nr=5, wfi_use_prev_wf_method_nr=6,&
       wfi_ps_method_nr=7, wfi_frozen_method_nr=8, wfi_aspc=9

  PUBLIC :: qs_wf_snapshot_type, qs_wf_snapshot_p_type, &
       qs_wf_history_type
  PUBLIC :: wfs_retain, wfs_release, wfi_retain, wfi_release, wfi_get_snapshot
!!***
!****************************************************************************

!!****s* qs_wf_history_types/qs_wf_snapshot_type [1.0] *
!!
!!   NAME
!!     qs_wf_snapshot_type
!!
!!   FUNCTION
!!     represent a past snapshot of the wavefunction.
!!     some elements might not be associated (to spare memory)
!!     depending on how the snapshot was taken
!!
!!   NOTES
!!     keep trak also of occupation numbers and energies?
!!
!!   ATTRIBUTES
!!     - wf: the wavefunctions
!!     - wf_mol: the wavefunctions in a block diagonal form
!!     - rho_r: the density in r space
!!     - rho_g: the density in g space
!!     - rho_ao: the density in ao space
!!     - overlap: the overlap matrix
!!     - rho_frozen: the frozen density structure
!!     - dt: the time of the snapshot (wrf to te previous snapshot!)
!!     - id_nr: unique identification number
!!     - ref_count: reference count (see doc/ReferenceCounting.html)
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!     02.2005 wf_mol added [MI]
!!
!!   SOURCE
!***************************************************************************
  TYPE qs_wf_snapshot_type
     TYPE(cp_fm_p_type), DIMENSION(:), POINTER :: wf
     TYPE(kg_fm_p_type),  POINTER              :: wf_mol 
     TYPE(pw_p_type), DIMENSION(:), POINTER :: rho_r
     TYPE(pw_p_type), DIMENSION(:), POINTER :: rho_g
     TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: rho_ao
     TYPE(real_matrix_type), POINTER :: overlap
     TYPE(qs_rho_type), POINTER :: rho_frozen
     REAL(KIND = dp) :: dt
     INTEGER :: id_nr, ref_count
  END TYPE qs_wf_snapshot_type
!!***
!****************************************************************************

!!****s* qs_wf_history_types/qs_wf_snapshot_p_type [1.0] *
!!
!!   NAME
!!     qs_wf_snapshot_p_type
!!
!!   FUNCTION
!!     pointer to a snapshot
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - snapshot: the pointer to the snapshot
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!   SOURCE
!***************************************************************************
  TYPE qs_wf_snapshot_p_type
     TYPE(qs_wf_snapshot_type), POINTER :: snapshot
  END TYPE qs_wf_snapshot_p_type
!!***
!****************************************************************************

!!****s* qs_wf_history_types/qs_wf_history_type [1.0] *
!!
!!   NAME
!!     qs_wf_history_type
!!
!!   FUNCTION
!!     keeps track of the previous wavefunctions and can extrapolate them
!!     for the next step of md
!!
!!   NOTES
!!     use a linked list for the past states ?
!!
!!   ATTRIBUTES
!!     - id_nr: unique identification number
!!     - ref_cont: reference count (see doc/ReferenceCounting.html)
!!     - memory_depth: how many snapshots should be stored
!!     - last_state_index: index of the latest snapshot
!!     - past_states: array with the past states (index starts at
!!       last_state_index)
!!     - interpolation_method_nr: the tag of the method used to
!!       extrapolate the new start state for qs
!!     - snapshot_count: number of snapshot taken so far (cumulative,
!!       can be bigger than the history depth)
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!   SOURCE
!***************************************************************************
  TYPE qs_wf_history_type
     INTEGER :: id_nr, ref_count, memory_depth, last_state_index, &
          interpolation_method_nr, snapshot_count
     LOGICAL :: store_wf, store_rho_r, store_rho_g, store_rho_ao,&
          store_overlap, store_frozen_density
     TYPE(qs_wf_snapshot_p_type), DIMENSION(:), POINTER :: past_states
  END TYPE qs_wf_history_type
!!***
!****************************************************************************

CONTAINS

!!****f* qs_wf_history_types/wfs_retain [1.0] *
!!
!!   NAME
!!     wfs_retain
!!
!!   SYNOPSIS
!!     Subroutine wfs_retain(snapshot, error)
!!       Type(qs_wf_snapshot_type), Pointer:: snapshot
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine wfs_retain
!!
!!   FUNCTION
!!     retains a snapshot of a wavefunction(see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - snapshot: the snapshot to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE wfs_retain(snapshot,error)
    TYPE(qs_wf_snapshot_type), POINTER       :: snapshot
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'wfs_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(snapshot),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     snapshot%ref_count=snapshot%ref_count+1
  END IF
END SUBROUTINE wfs_retain
!***************************************************************************

!!****f* qs_wf_history_types/wfs_release [1.0] *
!!
!!   NAME
!!     wfs_release
!!
!!   SYNOPSIS
!!     Subroutine wfs_release(snapshot, error)
!!       Type(qs_wf_snapshot_type), Pointer:: snapshot
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine wfs_release
!!
!!   FUNCTION
!!     releases a snapshot of a wavefunction (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - snapshot: the snapshot to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!     02.2005 wf_mol added [MI]
!!
!!*** **********************************************************************
SUBROUTINE wfs_release(snapshot,error)
    TYPE(qs_wf_snapshot_type), POINTER       :: snapshot
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'wfs_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, icount, stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (ASSOCIATED(snapshot)) THEN
     CPPreconditionNoFail(snapshot%ref_count>0,cp_failure_level,routineP,error)
     snapshot%ref_count=snapshot%ref_count-1
     IF (snapshot%ref_count==0) THEN
        IF (ASSOCIATED(snapshot%wf)) THEN
           DO i=1,SIZE(snapshot%wf)
              CALL cp_fm_release(snapshot%wf(i)%matrix,error=error)
           END DO
           DEALLOCATE(snapshot%wf,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(snapshot%wf_mol)) THEN
           CALL kg_fm_mol_set_release(snapshot%wf_mol,error)
        END IF
        IF (ASSOCIATED(snapshot%rho_ao)) THEN
           CALL deallocate_matrix_set(snapshot%rho_ao)
        END IF
        IF (ASSOCIATED(snapshot%overlap)) THEN
           CALL deallocate_matrix(snapshot%overlap)
        END IF
        IF (ASSOCIATED(snapshot%rho_frozen)) THEN
           CALL qs_rho_release(snapshot%rho_frozen)
        END IF
        DEALLOCATE(snapshot,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(snapshot)
END SUBROUTINE wfs_release
!***************************************************************************

!!****f* qs_wf_history_types/wfi_retain [1.0] *
!!
!!   NAME
!!     wfi_retain
!!
!!   SYNOPSIS
!!     Subroutine wfi_retain(wf_history, error)
!!       Type(qs_wf_history_type), Pointer:: wf_history
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine wfi_retain
!!
!!   FUNCTION
!!     retains a wf history (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - wf_history: the wf_history to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE wfi_retain(wf_history,error)
    TYPE(qs_wf_history_type), POINTER        :: wf_history
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'wfi_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(wf_history),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     wf_history%ref_count=wf_history%ref_count+1
  END IF
END SUBROUTINE wfi_retain
!***************************************************************************

!!****f* qs_wf_history_types/wfi_release [1.0] *
!!
!!   NAME
!!     wfi_release
!!
!!   SYNOPSIS
!!     Subroutine wfi_release(wf_history, error)
!!       Type(qs_wf_history_type), Pointer:: wf_history
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine wfi_release
!!
!!   FUNCTION
!!     releases a wf_history of a wavefunction 
!!     (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - wf_history: the wf_history to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE wfi_release(wf_history,error)
    TYPE(qs_wf_history_type), POINTER        :: wf_history
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'wfi_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (ASSOCIATED(wf_history)) THEN
     CPPreconditionNoFail(wf_history%ref_count>0,cp_failure_level,routineP,error)
     wf_history%ref_count=wf_history%ref_count-1
     IF (wf_history%ref_count==0) THEN
        IF (ASSOCIATED(wf_history%past_states)) THEN
           DO i=1,SIZE(wf_history%past_states)
              CALL wfs_release(wf_history%past_states(i)%snapshot,&
                   error=error)
           END DO
           DEALLOCATE(wf_history%past_states,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        DEALLOCATE(wf_history,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(wf_history)
END SUBROUTINE wfi_release
!***************************************************************************

!!****f* qs_wf_history_types/wfi_get_snapshot [1.0] *
!!
!!   NAME
!!     wfi_get_snapshot
!!
!!   SYNOPSIS
!!     Function wfi_get_snapshot(wf_history, index, error) Result(res)
!!       Type(qs_wf_history_type), Pointer:: wf_history
!!       Integer, Intent (IN):: index
!!       Type(qs_wf_snapshot_type), Pointer:: res
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function wfi_get_snapshot
!!
!!   FUNCTION
!!     returns a snapshot, the first being the latest snapshot
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - wf_history: the plage where to get the snapshot
!!     - index: the index of the snapshot you want
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     12.2002 created [fawzi]
!!
!!*** **********************************************************************
FUNCTION wfi_get_snapshot(wf_history, index, error) RESULT(res)
    TYPE(qs_wf_history_type), POINTER        :: wf_history
    INTEGER, INTENT(in)                      :: index
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    TYPE(qs_wf_snapshot_type), POINTER       :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'wfi_get_snapshot', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  NULLIFY(res)
  
  CPPrecondition(ASSOCIATED(wf_history),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(ASSOCIATED(wf_history%past_states),cp_failure_level,routineP,error,failure)
     IF (index>wf_history%memory_depth.OR.index>wf_history%snapshot_count) THEN
        CPPrecondition(.FALSE.,cp_warning_level,routineP,error,failure)
     END IF
  END IF
  IF (.not.failure) THEN
     res => wf_history%past_states(&
          MODULO(wf_history%snapshot_count+1-index,&
          wf_history%memory_depth)+1)%snapshot
  END IF
END FUNCTION wfi_get_snapshot
!***************************************************************************

!!****f* qs_wf_history_types/wfi_clear_history [1.0] *
!!
!!   NAME
!!     wfi_clear_history
!!
!!   FUNCTION
!!     clears the history taken so far
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - wf_history: the history to clear
!!     - flush_buffer: if the memory of the buffer should be released
!!       defaults to false
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     05.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE wfi_clear_history(wf_history,flush_buffer,error)
    TYPE(qs_wf_history_type), POINTER        :: wf_history
    LOGICAL, INTENT(in), OPTIONAL            :: flush_buffer
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'wfi_clear_history', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure, my_flush_buffer

  failure=.FALSE.
  my_flush_buffer=.FALSE.
  IF (PRESENT(flush_buffer)) my_flush_buffer=flush_buffer
  
  CPPrecondition(ASSOCIATED(wf_history),cp_failure_level,routineP,error,failure)
  CPPrecondition(wf_history%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     wf_history%snapshot_count=0
     IF (my_flush_buffer) THEN
        DO i=1,SIZE(wf_history%past_states)
           CALL wfs_release(wf_history%past_states(i)%snapshot,error=error)
        END DO
     END IF
  END IF
END SUBROUTINE wfi_clear_history
!***************************************************************************

END MODULE qs_wf_history_types
