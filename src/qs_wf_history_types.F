!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief interpolate the wavefunctions to speed up the convergence when
!>      doing MD
!> \par History
!>      12.2002 created [fawzi]
!>      02.2005 wf_mol added [MI]
!> \author fawzi
! *****************************************************************************
MODULE qs_wf_history_types
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_deallocate_matrix,&
                                             cp_dbcsr_deallocate_matrix_set,&
                                             cp_dbcsr_p_type,&
                                             cp_dbcsr_type
  USE cp_fm_types,                     ONLY: cp_fm_p_type,&
                                             cp_fm_release
  USE kinds,                           ONLY: dp
  USE pw_types,                        ONLY: pw_p_type
  USE qs_rho_types,                    ONLY: qs_rho_release,&
                                             qs_rho_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_wf_history_types'

  PUBLIC :: qs_wf_snapshot_type,&
            qs_wf_history_type, qs_wf_history_p_type
  PUBLIC :: wfi_retain, wfi_release, wfi_get_snapshot

! *****************************************************************************
!> \brief represent a past snapshot of the wavefunction.
!>      some elements might not be associated (to spare memory)
!>      depending on how the snapshot was taken
!> \param wf the wavefunctions
!> \param rho_r the density in r space
!> \param rho_g the density in g space
!> \param rho_ao the density in ao space
!> \param overlap the overlap matrix
!> \param rho_frozen the frozen density structure
!> \param dt the time of the snapshot (wrf to te previous snapshot!)
!> \param id_nr unique identification number
!> \param ref_count reference count (see doc/ReferenceCounting.html)
!> \note
!>      keep trak also of occupation numbers and energies?
!> \par History
!>      02.2003 created [fawzi]
!>      02.2005 wf_mol added [MI]
!> \author fawzi
! *****************************************************************************
  TYPE qs_wf_snapshot_type
     TYPE(cp_fm_p_type), DIMENSION(:), POINTER :: wf
     TYPE(pw_p_type), DIMENSION(:), POINTER :: rho_r
     TYPE(pw_p_type), DIMENSION(:), POINTER :: rho_g
     TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER :: rho_ao
     TYPE(cp_dbcsr_p_type), DIMENSION(:,:), POINTER :: rho_ao_kp
     TYPE(cp_dbcsr_type), POINTER :: overlap
     TYPE(qs_rho_type), POINTER :: rho_frozen
     REAL(KIND = dp) :: dt
     INTEGER :: id_nr, ref_count
  END TYPE qs_wf_snapshot_type

! *****************************************************************************
!> \brief pointer to a snapshot
!> \param snapshot the pointer to the snapshot
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
  TYPE qs_wf_snapshot_p_type
     TYPE(qs_wf_snapshot_type), POINTER :: snapshot
  END TYPE qs_wf_snapshot_p_type

! *****************************************************************************
!> \brief keeps track of the previous wavefunctions and can extrapolate them
!>      for the next step of md
!> \param id_nr unique identification number
!> \param ref_cont reference count (see doc/ReferenceCounting.html)
!> \param memory_depth how many snapshots should be stored
!> \param last_state_index index of the latest snapshot
!> \param past_states array with the past states (index starts at
!>        last_state_index)
!> \param interpolation_method_nr the tag of the method used to
!>        extrapolate the new start state for qs
!> \param snapshot_count number of snapshot taken so far (cumulative,
!>        can be bigger than the history depth)
!> \note
!>      use a linked list for the past states ?
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
  TYPE qs_wf_history_type
     INTEGER :: id_nr, ref_count, memory_depth, last_state_index, &
          interpolation_method_nr, snapshot_count
     LOGICAL :: store_wf, store_rho_r, store_rho_g, store_rho_ao, store_rho_ao_kp, &
          store_overlap, store_frozen_density
     TYPE(qs_wf_snapshot_p_type), DIMENSION(:), POINTER :: past_states
  END TYPE qs_wf_history_type

! *****************************************************************************
!> \brief to create arrays of pointers to qs_wf_history_type
!> \param wf_hist the pointer to the wf history
!> \author fawzi
! *****************************************************************************
  TYPE qs_wf_history_p_type
     TYPE(qs_wf_history_type), POINTER :: wf_history
  END TYPE qs_wf_history_p_type

CONTAINS

! *****************************************************************************
!> \brief releases a snapshot of a wavefunction (see doc/ReferenceCounting.html)
!> \param snapshot the snapshot to release
!> \par History
!>      02.2003 created [fawzi]
!>      02.2005 wf_mol added [MI]
!> \author fawzi
! *****************************************************************************
SUBROUTINE wfs_release(snapshot)
    TYPE(qs_wf_snapshot_type), POINTER       :: snapshot

    CHARACTER(len=*), PARAMETER :: routineN = 'wfs_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

  failure=.FALSE.

  IF (ASSOCIATED(snapshot)) THEN
     CPPreconditionNoFail(snapshot%ref_count>0,cp_failure_level,routineP)
     snapshot%ref_count=snapshot%ref_count-1
     IF (snapshot%ref_count==0) THEN
        IF (ASSOCIATED(snapshot%wf)) THEN
           DO i=1,SIZE(snapshot%wf)
              CALL cp_fm_release(snapshot%wf(i)%matrix)
           END DO
           DEALLOCATE(snapshot%wf,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP)
        END IF
        ! snapshot%rho_r & snapshot%rho_g is deallocated in wfs_update
        ! of qs_wf_history_methods, in case you wonder about it.
        IF (ASSOCIATED(snapshot%rho_ao)) THEN
           CALL cp_dbcsr_deallocate_matrix_set(snapshot%rho_ao)
        END IF
        IF (ASSOCIATED(snapshot%rho_ao_kp)) THEN
           CALL cp_dbcsr_deallocate_matrix_set(snapshot%rho_ao_kp)
        END IF
        IF (ASSOCIATED(snapshot%overlap)) THEN
           CALL cp_dbcsr_deallocate_matrix(snapshot%overlap)
        END IF
        IF (ASSOCIATED(snapshot%rho_frozen)) THEN
           CALL qs_rho_release(snapshot%rho_frozen)
        END IF
        DEALLOCATE(snapshot,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP)
     END IF
  END IF
  NULLIFY(snapshot)
END SUBROUTINE wfs_release

! *****************************************************************************
!> \brief retains a wf history (see doc/ReferenceCounting.html)
!> \param wf_history the wf_history to retain
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
SUBROUTINE wfi_retain(wf_history)
    TYPE(qs_wf_history_type), POINTER        :: wf_history

    CHARACTER(len=*), PARAMETER :: routineN = 'wfi_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(wf_history),cp_failure_level,routineP,failure)
  wf_history%ref_count=wf_history%ref_count+1

END SUBROUTINE wfi_retain

! *****************************************************************************
!> \brief releases a wf_history of a wavefunction
!>      (see doc/ReferenceCounting.html)
!> \param wf_history the wf_history to release
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
SUBROUTINE wfi_release(wf_history)
    TYPE(qs_wf_history_type), POINTER        :: wf_history

    CHARACTER(len=*), PARAMETER :: routineN = 'wfi_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

  failure=.FALSE.

  IF (ASSOCIATED(wf_history)) THEN
     CPPreconditionNoFail(wf_history%ref_count>0,cp_failure_level,routineP)
     wf_history%ref_count=wf_history%ref_count-1
     IF (wf_history%ref_count==0) THEN
        IF (ASSOCIATED(wf_history%past_states)) THEN
           DO i=1,SIZE(wf_history%past_states)
              CALL wfs_release(wf_history%past_states(i)%snapshot)
           END DO
           DEALLOCATE(wf_history%past_states,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP)
        END IF
        DEALLOCATE(wf_history,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP)
     END IF
  END IF
  NULLIFY(wf_history)
END SUBROUTINE wfi_release

! *****************************************************************************
!> \brief returns a snapshot, the first being the latest snapshot
!> \param wf_history the plage where to get the snapshot
!> \param index the index of the snapshot you want
!> \retval res ...
!> \par History
!>      12.2002 created [fawzi]
!> \author fawzi
! *****************************************************************************
FUNCTION wfi_get_snapshot(wf_history, index) RESULT(res)
    TYPE(qs_wf_history_type), POINTER        :: wf_history
    INTEGER, INTENT(in)                      :: index
    TYPE(qs_wf_snapshot_type), POINTER       :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'wfi_get_snapshot', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  NULLIFY(res)

  CPPrecondition(ASSOCIATED(wf_history),cp_failure_level,routineP,failure)
  CPPrecondition(ASSOCIATED(wf_history%past_states),cp_failure_level,routineP,failure)
  IF (index>wf_history%memory_depth.OR.index>wf_history%snapshot_count) THEN
     CPPrecondition(.FALSE.,cp_warning_level,routineP,failure)
  END IF
  res => wf_history%past_states(&
       MODULO(wf_history%snapshot_count+1-index,&
       wf_history%memory_depth)+1)%snapshot
END FUNCTION wfi_get_snapshot

END MODULE qs_wf_history_types
