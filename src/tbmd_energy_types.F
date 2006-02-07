!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2006  CP2K developers group                            !
!-----------------------------------------------------------------------------!

!!****h* cp2k/tbmd_energy_types [1.0] *
!!
!!   NAME
!!     tbmd_energy_types
!!
!!   FUNCTION
!!     The tight-binding energy type
!!
!!   NOTES
!!     -
!!     
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!   SOURCE
!****************************************************************************
MODULE tbmd_energy_types
  USE kinds,                           ONLY: dp
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tbmd_energy_types'
  INTEGER, PRIVATE, SAVE :: last_tbmd_id=0

  ! *** Public data types ***
  PUBLIC :: tbmd_energy_type, &
            tbmd_energy_p_type

  ! *** Public subroutines ***
  PUBLIC :: tbmd_energy_retain, &
            tbmd_energy_release, &
            tbmd_energy_set,&
            tbmd_energy_get, &
            tbmd_energy_clear, &
            tbmd_energy_create
!!***
!****************************************************************************

!!****s* tbmd_energy_types/tbmd_energy_type [1.0] *
!!
!!   NAME
!!     tbmd_energy_type
!!
!!   FUNCTION
!!     The TBMD energy type
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - id_nr: unique identifier
!!     - ref_count: reference counter (see doc/ReferenceCounting.html)
!!     - E_kinetic: kinetic energy
!!     - E_bandstructure: band-structure energy
!!     - E_repulsion: pair-repulsion energy
!!     - E_correction: (bond) correction energy
!!     - E_tb: total tight-binding energy
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!   SOURCE
!****************************************************************************
  TYPE tbmd_energy_type
     INTEGER                                         :: id_nr, &
                                                        ref_count
     REAL(KIND=dp)                                   :: E_kinetic, &
                                                        E_bandstructure, &
                                                        E_repulsion, &
                                                        E_correction, &
                                                        E_tb

  END TYPE tbmd_energy_type
!!***
!****************************************************************************

!!****s* tbmd_energy_types/tbmd_energy_p_type [1.0] *
!!
!!   NAME
!!     tbmd_energy_p_type
!!
!!   FUNCTION
!!     Type which allows the use of arrays of pointers
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - tbmd_energy_ptr: Pointer to the actual tbmd energy type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!   SOURCE
!***************************************************************************
  TYPE tbmd_energy_p_type
     TYPE(tbmd_energy_type), POINTER :: tbmd_energy_ptr
  END TYPE tbmd_energy_p_type
!!***
!****************************************************************************

CONTAINS

!!****f* tbmd_energy_types/tbmd_energy_retain [1.0] *
!!
!!   NAME
!!     tbmd_energy_retain
!!
!!   FUNCTION
!!     Retains a tbmd force type (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - tbmd_energy: The tbmd energy type to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE tbmd_energy_retain(tbmd_energy, error)

    TYPE(tbmd_energy_type), POINTER          :: tbmd_energy
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_energy_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(tbmd_energy), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
      CPPrecondition(tbmd_energy%ref_count>0, cp_failure_level, routineP, error, failure)
      tbmd_energy%ref_count = tbmd_energy%ref_count+1
    END IF
  END SUBROUTINE tbmd_energy_retain
!***************************************************************************

!!****f* tbmd_energy_types/tbmd_energy_release [1.0] *
!!
!!   NAME
!!     tbmd_energy_release
!!
!!   FUNCTION
!!     Releases the given tbmd energy type (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - tbmd_energy: The tbmd energy type to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE tbmd_energy_release(tbmd_energy, error)

    TYPE(tbmd_energy_type), POINTER          :: tbmd_energy
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_energy_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    IF (ASSOCIATED(tbmd_energy)) THEN
      CPPrecondition(tbmd_energy%ref_count>0, cp_failure_level, routineP, error, failure)
      tbmd_energy%ref_count = tbmd_energy%ref_count-1
      IF (tbmd_energy%ref_count<1) THEN
        DEALLOCATE(tbmd_energy, stat=stat)
        CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)
      END IF
    END IF
    NULLIFY(tbmd_energy)
  END SUBROUTINE tbmd_energy_release
!***************************************************************************

!!****f* tbmd_energy_types/tbmd_energy_get [1.0] *
!!
!!   NAME
!!     tbmd_energy_get
!!
!!   FUNCTION
!!     Returns various attributes of the tbmd type
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - tbmd_energy: The enquired tbmd energy type
!!     - id_nr: unique identifier
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!     - E_kinetic: kinetic energy
!!     - E_bandstructure: band-structure energy
!!     - E_repulsion: pair-repulsion energy
!!     - E_correction: (bond) correction energy
!!     - E_tb: total tight-binding energy
!!
!!     For possible missing arguments see the attributes of tbmd_energy_type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE tbmd_energy_get(tbmd_energy, id_nr, E_kinetic, &
                             E_bandstructure, E_repulsion, &
                             E_correction, E_tb, error)

    TYPE(tbmd_energy_type), POINTER          :: tbmd_energy
    INTEGER, INTENT(OUT), OPTIONAL           :: id_nr
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: E_kinetic, E_bandstructure, &
                                                E_repulsion, E_correction, &
                                                E_tb
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_energy_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(tbmd_energy), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
      CPPrecondition(tbmd_energy%ref_count>0, cp_failure_level, routineP, error, failure)
    END IF

    IF (.NOT. failure) THEN
      IF (PRESENT(id_nr)) id_nr = tbmd_energy%id_nr
      IF (PRESENT(E_kinetic)) E_kinetic = tbmd_energy%E_kinetic
      IF (PRESENT(E_bandstructure)) E_bandstructure = tbmd_energy%E_bandstructure
      IF (PRESENT(E_repulsion)) E_repulsion = tbmd_energy%E_repulsion
      IF (PRESENT(E_correction)) E_correction = tbmd_energy%E_correction
      IF (PRESENT(E_tb)) E_tb = tbmd_energy%E_tb
    END IF
  END SUBROUTINE tbmd_energy_get
!***************************************************************************

!!****f* tbmd_energy_types/tbmd_energy_set [1.0] *
!!
!!   NAME
!!     tbmd_energy_set
!!
!!   FUNCTION
!!     Sets various attributes of the tbmd energy type
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - tbmd_energy: The enquired tbmd energy type
!!     - id_nr: unique identifier
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!     - E_kinetic: kinetic energy
!!     - E_bandstructure: band-structure energy
!!     - E_repulsion: pair-repulsion energy
!!     - E_correction: (bond) correction energy
!!     - E_tb: total tight-binding energy
!!
!!     For possible missing arguments see the attributes of tbmd_energy_type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE tbmd_energy_set(tbmd_energy, id_nr, E_kinetic, &
                             E_bandstructure, E_repulsion, &
                             E_correction, E_tb, error)

    TYPE(tbmd_energy_type), POINTER          :: tbmd_energy
    INTEGER, INTENT(IN), OPTIONAL            :: id_nr
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: E_kinetic, E_bandstructure, &
                                                E_repulsion, E_correction, &
                                                E_tb
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_energy_set', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(tbmd_energy), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
      CPPrecondition(tbmd_energy%ref_count>0, cp_failure_level, routineP, error, failure)
    END IF

    IF (.NOT. failure) THEN
      IF (PRESENT(id_nr)) tbmd_energy%id_nr = id_nr
      IF (PRESENT(E_kinetic)) tbmd_energy%E_kinetic = E_kinetic
      IF (PRESENT(E_bandstructure)) tbmd_energy%E_bandstructure = E_bandstructure
      IF (PRESENT(E_repulsion)) tbmd_energy%E_repulsion = E_repulsion
      IF (PRESENT(E_correction)) tbmd_energy%E_correction = E_correction
      IF (PRESENT(E_tb)) tbmd_energy%E_tb = E_tb
    END IF
  END SUBROUTINE tbmd_energy_set
!***************************************************************************

!!****f* tbmd_energy_types/tbmd_energy_clear [1.0] *
!!
!!   NAME
!!     tbmd_energy_clear
!!
!!   FUNCTION
!!     Reinitializes the tbmd energy type
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - tbmd_energy: The tbmd energy type to be reinitialized
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE tbmd_energy_clear(tbmd_energy, error)

    TYPE(tbmd_energy_type), POINTER          :: tbmd_energy
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_energy_clear', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(tbmd_energy), cp_failure_level, routineP, error, failure)
 
    IF (.NOT. failure) THEN
      tbmd_energy%E_kinetic       = 0.0_dp
      tbmd_energy%E_bandstructure = 0.0_dp
      tbmd_energy%E_repulsion     = 0.0_dp
      tbmd_energy%E_correction    = 0.0_dp
      tbmd_energy%E_tb            = 0.0_dp
    ELSE
      CALL tbmd_energy_create(tbmd_energy, error=error)
    END IF
  END SUBROUTINE tbmd_energy_clear
!***************************************************************************

!!****f* tbmd_energy_types/tbmd_energy_create [1.0] *
!!
!!   NAME
!!     tbmd_energy_create
!!
!!   FUNCTION
!!     Creates a tbmd energy type instance
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - tbmd_energy: The tbmd energy type to be created
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE tbmd_energy_create(tbmd_energy, error)

    TYPE(tbmd_energy_type), POINTER          :: tbmd_energy
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_energy_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    ALLOCATE(tbmd_energy, stat=stat)
    CPPostcondition(stat == 0, cp_failure_level, routineP, error, failure)
   
    IF (.NOT. failure) THEN
      tbmd_energy%ref_count = 1
      last_tbmd_id=last_tbmd_id+1
      tbmd_energy%id_nr = last_tbmd_id
      CALL tbmd_energy_clear(tbmd_energy, error=error)
    END IF
  END SUBROUTINE tbmd_energy_create
!***************************************************************************

END MODULE tbmd_energy_types
