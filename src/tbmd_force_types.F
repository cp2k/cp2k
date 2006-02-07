!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2006  CP2K developers group                            !
!-----------------------------------------------------------------------------!

!!****h* cp2k/tbmd_force_types [1.0] *
!!
!!   NAME
!!     tbmd_force_types
!!
!!   FUNCTION
!!     The tight-binding force type.
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
MODULE tbmd_force_types
  USE kinds,                           ONLY: dp
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tbmd_force_types'
  INTEGER, PRIVATE, SAVE :: last_tbmd_id=0

  ! *** Public data types ***
  PUBLIC :: tbmd_force_type, &
            tbmd_force_p_type

  ! *** Public subroutines ***
  PUBLIC :: tbmd_force_retain, &
            tbmd_force_release, &
            tbmd_force_set,&
            tbmd_force_get, &
            tbmd_force_clear, &
            tbmd_force_create
!!***
!****************************************************************************

!!****s* tbmd_force_types/tbmd_force_type [1.0] *
!!
!!   NAME
!!     tbmd_force_type
!!
!!   FUNCTION
!!     The TBMD force type
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - id_nr: unique identifier
!!     - ref_count: reference counter (see doc/ReferenceCounting.html)
!!     - f_tb: The total tight-binding forces
!!     - f_kinetic: Kinetic forces
!!     - f_bandstruct: Forces of the band-structure term
!!     - f_repulsion: Forces due to the repulsive pair-potential
!!     - f_total: The sum of the total forces of all involved atoms
!!         In principle this sum should be zero, due to translation invariance
!!         It can therefore be used as a criterion for the quality of the forces
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!   SOURCE
!****************************************************************************
  TYPE tbmd_force_type
    INTEGER                                          :: id_nr, &
                                                        ref_count
    REAL(KIND=dp), DIMENSION(:,:), POINTER           :: f_tb, &
                                                        f_kinetic, &
                                                        f_bandstruct, &
                                                        f_repulsion
    REAL(KIND=dp), DIMENSION(:), POINTER             :: f_total

  END TYPE tbmd_force_type
!!***
!****************************************************************************

!!****s* tbmd_force_types/tbmd_force_p_type [1.0] *
!!
!!   NAME
!!     tbmd_force_p_type
!!
!!   FUNCTION
!!     Type which allows the use of arrays of pointers
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - tbmd_force_ptr: Pointer to the actual tbmd force type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!   SOURCE
!***************************************************************************
  TYPE tbmd_force_p_type
     TYPE(tbmd_force_type), POINTER :: tbmd_force_ptr
  END TYPE tbmd_force_p_type
!!***
!****************************************************************************

CONTAINS

!!****f* tbmd_force_types/tbmd_force_retain [1.0] *
!!
!!   NAME
!!     tbmd_force_retain
!!
!!   FUNCTION
!!     Retains a tbmd force type (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - tbmd_force: The tbmd force type to retain
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
  SUBROUTINE tbmd_force_retain(tbmd_force, error)

    TYPE(tbmd_force_type), POINTER           :: tbmd_force
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_force_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(tbmd_force), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
      CPPrecondition(tbmd_force%ref_count>0, cp_failure_level, routineP, error, failure)
      tbmd_force%ref_count = tbmd_force%ref_count+1
    END IF
  END SUBROUTINE tbmd_force_retain
!***************************************************************************

!!****f* tbmd_force_types/tbmd_force_release [1.0] *
!!
!!   NAME
!!     tbmd_force_release
!!
!!   FUNCTION
!!     Releases the given tbmd force type (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - tbmd_force: The tbmd force type to release
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
  SUBROUTINE tbmd_force_release(tbmd_force, error)

    TYPE(tbmd_force_type), POINTER           :: tbmd_force
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_force_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.

    
    IF (ASSOCIATED(tbmd_force)) THEN
      CPPrecondition(tbmd_force%ref_count>0, cp_failure_level, routineP, error, failure)
      tbmd_force%ref_count = tbmd_force%ref_count-1
      IF (tbmd_force%ref_count<1) THEN
        IF (ASSOCIATED(tbmd_force%f_tb)) THEN
          DEALLOCATE(tbmd_force%f_tb, stat=stat)
          CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)
        END IF
        IF (ASSOCIATED(tbmd_force%f_kinetic)) THEN
          DEALLOCATE(tbmd_force%f_kinetic, stat=stat)
          CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)
        END IF
        IF (ASSOCIATED(tbmd_force%f_bandstruct)) THEN
          DEALLOCATE(tbmd_force%f_bandstruct, stat=stat)
          CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)
        END IF
        IF (ASSOCIATED(tbmd_force%f_repulsion)) THEN
          DEALLOCATE(tbmd_force%f_repulsion, stat=stat)
          CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)
        END IF
        IF (ASSOCIATED(tbmd_force%f_total)) THEN
          DEALLOCATE(tbmd_force%f_total, stat=stat)
          CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)
        END IF
        DEALLOCATE(tbmd_force, stat=stat)
        CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)
      END IF
    END IF
    NULLIFY(tbmd_force)

  END SUBROUTINE tbmd_force_release
!***************************************************************************

!!****f* tbmd_force_types/tbmd_force_get [1.0] *
!!
!!   NAME
!!     tbmd_force_get
!!
!!   FUNCTION
!!     Returns various attributes of the tbmd force type
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - tbmd_force: The enquired tbmd force type
!!     - id_nr: unique identifier
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!     - f_kinetic: The forces due to the kinetic energy
!!     - f_bandstruct: Forces of the band-structure term
!!     - f_repulsion: Forces due to the repulsive pair-potential
!!     - f_tb: The total tight-binding forces
!!     - f_total: The sum of the total forces of all involved atoms
!!         In principle this sum should be zero, due to translation invariance
!!         It can therefore be used as a criterion for the quality of the forces
!!
!!     For possible missing arguments see the attributes of tbmd_force_type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE tbmd_force_get(tbmd_force, id_nr, f_kinetic, f_bandstruct, &
                            f_repulsion, f_tb, f_total, error)

    TYPE(tbmd_force_type), POINTER           :: tbmd_force
    INTEGER, INTENT(OUT), OPTIONAL           :: id_nr
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: f_kinetic, f_bandstruct, &
                                                f_repulsion, f_tb
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: f_total
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_force_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(tbmd_force), cp_failure_level, routineP, error, failure)
    IF (.NOT.failure) THEN
      CPPrecondition(tbmd_force%ref_count>0, cp_failure_level, routineP, error, failure)
    END IF

    IF (.NOT. failure) THEN
      IF (PRESENT(id_nr)) id_nr = tbmd_force%id_nr
      IF (PRESENT(f_kinetic)) f_kinetic = tbmd_force%f_kinetic
      IF (PRESENT(f_bandstruct)) f_bandstruct = tbmd_force%f_bandstruct
      IF (PRESENT(f_repulsion)) f_repulsion = tbmd_force%f_repulsion
      IF (PRESENT(f_tb)) f_tb = tbmd_force%f_tb
      IF (PRESENT(f_total)) f_total = tbmd_force%f_total
    END IF
  END SUBROUTINE tbmd_force_get
!***************************************************************************

!!****f* tbmd_force_types/tbmd_force_set [1.0] *
!!
!!   NAME
!!     tbmd_force_set
!!
!!   FUNCTION
!!     Sets various attributes of the tbmd force type
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - tbmd_force: The enquired tbmd force type
!!     - id_nr: unique identifier
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!     - f_kinetic: The forces due to the kinetic energy
!!     - f_bandstruct: Forces of the band-structure term
!!     - f_repulsion: Forces due to the repulsive pair-potential
!!     - f_tb: The total tight-binding forces
!!     - f_total: The sum of the total forces of all involved atoms
!!         In principle this sum should be zero, due to translation invariance
!!         It can therefore be used as a criterion for the quality of the forces
!!
!!     For possible missing arguments see the attributes of tbmd_force_type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE tbmd_force_set(tbmd_force, id_nr, f_kinetic, f_bandstruct, &
                            f_repulsion, f_tb, f_total, error)

    TYPE(tbmd_force_type), POINTER           :: tbmd_force
    INTEGER, INTENT(IN), OPTIONAL            :: id_nr
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: f_kinetic, f_bandstruct, &
                                                f_repulsion, f_tb
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: f_total
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_force_set', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(tbmd_force), cp_failure_level, routineP, error, failure)
    IF (.NOT.failure) THEN
      CPPrecondition(tbmd_force%ref_count>0, cp_failure_level, routineP, error, failure)
    END IF

    IF (.NOT. failure) THEN
      IF (PRESENT(id_nr)) tbmd_force%id_nr = id_nr
      IF (PRESENT(f_kinetic)) tbmd_force%f_kinetic = f_kinetic
      IF (PRESENT(f_bandstruct)) tbmd_force%f_bandstruct = f_bandstruct
      IF (PRESENT(f_repulsion)) tbmd_force%f_repulsion = f_repulsion
      IF (PRESENT(f_tb)) tbmd_force%f_tb = f_tb
      IF (PRESENT(f_total)) tbmd_force%f_total = f_total
    END IF
  END SUBROUTINE tbmd_force_set
!***************************************************************************

!!****f* tbmd_force_types/tbmd_force_clear [1.0] *
!!
!!   NAME
!!     tbmd_force_clear
!!
!!   FUNCTION
!!     Reinitializes the tbmd force type
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - tbmd_force: The tbmd force type to be cleared
!!     - natom: The number of atoms
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
  SUBROUTINE tbmd_force_clear(tbmd_force, natom, error)

    TYPE(tbmd_force_type), POINTER           :: tbmd_force
    INTEGER                                  :: natom
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_force_clear', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j
    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(tbmd_force), cp_failure_level, routineP, error, failure)

    IF (.NOT. failure) THEN
      DO i = 1, natom
        DO j = 1, 3
          tbmd_force%f_kinetic(i,j)    = 0.0_dp
          tbmd_force%f_bandstruct(i,j) = 0.0_dp
          tbmd_force%f_repulsion(i,j)  = 0.0_dp
          tbmd_force%f_tb(i,j)         = 0.0_dp
        END DO
      END DO
      tbmd_force%f_total(:) = 0.0_dp
    ELSE
      CALL tbmd_force_create(tbmd_force, natom, error=error)
    END IF
  END SUBROUTINE tbmd_force_clear
!***************************************************************************

!!****f* tbmd_force_types/tbmd_force_create [1.0] *
!!
!!   NAME
!!     tbmd_force_create
!!
!!   FUNCTION
!!     Creates a tbmd force type instance
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - tbmd_force: The tbmd force type to be created
!!     - natom: The number of atoms
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
  SUBROUTINE tbmd_force_create(tbmd_force, natom, error)

    TYPE(tbmd_force_type), POINTER           :: tbmd_force
    INTEGER                                  :: natom
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_force_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.

    ALLOCATE(tbmd_force, stat=stat)
    CPPostcondition(stat == 0, cp_failure_level, routineP, error, failure)
    ALLOCATE(tbmd_force%f_kinetic(natom,3), stat=stat)
    CPPostcondition(stat == 0, cp_failure_level, routineP, error, failure)

    ALLOCATE(tbmd_force%f_bandstruct(natom,3), stat=stat)
    CPPostcondition(stat == 0, cp_failure_level, routineP, error, failure)
    ALLOCATE(tbmd_force%f_repulsion(natom,3), stat=stat)
    CPPostcondition(stat == 0, cp_failure_level, routineP, error, failure)
    ALLOCATE(tbmd_force%f_tb(natom,3), stat=stat)
    CPPostcondition(stat == 0, cp_failure_level, routineP, error, failure)

    ALLOCATE(tbmd_force%f_total(3), stat=stat)
    CPPostcondition(stat == 0, cp_failure_level, routineP, error, failure)

    IF (.NOT. failure) THEN
      tbmd_force%ref_count = 1
      last_tbmd_id=last_tbmd_id+1
      tbmd_force%id_nr = last_tbmd_id
      CALL tbmd_force_clear(tbmd_force, natom, error=error)
    END IF
  END SUBROUTINE tbmd_force_create
!***************************************************************************

END MODULE tbmd_force_types
