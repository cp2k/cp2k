!-----------------------------------------------------------------------------! 
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2005  CP2K developers group                            !
!-----------------------------------------------------------------------------!


!!****** cp2k/harris_env_types [1.0] *
!!
!!   NAME
!!     harris_energy_types
!!
!!   FUNCTION
!!     The Harris environment, which contains beside the harris energy and force 
!!     types all the informations needed to perform and analyse the Harris 
!!     functional energy and force corrections.
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     07.2005 initial create [tdk]
!!
!!   SOURCE
!******************************************************************************
MODULE harris_env_types

  USE harris_energy_types,             ONLY: harris_energy_clear,&
                                             harris_energy_release,&
                                             harris_energy_type
!harris_energy_create
  USE harris_force_types,              ONLY: harris_force_clear,&
                                             harris_force_release,&
                                             harris_force_type
!harris_force_create
  USE qs_rho_types,                    ONLY: qs_rho_release,&
                                             qs_rho_type
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'harris_env_types'
  INTEGER, PRIVATE, SAVE :: last_harris_id=0

  ! *** Public data types ***
  PUBLIC :: harris_env_type, &
            harris_env_p_type

  ! *** Public subroutines ***
  PUBLIC :: harris_env_retain, &
            harris_env_release, &
            harris_env_set,&
            harris_env_get, &
            harris_env_clear, &
            harris_env_create
!!***
!****************************************************************************

!!****s* harris_env_types/harris_env_type [1.0] *
!!
!!   NAME
!!     harris_env_type
!!
!!   FUNCTION
!!     The Harris Environment
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - id_nr: unique identifier
!!     - ref_count: reference counter (see doc/ReferenceCounting.html)
!!     - harris_energy: The harris energy type
!!     - harris_force: The harris force type
!!     - rho: The rho data structure in which diverse representations of 
!!            the density are stored
!!     - rho_diff: rho_out - rho_in
!!     - natom: The number of involved atoms
!!     - nspins: The number of spins
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     07.2005 initial create [tdk]
!!
!!   SOURCE
!****************************************************************************
  TYPE harris_env_type
    INTEGER                                          :: id_nr, &
                                                        ref_count
    TYPE(harris_energy_type), POINTER                :: harris_energy
    TYPE(harris_force_type), POINTER                 :: harris_force
    TYPE(qs_rho_type), POINTER                       :: rho, &
                                                        rho_diff
    INTEGER                                          :: natom, &
                                                        nspins

  END TYPE harris_env_type
!!***
!****************************************************************************

!!****s* harris_env_types/harris_env_p_type [1.0] *
!!
!!   NAME
!!     harris_env_p_type
!!
!!   FUNCTION
!!     Type which allows the use of arrays of pointers
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - harris_env_ptr: Pointer to the actual harris environment
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     07.2005 initial create [tdk]
!!
!!   SOURCE
!***************************************************************************
  TYPE harris_env_p_type
    TYPE(harris_env_type), POINTER :: harris_env_ptr
  END TYPE harris_env_p_type
!!***
!****************************************************************************

CONTAINS

!!****f* harris_env_types/harris_env_retain [1.0] *
!!
!!   NAME
!!     harris_env_retain
!!
!!   SYNOPSIS
!!     SUBROUTINE harris_env_retain(harris_env, error)
!!       TYPE(harris_env_type), POINTER :: harris_env
!!       TYPE(cp_error_type), OPTIONAL, INTENT(INOUT):: error
!!     END SUBROUTINE harris_env_retain
!!
!!   FUNCTION
!!     Retains an harris environment (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - harris_env: The harris environment to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!     For possible missing arguments see the attributes of harris_env_type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     07.2005 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE harris_env_retain(harris_env, error)

    TYPE(harris_env_type), POINTER           :: harris_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(harris_env), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
      CPPrecondition(harris_env%ref_count>0, cp_failure_level, routineP, error, failure)
      harris_env%ref_count = harris_env%ref_count+1
    END IF
  END SUBROUTINE harris_env_retain
!***************************************************************************

!!****f* harris_env_types/harris_env_release [1.0] *
!!
!!   NAME
!!     harris_env_release
!!
!!   FUNCTION
!!     Releases the given harris environment (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - harris_env: The harris environment to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!     For possible missing arguments see the attributes of harris_env_type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     07.2005 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE harris_env_release(harris_env, error)

    TYPE(harris_env_type), POINTER           :: harris_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    IF (ASSOCIATED(harris_env)) THEN
      CPPrecondition(harris_env%ref_count>0, cp_failure_level, routineP, error, failure)
      harris_env%ref_count = harris_env%ref_count-1
      IF (harris_env%ref_count<1) THEN
        IF (ASSOCIATED(harris_env%harris_energy)) THEN
          CALL harris_energy_release(harris_env%harris_energy, error=error)
        END IF
        IF (ASSOCIATED(harris_env%harris_force)) THEN
          CALL harris_force_release(harris_env%harris_force, error=error)
        END IF
        IF (ASSOCIATED(harris_env%rho)) THEN
          CALL qs_rho_release(harris_env%rho, error=error)
        END IF
        IF (ASSOCIATED(harris_env%rho_diff)) THEN
          CALL qs_rho_release(harris_env%rho_diff, error=error)
        END IF
        DEALLOCATE(harris_env, stat=stat)
        CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)
      END IF
    END IF
    NULLIFY(harris_env)
  END SUBROUTINE harris_env_release
!***************************************************************************

!!****f* harris_env_types/harris_env_get [1.0] *
!!
!!   NAME
!!     harris_env_get
!!
!!   FUNCTION
!!     Returns various attributes of the harris environment
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - harris_env: The enquired harris environment
!!     - id_nr: unique identifier
!!     - harris_energy: The harris energy type
!!     - harris_force: The harris force type
!!       see module cp_error_handling
!!     - rho: The rho data structure in which diverse representations of 
!!            the density are stored
!!     - rho_diff: rho_out - rho_in
!!     - error: variable to control error logging, stopping,... 
!!
!!     For possible missing arguments see the attributes of harris_env_type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     07.2005 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE harris_env_get(harris_env, id_nr, harris_energy, harris_force, &
                            rho, rho_diff, error)

    TYPE(harris_env_type), POINTER           :: harris_env
    INTEGER, INTENT(OUT), OPTIONAL           :: id_nr
    TYPE(harris_energy_type), OPTIONAL, &
      POINTER                                :: harris_energy
    TYPE(harris_force_type), OPTIONAL, &
      POINTER                                :: harris_force
    TYPE(qs_rho_type), OPTIONAL, POINTER     :: rho, rho_diff
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_env_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(harris_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(harris_env%ref_count>0, cp_failure_level, routineP, error, failure)

    IF (.NOT. failure) THEN
      IF (PRESENT(id_nr)) id_nr = harris_env%id_nr
      IF (PRESENT(harris_energy)) harris_energy => harris_env%harris_energy
      IF (PRESENT(harris_force)) harris_force => harris_env%harris_force
      IF (PRESENT(rho)) rho => harris_env%rho
      IF (PRESENT(rho_diff)) rho_diff => harris_env%rho_diff
    END IF

  END SUBROUTINE harris_env_get
!***************************************************************************

!!****f* harris_env_types/harris_env_set [1.0] *
!!
!!   NAME
!!     harris_env_set
!!
!!   FUNCTION
!!     Sets various attributes of the harris environment
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - harris_env: The enquired harris environment
!!     - id_nr: unique identifier
!!     - harris_energy: The harris energy type
!!     - harris_force: The harris force type
!!     - rho: The rho data structure in which diverse representations of 
!!            the density are stored
!!     - rho_diff: rho_out - rho_in
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!     For possible missing arguments see the attributes of harris_env_type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     07.2005 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE harris_env_set(harris_env, id_nr, harris_energy, harris_force, &
                            rho, rho_diff, error)

    TYPE(harris_env_type), POINTER           :: harris_env
    INTEGER, INTENT(IN), OPTIONAL            :: id_nr
    TYPE(harris_energy_type), OPTIONAL, &
      POINTER                                :: harris_energy
    TYPE(harris_force_type), OPTIONAL, &
      POINTER                                :: harris_force
    TYPE(qs_rho_type), OPTIONAL, POINTER     :: rho, rho_diff
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_env_set', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(harris_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(harris_env%ref_count>0, cp_failure_level, routineP, error, failure)

    IF (.NOT. failure) THEN
      IF (PRESENT(id_nr)) harris_env%id_nr = id_nr
      IF (PRESENT(harris_energy)) harris_env%harris_energy => harris_energy
      IF (PRESENT(harris_force)) harris_env%harris_force => harris_force
      IF (PRESENT(rho)) harris_env%rho => rho
      IF (PRESENT(rho_diff)) harris_env%rho_diff => rho_diff
    END IF
  END SUBROUTINE harris_env_set
!***************************************************************************

!!****f* harris_env_types/harris_env_clear [1.0] *
!!
!!   NAME
!!     harris_env_clear
!!
!!   FUNCTION
!!     Reinitializes the harris environment
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - harris_env: The harris energy environment to be reinitialized
!!     - natom: Number of atoms
!!     - nspins: Number of spins
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!     For possible missing arguments see the attributes of harris_env_type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     07.2005 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE harris_env_clear(harris_env, natom, nspins, error)

    TYPE(harris_env_type), POINTER           :: harris_env
    INTEGER, INTENT(IN)                      :: natom, nspins
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_env_clear', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(harris_env), cp_failure_level, routineP, error, failure)
 
    IF (.NOT. failure) THEN
      harris_env%natom  = natom
      harris_env%nspins = nspins

      CALL harris_energy_clear(harris_env%harris_energy, error=error)
      CALL harris_force_clear(harris_env%harris_force, natom, error=error)

      !CALL coeff_zero(harris_env%rho_core)
      
      IF (ASSOCIATED(harris_env%rho)) THEN
        CALL qs_rho_release(harris_env%rho, error=error)
      END IF
      IF (ASSOCIATED(harris_env%rho_diff)) THEN
        CALL qs_rho_release(harris_env%rho_diff, error=error)
      END IF
    ELSE
      CALL harris_env_create(harris_env, natom, nspins, error=error)
    END IF
  END SUBROUTINE harris_env_clear
!***************************************************************************

!!****f* harris_env_types/harris_env_create [1.0] *
!!
!!   NAME
!!     harris_env_create
!!
!!   FUNCTION
!!     Creates an harris environment
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - harris_env: The harris environment to be created
!!     - natom: Number of atoms
!!     - nspins: Number of spins
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!     For possible missing arguments see the attributes of harris_env_type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     07.2005 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE harris_env_create(harris_env, natom, nspins, error)

    TYPE(harris_env_type), POINTER           :: harris_env
    INTEGER, INTENT(IN)                      :: natom, nspins
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    ALLOCATE(harris_env, stat=stat)
    CPPostcondition(stat == 0, cp_failure_level, routineP, error, failure)

    !IF (.NOT. failure) THEN
    !  ALLOCATE(harris_env%harris_energy, stat=stat)
    !  CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    !END IF

    !IF (.NOT. failure) THEN
    !  ALLOCATE(harris_env%harris_force, stat=stat)
    !  CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    !END IF

    !IF (.NOT. failure) THEN
    !  ALLOCATE(harris_env%rho, stat=stat)
    !  CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    !END IF
   
    IF (.NOT. failure) THEN
      harris_env%ref_count = 1
      last_harris_id=last_harris_id+1
      harris_env%id_nr = last_harris_id
      !CALL harris_env_clear(harris_env, natom=natom, nspins=nspins, error=error)
      NULLIFY(harris_env%rho)
      NULLIFY(harris_env%rho_diff)
      harris_env%natom = natom
      harris_env%nspins = nspins

    END IF
  END SUBROUTINE harris_env_create
!***************************************************************************

END MODULE harris_env_types
