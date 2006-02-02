!-----------------------------------------------------------------------------!  
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2005  CP2K developers group                            !
!-----------------------------------------------------------------------------!


!!****** cp2k/harris_energy_types [1.0] *
!!
!!   NAME
!!     harris_energy_types
!!
!!   FUNCTION
!!     The Harris energy type, which contains all the informations 
!!     needed to perform and analyse the Harris functional energy correction.
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     06.2005 initial create [tdk]
!!
!!   SOURCE
!******************************************************************************
MODULE harris_energy_types

  USE kinds,                           ONLY: dp
#include "cp_common_uses.h"
 
  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='harris_energy_types'
  INTEGER, PRIVATE, save :: last_harris_id=0

  ! *** Public data types ***
  PUBLIC :: harris_energy_type, &
            harris_energy_p_type

  ! *** Public subroutines ***
  PUBLIC :: harris_energy_retain, &
            harris_energy_release, &
            harris_energy_set,&
            harris_energy_get, &
            harris_energy_clear, &
            harris_energy_create
!!***
!****************************************************************************

!!****s* harris_energy_types/harris_energy_type [1.0] *
!!
!!   NAME
!!     harris_energy_type
!!
!!   FUNCTION
!!     The Harris energy type
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - id_nr: unique identifier
!!     - ref_count: reference counter (see doc/ReferenceCounting.html)
!!     - Ehartree_elec: The electronic hartree energy
!!     - Ehartree_core: The hartree energy computed with rho_core
!!     - Exc: The exchange-correlation energy
!!     - integral_vxc: The xc potential integral
!!     - Ecore: The core hamiltonian energy
!!     - Ecore_overlap: The overalp energy
!!     - Ecore_self: The self-interaction energy
!!     - EII: The ion-ion Energy
!!     - Eharris_correction: The energy correction due to the harris functional
!!     - sum_of_eigenvalues: The sum of the eigenvalues
!!     - Eharris: The actual harris functional energy
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     06.2005 initial create [tdk]
!!
!!   SOURCE
!****************************************************************************
  TYPE harris_energy_type
     INTEGER                                         :: id_nr, &
                                                        ref_count
     REAL(KIND=dp)                                   :: Ehartree_elec, &
                                                        Ehartree_core, &
                                                        !Ehartree_tot, &
                                                        !Ehartree_elec_core, &
                                                        Exc, &
                                                        integral_vxc, &
                                                        Ecore, &
                                                        Ecore_overlap, &
                                                        Ecore_self, &
                                                        EII, &
                                                        Eharris_correction, &
                                                        sum_of_eigenvalues, &
                                                        Eharris

  END TYPE harris_energy_type
!!***
!****************************************************************************

!!****s* harris_energy_types/harris_energy_p_type [1.0] *
!!
!!   NAME
!!     harris_energy_p_type
!!
!!   FUNCTION
!!     Type which allows the use of arrays of pointers
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - harris_energy_ptr: Pointer to the actual harris energy type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     06.2005 initial create [tdk]
!!
!!   SOURCE
!***************************************************************************
  TYPE harris_energy_p_type
     TYPE(harris_energy_type), POINTER :: harris_energy_ptr
  END TYPE harris_energy_p_type
!!***
!****************************************************************************

CONTAINS

!!****f* harris_energy_types/harris_energy_retain [1.0] *
!!
!!   NAME
!!     harris_energy_retain
!!
!!   SYNOPSIS
!!     SUBROUTINE harris_energy_retain(harris_energy, error)
!!       TYPE(harris_energy_type), POINTER :: harris_energy
!!       TYPE(cp_error_type), OPTIONAL, INTENT(INOUT):: error
!!     END SUBROUTINE harris_energy_retain
!!
!!   FUNCTION
!!     Retains a harris force type (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - harris_energy: The harris energy type to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!     For possible missing arguments see the attributes of harris_energy_type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     06.2005 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE harris_energy_retain(harris_energy, error)

    TYPE(harris_energy_type), POINTER        :: harris_energy
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_energy_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------
   
    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(harris_energy), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
      CPPrecondition(harris_energy%ref_count>0, cp_failure_level, routineP, error, failure)
      harris_energy%ref_count = harris_energy%ref_count+1
    END IF
  END SUBROUTINE harris_energy_retain
!***************************************************************************

!!****f* harris_energy_types/harris_energy_release [1.0] *
!!
!!   NAME
!!     harris_energy_release
!!
!!   SYNOPSIS
!!     SUBROUTINE harris_energy_release(harris_energy, error)
!!       TYPE(harris_energy_type), POINTER :: harris_energy
!!       TYPE(cp_error_type), OPTIONAL, INTENT(INOUT) :: error
!!     END SUBROUTINE harris_energy_release
!!
!!   FUNCTION
!!     Releases the given harris energy type (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - harris_energy: The harris energy type to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!     For possible missing arguments see the attributes of harris_energy_type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     06.2005 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE harris_energy_release(harris_energy, error)

    TYPE(harris_energy_type), POINTER        :: harris_energy
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_energy_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    IF (ASSOCIATED(harris_energy)) THEN
      CPPrecondition(harris_energy%ref_count>0, cp_failure_level, routineP, error, failure)
      harris_energy%ref_count = harris_energy%ref_count-1
      IF (harris_energy%ref_count<1) THEN
        DEALLOCATE(harris_energy, stat=stat)
        CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)
      END IF
    END IF
    NULLIFY(harris_energy)
  END SUBROUTINE harris_energy_release
!***************************************************************************

!!****f* harris_energy_types/harris_energy_get [1.0] *
!!
!!   NAME
!!     harris_energy_get
!!
!!   SYNOPSIS
!!     SUBROUTINE harris_energy_get(harris_energy, id_nr, Ehartree_elec, 
!!                                  Ehartree_core, Exc, integral_vxc, Ecore, 
!!                                  Ecore_overlap, Ecore_self, EII, 
!!                                  Eharris_correction, sum_of_eigenvalues, 
!!                                  Eharris, error)
!!       TYPE(harris_energy_type), POINTER    :: harris_energy
!!       INTEGER, INTENT(OUT), OPTIONAL       :: id_nr
!!       REAL(KIND=dp), INTENT(OUT), OPTIONAL :: Ehartreee_elec, 
!!                                               Ehartree_core, Exc,
!!                                               integral_vxc, Ecore,
!!                                               Ecore_overlap, Ecore_self
!!                                               EII, Eharris_correction,
!!                                               sum_of_eigenvalues, Eharris
!!       TYPE(cp_error_type), OPTIONAL, INTENT(INOUT) :: error
!!     END SUBROUTINE harris_energy_get
!!
!!   FUNCTION
!!     Returns various attributes of the harris type
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - harris_energy: The enquired harris energy type
!!     - id_nr: unique identifier
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!     - Ehartree_elec: The electronic hartree energy
!!     - Ehartree_core: The hartree energy computed with rho_core
!!     - Exc: The exchange-correlation energy
!!     - integral_vxc: The xc potential integral
!!     - Ecore: The core hamiltonian energy
!!     - Ecore_overlap: The overalp energy
!!     - Ecore_self: The self-interaction energy
!!     - EII: The ion-ion Energy
!!     - Eharris_correction: The energy correction due to the harris functional
!!     - sum_of_eigenvalues: The sum of the eigenvalues
!!     - Eharris: The actual harris functional energy
!!
!!     For possible missing arguments see the attributes of harris_energy_type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     06.2005 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE harris_energy_get(harris_energy, id_nr, Ehartree_elec, &
                               Ehartree_core, Exc, integral_vxc, Ecore, &
                               Ecore_overlap, Ecore_self, EII, &
                               Eharris_correction, sum_of_eigenvalues, &
                               Eharris, error)

    TYPE(harris_energy_type), POINTER          :: harris_energy
    INTEGER, INTENT(OUT), OPTIONAL             :: id_nr
    REAL(KIND=dp), INTENT(OUT), OPTIONAL       :: Ehartree_elec, &
                                                  Ehartree_core, &
                                                  Exc, &
                                                  integral_vxc, &
                                                  Ecore, &
                                                  Ecore_overlap, &
                                                  Ecore_self, &
                                                  EII, &
                                                  Eharris_correction, &
                                                  sum_of_eigenvalues, &
                                                  Eharris
    TYPE(cp_error_type), INTENT(INOUT), &
         OPTIONAL                              :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_energy_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                    :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(harris_energy), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
      CPPrecondition(harris_energy%ref_count>0, cp_failure_level, routineP, error, failure)
    END IF

    IF (.NOT. failure) THEN
      IF (PRESENT(id_nr)) id_nr = harris_energy%id_nr
      IF (PRESENT(Ehartree_elec)) Ehartree_elec = harris_energy%Ehartree_elec
      IF (PRESENT(Ehartree_core)) Ehartree_core = harris_energy%Ehartree_core
      IF (PRESENT(Exc)) Exc = harris_energy%Exc
      IF (PRESENT(integral_vxc)) integral_vxc = harris_energy%integral_vxc
      IF (PRESENT(Ecore)) Ecore = harris_energy%Ecore
      IF (PRESENT(Ecore_overlap)) Ecore_overlap = harris_energy%Ecore_overlap
      IF (PRESENT(Ecore_self)) Ecore_self = harris_energy%Ecore_self
      IF (PRESENT(EII)) EII = harris_energy%EII
      IF (PRESENT(Eharris_correction)) THEN
        Eharris_correction = harris_energy%Eharris_correction
      END IF
      IF (PRESENT(sum_of_eigenvalues)) THEN
        sum_of_eigenvalues = harris_energy%sum_of_eigenvalues
      END IF
      IF (PRESENT(Eharris)) Eharris = harris_energy%Eharris
    END IF
  END SUBROUTINE harris_energy_get
!***************************************************************************

!!****f* harris_energy_types/harris_energy_set [1.0] *
!!
!!   NAME
!!     harris_energy_set
!!
!!   SYNOPSIS
!!     SUBROUTINE harris_energy_set(harris_energy, id_nr, Ehartree_elec, 
!!                                  Ehartree_core, Exc, integral_vxc, Ecore, 
!!                                  Ecore_overlap, Ecore_self, EII, 
!!                                  Eharris_correction, sum_of_eigenvalues, 
!!                                  Eharris, error)
!!       TYPE(harris_energy_type), POINTER    :: harris_energy
!!       INTEGER, INTENT(OUT), OPTIONAL       :: id_nr
!!       REAL(KIND=dp), INTENT(OUT), OPTIONAL :: Ehartreee_elec, 
!!                                               Ehartree_core, Exc,
!!                                               integral_vxc, Ecore,
!!                                               Ecore_overlap, Ecore_self,
!!                                               EII, Eharris_correction,
!!                                               sum_of_eigenvalues, Eharris
!!       TYPE(cp_error_type), OPTIONAL, INTENT(INOUT) :: error
!!     END SUBROUTINE harris_energy_set
!!
!!   FUNCTION
!!     Sets various attributes of the harris energy type
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - harris_energy: The enquired harris energy type
!!     - id_nr: unique identifier
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!     - Ehartree_elec: The electronic hartree energy
!!     - Ehartree_core: The hartree energy computed with rho_core
!!     - Exc: The exchange-correlation energy
!!     - integral_vxc: The xc potential integral
!!     - Ecore: The core hamiltonian energy
!!     - Ecore_overlap: The overalp energy
!!     - Ecore_self: The self-interaction energy
!!     - EII: The ion-ion Energy
!!     - Eharris_correction: The energy correction due to the harris functional
!!     - sum_of_eigenvalues: The sum of the eigenvalues
!!     - Eharris: The actual harris functional energy
!!
!!     For possible missing arguments see the attributes of harris_energy_type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     06.2005 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE harris_energy_set(harris_energy, id_nr, Ehartree_elec, &
                               Ehartree_core, Exc, integral_vxc, Ecore, &
                               Ecore_overlap, Ecore_self, EII, &
                               Eharris_correction, sum_of_eigenvalues, &
                               Eharris, error)

    TYPE(harris_energy_type), POINTER          :: harris_energy
    INTEGER, INTENT(IN), OPTIONAL              :: id_nr
    REAL(KIND=dp), INTENT(IN), OPTIONAL        :: Ehartree_elec, &
                                                  Ehartree_core, &
                                                  Exc, &
                                                  integral_vxc, &
                                                  Ecore, &
                                                  Ecore_overlap, &
                                                  Ecore_self, &
                                                  EII, &
                                                  Eharris_correction, &
                                                  sum_of_eigenvalues, &
                                                  Eharris
    TYPE(cp_error_type), INTENT(INOUT), &
         OPTIONAL                              :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_energy_set', &
      routineP = moduleN//':'//routineN

    LOGICAL                                    :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(harris_energy), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
      CPPrecondition(harris_energy%ref_count>0, cp_failure_level, routineP, error, failure)
    END IF

    IF (.NOT. failure) THEN
      IF (PRESENT(id_nr)) harris_energy%id_nr = id_nr
      IF (PRESENT(Ehartree_elec)) harris_energy%Ehartree_elec = Ehartree_elec
      IF (PRESENT(Ehartree_core)) harris_energy%Ehartree_core = Ehartree_core
      IF (PRESENT(Exc)) harris_energy%Exc = Exc
      IF (PRESENT(integral_vxc)) harris_energy%integral_vxc = integral_vxc
      IF (PRESENT(Ecore)) harris_energy%Ecore = Ecore
      IF (PRESENT(Ecore_overlap)) harris_energy%Ecore_overlap = Ecore_overlap
      IF (PRESENT(Ecore_self)) harris_energy%Ecore_self = Ecore_self
      IF (PRESENT(EII)) harris_energy%EII = EII
      IF (PRESENT(Eharris_correction)) THEN
        harris_energy%Eharris_correction = Eharris_correction
      END IF
      IF (PRESENT(sum_of_eigenvalues)) THEN
        harris_energy%sum_of_eigenvalues = sum_of_eigenvalues
      END IF
      IF (PRESENT(Eharris)) harris_energy%Eharris = Eharris
    END IF
  END SUBROUTINE harris_energy_set
!***************************************************************************

!!****f* harris_energy_types/harris_energy_clear [1.0] *
!!
!!   NAME
!!     harris_energy_clear
!!
!!   SYNOPSIS
!!     SUBROUTINE harris_energy_clear(harris_energy, error)
!!       TYPE(harris_energy_type), POINTER :: harris_energy
!!       TYPE(cp_error_type), OPTIONAL, INTENT(INOUT) :: error
!!     END SUBROUTINE harris_energy_clear
!!
!!   FUNCTION
!!     Reinitializes the harris energy type
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - harris_energy: The harris energy type to be reinitialized
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!     For possible missing arguments see the attributes of harris_energy_type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     06.2005 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE harris_energy_clear(harris_energy, error)

    TYPE(harris_energy_type), POINTER        :: harris_energy
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_energy_clear', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(harris_energy), cp_failure_level, routineP, error, failure)
 
    IF (.NOT. failure) THEN
      harris_energy%Ehartree_elec      = 0.0_dp
      harris_energy%Ehartree_core      = 0.0_dp
      harris_energy%Exc                = 0.0_dp
      harris_energy%integral_vxc       = 0.0_dp
      harris_energy%Ecore              = 0.0_dp
      harris_energy%Ecore_overlap      = 0.0_dp
      harris_energy%Ecore_self         = 0.0_dp
      harris_energy%EII                = 0.0_dp
      harris_energy%Eharris_correction = 0.0_dp
      harris_energy%sum_of_eigenvalues = 0.0_dp
      harris_energy%Eharris            = 0.0_dp
    ELSE
      CALL harris_energy_create(harris_energy, error=error)
    END IF
  END SUBROUTINE harris_energy_clear
!***************************************************************************

!!****f* harris_energy_types/harris_energy_create [1.0] *
!!
!!   NAME
!!     harris_energy_create
!!
!!   SYNOPSIS
!!     SUBROUTINE harris_energy_create(harris_energy, error)
!!       TYPE(harris_energy_type), POINTER :: harris_energy
!!       TYPE(cp_error_type), OPTIONAL, INTENT(INOUT) :: error
!!     END SUBROUTINE harris_energy_create
!!
!!   FUNCTION
!!     Creates a harris energy type instance
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - harris_energy: The harris energy type to be created
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!     For possible missing arguments see the attributes of harris_energy_type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     06.2005 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE harris_energy_create(harris_energy, error)

    TYPE(harris_energy_type), POINTER        :: harris_energy
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_energy_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    ALLOCATE(harris_energy, stat=stat)
    CPPostcondition(stat == 0, cp_failure_level, routineP, error, failure)
   
    IF (.NOT. failure) THEN
      harris_energy%ref_count = 1
      last_harris_id=last_harris_id+1
      harris_energy%id_nr = last_harris_id
      CALL harris_energy_clear(harris_energy, error=error)
    END IF
  END SUBROUTINE harris_energy_create
!***************************************************************************

END MODULE harris_energy_types
