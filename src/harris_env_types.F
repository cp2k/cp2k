!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2005  CP2K developers group                            !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****** cp2k/harris_env_types [1.0] *
!!
!!   NAME
!!     harris_env_types
!!
!!   FUNCTION
!!     The Harris environment type, which contains all the informations 
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
MODULE harris_env_types
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE kinds,                           ONLY: dp
 
  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='harris_env_types'
  INTEGER, PRIVATE, save :: last_harris_id=0

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
  TYPE harris_env_type
     INTEGER            :: id_nr, ref_count
     REAL(KIND=dp)      :: Ehartree_elec, &
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
!!     - harris_env_ptr: Pointer to the actual harris environment type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     06.2005 initial create [tdk]
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
!!       TYPE(harris_env_type), POINTER:: harris_env
!!       TYPE(cp_error_type), OPTIONAL, INTENT(INOUT):: error
!!     END SUBROUTINE harris_env_retain
!!
!!   FUNCTION
!!     Retains a harris environment type (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - harris_env_type: The harris environment type to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!     For possible missing arguments see the attributes of harris_env_type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     06.2005 initial create [tdk]
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
!!   SYNOPSIS
!!     SUBROUTINE harris_env_release(harris_env, error)
!!       TYPE(harris_env_type), POINTER:: harris_env
!!       TYPE(cp_error_type), OPTIONAL, INTENT(INOUT):: error
!!     END SUBROUTINE harris_env_release
!!
!!   FUNCTION
!!     Releases the given harris environment type 
!!     (see doc/ReferenceCounting.html)
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
!!     06.2005 initial create [tdk]
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
!!   SYNOPSIS
!!     SUBROUTINE ep_env_get(harris_env, id_nr, Ehartree_elec, Ehartree_core, 
!!                           Exc, integral_vxc, Ecore, Ecore_overlap,
!!                           Ecore_self, EII, Eharris_correction, 
!!                           sum_of_eigenvalues, Eharris, error)
!!       TYPE(harris_env_type), POINTER:: harris_env
!!       INTEGER, INTENT(OUT), OPTIONAL:: id_nr
!!       REAL(KIND=dp), INTENT(OUT), OPTIONAL :: Ehartreee_elec, 
!!                                               Ehartree_core, Exc,
!!                                               integral_vxc, Ecore,
!!                                               Ecore_overlap, Ecore_self
!!                                               EII, Eharris_correction,
!!                                               sum_of_eigenvalues, Eharris
!!       TYPE(cp_error_type), OPTIONAL, INTENT(INOUT):: error
!!     END SUBROUTINE ep_env_get
!!
!!   FUNCTION
!!     Returns various attributes of the harris environment.
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - harris_env: The enquired harris environment
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
!!     For possible missing arguments see the attributes of harris_env_type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     06.2005 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE harris_env_get(harris_env, id_nr, Ehartree_elec, Ehartree_core, &
                            Exc, integral_vxc, Ecore, Ecore_overlap, &
                            Ecore_self, EII, Eharris_correction, &
                            sum_of_eigenvalues, Eharris, error)
    TYPE(harris_env_type), POINTER           :: harris_env
    INTEGER, INTENT(OUT), OPTIONAL           :: id_nr
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: Ehartree_elec, &
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
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_env_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(harris_env), cp_failure_level, routineP, error, failure)
    IF (.not.failure) THEN
      CPPrecondition(harris_env%ref_count>0, cp_failure_level, routineP, error, failure)
    END IF

    IF (.NOT. failure) THEN
      IF (PRESENT(id_nr)) id_nr = harris_env%id_nr
      IF (PRESENT(Ehartree_elec)) Ehartree_elec = harris_env%Ehartree_elec
      IF (PRESENT(Ehartree_core)) Ehartree_core = harris_env%Ehartree_core
      IF (PRESENT(Exc)) Exc = harris_env%Exc
      IF (PRESENT(integral_vxc)) integral_vxc = harris_env%integral_vxc
      IF (PRESENT(Ecore)) Ecore = harris_env%Ecore
      IF (PRESENT(Ecore_overlap)) Ecore_overlap = harris_env%Ecore_overlap
      IF (PRESENT(Ecore_self)) Ecore_self = harris_env%Ecore_self
      IF (PRESENT(EII)) EII = harris_env%EII
      IF (PRESENT(Eharris_correction)) THEN
        Eharris_correction = harris_env%Eharris_correction
      END IF
      IF (PRESENT(sum_of_eigenvalues)) THEN
        sum_of_eigenvalues = harris_env%sum_of_eigenvalues
      END IF
      IF (PRESENT(Eharris)) Eharris = harris_env%Eharris
    END IF
  END SUBROUTINE harris_env_get
!***************************************************************************

!!****f* harris_env_types/harris_env_set [1.0] *
!!
!!   NAME
!!     harris_env_set
!!
!!   SYNOPSIS
!!     SUBROUTINE ep_env_set(harris_env, id_nr, Ehartree_elec, Ehartree_core, 
!!                           Exc, integral_vxc, Ecore, Ecore_overlap, 
!!                           Ecore_self, EII, Eharris_correction, 
!!                           sum_of_eigenvalues, Eharris, error)
!!       TYPE(harris_env_type), POINTER       :: harris_env
!!       INTEGER, INTENT(OUT), OPTIONAL       :: id_nr
!!       REAL(KIND=dp), INTENT(OUT), OPTIONAL :: Ehartreee_elec, 
!!                                               Ehartree_core, Exc,
!!                                               integral_vxc, Ecore,
!!                                               Ecore_overlap, Ecore_self,
!!                                               EII, Eharris_correction,
!!                                               sum_of_eigenvalues, Eharris
!!       TYPE(cp_error_type), OPTIONAL, INTENT(INOUT):: error
!!     END SUBROUTINE ep_env_set
!!
!!   FUNCTION
!!     Sets various attributes of the harris environment.
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - harris_env: The enquired harris environment
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
!!     For possible missing arguments see the attributes of harris_env_type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     06.2005 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE harris_env_set(harris_env, id_nr, Ehartree_elec, Ehartree_core, &
                            Exc, integral_vxc, Ecore, Ecore_overlap, &
                            Ecore_self, EII, Eharris_correction, &
                            sum_of_eigenvalues, Eharris, error)
    TYPE(harris_env_type), POINTER           :: harris_env
    INTEGER, INTENT(OUT), OPTIONAL           :: id_nr
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: Ehartree_elec, &
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
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_env_set', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(harris_env), cp_failure_level, routineP, error, failure)
    IF (.not.failure) THEN
      CPPrecondition(harris_env%ref_count>0, cp_failure_level, routineP, error, failure)
    END IF

    IF (.NOT. failure) THEN
      IF (PRESENT(id_nr)) harris_env%id_nr = id_nr
      IF (PRESENT(Ehartree_elec)) harris_env%Ehartree_elec = Ehartree_elec
      IF (PRESENT(Ehartree_core)) harris_env%Ehartree_core = Ehartree_core
      IF (PRESENT(Exc)) harris_env%Exc = Exc
      IF (PRESENT(integral_vxc)) harris_env%integral_vxc = integral_vxc
      IF (PRESENT(Ecore)) harris_env%Ecore = Ecore
      IF (PRESENT(Ecore_overlap)) harris_env%Ecore_overlap = Ecore_overlap
      IF (PRESENT(Ecore_self)) harris_env%Ecore_self = Ecore_self
      IF (PRESENT(EII)) harris_env%EII = EII
      IF (PRESENT(Eharris_correction)) THEN
        harris_env%Eharris_correction = Eharris_correction
      END IF
      IF (PRESENT(sum_of_eigenvalues)) THEN
        harris_env%sum_of_eigenvalues = sum_of_eigenvalues
      END IF
      If (PRESENT(Eharris)) harris_env%Eharris = Eharris
    END IF
  END SUBROUTINE harris_env_set
!***************************************************************************

!!****f* harris_env_types/harris_env_clear [1.0] *
!!
!!   NAME
!!     harris_env_clear
!!
!!   SYNOPSIS
!!     SUBROUTINE harris_env_clear(harris_env, error)
!!       TYPE(harris_env_type), POINTER:: harris_env
!!       TYPE(cp_error_type), OPTIONAL, INTENT(INOUT):: error
!!     END SUBROUTINE ep_env_init
!!
!!   FUNCTION
!!     Reinitializes the harris environment.
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - harris_env: The harris environment to be initialized
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!     For possible missing arguments see the attributes of harris_env_type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     06.2005 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE harris_env_clear(harris_env, error)
    TYPE(harris_env_type), POINTER           :: harris_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_env_clear', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(harris_env), cp_failure_level, routineP, error, failure)
 
    IF (.NOT. failure) THEN
      harris_env%Ehartree_elec = 0.0_dp
      harris_env%Ehartree_core = 0.0_dp
      harris_env%Exc = 0.0_dp
      harris_env%integral_vxc = 0.0_dp
      harris_env%Ecore = 0.0_dp
      harris_env%Ecore_overlap = 0.0_dp
      harris_env%Ecore_self = 0.0_dp
      harris_env%EII = 0.0_dp
      harris_env%Eharris_correction = 0.0_dp
      harris_env%sum_of_eigenvalues = 0.0_dp
      harris_env%Eharris = 0.0_dp
    ELSE
      CALL harris_env_create(harris_env, error=error)
    END IF
  END SUBROUTINE harris_env_clear
!***************************************************************************

!!****f* harris_env_types/harris_env_create [1.0] *
!!
!!   NAME
!!     harris_env_create
!!
!!   SYNOPSIS
!!     SUBROUTINE ep_env_create(harris_env, error)
!!       TYPE(harris_env_type), POINTER:: harris_env
!!       TYPE(cp_error_type), OPTIONAL, INTENT(INOUT):: error
!!     END SUBROUTINE ep_env_create
!!
!!   FUNCTION
!!     Creates a harris environment.
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - harris_env: The harris environment to be created
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!     For possible missing arguments see the attributes of harris_env_type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     06.2005 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE harris_env_create(harris_env, error)
    TYPE(harris_env_type), POINTER           :: harris_env
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
    IF (.NOT. failure) THEN
      harris_env%ref_count = 1
      last_harris_id=last_harris_id+1
      harris_env%id_nr = last_harris_id
      CALL harris_env_clear(harris_env, error=error)      
    END IF
  END SUBROUTINE harris_env_create
!***************************************************************************

END MODULE harris_env_types
