!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2006  CP2K developers group                            !
!-----------------------------------------------------------------------------!

!!****h* cp2k/tbmd_model_types [1.0] *
!!
!!   NAME
!!     tbmd_model_types
!!
!!   FUNCTION
!!     The tight-binding model type.
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
MODULE tbmd_model_types
  USE kinds,                           ONLY: dp

#include "cp_common_uses.h"
 
  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='tbmd_model_types'
  INTEGER, PRIVATE, save :: last_tbmd_id=0

  ! *** Public data types ***
  PUBLIC :: tbmd_model_type, &
            tbmd_model_p_type

  ! *** Public subroutines ***
  PUBLIC :: tbmd_model_retain, &
            tbmd_model_release, &
            tbmd_model_set,&
            tbmd_model_get, &
            tbmd_model_clear, &
            tbmd_model_create
!!***
!****************************************************************************

!!****s* tbmd_model_types/tbmd_model_type [1.0] *
!!
!!   NAME
!!     tbmd_model_type
!!
!!   FUNCTION
!!     The TBMD model type
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - id_nr: unique identifier
!!     - ref_count: reference counter (see doc/ReferenceCounting.html)
!!     - h_ssS: The hopping coefficient h(1)_ssSigma
!!     - h_spS: The hopping coefficient h(1)_spSigma
!!     - h_ppS: The hopping coefficient h(1)_ppSigma
!!     - h_ppP: The hopping coefficient h(1)_ppPi
!!     - n: The hopping exponent
!!     - r_0: Nearest neighbor distance 
!!     - Phi: The pair-potential parameter Phi(1)
!!     - m: The value of the pair-potnetial exponent
!!     - d_c: The pair-potential cutoff
!!     - m_c: Pair-potential exponent m_c
!!     - dEpsilon: s-p energy splitting
!!     - E_s: Energy of the s-Orbital
!!     - E_0: Energy shift constant
!!     - r_calpha: Individual cutoffs for each hopping parameter
!!     - n_calpha: Individual exponents for each hopping parameter
!!     - C14: Paramters C1-C4 of Kwon's 4th order polynomial for the 
!!            pair-repulsion
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!   SOURCE
!****************************************************************************
  TYPE tbmd_model_type
    INTEGER                                          :: id_nr, &
                                                        ref_count
    REAL(KIND=dp)                                    :: h_ssS, &
                                                        h_spS, &
                                                        h_ppS, &
                                                        h_ppP, &
                                                        n, &
                                                        r_0, &
                                                        Phi, &
                                                        m, &
                                                        d_c, &
                                                        m_c, &
                                                        dEpsilon, &
                                                        E_s, &
                                                        E_0
    REAL(KIND=dp), DIMENSION(4)                      :: r_calpha, &
                                                        n_calpha, &
                                                        C14

  END TYPE tbmd_model_type
!!***
!****************************************************************************

!!****s* tbmd_model_types/tbmd_model_p_type [1.0] *
!!
!!   NAME
!!     tbmd_model_p_type
!!
!!   FUNCTION
!!     Type which allows the use of arrays of pointers
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - tbmd_model_ptr: Pointer to the actual tbmd model type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!   SOURCE
!***************************************************************************
  TYPE tbmd_model_p_type
     TYPE(tbmd_model_type), POINTER :: tbmd_model_ptr
  END TYPE tbmd_model_p_type
!!***
!****************************************************************************
CONTAINS

!!****f* tbmd_model_types/tbmd_model_retain [1.0] *
!!
!!   NAME
!!     tbmd_model_retain
!!
!!   FUNCTION
!!     Retains a tbmd model type (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - tbmd_model: The tbmd model type to retain
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
  SUBROUTINE tbmd_model_retain(tbmd_model, error)

    TYPE(tbmd_model_type), POINTER           :: tbmd_model
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_model_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------
   
    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(tbmd_model), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
      CPPrecondition(tbmd_model%ref_count>0, cp_failure_level, routineP, error, failure)
      tbmd_model%ref_count = tbmd_model%ref_count+1
    END IF
  END SUBROUTINE tbmd_model_retain
!***************************************************************************

!!****f* tbmd_model_types/tbmd_model_release [1.0] *
!!
!!   NAME
!!     tbmd_model_release
!!
!!   FUNCTION
!!     Releases the given tbmd model type (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - tbmd_model: The tbmd model type to release
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
  SUBROUTINE tbmd_model_release(tbmd_model, error)

    TYPE(tbmd_model_type), POINTER                  :: tbmd_model
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                                      :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_model_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                         :: stat
    LOGICAL                                         :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.

    
    IF (ASSOCIATED(tbmd_model)) THEN
      CPPrecondition(tbmd_model%ref_count>0, cp_failure_level, routineP, error, failure)
      tbmd_model%ref_count = tbmd_model%ref_count-1
      IF (tbmd_model%ref_count<1) THEN
        !IF (ASSOCIATED(tbmd_model%f_total)) THEN
        !  DEALLOCATE(tbmd_model%f_total, stat=stat)
        !  CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)
        !END IF
        DEALLOCATE(tbmd_model, stat=stat)
        CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)
      END IF
    END IF
    NULLIFY(tbmd_model)

  END SUBROUTINE tbmd_model_release
!***************************************************************************

!!****f* tbmd_model_types/tbmd_model_get [1.0] *
!!
!!   NAME
!!     tbmd_model_get
!!
!!   FUNCTION
!!     Returns various attributes of the tbmd model type
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - tbmd_model: The enquired tbmd model type
!!     - id_nr: unique identifier
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handlin
!!     - h_ssS: The hopping coefficient h(1)_ssSigma
!!     - h_spS: The hopping coefficient h(1)_spSigma
!!     - h_ppS: The hopping coefficient h(1)_ppSigma
!!     - h_ppP: The hopping coefficient h(1)_ppPi
!!     - n: The hopping exponent
!!     - r_0: Nearest neighbor distance 
!!     - Phi: The pair-potential parameter Phi(1)
!!     - m: The value of the pair-potnetial exponent
!!     - d_c: The pair-potential cutoff
!!     - m_c: Pair-potential exponent m_c
!!     - dEpsilon: s-p energy splitting
!!     - E_s: Energy of the s-Orbital
!!     - E_0: Energy shift constant
!!     - r_calpha: Individual cutoffs for each hopping parameter
!!     - n_calpha: Individual exponents for each hopping parameter
!!     - C14: Paramters C1-C4 of Kwon's 4th order polynomial for the 
!!            pair-repulsion
!!
!!     For possible missing arguments see the attributes of tbmd_model_type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE tbmd_model_get(tbmd_model, id_nr, h_ssS, h_spS, h_ppS, h_ppP, &
                            n, r_0, Phi, m, d_c, m_c, dEpsilon, E_s, E_0, &
                            r_calpha, n_calpha, C14, error)

    TYPE(tbmd_model_type), POINTER                   :: tbmd_model
    INTEGER, INTENT(OUT), OPTIONAL                   :: id_nr
    REAL(KIND=dp)                                    :: h_ssS, &
                                                        h_spS, &
                                                        h_ppS, &
                                                        h_ppP, &
                                                        n, &
                                                        r_0, &
                                                        Phi, &
                                                        m, &
                                                        d_c, &
                                                        m_c, &
                                                        dEpsilon, &
                                                        E_s, &
                                                        E_0
    REAL(KIND=dp), DIMENSION(4)                      :: r_calpha, &
                                                        n_calpha, &
                                                        C14
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                                       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_model_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                          :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(tbmd_model), cp_failure_level, routineP, error, failure)
    IF (.NOT.failure) THEN
      CPPrecondition(tbmd_model%ref_count>0, cp_failure_level, routineP, error, failure)
    END IF

    !IF (.NOT. failure) THEN
    !  IF (PRESENT(f_total)) f_total = tbmd_model%f_total
    !END IF
  END SUBROUTINE tbmd_model_get
!***************************************************************************

!!****f* tbmd_model_types/tbmd_model_set [1.0] *
!!
!!   NAME
!!     tbmd_model_set
!!
!!   FUNCTION
!!     Sets various attributes of the tbmd model type
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - tbmd_model: The enquired tbmd model type
!!     - id_nr: unique identifier
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!     - h_ssS: The hopping coefficient h(1)_ssSigma
!!     - h_spS: The hopping coefficient h(1)_spSigma
!!     - h_ppS: The hopping coefficient h(1)_ppSigma
!!     - h_ppP: The hopping coefficient h(1)_ppPi
!!     - n: The hopping exponent
!!     - r_0: Nearest neighbor distance 
!!     - Phi: The pair-potential parameter Phi(1)
!!     - m: The value of the pair-potnetial exponent
!!     - d_c: The pair-potential cutoff
!!     - m_c: Pair-potential exponent m_c
!!     - dEpsilon: s-p energy splitting
!!     - E_s: Energy of the s-Orbital
!!     - E_0: Energy shift constant
!!     - r_calpha: Individual cutoffs for each hopping parameter
!!     - n_calpha: Individual exponents for each hopping parameter
!!     - C14: Paramters C1-C4 of Kwon's 4th order polynomial for the 
!!            pair-repulsion
!!
!!     For possible missing arguments see the attributes of tbmd_model_type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE tbmd_model_set(tbmd_model, id_nr, h_ssS, h_spS, h_ppS, h_ppP, &
                            n, r_0, Phi, m, d_c, m_c, dEpsilon, E_s, E_0, &
                            r_calpha, n_calpha, C14, error)

    TYPE(tbmd_model_type), POINTER                   :: tbmd_model
    INTEGER, INTENT(IN), OPTIONAL                    :: id_nr
    REAL(KIND=dp)                                    :: h_ssS, &
                                                        h_spS, &
                                                        h_ppS, &
                                                        h_ppP, &
                                                        n, &
                                                        r_0, &
                                                        Phi, &
                                                        m, &
                                                        d_c, &
                                                        m_c, &
                                                        dEpsilon, &
                                                        E_s, &
                                                        E_0
    REAL(KIND=dp), DIMENSION(4)                      :: r_calpha, &
                                                        n_calpha, &
                                                        C14
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                                       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_model_set', &
      routineP = moduleN//':'//routineN

    LOGICAL                                          :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(tbmd_model), cp_failure_level, routineP, error, failure)
    IF (.NOT.failure) THEN
      CPPrecondition(tbmd_model%ref_count>0, cp_failure_level, routineP, error, failure)
    END IF

    IF (.NOT. failure) THEN
      !IF (PRESENT(f_total)) tbmd_model%f_total = f_total
    END IF
  END SUBROUTINE tbmd_model_set
!***************************************************************************

!!****f* tbmd_model_types/tbmd_model_clear [1.0] *
!!
!!   NAME
!!     tbmd_model_clear
!!
!!   FUNCTION
!!     Reinitializes the tbmd model type
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - tbmd_model: The tbmd model type to be cleared
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
  SUBROUTINE tbmd_model_clear(tbmd_model, error)

    TYPE(tbmd_model_type), POINTER                  :: tbmd_model
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                                      :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_model_clear', &
      routineP = moduleN//':'//routineN

    LOGICAL                                         :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(tbmd_model), cp_failure_level, routineP, error, failure)

    IF (.NOT. failure) THEN
      !tbmd_tbmd%f_total(:) = 0.0_dp
    ELSE
      CALL tbmd_model_create(tbmd_model, error=error)
    END IF
  END SUBROUTINE tbmd_model_clear
!***************************************************************************

!!****f* tbmd_model_types/tbmd_model_create [1.0] *
!!
!!   NAME
!!     tbmd_model_create
!!
!!   FUNCTION
!!     Creates a tbmd model type instance
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - tbmd_model: The tbmd model type to be created
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
  SUBROUTINE tbmd_model_create(tbmd_model, error)

    TYPE(tbmd_model_type), POINTER                  :: tbmd_model
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                                      :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_model_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                         :: stat
    LOGICAL                                         :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.

    !ALLOCATE(tbmd_model%f_total(3), stat=stat)
    !CPPostcondition(stat == 0, cp_failure_level, routineP, error, failure)

    IF (.NOT. failure) THEN
      tbmd_model%ref_count = 1
      last_tbmd_id=last_tbmd_id+1
      tbmd_model%id_nr = last_tbmd_id
      CALL tbmd_model_clear(tbmd_model, error=error)
    END IF
  END SUBROUTINE tbmd_model_create
!***************************************************************************

END MODULE tbmd_model_types
