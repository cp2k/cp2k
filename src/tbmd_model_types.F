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

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tbmd_model_types'
  INTEGER, PRIVATE, SAVE :: last_tbmd_id=0

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
!!     E_tot = Sum_n[2*<Psi_n|H_TB|Psi_n> + E_rep + E_0*N_at]
!!     My general form for hopping parameters and pair potentials: 
!!     h_alpha(r) = h_alpha(r_0)(r_0/r)**n * exp{n*[-(r/r_calpha)**n_calpha + 
!!                                               (r_0/r_calpha)**n_calpha]}
!!     E_rep = Sum_i[f*Sum_j[Phi(r_ij)]], 
!!             with f(x)=C_1*x + C_2*x**2 + C_3*x**3 + C_4*x**4 and
!!     Phi(r) = Phi(1)*(r_0/r)**m * exp{m*[-(r/d_c)**m_c + (r_0/d_c)**m_c]}
!!
!!   ATTRIBUTES
!!     - id_nr: unique identifier
!!     - ref_count: reference counter (see doc/ReferenceCounting.html)
!!     - tb_model_name: Name of the tight-binding parametrization
!!     - noisy_tb: Logical which indicates if our linear scaling 
!!                 stoachstic tight-binding method is used.
!!     - h_ssS: The hopping coefficient h(1)_ssSigma
!!     - h_spS: The hopping coefficient h(1)_spSigma
!!     - h_ppS: The hopping coefficient h(1)_ppSigma
!!     - h_ppP: The hopping coefficient h(1)_ppPi
!!     - n: The hopping exponent
!!     - r_0: Nearest neighbor distance 
!!     - Phi: The pair-potential parameter Phi(1)
!!     - m: The value of the pair-potential exponent
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
                                                        ref_count, &
                                                        tb_model_name
    LOGICAL                                          :: noisy_tb
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

    TYPE(tbmd_model_type), POINTER           :: tbmd_model
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_model_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.

    
    IF (ASSOCIATED(tbmd_model)) THEN
      CPPrecondition(tbmd_model%ref_count>0, cp_failure_level, routineP, error, failure)
      tbmd_model%ref_count = tbmd_model%ref_count-1
      IF (tbmd_model%ref_count<1) THEN
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
!!     - tb_model_name: Name of the tight-binding parametrization
!!     - noisy_tb: Logical which indicates if our linear scaling 
!!                 stoachstic tight-binding method is used.
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
                            r_calpha, n_calpha, C14, tb_model_name, &
                            noisy_tb, error)

    TYPE(tbmd_model_type), POINTER           :: tbmd_model
    INTEGER, INTENT(OUT), OPTIONAL           :: id_nr
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: h_ssS, h_spS, h_ppS, h_ppP, &
                                                n, r_0, Phi, m, d_c, m_c, &
                                                dEpsilon, E_s, E_0
    REAL(KIND=dp), DIMENSION(4), &
      INTENT(OUT), OPTIONAL                  :: r_calpha, n_calpha, C14
    INTEGER, INTENT(OUT), OPTIONAL           :: tb_model_name
    LOGICAL, INTENT(OUT), OPTIONAL           :: noisy_tb
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_model_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(tbmd_model), cp_failure_level, routineP, error, failure)
    IF (.NOT.failure) THEN
      CPPrecondition(tbmd_model%ref_count>0, cp_failure_level, routineP, error, failure)
    END IF

    IF (.NOT. failure) THEN
      IF (PRESENT(id_nr)) id_nr = tbmd_model%id_nr
      IF (PRESENT(tb_model_name)) THEN 
        tb_model_name = tbmd_model%tb_model_name
      END IF
      IF (PRESENT(noisy_tb)) noisy_tb = tbmd_model%noisy_tb
      IF (PRESENT(h_ssS)) h_ssS = tbmd_model%h_ssS
      IF (PRESENT(h_spS)) h_spS = tbmd_model%h_spS
      IF (PRESENT(h_ppS)) h_ppS = tbmd_model%h_ppS
      IF (PRESENT(h_ppP)) h_ppP = tbmd_model%h_ppP
      IF (PRESENT(n)) n = tbmd_model%n
      IF (PRESENT(r_0)) r_0 = tbmd_model%r_0
      IF (PRESENT(Phi)) Phi = tbmd_model%Phi
      IF (PRESENT(m)) m = tbmd_model%m
      IF (PRESENT(d_c)) d_c = tbmd_model%d_c
      IF (PRESENT(m_c)) m_c = tbmd_model%m_c
      IF (PRESENT(dEpsilon)) dEpsilon = tbmd_model%dEpsilon
      IF (PRESENT(E_s)) E_s = tbmd_model%E_s
      IF (PRESENT(E_0)) E_0 = tbmd_model%E_0
      IF (PRESENT(r_calpha)) r_calpha = tbmd_model%r_calpha
      IF (PRESENT(n_calpha)) n_calpha = tbmd_model%n_calpha
      IF (PRESENT(C14)) C14 = tbmd_model%C14
    END IF
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
!!     - tb_model_name: Name of the tight-binding parametrization
!!     - noisy_tb: Logical which indicates if our linear scaling 
!!                 stoachstic tight-binding method is used.
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
                            r_calpha, n_calpha, C14, tb_model_name, &
                            noisy_tb, error)

    TYPE(tbmd_model_type), POINTER           :: tbmd_model
    INTEGER, INTENT(IN), OPTIONAL            :: id_nr
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: h_ssS, h_spS, h_ppS, h_ppP, &
                                                n, r_0, Phi, m, d_c, m_c, &
                                                dEpsilon, E_s, E_0
    REAL(KIND=dp), DIMENSION(4), &
      INTENT(IN), OPTIONAL                   :: r_calpha, n_calpha, C14
    INTEGER, INTENT(IN), OPTIONAL            :: tb_model_name
    LOGICAL, INTENT(IN), OPTIONAL            :: noisy_tb
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_model_set', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(tbmd_model), cp_failure_level, routineP, error, failure)
    IF (.NOT.failure) THEN
      CPPrecondition(tbmd_model%ref_count>0, cp_failure_level, routineP, error, failure)
    END IF

    IF (.NOT. failure) THEN
      IF (PRESENT(id_nr)) tbmd_model%id_nr = id_nr
      IF (PRESENT(tb_model_name)) THEN
        tbmd_model%tb_model_name = tb_model_name
      END IF
      IF (PRESENT(noisy_tb)) tbmd_model%noisy_tb = noisy_tb
      IF (PRESENT(h_ssS)) tbmd_model%h_ssS = h_ssS
      IF (PRESENT(h_spS)) tbmd_model%h_spS = h_spS
      IF (PRESENT(h_ppS)) tbmd_model%h_ppS = h_ppS
      IF (PRESENT(h_ppP)) tbmd_model%h_ppP = h_ppP
      IF (PRESENT(n)) tbmd_model%n = n
      IF (PRESENT(r_0)) tbmd_model%r_0 = r_0
      IF (PRESENT(Phi)) tbmd_model%Phi = Phi
      IF (PRESENT(m)) tbmd_model%m = m
      IF (PRESENT(d_c)) tbmd_model%d_c = d_c
      IF (PRESENT(m_c)) tbmd_model%m_c = m_c
      IF (PRESENT(dEpsilon)) tbmd_model%dEpsilon = dEpsilon
      IF (PRESENT(E_s)) tbmd_model%E_s = E_s
      IF (PRESENT(E_0)) tbmd_model%E_0 = E_0
      IF (PRESENT(r_calpha)) tbmd_model%r_calpha = r_calpha
      IF (PRESENT(n_calpha)) tbmd_model%n_calpha = n_calpha
      IF (PRESENT(C14)) tbmd_model%C14 = C14
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

    TYPE(tbmd_model_type), POINTER           :: tbmd_model
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_model_clear', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(tbmd_model), cp_failure_level, routineP, error, failure)

    IF (.NOT. failure) THEN
      tbmd_model%tb_model_name = 0
      tbmd_model%noisy_tb      = .FALSE.
      tbmd_model%h_ssS         = 0.0_dp
      tbmd_model%h_spS         = 0.0_dp
      tbmd_model%h_ppS         = 0.0_dp
      tbmd_model%h_ppP         = 0.0_dp
      tbmd_model%n             = 0.0_dp
      tbmd_model%r_0           = 0.0_dp
      tbmd_model%Phi           = 0.0_dp
      tbmd_model%m             = 0.0_dp
      tbmd_model%d_c           = 0.0_dp
      tbmd_model%m_c           = 0.0_dp
      tbmd_model%dEpsilon      = 0.0_dp
      tbmd_model%E_s           = 0.0_dp
      tbmd_model%E_0           = 0.0_dp
      tbmd_model%r_calpha(:)   = 0.0_dp
      tbmd_model%n_calpha(:)   = 0.0_dp
      tbmd_model%C14(:)        = 0.0_dp
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

    TYPE(tbmd_model_type), POINTER           :: tbmd_model
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_model_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.

    ALLOCATE(tbmd_model, stat=stat)
    CPPostcondition(stat == 0, cp_failure_level, routineP, error, failure)

    IF (.NOT. failure) THEN
      tbmd_model%ref_count = 1
      last_tbmd_id=last_tbmd_id+1
      tbmd_model%id_nr = last_tbmd_id
      CALL tbmd_model_clear(tbmd_model, error=error)
    END IF
  END SUBROUTINE tbmd_model_create
!***************************************************************************

END MODULE tbmd_model_types
