!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief The Harris energy type, which contains all the informations
!>      needed to perform and analyse the Harris functional energy correction.
!> \par History
!>      06.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
MODULE harris_energy_types

  USE f77_blas
  USE kinds,                           ONLY: dp
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'harris_energy_types'
  INTEGER, PRIVATE, SAVE :: last_harris_id=0

  ! *** Public data types ***
  PUBLIC :: harris_energy_type

  ! *** Public subroutines ***
  PUBLIC :: harris_energy_release,&
            harris_energy_clear,&
            harris_energy_create

! *****************************************************************************
!> \brief The Harris energy type
!> \param id_nr unique identifier
!> \param ref_count reference counter (see doc/ReferenceCounting.html)
!> \param Ehartree_elec The electronic hartree energy
!> \param Ehartree_core The hartree energy computed with rho_core
!> \param Exc The exchange-correlation energy
!> \param integral_vxc The xc potential integral
!> \param Ecore The core hamiltonian energy
!> \param Ecore_overlap The overalp energy
!> \param Ecore_self The self-interaction energy
!> \param EII The ion-ion Energy
!> \param Eharris_correction The energy correction due to the harris functional
!> \param sum_of_eigenvalues The sum of the eigenvalues
!> \param Eharris The actual harris functional energy
!> \par History
!>      06.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
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

! *****************************************************************************
!> \brief Type which allows the use of arrays of pointers
!> \param harris_energy_ptr Pointer to the actual harris energy type
!> \par History
!>      06.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  TYPE harris_energy_p_type
     TYPE(harris_energy_type), POINTER :: harris_energy_ptr
  END TYPE harris_energy_p_type

CONTAINS

! *****************************************************************************
!> \brief Retains a harris force type (see doc/ReferenceCounting.html)
!> \param harris_energy The harris energy type to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!>
!>      For possible missing arguments see the attributes of harris_energy_type
!> \par History
!>      06.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  SUBROUTINE harris_energy_retain(harris_energy, error)

    TYPE(harris_energy_type), POINTER        :: harris_energy
    TYPE(cp_error_type), INTENT(INOUT)       :: error

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

! *****************************************************************************
!> \brief Releases the given harris energy type (see doc/ReferenceCounting.html)
!> \param harris_energy The harris energy type to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!>
!>      For possible missing arguments see the attributes of harris_energy_type
!> \par History
!>      06.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  SUBROUTINE harris_energy_release(harris_energy, error)

    TYPE(harris_energy_type), POINTER        :: harris_energy
    TYPE(cp_error_type), INTENT(INOUT)       :: error

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

! *****************************************************************************
!> \brief Returns various attributes of the harris type
!> \param harris_energy The enquired harris energy type
!> \param id_nr unique identifier
!> \param Ehartree_elec The electronic hartree energy
!> \param Ehartree_core The hartree energy computed with rho_core
!> \param Exc The exchange-correlation energy
!> \param integral_vxc The xc potential integral
!> \param Ecore The core hamiltonian energy
!> \param Ecore_overlap The overalp energy
!> \param Ecore_self The self-interaction energy
!> \param EII The ion-ion Energy
!> \param Eharris_correction The energy correction due to the harris functional
!> \param sum_of_eigenvalues The sum of the eigenvalues
!> \param Eharris The actual harris functional energy
!>
!>      For possible missing arguments see the attributes of harris_energy_type
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      06.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  SUBROUTINE harris_energy_get(harris_energy, id_nr, Ehartree_elec, &
                               Ehartree_core, Exc, integral_vxc, Ecore, &
                               Ecore_overlap, Ecore_self, EII, &
                               Eharris_correction, sum_of_eigenvalues, &
                               Eharris, error)

    TYPE(harris_energy_type), POINTER        :: harris_energy
    INTEGER, INTENT(OUT), OPTIONAL           :: id_nr
    REAL(KIND=dp), INTENT(OUT), OPTIONAL :: Ehartree_elec, Ehartree_core, &
      Exc, integral_vxc, Ecore, Ecore_overlap, Ecore_self, EII, &
      Eharris_correction, sum_of_eigenvalues, Eharris
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_energy_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

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

! *****************************************************************************
!> \brief Sets various attributes of the harris energy type
!> \param harris_energy The enquired harris energy type
!> \param id_nr unique identifier
!> \param Ehartree_elec The electronic hartree energy
!> \param Ehartree_core The hartree energy computed with rho_core
!> \param Exc The exchange-correlation energy
!> \param integral_vxc The xc potential integral
!> \param Ecore The core hamiltonian energy
!> \param Ecore_overlap The overalp energy
!> \param Ecore_self The self-interaction energy
!> \param EII The ion-ion Energy
!> \param Eharris_correction The energy correction due to the harris functional
!> \param sum_of_eigenvalues The sum of the eigenvalues
!> \param Eharris The actual harris functional energy
!>
!>      For possible missing arguments see the attributes of harris_energy_type
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      06.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  SUBROUTINE harris_energy_set(harris_energy, id_nr, Ehartree_elec, &
                               Ehartree_core, Exc, integral_vxc, Ecore, &
                               Ecore_overlap, Ecore_self, EII, &
                               Eharris_correction, sum_of_eigenvalues, &
                               Eharris, error)

    TYPE(harris_energy_type), POINTER        :: harris_energy
    INTEGER, INTENT(IN), OPTIONAL            :: id_nr
    REAL(KIND=dp), INTENT(IN), OPTIONAL :: Ehartree_elec, Ehartree_core, Exc, &
      integral_vxc, Ecore, Ecore_overlap, Ecore_self, EII, &
      Eharris_correction, sum_of_eigenvalues, Eharris
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_energy_set', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

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

! *****************************************************************************
!> \brief Reinitializes the harris energy type
!> \param harris_energy The harris energy type to be reinitialized
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!>
!>      For possible missing arguments see the attributes of harris_energy_type
!> \par History
!>      06.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  SUBROUTINE harris_energy_clear(harris_energy, error)

    TYPE(harris_energy_type), POINTER        :: harris_energy
    TYPE(cp_error_type), INTENT(INOUT)       :: error

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

! *****************************************************************************
!> \brief Creates a harris energy type instance
!> \param harris_energy The harris energy type to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!>
!>      For possible missing arguments see the attributes of harris_energy_type
!> \par History
!>      06.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  SUBROUTINE harris_energy_create(harris_energy, error)

    TYPE(harris_energy_type), POINTER        :: harris_energy
    TYPE(cp_error_type), INTENT(INOUT)       :: error

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

END MODULE harris_energy_types
