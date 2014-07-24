!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Definition of the semi empirical multipole integral expansions types
!> \author Teodoro Laino [tlaino] - 08.2008 Zurich University
! *****************************************************************************
MODULE semi_empirical_mpole_types
  USE f77_blas
  USE kinds,                           ONLY: dp
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'semi_empirical_mpole_types'

! *****************************************************************************
!> \brief Semi-empirical integral multipole expansion type
!> \author Teodoro Laino [tlaino] - 08.2008 Zurich University
! *****************************************************************************
  TYPE semi_empirical_mpole_type
     LOGICAL, DIMENSION(3)                    :: task
     INTEGER                                  :: indi, indj
     REAL(KIND=dp)                            :: c
     REAL(KIND=dp), DIMENSION(3)              :: d
     REAL(KIND=dp), DIMENSION(3,3)            :: qc ! quadrupole cartesian
     REAL(KIND=dp), DIMENSION(5)              :: qs ! quadrupole spherical
     ! alternative definition used in GKS integral routines
     REAL(KIND=dp)                            :: cs
     REAL(KIND=dp), DIMENSION(3)              :: ds
     REAL(KIND=dp), DIMENSION(3,3)            :: qq ! quadrupole cartesian
  END TYPE semi_empirical_mpole_type

! *****************************************************************************
!> \brief Semi-empirical integral multipole expansion type - pointer type
!> \author Teodoro Laino [tlaino] - 08.2008 Zurich University
! *****************************************************************************
  TYPE semi_empirical_mpole_p_type
     TYPE(semi_empirical_mpole_type), POINTER :: mpole
  END TYPE semi_empirical_mpole_p_type

! *****************************************************************************
!> \brief Global Multipolar NDDO information type
!> \author Teodoro Laino [tlaino] - 08.2008 Zurich University
! *****************************************************************************
  TYPE nddo_mpole_type
     REAL(KIND=dp), DIMENSION(:), POINTER     :: charge,     efield0
     REAL(KIND=dp), DIMENSION(:,:), POINTER   :: dipole,     efield1, efield2
     REAL(KIND=dp), DIMENSION(:,:,:), POINTER :: quadrupole
  END TYPE nddo_mpole_type


  PUBLIC :: semi_empirical_mpole_type,&
            semi_empirical_mpole_p_type,&
            semi_empirical_mpole_p_create,&
            semi_empirical_mpole_p_release,&
            nddo_mpole_type,&
            nddo_mpole_create,&
            nddo_mpole_release

CONTAINS

! *****************************************************************************
!> \brief Allocate semi-empirical mpole type
!> \param mpole ...
!> \param ndim ...
!> \param error ...
!> \author Teodoro Laino [tlaino] - 08.2008 Zurich University
! *****************************************************************************
  SUBROUTINE semi_empirical_mpole_p_create(mpole, ndim, error)
    TYPE(semi_empirical_mpole_p_type), &
      DIMENSION(:), POINTER                  :: mpole
    INTEGER, INTENT(IN)                      :: ndim
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'semi_empirical_mpole_p_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(.NOT.ASSOCIATED(mpole),cp_failure_level,routineP,error,failure)
    ALLOCATE (mpole(ndim),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DO i = 1, ndim
       NULLIFY(mpole(i)%mpole)
       CALL semi_empirical_mpole_create(mpole(i)%mpole, error=error)
    END DO

  END SUBROUTINE semi_empirical_mpole_p_create

! *****************************************************************************
!> \brief Deallocate the semi-empirical mpole type
!> \param mpole ...
!> \param error ...
!> \author Teodoro Laino [tlaino] - 08.2008 Zurich University
! *****************************************************************************
  SUBROUTINE semi_empirical_mpole_p_release(mpole, error)
    TYPE(semi_empirical_mpole_p_type), &
      DIMENSION(:), POINTER                  :: mpole
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'semi_empirical_mpole_p_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(mpole)) THEN
       DO i = 1, SIZE(mpole)
          CALL semi_empirical_mpole_release(mpole(i)%mpole,error=error)
       END DO
       DEALLOCATE (mpole, stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE semi_empirical_mpole_p_release

! *****************************************************************************
!> \brief Allocate semi-empirical mpole type
!> \param mpole ...
!> \param error ...
!> \author Teodoro Laino [tlaino] - 08.2008 Zurich University
! *****************************************************************************
  SUBROUTINE semi_empirical_mpole_create(mpole, error)
    TYPE(semi_empirical_mpole_type), POINTER :: mpole
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'semi_empirical_mpole_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(.NOT.ASSOCIATED(mpole),cp_failure_level,routineP,error,failure)
    ALLOCATE (mpole,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    mpole%task = .FALSE.
    mpole%indi = 0
    mpole%indj = 0
    mpole%c    = HUGE(0.0_dp)
    mpole%d    = HUGE(0.0_dp)
    mpole%qc   = HUGE(0.0_dp)
    mpole%qs   = HUGE(0.0_dp)
    mpole%cs   = HUGE(0.0_dp)
    mpole%ds   = HUGE(0.0_dp)
    mpole%qq   = HUGE(0.0_dp)
  END SUBROUTINE semi_empirical_mpole_create

! *****************************************************************************
!> \brief Deallocate the semi-empirical mpole type
!> \param mpole ...
!> \param error ...
!> \author Teodoro Laino [tlaino] - 08.2008 Zurich University
! *****************************************************************************
  SUBROUTINE semi_empirical_mpole_release(mpole, error)
    TYPE(semi_empirical_mpole_type), POINTER :: mpole
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'semi_empirical_mpole_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(mpole)) THEN
       DEALLOCATE (mpole, stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE semi_empirical_mpole_release

! *****************************************************************************
!> \brief Allocate NDDO multipole type
!> \param nddo_mpole ...
!> \param error ...
!> \author Teodoro Laino [tlaino] - 08.2008 Zurich University
! *****************************************************************************
  SUBROUTINE nddo_mpole_create(nddo_mpole, error)
    TYPE(nddo_mpole_type), POINTER           :: nddo_mpole
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'nddo_mpole_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(.NOT.ASSOCIATED(nddo_mpole),cp_failure_level,routineP,error,failure)
    ALLOCATE (nddo_mpole,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    NULLIFY(nddo_mpole%charge)
    NULLIFY(nddo_mpole%dipole)
    NULLIFY(nddo_mpole%quadrupole)
    NULLIFY(nddo_mpole%efield0)
    NULLIFY(nddo_mpole%efield1)
    NULLIFY(nddo_mpole%efield2)
  END SUBROUTINE nddo_mpole_create

! *****************************************************************************
!> \brief Deallocate NDDO multipole type
!> \param nddo_mpole ...
!> \param error ...
!> \author Teodoro Laino [tlaino] - 08.2008 Zurich University
! *****************************************************************************
  SUBROUTINE nddo_mpole_release(nddo_mpole, error)
    TYPE(nddo_mpole_type), POINTER           :: nddo_mpole
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'nddo_mpole_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(nddo_mpole)) THEN
       IF (ASSOCIATED(nddo_mpole%charge)) THEN
          DEALLOCATE(nddo_mpole%charge,stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(nddo_mpole%dipole)) THEN
          DEALLOCATE(nddo_mpole%dipole,stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(nddo_mpole%quadrupole)) THEN
          DEALLOCATE(nddo_mpole%quadrupole,stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(nddo_mpole%efield0)) THEN
          DEALLOCATE(nddo_mpole%efield0,stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(nddo_mpole%efield1)) THEN
          DEALLOCATE(nddo_mpole%efield1,stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(nddo_mpole%efield2)) THEN
          DEALLOCATE(nddo_mpole%efield2,stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       DEALLOCATE (nddo_mpole, stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE nddo_mpole_release

END MODULE semi_empirical_mpole_types
