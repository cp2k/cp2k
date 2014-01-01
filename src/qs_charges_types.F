!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief container for information about total charges on the grids
!> \par History
!>      10.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE qs_charges_types
  USE f77_blas
  USE kinds,                           ONLY: dp
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_charges_types'

  PUBLIC :: qs_charges_type
  PUBLIC :: qs_charges_create, qs_charges_retain, qs_charges_release
!***

! *****************************************************************************
!> \brief Container for information about total charges on the grids
!> \param total_rho_core_rspace total charge on the rho_core grid
!> \param total_rho_rspace total charge in the real space
!> \param total_rho_gspace total charge in the g space
!> \note
!>      this type is loosing the reason to exist...
!> \par History
!>      10.2002 created [fawzi]
!>      11.2002 moved total_rho_elec_rspace to qs_rho_type
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE qs_charges_type
     INTEGER :: ref_count
     REAL(KIND = dp) :: total_rho_core_rspace, total_rho_gspace
     REAL(KIND = dp) :: total_rho0_soft_rspace, total_rho0_hard_lebedev
     REAL(KIND = dp) :: total_rho_soft_gspace
     REAL(KIND = dp), DIMENSION(:), POINTER  :: total_rho1_hard,&
                                                total_rho1_soft
     REAL(KIND = dp) :: background
  END TYPE qs_charges_type

CONTAINS

! *****************************************************************************
!> \brief creates a charges object
!> \param qs_charges the charges object to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE qs_charges_create(qs_charges,nspins,total_rho_core_rspace, &
     total_rho_gspace, error)
    TYPE(qs_charges_type), POINTER           :: qs_charges
    INTEGER, INTENT(in)                      :: nspins
    REAL(KIND=dp), INTENT(in), OPTIONAL      :: total_rho_core_rspace, &
                                                total_rho_gspace
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_charges_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.

  ALLOCATE(qs_charges, stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     qs_charges%total_rho_core_rspace=0.0_dp
     IF (PRESENT(total_rho_core_rspace)) &
          qs_charges%total_rho_core_rspace=total_rho_core_rspace
     qs_charges%total_rho_gspace=0.0_dp
     IF (PRESENT(total_rho_gspace)) &
          qs_charges%total_rho_gspace=total_rho_gspace
     qs_charges%total_rho_soft_gspace = 0.0_dp
     qs_charges%total_rho0_hard_lebedev = 0.0_dp
     qs_charges%total_rho_soft_gspace = 0.0_dp
     qs_charges%background = 0.0_dp
     ALLOCATE(qs_charges%total_rho1_hard(nspins),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     qs_charges%total_rho1_hard(:) = 0.0_dp
     ALLOCATE(qs_charges%total_rho1_soft(nspins),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     qs_charges%total_rho1_soft(:) = 0.0_dp
     qs_charges%ref_count=1
  END IF
END SUBROUTINE qs_charges_create

! *****************************************************************************
!> \brief retains the given qs_charges (see cp2k/doc/ReferenceCounting.html)
!> \param qs_charges the object to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE qs_charges_retain(qs_charges, error)
    TYPE(qs_charges_type), POINTER           :: qs_charges
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_charges_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(qs_charges),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(qs_charges%ref_count>0,cp_failure_level,routineP,error,failure)
     qs_charges%ref_count=qs_charges%ref_count+1
  END IF
END SUBROUTINE qs_charges_retain

! *****************************************************************************
!> \brief releases the charges object (see cp2k/doc/ReferenceCounting.html)
!> \param qs_charges the object to be released
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE qs_charges_release(qs_charges, error)
    TYPE(qs_charges_type), POINTER           :: qs_charges
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_charges_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.

  IF (ASSOCIATED(qs_charges)) THEN
     CPPrecondition(qs_charges%ref_count>0,cp_failure_level,routineP,error,failure)
     qs_charges%ref_count=qs_charges%ref_count-1
     IF (qs_charges%ref_count<1) THEN
        DEALLOCATE(qs_charges%total_rho1_hard,stat=stat)
        DEALLOCATE(qs_charges%total_rho1_soft,stat=stat)
        DEALLOCATE(qs_charges,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(qs_charges)
END SUBROUTINE qs_charges_release

END MODULE qs_charges_types
