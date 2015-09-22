!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief container for information about total charges on the grids
!> \par History
!>      10.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE qs_charges_types
  
  USE kinds,                           ONLY: dp
#include "./common/cp_common_uses.f90"

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
!> \param nspins ...
!> \param total_rho_core_rspace ...
!> \param total_rho_gspace ...
!> \par History
!>      10.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE qs_charges_create(qs_charges,nspins,total_rho_core_rspace, &
     total_rho_gspace)
    TYPE(qs_charges_type), POINTER           :: qs_charges
    INTEGER, INTENT(in)                      :: nspins
    REAL(KIND=dp), INTENT(in), OPTIONAL      :: total_rho_core_rspace, &
                                                total_rho_gspace

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_charges_create', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  ALLOCATE(qs_charges)
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
  ALLOCATE(qs_charges%total_rho1_hard(nspins))
  qs_charges%total_rho1_hard(:) = 0.0_dp
  ALLOCATE(qs_charges%total_rho1_soft(nspins))
  qs_charges%total_rho1_soft(:) = 0.0_dp
  qs_charges%ref_count=1
END SUBROUTINE qs_charges_create

! *****************************************************************************
!> \brief retains the given qs_charges (see cp2k/doc/ReferenceCounting.html)
!> \param qs_charges the object to retain
!> \par History
!>      10.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE qs_charges_retain(qs_charges)
    TYPE(qs_charges_type), POINTER           :: qs_charges

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_charges_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPASSERT(ASSOCIATED(qs_charges))
  CPASSERT(qs_charges%ref_count>0)
  qs_charges%ref_count=qs_charges%ref_count+1
END SUBROUTINE qs_charges_retain

! *****************************************************************************
!> \brief releases the charges object (see cp2k/doc/ReferenceCounting.html)
!> \param qs_charges the object to be released
!> \par History
!>      10.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE qs_charges_release(qs_charges)
    TYPE(qs_charges_type), POINTER           :: qs_charges

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_charges_release', &
      routineP = moduleN//':'//routineN

  IF (ASSOCIATED(qs_charges)) THEN
     CPASSERT(qs_charges%ref_count>0)
     qs_charges%ref_count=qs_charges%ref_count-1
     IF (qs_charges%ref_count<1) THEN
        DEALLOCATE(qs_charges%total_rho1_hard)
        DEALLOCATE(qs_charges%total_rho1_soft)
        DEALLOCATE(qs_charges)
     END IF
  END IF
  NULLIFY(qs_charges)
END SUBROUTINE qs_charges_release

END MODULE qs_charges_types
