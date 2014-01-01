!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      none
!> \author CJM
! *****************************************************************************
MODULE scp_force_types
  USE f77_blas
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_size
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'scp_force_types'

  TYPE scp_force_type
    REAL(KIND=dp), DIMENSION(:,:), POINTER :: f_scp, f_scp_disp, f_total
  END TYPE scp_force_type

  PUBLIC :: scp_force_type

  PUBLIC :: allocate_scp_force,&
            deallocate_scp_force,&
            init_scp_force

CONTAINS

! *****************************************************************************
!> \brief   Allocate a scp force data structure.
!> \version 1.0
! *****************************************************************************
  SUBROUTINE allocate_scp_force(scp_forces,natom_of_kind)

    TYPE(scp_force_type), DIMENSION(:), &
      POINTER                                :: scp_forces
    INTEGER, DIMENSION(:), INTENT(IN)        :: natom_of_kind

    CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_scp_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, istat, n, nkind

    IF (ASSOCIATED(scp_forces)) THEN
      CALL deallocate_scp_force(scp_forces)
    END IF

    nkind = SIZE(natom_of_kind)

    ALLOCATE (scp_forces(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "scp_forces",int_size*nkind)

    DO ikind=1,nkind

      n = natom_of_kind(ikind)

      ALLOCATE (scp_forces(ikind)%f_scp (3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,&
                         "scp_forces(ikind)%f_scp",dp_size*3*n)
      END IF

      ALLOCATE (scp_forces(ikind)%f_scp_disp(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,&
                         "scp_forces(ikind)%f_scp_disp",dp_size*3*n)
      END IF

      ALLOCATE (scp_forces(ikind)%f_total(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,&
                         "scp_forces(ikind)%f_total",dp_size*3*n)
      END IF

    END DO

  END SUBROUTINE allocate_scp_force

! *****************************************************************************
!> \brief   Deallocate a scp force data structure.
!> \version 1.0
! *****************************************************************************
  SUBROUTINE deallocate_scp_force(scp_forces)

    TYPE(scp_force_type), DIMENSION(:), &
      POINTER                                :: scp_forces

    CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_scp_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, istat, nkind

    IF (.NOT.ASSOCIATED(scp_forces)) THEN
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The scp_forces pointer is not associated "//&
                        "and cannot be deallocated")
    END IF

    nkind = SIZE(scp_forces)

    DO ikind=1,nkind

      IF (ASSOCIATED(scp_forces(ikind)%f_scp)) THEN
        DEALLOCATE (scp_forces(ikind)%f_scp,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,__LINE__,"scp_forces%f_scp")
        END IF
      END IF

      IF (ASSOCIATED(scp_forces(ikind)%f_scp_disp)) THEN
        DEALLOCATE (scp_forces(ikind)%f_scp_disp,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,__LINE__,"scp_forces%f_scp_disp")
        END IF
      END IF


      IF (ASSOCIATED(scp_forces(ikind)%f_total)) THEN
        DEALLOCATE (scp_forces(ikind)%f_total,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,__LINE__,"scp_forces%f_total")
        END IF
      END IF

    END DO

    DEALLOCATE (scp_forces,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"scp_forces")

  END SUBROUTINE deallocate_scp_force

! *****************************************************************************
!> \brief   Initialize a scp force data structure.
!> \version 1.0
! *****************************************************************************
  SUBROUTINE init_scp_force(scp_forces)

    TYPE(scp_force_type), DIMENSION(:), &
      POINTER                                :: scp_forces

    CHARACTER(LEN=*), PARAMETER :: routineN = 'init_scp_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind

    IF (.NOT.ASSOCIATED(scp_forces)) THEN
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The scp_forces pointer is not associated "//&
                        "and cannot be initialized")
    END IF

    DO ikind=1,SIZE(scp_forces)
      scp_forces(ikind)%f_scp(:,:) = 0.0_dp
      scp_forces(ikind)%f_total(:,:) = 0.0_dp
      scp_forces(ikind)%f_scp_disp(:,:) = 0.0_dp
    END DO

  END SUBROUTINE init_scp_force

END MODULE scp_force_types

