!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      none
!> \author gloria (01.10.2002)
! *****************************************************************************
MODULE kg_force_types

  USE f77_blas
  USE kinds,                           ONLY: dp
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
#include "cp_common_uses.h"

  IMPLICIT NONE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kg_force_types'
  PRIVATE

! *****************************************************************************
  TYPE kg_force_type
    REAL(KIND=dp), DIMENSION(:,:), POINTER :: f_rspace_core,&
                                         f_hartree_core,&
                                         f_rho,&
                                         f_xc,&
                                         f_bc,&
                                         f_hartree,&
                                         f_ppl,&
                                         f_total
  END TYPE kg_force_type

! *** Public data types ***

  PUBLIC :: kg_force_type

! *** Public subroutines ***

  PUBLIC :: allocate_kg_force,&
            deallocate_kg_force,&
            init_kg_force

CONTAINS

! *****************************************************************************
!> \brief   Allocate a kg force data structure.
!> \version 1.0
! *****************************************************************************
  SUBROUTINE allocate_kg_force(kg_forces,natom_of_kind)
    TYPE(kg_force_type), DIMENSION(:), &
      POINTER                                :: kg_forces
    INTEGER, DIMENSION(:), INTENT(IN)        :: natom_of_kind

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_kg_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, istat, n, nkind

    IF (ASSOCIATED(kg_forces)) THEN
      CALL deallocate_kg_force(kg_forces)
    END IF

    nkind = SIZE(natom_of_kind)

    ALLOCATE (kg_forces(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"kg_forces",nkind)

    DO ikind=1,nkind

      n = natom_of_kind(ikind)

      ALLOCATE (kg_forces(ikind)%f_rspace_core(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineP,"kg_forces(ikind)%f_rspace_core",3*n)
      END IF

      ALLOCATE (kg_forces(ikind)%f_hartree_core(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineP,"kg_forces(ikind)%f_hartree_core",3*n)
      END IF

      ALLOCATE (kg_forces(ikind)%f_rho(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineP,"kg_forces(ikind)%f_rho",3*n)
      END IF

      ALLOCATE (kg_forces(ikind)%f_hartree(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineP,"kg_forces(ikind)%f_hartree",3*n)
      END IF

      ALLOCATE (kg_forces(ikind)%f_xc(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineP,"kg_forces(ikind)%f_xc",3*n)
      END IF

      ALLOCATE (kg_forces(ikind)%f_bc(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineP,"kg_forces(ikind)%f_bc",3*n)
      END IF

      ALLOCATE (kg_forces(ikind)%f_ppl(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineP,"kg_forces(ikind)%f_ppl",3*n)
      END IF

      ALLOCATE (kg_forces(ikind)%f_total(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineP,"kg_forces(ikind)%f_total",3*n)
      END IF

    END DO

  END SUBROUTINE allocate_kg_force 

! *****************************************************************************
!> \brief   Deallocate a kg force data structure.
!> \author  gt
!> \date    30.09.2002
!> \version 1.0
! *****************************************************************************
  SUBROUTINE deallocate_kg_force(kg_forces)
    TYPE(kg_force_type), DIMENSION(:), &
      POINTER                                :: kg_forces

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_kg_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, istat, nkind

    IF (.NOT.ASSOCIATED(kg_forces)) THEN
      CALL stop_program(routineP,"The kg_forces pointer is not associated "//&
                                "and cannot be deallocated")
    END IF

    nkind = SIZE(kg_forces)

    DO ikind=1,nkind

      IF (ASSOCIATED(kg_forces(ikind)%f_rspace_core)) THEN
        DEALLOCATE (kg_forces(ikind)%f_rspace_core,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineP,"kg_forces%f_rspace_core")
        END IF
      END IF

      IF (ASSOCIATED(kg_forces(ikind)%f_hartree_core)) THEN
        DEALLOCATE (kg_forces(ikind)%f_hartree_core,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineP,"kg_forces%f_hartree_core")
        END IF
      END IF

      IF (ASSOCIATED(kg_forces(ikind)%f_rho)) THEN
        DEALLOCATE (kg_forces(ikind)%f_rho,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineP,"kg_forces%f_rho")
        END IF
      END IF

      IF (ASSOCIATED(kg_forces(ikind)%f_ppl)) THEN
        DEALLOCATE (kg_forces(ikind)%f_ppl,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineP,"kg_forces%f_ppl")
        END IF
      END IF

      IF (ASSOCIATED(kg_forces(ikind)%f_hartree)) THEN
        DEALLOCATE (kg_forces(ikind)%f_hartree,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineP,"kg_forces%f_hartree")
        END IF
      END IF

      IF (ASSOCIATED(kg_forces(ikind)%f_xc)) THEN
        DEALLOCATE (kg_forces(ikind)%f_xc,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineP,"kg_forces%f_xc")
        END IF
      END IF

      IF (ASSOCIATED(kg_forces(ikind)%f_bc)) THEN
        DEALLOCATE (kg_forces(ikind)%f_bc,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineP,"kg_forces%f_bc")
        END IF
      END IF

      IF (ASSOCIATED(kg_forces(ikind)%f_total)) THEN
        DEALLOCATE (kg_forces(ikind)%f_total,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineP,"kg_forces%f_total")
        END IF
      END IF

    END DO

    DEALLOCATE (kg_forces,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"kg_forces")

  END SUBROUTINE deallocate_kg_force

! *****************************************************************************
!> \brief   Initialize a kg force data structure.
!> \author  gt
!> \date    01.10.2002
!> \version 1.0
! *****************************************************************************
  SUBROUTINE init_kg_force(kg_forces)
    TYPE(kg_force_type), DIMENSION(:), &
      POINTER                                :: kg_forces

    CHARACTER(len=*), PARAMETER :: routineN = 'init_kg_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind

    IF (.NOT.ASSOCIATED(kg_forces)) THEN
      CALL stop_program(routineP,"The kg_forces pointer is not associated "//&
                                "and cannot be initialized")
    END IF

    DO ikind=1,SIZE(kg_forces)
      kg_forces(ikind)%f_rspace_core(:,:) = 0.0_dp
      kg_forces(ikind)%f_hartree_core(:,:) = 0.0_dp
      kg_forces(ikind)%f_ppl(:,:) = 0.0_dp
      kg_forces(ikind)%f_rho(:,:) = 0.0_dp
      kg_forces(ikind)%f_hartree(:,:) = 0.0_dp
      kg_forces(ikind)%f_xc(:,:) = 0.0_dp
      kg_forces(ikind)%f_bc(:,:) = 0.0_dp
      kg_forces(ikind)%f_total(:,:) = 0.0_dp
    END DO

  END SUBROUTINE init_kg_force

END MODULE kg_force_types

