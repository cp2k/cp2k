!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/kg_force_types [1.0] *
!!
!!   NAME
!!     kg_force_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     gloria (01.10.2002)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE kg_force_types

  USE kinds,                           ONLY: dp
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE

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

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE allocate_kg_force(kg_forces,natom_of_kind)

!   Purpose: Allocate a kg force data structure.

!   ***************************************************************************
    

    TYPE(kg_force_type), DIMENSION(:), &
      POINTER                                :: kg_forces
    INTEGER, DIMENSION(:), INTENT(IN)        :: natom_of_kind

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE allocate_kg_force (MODULE kg_force_types)"

    INTEGER                                  :: ikind, istat, n, nkind

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(kg_forces)) THEN
      CALL deallocate_kg_force(kg_forces)
    END IF

    nkind = SIZE(natom_of_kind)

    ALLOCATE (kg_forces(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"kg_forces",nkind)

    DO ikind=1,nkind

      n = natom_of_kind(ikind)

      ALLOCATE (kg_forces(ikind)%f_rspace_core(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"kg_forces(ikind)%f_rspace_core",3*n)
      END IF

      ALLOCATE (kg_forces(ikind)%f_hartree_core(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"kg_forces(ikind)%f_hartree_core",3*n)
      END IF

      ALLOCATE (kg_forces(ikind)%f_rho(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"kg_forces(ikind)%f_rho",3*n)
      END IF

      ALLOCATE (kg_forces(ikind)%f_hartree(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"kg_forces(ikind)%f_hartree",3*n)
      END IF

      ALLOCATE (kg_forces(ikind)%f_xc(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"kg_forces(ikind)%f_xc",3*n)
      END IF

      ALLOCATE (kg_forces(ikind)%f_bc(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"kg_forces(ikind)%f_bc",3*n)
      END IF

      ALLOCATE (kg_forces(ikind)%f_ppl(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"kg_forces(ikind)%f_ppl",3*n)
      END IF

      ALLOCATE (kg_forces(ikind)%f_total(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"kg_forces(ikind)%f_total",3*n)
      END IF

    END DO

  END SUBROUTINE allocate_kg_force 

! *****************************************************************************

  SUBROUTINE deallocate_kg_force(kg_forces)

!   Purpose: Deallocate a kg force data structure.

!   History: - Creation (30.09.2002,gt)

!   ***************************************************************************

    TYPE(kg_force_type), DIMENSION(:), &
      POINTER                                :: kg_forces

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE deallocate_kg_force (MODULE kg_force_types)"

    INTEGER                                  :: ikind, istat, nkind

!   ---------------------------------------------------------------------------

    IF (.NOT.ASSOCIATED(kg_forces)) THEN
      CALL stop_program(routine,"The kg_forces pointer is not associated "//&
                                "and cannot be deallocated")
    END IF

    nkind = SIZE(kg_forces)

    DO ikind=1,nkind

      IF (ASSOCIATED(kg_forces(ikind)%f_rspace_core)) THEN
        DEALLOCATE (kg_forces(ikind)%f_rspace_core,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routine,"kg_forces%f_rspace_core")
        END IF
      END IF

      IF (ASSOCIATED(kg_forces(ikind)%f_hartree_core)) THEN
        DEALLOCATE (kg_forces(ikind)%f_hartree_core,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routine,"kg_forces%f_hartree_core")
        END IF
      END IF

      IF (ASSOCIATED(kg_forces(ikind)%f_rho)) THEN
        DEALLOCATE (kg_forces(ikind)%f_rho,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routine,"kg_forces%f_rho")
        END IF
      END IF

      IF (ASSOCIATED(kg_forces(ikind)%f_ppl)) THEN
        DEALLOCATE (kg_forces(ikind)%f_ppl,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routine,"kg_forces%f_ppl")
        END IF
      END IF

      IF (ASSOCIATED(kg_forces(ikind)%f_hartree)) THEN
        DEALLOCATE (kg_forces(ikind)%f_hartree,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routine,"kg_forces%f_hartree")
        END IF
      END IF

      IF (ASSOCIATED(kg_forces(ikind)%f_xc)) THEN
        DEALLOCATE (kg_forces(ikind)%f_xc,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routine,"kg_forces%f_xc")
        END IF
      END IF

      IF (ASSOCIATED(kg_forces(ikind)%f_bc)) THEN
        DEALLOCATE (kg_forces(ikind)%f_bc,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routine,"kg_forces%f_bc")
        END IF
      END IF

      IF (ASSOCIATED(kg_forces(ikind)%f_total)) THEN
        DEALLOCATE (kg_forces(ikind)%f_total,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routine,"kg_forces%f_total")
        END IF
      END IF

    END DO

    DEALLOCATE (kg_forces,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"kg_forces")

  END SUBROUTINE deallocate_kg_force

! *****************************************************************************

  SUBROUTINE init_kg_force(kg_forces)

!   Purpose: Initialize a kg force data structure.

!   History: - Creation (01.10.2002,gt)

!   ***************************************************************************


    TYPE(kg_force_type), DIMENSION(:), &
      POINTER                                :: kg_forces

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE init_kg_force (MODULE kg_force_types)"

    INTEGER                                  :: ikind

!   ---------------------------------------------------------------------------

    IF (.NOT.ASSOCIATED(kg_forces)) THEN
      CALL stop_program(routine,"The kg_forces pointer is not associated "//&
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

! *****************************************************************************

END MODULE kg_force_types

! *****************************************************************************
