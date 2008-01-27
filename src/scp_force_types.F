!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/scp_force_types [1.0] *
!!
!!   NAME
!!     scp_force_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE scp_force_types

  USE kinds,                           ONLY: dp
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE

  TYPE scp_force_type
    REAL(KIND=dp), DIMENSION(:,:), POINTER :: f_scp, f_total
  END TYPE scp_force_type

! *** Public data types ***

  PUBLIC :: scp_force_type

! *** Public subroutines ***

  PUBLIC :: allocate_scp_force,&
            deallocate_scp_force,&
            init_scp_force

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE allocate_scp_force(scp_forces,natom_of_kind)

!   Purpose: Allocate a scp force data structure.

!   ***************************************************************************
    

    TYPE(scp_force_type), DIMENSION(:), &
      POINTER                                :: scp_forces
    INTEGER, DIMENSION(:), INTENT(IN)        :: natom_of_kind

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE allocate_scp_force (MODULE scp_force_types)"

    INTEGER                                  :: ikind, istat, n, nkind

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(scp_forces)) THEN
      CALL deallocate_scp_force(scp_forces)
    END IF

    nkind = SIZE(natom_of_kind)

    ALLOCATE (scp_forces(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"scp_forces",nkind)

    DO ikind=1,nkind

      n = natom_of_kind(ikind)

      ALLOCATE (scp_forces(ikind)%f_scp (3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"scp_forces(ikind)%f_scp",3*n)
      END IF

      ALLOCATE (scp_forces(ikind)%f_total(3,n),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"scp_forces(ikind)%f_total",3*n)
      END IF

    END DO

  END SUBROUTINE allocate_scp_force 

! *****************************************************************************

  SUBROUTINE deallocate_scp_force(scp_forces)

!   Purpose: Deallocate a scp force data structure.


!   ***************************************************************************

    TYPE(scp_force_type), DIMENSION(:), &
      POINTER                                :: scp_forces

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE deallocate_scp_force (MODULE scp_force_types)"

    INTEGER                                  :: ikind, istat, nkind

!   ---------------------------------------------------------------------------

    IF (.NOT.ASSOCIATED(scp_forces)) THEN
      CALL stop_program(routine,"The scp_forces pointer is not associated "//&
                                "and cannot be deallocated")
    END IF

    nkind = SIZE(scp_forces)

    DO ikind=1,nkind

      IF (ASSOCIATED(scp_forces(ikind)%f_scp)) THEN
        DEALLOCATE (scp_forces(ikind)%f_scp,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routine,"scp_forces%f_scp")
        END IF
      END IF

      IF (ASSOCIATED(scp_forces(ikind)%f_total)) THEN
        DEALLOCATE (scp_forces(ikind)%f_total,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routine,"scp_forces%f_total")
        END IF
      END IF

    END DO

    DEALLOCATE (scp_forces,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"scp_forces")

  END SUBROUTINE deallocate_scp_force

! *****************************************************************************

  SUBROUTINE init_scp_force(scp_forces)

!   Purpose: Initialize a scp force data structure.


!   ***************************************************************************


    TYPE(scp_force_type), DIMENSION(:), &
      POINTER                                :: scp_forces

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE init_scp_force (MODULE scp_force_types)"

    INTEGER                                  :: ikind

!   ---------------------------------------------------------------------------

    IF (.NOT.ASSOCIATED(scp_forces)) THEN
      CALL stop_program(routine,"The scp_forces pointer is not associated "//&
                                "and cannot be initialized")
    END IF

    DO ikind=1,SIZE(scp_forces)
      scp_forces(ikind)%f_scp(:,:) = 0.0_dp
      scp_forces(ikind)%f_total(:,:) = 0.0_dp
    END DO

  END SUBROUTINE init_scp_force

! *****************************************************************************

END MODULE scp_force_types

! *****************************************************************************
