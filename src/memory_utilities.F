!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/memory_utilities [1.0] *
!!
!!   NAME
!!     memory_utilities
!!
!!   FUNCTION
!!     Utility routines for the memory handling.
!!
!!   AUTHOR
!!     Matthias Krack (25.06.1999)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE memory_utilities

! *****************************************************************************

  USE kinds,  ONLY: int_size,&
                    wp => dp,&
                    wp_size => dp_size

  USE termination, ONLY: stop_memory

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: reallocate

! *****************************************************************************

  INTERFACE reallocate
    MODULE PROCEDURE reallocate_i1,reallocate_i2,reallocate_i3,&
                     reallocate_r1,reallocate_r2,reallocate_r3
  END INTERFACE

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE reallocate_i1(p,lb1_new,ub1_new)

!   Purpose: (Re)Allocate an integer vector with a new dimension.

!   History: - Creation (18.07.2002,MK)

!   ***************************************************************************

    INTEGER, DIMENSION(:), POINTER :: p
    INTEGER, INTENT(IN)            :: lb1_new,ub1_new

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE reallocate_i1 (MODULE memory_utilities)"

!   *** Local variables ***

    INTEGER :: istat,lb1,lb1_old,new_size,old_size,ub1,ub1_old

    INTEGER, DIMENSION(:), ALLOCATABLE :: work

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(p)) THEN
      lb1_old = LBOUND(p,1)
      ub1_old = UBOUND(p,1)
      ALLOCATE (work(lb1_old:ub1_old),STAT=istat)
      IF (istat /= 0) THEN
        old_size = (ub1_old - lb1_old + 1)*int_size
        CALL stop_memory(routine,"work",old_size)
      END IF
      work(:) = p(:)
      DEALLOCATE (p,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"p")
    END IF

    ALLOCATE (p(lb1_new:ub1_new),STAT=istat)
    IF (istat /= 0) THEN
      new_size = (ub1_new - lb1_new + 1)*int_size
      CALL stop_memory(routine,"p",new_size)
    END IF
    p(:) = 0

    IF (ALLOCATED(work)) THEN
      lb1 = MAX(lb1_new,lb1_old)
      ub1 = MIN(ub1_new,ub1_old)
      p(lb1:ub1) = work(lb1:ub1)
      DEALLOCATE (work,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"work")
    END IF

  END SUBROUTINE reallocate_i1

! *****************************************************************************

  SUBROUTINE reallocate_i2(p,lb1_new,ub1_new,lb2_new,ub2_new)

!   Purpose: (Re)Allocate an integer array with new dimensions.

!   History: - Creation (18.07.2002,MK)

!   ***************************************************************************

    INTEGER, DIMENSION(:,:), POINTER :: p
    INTEGER, INTENT(IN)              :: lb1_new,lb2_new,ub1_new,ub2_new

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE reallocate_i2 (MODULE memory_utilities)"

!   *** Local variables ***

    INTEGER :: istat,lb1,lb1_old,lb2,lb2_old,new_size,old_size,ub1,ub1_old,&
               ub2,ub2_old

    INTEGER, DIMENSION(:,:), ALLOCATABLE :: work

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(p)) THEN
      lb1_old = LBOUND(p,1)
      ub1_old = UBOUND(p,1)
      lb2_old = LBOUND(p,2)
      ub2_old = UBOUND(p,2)
      ALLOCATE (work(lb1_old:ub1_old,lb2_old:ub2_old),STAT=istat)
      IF (istat /= 0) THEN
        old_size = (ub1_old - lb1_old + 1)*(ub2_old - lb2_old + 1)*int_size
        CALL stop_memory(routine,"work",old_size)
      END IF
      work(:,:) = p(:,:)
      DEALLOCATE (p,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"p")
    END IF

    ALLOCATE (p(lb1_new:ub1_new,lb2_new:ub2_new),STAT=istat)
    IF (istat /= 0) THEN
      new_size = (ub1_new - lb1_new + 1)*(ub2_new - lb2_new + 1)*int_size
      CALL stop_memory(routine,"p",new_size)
    END IF
    p(:,:) = 0

    IF (ALLOCATED(work)) THEN
      lb1 = MAX(lb1_new,lb1_old)
      ub1 = MIN(ub1_new,ub1_old)
      lb2 = MAX(lb2_new,lb2_old)
      ub2 = MIN(ub2_new,ub2_old)
      p(lb1:ub1,lb2:ub2) = work(lb1:ub1,lb2:ub2)
      DEALLOCATE (work,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"work")
    END IF

  END SUBROUTINE reallocate_i2

! *****************************************************************************

  SUBROUTINE reallocate_i3(p,lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new)

!   Purpose: (Re)Allocate an integer array with new dimensions.

!   History: - Creation (18.07.2002,MK)

!   ***************************************************************************

    INTEGER, DIMENSION(:,:,:), POINTER :: p
    INTEGER, INTENT(IN)                :: lb1_new,lb2_new,lb3_new,&
                                          ub1_new,ub2_new,ub3_new

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE reallocate_i3 (MODULE memory_utilities)"

!   *** Local variables ***

    INTEGER :: istat,lb1,lb1_old,lb2,lb2_old,lb3,lb3_old,new_size,old_size,&
               ub1,ub1_old,ub2,ub2_old,ub3,ub3_old

    INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: work

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(p)) THEN
      lb1_old = LBOUND(p,1)
      ub1_old = UBOUND(p,1)
      lb2_old = LBOUND(p,2)
      ub2_old = UBOUND(p,2)
      lb3_old = LBOUND(p,3)
      ub3_old = UBOUND(p,3)
      ALLOCATE (work(lb1_old:ub1_old,&
                     lb2_old:ub2_old,&
                     lb3_old:ub3_old),STAT=istat)
      IF (istat /= 0) THEN
        old_size = (ub1_old - lb1_old + 1)*&
                   (ub2_old - lb2_old + 1)*&
                   (ub3_old - lb3_old + 1)*int_size
        CALL stop_memory(routine,"work",old_size)
      END IF
      work(:,:,:) = p(:,:,:)
      DEALLOCATE (p,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"p")
    END IF

    ALLOCATE (p(lb1_new:ub1_new,lb2_new:ub2_new,lb3_new:ub3_new),STAT=istat)
    IF (istat /= 0) THEN
      new_size = (ub1_new - lb1_new + 1)*&
                 (ub2_new - lb2_new + 1)*&
                 (ub3_new - lb3_new + 1)*int_size
      CALL stop_memory(routine,"p",new_size)
    END IF
    p(:,:,:) = 0

    IF (ALLOCATED(work)) THEN
      lb1 = MAX(lb1_new,lb1_old)
      ub1 = MIN(ub1_new,ub1_old)
      lb2 = MAX(lb2_new,lb2_old)
      ub2 = MIN(ub2_new,ub2_old)
      lb3 = MAX(lb3_new,lb3_old)
      ub3 = MIN(ub3_new,ub3_old)
      p(lb1:ub1,lb2:ub2,lb3:ub3) = work(lb1:ub1,lb2:ub2,lb3:ub3)
      DEALLOCATE (work,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"work")
    END IF

  END SUBROUTINE reallocate_i3

! *****************************************************************************

  SUBROUTINE reallocate_r1(p,lb1_new,ub1_new)

!   Purpose: (Re)Allocate an real vector with a new dimension.

!   History: - Creation (18.07.2002,MK)

!   ***************************************************************************

    REAL(wp), DIMENSION(:), POINTER :: p
    INTEGER, INTENT(IN)             :: lb1_new,ub1_new

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE reallocate_r1 (MODULE memory_utilities)"

!   *** Local variables ***

    INTEGER :: istat,lb1,lb1_old,new_size,old_size,ub1,ub1_old

    REAL(wp), DIMENSION(:), ALLOCATABLE :: work

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(p)) THEN
      lb1_old = LBOUND(p,1)
      ub1_old = UBOUND(p,1)
      ALLOCATE (work(lb1_old:ub1_old),STAT=istat)
      IF (istat /= 0) THEN
        old_size = (ub1_old - lb1_old + 1)*wp_size
        CALL stop_memory(routine,"work",old_size)
      END IF
      work(:) = p(:)
      DEALLOCATE (p,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"p")
    END IF

    ALLOCATE (p(lb1_new:ub1_new),STAT=istat)
    IF (istat /= 0) THEN
      new_size = (ub1_new - lb1_new + 1)*wp_size
      CALL stop_memory(routine,"p",new_size)
    END IF
    p(:) = 0.0_wp

    IF (ALLOCATED(work)) THEN
      lb1 = MAX(lb1_new,lb1_old)
      ub1 = MIN(ub1_new,ub1_old)
      p(lb1:ub1) = work(lb1:ub1)
      DEALLOCATE (work,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"work")
    END IF

  END SUBROUTINE reallocate_r1

! *****************************************************************************

  SUBROUTINE reallocate_r2(p,lb1_new,ub1_new,lb2_new,ub2_new)

!   Purpose: (Re)Allocate an real array with new dimensions.

!   History: - Creation (18.07.2002,MK)

!   ***************************************************************************

    REAL(wp), DIMENSION(:,:), POINTER :: p
    INTEGER, INTENT(IN)               :: lb1_new,lb2_new,ub1_new,ub2_new

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE reallocate_r2 (MODULE memory_utilities)"

!   *** Local variables ***

    INTEGER :: istat,lb1,lb1_old,lb2,lb2_old,new_size,old_size,ub1,ub1_old,&
               ub2,ub2_old

    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: work

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(p)) THEN
      lb1_old = LBOUND(p,1)
      ub1_old = UBOUND(p,1)
      lb2_old = LBOUND(p,2)
      ub2_old = UBOUND(p,2)
      ALLOCATE (work(lb1_old:ub1_old,lb2_old:ub2_old),STAT=istat)
      IF (istat /= 0) THEN
        old_size = (ub1_old - lb1_old + 1)*(ub2_old - lb2_old + 1)*wp_size
        CALL stop_memory(routine,"work",old_size)
      END IF
      work(:,:) = p(:,:)
      DEALLOCATE (p,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"p")
    END IF

    ALLOCATE (p(lb1_new:ub1_new,lb2_new:ub2_new),STAT=istat)
    IF (istat /= 0) THEN
      new_size = (ub1_new - lb1_new + 1)*(ub2_new - lb2_new + 1)*wp_size
      CALL stop_memory(routine,"p",new_size)
    END IF
    p(:,:) = 0.0_wp

    IF (ALLOCATED(work)) THEN
      lb1 = MAX(lb1_new,lb1_old)
      ub1 = MIN(ub1_new,ub1_old)
      lb2 = MAX(lb2_new,lb2_old)
      ub2 = MIN(ub2_new,ub2_old)
      p(lb1:ub1,lb2:ub2) = work(lb1:ub1,lb2:ub2)
      DEALLOCATE (work,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"work")
    END IF

  END SUBROUTINE reallocate_r2

! *****************************************************************************

  SUBROUTINE reallocate_r3(p,lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new)

!   Purpose: (Re)Allocate an real array with new dimensions.

!   History: - Creation (18.07.2002,MK)

!   ***************************************************************************

    REAL(wp), DIMENSION(:,:,:), POINTER :: p
    INTEGER, INTENT(IN)                 :: lb1_new,lb2_new,lb3_new,&
                                           ub1_new,ub2_new,ub3_new

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE reallocate_r3 (MODULE memory_utilities)"

!   *** Local variables ***

    INTEGER :: istat,lb1,lb1_old,lb2,lb2_old,lb3,lb3_old,new_size,old_size,&
               ub1,ub1_old,ub2,ub2_old,ub3,ub3_old

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: work

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(p)) THEN
      lb1_old = LBOUND(p,1)
      ub1_old = UBOUND(p,1)
      lb2_old = LBOUND(p,2)
      ub2_old = UBOUND(p,2)
      lb3_old = LBOUND(p,3)
      ub3_old = UBOUND(p,3)
      ALLOCATE (work(lb1_old:ub1_old,&
                     lb2_old:ub2_old,&
                     lb3_old:ub3_old),STAT=istat)
      IF (istat /= 0) THEN
        old_size = (ub1_old - lb1_old + 1)*&
                   (ub2_old - lb2_old + 1)*&
                   (ub3_old - lb3_old + 1)*wp_size
        CALL stop_memory(routine,"work",old_size)
      END IF
      work(:,:,:) = p(:,:,:)
      DEALLOCATE (p,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"p")
    END IF

    ALLOCATE (p(lb1_new:ub1_new,lb2_new:ub2_new,lb3_new:ub3_new),STAT=istat)
    IF (istat /= 0) THEN
      new_size = (ub1_new - lb1_new + 1)*&
                 (ub2_new - lb2_new + 1)*&
                 (ub3_new - lb3_new + 1)*wp_size
      CALL stop_memory(routine,"p",new_size)
    END IF
    p(:,:,:) = 0.0_wp

    IF (ALLOCATED(work)) THEN
      lb1 = MAX(lb1_new,lb1_old)
      ub1 = MIN(ub1_new,ub1_old)
      lb2 = MAX(lb2_new,lb2_old)
      ub2 = MIN(ub2_new,ub2_old)
      lb3 = MAX(lb3_new,lb3_old)
      ub3 = MIN(ub3_new,ub3_old)
      p(lb1:ub1,lb2:ub2,lb3:ub3) = work(lb1:ub1,lb2:ub2,lb3:ub3)
      DEALLOCATE (work,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"work")
    END IF

  END SUBROUTINE reallocate_r3

! *****************************************************************************

END MODULE memory_utilities
