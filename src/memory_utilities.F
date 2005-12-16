!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2002  CP2K developers group
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

  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             dp_size,&
                                             int_size
  USE machine,                         ONLY: print_memory
  USE message_passing,                 ONLY: mp_sync
  USE termination,                     ONLY: stop_memory

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "memory_utilities"

  PUBLIC :: reallocate,&
            write_memory

! *****************************************************************************

  INTERFACE reallocate
    MODULE PROCEDURE reallocate_c1,reallocate_c2,reallocate_c3,reallocate_c4,&
                     reallocate_i1,reallocate_i2,reallocate_i3,reallocate_i4,&
                     reallocate_r1,reallocate_r2,reallocate_r3,reallocate_r4,&
                     reallocate_r5,reallocate_s1
  END INTERFACE

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE reallocate_c1(p,lb1_new,ub1_new)

    ! (Re)Allocate a complex vector with a new dimension (28.11.2005,MK)

    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "reallocate_c1"

    INTEGER                                  :: istat, lb1, lb1_old, &
                                                new_size, old_size, ub1, &
                                                ub1_old

    COMPLEX(KIND=dp), ALLOCATABLE,&
      DIMENSION(:)                           :: work

    ! -------------------------------------------------------------------------

    IF (ASSOCIATED(p)) THEN
      lb1_old = LBOUND(p,1)
      ub1_old = UBOUND(p,1)
      lb1 = MAX(lb1_new,lb1_old)
      ub1 = MIN(ub1_new,ub1_old)
      ALLOCATE (work(lb1:ub1),STAT=istat)
      IF (istat /= 0) THEN
        old_size = (ub1_old - lb1_old + 1)*dp_size
        CALL stop_memory(routine_name,module_name,__LINE__,"work",old_size)
      END IF
      work(lb1:ub1) = p(lb1:ub1)
      DEALLOCATE (p,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"p")
    END IF

    ALLOCATE (p(lb1_new:ub1_new),STAT=istat)
    IF (istat /= 0) THEN
      new_size = (ub1_new - lb1_new + 1)*dp_size
      CALL stop_memory(routine_name,module_name,__LINE__,"p",new_size)
    END IF
    p(:) = (0.0_dp,0.0_dp)

    IF (ALLOCATED(work)) THEN
      lb1 = MAX(lb1_new,lb1_old)
      ub1 = MIN(ub1_new,ub1_old)
      p(lb1:ub1) = work(lb1:ub1)
      DEALLOCATE (work,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine_name,module_name,__LINE__,"work")
      END IF
    END IF

  END SUBROUTINE reallocate_c1

  ! ***************************************************************************

  SUBROUTINE reallocate_c2(p,lb1_new,ub1_new,lb2_new,ub2_new)

    ! (Re)Allocate a complex array with new dimensions (28.11.2005,MK)

    COMPLEX(KIND=dp), DIMENSION(:, :),&
      POINTER                                :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new, lb2_new, &
                                                ub2_new

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "reallocate_c2"

    INTEGER                                  :: istat, lb1, lb1_old, lb2, &
                                                lb2_old, new_size, old_size, &
                                                ub1, ub1_old, ub2, ub2_old
    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: work

    ! -------------------------------------------------------------------------

    IF (ASSOCIATED(p)) THEN
      lb1_old = LBOUND(p,1)
      ub1_old = UBOUND(p,1)
      lb2_old = LBOUND(p,2)
      ub2_old = UBOUND(p,2)
      ALLOCATE (work(lb1_old:ub1_old,lb2_old:ub2_old),STAT=istat)
      IF (istat /= 0) THEN
        old_size = (ub1_old - lb1_old + 1)*(ub2_old - lb2_old + 1)*dp_size
        CALL stop_memory(routine_name,module_name,__LINE__,"work",old_size)
      END IF
      work(:,:) = p(:,:)
      DEALLOCATE (p,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"p")
    END IF

    ALLOCATE (p(lb1_new:ub1_new,lb2_new:ub2_new),STAT=istat)
    IF (istat /= 0) THEN
      new_size = (ub1_new - lb1_new + 1)*(ub2_new - lb2_new + 1)*dp_size
      CALL stop_memory(routine_name,module_name,__LINE__,"p",new_size)
    END IF
    p(:,:) = (0.0_dp,0.0_dp)

    IF (ALLOCATED(work)) THEN
      lb1 = MAX(lb1_new,lb1_old)
      ub1 = MIN(ub1_new,ub1_old)
      lb2 = MAX(lb2_new,lb2_old)
      ub2 = MIN(ub2_new,ub2_old)
      p(lb1:ub1,lb2:ub2) = work(lb1:ub1,lb2:ub2)
      DEALLOCATE (work,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine_name,module_name,__LINE__,"work")
      END IF
    END IF

  END SUBROUTINE reallocate_c2

  ! ***************************************************************************

  SUBROUTINE reallocate_c3(p,lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new)

    ! (Re)Allocate an complex array with new dimensions (28.11.2005,MK)

    COMPLEX(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new, lb2_new, &
                                                ub2_new, lb3_new, ub3_new

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "reallocate_c3"

    INTEGER :: istat, lb1, lb1_old, lb2, lb2_old, lb3, lb3_old, new_size, &
      old_size, ub1, ub1_old, ub2, ub2_old, ub3, ub3_old
    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: work

    ! -------------------------------------------------------------------------

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
                   (ub3_old - lb3_old + 1)*dp_size
        CALL stop_memory(routine_name,module_name,__LINE__,"work",old_size)
      END IF
      work(:,:,:) = p(:,:,:)
      DEALLOCATE (p,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"p")
    END IF

    ALLOCATE (p(lb1_new:ub1_new,lb2_new:ub2_new,lb3_new:ub3_new),STAT=istat)
    IF (istat /= 0) THEN
      new_size = (ub1_new - lb1_new + 1)*&
                 (ub2_new - lb2_new + 1)*&
                 (ub3_new - lb3_new + 1)*dp_size
      CALL stop_memory(routine_name,module_name,__LINE__,"p",new_size)
    END IF
    p(:,:,:) = (0.0_dp,0.0_dp)

    IF (ALLOCATED(work)) THEN
      lb1 = MAX(lb1_new,lb1_old)
      ub1 = MIN(ub1_new,ub1_old)
      lb2 = MAX(lb2_new,lb2_old)
      ub2 = MIN(ub2_new,ub2_old)
      lb3 = MAX(lb3_new,lb3_old)
      ub3 = MIN(ub3_new,ub3_old)
      p(lb1:ub1,lb2:ub2,lb3:ub3) = work(lb1:ub1,lb2:ub2,lb3:ub3)
      DEALLOCATE (work,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine_name,module_name,__LINE__,"work")
      END IF
    END IF

  END SUBROUTINE reallocate_c3

  ! ***************************************************************************

  SUBROUTINE reallocate_c4(p,lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new,&
                             lb4_new,ub4_new)

    ! (Re)Allocate an complex array with new dimensions (28.11.2005,MK)

    COMPLEX(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new, lb2_new, &
                                                ub2_new, lb3_new, ub3_new, &
                                                lb4_new, ub4_new

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "reallocate_c4"

    INTEGER :: istat, lb1, lb1_old, lb2, lb2_old, lb3, lb3_old, lb4, lb4_old, &
      new_size, old_size, ub1, ub1_old, ub2, ub2_old, ub3, ub3_old, ub4, &
      ub4_old
    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: work

    ! -------------------------------------------------------------------------

    IF (ASSOCIATED(p)) THEN
      lb1_old = LBOUND(p,1)
      ub1_old = UBOUND(p,1)
      lb2_old = LBOUND(p,2)
      ub2_old = UBOUND(p,2)
      lb3_old = LBOUND(p,3)
      ub3_old = UBOUND(p,3)
      lb4_old = LBOUND(p,4)
      ub4_old = UBOUND(p,4)
      ALLOCATE (work(lb1_old:ub1_old,&
                     lb2_old:ub2_old,&
                     lb3_old:ub3_old,&
                     lb4_old:ub4_old),STAT=istat)
      IF (istat /= 0) THEN
        old_size = (ub1_old - lb1_old + 1)*&
                   (ub2_old - lb2_old + 1)*&
                   (ub3_old - lb3_old + 1)*&
                   (ub4_old - lb4_old + 1)*dp_size
        CALL stop_memory(routine_name,module_name,__LINE__,"work",old_size)
      END IF
      work(:,:,:,:) = p(:,:,:,:)
      DEALLOCATE (p,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"p")
    END IF

    ALLOCATE (p(lb1_new:ub1_new,lb2_new:ub2_new,lb3_new:ub3_new,&
                lb4_new:ub4_new),STAT=istat)
    IF (istat /= 0) THEN
      new_size = (ub1_new - lb1_new + 1)*&
                 (ub2_new - lb2_new + 1)*&
                 (ub3_new - lb3_new + 1)*&
                 (ub4_new - lb4_new + 1)*dp_size
      CALL stop_memory(routine_name,module_name,__LINE__,"p",new_size)
    END IF
    p(:,:,:,:) = (0.0_dp,0.0_dp)

    IF (ALLOCATED(work)) THEN
      lb1 = MAX(lb1_new,lb1_old)
      ub1 = MIN(ub1_new,ub1_old)
      lb2 = MAX(lb2_new,lb2_old)
      ub2 = MIN(ub2_new,ub2_old)
      lb3 = MAX(lb3_new,lb3_old)
      ub3 = MIN(ub3_new,ub3_old)
      lb4 = MAX(lb4_new,lb4_old)
      ub4 = MIN(ub4_new,ub4_old)
      p(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4) = work(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4)
      DEALLOCATE (work,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine_name,module_name,__LINE__,"work")
      END IF
    END IF

  END SUBROUTINE reallocate_c4

! *****************************************************************************

  SUBROUTINE reallocate_i1(p,lb1_new,ub1_new)

!   Purpose: (Re)Allocate an integer vector with a new dimension.

!   History: - Creation (18.07.2002,MK)

!   ***************************************************************************

    INTEGER, DIMENSION(:), POINTER           :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "reallocate_i1"

    INTEGER                                  :: istat, lb1, lb1_old, &
                                                new_size, old_size, ub1, &
                                                ub1_old
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: work

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(p)) THEN
      lb1_old = LBOUND(p,1)
      ub1_old = UBOUND(p,1)
      lb1 = MAX(lb1_new,lb1_old)
      ub1 = MIN(ub1_new,ub1_old)
      ALLOCATE (work(lb1:ub1),STAT=istat)
      IF (istat /= 0) THEN
        old_size = (ub1 - lb1 + 1)*int_size
        CALL stop_memory(routine_name,module_name,__LINE__,"work",old_size)
      END IF
      work(lb1:ub1) = p(lb1:ub1)
      DEALLOCATE (p,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"p")
    END IF

    ALLOCATE (p(lb1_new:ub1_new),STAT=istat)
    IF (istat /= 0) THEN
      new_size = (ub1_new - lb1_new + 1)*int_size
      CALL stop_memory(routine_name,module_name,__LINE__,"p",new_size)
    END IF
    p(:) = 0

    IF (ALLOCATED(work)) THEN
      lb1 = MAX(lb1_new,lb1_old)
      ub1 = MIN(ub1_new,ub1_old)
      p(lb1:ub1) = work(lb1:ub1)
      DEALLOCATE (work,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine_name,module_name,__LINE__,"work")
      END IF
    END IF

  END SUBROUTINE reallocate_i1

! *****************************************************************************

  SUBROUTINE reallocate_i2(p,lb1_new,ub1_new,lb2_new,ub2_new)

!   Purpose: (Re)Allocate an integer array with new dimensions.

!   History: - Creation (18.07.2002,MK)

!   ***************************************************************************

    INTEGER, DIMENSION(:, :), POINTER        :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new, lb2_new, &
                                                ub2_new

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "reallocate_i2"

    INTEGER                                  :: istat, lb1, lb1_old, lb2, &
                                                lb2_old, new_size, old_size, &
                                                ub1, ub1_old, ub2, ub2_old
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: work

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(p)) THEN
      lb1_old = LBOUND(p,1)
      ub1_old = UBOUND(p,1)
      lb2_old = LBOUND(p,2)
      ub2_old = UBOUND(p,2)
      lb1 = MAX(lb1_new,lb1_old)
      ub1 = MIN(ub1_new,ub1_old)
      lb2 = MAX(lb2_new,lb2_old)
      ub2 = MIN(ub2_new,ub2_old)
      ALLOCATE (work(lb1:ub1,lb2:ub2),STAT=istat)
      IF (istat /= 0) THEN
        old_size = (ub1 - lb1 + 1)*(ub2 - lb2 + 1)*int_size
        CALL stop_memory(routine_name,module_name,__LINE__,"work",old_size)
      END IF
      work(lb1:ub1,lb2:ub2) = p(lb1:ub1,lb2:ub2)
      DEALLOCATE (p,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"p")
    END IF

    ALLOCATE (p(lb1_new:ub1_new,lb2_new:ub2_new),STAT=istat)
    IF (istat /= 0) THEN
      new_size = (ub1_new - lb1_new + 1)*(ub2_new - lb2_new + 1)*int_size
      CALL stop_memory(routine_name,module_name,__LINE__,"p",new_size)
    END IF
    p(:,:) = 0

    IF (ALLOCATED(work)) THEN
      p(lb1:ub1,lb2:ub2) = work(lb1:ub1,lb2:ub2)
      DEALLOCATE (work,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine_name,module_name,__LINE__,"work")
      END IF
    END IF

  END SUBROUTINE reallocate_i2

! *****************************************************************************

  SUBROUTINE reallocate_i3(p,lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new)

!   Purpose: (Re)Allocate an integer array with new dimensions.

!   History: - Creation (18.07.2002,MK)

!   ***************************************************************************

    INTEGER, DIMENSION(:, :, :), POINTER     :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new, lb2_new, &
                                                ub2_new, lb3_new, ub3_new

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "reallocate_i3"

    INTEGER :: istat, lb1, lb1_old, lb2, lb2_old, lb3, lb3_old, new_size, &
      old_size, ub1, ub1_old, ub2, ub2_old, ub3, ub3_old
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: work

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
        CALL stop_memory(routine_name,module_name,__LINE__,"work",old_size)
      END IF
      work(:,:,:) = p(:,:,:)
      DEALLOCATE (p,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"p")
    END IF

    ALLOCATE (p(lb1_new:ub1_new,lb2_new:ub2_new,lb3_new:ub3_new),STAT=istat)
    IF (istat /= 0) THEN
      new_size = (ub1_new - lb1_new + 1)*&
                 (ub2_new - lb2_new + 1)*&
                 (ub3_new - lb3_new + 1)*int_size
      CALL stop_memory(routine_name,module_name,__LINE__,"p",new_size)
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
      IF (istat /= 0) THEN
        CALL stop_memory(routine_name,module_name,__LINE__,"work")
      END IF
    END IF

  END SUBROUTINE reallocate_i3

! *****************************************************************************

  SUBROUTINE reallocate_i4(p,lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new,&
                             lb4_new,ub4_new)

!   Purpose: (Re)Allocate an integer array with new dimensions.

!   History: - Creation (04.10.2002,MK)

!   ***************************************************************************

    INTEGER, DIMENSION(:, :, :, :), POINTER  :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new, lb2_new, &
                                                ub2_new, lb3_new, ub3_new, &
                                                lb4_new, ub4_new

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "reallocate_i4"

    INTEGER :: istat, lb1, lb1_old, lb2, lb2_old, lb3, lb3_old, lb4, lb4_old, &
      new_size, old_size, ub1, ub1_old, ub2, ub2_old, ub3, ub3_old, ub4, &
      ub4_old
    INTEGER, ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: work

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(p)) THEN
      lb1_old = LBOUND(p,1)
      ub1_old = UBOUND(p,1)
      lb2_old = LBOUND(p,2)
      ub2_old = UBOUND(p,2)
      lb3_old = LBOUND(p,3)
      ub3_old = UBOUND(p,3)
      lb4_old = LBOUND(p,4)
      ub4_old = UBOUND(p,4)
      ALLOCATE (work(lb1_old:ub1_old,&
                     lb2_old:ub2_old,&
                     lb3_old:ub3_old,&
                     lb4_old:ub4_old),STAT=istat)
      IF (istat /= 0) THEN
        old_size = (ub1_old - lb1_old + 1)*&
                   (ub2_old - lb2_old + 1)*&
                   (ub3_old - lb3_old + 1)*&
                   (ub4_old - lb4_old + 1)*int_size
        CALL stop_memory(routine_name,module_name,__LINE__,"work",old_size)
      END IF
      work(:,:,:,:) = p(:,:,:,:)
      DEALLOCATE (p,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"p")
    END IF

    ALLOCATE (p(lb1_new:ub1_new,lb2_new:ub2_new,lb3_new:ub3_new,&
                lb4_new:ub4_new),STAT=istat)
    IF (istat /= 0) THEN
      new_size = (ub1_new - lb1_new + 1)*&
                 (ub2_new - lb2_new + 1)*&
                 (ub3_new - lb3_new + 1)*&
                 (ub4_new - lb4_new + 1)*int_size
      CALL stop_memory(routine_name,module_name,__LINE__,"p",new_size)
    END IF
    p(:,:,:,:) = 0

    IF (ALLOCATED(work)) THEN
      lb1 = MAX(lb1_new,lb1_old)
      ub1 = MIN(ub1_new,ub1_old)
      lb2 = MAX(lb2_new,lb2_old)
      ub2 = MIN(ub2_new,ub2_old)
      lb3 = MAX(lb3_new,lb3_old)
      ub3 = MIN(ub3_new,ub3_old)
      lb4 = MAX(lb4_new,lb4_old)
      ub4 = MIN(ub4_new,ub4_old)
      p(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4) = work(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4)
      DEALLOCATE (work,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine_name,module_name,__LINE__,"work")
      END IF
    END IF

  END SUBROUTINE reallocate_i4

! *****************************************************************************

  SUBROUTINE reallocate_r1(p,lb1_new,ub1_new)

!   Purpose: (Re)Allocate an real vector with a new dimension.

!   History: - Creation (18.07.2002,MK)

!   ***************************************************************************

    REAL(KIND=dp), DIMENSION(:), POINTER     :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "reallocate_r1"

    INTEGER                                  :: istat, lb1, lb1_old, &
                                                new_size, old_size, ub1, &
                                                ub1_old
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: work

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(p)) THEN
      lb1_old = LBOUND(p,1)
      ub1_old = UBOUND(p,1)
      lb1 = MAX(lb1_new,lb1_old)
      ub1 = MIN(ub1_new,ub1_old)
      ALLOCATE (work(lb1:ub1),STAT=istat)
      IF (istat /= 0) THEN
        old_size = (ub1_old - lb1_old + 1)*dp_size
        CALL stop_memory(routine_name,module_name,__LINE__,"work",old_size)
      END IF
      work(lb1:ub1) = p(lb1:ub1)
      DEALLOCATE (p,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"p")
    END IF

    ALLOCATE (p(lb1_new:ub1_new),STAT=istat)
    IF (istat /= 0) THEN
      new_size = (ub1_new - lb1_new + 1)*dp_size
      CALL stop_memory(routine_name,module_name,__LINE__,"p",new_size)
    END IF
    p(:) = 0.0_dp

    IF (ALLOCATED(work)) THEN
      lb1 = MAX(lb1_new,lb1_old)
      ub1 = MIN(ub1_new,ub1_old)
      p(lb1:ub1) = work(lb1:ub1)
      DEALLOCATE (work,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine_name,module_name,__LINE__,"work")
      END IF
    END IF

  END SUBROUTINE reallocate_r1

! *****************************************************************************

  SUBROUTINE reallocate_r2(p,lb1_new,ub1_new,lb2_new,ub2_new)

!   Purpose: (Re)Allocate an real array with new dimensions.

!   History: - Creation (18.07.2002,MK)

!   ***************************************************************************

    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new, lb2_new, &
                                                ub2_new

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "reallocate_r2"

    INTEGER                                  :: istat, lb1, lb1_old, lb2, &
                                                lb2_old, new_size, old_size, &
                                                ub1, ub1_old, ub2, ub2_old
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: work

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(p)) THEN
      lb1_old = LBOUND(p,1)
      ub1_old = UBOUND(p,1)
      lb2_old = LBOUND(p,2)
      ub2_old = UBOUND(p,2)
      ALLOCATE (work(lb1_old:ub1_old,lb2_old:ub2_old),STAT=istat)
      IF (istat /= 0) THEN
        old_size = (ub1_old - lb1_old + 1)*(ub2_old - lb2_old + 1)*dp_size
        CALL stop_memory(routine_name,module_name,__LINE__,"work",old_size)
      END IF
      work(:,:) = p(:,:)
      DEALLOCATE (p,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"p")
    END IF

    ALLOCATE (p(lb1_new:ub1_new,lb2_new:ub2_new),STAT=istat)
    IF (istat /= 0) THEN
      new_size = (ub1_new - lb1_new + 1)*(ub2_new - lb2_new + 1)*dp_size
      CALL stop_memory(routine_name,module_name,__LINE__,"p",new_size)
    END IF
    p(:,:) = 0.0_dp

    IF (ALLOCATED(work)) THEN
      lb1 = MAX(lb1_new,lb1_old)
      ub1 = MIN(ub1_new,ub1_old)
      lb2 = MAX(lb2_new,lb2_old)
      ub2 = MIN(ub2_new,ub2_old)
      p(lb1:ub1,lb2:ub2) = work(lb1:ub1,lb2:ub2)
      DEALLOCATE (work,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine_name,module_name,__LINE__,"work")
      END IF
    END IF

  END SUBROUTINE reallocate_r2

! *****************************************************************************

  SUBROUTINE reallocate_r3(p,lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new)

!   Purpose: (Re)Allocate an real array with new dimensions.

!   History: - Creation (18.07.2002,MK)

!   ***************************************************************************

    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new, lb2_new, &
                                                ub2_new, lb3_new, ub3_new

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "reallocate_r3"

    INTEGER :: istat, lb1, lb1_old, lb2, lb2_old, lb3, lb3_old, new_size, &
      old_size, ub1, ub1_old, ub2, ub2_old, ub3, ub3_old
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: work

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
                   (ub3_old - lb3_old + 1)*dp_size
        CALL stop_memory(routine_name,module_name,__LINE__,"work",old_size)
      END IF
      work(:,:,:) = p(:,:,:)
      DEALLOCATE (p,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"p")
    END IF

    ALLOCATE (p(lb1_new:ub1_new,lb2_new:ub2_new,lb3_new:ub3_new),STAT=istat)
    IF (istat /= 0) THEN
      new_size = (ub1_new - lb1_new + 1)*&
                 (ub2_new - lb2_new + 1)*&
                 (ub3_new - lb3_new + 1)*dp_size
      CALL stop_memory(routine_name,module_name,__LINE__,"p",new_size)
    END IF
    p(:,:,:) = 0.0_dp

    IF (ALLOCATED(work)) THEN
      lb1 = MAX(lb1_new,lb1_old)
      ub1 = MIN(ub1_new,ub1_old)
      lb2 = MAX(lb2_new,lb2_old)
      ub2 = MIN(ub2_new,ub2_old)
      lb3 = MAX(lb3_new,lb3_old)
      ub3 = MIN(ub3_new,ub3_old)
      p(lb1:ub1,lb2:ub2,lb3:ub3) = work(lb1:ub1,lb2:ub2,lb3:ub3)
      DEALLOCATE (work,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine_name,module_name,__LINE__,"work")
      END IF
    END IF

  END SUBROUTINE reallocate_r3

! *****************************************************************************

  SUBROUTINE reallocate_r4(p,lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new,&
                             lb4_new,ub4_new)

!   Purpose: (Re)Allocate an real array with new dimensions.

!   History: - Creation (04.10.2002,MK)

!   ***************************************************************************

    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new, lb2_new, &
                                                ub2_new, lb3_new, ub3_new, &
                                                lb4_new, ub4_new

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "reallocate_r4"

    INTEGER :: istat, lb1, lb1_old, lb2, lb2_old, lb3, lb3_old, lb4, lb4_old, &
      new_size, old_size, ub1, ub1_old, ub2, ub2_old, ub3, ub3_old, ub4, &
      ub4_old
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: work

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(p)) THEN
      lb1_old = LBOUND(p,1)
      ub1_old = UBOUND(p,1)
      lb2_old = LBOUND(p,2)
      ub2_old = UBOUND(p,2)
      lb3_old = LBOUND(p,3)
      ub3_old = UBOUND(p,3)
      lb4_old = LBOUND(p,4)
      ub4_old = UBOUND(p,4)
      ALLOCATE (work(lb1_old:ub1_old,&
                     lb2_old:ub2_old,&
                     lb3_old:ub3_old,&
                     lb4_old:ub4_old),STAT=istat)
      IF (istat /= 0) THEN
        old_size = (ub1_old - lb1_old + 1)*&
                   (ub2_old - lb2_old + 1)*&
                   (ub3_old - lb3_old + 1)*&
                   (ub4_old - lb4_old + 1)*dp_size
        CALL stop_memory(routine_name,module_name,__LINE__,"work",old_size)
      END IF
      work(:,:,:,:) = p(:,:,:,:)
      DEALLOCATE (p,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"p")
    END IF

    ALLOCATE (p(lb1_new:ub1_new,lb2_new:ub2_new,lb3_new:ub3_new,&
                lb4_new:ub4_new),STAT=istat)
    IF (istat /= 0) THEN
      new_size = (ub1_new - lb1_new + 1)*&
                 (ub2_new - lb2_new + 1)*&
                 (ub3_new - lb3_new + 1)*&
                 (ub4_new - lb4_new + 1)*dp_size
      CALL stop_memory(routine_name,module_name,__LINE__,"p",new_size)
    END IF
    p(:,:,:,:) = 0.0_dp

    IF (ALLOCATED(work)) THEN
      lb1 = MAX(lb1_new,lb1_old)
      ub1 = MIN(ub1_new,ub1_old)
      lb2 = MAX(lb2_new,lb2_old)
      ub2 = MIN(ub2_new,ub2_old)
      lb3 = MAX(lb3_new,lb3_old)
      ub3 = MIN(ub3_new,ub3_old)
      lb4 = MAX(lb4_new,lb4_old)
      ub4 = MIN(ub4_new,ub4_old)
      p(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4) = work(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4)
      DEALLOCATE (work,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine_name,module_name,__LINE__,"work")
      END IF
    END IF

  END SUBROUTINE reallocate_r4

! *****************************************************************************

  SUBROUTINE reallocate_r5(p,lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new,&
                             lb4_new,ub4_new,lb5_new,ub5_new)

!   Purpose: (Re)Allocate an real array with new dimensions.

!   History: - Creation (04.10.2002,MK)

!   ***************************************************************************

    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :), POINTER      :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new, lb2_new, &
                                                ub2_new, lb3_new, ub3_new, &
                                                lb4_new, ub4_new, lb5_new, &
                                                ub5_new

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "reallocate_r5"

    INTEGER :: istat, lb1, lb1_old, lb2, lb2_old, lb3, lb3_old, lb4, lb4_old, &
      lb5, lb5_old, new_size, old_size, ub1, ub1_old, ub2, ub2_old, ub3, &
      ub3_old, ub4, ub4_old, ub5, ub5_old
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :, :)               :: work

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(p)) THEN
      lb1_old = LBOUND(p,1)
      ub1_old = UBOUND(p,1)
      lb2_old = LBOUND(p,2)
      ub2_old = UBOUND(p,2)
      lb3_old = LBOUND(p,3)
      ub3_old = UBOUND(p,3)
      lb4_old = LBOUND(p,4)
      ub4_old = UBOUND(p,4)
      lb5_old = LBOUND(p,5)
      ub5_old = UBOUND(p,5)
      ALLOCATE (work(lb1_old:ub1_old,&
                     lb2_old:ub2_old,&
                     lb3_old:ub3_old,&
                     lb4_old:ub4_old,&
                     lb5_old:ub5_old),STAT=istat)
      IF (istat /= 0) THEN
        old_size = (ub1_old - lb1_old + 1)*&
                   (ub2_old - lb2_old + 1)*&
                   (ub3_old - lb3_old + 1)*&
                   (ub4_old - lb4_old + 1)*&
                   (ub5_old - lb5_old + 1)*dp_size
        CALL stop_memory(routine_name,module_name,__LINE__,"work",old_size)
      END IF
      work(:,:,:,:,:) = p(:,:,:,:,:)
      DEALLOCATE (p,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"p")
    END IF

    ALLOCATE (p(lb1_new:ub1_new,lb2_new:ub2_new,lb3_new:ub3_new,&
                lb4_new:ub4_new,lb5_new:ub5_new),STAT=istat)
    IF (istat /= 0) THEN
      new_size = (ub1_new - lb1_new + 1)*&
                 (ub2_new - lb2_new + 1)*&
                 (ub3_new - lb3_new + 1)*&
                 (ub4_new - lb4_new + 1)*&
                 (ub5_new - lb5_new + 1)*dp_size
      CALL stop_memory(routine_name,module_name,__LINE__,"p",new_size)
    END IF
    p(:,:,:,:,:) = 0.0_dp

    IF (ALLOCATED(work)) THEN
      lb1 = MAX(lb1_new,lb1_old)
      ub1 = MIN(ub1_new,ub1_old)
      lb2 = MAX(lb2_new,lb2_old)
      ub2 = MIN(ub2_new,ub2_old)
      lb3 = MAX(lb3_new,lb3_old)
      ub3 = MIN(ub3_new,ub3_old)
      lb4 = MAX(lb4_new,lb4_old)
      ub4 = MIN(ub4_new,ub4_old)
      lb5 = MAX(lb5_new,lb5_old)
      ub5 = MIN(ub5_new,ub5_old)
      p(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4,lb5:ub5) =&
        work(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4,lb5:ub5)
      DEALLOCATE (work,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine_name,module_name,__LINE__,"work")
      END IF
    END IF

  END SUBROUTINE reallocate_r5

! *****************************************************************************

  SUBROUTINE reallocate_s1(p,lb1_new,ub1_new)

!   Purpose: (Re)Allocate a vector of string variables of default string
!            length.

!TC This routine doesnt work on SUN/Solaris!!! It should probably not be used.

!   ***************************************************************************

    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: p
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "reallocate_s1"

    CHARACTER(LEN=default_string_length), &
      ALLOCATABLE, DIMENSION(:)              :: work
    INTEGER                                  :: istat, lb1, lb1_old, &
                                                new_size, old_size, ub1, &
                                                ub1_old

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(p)) THEN
      lb1_old = LBOUND(p,1)
      ub1_old = UBOUND(p,1)
      ALLOCATE (work(lb1_old:ub1_old),STAT=istat)
      IF (istat /= 0) THEN
        old_size = (ub1_old - lb1_old + 1)*1
        CALL stop_memory(routine_name,module_name,__LINE__,"work",old_size)
      END IF
      work(:) = p(:)
      DEALLOCATE (p,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"p")
    END IF

    ALLOCATE (p(lb1_new:ub1_new),STAT=istat)
    IF (istat /= 0) THEN
      new_size = (ub1_new - lb1_new + 1)*1
      CALL stop_memory(routine_name,module_name,__LINE__,"p",new_size)
    END IF
    p(:) = ""

    IF (ALLOCATED(work)) THEN
      lb1 = MAX(lb1_new,lb1_old)
      ub1 = MIN(ub1_new,ub1_old)
      p(lb1:ub1) = work(lb1:ub1)
      DEALLOCATE (work,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine_name,module_name,__LINE__,"work")
      END IF
    END IF

  END SUBROUTINE reallocate_s1

! *****************************************************************************

  SUBROUTINE write_memory(checkpoint,globenv)

!   Purpose: Write the current memory size of all processe in KBytes.

!   History: - Creation (24.09.2002,MK)

!   ***************************************************************************

    CHARACTER(LEN=*), INTENT(IN)             :: checkpoint
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=LEN_TRIM(checkpoint)+7)    :: string
    INTEGER                                  :: group, ipe, mype, npe

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%print%memory) RETURN

    group = globenv%group
    npe = globenv%num_pe
    mype = globenv%mepos

    DO ipe=0,npe-1
      WRITE (UNIT=string,FMT="(I6,A)") ipe,":"//TRIM(checkpoint)
      string = ADJUSTL(string)
      IF (ipe == mype) CALL print_memory(string)
      CALL mp_sync(group)
    END DO

  END SUBROUTINE write_memory

! *****************************************************************************

END MODULE memory_utilities
