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

  USE kinds,  ONLY: wp => dp,&
                    wp_size => dp_size

  USE termination, ONLY: stop_memory

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: reallocate

! *****************************************************************************

  INTERFACE reallocate
    MODULE PROCEDURE reallocate_i0,reallocate_i1,reallocate_i2,reallocate_i3,&
                     reallocate_i4,&
                     reallocate_r0,reallocate_r1,reallocate_r2,reallocate_r3,&
                     reallocate_r4
  END INTERFACE

! *****************************************************************************

CONTAINS

! *****************************************************************************

  FUNCTION reallocate_i0(p_old) RESULT(p_new)

!   Purpose: (Re)Allocate a pointer to an integer number.

!   History: - Creation (20.08.1999, Matthias Krack)

!   ***************************************************************************

    INTEGER, POINTER :: p_old

    INTEGER, POINTER :: p_new

!   *** Local variables ***

    INTEGER :: istat

!   ---------------------------------------------------------------------------

    ALLOCATE (p_new,STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory("FUNCTION reallocate_i0 (MODULE memory_utilities)",&
                       "p_new",BIT_SIZE(0)/8)
    END IF

    p_new = 0

    IF (ASSOCIATED(p_old)) THEN
      p_new = p_old
      DEALLOCATE (p_old,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory("FUNCTION reallocate_i0 (MODULE memory_utilities)",&
                         "p_old")
      END IF
    END IF

  END FUNCTION reallocate_i0

! *****************************************************************************

  FUNCTION reallocate_i1(p_old,n1_new,n2_new) RESULT(p_new)

!   Purpose: (Re)Allocate an integer vector with a new dimension.

!   History: - Creation (20.08.1999, Matthias Krack)

!   ***************************************************************************

    INTEGER, INTENT(IN) :: n1_new,n2_new

    INTEGER, DIMENSION(:), POINTER :: p_old

    INTEGER, DIMENSION(:), POINTER :: p_new

!   *** Local variables ***

    INTEGER :: istat,n1,n1_old,n2,n2_old,new_size

!   ---------------------------------------------------------------------------

    ALLOCATE (p_new(n1_new:n2_new),STAT=istat)
    IF (istat /= 0) THEN
      new_size = (n2_new - n1_new + 1)*BIT_SIZE(0)/8
      CALL stop_memory("FUNCTION reallocate_i1 (MODULE memory_utilities)",&
                       "p_new",new_size)
    END IF

    p_new(:) = 0

    IF (ASSOCIATED(p_old)) THEN
      n1_old = LBOUND(p_old,1)
      n2_old = UBOUND(p_old,1)
      n1 = MAX(n1_new,n1_old)
      n2 = MIN(n2_new,n2_old)
      p_new(n1:n2) = p_old(n1:n2)
      DEALLOCATE (p_old,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory("FUNCTION reallocate_i1 (MODULE memory_utilities)",&
                         "p_old")
      END IF
    END IF

  END FUNCTION reallocate_i1

! *****************************************************************************

  FUNCTION reallocate_i2(p_old,m1_new,m2_new,&
                               n1_new,n2_new) RESULT(p_new)

!   Purpose: (Re)Allocate an integer array with new dimensions.

!   History: - Creation (25.06.1999, Matthias Krack)

!   ***************************************************************************

    INTEGER, INTENT(IN) :: m1_new,m2_new,n1_new,n2_new

    INTEGER, DIMENSION(:,:), POINTER :: p_old

    INTEGER, DIMENSION(:,:), POINTER :: p_new

!   *** Local variables ***

    INTEGER :: istat,m1,m1_old,m2,m2_old,n1,n1_old,n2,n2_old,new_size

!   ---------------------------------------------------------------------------

    ALLOCATE (p_new(m1_new:m2_new,n1_new:n2_new),STAT=istat)
    IF (istat /= 0) THEN
      new_size = (m2_new - m1_new + 1)*&
                 (n2_new - n1_new + 1)*BIT_SIZE(0)/8
      CALL stop_memory("FUNCTION reallocate_i2 (MODULE memory_utilities)",&
                       "p_new",new_size)
    END IF

    p_new(:,:) = 0

    IF (ASSOCIATED(p_old)) THEN
      m1_old = LBOUND(p_old,1)
      m2_old = UBOUND(p_old,1)
      n1_old = LBOUND(p_old,2)
      n2_old = UBOUND(p_old,2)
      m1 = MAX(m1_new,m1_old)
      m2 = MIN(m2_new,m2_old)
      n1 = MAX(n1_new,n1_old)
      n2 = MIN(n2_new,n2_old)
      p_new(m1:m2,n1:n2) = p_old(m1:m2,n1:n2)
      DEALLOCATE (p_old,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory("FUNCTION reallocate_i2 (MODULE memory_utilities)",&
                         "p_old")
      END IF
    END IF

  END FUNCTION reallocate_i2

! *****************************************************************************

  FUNCTION reallocate_i3(p_old,l1_new,l2_new,&
                               m1_new,m2_new,&
                               n1_new,n2_new) RESULT(p_new)

!   Purpose: (Re)Allocate an integer array with new dimensions.

!   History: - Creation (25.06.1999, Matthias Krack)

!   ***************************************************************************

    INTEGER, INTENT(IN) :: l1_new,l2_new,m1_new,m2_new,n1_new,n2_new

    INTEGER, DIMENSION(:,:,:), POINTER :: p_old

    INTEGER, DIMENSION(:,:,:), POINTER :: p_new

!   *** Local variables ***

    INTEGER :: istat,l1,l1_old,l2,l2_old,m1,m1_old,m2,m2_old,n1,n1_old,n2,&
               n2_old,new_size

!   ---------------------------------------------------------------------------

    ALLOCATE (p_new(l1_new:l2_new,m1_new:m2_new,n1_new:n2_new),STAT=istat)
    IF (istat /= 0) THEN
      new_size = (l2_new - l1_new + 1)*&
                 (m2_new - m1_new + 1)*&
                 (n2_new - n1_new + 1)*BIT_SIZE(0)/8
      CALL stop_memory("FUNCTION reallocate_i3 (MODULE memory_utilities)",&
                       "p_new",new_size)
    END IF

    p_new(:,:,:) = 0

    IF (ASSOCIATED(p_old)) THEN
      l1_old = LBOUND(p_old,1)
      l2_old = UBOUND(p_old,1)
      m1_old = LBOUND(p_old,2)
      m2_old = UBOUND(p_old,2)
      n1_old = LBOUND(p_old,3)
      n2_old = UBOUND(p_old,3)
      l1 = MAX(l1_new,l1_old)
      l2 = MIN(l2_new,l2_old)
      m1 = MAX(m1_new,m1_old)
      m2 = MIN(m2_new,m2_old)
      n1 = MAX(n1_new,n1_old)
      n2 = MIN(n2_new,n2_old)
      p_new(l1:l2,m1:m2,n1:n2) = p_old(l1:l2,m1:m2,n1:n2)
      DEALLOCATE (p_old,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory("FUNCTION reallocate_i3 (MODULE memory_utilities)",&
                         "p_old")
      END IF
    END IF

  END FUNCTION reallocate_i3

! *****************************************************************************

  FUNCTION reallocate_i4(p_old,k1_new,k2_new,&
                               l1_new,l2_new,&
                               m1_new,m2_new,&
                               n1_new,n2_new) RESULT(p_new)

!   Purpose: (Re)Allocate an integer array with new dimensions.

!   History: - Creation (25.06.1999, Matthias Krack)

!   ***************************************************************************

    INTEGER, INTENT(IN) :: k1_new,k2_new,l1_new,l2_new,m1_new,m2_new,n1_new,&
                           n2_new

    INTEGER, DIMENSION(:,:,:,:), POINTER :: p_old

    INTEGER, DIMENSION(:,:,:,:), POINTER :: p_new

!   *** Local variables ***

    INTEGER :: istat,k1,k1_old,k2,k2_old,l1,l1_old,l2,l2_old,m1,m1_old,m2,&
               m2_old,n1,n1_old,n2,n2_old,new_size

!   ---------------------------------------------------------------------------

    ALLOCATE (p_new(k1_new:k2_new,l1_new:l2_new,m1_new:m2_new,n1_new:n2_new),&
              STAT=istat)
    IF (istat /= 0) THEN
      new_size = (k2_new - k1_new + 1)*&
                 (l2_new - l1_new + 1)*&
                 (m2_new - m1_new + 1)*&
                 (n2_new - n1_new + 1)*BIT_SIZE(0)/8
      CALL stop_memory("FUNCTION reallocate_i4 (MODULE memory_utilities)",&
                       "p_new",new_size)
    END IF

    p_new(:,:,:,:) = 0

    IF (ASSOCIATED(p_old)) THEN
      k1_old = LBOUND(p_old,1)
      k2_old = UBOUND(p_old,1)
      l1_old = LBOUND(p_old,2)
      l2_old = UBOUND(p_old,2)
      m1_old = LBOUND(p_old,3)
      m2_old = UBOUND(p_old,3)
      n1_old = LBOUND(p_old,4)
      n2_old = UBOUND(p_old,4)
      k1 = MAX(k1_new,k1_old)
      k2 = MIN(k2_new,k2_old)
      l1 = MAX(l1_new,l1_old)
      l2 = MIN(l2_new,l2_old)
      m1 = MAX(m1_new,m1_old)
      m2 = MIN(m2_new,m2_old)
      n1 = MAX(n1_new,n1_old)
      n2 = MIN(n2_new,n2_old)
      p_new(k1:k2,l1:l2,m1:m2,n1:n2) = p_old(k1:k2,l1:l2,m1:m2,n1:n2)
      DEALLOCATE (p_old,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory("FUNCTION reallocate_i4 (MODULE memory_utilities)",&
                         "p_old")
      END IF
    END IF

  END FUNCTION reallocate_i4

! *****************************************************************************

  FUNCTION reallocate_r0(p_old) RESULT(p_new)

!   Purpose: (Re)Allocate a pointer to a real number.

!   History: - Creation (10.05.2000, Matthias Krack)

!   ***************************************************************************

    REAL(wp), POINTER :: p_old

    REAL(wp), POINTER :: p_new

!   *** Local variables ***

    INTEGER :: istat

!   ---------------------------------------------------------------------------

    ALLOCATE (p_new,STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory("FUNCTION reallocate_r0 (MODULE memory_utilities)",&
                       "p_new",wp_size)
    END IF

    p_new = 0.0_wp

    IF (ASSOCIATED(p_old)) THEN
      p_new = p_old
      DEALLOCATE (p_old,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory("FUNCTION reallocate_r0 (MODULE memory_utilities)",&
                         "p_old")
      END IF
    END IF

  END FUNCTION reallocate_r0

! *****************************************************************************

  FUNCTION reallocate_r1(p_old,n1_new,n2_new) RESULT(p_new)

!   Purpose: (Re)Allocate a real vector of with a new dimension.

!   History: - Creation (20.08.1999, Matthias Krack)

!   ***************************************************************************

    INTEGER, INTENT(IN) :: n1_new,n2_new

    REAL(wp), DIMENSION(:), POINTER :: p_old

    REAL(wp), DIMENSION(:), POINTER :: p_new

!   *** Local variables ***

    INTEGER :: istat,n1,n1_old,n2,n2_old,new_size

!   ---------------------------------------------------------------------------

    ALLOCATE (p_new(n1_new:n2_new),STAT=istat)
    IF (istat /= 0) THEN
      new_size = (n2_new - n1_new + 1)*wp_size
      CALL stop_memory("FUNCTION reallocate_r1 (MODULE memory_utilities)",&
                       "p_new",new_size)
    END IF

    p_new(:) = 0.0_wp

    IF (ASSOCIATED(p_old)) THEN
      n1_old = LBOUND(p_old,1)
      n2_old = UBOUND(p_old,1)
      n1 = MAX(n1_new,n1_old)
      n2 = MIN(n2_new,n2_old)
      p_new(n1:n2) = p_old(n1:n2)
      DEALLOCATE (p_old,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory("FUNCTION reallocate_r1 (MODULE memory_utilities)",&
                         "p_old")
      END IF
    END IF

  END FUNCTION reallocate_r1

! *****************************************************************************

  FUNCTION reallocate_r2(p_old,m1_new,m2_new,&
                               n1_new,n2_new) RESULT(p_new)

!   Purpose: (Re)Allocate a real array with new dimensions.

!   History: - Creation (25.06.1999, Matthias Krack)

!   ***************************************************************************

    INTEGER, INTENT(IN) :: m1_new,m2_new,n1_new,n2_new

    REAL(wp), DIMENSION(:,:), POINTER :: p_old

    REAL(wp), DIMENSION(:,:), POINTER :: p_new

!   *** Local variables ***

    INTEGER :: istat,m1,m1_old,m2,m2_old,n1,n1_old,n2,n2_old,new_size

!   ---------------------------------------------------------------------------

    ALLOCATE (p_new(m1_new:m2_new,n1_new:n2_new),STAT=istat)
    IF (istat /= 0) THEN
      new_size = (m2_new - m1_new + 1)*&
                 (n2_new - n1_new + 1)*wp_size
      CALL stop_memory("FUNCTION reallocate_r2 (MODULE memory_utilities)",&
                       "p_new",new_size)
    END IF

    p_new(:,:) = 0.0_wp

    IF (ASSOCIATED(p_old)) THEN
      m1_old = LBOUND(p_old,1)
      m2_old = UBOUND(p_old,1)
      n1_old = LBOUND(p_old,2)
      n2_old = UBOUND(p_old,2)
      m1 = MAX(m1_new,m1_old)
      m2 = MIN(m2_new,m2_old)
      n1 = MAX(n1_new,n1_old)
      n2 = MIN(n2_new,n2_old)
      p_new(m1:m2,n1:n2) = p_old(m1:m2,n1:n2)
      DEALLOCATE (p_old,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory("FUNCTION reallocate_r2 (MODULE memory_utilities)",&
                         "p_old")
      END IF
    END IF

  END FUNCTION reallocate_r2

! *****************************************************************************

  FUNCTION reallocate_r3(p_old,l1_new,l2_new,&
                               m1_new,m2_new,&
                               n1_new,n2_new) RESULT(p_new)

!   Purpose: (Re)Allocate a real array with new dimensions.

!   History: - Creation (25.06.1999, Matthias Krack)

!   ***************************************************************************

    INTEGER, INTENT(IN) :: l1_new,l2_new,m1_new,m2_new,n1_new,n2_new

    REAL(wp), DIMENSION(:,:,:), POINTER :: p_old

    REAL(wp), DIMENSION(:,:,:), POINTER :: p_new

!   *** Local variables ***

    INTEGER :: istat,l1,l1_old,l2,l2_old,m1,m1_old,m2,m2_old,n1,n1_old,n2,&
               n2_old,new_size

!   ---------------------------------------------------------------------------

    ALLOCATE (p_new(l1_new:l2_new,m1_new:m2_new,n1_new:n2_new),STAT=istat)
    IF (istat /= 0) THEN
      new_size = (l2_new - l1_new + 1)*&
                 (m2_new - m1_new + 1)*&
                 (n2_new - n1_new + 1)*wp_size
      CALL stop_memory("FUNCTION reallocate_r3 (MODULE memory_utilities)",&
                       "p_new",new_size)
    END IF

    p_new(:,:,:) = 0.0_wp

    IF (ASSOCIATED(p_old)) THEN
      l1_old = LBOUND(p_old,1)
      l2_old = UBOUND(p_old,1)
      m1_old = LBOUND(p_old,2)
      m2_old = UBOUND(p_old,2)
      n1_old = LBOUND(p_old,3)
      n2_old = UBOUND(p_old,3)
      l1 = MAX(l1_new,l1_old)
      l2 = MIN(l2_new,l2_old)
      m1 = MAX(m1_new,m1_old)
      m2 = MIN(m2_new,m2_old)
      n1 = MAX(n1_new,n1_old)
      n2 = MIN(n2_new,n2_old)
      p_new(l1:l2,m1:m2,n1:n2) = p_old(l1:l2,m1:m2,n1:n2)
      DEALLOCATE (p_old,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory("FUNCTION reallocate_r3 (MODULE memory_utilities)",&
                         "p_old")
      END IF
    END IF

  END FUNCTION reallocate_r3

! *****************************************************************************

  FUNCTION reallocate_r4(p_old,k1_new,k2_new,&
                               l1_new,l2_new,&
                               m1_new,m2_new,&
                               n1_new,n2_new) RESULT(p_new)

!   Purpose: (Re)Allocate a real array with new dimensions.

!   History: - Creation (25.06.1999, Matthias Krack)

!   ***************************************************************************

    INTEGER, INTENT(IN) :: k1_new,k2_new,l1_new,l2_new,m1_new,m2_new,n1_new,&
                           n2_new

    REAL(wp), DIMENSION(:,:,:,:), POINTER :: p_old

    REAL(wp), DIMENSION(:,:,:,:), POINTER :: p_new

!   *** Local variables ***

    INTEGER :: istat,k1,k1_old,k2,k2_old,l1,l1_old,l2,l2_old,m1,m1_old,m2,&
               m2_old,n1,n1_old,n2,n2_old,new_size

!   ---------------------------------------------------------------------------

    ALLOCATE (p_new(k1_new:k2_new,l1_new:l2_new,m1_new:m2_new,n1_new:n2_new),&
              STAT=istat)
    IF (istat /= 0) THEN
      new_size = (k2_new - k1_new + 1)*&
                 (l2_new - l1_new + 1)*&
                 (m2_new - m1_new + 1)*&
                 (n2_new - n1_new + 1)*wp_size
      CALL stop_memory("FUNCTION reallocate_r4 (MODULE memory_utilities)",&
                       "p_new",new_size)
    END IF

    p_new(:,:,:,:) = 0.0_wp

    IF (ASSOCIATED(p_old)) THEN
      k1_old = LBOUND(p_old,1)
      k2_old = UBOUND(p_old,1)
      l1_old = LBOUND(p_old,2)
      l2_old = UBOUND(p_old,2)
      m1_old = LBOUND(p_old,3)
      m2_old = UBOUND(p_old,3)
      n1_old = LBOUND(p_old,4)
      n2_old = UBOUND(p_old,4)
      k1 = MAX(k1_new,k1_old)
      k2 = MIN(k2_new,k2_old)
      l1 = MAX(l1_new,l1_old)
      l2 = MIN(l2_new,l2_old)
      m1 = MAX(m1_new,m1_old)
      m2 = MIN(m2_new,m2_old)
      n1 = MAX(n1_new,n1_old)
      n2 = MIN(n2_new,n2_old)
      p_new(k1:k2,l1:l2,m1:m2,n1:n2) = p_old(k1:k2,l1:l2,m1:m2,n1:n2)
      DEALLOCATE (p_old,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory("FUNCTION reallocate_r4 (MODULE memory_utilities)",&
                         "p_old")
      END IF
    END IF

  END FUNCTION reallocate_r4

! *****************************************************************************

END MODULE memory_utilities
