!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2003  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/orbital_pointers [1.0] *
!!
!!   NAME
!!     orbital_pointers
!!
!!   FUNCTION
!!      Define all orbital pointers.
!!
!!   AUTHOR
!!     Matthias Krack (07.06.2000)
!!
!!   MODIFICATION HISTORY
!!     reallocate eliminated (17.07.2002,MK)
!!
!!   SOURCE
!******************************************************************************

MODULE orbital_pointers

! *****************************************************************************

! co    : Cartesian orbital pointer for a orbital shell.
! coset : Cartesian orbital pointer for a set of orbitals.
! nco   : Number of Cartesian orbitals for the angular momentum quantum
!         number l.
! ncoset: Number of Cartesian orbitals up to the angular momentum quantum
!         number l.
! nso   : Number of spherical orbitals for the angular momentum quantum
!         number l.
! nsoset: Number of spherical orbitals up to the angular momentum quantum
!         number l.

! *****************************************************************************

  USE kinds,                           ONLY: int_size
  USE termination,                     ONLY: stop_memory

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "orbital_pointers"

  INTEGER, SAVE :: current_maxl = -1

  INTEGER, DIMENSION(:), ALLOCATABLE     :: nco,ncoset,nso,nsoset
  INTEGER, DIMENSION(:,:), ALLOCATABLE   :: indco
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: co,coset

! *** Public subroutines ***

  PUBLIC :: init_orbital_pointers

! *** Public variables ***

  PUBLIC :: co,&
            coset,&
            current_maxl,&
            indco,&
            nco,&
            ncoset,&
            nso

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE init_orbital_pointers(maxl)

!   Purpose: Allocate and initialize the orbital pointers.

!   History: - Creation (07.06.2000, Matthias Krack)

!   ***************************************************************************

    INTEGER, INTENT(IN)                      :: maxl

    CHARACTER(LEN=*), PARAMETER :: routine_name = "init_orbital_pointers"

    INTEGER                                  :: istat, l, lx, ly, lz

!   ---------------------------------------------------------------------------
!   *** Quick return, if the current initialization is sufficient ***

    IF (maxl <= current_maxl) RETURN

!   *** Deallocate all arrays, if it is not the first call ***

    IF (current_maxl > -1) THEN
      DEALLOCATE (co,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine_name,module_name,__LINE__,"co")
      END IF
      DEALLOCATE (coset,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine_name,module_name,__LINE__,"coset")
      END IF
      DEALLOCATE (indco,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine_name,module_name,__LINE__,"indco")
      END IF
      DEALLOCATE (nco,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine_name,module_name,__LINE__,"nco")
      END IF
      DEALLOCATE (ncoset,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine_name,module_name,__LINE__,"ncoset")
      END IF
      DEALLOCATE (nso,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine_name,module_name,__LINE__,"nso")
      END IF
      DEALLOCATE (nsoset,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine_name,module_name,__LINE__,"nsoset")
      END IF
    END IF

!   *** Number of Cartesian orbitals for each l ***

    ALLOCATE (nco(-1:maxl),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine_name,module_name,__LINE__,&
                       "nco",(maxl+2)*int_size)
    END IF

    nco(-1) = 0

    DO l=0,maxl
      nco(l) = (l + 1)*(l + 2)/2
    END DO

!   *** Number of Cartesian orbitals up to l ***

    ALLOCATE (ncoset(-1:maxl),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine_name,module_name,__LINE__,&
                       "ncoset",(maxl+2)*int_size)
    END IF

    ncoset(-1) = 0

    DO l=0,maxl
      ncoset(l) = ncoset(l-1) + nco(l)
    END DO

!   *** Build the Cartesian orbital pointer and the shell orbital pointer ***

    ALLOCATE (co(0:maxl,0:maxl,0:maxl),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine_name,module_name,__LINE__,&
                       "co",(maxl+1)**3*int_size)
    END IF

    co(:,:,:) = 0

    ALLOCATE (coset(-1:maxl,-1:maxl,-1:maxl),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine_name,module_name,__LINE__,&
                       "coset",(maxl+2)**3*int_size)
    END IF

    coset(:,:,:) = 0

    coset(-1,:,:) = 1
    coset(:,-1,:) = 1
    coset(:,:,-1) = 1

    DO lx=0,maxl
      DO ly=0,maxl
        DO lz=0,maxl
          l = lx + ly + lz
          IF (l > maxl) CYCLE
          co(lx,ly,lz) = 1 + (l - lx)*(l - lx + 1)/2 + lz
          coset(lx,ly,lz) = ncoset(l-1) + co(lx,ly,lz)
        END DO
      END DO
    END DO

    ALLOCATE (indco(3,ncoset(maxl)),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine_name,module_name,__LINE__,&
                       "indco",3*ncoset(maxl)*int_size)
    END IF

    indco(:,:) = 0

    DO l=0,maxl
      DO lx=0,l
        DO ly=0,l-lx
          lz = l - lx - ly
          indco(1:3,coset(lx,ly,lz)) = (/lx,ly,lz/)
        END DO
      END DO
    END DO

!   *** Number of spherical orbitals for each l ***

    ALLOCATE (nso(-1:maxl),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine_name,module_name,__LINE__,&
                       "nso",(maxl+2)*int_size)
    END IF

    nso(-1) = 0

    DO l=0,maxl
      nso(l) = 2*l + 1
    END DO

!   *** Number of spherical orbitals up to l ***

    ALLOCATE (nsoset(-1:maxl),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine_name,module_name,__LINE__,&
                       "nsoset",(maxl+2)*int_size)
    END IF

    nsoset(-1) = 0

    DO l=0,maxl
      nsoset(l) = nsoset(l-1) + nso(l)
    END DO

!   *** Save the new initialization status ***

    current_maxl = maxl

  END SUBROUTINE init_orbital_pointers

! *****************************************************************************

END MODULE orbital_pointers
