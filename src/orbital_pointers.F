!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
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

  USE kinds, ONLY: int_size

  USE termination, ONLY: stop_memory

  IMPLICIT NONE

  PRIVATE

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

    INTEGER, INTENT(IN) :: maxl

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE init_orbital_pointers (MODULE orbital_pointers)"

!   *** Local variables ***

    INTEGER :: istat,l,lx,ly,lz

!   ---------------------------------------------------------------------------

!   *** Quick return, if the current initialization is sufficient ***

    IF (maxl <= current_maxl) RETURN

!   *** Initialize pointers in the first call ***

    IF (current_maxl > -1) THEN
      DEALLOCATE (co,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"co")
      DEALLOCATE (coset,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"coset")
      DEALLOCATE (indco,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"indco")
      DEALLOCATE (nco,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"nco")
      DEALLOCATE (ncoset,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"ncoset")
      DEALLOCATE (nso,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"nso")
      DEALLOCATE (nsoset,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"nsoset")
    END IF

!   *** Number of Cartesian orbitals for each l ***

    ALLOCATE (nco(-1:maxl),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"nco",(maxl+2)*int_size)

    nco(-1) = 0
    DO l=0,maxl
      nco(l) = (l + 1)*(l + 2)/2
    END DO

!   *** Number of Cartesian orbitals up to l ***

    ALLOCATE (ncoset(-1:maxl),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"ncoset",(maxl+2)*int_size)

    ncoset(-1) = 0
    DO l=0,maxl
      ncoset(l) = ncoset(l-1) + nco(l)
    END DO

!   *** Build the Cartesian orbital pointer and the shell orbital pointer ***

    ALLOCATE (co(0:maxl,0:maxl,0:maxl),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"co",(maxl+1)**3*int_size)
    co(:,:,:) = 0

    ALLOCATE (coset(-1:maxl,-1:maxl,-1:maxl),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"coset",(maxl+2)**3*int_size)
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
    IF (istat /= 0) CALL stop_memory(routine,"indco",3*ncoset(maxl)*int_size)
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
    IF (istat /= 0) CALL stop_memory(routine,"nso",(maxl+2)*int_size)

    nso(-1) = 0
    DO l=0,maxl
      nso(l) = 2*l + 1
    END DO

!   *** Number of spherical orbitals up to l ***

    ALLOCATE (nsoset(-1:maxl),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"nsoset",(maxl+2)*int_size)

    nsoset(-1) = 0
    DO l=0,maxl
      nsoset(l) = nsoset(l-1) + nso(l)
    END DO

!   *** Save the new initialization status ***

    current_maxl = maxl

  END SUBROUTINE init_orbital_pointers

! *****************************************************************************

END MODULE orbital_pointers
