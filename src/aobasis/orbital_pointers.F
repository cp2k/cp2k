!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Provides Cartesian and spherical orbital pointers and indices
!> \par History
!>      - reallocate eliminated (17.07.2002,MK)
!>      - restructured and cleaned (20.05.2004,MK)
!> \author Matthias Krack (07.06.2000)
! *****************************************************************************
MODULE orbital_pointers

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
! nsoset_pm: Number of spherical orbitals up to the angular momentum quantum
!         number l for positive m only, i.e. m >= 0

  
  USE kinds,                           ONLY: int_size
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'orbital_pointers'

  INTEGER, SAVE :: current_maxl = -1

  INTEGER, DIMENSION(:), ALLOCATABLE     :: nco,ncoset,nso,nsoset,nsoset_pm
  INTEGER, DIMENSION(:,:), ALLOCATABLE   :: indco,indso,indso_inv
  INTEGER, DIMENSION(:,:), ALLOCATABLE   :: so,soset
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: co,coset

! *** Public subroutines ***

  PUBLIC :: deallocate_orbital_pointers,&
            init_orbital_pointers

! *** Public variables ***

  PUBLIC :: co,&
            coset,&
            current_maxl,&
            indco,&
            indso,&
            indso_inv,&
            nco,&
            ncoset,&
            nso,&
            nsoset,&
            nsoset_pm,&
            soset

CONTAINS

! *****************************************************************************
!> \brief  Allocate and initialize the orbital pointers.
!> \param maxl ...
!> \date   20.05.2004
!> \author MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE create_orbital_pointers(maxl)
    INTEGER, INTENT(IN)                      :: maxl

    CHARACTER(len=*), PARAMETER :: routineN = 'create_orbital_pointers', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iso, istat, l, lx, ly, lz, m

    IF (current_maxl > -1) THEN
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "Orbital pointers are already allocated. "//&
                        "Use the init routine for an update")
    END IF

    IF (maxl < 0) THEN
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "A negative maximum angular momentum quantum "//&
                        "number is invalid")
    END IF

!   *** Number of Cartesian orbitals for each l ***

    ALLOCATE (nco(-1:maxl),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routineN,moduleN,__LINE__,&
                       "nco",(maxl+2)*int_size)
    END IF

    nco(-1) = 0

    DO l=0,maxl
      nco(l) = (l + 1)*(l + 2)/2
    END DO

!   *** Number of Cartesian orbitals up to l ***

    ALLOCATE (ncoset(-1:maxl),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routineN,moduleN,__LINE__,&
                       "ncoset",(maxl+2)*int_size)
    END IF

    ncoset(-1) = 0

    DO l=0,maxl
      ncoset(l) = ncoset(l-1) + nco(l)
    END DO

!   *** Build the Cartesian orbital pointer and the shell orbital pointer ***

    ALLOCATE (co(0:maxl,0:maxl,0:maxl),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routineN,moduleN,__LINE__,&
                       "co",(maxl+1)**3*int_size)
    END IF

    co(:,:,:) = 0

    ALLOCATE (coset(-1:maxl,-1:maxl,-1:maxl),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routineN,moduleN,__LINE__,&
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
      CALL stop_memory(routineN,moduleN,__LINE__,&
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
      CALL stop_memory(routineN,moduleN,__LINE__,&
                       "nso",(maxl+2)*int_size)
    END IF

    nso(-1) = 0

    DO l=0,maxl
      nso(l) = 2*l + 1
    END DO

!   *** Number of spherical orbitals up to l ***

    ALLOCATE (nsoset(-1:maxl),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routineN,moduleN,__LINE__,&
                       "nsoset",(maxl+2)*int_size)
    END IF
    ALLOCATE (nsoset_pm(-1:maxl),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routineN,moduleN,__LINE__,&
                       "nsoset_pm",(maxl+2)*int_size)
    END IF

    nsoset(-1) = 0
    nsoset_pm(-1) = 0

    DO l=0,maxl
      nsoset(l) = nsoset(l-1) + nso(l)
      nsoset_pm(l) = nsoset_pm(l-1) + l + 1
    END DO

    ALLOCATE (indso(2,nsoset(maxl)),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routineN,moduleN,__LINE__,&
                       "indso",2*nsoset(maxl)*int_size)
    END IF
    ! indso_inv: inverse to indso
    ALLOCATE (indso_inv(0:maxl,-maxl:maxl),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routineN,moduleN,__LINE__,&
                    "indso_inv",2*(1+maxl)*(2*maxl+1)*int_size)
    END IF

    indso(:,:) = 0
    indso_inv(:,:) = 0

    iso = 0
    DO l=0,maxl
      DO m=-l,l
        iso = iso + 1
        indso(1:2,iso) = (/l,m/)
        indso_inv(l,m) = iso
      END DO
    END DO

    ALLOCATE (so(0:maxl,-maxl:maxl),soset(0:maxl,-maxl:maxl),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routineN,moduleN,__LINE__,&
                    "so,soset",2*(1+maxl)*(2*maxl+1)*int_size)
    END IF

    soset(:,:) = 0
    DO l = 0,maxl
      DO m = -l,l
        so(l,m) = nso(l) - (l-m)
        soset(l,m) = nsoset(l-1) + nso(l) - (l-m)
      END DO
    END DO

!   *** Save initialization status ***

    current_maxl = maxl

  END SUBROUTINE create_orbital_pointers

! *****************************************************************************
!> \brief   Deallocate the orbital pointers.
!> \date    20.05.2005
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE deallocate_orbital_pointers()

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_orbital_pointers', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat

    IF (current_maxl > -1) THEN

      DEALLOCATE (co,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,"co")
      END IF

      DEALLOCATE (coset,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,"coset")
      END IF

      DEALLOCATE (indco,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,"indco")
      END IF

      DEALLOCATE (indso,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,"indso")
      END IF

      DEALLOCATE (indso_inv,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,"indso_inv")
      END IF

      DEALLOCATE (nco,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,"nco")
      END IF

      DEALLOCATE (ncoset,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,"ncoset")
      END IF

      DEALLOCATE (nso,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,"nso")
      END IF

      DEALLOCATE (nsoset,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,"nsoset")
      END IF

      DEALLOCATE (nsoset_pm,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,"nsoset_pm")
      END IF

      DEALLOCATE (so,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,"nsoset")
      END IF

      DEALLOCATE (soset,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,"nsoset")
      END IF

      current_maxl = -1

    END IF

  END SUBROUTINE deallocate_orbital_pointers

! *****************************************************************************
!> \brief   Initialize or update the orbital pointers.
!> \param maxl ...
!> \date    07.06.2000
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE init_orbital_pointers(maxl)
    INTEGER, INTENT(IN)                      :: maxl

    CHARACTER(len=*), PARAMETER :: routineN = 'init_orbital_pointers', &
      routineP = moduleN//':'//routineN

    IF (maxl < 0) THEN
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "A negative maximum angular momentum quantum "//&
                        "number is invalid")
    END IF

!   *** Check, if the current initialization is sufficient ***

    IF (maxl > current_maxl) THEN
      CALL deallocate_orbital_pointers()
      CALL create_orbital_pointers(maxl)
    END IF

  END SUBROUTINE init_orbital_pointers

END MODULE orbital_pointers
