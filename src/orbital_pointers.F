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
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE orbital_pointers

! *****************************************************************************

! Index:

! SUBROUTINE init_orbital_pointers(maxl)

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

  IMPLICIT NONE

  PRIVATE

  INTEGER :: current_maxl = -1

  INTEGER, DIMENSION(:), POINTER     :: nco,ncoset,nso,nsoset
  INTEGER, DIMENSION(:,:), POINTER   :: indco
  INTEGER, DIMENSION(:,:,:), POINTER :: co,coset

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

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE init_orbital_pointers(maxl)

!   Purpose: Allocate and initialize the orbital pointers.

!   History: - Creation (07.06.2000, Matthias Krack)

!   ***************************************************************************

    USE memory_utilities, ONLY: reallocate

    INTEGER, INTENT(IN) :: maxl

!   *** Local variables ***

    INTEGER :: l,lx,ly,lz

!   ---------------------------------------------------------------------------

!   *** Quick return, if the current initialization is sufficient ***

    IF (maxl <= current_maxl) RETURN

!   *** Initialize pointers in the first call ***

    IF (current_maxl == -1) NULLIFY(co,coset,indco,nco,ncoset,nso,nsoset)

!   *** Number of Cartesian orbitals for each l ***

    nco => reallocate(nco,-1,maxl)

    DO l=0,maxl
      nco(l) = (l + 1)*(l + 2)/2
    END DO

!   *** Number of Cartesian orbitals up to l ***

    ncoset => reallocate(ncoset,-1,maxl)

    DO l=0,maxl
      ncoset(l) = ncoset(l-1) + nco(l)
    END DO

!   *** Build the Cartesian orbital pointer and the shell orbital pointer ***

    co => reallocate(co,0,maxl,0,maxl,0,maxl)

    coset => reallocate(coset,0,maxl,0,maxl,0,maxl)

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

    indco => reallocate(indco,1,3,1,ncoset(maxl))

    DO l=0,maxl
      DO lx=0,l
        DO ly=0,l-lx
          lz = l - lx - ly
          indco(1:3,coset(lx,ly,lz)) = (/lx,ly,lz/)
        END DO
      END DO
    END DO

!   *** Number of spherical orbitals for each l ***

    nso => reallocate(nso,-1,maxl)

    DO l=0,maxl
      nso(l) = 2*l + 1
    END DO

!   *** Number of spherical orbitals up to l ***

    nsoset => reallocate(nsoset,-1,maxl)

    DO l=0,maxl
      nsoset(l) = nsoset(l-1) + nso(l)
    END DO

!   *** Save the new initialization status ***

    current_maxl = maxl

  END SUBROUTINE init_orbital_pointers

! *****************************************************************************

END MODULE orbital_pointers
