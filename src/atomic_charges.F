!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2006  CP2K group                                            !
!-----------------------------------------------------------------------------!


!!****** cp2k/atomic_charges [1.0] *
!!
!!   NAME
!!     atomic_charges 
!!
!!   FUNCTION
!!     simple routine to print charges for all atomic charge methods
!!     (currently mulliken, lowdin and ddapc)
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     Joost VandeVondele [2006.03]         
!!
!!   SOURCE
!******************************************************************************

MODULE atomic_charges
  USE atomic_kind_types,               ONLY: get_atomic_kind
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atomic_charges'

! *** Public subroutines ***

  PUBLIC :: print_atomic_charges

! *****************************************************************************

CONTAINS

!!****f* atomic_charges/print_atomic_charges*
!!
!!   NAME
!!     print_atomic_charges
!!
!!   FUNCTION
!!     generates a unified output format for atomic charges
!!
!!   NOTES
!!     charges are computed per spin in the LSD case
!!
!!   INPUTS
!!    - electronic_charges : (natom,nspin), the number of electrons of (so positive) per spin 
!!                           if (nspin==1) it is the sum of alpha and beta electrons
!!    - atomic_charges     : truely the atomic charge (taking Z into account, atoms negative, no spin)
!!
!!   MODIFICATION HISTORY
!!     03.2006 created [Joost VandeVondele]
!!
!!   SOURCE
!!*** **********************************************************************
  SUBROUTINE print_atomic_charges(particle_set, scr, title, electronic_charges, atomic_charges)

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER                                  :: scr
    CHARACTER(LEN=*)                         :: title
    REAL(KIND=dp), DIMENSION(:, :), OPTIONAL :: electronic_charges
    REAL(KIND=dp), DIMENSION(:), OPTIONAL    :: atomic_charges

    CHARACTER(LEN=*), PARAMETER :: routineN = 'print_atomic_charges'

    CHARACTER(LEN=2)                         :: element_symbol
    INTEGER                                  :: handle, iatom, natom, nspin
    REAL(KIND=dp)                            :: total_charge, zeff

    CALL timeset(routineN,handle)
    
    IF (PRESENT(electronic_charges)) THEN
      nspin=SIZE(electronic_charges,2)
      natom=SIZE(electronic_charges,1)
    ELSE
      natom=SIZE(atomic_charges,1)
      nspin=0
    ENDIF

    IF (scr>0) THEN
       WRITE(scr,'(T2,A)') title
       SELECT CASE (nspin)
       CASE(0,1)
       WRITE(scr,'(A)') "  Atom     |    Charge"
       CASE DEFAULT
       WRITE(scr,'(A)') "  Atom     |    Charge | Spin diff charge"
       END SELECT
       total_charge = 0.0_dp
       IF (SIZE(particle_set) .NE. natom) &
           CALL stop_program(routineN,"unexpected number of atoms/charges")
       WRITE(scr,'(A)') ""
       DO iatom=1,natom
          CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                             element_symbol=element_symbol, zeff=zeff)

          SELECT CASE (nspin)
          CASE(0)
             WRITE(scr,'(I6,A2,A2,F12.6)') iatom,"  ",element_symbol,atomic_charges(iatom)
             total_charge=total_charge+atomic_charges(iatom)
          CASE(1)
             WRITE(scr,'(I6,A2,A2,F12.6)') iatom,"  ",element_symbol,zeff-electronic_charges(iatom,1)
             total_charge=total_charge+zeff-electronic_charges(iatom,1)
          CASE DEFAULT
             WRITE(scr,'(I6,A2,A2,2F12.6)') iatom,"  ",element_symbol, &
                      zeff-(electronic_charges(iatom,1)+electronic_charges(iatom,2)), &
                           (electronic_charges(iatom,1)-electronic_charges(iatom,2))
             total_charge=total_charge+ zeff-(electronic_charges(iatom,1)+electronic_charges(iatom,2))
          END SELECT
       ENDDO
       WRITE(scr,'(A,F10.6)') "  Total     ",total_charge
       WRITE(scr,'(A)') ""
    ENDIF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE print_atomic_charges

END MODULE atomic_charges
