!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief simple routine to print charges for all atomic charge methods
!>      (currently mulliken, lowdin and ddapc)
!> \par History
!>      Joost VandeVondele [2006.03]
! *****************************************************************************
MODULE atomic_charges
  USE atomic_kind_types,               ONLY: get_atomic_kind
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             qs_kind_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atomic_charges'

  PUBLIC :: print_atomic_charges, print_multipoles

CONTAINS

! *****************************************************************************
!> \brief generates a unified output format for atomic charges
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param scr ...
!> \param title ...
!> \param electronic_charges (natom,nspin), the number of electrons of (so positive) per spin
!>                            if (nspin==1) it is the sum of alpha and beta electrons
!> \param atomic_charges truely the atomic charge (taking Z into account, atoms negative, no spin)
!> \par History
!>      03.2006 created [Joost VandeVondele]
!> \note
!>      charges are computed per spin in the LSD case
! *****************************************************************************
  SUBROUTINE print_atomic_charges(particle_set, qs_kind_set, scr, title, electronic_charges,&
                                  atomic_charges)

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    INTEGER                                  :: scr
    CHARACTER(LEN=*)                         :: title
    REAL(KIND=dp), DIMENSION(:, :), OPTIONAL :: electronic_charges
    REAL(KIND=dp), DIMENSION(:), OPTIONAL    :: atomic_charges

    CHARACTER(len=*), PARAMETER :: routineN = 'print_atomic_charges', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: element_symbol
    INTEGER                                  :: handle, iatom, ikind, natom, &
                                                nspin
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
       IF(title=="RESP charges:") THEN
        WRITE(scr,'(A)') "  Type |   Atom   |    Charge"
       ELSE
        WRITE(scr,'(A)') "  Atom     |    Charge"
       ENDIF
       CASE DEFAULT
       WRITE(scr,'(A)') "  Atom     |    Charge | Spin diff charge"
       END SELECT
       total_charge = 0.0_dp
       IF (SIZE(particle_set) /= natom) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "Unexpected number of atoms/charges")
       END IF
       WRITE(scr,'(A)') ""
       DO iatom=1,natom
          CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                             element_symbol=element_symbol, kind_number=ikind)
          CALL get_qs_kind(qs_kind_set(ikind), zeff=zeff)

          SELECT CASE (nspin)
          CASE(0)
             IF(title=="RESP charges:") THEN
              WRITE(scr,'(T3,A4,2X,I6,A2,A2,F12.6)') "RESP",iatom,"  ",element_symbol,atomic_charges(iatom)
              total_charge=total_charge+atomic_charges(iatom)
             ELSE
              WRITE(scr,'(I6,A2,A2,F12.6)') iatom,"  ",element_symbol,atomic_charges(iatom)
              total_charge=total_charge+atomic_charges(iatom)
             ENDIF
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
       IF(title=="RESP charges:") THEN
          WRITE(scr,'(A,F10.6)') "  Total             ",total_charge
       ELSE
          WRITE(scr,'(A,F10.6)') "  Total     ",total_charge
       ENDIF
       WRITE(scr,'(A)') ""
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE print_atomic_charges

! *****************************************************************************
!> \brief ...
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param scr ...
!> \param title ...
!> \param charge ...
!> \param dipole ...
!> \param quadrupole ...
! *****************************************************************************
  SUBROUTINE print_multipoles(particle_set, qs_kind_set, scr, title, charge, dipole, quadrupole)

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    INTEGER                                  :: scr
    CHARACTER(LEN=*)                         :: title
    REAL(KIND=dp), DIMENSION(:), OPTIONAL    :: charge
    REAL(KIND=dp), DIMENSION(:, :), OPTIONAL :: dipole
    REAL(KIND=dp), DIMENSION(:, :, :), &
      OPTIONAL                               :: quadrupole

    CHARACTER(len=*), PARAMETER :: routineN = 'print_multipoles', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: element_symbol
    INTEGER                                  :: handle, i, iatom, ikind, natom
    REAL(KIND=dp)                            :: zeff

    CALL timeset(routineN,handle)

    natom=0
    IF (PRESENT(charge)) THEN
      natom=SIZE(charge)
    ENDIF

    IF (scr>0) THEN

       WRITE(scr,'(T2,A)') 'multipoles:'

       DO iatom=1,natom
          CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                             element_symbol=element_symbol, kind_number=ikind)
          CALL get_qs_kind(qs_kind_set(ikind), zeff=zeff)

          WRITE(scr,'(a,i5)') ' iatom= ',iatom
          WRITE(scr,'(a,a2)') ' element_symbol= ',element_symbol
          WRITE(scr,'(a,f20.10)') ' zeff= ',zeff

          WRITE(scr,'(a, f20.10)') 'charge =     ',charge(iatom)
          WRITE(scr,'(a,3f20.10)') 'dipole =     ',dipole(:,iatom)
          WRITE(scr,'(a)') 'quadrupole = '
          DO i=1,3
            WRITE(scr,'(3f20.10)') quadrupole(i,:,iatom)
          ENDDO

       ENDDO
       WRITE(scr,'(A)') ""
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE print_multipoles

END MODULE atomic_charges
