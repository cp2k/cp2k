!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Types needed for a for a Harris model calculation
!> \par History
!>       2024.07 created
!> \author JGH
! **************************************************************************************************
MODULE qs_harris_types
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE pw_types,                        ONLY: pw_r3d_rs_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_harris_types'

! *****************************************************************************
   TYPE rho_vec_type
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)       :: rvecs
   END TYPE rho_vec_type

   TYPE harris_rhoin_type
      CHARACTER(LEN=default_string_length)             :: basis_type = "NDef"
      TYPE(rho_vec_type), ALLOCATABLE, DIMENSION(:, :) :: rhovec
      TYPE(rho_vec_type), ALLOCATABLE, DIMENSION(:, :) :: intvec
      INTEGER                                          :: nspin = 0
      INTEGER                                          :: nbas = 0
      INTEGER, ALLOCATABLE, DIMENSION(:, :)            :: basptr
      LOGICAL                                          :: frozen = .FALSE.
   END TYPE harris_rhoin_type

   TYPE harris_energy_type
      REAL(KIND=dp)                                    :: eharris = 0.0_dp
      REAL(KIND=dp)                                    :: eband = 0.0_dp
      REAL(KIND=dp)                                    :: exc_correction = 0.0_dp
      REAL(KIND=dp)                                    :: eh_correction = 0.0_dp
      REAL(KIND=dp)                                    :: ewald_correction = 0.0_dp
      REAL(KIND=dp)                                    :: dispersion = 0.0_dp
   END TYPE harris_energy_type

! *****************************************************************************
!> \brief Contains information on the Harris method
!> \par History
!>       07.2024 created
!> \author JGH
! *****************************************************************************
   TYPE harris_type
      INTEGER                                          :: energy_functional = 0
      INTEGER                                          :: density_source = 0
      INTEGER                                          :: orbital_basis = 0
      !
      TYPE(harris_energy_type)                         :: energy
      !
      TYPE(harris_rhoin_type)                          :: rhoin
      !
      TYPE(pw_r3d_rs_type)                             :: vh_rspace = pw_r3d_rs_type()
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER      :: vxc_rspace => Null()

      !
      LOGICAL                                          :: debug_forces = .FALSE.
      LOGICAL                                          :: debug_stress = .FALSE.
   END TYPE harris_type
! **************************************************************************************************

   PUBLIC :: harris_type, harris_energy_type, harris_env_release, &
             harris_print_energy, harris_rhoin_type, harris_rhoin_init

! **************************************************************************************************

CONTAINS

! **************************************************************************************************

! **************************************************************************************************
!> \brief ...
!> \param iounit ...
!> \param energy ...
! **************************************************************************************************
   SUBROUTINE harris_print_energy(iounit, energy)
      INTEGER, INTENT(IN)                                :: iounit
      TYPE(harris_energy_type)                           :: energy

      IF (iounit > 0) THEN
         WRITE (UNIT=iounit, FMT="(/,(T2,A))") "HARRIS MODEL ENERGY INFORMATION"
         WRITE (UNIT=iounit, FMT="((T3,A,T56,F25.14))") &
            "Harris model energy:                           ", energy%eharris, &
            "Band energy:                                   ", energy%eband, &
            "Hartree correction energy:                     ", energy%eh_correction, &
            "XC correction energy:                          ", energy%exc_correction, &
            "Ewald sum correction energy:                   ", energy%ewald_correction, &
            "Dispersion energy (pair potential):            ", energy%dispersion
      END IF

   END SUBROUTINE harris_print_energy

! **************************************************************************************************
!> \brief ...
!> \param rhoin ...
!> \param basis_type ...
!> \param qs_kind_set ...
!> \param atomic_kind_set ...
!> \param local_particles ...
!> \param nspin ...
! **************************************************************************************************
   SUBROUTINE harris_rhoin_init(rhoin, basis_type, qs_kind_set, atomic_kind_set, &
                                local_particles, nspin)
      TYPE(harris_rhoin_type)                            :: rhoin
      CHARACTER(LEN=*)                                   :: basis_type
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(distribution_1d_type), POINTER                :: local_particles
      INTEGER, INTENT(IN)                                :: nspin

      INTEGER                                            :: iatom, ikind, iptr, ispin, natom, nkind, &
                                                            nparticle_local, nsgf
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of, nbasf
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind

      CALL harris_rhoin_release(rhoin)

      rhoin%basis_type = basis_type
      rhoin%nspin = nspin

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                               atom_of_kind=atom_of_kind, kind_of=kind_of)
      natom = SIZE(atom_of_kind)
      nkind = SIZE(qs_kind_set)

      ALLOCATE (nbasf(nkind))
      DO ikind = 1, nkind
         qs_kind => qs_kind_set(ikind)
         CALL get_qs_kind(qs_kind, basis_set=basis_set, basis_type=basis_type)
         CALL get_gto_basis_set(basis_set, nsgf=nsgf)
         nbasf(ikind) = nsgf
      END DO

      ALLOCATE (rhoin%basptr(natom, 2))
      iptr = 1
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         rhoin%basptr(iatom, 1) = iptr
         iptr = iptr + nbasf(ikind)
         rhoin%basptr(iatom, 2) = iptr - 1
      END DO
      rhoin%nbas = iptr - 1

      ALLOCATE (rhoin%rhovec(nkind, nspin))
      DO ispin = 1, nspin
         DO ikind = 1, nkind
            nsgf = nbasf(ikind)
            nparticle_local = local_particles%n_el(ikind)
            ALLOCATE (rhoin%rhovec(ikind, ispin)%rvecs(nsgf, nparticle_local))
         END DO
      END DO

      ALLOCATE (rhoin%intvec(nkind, nspin))
      DO ispin = 1, nspin
         DO ikind = 1, nkind
            nsgf = nbasf(ikind)
            nparticle_local = local_particles%n_el(ikind)
            ALLOCATE (rhoin%intvec(ikind, ispin)%rvecs(nsgf, nparticle_local))
         END DO
      END DO

      DEALLOCATE (nbasf)

   END SUBROUTINE harris_rhoin_init

! **************************************************************************************************
!> \brief ...
!> \param harris_env ...
! **************************************************************************************************
   SUBROUTINE harris_env_release(harris_env)
      TYPE(harris_type), POINTER                         :: harris_env

      INTEGER                                            :: iab

      IF (ASSOCIATED(harris_env)) THEN
         !
         CALL harris_rhoin_release(harris_env%rhoin)
         !
         IF (ASSOCIATED(harris_env%vh_rspace%pw_grid)) THEN
            CALL harris_env%vh_rspace%release()
         END IF
         IF (ASSOCIATED(harris_env%vxc_rspace)) THEN
            DO iab = 1, SIZE(harris_env%vxc_rspace)
               CALL harris_env%vxc_rspace(iab)%release()
            END DO
            DEALLOCATE (harris_env%vxc_rspace)
         END IF
         !
         DEALLOCATE (harris_env)
      END IF

      NULLIFY (harris_env)

   END SUBROUTINE harris_env_release

! **************************************************************************************************
!> \brief ...
!> \param rhoin ...
! **************************************************************************************************
   SUBROUTINE harris_rhoin_release(rhoin)
      TYPE(harris_rhoin_type)                            :: rhoin

      INTEGER                                            :: i, j

      IF (ALLOCATED(rhoin%rhovec)) THEN
         DO i = 1, SIZE(rhoin%rhovec, 2)
            DO j = 1, SIZE(rhoin%rhovec, 1)
               IF (ALLOCATED(rhoin%rhovec(j, i)%rvecs)) THEN
                  DEALLOCATE (rhoin%rhovec(j, i)%rvecs)
               END IF
            END DO
         END DO
         DEALLOCATE (rhoin%rhovec)
      END IF
      IF (ALLOCATED(rhoin%intvec)) THEN
         DO i = 1, SIZE(rhoin%intvec, 2)
            DO j = 1, SIZE(rhoin%intvec, 1)
               IF (ALLOCATED(rhoin%intvec(j, i)%rvecs)) THEN
                  DEALLOCATE (rhoin%intvec(j, i)%rvecs)
               END IF
            END DO
         END DO
         DEALLOCATE (rhoin%intvec)
      END IF
      IF (ALLOCATED(rhoin%basptr)) THEN
         DEALLOCATE (rhoin%basptr)
      END IF
      rhoin%basis_type = "NDef"
      rhoin%nspin = 0
      rhoin%nbas = 0
      rhoin%frozen = .FALSE.

   END SUBROUTINE harris_rhoin_release

END MODULE qs_harris_types
