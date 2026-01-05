!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief A wrapper around pw_to_openpmd() which accepts particle_list_type
!> \author Ole Schuett, Franz Poeschel
! **************************************************************************************************
MODULE cp_realspace_grid_openpmd
   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE kinds,                           ONLY: dp
   USE particle_list_types,             ONLY: particle_list_type
   USE pw_types,                        ONLY: pw_r3d_rs_type
   USE realspace_grid_openpmd,          ONLY: pw_to_openpmd
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: cp_pw_to_openpmd

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_realspace_grid_openpmd'

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param pw ...
!> \param unit_nr ...
!> \param title ...
!> \param particles ...
!> \param zeff ...
!> \param stride ...
!> \param zero_tails ...
!> \param silent ...
!> \param mpi_io ...
! **************************************************************************************************
   SUBROUTINE cp_pw_to_openpmd( &
      pw, &
      unit_nr, &
      title, &
      particles, &
      zeff, &
      stride, &
      zero_tails, &
      silent, &
      mpi_io &
      )
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: pw
      INTEGER, INTENT(IN)                                :: unit_nr
      CHARACTER(*), INTENT(IN), OPTIONAL                 :: title
      TYPE(particle_list_type), POINTER                  :: particles
      REAL(KIND=dp), DIMENSION(:), OPTIONAL              :: zeff
      INTEGER, DIMENSION(:), OPTIONAL, POINTER           :: stride
      LOGICAL, INTENT(IN), OPTIONAL                      :: zero_tails, silent, mpi_io

      INTEGER                                            :: i, n
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: particles_z
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: particles_r
      TYPE(particle_list_type), POINTER                  :: my_particles

      NULLIFY (my_particles)
      my_particles => particles
      IF (ASSOCIATED(my_particles)) THEN
         n = my_particles%n_els
         ALLOCATE (particles_z(n))
         ALLOCATE (particles_r(3, n))
         DO i = 1, n
            CALL get_atomic_kind(my_particles%els(i)%atomic_kind, z=particles_z(i))
            particles_r(:, i) = my_particles%els(i)%r(:)
         END DO

         CALL pw_to_openpmd(pw=pw, unit_nr=unit_nr, title=title, &
                            particles_z=particles_z, particles_r=particles_r, &
                            particles_zeff=zeff, &
                            stride=stride, zero_tails=zero_tails, &
                            silent=silent, mpi_io=mpi_io)
      ELSE
         CALL pw_to_openpmd(pw=pw, unit_nr=unit_nr, title=title, &
                            stride=stride, zero_tails=zero_tails, &
                            silent=silent, mpi_io=mpi_io)
      END IF

   END SUBROUTINE cp_pw_to_openpmd

END MODULE cp_realspace_grid_openpmd
