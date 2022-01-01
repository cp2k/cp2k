!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief A wrapper around pw_to_cube() which accepts particle_list_type
!> \author Ole Schuett
! **************************************************************************************************
MODULE cp_realspace_grid_cube
   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE cp_output_handling,              ONLY: cp_mpi_io_get
   USE kinds,                           ONLY: dp
   USE particle_list_types,             ONLY: particle_list_type
   USE pw_types,                        ONLY: pw_type
   USE realspace_grid_cube,             ONLY: cube_to_pw,&
                                              pw_to_cube,&
                                              pw_to_simple_volumetric
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: cp_pw_to_cube, cp_pw_to_simple_volumetric, cp_cube_to_pw

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_realspace_grid_cube'

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param pw ...
!> \param unit_nr ...
!> \param title ...
!> \param particles ...
!> \param stride ...
!> \param zero_tails ...
!> \param silent minimal I/O
!> \param mpi_io True if cube should be written in parallel using MPI
! **************************************************************************************************
   SUBROUTINE cp_pw_to_cube(pw, unit_nr, title, particles, stride, zero_tails, silent, mpi_io)
      TYPE(pw_type), POINTER                             :: pw
      INTEGER, INTENT(IN)                                :: unit_nr
      CHARACTER(*), INTENT(IN), OPTIONAL                 :: title
      TYPE(particle_list_type), POINTER                  :: particles
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

         CALL pw_to_cube(pw=pw, unit_nr=unit_nr, title=title, &
                         particles_z=particles_z, particles_r=particles_r, &
                         stride=stride, zero_tails=zero_tails, &
                         silent=silent, mpi_io=mpi_io)
      ELSE
         CALL pw_to_cube(pw=pw, unit_nr=unit_nr, title=title, &
                         stride=stride, zero_tails=zero_tails, &
                         silent=silent, mpi_io=mpi_io)
      END IF

   END SUBROUTINE cp_pw_to_cube

! **************************************************************************************************
!> \brief Prints grid in a simple format: X Y Z value
!> \param pw ...
!> \param unit_nr ...
!> \param stride ...
!> \param pw2 ...
!> \par History
!>      Created [Vladimir Rybkin] (08.2017)
! **************************************************************************************************
   SUBROUTINE cp_pw_to_simple_volumetric(pw, unit_nr, stride, pw2)
      TYPE(pw_type), POINTER                             :: pw
      INTEGER, INTENT(IN)                                :: unit_nr
      INTEGER, DIMENSION(:), OPTIONAL, POINTER           :: stride
      TYPE(pw_type), OPTIONAL, POINTER                   :: pw2

      IF (.NOT. PRESENT(pw2)) THEN
         CALL pw_to_simple_volumetric(pw, unit_nr, stride)
      ELSE
         CALL pw_to_simple_volumetric(pw, unit_nr, stride, pw2)
      END IF

   END SUBROUTINE cp_pw_to_simple_volumetric

! **************************************************************************************************
!> \brief Thin wrapper around routine cube_to_pw
!> \param grid     pw to read from cube file
!> \param filename name of cube file
!> \param scaling  scale values before storing
!> \param silent   minimal I/O
!> \par History
!>      Created [Nico Holmberg] (09.2018)
! **************************************************************************************************
   SUBROUTINE cp_cube_to_pw(grid, filename, scaling, silent)
      TYPE(pw_type), POINTER                             :: grid
      CHARACTER(len=*), INTENT(in)                       :: filename
      REAL(kind=dp), INTENT(in)                          :: scaling
      LOGICAL, INTENT(in), OPTIONAL                      :: silent

      LOGICAL                                            :: parallel_read

      ! Determine whether to use MPI I/O for reading cube filename
      parallel_read = .TRUE.
      ! Parallel routine falls back to stream read in serial mode,
      ! but it has slight overhead compared to sequential read
      ! Therefore, we use sequential version in serial mode
      IF (grid%pw_grid%para%group_size == 1) parallel_read = .FALSE.
      ! Check if MPI I/O was disabled in GLOBAL section
      IF (.NOT. cp_mpi_io_get()) parallel_read = .FALSE.

      CALL cube_to_pw(grid, filename, scaling, parallel_read, silent=silent)

   END SUBROUTINE cp_cube_to_pw

END MODULE cp_realspace_grid_cube
