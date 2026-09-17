!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Main setup driver for RI-RS grids.
!> \par History
!>      09.2026 created
! **************************************************************************************************
MODULE gw_ri_rs_grid_setup_main
   USE gw_ri_rs_grid_from_file,         ONLY: read_ri_rs_grid_from_file
   USE gw_ri_rs_grid_optimization,      ONLY: optimize_ri_rs_grid
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE particle_types,                  ONLY: particle_type
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE util,                            ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: setup_ri_rs_grid

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_ri_rs_grid_setup_main'

CONTAINS

! **************************************************************************************************
!> \brief Get RI-RS grid points {r_l}, either by on-the-fly optimization or
!>        reading pretabulated atomic grids
!> \param bs_env Band-structure environment containing GW RI-RS parameters.
!> \param grid_points x,y,z RI-RS grid coordinates, size (3, ngrid)
! **************************************************************************************************
   SUBROUTINE setup_ri_rs_grid(bs_env, grid_points)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)            :: grid_points(:, :)

      IF (bs_env%ri_rs%grid_opt%enabled) THEN

         ! Initialize Lebedev grids for every element from CP2K routines and then
         ! optimize the coordinates of these grid points to minimize RIRS error
         !
         ! [(μν|P) - Σ_l ϕ_μ(r_l) ϕ_ν(r_l) Z_lP^(A)]²
         !
         CALL optimize_ri_rs_grid(bs_env)

      ELSE

         ! Read pretabulated atom-relative RI-RS grids from the data files
         CALL read_ri_rs_grid_from_file(bs_env)

      END IF

      ! Move atom-relative grids to their atomic centres and collect the global grid.
      CALL assemble_ri_rs_grid(bs_env, grid_points)

   END SUBROUTINE setup_ri_rs_grid

! **************************************************************************************************
!> \brief Move atom-relative RI-RS grids to their atomic centres and assemble the global grid.
!> \param bs_env ...
!> \param ri_rs_grid_points ...
! **************************************************************************************************
   SUBROUTINE assemble_ri_rs_grid(bs_env, ri_rs_grid_points)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)            :: ri_rs_grid_points(:, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'assemble_ri_rs_grid'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CPASSERT(ALLOCATED(bs_env%ri_rs%atomic_grids))
      CALL assemble_ri_rs_grid_points(bs_env, ri_rs_grid_points)
      CALL release_atomic_grids(bs_env)

      CALL timestop(handle)

   END SUBROUTINE assemble_ri_rs_grid

! **************************************************************************************************
!> \brief Assemble the global RI-RS grid in spatial atom order from atom-relative grids.
!> \param bs_env ...
!> \param ri_rs_grid_points ...
! **************************************************************************************************
   SUBROUTINE assemble_ri_rs_grid_points(bs_env, ri_rs_grid_points)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)            :: ri_rs_grid_points(:, :)

      INTEGER                                            :: atom_grid_end, atom_grid_start, iatom, &
                                                            ilayout, natom
      INTEGER, ALLOCATABLE                               :: atom_grid_offsets(:), atom_order(:)
      REAL(KIND=dp)                                      :: atom_center(3)

      natom = bs_env%n_atom
      CALL spatial_atom_order(bs_env%ri_rs%particle_set, atom_order)

      ALLOCATE (bs_env%ri_rs%grid_atom_boundaries(natom + 1), atom_grid_offsets(natom))
      bs_env%ri_rs%n_grid_points = 0
      DO ilayout = 1, natom
         iatom = atom_order(ilayout)
         atom_grid_offsets(iatom) = bs_env%ri_rs%n_grid_points + 1
         bs_env%ri_rs%grid_atom_boundaries(ilayout) = bs_env%ri_rs%n_grid_points + 1
         bs_env%ri_rs%n_grid_points = bs_env%ri_rs%n_grid_points + &
                                      bs_env%ri_rs%atomic_grids(iatom)%npts
      END DO
      bs_env%ri_rs%grid_atom_boundaries(natom + 1) = bs_env%ri_rs%n_grid_points + 1

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, FMT="(T2,A,T69,I12)") &
            'Total grid points used for RI-RS:', bs_env%ri_rs%n_grid_points
         WRITE (bs_env%unit_nr, "(A)") ' '
      END IF

      ALLOCATE (ri_rs_grid_points(3, bs_env%ri_rs%n_grid_points))
      !$OMP PARALLEL DO DEFAULT(NONE) &
      !$OMP SHARED(ri_rs_grid_points, atom_grid_offsets, bs_env, natom) &
      !$OMP PRIVATE(iatom, atom_center, atom_grid_start, atom_grid_end) &
      !$OMP SCHEDULE(DYNAMIC, 1)
      DO iatom = 1, natom
         atom_center(:) = bs_env%ri_rs%particle_set(iatom)%r(:)
         atom_grid_start = atom_grid_offsets(iatom)
         atom_grid_end = atom_grid_start + bs_env%ri_rs%atomic_grids(iatom)%npts - 1

         ri_rs_grid_points(1, atom_grid_start:atom_grid_end) = &
            bs_env%ri_rs%atomic_grids(iatom)%raw_points(1, :) + atom_center(1)
         ri_rs_grid_points(2, atom_grid_start:atom_grid_end) = &
            bs_env%ri_rs%atomic_grids(iatom)%raw_points(2, :) + atom_center(2)
         ri_rs_grid_points(3, atom_grid_start:atom_grid_end) = &
            bs_env%ri_rs%atomic_grids(iatom)%raw_points(3, :) + atom_center(3)
      END DO
      !$OMP END PARALLEL DO
   END SUBROUTINE assemble_ri_rs_grid_points

! **************************************************************************************************
!> \brief Release the atom-relative RI-RS grids after assembling the global molecular grid.
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE release_atomic_grids(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      INTEGER                                            :: iatom

      DO iatom = 1, SIZE(bs_env%ri_rs%atomic_grids)
         DEALLOCATE (bs_env%ri_rs%atomic_grids(iatom)%raw_points)
      END DO
      DEALLOCATE (bs_env%ri_rs%atomic_grids)
   END SUBROUTINE release_atomic_grids

! **************************************************************************************************
!> \brief Order atoms by a three-dimensional Morton code so consecutive grid runs remain local.
!> \param particle_set ...
!> \param order ...
! **************************************************************************************************
   SUBROUTINE spatial_atom_order(particle_set, order)
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: order

      CHARACTER(LEN=*), PARAMETER :: routineN = 'spatial_atom_order'
      INTEGER, PARAMETER                                 :: nbits = 21

      INTEGER                                            :: handle, iatom, idimension, natom
      INTEGER(KIND=int_8)                                :: cmax, integer_coordinate(3), m1, m2, m3
      INTEGER(KIND=int_8), ALLOCATABLE                   :: morton_code(:)
      REAL(KIND=dp)                                      :: hi(3), lo(3), span(3)

      CALL timeset(routineN, handle)

      natom = SIZE(particle_set)
      ALLOCATE (order(natom), morton_code(natom))
      cmax = ISHFT(1_int_8, nbits) - 1_int_8

      lo(:) = HUGE(1.0_dp)
      hi(:) = -HUGE(1.0_dp)
      DO iatom = 1, natom
         DO idimension = 1, 3
            lo(idimension) = MIN(lo(idimension), particle_set(iatom)%r(idimension))
            hi(idimension) = MAX(hi(idimension), particle_set(iatom)%r(idimension))
         END DO
      END DO
      span(:) = hi(:) - lo(:)
      DO idimension = 1, 3
         IF (span(idimension) <= 0.0_dp) span(idimension) = 1.0_dp
      END DO

      DO iatom = 1, natom
         DO idimension = 1, 3
            integer_coordinate(idimension) = &
               INT(((particle_set(iatom)%r(idimension) - lo(idimension))/span(idimension))* &
                   REAL(cmax, dp), int_8)
            integer_coordinate(idimension) = &
               MIN(cmax, MAX(0_int_8, integer_coordinate(idimension)))
         END DO
         CALL morton_split3(integer_coordinate(1), m1)
         CALL morton_split3(integer_coordinate(2), m2)
         CALL morton_split3(integer_coordinate(3), m3)
         morton_code(iatom) = IOR(IOR(m1, ISHFT(m2, 1)), ISHFT(m3, 2))
      END DO

      CALL sort(morton_code, natom, order)
      DEALLOCATE (morton_code)

      CALL timestop(handle)

   END SUBROUTINE spatial_atom_order

! **************************************************************************************************
!> \brief Spread the low 21 bits of an integer over every third bit of a Morton code.
!> \param input_integer ...
!> \param spread_integer ...
! **************************************************************************************************
   SUBROUTINE morton_split3(input_integer, spread_integer)
      INTEGER(KIND=int_8), INTENT(IN)                    :: input_integer
      INTEGER(KIND=int_8), INTENT(OUT)                   :: spread_integer

      spread_integer = IAND(input_integer, INT(z'1FFFFF', int_8))
      spread_integer = IAND(IOR(spread_integer, ISHFT(spread_integer, 32)), &
                            INT(z'1F00000000FFFF', int_8))
      spread_integer = IAND(IOR(spread_integer, ISHFT(spread_integer, 16)), &
                            INT(z'1F0000FF0000FF', int_8))
      spread_integer = IAND(IOR(spread_integer, ISHFT(spread_integer, 8)), &
                            INT(z'100F00F00F00F00F', int_8))
      spread_integer = IAND(IOR(spread_integer, ISHFT(spread_integer, 4)), &
                            INT(z'10C30C30C30C30C3', int_8))
      spread_integer = IAND(IOR(spread_integer, ISHFT(spread_integer, 2)), &
                            INT(z'1249249249249249', int_8))
   END SUBROUTINE morton_split3

END MODULE gw_ri_rs_grid_setup_main
