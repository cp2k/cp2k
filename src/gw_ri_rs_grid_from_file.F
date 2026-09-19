!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Read pretabulated atom-relative RI-RS grids from data files.
!> \par History
!>      09.2026 created
! **************************************************************************************************
MODULE gw_ri_rs_grid_from_file
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE particle_types,                  ONLY: particle_type
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type,&
                                              rirs_grid_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: read_ri_rs_grid_from_file

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_ri_rs_grid_from_file'

CONTAINS

! **************************************************************************************************
!> \brief Read pretabulated atom-relative RI-RS grids from the CP2K data files.
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE read_ri_rs_grid_from_file(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'read_ri_rs_grid_from_file'

      CHARACTER(LEN=default_path_length)                 :: filepath
      CHARACTER(LEN=default_string_length)               :: element_symbol
      INTEGER                                            :: handle, iatom, ikind, natom, nkind
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(rirs_grid_type), ALLOCATABLE                  :: kind_grids(:)

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(bs_env%ri_rs%atomic_kind_set))
      nkind = SIZE(bs_env%ri_rs%atomic_kind_set)
      particle_set => bs_env%ri_rs%particle_set
      CPASSERT(ASSOCIATED(particle_set))
      natom = bs_env%n_atom
      ALLOCATE (kind_grids(nkind))

      DO ikind = 1, nkind
         element_symbol = TRIM(bs_env%ri_rs%atomic_kind_set(ikind)%element_symbol)
         CALL get_ri_rs_grid_filepath(element_symbol, bs_env%ri_rs%grid_select, &
                                      bs_env%ri_rs%grid_file_suffix, filepath)
         CALL read_ri_rs_grid_file(TRIM(filepath), kind_grids(ikind)%raw_points)
         kind_grids(ikind)%npts = SIZE(kind_grids(ikind)%raw_points, 2)
      END DO

      IF (ALLOCATED(bs_env%ri_rs%atomic_grids)) DEALLOCATE (bs_env%ri_rs%atomic_grids)
      ALLOCATE (bs_env%ri_rs%atomic_grids(natom))
      DO iatom = 1, natom
         ikind = particle_set(iatom)%atomic_kind%kind_number
         bs_env%ri_rs%atomic_grids(iatom) = kind_grids(ikind)
      END DO
      DEALLOCATE (kind_grids)

      CALL timestop(handle)

   END SUBROUTINE read_ri_rs_grid_from_file

! **************************************************************************************************
!> \brief Construct the path of a pretabulated RI-RS grid file.
!> \param element_symbol ...
!> \param grid_select ...
!> \param grid_file_suffix ...
!> \param filepath ...
! **************************************************************************************************
   SUBROUTINE get_ri_rs_grid_filepath(element_symbol, grid_select, grid_file_suffix, filepath)
      CHARACTER(LEN=*), INTENT(IN)                       :: element_symbol
      INTEGER, INTENT(IN)                                :: grid_select
      CHARACTER(LEN=*), INTENT(IN)                       :: grid_file_suffix
      CHARACTER(LEN=*), INTENT(OUT)                      :: filepath

      CHARACTER(LEN=default_string_length)               :: suffix

      SELECT CASE (grid_select)
      CASE (1)
         suffix = "_def2-tzvp-rs.ion"
      CASE (2)
         suffix = "_cc-pvtz-rs.ion"
      CASE (3)
         IF (LEN_TRIM(grid_file_suffix) > 0) THEN
            suffix = TRIM(grid_file_suffix)
         ELSE
            suffix = "_rirs.ion"
         END IF
      CASE DEFAULT
         CPABORT("Unknown grid_select (1=def2-TZVPP, 2=cc-pVTZ, 3=user-provided).")
      END SELECT
      filepath = "ri_rs_grid/"//TRIM(element_symbol)//TRIM(suffix)
   END SUBROUTINE get_ri_rs_grid_filepath

! **************************************************************************************************
!> \brief Read Cartesian grid points from the existing CP2K RI-RS .ion format.
!> \param filename Complete input filename; no suffix or element-name construction is performed.
!> \param points Grid points in Bohr, indexed (alpha,l).
! **************************************************************************************************
   SUBROUTINE read_ri_rs_grid_file(filename, points)
      CHARACTER(LEN=*), INTENT(IN)                       :: filename
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: points

      CHARACTER(LEN=*), PARAMETER :: routineN = 'read_ri_rs_grid_file'

      CHARACTER(LEN=default_path_length)                 :: line
      INTEGER                                            :: handle, ierr, ipoint, iunit, npoints
      LOGICAL                                            :: found_points

      CALL timeset(routineN, handle)
      CALL open_file(file_name=TRIM(filename), unit_number=iunit, &
                     file_action='READ', file_status='OLD')

      CALL read_ri_rs_grid_header(iunit, filename, npoints)

      REWIND (iunit)
      found_points = .FALSE.
      DO
         READ (iunit, '(A)', IOSTAT=ierr) line
         IF (ierr /= 0) EXIT
         IF (INDEX(line, '<grid_points>') > 0) THEN
            found_points = .TRUE.
            EXIT
         END IF
      END DO
      IF (.NOT. found_points) THEN
         CPABORT('RI-RS .ion file has no <grid_points> block: '//TRIM(filename))
      END IF

      ALLOCATE (points(3, npoints))
      DO ipoint = 1, npoints
         READ (iunit, *, IOSTAT=ierr) points(:, ipoint)
         IF (ierr /= 0) CPABORT('Invalid grid point in RI-RS .ion file: '//TRIM(filename))
      END DO
      CALL close_file(unit_number=iunit)
      CALL timestop(handle)
   END SUBROUTINE read_ri_rs_grid_file

! **************************************************************************************************
!> \brief Read and validate the point count in an open RI-RS grid file.
!> \param iunit ...
!> \param filename ...
!> \param npoints ...
! **************************************************************************************************
   SUBROUTINE read_ri_rs_grid_header(iunit, filename, npoints)
      INTEGER, INTENT(IN)                                :: iunit
      CHARACTER(LEN=*), INTENT(IN)                       :: filename
      INTEGER, INTENT(OUT)                               :: npoints

      CHARACTER(LEN=default_path_length)                 :: line
      INTEGER                                            :: colon, ierr
      LOGICAL                                            :: found_size

      found_size = .FALSE.
      npoints = 0
      DO
         READ (iunit, '(A)', IOSTAT=ierr) line
         IF (ierr /= 0) EXIT
         IF (INDEX(line, 'n points') > 0) THEN
            colon = INDEX(line, ':')
            IF (colon > 0) THEN
               READ (line(colon + 1:), *, IOSTAT=ierr) npoints
               found_size = ierr == 0 .AND. npoints > 0
            END IF
            EXIT
         END IF
      END DO
      IF (.NOT. found_size) THEN
         CPABORT('RI-RS .ion file has no valid n points field: '//TRIM(filename))
      END IF
   END SUBROUTINE read_ri_rs_grid_header

END MODULE gw_ri_rs_grid_from_file
