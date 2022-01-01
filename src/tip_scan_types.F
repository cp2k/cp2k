!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE tip_scan_types
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE pw_types,                        ONLY: pw_release,&
                                              pw_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tip_scan_types'

   PUBLIC :: scanning_type, release_scanning_type, read_scanning_section

! **************************************************************************************************
   TYPE scanning_type
      INTEGER                              :: num_scan_points
      REAL(KIND=dp), DIMENSION(3)          :: ref_point
      REAL(KIND=dp), DIMENSION(:, :), &
         ALLOCATABLE                       :: tip_pos
      CHARACTER(LEN=default_string_length) :: tip_cube_file
      TYPE(pw_type), POINTER               :: tip_pw_r
      TYPE(pw_type), POINTER               :: tip_pw_g
   END TYPE scanning_type
! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param scan_info ...
!> \param input_section ...
! **************************************************************************************************
   SUBROUTINE read_scanning_section(scan_info, input_section)

      TYPE(scanning_type), INTENT(INOUT)                 :: scan_info
      TYPE(section_vals_type), POINTER                   :: input_section

      CHARACTER(LEN=default_string_length)               :: schar
      INTEGER                                            :: ii, ix, iy, iz, nx, ny, nz
      INTEGER, DIMENSION(:), POINTER                     :: ilist
      REAL(KIND=dp)                                      :: dx, dy, dz
      REAL(KIND=dp), DIMENSION(:), POINTER               :: rlist, rpoint

      CALL section_vals_val_get(input_section, "SCAN_DIRECTION", c_val=schar)
      CALL section_vals_val_get(input_section, "REFERENCE_POINT", r_vals=rpoint)
      CALL section_vals_val_get(input_section, "SCAN_POINTS", i_vals=ilist)
      CALL section_vals_val_get(input_section, "SCAN_STEP", r_vals=rlist)

      nx = 1
      ny = 1
      nz = 1

      dx = 0.0_dp
      dy = 0.0_dp
      dz = 0.0_dp

      SELECT CASE (schar)
      CASE ("X")
         CPASSERT(SIZE(ilist) >= 1)
         CPASSERT(SIZE(rlist) >= 1)
         nx = ilist(1)
         dx = rlist(1)
      CASE ("Y")
         CPASSERT(SIZE(ilist) >= 1)
         CPASSERT(SIZE(rlist) >= 1)
         ny = ilist(1)
         dy = rlist(1)
      CASE ("Z")
         CPASSERT(SIZE(ilist) >= 1)
         CPASSERT(SIZE(rlist) >= 1)
         nz = ilist(1)
         dz = rlist(1)
      CASE ("XY")
         CPASSERT(SIZE(ilist) >= 2)
         CPASSERT(SIZE(rlist) >= 2)
         nx = ilist(1)
         ny = ilist(2)
         dx = rlist(1)
         dy = rlist(2)
      CASE ("XZ")
         CPASSERT(SIZE(ilist) >= 2)
         CPASSERT(SIZE(rlist) >= 2)
         nx = ilist(1)
         nz = ilist(2)
         dx = rlist(1)
         dz = rlist(2)
      CASE ("YZ")
         CPASSERT(SIZE(ilist) >= 2)
         CPASSERT(SIZE(rlist) >= 2)
         ny = ilist(1)
         nz = ilist(2)
         dy = rlist(1)
         dz = rlist(2)
      CASE ("XYZ")
         CPASSERT(SIZE(ilist) >= 3)
         CPASSERT(SIZE(rlist) >= 3)
         nx = ilist(1)
         ny = ilist(2)
         nz = ilist(3)
         dx = rlist(1)
         dy = rlist(2)
         dz = rlist(3)
      CASE DEFAULT
         CPABORT("Invalid Scan Type")
      END SELECT

      scan_info%ref_point(1:3) = rpoint(1:3)
      scan_info%num_scan_points = nx*ny*nz
      ALLOCATE (scan_info%tip_pos(3, nx*ny*nz))
      rpoint(1) = rpoint(1) - 0.5_dp*(nx - 1)*dx
      rpoint(2) = rpoint(2) - 0.5_dp*(ny - 1)*dy
      rpoint(3) = rpoint(3) - 0.5_dp*(nz - 1)*dz

      ii = 0
      DO iz = 1, nz
         DO iy = 1, ny
            DO ix = 1, nx
               ii = ii + 1
               scan_info%tip_pos(1, ii) = rpoint(1) + (ix - 1)*dx
               scan_info%tip_pos(2, ii) = rpoint(2) + (iy - 1)*dy
               scan_info%tip_pos(3, ii) = rpoint(3) + (iz - 1)*dz
            END DO
         END DO
      END DO

      ! tip potential file name
      CALL section_vals_val_get(input_section, "TIP_FILENAME", c_val=schar)
      scan_info%tip_cube_file = schar

      NULLIFY (scan_info%tip_pw_r)
      NULLIFY (scan_info%tip_pw_g)

   END SUBROUTINE read_scanning_section

! **************************************************************************************************
!> \brief ...
!> \param scan_info ...
! **************************************************************************************************
   SUBROUTINE release_scanning_type(scan_info)

      TYPE(scanning_type), INTENT(INOUT)                 :: scan_info

      scan_info%num_scan_points = 0
      scan_info%ref_point = 0.0_dp
      IF (ALLOCATED(scan_info%tip_pos)) THEN
         DEALLOCATE (scan_info%tip_pos)
      END IF
      IF (ASSOCIATED(scan_info%tip_pw_r)) THEN
         CALL pw_release(scan_info%tip_pw_r)
      END IF
      IF (ASSOCIATED(scan_info%tip_pw_g)) THEN
         CALL pw_release(scan_info%tip_pw_g)
      END IF

   END SUBROUTINE release_scanning_type

END MODULE tip_scan_types
