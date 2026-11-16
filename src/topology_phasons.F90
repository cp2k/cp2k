!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Native, bounded-memory snapshot/link postprocessor for closed 2D/4D parameter meshes.
!> \note No SCF solver or external topology package is invoked. Run one process with OpenMP threads.
! **************************************************************************************************
PROGRAM topology_phasons
   USE iso_fortran_env,                 ONLY: output_unit
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE topology_phason_analysis,        ONLY: analyze_mesh,&
                                              compare_refinements,&
                                              phason_result_type
   USE topology_snapshot,               ONLY: check_snapshot,&
                                              read_snapshot,&
                                              snapshot_overlap,&
                                              snapshot_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   CHARACTER(LEN=default_path_length) :: mode, filename, arg, filename_b
   TYPE(snapshot_type) :: left, right
   TYPE(phason_result_type) :: coarse, fine
   COMPLEX(KIND=dp), ALLOCATABLE :: overlap(:, :)
   INTEGER :: ia, ib, shift(3), i, j, ios, io_unit
   REAL(KIND=dp) :: metric, gap

   io_unit = output_unit
   CALL get_command_argument(1, mode)
   SELECT CASE (TRIM(mode))
   CASE ('overlap')
      IF (COMMAND_ARGUMENT_COUNT() /= 8) THEN
         CPABORT('Usage: topology_phasons overlap left point right point gx gy gz')
      END IF
      CALL get_command_argument(2, filename)
      CALL get_command_argument(3, arg)
      READ (arg, *, IOSTAT=ios) ia
      CPASSERT(ios == 0)
      CALL get_command_argument(4, filename_b)
      CALL get_command_argument(5, arg)
      READ (arg, *, IOSTAT=ios) ib
      CPASSERT(ios == 0)
      DO i = 1, 3
         CALL get_command_argument(5 + i, arg)
         READ (arg, *, IOSTAT=ios) shift(i)
         CPASSERT(ios == 0)
      END DO
      CALL read_snapshot(TRIM(filename), ia, left)
      CALL read_snapshot(TRIM(filename_b), ib, right)
      CALL check_snapshot(left, 1.e-7_dp, 1.e-7_dp, metric, gap)
      CALL check_snapshot(right, 1.e-7_dp, 1.e-7_dp, metric, gap)
      ALLOCATE (overlap(left%rank, right%rank))
      CALL snapshot_overlap(left, right, left%k, right%k + REAL(shift, dp), overlap)
      WRITE (io_unit, '(A,I0)') 'OVERLAP ', left%rank
      DO j = 1, right%rank
         DO i = 1, left%rank
            WRITE (io_unit, '(2ES26.17)') REAL(overlap(i, j), dp), AIMAG(overlap(i, j))
         END DO
      END DO
   CASE ('mesh', 'links')
      IF (COMMAND_ARGUMENT_COUNT() /= 2) THEN
         CPABORT('Usage: topology_phasons mesh|links manifest')
      END IF
      CALL get_command_argument(2, filename)
      CALL analyze_mesh(TRIM(filename), TRIM(mode) == 'mesh', fine, io_unit)
   CASE ('converge', 'converge-links')
      IF (COMMAND_ARGUMENT_COUNT() /= 4) THEN
         CPABORT('Usage: topology_phasons converge coarse fine tolerance')
      END IF
      CALL get_command_argument(2, filename)
      CALL get_command_argument(3, filename_b)
      CALL get_command_argument(4, arg)
      READ (arg, *, IOSTAT=ios) metric
      CPASSERT(ios == 0)
      CALL analyze_mesh(TRIM(filename), TRIM(mode) == 'converge', coarse, io_unit)
      CALL analyze_mesh(TRIM(filename_b), TRIM(mode) == 'converge', fine, io_unit)
      CALL compare_refinements(coarse, fine, metric, io_unit)
   CASE DEFAULT
      WRITE (io_unit, '(A)') 'Usage: topology_phasons mesh manifest | links file | overlap left point right point gx gy gz'
      CPABORT('Unknown or missing topology_phasons mode')
   END SELECT
END PROGRAM topology_phasons
