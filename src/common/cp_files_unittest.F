!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM cp_files_unittest
   USE base_hooks,                      ONLY: cp_abort_hook
   USE cp_files,                        ONLY: close_file,&
                                              get_unit_number,&
                                              init_preconnection_list,&
                                              open_file
   USE machine,                         ONLY: default_input_unit,&
                                              default_output_unit
   USE message_passing,                 ONLY: mp_abort,&
                                              mp_comm_type,&
                                              mp_world_finalize,&
                                              mp_world_init
#include "../base/base_uses.f90"

   IMPLICIT NONE

   TYPE(mp_comm_type) :: mpi_comm

   CALL mp_world_init(mpi_comm)
   cp_abort_hook => abort_test
   ! All ranks share the test directory; only the source rank may touch named files.
   IF (mpi_comm%is_source()) CALL test_files()
   CALL mp_world_finalize()
   NULLIFY (cp_abort_hook)

CONTAINS

! **************************************************************************************************
!> \brief Report a test failure and terminate all ranks, including those waiting in finalization.
!> \param location Source location of the failure
!> \param message Failure description
! **************************************************************************************************
   SUBROUTINE abort_test(location, message)
      CHARACTER(LEN=*), INTENT(IN)                       :: location, message

      WRITE (UNIT=default_output_unit, FMT='(A)') TRIM(location)//': '//TRIM(message)
      FLUSH (default_output_unit)
      CALL mp_abort()
   END SUBROUTINE abort_test

! **************************************************************************************************
!> \brief Exercise file lifetimes, retained connections and file positions on one MPI rank.
! **************************************************************************************************
   SUBROUTINE test_files()
      CHARACTER(LEN=64)                                  :: name
      INTEGER                                            :: i, ios, saved, u, units(16), v, value
      LOGICAL                                            :: exists, named, opened

      CALL init_preconnection_list()
      CALL open_file('cp_files_test.dat', unit_number=u, file_status='REPLACE', file_action='READWRITE')
      IF (u <= 0 .OR. u == default_input_unit .OR. u == default_output_unit) THEN
         CPABORT('Invalid unit')
      END IF
      WRITE (u, *) 11
      saved = u
      CALL close_file(u, keep_preconnection=.TRUE.)
      INQUIRE (UNIT=u, OPENED=opened)
      IF (.NOT. opened) CPABORT('Kept connection closed')
      IF (get_unit_number('cp_files_test.dat') /= u) CPABORT('Named lookup failed')
      CALL open_file('cp_files_other.dat', unit_number=v, file_status='REPLACE', file_action='READWRITE')
      IF (v == saved) CPABORT('Kept unit reused')
      CALL open_file('cp_files_test.dat', unit_number=u, file_action='READWRITE', file_position='APPEND')
      IF (u /= saved) CPABORT('Reopen changed unit')
      WRITE (u, *) 22
      CALL close_file(u, keep_preconnection=.TRUE.)
      CALL open_file('cp_files_test.dat', unit_number=u, file_action='READWRITE')
      READ (u, *) value
      IF (value /= 11) CPABORT('Reopen did not rewind')
      READ (u, *) value
      IF (value /= 22) CPABORT('Append lost data')
      CALL close_file(u)
      INQUIRE (UNIT=u, OPENED=opened)
      IF (opened) CPABORT('Final close failed')
      CALL open_file('cp_files_test.dat', unit_number=u, file_action='READWRITE', file_position='APPEND')
      WRITE (u, *) 33
      REWIND (u)
      READ (u, *) value
      READ (u, *) value
      READ (u, *) value
      IF (value /= 33) CPABORT('Closed-file append failed')
      CALL close_file(u, file_status='DELETE')
      INQUIRE (FILE='cp_files_test.dat', EXIST=exists)
      IF (exists) CPABORT('Delete failed')
      CALL close_file(v, file_status='DELETE')

      CALL open_file('', unit_number=u, file_status='SCRATCH', file_action='READWRITE')
      INQUIRE (UNIT=u, NAMED=named)
      IF (named) CPABORT('Scratch is named')
      WRITE (u, *) 44
      REWIND (u)
      READ (u, *) value
      IF (value /= 44) CPABORT('Scratch roundtrip failed')
      CALL close_file(u)
      INQUIRE (UNIT=u, OPENED=opened)
      IF (opened) CPABORT('Scratch close failed')

      u = 91
      CALL open_file('cp_files_explicit.dat', unit_number=u, file_status='REPLACE', &
                     file_action='WRITE', skip_get_unit_number=.TRUE.)
      IF (u /= 91) CPABORT('Explicit unit changed')
      CALL close_file(u, file_status='DELETE')

      ! Binary restart readers resume sequential records without reconnecting or rewinding.
      CALL open_file('cp_files_binary.dat', unit_number=u, file_status='REPLACE', &
                     file_form='UNFORMATTED', file_action='READWRITE')
      WRITE (u) 11
      WRITE (u) 22
      REWIND (u)
      READ (u) value
      IF (value /= 11) CPABORT('Binary first record failed')
      saved = u
      CALL close_file(u, keep_preconnection=.TRUE.)
      CALL open_file('cp_files_binary.dat', unit_number=u, file_form='UNFORMATTED', &
                     file_action='READWRITE', file_position='ASIS')
      IF (u /= saved) CPABORT('Binary reopen changed unit')
      READ (u) value
      IF (value /= 22) CPABORT('Binary ASIS lost position')
      CALL close_file(u, file_status='DELETE')

      DO i = 1, SIZE(units)
         WRITE (name, '(A,I0)') 'cp_files_multi_', i
         CALL open_file(TRIM(name), unit_number=units(i), file_status='REPLACE', file_action='WRITE')
         IF (ANY(units(:i - 1) == units(i))) CPABORT('Multiple files share unit')
      END DO
      DO i = 1, SIZE(units)
         CALL close_file(units(i), file_status='DELETE')
         INQUIRE (UNIT=units(i), OPENED=opened, IOSTAT=ios)
         IF (ios /= 0 .OR. opened) CPABORT('Multiple file close failed')
      END DO
      ! Exceed the former ten-slot registry (two slots were reserved for standard I/O).
      DO i = 1, SIZE(units)
         WRITE (name, '(A,I0)') 'cp_files_kept_', i
         CALL open_file(TRIM(name), unit_number=units(i), file_status='REPLACE', file_action='WRITE')
         IF (get_unit_number(TRIM(name)) /= units(i)) CPABORT('Open-file lookup failed')
         CALL close_file(units(i), keep_preconnection=.TRUE.)
         IF (ANY(units(:i - 1) == units(i))) CPABORT('Kept files share unit')
      END DO
      CALL init_preconnection_list()
      DO i = 1, SIZE(units)
         WRITE (name, '(A,I0)') 'cp_files_kept_', i
         CALL open_file(TRIM(name), unit_number=u, file_action='WRITE')
         IF (u /= units(i)) CPABORT('Kept-file lookup failed')
         CALL close_file(u, file_status='DELETE')
         INQUIRE (UNIT=u, OPENED=opened)
         IF (opened) CPABORT('Kept-file release failed')
      END DO
      WRITE (UNIT=default_output_unit, FMT='(A)') 'cp_files_unittest: PASS'
   END SUBROUTINE test_files
END PROGRAM cp_files_unittest
