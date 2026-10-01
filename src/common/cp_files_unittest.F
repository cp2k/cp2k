!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM cp_files_unittest
   USE cp_files,                        ONLY: close_file,&
                                              get_unit_number,&
                                              init_preconnection_list,&
                                              open_file
   USE machine,                         ONLY: default_input_unit,&
                                              default_output_unit

   IMPLICIT NONE
   INTEGER :: u, v, saved, value, ios, i
   INTEGER :: units(16)
   CHARACTER(LEN=64) :: name
   LOGICAL :: opened, named, exists

   CALL init_preconnection_list()
   CALL open_file('cp_files_test.dat', unit_number=u, file_status='REPLACE', file_action='READWRITE')
   IF (u <= 0 .OR. u == default_input_unit .OR. u == default_output_unit) ERROR STOP 'Invalid unit'
   WRITE (u, *) 11
   saved = u
   CALL close_file(u, keep_preconnection=.TRUE.)
   INQUIRE (UNIT=u, OPENED=opened)
   IF (.NOT. opened) ERROR STOP 'Kept connection closed'
   IF (get_unit_number('cp_files_test.dat') /= u) ERROR STOP 'Named lookup failed'
   CALL open_file('cp_files_other.dat', unit_number=v, file_status='REPLACE', file_action='READWRITE')
   IF (v == saved) ERROR STOP 'Kept unit reused'
   CALL open_file('cp_files_test.dat', unit_number=u, file_action='READWRITE', file_position='APPEND')
   IF (u /= saved) ERROR STOP 'Reopen changed unit'
   WRITE (u, *) 22
   CALL close_file(u, keep_preconnection=.TRUE.)
   CALL open_file('cp_files_test.dat', unit_number=u, file_action='READWRITE')
   READ (u, *) value
   IF (value /= 11) ERROR STOP 'Reopen did not rewind'
   READ (u, *) value
   IF (value /= 22) ERROR STOP 'Append lost data'
   CALL close_file(u)
   INQUIRE (UNIT=u, OPENED=opened)
   IF (opened) ERROR STOP 'Final close failed'
   CALL open_file('cp_files_test.dat', unit_number=u, file_action='READWRITE', file_position='APPEND')
   WRITE (u, *) 33
   REWIND (u)
   READ (u, *) value
   READ (u, *) value
   READ (u, *) value
   IF (value /= 33) ERROR STOP 'Closed-file append failed'
   CALL close_file(u, file_status='DELETE')
   INQUIRE (FILE='cp_files_test.dat', EXIST=exists)
   IF (exists) ERROR STOP 'Delete failed'
   CALL close_file(v, file_status='DELETE')

   CALL open_file('', unit_number=u, file_status='SCRATCH', file_action='READWRITE')
   INQUIRE (UNIT=u, NAMED=named)
   IF (named) ERROR STOP 'Scratch is named'
   WRITE (u, *) 44
   REWIND (u)
   READ (u, *) value
   IF (value /= 44) ERROR STOP 'Scratch roundtrip failed'
   CALL close_file(u)
   INQUIRE (UNIT=u, OPENED=opened)
   IF (opened) ERROR STOP 'Scratch close failed'

   u = 91
   CALL open_file('cp_files_explicit.dat', unit_number=u, file_status='REPLACE', &
                  file_action='WRITE', skip_get_unit_number=.TRUE.)
   IF (u /= 91) ERROR STOP 'Explicit unit changed'
   CALL close_file(u, file_status='DELETE')

   ! Binary restart readers resume sequential records without reconnecting or rewinding.
   CALL open_file('cp_files_binary.dat', unit_number=u, file_status='REPLACE', &
                  file_form='UNFORMATTED', file_action='READWRITE')
   WRITE (u) 11
   WRITE (u) 22
   REWIND (u)
   READ (u) value
   IF (value /= 11) ERROR STOP 'Binary first record failed'
   saved = u
   CALL close_file(u, keep_preconnection=.TRUE.)
   CALL open_file('cp_files_binary.dat', unit_number=u, file_form='UNFORMATTED', &
                  file_action='READWRITE', file_position='ASIS')
   IF (u /= saved) ERROR STOP 'Binary reopen changed unit'
   READ (u) value
   IF (value /= 22) ERROR STOP 'Binary ASIS lost position'
   CALL close_file(u, file_status='DELETE')

   DO i = 1, SIZE(units)
      WRITE (name, '(A,I0)') 'cp_files_multi_', i
      CALL open_file(TRIM(name), unit_number=units(i), file_status='REPLACE', file_action='WRITE')
      IF (ANY(units(:i - 1) == units(i))) ERROR STOP 'Multiple files share unit'
   END DO
   DO i = 1, SIZE(units)
      CALL close_file(units(i), file_status='DELETE')
      INQUIRE (UNIT=units(i), OPENED=opened, IOSTAT=ios)
      IF (ios /= 0 .OR. opened) ERROR STOP 'Multiple file close failed'
   END DO
   ! Exceed the former ten-slot registry (two slots were reserved for standard I/O).
   DO i = 1, SIZE(units)
      WRITE (name, '(A,I0)') 'cp_files_kept_', i
      CALL open_file(TRIM(name), unit_number=units(i), file_status='REPLACE', file_action='WRITE')
      IF (get_unit_number(TRIM(name)) /= units(i)) ERROR STOP 'Open-file lookup failed'
      CALL close_file(units(i), keep_preconnection=.TRUE.)
      IF (ANY(units(:i - 1) == units(i))) ERROR STOP 'Kept files share unit'
   END DO
   CALL init_preconnection_list()
   DO i = 1, SIZE(units)
      WRITE (name, '(A,I0)') 'cp_files_kept_', i
      CALL open_file(TRIM(name), unit_number=u, file_action='WRITE')
      IF (u /= units(i)) ERROR STOP 'Kept-file lookup failed'
      CALL close_file(u, file_status='DELETE')
      INQUIRE (UNIT=u, OPENED=opened)
      IF (opened) ERROR STOP 'Kept-file release failed'
   END DO
   PRINT *, 'cp_files_unittest: PASS'
END PROGRAM cp_files_unittest
