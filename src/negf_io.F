!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for reading and writing NEGF restart files.
!> \author Dmitry Ryndyk (12.2025)
! **************************************************************************************************
MODULE negf_io

   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_log_handling,                 ONLY: cp_logger_generate_filename,&
                                              cp_logger_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'negf_io'

   PUBLIC :: negf_restart_file_name, &
             negf_print_matrix_to_file, &
             negf_read_matrix_from_file

CONTAINS

! **************************************************************************************************
!> \brief Checks if the restart file exists and returns the filename.
!> \param filename  ...
!> \param exist     ...
!> \param negf_section ...
!> \param logger    ...
!> \param icontact  ...
!> \param ispin     ...
!> \param h00 ...
!> \param h01 ...
!> \param s00 ...
!> \param s01 ...
!> \par History
!>    * 12.2025 created  [Dmitry Ryndyk]
! **************************************************************************************************
   SUBROUTINE negf_restart_file_name(filename, exist, negf_section, logger, icontact, ispin, h00, h01, s00, s01)
      CHARACTER(LEN=default_path_length), INTENT(OUT)    :: filename
      LOGICAL, INTENT(OUT)                               :: exist
      TYPE(section_vals_type), POINTER                   :: negf_section
      TYPE(cp_logger_type), POINTER                      :: logger
      INTEGER, INTENT(IN), OPTIONAL                      :: icontact, ispin
      LOGICAL, INTENT(IN), OPTIONAL                      :: h00, h01, s00, s01

      CHARACTER(len=default_string_length)               :: middle_name, string1, string2
      LOGICAL                                            :: my_h00, my_h01, my_s00, my_s01
      TYPE(section_vals_type), POINTER                   :: contact_section, print_key

      my_h00 = .FALSE.
      IF (PRESENT(h00)) my_h00 = h00
      my_h01 = .FALSE.
      IF (PRESENT(h01)) my_h01 = h01
      my_s00 = .FALSE.
      IF (PRESENT(s00)) my_s00 = s00
      my_s01 = .FALSE.
      IF (PRESENT(s01)) my_s01 = s01

      exist = .FALSE.

      ! try to read from the filename that is generated automatically from the printkey
      contact_section => section_vals_get_subs_vals(negf_section, "CONTACT")
      print_key => section_vals_get_subs_vals(contact_section, "RESTART", i_rep_section=icontact)

      IF (my_h00) THEN

         WRITE (string1, *) icontact
         WRITE (string2, *) ispin
         IF (ispin == 0) THEN
            middle_name = "N"//TRIM(string1)//"-H00"
         ELSE
            middle_name = "N"//TRIM(string1)//"-H00-S"//TRIM(string2)
         END IF
         filename = negf_elecrodes_generate_filename(logger, print_key, middle_name=middle_name, &
                                                     extension=".hs", my_local=.FALSE.)
      END IF

      IF (my_h01) THEN
         WRITE (string1, *) icontact
         WRITE (string2, *) ispin
         IF (ispin == 0) THEN
            middle_name = "N"//TRIM(string1)//"-H01"
         ELSE
            middle_name = "N"//TRIM(string1)//"-H01-S"//TRIM(string2)
         END IF
         filename = negf_elecrodes_generate_filename(logger, print_key, middle_name=middle_name, &
                                                     extension=".hs", my_local=.FALSE.)
      END IF

      IF (my_s00) THEN
         WRITE (string1, *) icontact
         middle_name = "N"//TRIM(string1)//"-S00"
         filename = negf_elecrodes_generate_filename(logger, print_key, middle_name=middle_name, &
                                                     extension=".hs", my_local=.FALSE.)
      END IF

      IF (my_s01) THEN
         WRITE (string1, *) icontact
         middle_name = "N"//TRIM(string1)//"-S01"
         filename = negf_elecrodes_generate_filename(logger, print_key, middle_name=middle_name, &
                                                     extension=".hs", my_local=.FALSE.)
      END IF

      INQUIRE (FILE=filename, exist=exist)

   END SUBROUTINE negf_restart_file_name

! **************************************************************************************************
!> \brief ...
!> \param logger the logger for the parallel environment, iteration info
!>        and filename generation
!> \param print_key ...
!> \param middle_name name to be added to the generated filename, useful when
!>        print_key activates different distinct outputs, to be able to
!>        distinguish them
!> \param extension extension to be applied to the filename (including the ".")
!> \param my_local if the unit should be local to this task, or global to the
!>        program (defaults to false).
!> \return ...
!> \par History
!>    * 12.2025 created  [Dmitry Ryndyk]
! **************************************************************************************************
   FUNCTION negf_elecrodes_generate_filename(logger, print_key, middle_name, extension, &
                                             my_local) RESULT(filename)
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_key
      CHARACTER(len=*), INTENT(IN), OPTIONAL             :: middle_name
      CHARACTER(len=*), INTENT(IN)                       :: extension
      LOGICAL, INTENT(IN)                                :: my_local
      CHARACTER(len=default_path_length)                 :: filename

      CHARACTER(len=default_path_length)                 :: outName, outPath, postfix, root
      CHARACTER(len=default_string_length)               :: my_middle_name
      INTEGER                                            :: my_ind1, my_ind2
      LOGICAL                                            :: has_root

      CALL section_vals_val_get(print_key, "FILENAME", c_val=outPath)
      IF (outPath(1:1) == '=') THEN
         CPASSERT(LEN(outPath) - 1 <= LEN(filename))
         filename = outPath(2:)
         RETURN
      END IF
      IF (outPath == "__STD_OUT__") outPath = ""
      outName = outPath
      has_root = .FALSE.
      my_ind1 = INDEX(outPath, "/")
      my_ind2 = LEN_TRIM(outPath)
      IF (my_ind1 /= 0) THEN
         has_root = .TRUE.
         DO WHILE (INDEX(outPath(my_ind1 + 1:my_ind2), "/") /= 0)
            my_ind1 = INDEX(outPath(my_ind1 + 1:my_ind2), "/") + my_ind1
         END DO
         IF (my_ind1 == my_ind2) THEN
            outName = ""
         ELSE
            outName = outPath(my_ind1 + 1:my_ind2)
         END IF
      END IF

      IF (PRESENT(middle_name)) THEN
         IF (outName /= "") THEN
            my_middle_name = "-"//TRIM(outName)//"-"//middle_name
         ELSE
            my_middle_name = "-"//middle_name
         END IF
      ELSE
         IF (outName /= "") THEN
            my_middle_name = "-"//TRIM(outName)
         ELSE
            my_middle_name = ""
         END IF
      END IF

      IF (.NOT. has_root) THEN
         root = TRIM(logger%iter_info%project_name)//TRIM(my_middle_name)
      ELSE IF (outName == "") THEN
         root = outPath(1:my_ind1)//TRIM(logger%iter_info%project_name)//TRIM(my_middle_name)
      ELSE
         root = outPath(1:my_ind1)//my_middle_name(2:LEN_TRIM(my_middle_name))
      END IF

      postfix = extension

      ! and let the logger generate the filename
      CALL cp_logger_generate_filename(logger, res=filename, &
                                       root=root, postfix=postfix, local=my_local)

   END FUNCTION negf_elecrodes_generate_filename

! **************************************************************************************************
!> \brief Prints full matrix to a file.
!> \param filename ...
!> \param matrix   ...
!> \par History
!>    * 12.2025 created  [Dmitry Ryndyk]
! **************************************************************************************************
   SUBROUTINE negf_print_matrix_to_file(filename, matrix)
      CHARACTER(LEN=default_path_length), INTENT(IN)     :: filename
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: matrix

      CHARACTER(len=100)                                 :: sfmt
      INTEGER                                            :: i, j, ncol, nrow, print_unit

      CALL open_file(file_name=filename, file_status="REPLACE", &
                     file_form="FORMATTED", file_action="WRITE", &
                     file_position="REWIND", unit_number=print_unit)

      nrow = SIZE(matrix, 1)
      ncol = SIZE(matrix, 2)
      WRITE (sfmt, "('(',i0,'(E15.5))')") ncol
      WRITE (print_unit, *) nrow, ncol
      DO i = 1, nrow
         WRITE (print_unit, sfmt) (matrix(i, j), j=1, ncol)
      END DO

      CALL close_file(print_unit)

   END SUBROUTINE negf_print_matrix_to_file

! **************************************************************************************************
!> \brief Reads full matrix from a file.
!> \param filename ...
!> \param matrix   ...
!> \par History
!>    * 12.2025 created  [Dmitry Ryndyk]
! **************************************************************************************************
   SUBROUTINE negf_read_matrix_from_file(filename, matrix)
      CHARACTER(LEN=default_path_length), INTENT(IN)     :: filename
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: matrix

      INTEGER                                            :: i, j, ncol, nrow, print_unit

      CALL open_file(file_name=filename, file_status="OLD", &
                     file_form="FORMATTED", file_action="READ", &
                     file_position="REWIND", unit_number=print_unit)

      READ (print_unit, *) nrow, ncol
      DO i = 1, nrow
         READ (print_unit, *) (matrix(i, j), j=1, ncol)
      END DO

      CALL close_file(print_unit)

   END SUBROUTINE negf_read_matrix_from_file

END MODULE negf_io
