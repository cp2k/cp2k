!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utility routines to open and close files.
!> \par History
!>      - Creation CP2K_WORKSHOP 1.0 TEAM
!>      - Revised (18.02.2011, MK)
!>      - Enhanced error checking (22.02.2011, MK)
!>      - Preconnected units removed (03.10.2026, MK)
!> \note Positive unit numbers are used instead of OPEN(NEWUNIT=...), because NEWUNIT returns
!>       negative numbers and CP2K treats a unit number > 0 as "output enabled".
!>       open_file is thread-safe: finding a free unit and connecting the file are done within
!>       the OpenMP critical region cp_files_connect. Concurrent use of the same unit (I/O,
!>       REWIND, close_file) must be avoided by the caller.
!> \author Matthias Krack (MK)
! **************************************************************************************************
MODULE cp_files
   USE ISO_C_BINDING,                   ONLY: C_CHAR,&
                                              C_F_POINTER,&
                                              C_NULL_CHAR,&
                                              C_PTR
   USE kinds,                           ONLY: default_path_length
   USE machine,                         ONLY: default_error_unit,&
                                              default_input_unit,&
                                              default_output_unit,&
                                              m_getcwd
   USE string_utilities,                ONLY: uppercase
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: close_file, &
             open_file, &
             file_exists, &
             get_data_dir, &
             discover_file, &
             join_paths

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_files'

   ! Upper bound for the search of a free logical unit number
   INTEGER, PARAMETER :: max_unit_number = 1024

CONTAINS

! **************************************************************************************************
!> \brief Closes an open file given by its logical unit number.
!> \param unit_number Logical unit number of the file to be closed
!> \param file_status File status after closing: KEEP or DELETE (default: KEEP for named,
!>                    DELETE for unnamed SCRATCH files)
!> \author Matthias Krack (MK)
!> \note  The program is aborted, if the unit is one of the default units (stdin, stdout,
!>        stderr). A unit which is not open is ignored silently.
! **************************************************************************************************
   SUBROUTINE close_file(unit_number, file_status)

      INTEGER, INTENT(IN)                                :: unit_number
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: file_status

      CHARACTER(LEN=11)                                  :: status_string
      CHARACTER(LEN=2*default_path_length)               :: message
      CHARACTER(LEN=default_path_length)                 :: iomsg_string
      INTEGER                                            :: istat
      LOGICAL                                            :: exists, is_named, is_open

      INQUIRE (UNIT=unit_number, EXIST=exists, NAMED=is_named, OPENED=is_open, IOSTAT=istat)

      IF (istat /= 0) THEN
         WRITE (UNIT=message, FMT="(A,I0,A,I0,A)") &
            "An error occurred inquiring the unit with the number ", unit_number, &
            " (IOSTAT = ", istat, ")"
         CALL cp_abort(__LOCATION__, TRIM(message))
      ELSE IF (.NOT. exists) THEN
         WRITE (UNIT=message, FMT="(A,I0,A)") &
            "The specified unit number ", unit_number, &
            " cannot be closed, because it does not exist."
         CALL cp_abort(__LOCATION__, TRIM(message))
      END IF

      ! Close the specified file

      IF (is_open) THEN

         ! Refuse to close any default unit
         SELECT CASE (unit_number)
         CASE (default_error_unit)
            WRITE (UNIT=message, FMT="(A,I0)") &
               "Attempt to close the default error unit number ", unit_number
            CALL cp_abort(__LOCATION__, TRIM(message))
         CASE (default_input_unit)
            WRITE (UNIT=message, FMT="(A,I0)") &
               "Attempt to close the default input unit number ", unit_number
            CALL cp_abort(__LOCATION__, TRIM(message))
         CASE (default_output_unit)
            WRITE (UNIT=message, FMT="(A,I0)") &
               "Attempt to close the default output unit number ", unit_number
            CALL cp_abort(__LOCATION__, TRIM(message))
         END SELECT

         ! Define status after closing the file (KEEP for named, DELETE for unnamed SCRATCH files)
         IF (PRESENT(file_status)) THEN
            IF (LEN_TRIM(ADJUSTL(file_status)) > LEN(status_string)) THEN
               CALL cp_abort(__LOCATION__, &
                             "An invalid file STATUS="//TRIM(ADJUSTL(file_status))// &
                             " was specified.")
            END IF
            status_string = TRIM(ADJUSTL(file_status))
            CALL uppercase(status_string)
         ELSE IF (.NOT. is_named) THEN
            status_string = "DELETE"
         ELSE
            status_string = "KEEP"
         END IF

         SELECT CASE (status_string)
         CASE ("DELETE", "KEEP")
            CLOSE (UNIT=unit_number, IOSTAT=istat, IOMSG=iomsg_string, STATUS=TRIM(status_string))
            IF (istat /= 0) THEN
               CALL cp_abort(__LOCATION__, &
                             "An error occurred closing the file with the logical unit number "// &
                             TRIM(int_to_string(unit_number))//" (IOSTAT="// &
                             TRIM(int_to_string(istat))//"): "//TRIM(iomsg_string))
            END IF
         CASE DEFAULT
            CALL cp_abort(__LOCATION__, &
                          "An invalid file STATUS="//TRIM(status_string)//" was specified.")
         END SELECT

      END IF

      ! A unit which is not open is ignored silently. close_file is called by cp_logger_release,
      ! e.g. for a unit shared by several loggers and already closed by one of them. No warning can
      ! be issued at that point, since the logger being released might still be the default logger.

   END SUBROUTINE close_file

! **************************************************************************************************
!> \brief Returns the first free positive logical unit number.
!> \return First free unit in 1..max_unit_number (excluding the default units), or -1 if none
!>         is free
!> \author Matthias Krack (MK)
!> \note  Must be called within the critical region cp_files_connect together with the OPEN
!>        statement, so that no other thread can obtain the same unit number in between.
! **************************************************************************************************
   FUNCTION get_unit_number() RESULT(unit_number)

      INTEGER                                            :: unit_number

      INTEGER                                            :: istat
      LOGICAL                                            :: exists, is_open

      DO unit_number = 1, max_unit_number
         SELECT CASE (unit_number)
         CASE (default_error_unit, default_input_unit, default_output_unit)
            CYCLE
         CASE DEFAULT
            INQUIRE (UNIT=unit_number, EXIST=exists, OPENED=is_open, IOSTAT=istat)
            IF (exists .AND. (.NOT. is_open) .AND. (istat == 0)) RETURN
         END SELECT
      END DO

      unit_number = -1

   END FUNCTION get_unit_number

! **************************************************************************************************
!> \brief Opens a file and returns its logical unit number.
!> \param file_name File name; must be omitted for an anonymous SCRATCH file
!> \param unit_number Logical unit number of the opened file; input if skip_get_unit_number is true
!> \param file_status Fortran file status OLD, NEW, REPLACE, UNKNOWN, or SCRATCH for an
!>                    anonymous temporary file (default: OLD)
!> \param file_form FORMATTED or UNFORMATTED (default: FORMATTED)
!> \param file_action READ, WRITE, or READWRITE (default: WRITE for STATUS=NEW, REPLACE, and
!>                    SCRATCH, otherwise READ)
!> \param file_position ASIS, REWIND, or APPEND (default: REWIND)
!> \param file_pad YES or NO (default: YES); only allowed for FORMATTED files
!> \param file_access SEQUENTIAL or STREAM (default: SEQUENTIAL); DIRECT is not supported
!> \param skip_get_unit_number Use the unit number passed in instead of finding a free one
!> \param debug Print the connection properties to the default output unit
!> \author Matthias Krack (MK)
!> \note  - If the file is already open, its unit number is returned and no new connection is
!>          made. ACTION, FORM, and ACCESS must match the existing connection; with
!>          POSITION=REWIND the file is rewound.
!>        - For STATUS=OLD and ACTION=READ, the file is also looked up in the data directory
!>          (see discover_file). Files opened for writing are never taken from there.
!>        - Any error aborts the program; there is no IOSTAT return.
! **************************************************************************************************
   SUBROUTINE open_file(file_name, unit_number, file_status, file_form, file_action, &
                        file_position, file_pad, file_access, skip_get_unit_number, debug)

      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: file_name
      INTEGER, INTENT(INOUT)                             :: unit_number
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: file_status, file_form, file_action, &
                                                            file_position, file_pad, file_access
      LOGICAL, INTENT(IN), OPTIONAL                      :: skip_get_unit_number, debug

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'open_file'

      CHARACTER(LEN=11) :: access_string, action_string, current_access, current_action, &
         current_form, form_string, pad_string, position_string, status_string
      CHARACTER(LEN=2*default_path_length)               :: message
      CHARACTER(LEN=default_path_length)                 :: cwd, iomsg_string, real_file_name
      INTEGER                                            :: inquire_stat, open_stat, &
                                                            requested_unit_number, unit_stat
      LOGICAL                                            :: connect, exists, is_named, is_open, &
                                                            skip_get, unit_in_use

      real_file_name = ""

      IF (PRESENT(skip_get_unit_number)) THEN
         skip_get = skip_get_unit_number
      ELSE
         skip_get = .FALSE.
      END IF

      IF (.NOT. skip_get) unit_number = -1

      IF (PRESENT(file_access)) THEN
         access_string = TRIM(ADJUSTL(file_access))
         CALL uppercase(access_string)
         IF (access_string == "DIRECT") THEN
            CALL cp_abort(__LOCATION__, &
                          "ACCESS=DIRECT is not supported.")
         END IF
      ELSE
         access_string = "SEQUENTIAL"
      END IF

      IF (PRESENT(file_status)) THEN
         status_string = TRIM(ADJUSTL(file_status))
         CALL uppercase(status_string)
      ELSE
         status_string = "OLD"
      END IF

      IF (PRESENT(file_form)) THEN
         form_string = TRIM(ADJUSTL(file_form))
         CALL uppercase(form_string)
      ELSE
         form_string = "FORMATTED"
      END IF

      IF (PRESENT(file_pad)) THEN
         pad_string = TRIM(ADJUSTL(file_pad))
         CALL uppercase(pad_string)
         IF (form_string == "UNFORMATTED") THEN
            CALL cp_abort(__LOCATION__, &
                          "The PAD specifier is not allowed for an UNFORMATTED file.")
         END IF
      ELSE
         pad_string = "YES"
      END IF

      IF (PRESENT(file_action)) THEN
         action_string = TRIM(ADJUSTL(file_action))
         CALL uppercase(action_string)
      ELSE
         ! A new or emptied file cannot be read, thus WRITE is the only sensible default
         SELECT CASE (status_string)
         CASE ("NEW", "REPLACE", "SCRATCH")
            action_string = "WRITE"
         CASE DEFAULT
            action_string = "READ"
         END SELECT
      END IF

      IF (PRESENT(file_position)) THEN
         position_string = TRIM(ADJUSTL(file_position))
         CALL uppercase(position_string)
      ELSE
         position_string = "REWIND"
      END IF

      IF (PRESENT(file_name)) THEN
         IF (LEN_TRIM(file_name) == 0) THEN
            CALL cp_abort(__LOCATION__, &
                          "The specified FILE name is empty.")
         END IF
         IF (LEN_TRIM(file_name) > LEN(real_file_name)) THEN
            CALL cp_abort(__LOCATION__, &
                          "The FILE name <"//TRIM(file_name)//"> exceeds the maximum length of "// &
                          TRIM(int_to_string(LEN(real_file_name)))//" characters.")
         END IF
         IF (file_name(1:1) == " ") THEN
            CALL cp_abort(__LOCATION__, &
                          "The FILE name <"//TRIM(file_name)//"> has leading blanks.")
         END IF
         IF (status_string == "SCRATCH") THEN
            CALL cp_abort(__LOCATION__, &
                          "A FILE name specifier is not allowed for SCRATCH files.")
         END IF
         ! Resolve the file name. The data directory is only searched for existing files
         ! opened read-only, so that installed data files are never modified.
         IF ((status_string == "OLD") .AND. (action_string == "READ")) THEN
            real_file_name = discover_file(file_name)
         ELSE
            real_file_name = TRIM(ADJUSTL(file_name))
         END IF
      END IF

      SELECT CASE (status_string)
      CASE ("SCRATCH")
         ! No file name required
      CASE ("NEW", "OLD", "REPLACE", "UNKNOWN")
         IF (LEN_TRIM(real_file_name) == 0) THEN
            CALL cp_abort(__LOCATION__, &
                          "A FILE name length of zero for a file with STATUS="// &
                          TRIM(status_string)//" is not allowed.")
         END IF
      CASE DEFAULT
         CALL cp_abort(__LOCATION__, "A file STATUS="//TRIM(status_string)//" is not supported.")
      END SELECT

      ! Inquire the file, find a free unit, and connect the file atomically with respect to other
      ! threads. All checks and error messages follow outside the critical region, because
      ! cp_abort and cp_warn may call open_file again via the logger, which would deadlock.
      requested_unit_number = unit_number
      connect = .FALSE.
      exists = .FALSE.
      is_open = .FALSE.
      unit_in_use = .FALSE.
      inquire_stat = 0
      open_stat = 0
      unit_stat = 0
      iomsg_string = ""
!$OMP CRITICAL(cp_files_connect)
      IF (status_string == "SCRATCH") THEN
         connect = .TRUE.
      ELSE
         INQUIRE (FILE=TRIM(real_file_name), EXIST=exists, OPENED=is_open, IOSTAT=inquire_stat)
         IF (inquire_stat == 0) THEN
            IF (is_open) THEN
               INQUIRE (FILE=TRIM(real_file_name), NUMBER=unit_number, &
                        ACTION=current_action, FORM=current_form, ACCESS=current_access)
            ELSE
               SELECT CASE (status_string)
               CASE ("NEW")
                  connect = (.NOT. exists)
               CASE ("OLD")
                  connect = exists
               CASE DEFAULT
                  connect = .TRUE.
               END SELECT
            END IF
         END IF
      END IF
      IF (connect) THEN
         ! Retrieve unit number automatically if not suppressed explicitly. A requested unit must
         ! not be connected to another file (which OPEN would close silently) except for the
         ! default output unit, which is redirected on purpose (e.g. by the command line flag -o)
         IF (skip_get) THEN
            IF ((unit_number > 0) .AND. (unit_number /= default_output_unit)) THEN
               INQUIRE (UNIT=unit_number, OPENED=unit_in_use, IOSTAT=unit_stat)
            END IF
         ELSE
            unit_number = get_unit_number()
         END IF
         IF ((unit_number > 0) .AND. (.NOT. unit_in_use) .AND. (unit_stat == 0)) THEN
            IF (status_string == "SCRATCH") THEN
               IF (form_string == "FORMATTED") THEN
                  OPEN (UNIT=unit_number, &
                        STATUS=TRIM(status_string), &
                        ACCESS=TRIM(access_string), &
                        FORM=TRIM(form_string), &
                        POSITION=TRIM(position_string), &
                        ACTION=TRIM(action_string), &
                        PAD=TRIM(pad_string), &
                        IOMSG=iomsg_string, &
                        IOSTAT=open_stat)
               ELSE
                  OPEN (UNIT=unit_number, &
                        STATUS=TRIM(status_string), &
                        ACCESS=TRIM(access_string), &
                        FORM=TRIM(form_string), &
                        POSITION=TRIM(position_string), &
                        ACTION=TRIM(action_string), &
                        IOMSG=iomsg_string, &
                        IOSTAT=open_stat)
               END IF
            ELSE
               IF (form_string == "FORMATTED") THEN
                  OPEN (UNIT=unit_number, &
                        FILE=TRIM(real_file_name), &
                        STATUS=TRIM(status_string), &
                        ACCESS=TRIM(access_string), &
                        FORM=TRIM(form_string), &
                        POSITION=TRIM(position_string), &
                        ACTION=TRIM(action_string), &
                        PAD=TRIM(pad_string), &
                        IOMSG=iomsg_string, &
                        IOSTAT=open_stat)
               ELSE
                  OPEN (UNIT=unit_number, &
                        FILE=TRIM(real_file_name), &
                        STATUS=TRIM(status_string), &
                        ACCESS=TRIM(access_string), &
                        FORM=TRIM(form_string), &
                        POSITION=TRIM(position_string), &
                        ACTION=TRIM(action_string), &
                        IOMSG=iomsg_string, &
                        IOSTAT=open_stat)
               END IF
            END IF
         END IF
      END IF
!$OMP END CRITICAL(cp_files_connect)

      IF (inquire_stat /= 0) THEN
         CALL cp_abort(__LOCATION__, &
                       "An error occurred inquiring the file <"//TRIM(real_file_name)// &
                       "> (IOSTAT="//TRIM(int_to_string(inquire_stat))//").")
      END IF
      IF (unit_stat /= 0) THEN
         CALL cp_abort(__LOCATION__, &
                       "An error occurred inquiring the requested unit "// &
                       TRIM(int_to_string(unit_number))//" (IOSTAT="// &
                       TRIM(int_to_string(unit_stat))//").")
      END IF

      IF (is_open) THEN
         ! The file is already connected: check the connection and reuse it
         IF (skip_get .AND. (unit_number /= requested_unit_number)) THEN
            CALL cp_abort(__LOCATION__, &
                          "The file <"//TRIM(real_file_name)//"> is already connected to "// &
                          "unit "//TRIM(int_to_string(unit_number))//" and cannot be opened "// &
                          "on the requested unit "// &
                          TRIM(int_to_string(requested_unit_number))//".")
         END IF
         CALL uppercase(current_action)
         IF (TRIM(current_action) /= TRIM(action_string)) THEN
            CALL cp_abort(__LOCATION__, &
                          "Attempt to re-open the existing file <"// &
                          TRIM(real_file_name)//"> with the modified ACTION attribute "// &
                          TRIM(action_string)//". The current ACTION attribute is "// &
                          TRIM(current_action)//".")
         END IF
         CALL uppercase(current_form)
         IF (TRIM(current_form) /= TRIM(form_string)) THEN
            CALL cp_abort(__LOCATION__, &
                          "Attempt to re-open the existing "// &
                          TRIM(current_form)//" file <"//TRIM(real_file_name)// &
                          "> as "//TRIM(form_string)//" file.")
         END IF
         CALL uppercase(current_access)
         IF (TRIM(current_access) /= TRIM(access_string)) THEN
            CALL cp_abort(__LOCATION__, &
                          "Attempt to re-open the existing "// &
                          TRIM(current_access)//" file <"//TRIM(real_file_name)// &
                          "> with ACCESS="//TRIM(access_string)//".")
         END IF
         IF (status_string == "NEW") THEN
            CALL cp_abort(__LOCATION__, &
                          "Attempt to re-open the existing file <"// &
                          TRIM(real_file_name)//"> with STATUS=NEW.")
         END IF
         IF (position_string == "REWIND") REWIND (UNIT=unit_number)
      ELSE
         ! Check the outcome of the connection attempt
         SELECT CASE (status_string)
         CASE ("NEW")
            IF (exists) THEN
               CALL cp_abort(__LOCATION__, &
                             "The file <"//TRIM(real_file_name)// &
                             "> cannot be opened with STATUS=NEW, because it already exists.")
            END IF
         CASE ("OLD")
            IF (.NOT. exists) THEN
               IF (action_string == "READ") THEN
                  CALL cp_abort(__LOCATION__, &
                                "The specified file <"//TRIM(real_file_name)// &
                                "> cannot be opened with STATUS=OLD. It does not exist. "// &
                                "Data directory path: "//TRIM(get_data_dir()))
               ELSE
                  CALL cp_abort(__LOCATION__, &
                                "The specified file <"//TRIM(real_file_name)// &
                                "> cannot be opened with STATUS=OLD and ACTION="// &
                                TRIM(action_string)//". It does not exist.")
               END IF
            END IF
         END SELECT
         IF (unit_number < 1) THEN
            IF (skip_get) THEN
               CALL cp_abort(__LOCATION__, &
                             "Invalid logical unit number "//TRIM(int_to_string(unit_number))// &
                             " specified with skip_get_unit_number.")
            ELSE IF (status_string == "SCRATCH") THEN
               CALL cp_abort(__LOCATION__, &
                             "Could not obtain a unit number to open a SCRATCH file.")
            ELSE
               CALL cp_abort(__LOCATION__, &
                             "Cannot open the file <"//TRIM(real_file_name)// &
                             ">, because no unused logical unit number could be obtained.")
            END IF
         END IF
         IF (unit_in_use) THEN
            CALL cp_abort(__LOCATION__, &
                          "The requested logical unit number "//TRIM(int_to_string(unit_number))// &
                          " is already connected to another file.")
         END IF
         IF (open_stat /= 0) THEN
            IF (status_string == "SCRATCH") THEN
               CALL cp_abort(__LOCATION__, &
                             "An error occurred opening a SCRATCH file (UNIT="// &
                             TRIM(int_to_string(unit_number))//", IOSTAT="// &
                             TRIM(int_to_string(open_stat))//"): "//TRIM(iomsg_string))
            ELSE
               CALL m_getcwd(cwd)
               CALL cp_abort(__LOCATION__, &
                             "An error occurred opening the file <"//TRIM(real_file_name)// &
                             "> (UNIT="//TRIM(int_to_string(unit_number))// &
                             ", IOSTAT="//TRIM(int_to_string(open_stat))//"): "// &
                             TRIM(iomsg_string)//". Current working directory: "//TRIM(cwd))
            END IF
         END IF
      END IF

      IF (PRESENT(debug)) THEN
         IF (debug) THEN
            INQUIRE (UNIT=unit_number, &
                     OPENED=is_open, &
                     NAMED=is_named, &
                     POSITION=position_string, &
                     ACCESS=access_string, &
                     FORM=form_string, &
                     ACTION=action_string)
            ! The NAME specifier is undefined for unnamed (SCRATCH) files
            IF (is_named) THEN
               INQUIRE (UNIT=unit_number, NAME=message)
            ELSE
               message = "<unnamed>"
            END IF
            WRITE (UNIT=default_output_unit, FMT="(T2,A)") "BEGIN DEBUG "//TRIM(routineN)
            WRITE (UNIT=default_output_unit, FMT="(T3,A,I0)") "NUMBER  : ", unit_number
            WRITE (UNIT=default_output_unit, FMT="(T3,A,L1)") "OPENED  : ", is_open
            WRITE (UNIT=default_output_unit, FMT="(T3,A)") &
               "NAME    : "//TRIM(message), &
               "POSITION: "//TRIM(position_string), &
               "ACCESS  : "//TRIM(access_string), &
               "FORM    : "//TRIM(form_string), &
               "ACTION  : "//TRIM(action_string), &
               "END DEBUG "//TRIM(routineN)
         END IF
      END IF

   END SUBROUTINE open_file

! **************************************************************************************************
!> \brief Checks if a file exists, considering also the file discovery mechanism.
!> \param file_name Name of the file
!> \return True, if the file exists in the current or the data directory
!> \author Ole Schuett
! **************************************************************************************************
   FUNCTION file_exists(file_name) RESULT(exist)
      CHARACTER(LEN=*), INTENT(IN)                       :: file_name
      LOGICAL                                            :: exist

      CHARACTER(LEN=default_path_length)                 :: real_file_name
      INTEGER                                            :: istat

      real_file_name = discover_file(file_name)
      INQUIRE (FILE=TRIM(real_file_name), EXIST=exist, IOSTAT=istat)
      IF (istat /= 0) exist = .FALSE.

   END FUNCTION file_exists

! **************************************************************************************************
!> \brief Searches for a file in the current directory and then in the data directory.
!> \param file_name Name of the file
!> \return File name as given, if it exists, otherwise its path in the data directory, if
!>         it exists there, else the file name stripped of blanks
!> \author Ole Schuett
! **************************************************************************************************
   FUNCTION discover_file(file_name) RESULT(real_file_name)
      CHARACTER(LEN=*), INTENT(IN)                       :: file_name
      CHARACTER(LEN=default_path_length)                 :: real_file_name

      CHARACTER(LEN=default_path_length)                 :: candidate, data_dir
      INTEGER                                            :: stat
      LOGICAL                                            :: exists

      ! Strip leading and trailing blanks from file name
      real_file_name = TRIM(ADJUSTL(file_name))

      IF (LEN_TRIM(real_file_name) == 0) THEN
         CALL cp_abort(__LOCATION__, "A file name length of zero for an existing file is invalid.")
      END IF

      ! First try file name directly
      INQUIRE (FILE=TRIM(real_file_name), EXIST=exists, IOSTAT=stat)
      IF (stat == 0 .AND. exists) RETURN

      ! Then try the data directory
      data_dir = get_data_dir()
      IF (LEN_TRIM(data_dir) > 0) THEN
         candidate = join_paths(data_dir, real_file_name)
         INQUIRE (FILE=TRIM(candidate), EXIST=exists, IOSTAT=stat)
         IF (stat == 0 .AND. exists) THEN
            real_file_name = candidate
            RETURN
         END IF
      END IF

   END FUNCTION discover_file

! **************************************************************************************************
!> \brief Returns path of data directory if set, otherwise an empty string
!> \return Path of the data directory or an empty string
!> \author Ole Schuett
! **************************************************************************************************
   FUNCTION get_data_dir() RESULT(res)
      CHARACTER(len=default_path_length)                 :: res

      CHARACTER(LEN=1, KIND=C_CHAR), DIMENSION(:), &
         POINTER                                         :: path_f
      INTEGER                                            :: i
      TYPE(C_PTR)                                        :: path_c
      INTERFACE
         FUNCTION get_data_dir_c() BIND(C, name="get_data_dir")
            IMPORT :: C_PTR
            TYPE(C_PTR)                               :: get_data_dir_c
         END FUNCTION get_data_dir_c
      END INTERFACE

      path_c = get_data_dir_c()
      CALL C_F_POINTER(path_c, path_f, shape=[default_path_length])

      res = ""
      DO i = 1, default_path_length
         IF (path_f(i) == C_NULL_CHAR) RETURN
         res(i:i) = path_f(i)
      END DO

      CALL cp_abort(__LOCATION__, "The CP2K_DATA_DIR path is too long")

   END FUNCTION get_data_dir

! **************************************************************************************************
!> \brief Joins two file paths, inserting '/' as needed.
!> \param path1 Leading path, e.g. a directory
!> \param path2 Trailing path; returned unchanged if it is an absolute path
!> \return Joined path; the program is aborted, if it exceeds default_path_length
!> \author Ole Schuett
! **************************************************************************************************
   FUNCTION join_paths(path1, path2) RESULT(joined_path)
      CHARACTER(LEN=*), INTENT(IN)                       :: path1, path2
      CHARACTER(LEN=default_path_length)                 :: joined_path

      INTEGER                                            :: joined_length, n1, n2
      LOGICAL                                            :: absolute

      n1 = LEN_TRIM(path1)
      n2 = LEN_TRIM(path2)

      absolute = .FALSE.
      IF (n2 > 0) absolute = (path2(1:1) == '/')

      ! Determine the length of the joined path before assigning it to avoid silent truncation
      IF (absolute .OR. (n1 == 0)) THEN
         joined_length = n2
      ELSE IF (path1(n1:n1) == '/') THEN
         joined_length = n1 + n2
      ELSE
         joined_length = n1 + 1 + n2
      END IF
      IF (joined_length > LEN(joined_path)) THEN
         CALL cp_abort(__LOCATION__, &
                       "The joined path of <"//path1(1:n1)//"> and <"//path2(1:n2)// &
                       "> exceeds the maximum length of "// &
                       TRIM(int_to_string(LEN(joined_path)))//" characters.")
      END IF

      IF (absolute .OR. (n1 == 0)) THEN
         joined_path = path2(1:n2)
      ELSE IF (path1(n1:n1) == '/') THEN
         joined_path = path1(1:n1)//path2(1:n2)
      ELSE
         joined_path = path1(1:n1)//'/'//path2(1:n2)
      END IF

   END FUNCTION join_paths

! **************************************************************************************************
!> \brief Converts an integer to a left-adjusted string.
!> \param i Integer value
!> \return String representation of i
!> \author Matthias Krack (MK)
! **************************************************************************************************
   FUNCTION int_to_string(i) RESULT(string)
      INTEGER, INTENT(IN)                                :: i
      CHARACTER(LEN=11)                                  :: string

      WRITE (UNIT=string, FMT="(I0)") i

   END FUNCTION int_to_string

END MODULE cp_files
