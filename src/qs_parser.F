!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_parser [1.0] *
!!
!!   NAME
!!     qs_parser
!!
!!   FUNCTION
!!     Utility routines to read data from files.
!!
!!   AUTHOR
!!     Matthias Krack (22.11.1999)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE qs_parser

  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_sync
  USE string_utilities,                ONLY: compress,&
                                             uppercase
  USE termination,                     ONLY: get_error_unit,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "qs_parser"

  INTEGER, PARAMETER :: max_file_name_length = 200,&
                        max_line_length = 240

! *** Global variables ***

  CHARACTER(LEN=4),PARAMETER :: default_separators = ",:;="
  CHARACTER(LEN=3),PARAMETER :: default_end_section_label = "END"
  CHARACTER(LEN=1),PARAMETER :: blank_character = " ",&
                      default_comment_character = "#",&
                      default_continuation_character = CHAR(92),& ! backslash
                      default_section_character = "&",&
                      horizontal_tab = CHAR(9)
  INTEGER          :: input_unit = 0,&
                      max_message_length = 200,&
                      max_unit_number = 999

  CHARACTER(LEN=40) :: end_section,&
                       start_section
  CHARACTER(LEN=10) :: separators
  CHARACTER(LEN=1)  :: comment_character,&
                       continuation_character,&
                       section_character

  INTEGER, DIMENSION(2), PARAMETER :: reserved_unit_numbers = (/5,6/)

  CHARACTER(LEN=max_file_name_length) :: input_file_name
  CHARACTER(LEN=max_line_length)      :: input_line

! icol             : Number of the current column in the current input line
! icol1            : First column of the current input string
! icol2            : Last column of the current input string
! input_line_number: Number of the current input line read from the input file
! input_unit       : Logical unit number of the input file

  INTEGER :: group,icol,icol1,icol2,input_line_number,output_unit,source
  LOGICAL :: first_separator,ionode

! *** Public subroutines ***

  PUBLIC :: close_file,&
            finish_parser,&
            integer_object,&
            open_file,&
            read_object,&
            search_object,&
            start_parser,&
            stop_parser,&
            test_object

! *** Public functions ***

  PUBLIC :: get_unit_number

! *****************************************************************************

  INTERFACE read_object
    MODULE PROCEDURE read_integer_object,&
                     read_logical_object,&
                     read_real_object,&
                     read_string_object
  END INTERFACE

  INTERFACE stop_parser
    MODULE PROCEDURE stop_parser_new,stop_parser_old
  END INTERFACE

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE broadcast_input_information(ifail)

!   Purpose: Broadcast the input information.

!   History: - Creation (02.03.2001,MK)

!   ***************************************************************************

    INTEGER, INTENT(IN), OPTIONAL            :: ifail

!   ---------------------------------------------------------------------------

    IF (PRESENT(ifail)) CALL mp_bcast(ifail,source,group)
    CALL mp_bcast(input_line,source,group)
    CALL mp_bcast(input_line_number,source,group)
    CALL mp_bcast(icol,source,group)

  END SUBROUTINE broadcast_input_information

! *****************************************************************************

  SUBROUTINE close_file(unit_number,file_status)

!   Purpose: Close the current input file.

!   History: - Creation (14.02.2001,MK)

!   ***************************************************************************

    INTEGER, INTENT(IN)                      :: unit_number
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: file_status

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "close_file"

    CHARACTER(LEN=6)                         :: status_string
    CHARACTER(LEN=max_message_length)        :: message
    INTEGER                                  :: istat
    LOGICAL                                  :: exists, opened

!   ---------------------------------------------------------------------------
!   *** Check the specified input file name ***

    INQUIRE (UNIT=unit_number,EXIST=exists,OPENED=opened,IOSTAT=istat)

    IF (istat /= 0) THEN
      WRITE (UNIT=message,FMT="(A,I6,A,I6,A)")&
        "An error occurred inquiring the unit with the number ",unit_number,&
        " (IOSTAT = ",istat,")"
      CALL compress(message)
      CALL stop_program(routine_name,module_name,__LINE__,message)
    ELSE IF (.NOT.exists) THEN
      WRITE (UNIT=message,FMT="(A,I6,A,I6,A)")&
        "The specified unit number ",unit_number," does not exist"
      CALL compress(message)
      CALL stop_program(routine_name,module_name,__LINE__,message)
    END IF

!   *** Close the specified file ***

    IF (opened) THEN

      IF (PRESENT(file_status)) THEN
        status_string = file_status
        CALL uppercase(status_string)
      ELSE
        status_string = "KEEP"
      END IF

      CLOSE (UNIT=unit_number,IOSTAT=istat,STATUS=TRIM(status_string))

      IF (istat /= 0) THEN
        WRITE (UNIT=message,FMT="(A,I6,A,I6,A)")&
          "An error occurred closing the file <"//TRIM(input_file_name)//&
          "> with the unit number ",unit_number," (IOSTAT = ",istat,")"
        CALL compress(message)
        CALL stop_program(routine_name,module_name,__LINE__,message)
      END IF

    END IF

  END SUBROUTINE close_file

! *****************************************************************************

  SUBROUTINE finish_parser()

!   Purpose: Finish a parser run.

!   History: - Creation (14.02.2001,MK)

!   ***************************************************************************


    CHARACTER(LEN=*), PARAMETER              :: routine_name = "finish_parser"

!   ---------------------------------------------------------------------------

    IF (ionode) THEN
      IF (input_unit == 0) THEN
        CALL stop_program(routine_name,module_name,__LINE__,&
                          "Cannot finish an input session that was never "//&
                          " started")
      ELSE
        CALL close_file(unit_number=input_unit)
      END IF
    END IF

    CALL mp_sync(group)

    input_unit = 0

  END SUBROUTINE finish_parser

! *****************************************************************************

  SUBROUTINE get_next_input_line(nline,ifail)

!   Purpose: Read the next input line and broadcast the input information.
!            Skip (nline-1) lines and skip also all comment lines.

!   History: - Creation (22.11.1999,MK)

!   ***************************************************************************

    INTEGER, INTENT(IN)                      :: nline
    INTEGER, INTENT(OUT), OPTIONAL           :: ifail

!   ---------------------------------------------------------------------------

    IF (PRESENT(ifail)) THEN
      CALL read_next_input_line(nline,ifail)
      CALL broadcast_input_information(ifail)
    ELSE
      CALL read_next_input_line(nline)
      CALL broadcast_input_information()
    END IF

  END SUBROUTINE get_next_input_line

! *****************************************************************************

  SUBROUTINE get_next_input_string(string_length)

!   Purpose: Get the next input string from the input line.

!   History: - Creation (19.02.2001,MK)

!   ***************************************************************************

    INTEGER, INTENT(IN), OPTIONAL            :: string_length

    CHARACTER(LEN=*), PARAMETER :: routine_name = "get_next_input_string"

    INTEGER                                  :: i, length

!   ---------------------------------------------------------------------------

    IF (PRESENT(string_length)) THEN
      length = string_length
    ELSE
      length = 0
    END IF

    IF (length > 0) THEN

!     *** Fixed input string length ***

      icol1 = icol + 1
      icol2 = icol + length
      icol = icol2

      IF (icol > LEN(input_line)) THEN
        CALL stop_parser(routine_name,module_name,__LINE__,&
                         "Maximum input line length exceeded when reading "//&
                         "an input string using a fixed string length")
      END IF

!     *** Ignore blanks ***

      DO i=icol1,icol2
        IF (input_line(i:i) /= blank_character) EXIT
      END DO
      icol1 = MIN(i,icol2)

      DO i=icol2,icol1,-1
        IF (input_line(i:i) /= blank_character) EXIT
      END DO
      icol2 = MAX(i,icol1)

    ELSE

!     *** Variable input string length (automatic search) ***

!     *** Search for the beginning of the next input string ***

      DO

!       *** Increment the column counter ***

        icol = icol + 1

!       *** Quick return, if the end of line is found ***

        IF ((icol > LEN_TRIM(input_line)).OR.&
            (input_line(icol:icol) == comment_character)) THEN
          icol1 = 0
          icol2 = -1
          RETURN
        END IF

!       *** Check for input line continuation ***

        IF (input_line(icol:icol) == continuation_character) THEN
          DO i=LEN(input_line),1,-1
            IF ((input_line(i:i) /= blank_character).AND.&
                (input_line(i:i) /= horizontal_tab)) EXIT
          END DO
          IF (icol == i) THEN
            CALL get_next_input_line(1)
            CYCLE
          ELSE
            icol1 = i
            icol2 = i
            CALL stop_parser(routine_name,module_name,__LINE__,&
                             "Found non-blank tokens after the "//&
                             "line continuation character "//&
                             continuation_character)
          END IF
        END IF

!       *** Ignore all white space and accept only one     ***
!       *** separator token or a string in quotation marks ***

        IF ((input_line(icol:icol) == blank_character).OR.&
            (input_line(icol:icol) == horizontal_tab)) THEN
          CYCLE
        ELSE IF (INDEX(separators,input_line(icol:icol)) > 0) THEN
          IF (first_separator) THEN
            first_separator = .FALSE.
            CYCLE
          ELSE
            icol1 = icol
            icol2 = icol
            CALL stop_parser(routine_name,module_name,__LINE__,&
                             "Unexpected separator token <"//&
                             input_line(icol:icol)//"> found")
          END IF
        ELSE IF (input_line(icol:icol) == """") THEN
          first_separator = .TRUE.
          icol1 = icol + 1
          icol2 = icol + INDEX(input_line(icol1:),"""")
          IF (icol2 == icol) THEN
            icol1 = icol
            icol2 = icol
            CALL stop_parser(routine_name,module_name,__LINE__,&
                             "Unmatched quotation mark found")
          ELSE
            icol = icol2
            icol2 = icol2 - 1
            RETURN
          END IF
        ELSE
          first_separator = .TRUE.
          icol1 = icol
          EXIT
        END IF

      END DO

!     *** Search for the end of the next input string ***

      DO
        icol = icol + 1
        IF ((icol > LEN_TRIM(input_line)).OR.&
            (input_line(icol:icol) == blank_character).OR.&
            (input_line(icol:icol) == comment_character).OR.&
            (input_line(icol:icol) == continuation_character).OR.&
            (input_line(icol:icol) == horizontal_tab)) THEN
          EXIT
        ELSE IF (INDEX(separators,input_line(icol:icol)) > 0) THEN
          first_separator = .FALSE.
          EXIT
        END IF
      END DO

      icol2 = icol - 1

      IF (input_line(icol:icol) == continuation_character) icol = icol2

    END IF

  END SUBROUTINE get_next_input_string

! *****************************************************************************

  FUNCTION get_unit_number() RESULT(unit_number)

!   Purpose: Returns a unit number of a unit that exists and is not connected.

!   History: - Creation (13.02.2001,MK)

!   ***************************************************************************

    INTEGER                                  :: unit_number

    INTEGER                                  :: istat
    LOGICAL                                  :: exists, opened

!   ---------------------------------------------------------------------------

    DO unit_number=1,max_unit_number
      IF (ANY(unit_number == reserved_unit_numbers)) CYCLE
      INQUIRE (UNIT=unit_number,EXIST=exists,OPENED=opened,IOSTAT=istat)
      IF (exists.AND.(.NOT.opened).AND.(istat == 0)) RETURN
    END DO

    unit_number = -1

   END FUNCTION get_unit_number

! *****************************************************************************

  FUNCTION integer_object(string) RESULT(contains_integer_object)

!   Purpose: Check, if the string object contains an object of type integer.

!   History: - Creation (22.11.1999,MK)

!   ***************************************************************************

    CHARACTER(LEN=*), INTENT(IN)             :: string
    LOGICAL                                  :: contains_integer_object

    INTEGER                                  :: i, length

!   ---------------------------------------------------------------------------

    contains_integer_object = .TRUE.

    length = LEN_TRIM(string)

    IF (length == 0) THEN
      contains_integer_object = .FALSE.
      RETURN
    END IF

    IF ((INDEX("+-",string(1:1)) > 0).AND.(length == 1)) THEN
      contains_integer_object = .FALSE.
      RETURN
    END IF

    IF (INDEX("+-0123456789",string(1:1)) == 0) THEN
      contains_integer_object = .FALSE.
      RETURN
    END IF

    DO i=2,length
      IF (INDEX("0123456789",string(i:i)) == 0) THEN
        contains_integer_object = .FALSE.
        EXIT
      END IF
    END DO

  END FUNCTION integer_object

! *****************************************************************************

  SUBROUTINE open_file(file_name,file_status,file_form,file_action,&
                       file_position,unit_number)

!   Purpose: Open the file with the name "file_name".

!   History: - Creation (13.02.2001,MK)

!   ***************************************************************************

    CHARACTER(LEN=*), INTENT(IN)             :: file_name
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: file_status, file_form, &
                                                file_action, file_position
    INTEGER, INTENT(OUT)                     :: unit_number

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "open_file"

    CHARACTER(LEN=11)                        :: action_string, form_string, &
                                                position_string, status_string
    CHARACTER(LEN=max_message_length)        :: message
    INTEGER                                  :: istat
    LOGICAL                                  :: exists, opened

!   ---------------------------------------------------------------------------

    IF (PRESENT(file_status)) THEN
      status_string = file_status
      CALL uppercase(status_string)
    ELSE
      status_string = "OLD"
    END IF

    IF (PRESENT(file_form)) THEN
      form_string = file_form
      CALL uppercase(form_string)
    ELSE
      form_string = "FORMATTED"
    END IF

    IF (PRESENT(file_action)) THEN
      action_string = file_action
      CALL uppercase(action_string)
    ELSE
      action_string = "READ"
    END IF

    IF (PRESENT(file_position)) THEN
      position_string = file_position
      CALL uppercase(position_string)
    ELSE
      position_string = "REWIND"
    END IF

!   *** Check the specified input file name ***

    INQUIRE (FILE=TRIM(file_name),EXIST=exists,OPENED=opened,IOSTAT=istat)

    IF (istat /= 0) THEN
      WRITE (UNIT=message,FMT="(A,I6,A)")&
        "An error occurred inquiring the file <"//TRIM(file_name)//&
        "> (IOSTAT = ",istat,")"
      CALL stop_program(routine_name,module_name,__LINE__,message)
    ELSE IF (status_string == "OLD") THEN
      IF (.NOT.exists) THEN
        CALL stop_program(routine_name,module_name,__LINE__,&
                          "The specified file <"//TRIM(file_name)//&
                          "> does not exist")
      END IF
    END IF

!   *** Open the specified input file ***

    IF (opened) THEN

      INQUIRE (FILE=TRIM(file_name),NUMBER=unit_number)
      REWIND (UNIT=unit_number)

    ELSE

!     *** Find an unused unit number ***

      unit_number = get_unit_number()

      IF (unit_number < 0) THEN
        CALL stop_program(routine_name,module_name,__LINE__,&
                          "Could not find a logical unit number which "//&
                          "exists and is not connected")
      END IF

      OPEN (UNIT=unit_number,&
            FILE=TRIM(file_name),&
            STATUS=TRIM(status_string),&
            ACCESS="SEQUENTIAL",&
            FORM=TRIM(form_string),&
            POSITION=TRIM(position_string),&
            ACTION=TRIM(action_string),&
            IOSTAT=istat)

      IF (istat /= 0) THEN
        WRITE (UNIT=message,FMT="(A,I6,A,I6,A)")&
          "An error occurred opening the file <"//TRIM(file_name)//&
          "> with the unit number ",unit_number," (IOSTAT = ",istat,")"
        CALL compress(message)
        CALL stop_program(routine_name,module_name,__LINE__,message)
      END IF

    END IF

    input_file_name = TRIM(file_name)

  END SUBROUTINE open_file

! *****************************************************************************

  SUBROUTINE read_integer_object(object,lower_to_upper,newline,skip_lines,&
                                 string_length)

!   Purpose: Read an integer number.

!   History: - Creation (22.11.1999,MK)

!   ***************************************************************************

    INTEGER, INTENT(OUT)                     :: object
    LOGICAL, INTENT(IN), OPTIONAL            :: lower_to_upper, newline
    INTEGER, INTENT(IN), OPTIONAL            :: skip_lines, string_length

    CHARACTER(LEN=*), PARAMETER :: routine_name = "read_integer_object"

    INTEGER                                  :: nline

!   ---------------------------------------------------------------------------

    IF (PRESENT(skip_lines)) THEN
      nline = skip_lines
    ELSE
      nline = 0
    END IF

    IF (PRESENT(newline)) THEN
      IF (newline) nline = nline + 1
    END IF

    CALL get_next_input_line(nline)

    IF (PRESENT(string_length)) THEN
      CALL get_next_input_string(string_length)
    ELSE
      CALL get_next_input_string()
    END IF

    IF (icol1 > icol2) THEN
      icol1 = icol
      icol2 = icol
      CALL stop_parser(routine_name,module_name,__LINE__,&
                       "An integer type object was expected, "//&
                       "found end of line")
    END IF

    IF (integer_object(input_line(icol1:icol2))) THEN
      READ (UNIT=input_line(icol1:icol2),FMT=*) object
    ELSE
      CALL stop_parser(routine_name,module_name,__LINE__,&
                       "An integer type object was expected, found <"//&
                       input_line(icol1:icol2)//">")
    END IF

  END SUBROUTINE read_integer_object

! *****************************************************************************

  SUBROUTINE read_logical_object(object,newline,skip_lines,string_length)

!   Purpose: Read a string representing logical object.

!   History: - Creation (01.04.2003,FM)
!            - New version (08.07.2003,MK)

!   ***************************************************************************

    LOGICAL, INTENT(OUT)                     :: object
    LOGICAL, INTENT(IN), OPTIONAL            :: newline
    INTEGER, INTENT(IN), OPTIONAL            :: skip_lines, string_length

    CHARACTER(LEN=*), PARAMETER :: routine_name = "read_logical_object"

    CHARACTER(LEN=max_line_length)           :: input_string
    INTEGER                                  :: input_string_length, nline

!   ---------------------------------------------------------------------------

    IF (PRESENT(skip_lines)) THEN
      nline = skip_lines
    ELSE
      nline = 0
    END IF

    IF (PRESENT(newline)) THEN
      IF (newline) nline = nline + 1
    END IF

    CALL get_next_input_line(nline)

    IF (PRESENT(string_length)) THEN
      CALL get_next_input_string(string_length)
    ELSE
      CALL get_next_input_string()
    END IF

    input_string_length = icol2 - icol1 + 1

    IF (input_string_length == 0) THEN
      icol1 = icol
      icol2 = icol
      CALL stop_parser(routine_name,module_name,__LINE__,&
                       "A string representing a logical object was expected, "//&
                       "found end of line")
    ELSE
      input_string=""
      input_string(:input_string_length) = input_line(icol1:icol2)
    END IF

    CALL uppercase(input_string)

    SELECT CASE (TRIM(input_string))
    CASE ("0","F",".F.","FALSE",".FALSE.","N","NO","OFF")
      object = .FALSE.
    CASE ("1","T",".T.","TRUE",".TRUE.","Y","YES","ON")
      object = .TRUE.
    CASE DEFAULT
      CALL stop_parser(routine_name,module_name,__LINE__,&
                       "A string representing a logical object was expected, "//&
                       "found <"//TRIM(input_string)//">")
    END SELECT

  END SUBROUTINE read_logical_object

! *****************************************************************************

  SUBROUTINE read_next_input_line(nline,ifail)

!   Purpose: Read the next line from a logical unit "unit" (I/O node only).
!            Skip (nline-1) lines and skip also all comment lines.

!   History: - Creation (22.11.1999,MK)

!   ***************************************************************************

    INTEGER, INTENT(IN)                      :: nline
    INTEGER, INTENT(OUT), OPTIONAL           :: ifail

    CHARACTER(LEN=*), PARAMETER :: routine_name = "read_next_input_line"

    CHARACTER(LEN=max_message_length)        :: message
    INTEGER                                  :: icomment, iline, istat

!   ---------------------------------------------------------------------------

    IF (PRESENT(ifail)) ifail = 0

    IF (ionode) THEN

      iline = 0
      istat = 0

      DO WHILE (iline /= nline)

        iline = iline + 1

!       *** Increment the input line counter ***

        input_line_number = input_line_number + 1

!       *** Read next line from file ***

        READ (UNIT=input_unit,FMT="(A)",IOSTAT=istat) input_line

!       *** Handle read errors ***

        IF (istat /= 0) THEN
          IF (PRESENT(ifail)) THEN
            ifail = istat
            EXIT
          ELSE
            WRITE (UNIT=message,FMT="(A,I6,A,I6,A)")&
              "An error occurred reading the input file <"//&
              TRIM(input_file_name)//"> with the unit number ",input_unit,&
              " (IOSTAT = ",istat,")"
            CALL compress(message)
            CALL stop_program(routine_name,module_name,__LINE__,message)
          END IF
        END IF

!       *** Check for comments ***

        icomment = INDEX(input_line,comment_character)

!       *** Ignore empty lines ***

        IF (icomment > 0) THEN
!MK xlf 8.1.1 bug -> len_trim returns -1 for an empty string
!         IF (LEN_TRIM(input_line(:icomment-1)) == 0) iline = iline - 1
          IF (LEN_TRIM(input_line(:icomment-1)) <= 0) iline = iline - 1
        ELSE
!MK xlf 8.1.1 bug -> len_trim returns -1 for an empty string
!         IF (LEN_TRIM(input_line) == 0) iline = iline - 1
          IF (LEN_TRIM(input_line) <= 0) iline = iline - 1
        END IF

      END DO

!     *** Reset column pointer, if a new line was read ***

      IF (nline > 0) icol = 0

    END IF

  END SUBROUTINE read_next_input_line

! *****************************************************************************

  SUBROUTINE read_real_object(object,lower_to_upper,newline,skip_lines,&
                              string_length)

!   Purpose: Read a floating point number.

!   History: - Creation (22.11.1999,MK)

!   ***************************************************************************

    REAL(KIND=dp), INTENT(OUT)               :: object
    LOGICAL, INTENT(IN), OPTIONAL            :: lower_to_upper, newline
    INTEGER, INTENT(IN), OPTIONAL            :: skip_lines, string_length

    CHARACTER(LEN=*), PARAMETER :: routine_name = "read_real_object"

    CHARACTER(LEN=max_line_length), &
      DIMENSION(2)                           :: string
    INTEGER                                  :: idot_first, idot_last, &
                                                islash, istat, iz, nline, nz
    REAL(KIND=dp), DIMENSION(2)              :: z

!   ---------------------------------------------------------------------------

    IF (PRESENT(skip_lines)) THEN
      nline = skip_lines
    ELSE
      nline = 0
    END IF

    IF (PRESENT(newline)) THEN
      IF (newline) nline = nline + 1
    END IF

    CALL get_next_input_line(nline)

    IF (PRESENT(string_length)) THEN
      CALL get_next_input_string(string_length)
    ELSE
      CALL get_next_input_string()
    END IF

    IF (icol1 > icol2) THEN
      icol1 = icol
      icol2 = icol
      CALL stop_parser(routine_name,module_name,__LINE__,&
                       "A real type object was expected, found end of line")
    END IF

    islash = icol1 + INDEX(input_line(icol1:icol2),"/") - 1

    IF (islash > icol1) THEN
      nz = 2
      string(1) = input_line(icol1:islash-1)
      string(2) = input_line(islash+1:icol2)
    ELSE
      nz = 1
      string(1) = input_line(icol1:icol2)
      z(2) = 1.0_dp
    END IF

    DO iz=1,nz

      IF (LEN_TRIM(string(iz)) == 0) THEN
        CALL stop_parser(routine_name,module_name,__LINE__,&
                         "A real type object was expected, found end of line")
      END IF

      idot_first = INDEX(string(iz),".")

      IF (idot_first > 0) THEN
        idot_last = INDEX(string(iz),".",.TRUE.)
        istat = 0
        IF (idot_first == idot_last) THEN
          READ (UNIT=string(iz),FMT=*,IOSTAT=istat) z(iz)
        ELSE
          istat = 1
        END IF
        IF (istat /= 0) THEN
          CALL stop_parser(routine_name,module_name,__LINE__,&
                           "A real type object was expected, found <"//&
                           TRIM(string(iz))//">")
        END IF
      ELSE IF (integer_object(string(iz))) THEN
        READ (UNIT=string(iz),FMT=*) z(iz)
      ELSE
        CALL stop_parser(routine_name,module_name,__LINE__,&
                         "A real type object was expected, found <"//&
                         TRIM(string(iz))//">")
      END IF

    END DO

    object = z(1)/z(2)

  END SUBROUTINE read_real_object

! *****************************************************************************

  SUBROUTINE read_string_object(object,lower_to_upper,newline,skip_lines,&
                                string_length)

!   Purpose: Read a string.

!   History: - Creation (22.11.1999,MK)

!   ***************************************************************************

    CHARACTER(LEN=*), INTENT(OUT)            :: object
    LOGICAL, INTENT(IN), OPTIONAL            :: lower_to_upper, newline
    INTEGER, INTENT(IN), OPTIONAL            :: skip_lines, string_length

    CHARACTER(LEN=*), PARAMETER :: routine_name = "read_string_object"

    CHARACTER(LEN=max_message_length)        :: message
    INTEGER                                  :: input_string_length, nline

!   ---------------------------------------------------------------------------

    object = ""

    IF (PRESENT(skip_lines)) THEN
      nline = skip_lines
    ELSE
      nline = 0
    END IF

    IF (PRESENT(newline)) THEN
      IF (newline) nline = nline + 1
    END IF

    CALL get_next_input_line(nline)

    IF (PRESENT(string_length)) THEN
      CALL get_next_input_string(string_length)
    ELSE
      CALL get_next_input_string()
    END IF

    input_string_length = icol2 - icol1 + 1

    IF (input_string_length == 0) THEN
      icol1 = icol
      icol2 = icol
      CALL stop_parser(routine_name,module_name,__LINE__,&
                       "A string type object was expected, found end of line")
    ELSE IF (input_string_length > LEN(object)) THEN
      WRITE (UNIT=message,FMT="(A,I6,A)")&
        "The input string <"//input_line(icol1:icol2)//"> has more than ",&
        LEN(object)," characters and is therefore too long to fit in the "//&
        "specified variable"
      CALL compress(message)
      CALL stop_parser(routine_name,module_name,__LINE__,message)
    ELSE
      object(:input_string_length) = input_line(icol1:icol2)
    END IF

!   *** Convert lowercase to uppercase, if requested ***

    IF (PRESENT(lower_to_upper)) THEN
      IF (lower_to_upper) CALL uppercase(object)
    END IF

  END SUBROUTINE read_string_object

! *****************************************************************************

  SUBROUTINE search_object(string,ignore_case,found,line)

!   Purpose: Search a string pattern in a file defined by its logical unit
!            number "unit". A case sensitive search is performed, if
!            ignore_case is .FALSE..

!   History: - Creation (05.10.1999,MK)

!   ***************************************************************************

    CHARACTER(LEN=*), INTENT(IN)             :: string
    LOGICAL, INTENT(IN)                      :: ignore_case
    LOGICAL, INTENT(OUT)                     :: found
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: line

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "search_object"

    CHARACTER(LEN=LEN(string))               :: pattern
    CHARACTER(LEN=max_line_length+1)         :: current_line
    CHARACTER(LEN=max_message_length)        :: message
    INTEGER                                  :: ifail, ipattern

!   ---------------------------------------------------------------------------

    found = .FALSE.

    IF (PRESENT(line)) line = ""

!   *** Search for string pattern ***

    IF (ionode) THEN

      pattern = string

      IF (ignore_case) CALL uppercase(pattern)

      DO

        CALL read_next_input_line(1,ifail)

!       *** Exit loop, if the end of file is reached ***

        IF (ifail < 0) EXIT

!       *** Check the current line for string pattern ***

        current_line = input_line

        IF (ignore_case) CALL uppercase(current_line)

        ipattern = INDEX(current_line,pattern)

        IF (ipattern > 0) THEN
          found = .TRUE.
          icol = ipattern - 1
          IF (PRESENT(line)) THEN
            IF (LEN(line) < LEN_TRIM(input_line)) THEN
              WRITE (UNIT=message,FMT="(A,I6,A)")&
                "The returned input line has more than ",LEN(line),&
                " characters and is therefore too long to fit in the "//&
                "specified variable"
              CALL compress(message)
              CALL stop_program(routine_name,module_name,__LINE__,message)
            END IF
          END IF
          EXIT
        END IF

      END DO

    END IF

    CALL mp_bcast(found,source,group)

    CALL broadcast_input_information(ifail)

    IF (found) THEN
      IF (PRESENT(line)) line = input_line
      CALL get_next_input_string()
    END IF

  END SUBROUTINE search_object

! *****************************************************************************

  SUBROUTINE start_parser(file_name,globenv,end_section_label,&
                          separator_chars,start_section_label,&
                          start_section_name,comment_char,continuation_char,&
                          section_char,section_found)

!   Purpose: Start a parser run.

!   History: - Creation (14.02.2001,MK)

!   ***************************************************************************

    CHARACTER(LEN=*), INTENT(IN)             :: file_name
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: end_section_label, &
                                                separator_chars, &
                                                start_section_label, &
                                                start_section_name
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL   :: comment_char, &
                                                continuation_char, &
                                                section_char
    LOGICAL, INTENT(OUT), OPTIONAL           :: section_found

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "start_parser"

    CHARACTER(LEN=3)                         :: test_result
    CHARACTER(LEN=60)                        :: current_name, requested_name
    LOGICAL                                  :: found

!   ---------------------------------------------------------------------------

    ionode = globenv%ionode
    group = globenv%group
    source = globenv%source

    separators = ""
    comment_character = ""
    continuation_character = ""
    section_character = ""
    end_section = ""

!   *** Load the default values and overwrite them, if requested ***

    separators = default_separators
    IF (PRESENT(separator_chars)) separators = separator_chars
    comment_character = default_comment_character
    IF (PRESENT(comment_char)) comment_character = comment_char
    continuation_character = default_continuation_character
    IF (PRESENT(continuation_char)) continuation_character = continuation_char
    section_character = default_section_character
    IF (PRESENT(section_char)) section_character = section_char
    end_section = section_character//default_end_section_label
    IF (PRESENT(end_section_label)) THEN
      end_section = section_character//TRIM(end_section_label)
    END IF

    CALL compress(end_section)

!   *** Get the logical output unit number for error messages ***

    CALL get_error_unit(output_unit)

    IF (ionode) THEN
      IF (input_unit == 0) THEN
        CALL open_file(file_name=TRIM(file_name),unit_number=input_unit)
      ELSE
        WRITE (UNIT=output_unit,FMT="(/,T2,A,I6,/,T2,A)")&
          "Current input unit: ",input_unit,&
          "Current input file: """//TRIM(input_file_name)//""""
        CALL stop_program(routine_name,module_name,__LINE__,&
                          "Cannot start a new input session, because an "//&
                          "old input session was not finished")
      END IF
    END IF

    input_line_number = 0

    found = .FALSE.

    IF (PRESENT(start_section_label)) THEN
      start_section = section_character//TRIM(start_section_label)
      IF (PRESENT(start_section_name)) THEN
        requested_name = start_section_name
        CALL uppercase(requested_name)
        DO
          CALL search_object(TRIM(start_section),.TRUE.,found)
          IF (found) THEN
            test_result = test_object()
            IF (test_result == "EOL") THEN
              CYCLE
            ELSE
              CALL read_object(current_name,lower_to_upper=.TRUE.)
              IF (current_name == requested_name) EXIT
            END IF
          ELSE
            EXIT
          END IF
        END DO
        start_section = section_character//TRIM(start_section_label)//" "//&
                        TRIM(requested_name)
      ELSE
        CALL search_object(TRIM(start_section),.TRUE.,found)
      END IF
    ELSE
      start_section = section_character//"not specified"
    END IF

    CALL compress(start_section)

    IF (PRESENT(section_found)) section_found = found

  END SUBROUTINE start_parser

! *****************************************************************************

  SUBROUTINE stop_parser_old(routine,message)

!   Purpose: Print an input error message and stop the program run.

!            There are default error messages for common instances:

!            - END_OF_FILE or EOF
!            - INVALID_KEYWORD
!            - SECTION_NOT_FOUND

!   History: - Creation (15.02.2001,MK)

!   ***************************************************************************

    CHARACTER(LEN=*), INTENT(IN)             :: routine, message

    CHARACTER(LEN=max_message_length)        :: error_message

!   ---------------------------------------------------------------------------

    IF (ionode) THEN

      INQUIRE (UNIT=input_unit,NAME=input_file_name)

      SELECT CASE (message)
      CASE ("END_OF_FILE","EOF")
        error_message = "Unexpected end of file found while reading the "//&
                        "input section <"//TRIM(start_section)//&
                        "> from the input file <"//TRIM(input_file_name)//">"
      CASE ("SECTION_NOT_FOUND")
        error_message = "No input section <"//TRIM(start_section)//&
                        "> found while reading from the input file <"//&
                        TRIM(input_file_name)//">"
      CASE DEFAULT
        WRITE (UNIT=output_unit,FMT="(/,/,T2,A)")&
          "Current input file: "//TRIM(input_file_name)
        WRITE (UNIT=error_message,FMT="(A,I6,A,I6,A)")&
          "ERROR in line ",input_line_number,", column ",icol1,":"
        CALL compress(error_message)
        WRITE (UNIT=output_unit,FMT="(/,T2,A,/,/,T2,A,/,T2,A)")&
          TRIM(error_message),TRIM(input_line),&
          REPEAT(blank_character,MAX(0,icol1-1))//&
          REPEAT("^",MAX(1,icol2-icol1+1))
        SELECT CASE (message)
        CASE ("INVALID_KEYWORD")
          error_message = "Invalid keyword found in the input section <"//&
                          TRIM(start_section)//">"
        CASE DEFAULT
          error_message = message
        END SELECT
      END SELECT

      CALL stop_program(routine,error_message)

    END IF

    CALL mp_sync(group)

  END SUBROUTINE stop_parser_old

! *****************************************************************************

  SUBROUTINE stop_parser_new(routine_name,module_name,line_number,message)

!   Purpose: Print an input error message and stop the program run.

!            There are default error messages for common instances:

!            - END_OF_FILE or EOF
!            - INVALID_KEYWORD
!            - SECTION_NOT_FOUND

!   History: - Creation (15.02.2001,MK)
!            - Adapted to new stop_program version (24.07.2003,MK)

!   ***************************************************************************

    CHARACTER(LEN=*), INTENT(IN)             :: routine_name, module_name
    INTEGER, INTENT(IN)                      :: line_number
    CHARACTER(LEN=*), INTENT(IN)             :: message

    CHARACTER(LEN=max_line_length)           :: indicators
    CHARACTER(LEN=max_message_length)        :: error_message
    INTEGER                                  :: i

!   ---------------------------------------------------------------------------

    IF (ionode) THEN

      INQUIRE (UNIT=input_unit,NAME=input_file_name)

      SELECT CASE (message)
      CASE ("END_OF_FILE","EOF")
        error_message = "Unexpected end of file found while reading the "//&
                        "input section <"//TRIM(start_section)//&
                        "> from the input file <"//TRIM(input_file_name)//">"
      CASE ("SECTION_NOT_FOUND")
        error_message = "No input section <"//TRIM(start_section)//&
                        "> found while reading from the input file <"//&
                        TRIM(input_file_name)//">"
      CASE DEFAULT
        WRITE (UNIT=output_unit,FMT="(/,/,T2,A)")&
          "Current input file: "//TRIM(input_file_name)
        WRITE (UNIT=error_message,FMT="(A,I6,A,I6,A)")&
          "ERROR in line ",input_line_number,", column ",icol1,":"
        CALL compress(error_message)
        indicators = blank_character
        DO i=1,icol1-1
          IF (input_line(i:i) == horizontal_tab) THEN
            indicators(i:i) = horizontal_tab
          END IF
        END DO
        WRITE (UNIT=indicators(icol1:),FMT="(A)")&
          REPEAT("^",MAX(1,icol2-icol1+1))
        WRITE (UNIT=output_unit,FMT="(/,T2,A,/,/,T2,A,/,T2,A)")&
          TRIM(error_message),TRIM(input_line),TRIM(indicators)
        SELECT CASE (message)
        CASE ("INVALID_KEYWORD")
          error_message = "Invalid keyword found in the input section <"//&
                          TRIM(start_section)//">"
        CASE DEFAULT
          error_message = message
        END SELECT
      END SELECT

      CALL stop_program(routine_name,module_name,line_number,error_message)

    END IF

    CALL mp_sync(group)

  END SUBROUTINE stop_parser_new

! *****************************************************************************

  FUNCTION test_object(newline,skip_lines,string_length) RESULT(test_result)

!   Purpose: Test next input object.

!            test_result = "EOF": End of file
!            test_result = "EOL": End of line
!            test_result = "EOS": End of section
!            test_result = "FLT": Floating point number
!            test_result = "INT": Integer number
!            test_result = "STR": String

!   History: - Creation (23.11.1999,MK)

!   ***************************************************************************

    LOGICAL, INTENT(IN), OPTIONAL            :: newline
    INTEGER, INTENT(IN), OPTIONAL            :: skip_lines, string_length
    CHARACTER(LEN=3)                         :: test_result

    CHARACTER(LEN=max_line_length)           :: old_input_line
    CHARACTER(LEN=max_line_length), &
      DIMENSION(2)                           :: string
    INTEGER :: idot_first, idot_last, ifail, iline, islash, istat, iz, nline, &
      nz, old_icol, old_icol1, old_icol2, old_input_line_number
    REAL(KIND=dp)                            :: z

!   ---------------------------------------------------------------------------

    test_result = ""

!   *** Store current status ***

    old_input_line = input_line
    old_input_line_number = input_line_number
    old_icol = icol
    old_icol1 = icol1
    old_icol2 = icol2

    IF (PRESENT(skip_lines)) THEN
      nline = skip_lines
    ELSE
      nline = 0
    END IF

    IF (PRESENT(newline)) THEN
      IF (newline) nline = nline + 1
    END IF

    CALL get_next_input_line(nline,ifail)

    IF (ifail < 0) THEN
      test_result = "EOF"
      RETURN
    END IF

    IF (PRESENT(string_length)) THEN
      CALL get_next_input_string(string_length)
    ELSE
      CALL get_next_input_string()
    END IF

    IF (icol1 > icol2) THEN
      test_result = "EOL"
      RETURN
    END IF

    islash = icol1 + INDEX(input_line(icol1:icol2),"/") - 1

    IF (islash > icol1) THEN
      nz = 2
      string(1) = input_line(icol1:islash-1)
      string(2) = input_line(islash+1:icol2)
    ELSE
      nz = 1
      string(1) = input_line(icol1:icol2)
    END IF

    DO iz=1,nz

      IF (LEN_TRIM(string(iz)) == 0) THEN
        test_result = "EOL"
        EXIT
      END IF

      idot_first = INDEX(string(iz),".")

      IF (idot_first > 0) THEN
        idot_last = INDEX(string(iz),".",.TRUE.)
        istat = 0
        IF (idot_first == idot_last) THEN
          READ (UNIT=string(iz),FMT=*,IOSTAT=istat) z
        ELSE
          istat = 1
        END IF
        IF (istat /= 0) THEN
          test_result = "STR"
        ELSE
          test_result = "FLT"
        END IF
      ELSE IF (integer_object(string(iz))) THEN
        IF (nz == 1) THEN
          test_result = "INT"
        ELSE
          test_result = "FLT"
        END IF
      ELSE
        IF (string(iz) == end_section) THEN
          test_result = "EOS"
        ELSE
          test_result = "STR"
        END IF
      END IF

    END DO

!   *** Reset to old status ***

    IF (ionode) THEN
      IF (nline > 0) THEN
        DO iline=1,input_line_number-old_input_line_number
          BACKSPACE (input_unit)
        END DO
      END IF
    END IF

    input_line = old_input_line
    input_line_number = old_input_line_number
    icol = old_icol
    icol1 = old_icol1
    icol2 = old_icol2

  END FUNCTION test_object

! *****************************************************************************

END MODULE qs_parser
