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
!!     MK (22.11.1999)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE qs_parser

! *****************************************************************************

! icol             : Number of the current column in the current input line
! icol1            : First column of the current input string
! icol2            : Last column of the current input string
! input_line_number: Number of the current input line read from the input file
! input_unit       : Logical unit number of the input file

! *****************************************************************************

  USE kinds, ONLY: wp => dp

  USE global_types,     ONLY: global_environment_type
  USE message_passing,  ONLY: mp_bcast,&
                              mp_sync
  USE output_utilities, ONLY: print_warning
  USE string_utilities, ONLY: compress,&
                              uppercase
  USE termination,      ONLY: get_error_unit,&
                              stop_program

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "qs_parser"

  INTEGER, PARAMETER :: max_file_name_length = 200,&
                        max_line_length = 240

  CHARACTER(LEN=4) :: default_end_section_label = "END",&
                      default_separators = ",:;="
  CHARACTER(LEN=1) :: blank_character = " ",&
                      default_comment_character = "#",&
                      default_continuation_character = CHAR(92),& ! backslash
                      default_section_character = "&"
  INTEGER          :: input_unit = 0,&
                      max_keyword_length = 40,&
                      max_message_length = 200,&
                      max_unit_number = 999

  CHARACTER(LEN=40) :: end_section,start_section
  CHARACTER(LEN=10) :: separators
  CHARACTER(LEN=1)  :: comment_character,&
                       continuation_character,&
                       section_character

  INTEGER, DIMENSION(2), PARAMETER :: reserved_unit_numbers = (/5,6/)

  CHARACTER(LEN=max_file_name_length) :: input_file_name
  CHARACTER(LEN=max_line_length)      :: input_line

  INTEGER :: group,icol,icol1,icol2,input_line_number,output_unit,source
  LOGICAL :: first_separator,ionode

! *** Public subroutines ***

  PUBLIC :: close_file,&
            finish_parser,&
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
    MODULE PROCEDURE read_integer_object,read_real_object,read_string_object
  END INTERFACE

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE broadcast_input_information(ifail)

!   Purpose: Broadcast the input information.

!   History: - Creation (02.03.2001, Matthias Krack)

!   ***************************************************************************

    INTEGER, INTENT(IN), OPTIONAL :: ifail

!   ---------------------------------------------------------------------------

    IF (PRESENT(ifail)) CALL mp_bcast(ifail,source,group)
    CALL mp_bcast(input_line,source,group)
    CALL mp_bcast(input_line_number,source,group)
    CALL mp_bcast(icol,source,group)

  END SUBROUTINE broadcast_input_information

! *****************************************************************************

  SUBROUTINE close_file(unit_number,file_status)

!   Purpose: Close the current input file.

!   History: - Creation (14.02.2001, Matthias Krack)

!   ***************************************************************************

    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: file_status
    INTEGER, INTENT(IN)                    :: unit_number

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "close_file"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

!   *** Local variables ***

    CHARACTER(LEN=max_message_length) :: message
    CHARACTER(LEN=6)                  :: status_string
    INTEGER                           :: istat
    LOGICAL                           :: exists,opened

!   ---------------------------------------------------------------------------

!   *** Check the specified input file name ***

    INQUIRE (UNIT=unit_number,EXIST=exists,OPENED=opened,IOSTAT=istat)

    IF (istat /= 0) THEN
      WRITE (UNIT=message,FMT="(A,I6,A,I6,A)")&
        "An error occurred inquiring the unit with the number ",unit_number,&
        " (IOSTAT = ",istat,")"
      CALL compress(message)
      CALL stop_program(routine,message)
    ELSE IF (.NOT.exists) THEN
      WRITE (UNIT=message,FMT="(A,I6,A,I6,A)")&
        "The specified unit number ",unit_number," does not exist"
      CALL compress(message)
      CALL stop_program(routine,message)
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
        CALL stop_program(routine,message)
      END IF

    END IF

  END SUBROUTINE close_file

! *****************************************************************************

  SUBROUTINE finish_parser()

!   Purpose: Finish a parser run.

!   History: - Creation (14.02.2001, Matthias Krack)

!   ***************************************************************************

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "finish_parser"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

!   *** Local variables ***

    CHARACTER(LEN=max_message_length) :: message

!   ---------------------------------------------------------------------------

    IF (ionode) THEN
      IF (input_unit == 0) THEN
        CALL stop_program(routine,&
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

!   History: - Creation (22.11.1999, Matthias Krack)

!   ***************************************************************************

    INTEGER, INTENT(IN)            :: nline
    INTEGER, INTENT(OUT), OPTIONAL :: ifail

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

  SUBROUTINE get_next_input_string()

!   Purpose: Get the next input string from the input line.

!   History: - Creation (19.02.2001, Matthias Krack)

!   ***************************************************************************

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "get_next_input_string"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

!   *** Local variables ***

    CHARACTER(LEN=max_message_length) :: message

!   ---------------------------------------------------------------------------

!   *** Search for the beginning of the next input string ***

    DO

!     *** Increment the column counter ***

      icol = icol + 1

!     *** Quick return, if the end of line is found ***

      IF ((icol > LEN_TRIM(input_line)).OR.&
          (input_line(icol:icol) == comment_character)) THEN
        icol1 = 0
        icol2 = -1
        RETURN
      END IF

!     *** Check for input line continuation ***

      IF ((icol == LEN_TRIM(input_line)).AND.&
          (input_line(icol:icol) == continuation_character)) THEN
        CALL get_next_input_line(1)
        CYCLE
      END IF

!     *** Ignore all white space and accept only one     ***
!     *** separator token or a string in quotation marks ***

      IF (input_line(icol:icol) == blank_character) THEN
        CYCLE
      ELSE IF (INDEX(separators,input_line(icol:icol)) > 0) THEN
        IF (first_separator) THEN
          first_separator = .FALSE.
          CYCLE
        ELSE
          icol1 = icol
          icol2 = icol
          CALL stop_parser(routine,&
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
          CALL stop_parser(routine,"Unmatched quotation mark found")
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

!   *** Search for the end of the next input string ***

    DO

      icol = icol + 1

      IF ((icol > LEN_TRIM(input_line)).OR.&
          (input_line(icol:icol) == blank_character).OR.&
          (input_line(icol:icol) == comment_character).OR.&
          (input_line(icol:icol) == continuation_character)) THEN
        EXIT
      ELSE IF (INDEX(separators,input_line(icol:icol)) > 0) THEN
        first_separator = .FALSE.
        EXIT
      END IF

    END DO

    icol2 = icol - 1

    IF (input_line(icol:icol) == continuation_character) icol = icol2

  END SUBROUTINE get_next_input_string

! *****************************************************************************

  FUNCTION get_unit_number() RESULT(unit_number)

!   Purpose: Returns a unit number of a unit that exists and is not connected.

!   History: - Creation (13.02.2001, Matthias Krack)

!   ***************************************************************************

    INTEGER :: unit_number

!   *** Local variables ***

    INTEGER :: istat
    LOGICAL :: exists,opened

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

!   History: - Creation (22.11.1999, Matthias Krack)

!   ***************************************************************************

    CHARACTER(LEN=*), INTENT(IN) :: string

    LOGICAL :: contains_integer_object

!   *** Local variables ***

    INTEGER :: i,length

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

  SUBROUTINE open_file(file_name,file_status,file_form,file_action,unit_number)

!   Purpose: Open the file with the name "file_name".

!   History: - Creation (13.02.2001, Matthias Krack)

!   ***************************************************************************

    CHARACTER(LEN=*), INTENT(IN)           :: file_name
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: file_action,file_form,file_status
    INTEGER, INTENT(OUT)                   :: unit_number

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "open_file"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

!   *** Local variables ***

    CHARACTER(LEN=max_message_length) :: message
    CHARACTER(LEN=11)                 :: action_string,form_string,&
                                         status_string
    INTEGER                           :: istat
    LOGICAL                           :: exists,opened

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

!   *** Check the specified input file name ***

    INQUIRE (FILE=TRIM(file_name),EXIST=exists,OPENED=opened,IOSTAT=istat)

    IF (istat /= 0) THEN
      WRITE (UNIT=message,FMT="(A,I6,A)")&
        "An error occurred inquiring the file <"//TRIM(file_name)//&
        "> (IOSTAT = ",istat,")"
      CALL stop_program(routine,message)
    ELSE IF (status_string == "OLD") THEN
      IF (.NOT.exists) THEN
        CALL stop_program(routine,&
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
        CALL stop_program(routine,&
                          "Could not find a logical unit number which "//&
                          "exists and is not connected")
      END IF

      OPEN (UNIT=unit_number,FILE=TRIM(file_name),STATUS=TRIM(status_string),&
            ACCESS="SEQUENTIAL",FORM=TRIM(form_string),POSITION="REWIND",&
            ACTION=TRIM(action_string),IOSTAT=istat)

      IF (istat /= 0) THEN
        WRITE (UNIT=message,FMT="(A,I6,A,I6,A)")&
          "An error occurred opening the file <"//TRIM(file_name)//&
          "> with the unit number ",unit_number," (IOSTAT = ",istat,")"
        CALL compress(message)
        CALL stop_program(routine,message)
      END IF

    END IF

    input_file_name = TRIM(file_name)

  END SUBROUTINE open_file

! *****************************************************************************

  SUBROUTINE read_integer_object(object,lower_to_upper,newline,skip_lines)

!   Purpose: Read an integer number.

!   History: - Creation (22.11.1999, Matthias Krack)

!   ***************************************************************************

    INTEGER, INTENT(OUT)          :: object
    INTEGER, INTENT(IN), OPTIONAL :: skip_lines
    LOGICAL, INTENT(IN), OPTIONAL :: lower_to_upper,newline

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "read_integer_object"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

!   *** Local variables ***

    INTEGER :: nline

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

    CALL get_next_input_string()

    IF (icol1 > icol2) THEN
      icol1 = icol
      icol2 = icol
      CALL stop_parser(routine,&
                       "An integer type object was expected, "//&
                       "found end of line")
    END IF

    IF (integer_object(input_line(icol1:icol2))) THEN
      READ (UNIT=input_line(icol1:icol2),FMT=*) object
    ELSE
      CALL stop_parser(routine,&
                       "An integer type object was expected, found <"//&
                       input_line(icol1:icol2)//">")
    END IF

  END SUBROUTINE read_integer_object

! *****************************************************************************

  SUBROUTINE read_next_input_line(nline,ifail)

!   Purpose: Read the next line from a logical unit "unit" (I/O node only).
!            Skip (nline-1) lines and skip also all comment lines.

!   History: - Creation (22.11.1999, Matthias Krack)

!   ***************************************************************************

    INTEGER, INTENT(IN)            :: nline
    INTEGER, INTENT(OUT), OPTIONAL :: ifail

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "read_next_input_line"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

!   *** Local variables ***

    CHARACTER(LEN=max_message_length)   :: message
    INTEGER                             :: icomment,iline,istat

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
            CALL stop_program(routine,message)
          END IF
        END IF

!       *** Check for comments ***

        icomment = INDEX(input_line,comment_character)

!       *** Ignore empty lines ***

        IF (icomment > 0) THEN
          IF (LEN_TRIM(input_line(:icomment-1)) == 0) iline = iline - 1
        ELSE
          IF (LEN_TRIM(input_line) == 0) iline = iline - 1
        END IF

      END DO

!     *** Reset column pointer, if a new line was read ***

      IF (nline > 0) icol = 0

    END IF

  END SUBROUTINE read_next_input_line

! *****************************************************************************

  SUBROUTINE read_real_object(object,lower_to_upper,newline,skip_lines)

!   Purpose: Read a floating point number.

!   History: - Creation (22.11.1999, Matthias Krack)

!   ***************************************************************************

    REAL(wp), INTENT(OUT)         :: object
    INTEGER, INTENT(IN), OPTIONAL :: skip_lines
    LOGICAL, INTENT(IN), OPTIONAL :: lower_to_upper,newline

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "read_real_object"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

!   *** Local variables ***

    INTEGER :: idot_first,idot_last,islash,istat,iz,nline,nz

    CHARACTER(LEN=max_line_length), DIMENSION(2) :: string
    REAL(wp), DIMENSION(2)                       :: z

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

    CALL get_next_input_string()

    IF (icol1 > icol2) THEN
      icol1 = icol
      icol2 = icol
      CALL stop_parser(routine,&
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
      z(2) = 1.0_wp
    END IF

    DO iz=1,nz

      IF (LEN_TRIM(string(iz)) == 0) THEN
        CALL stop_parser(routine,&
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
          CALL stop_parser(routine,&
                           "A real type object was expected, found <"//&
                           TRIM(string(iz))//">")
        END IF
      ELSE IF (integer_object(string(iz))) THEN
        READ (UNIT=string(iz),FMT=*) z(iz)
      ELSE
        CALL stop_parser(routine,&
                         "A real type object was expected, found <"//&
                         TRIM(string(iz))//">")
      END IF

    END DO

    object = z(1)/z(2)

  END SUBROUTINE read_real_object

! *****************************************************************************

  SUBROUTINE read_string_object(object,lower_to_upper,newline,skip_lines)

!   Purpose: Read a string.

!   History: - Creation (22.11.1999, Matthias Krack)

!   ***************************************************************************

    CHARACTER(LEN=*), INTENT(OUT) :: object
    INTEGER, INTENT(IN), OPTIONAL :: skip_lines
    LOGICAL, INTENT(IN), OPTIONAL :: lower_to_upper,newline

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "read_string_object"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

!   *** Local variables ***

    CHARACTER(LEN=max_message_length) :: message
    INTEGER                           :: input_string_length,nline

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

    CALL get_next_input_string()

    input_string_length = icol2 - icol1 + 1

    IF (input_string_length == 0) THEN
      icol1 = icol
      icol2 = icol
      CALL stop_parser(routine,&
                       "A string type object was expected, found end of line")
    ELSE IF (input_string_length > LEN(object)) THEN
      WRITE (UNIT=message,FMT="(A,I6,A)")&
        "The input string <"//input_line(icol1:icol2)//"> has more than ",&
        LEN(object)," characters and is therefore too long to fit in the "//&
        "specified variable"
      CALL compress(message)
      CALL stop_parser(routine,message)
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

!   History: - Creation (05.10.1999, Matthias Krack)

!   ***************************************************************************

    CHARACTER(LEN=*), INTENT(IN)            :: string
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: line
    LOGICAL, INTENT(IN)                     :: ignore_case
    LOGICAL, INTENT(OUT)                    :: found

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "search_object"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

!   *** Local variables ***

    CHARACTER(LEN=max_line_length+1)  :: current_line
    CHARACTER(LEN=max_message_length) :: message
    CHARACTER(LEN=LEN(string))        :: pattern
    INTEGER                           :: ifail,ipattern

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
              CALL stop_program(routine,message)
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

!   History: - Creation (14.02.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    CHARACTER(LEN=*), INTENT(IN)              :: file_name
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: end_section_label,&
                                                 separator_chars,&
                                                 start_section_label,&
                                                 start_section_name
    CHARACTER(LEN=1), OPTIONAL, INTENT(IN)    :: comment_char,&
                                                 continuation_char,&
                                                 section_char
    LOGICAL, OPTIONAL, INTENT(OUT)            :: section_found

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "start_parser"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

!   *** Local variables ***

    CHARACTER(LEN=max_message_length) :: message
    CHARACTER(LEN=60)                 :: current_name,requested_name
    CHARACTER(LEN=3)                  :: test_result
    LOGICAL                           :: found

!   ---------------------------------------------------------------------------

    ionode = globenv%ionode
    group = globenv%group
    source = globenv%source

!   *** Load the default values ***

    separators = default_separators
    comment_character = default_comment_character
    continuation_character = default_continuation_character
    section_character = default_section_character
    end_section = section_character//default_end_section_label

!   *** Overwrite the default values, if requested ***

    IF (PRESENT(separator_chars)) separators = separator_chars
    IF (PRESENT(comment_char)) comment_character = comment_char
    IF (PRESENT(continuation_char)) continuation_character = continuation_char
    IF (PRESENT(section_char)) section_character = section_char
    IF (PRESENT(end_section_label)) THEN
      end_section = section_character//TRIM(end_section_label)
    END IF

!   *** Get the logical output unit number for error messages ***

    CALL get_error_unit(output_unit)

    IF (ionode) THEN
      IF (input_unit == 0) THEN
        CALL open_file(file_name=TRIM(file_name),unit_number=input_unit)
      ELSE
        WRITE (UNIT=output_unit,FMT="(/,T2,A,I6,/,T2,A)")&
          "Current input unit: ",input_unit,&
          "Current input file: """//TRIM(input_file_name)//""""
        CALL stop_program(routine,&
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
    END IF

    IF (PRESENT(section_found)) section_found = found

  END SUBROUTINE start_parser

! *****************************************************************************

  SUBROUTINE stop_parser(routine,message)

!   Purpose: Print an input error message and stop the program run.

!            There are default error messages for common instances:

!            - END_OF_FILE or EOF
!            - INVALID_KEYWORD
!            - SECTION_NOT_FOUND

!   History: - Creation (15.02.2001, Matthias Krack)

!   ***************************************************************************

    CHARACTER(LEN=*), INTENT(IN) :: message,routine

!   *** Local variables ***

    CHARACTER(LEN=max_message_length) :: error_message

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

  END SUBROUTINE stop_parser

! *****************************************************************************

  FUNCTION test_object(newline,skip_lines) RESULT(test_result)

!   Purpose: Test next input object.

!            test_result = "EOF": End of file
!            test_result = "EOL": End of line
!            test_result = "EOS": End of section
!            test_result = "FLT": Floating point number
!            test_result = "INT": Integer number
!            test_result = "STR": String

!   History: - Creation (23.11.1999, Matthias Krack)

!   ***************************************************************************

    INTEGER, INTENT(IN), OPTIONAL :: skip_lines
    LOGICAL, INTENT(IN), OPTIONAL :: newline

    CHARACTER(LEN=3) :: test_result

!   *** Local variables ***

    CHARACTER(LEN=max_line_length) :: old_input_line
    REAL(wp)                       :: z
    INTEGER                        :: idot_first,idot_last,ifail,iline,islash,&
                                      istat,iz,nline,nz,old_icol,old_icol1,&
                                      old_icol2,old_input_line_number

    CHARACTER(LEN=max_line_length), DIMENSION(2) :: string

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

    CALL get_next_input_string()

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
