!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004 CP2K developer group
!-----------------------------------------------------------------------------!
!!****** cp2k/cp_parser [1.0] *
!!
!!   NAME
!!     cp_parser
!!
!!   FUNCTION
!!     Utility routines to read data from files.
!!     Kept as close as possible to the old parser because 
!!       1. string handling is a weak point of fortran compilers, and it is
!!          easy to write correct things that do not work
!!       2. conversion of old code
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     22.11.1999 first version of the old parser (called qs_parser)
!!                Matthias Krack
!!     06.2004 removed module varibles, cp_parser_type, new module [fawzi]
!!
!!   SOURCE
!******************************************************************************

MODULE cp_parser
  USE cp_files,                        ONLY: open_file, close_file
  USE cp_para_env,                     ONLY: cp_para_env_create,&
                                             cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE message_passing,                 ONLY: MPI_COMM_SELF,&
                                             mp_bcast
  USE string_utilities,                ONLY: compress,&
                                             uppercase
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cp_parser_type
  PUBLIC :: parser_release, parser_retain, parser_create
  PUBLIC :: parser_next_token, parser_test_next_token
  PUBLIC :: parser_get_object, parser_get_real, parser_get_integer,&
       parser_get_logical, parser_get_string, parser_location,&
       parser_search_string,parser_get_next_line

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: moduleN = "cp_parser"

  INTEGER, PARAMETER, PUBLIC :: max_line_length = 240

! *** Global variables ***

  CHARACTER(LEN=4),PARAMETER :: default_separators = ",:;="
  CHARACTER(LEN=3),PARAMETER :: default_end_section_label = "END"
  CHARACTER(LEN=1),PARAMETER :: blank_character = " ",&
                      default_comment_character = "#",&
                      default_continuation_character = CHAR(92),& ! backslash
                      default_section_character = "&",&
                      horizontal_tab = CHAR(9)
  INTEGER,PARAMETER, PUBLIC :: max_unit_number = 999

  INTEGER, SAVE, PRIVATE :: last_parser_id=0

!!****s* cp_parser/cp_parser_type [1.0] *
!!
!!   NAME
!!     cp_parser_type
!!
!!   FUNCTION
!!     represent a parser
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     -
! icol             : Number of the current column in the current input line,
! -1 if at the end of the file
! icol1            : First column of the current input string
! icol2            : Last column of the current input string
! input_line_number: Number of the current input line read from the input file
! input_unit       : Logical unit number of the input file
!!
!!   AUTHOR
!!     fawzi
!!
!!   SOURCE
  !***************************************************************************
  TYPE cp_parser_type
     INTEGER :: id_nr,ref_count,input_unit    
     CHARACTER(LEN=default_string_length) :: end_section,&
          start_section
     CHARACTER(LEN=10) :: separators
     CHARACTER(LEN=1)  :: comment_character,&
          continuation_character,&
          section_character
     CHARACTER(LEN=default_path_length) :: input_file_name
     CHARACTER(LEN=max_line_length)      :: input_line
     INTEGER :: icol,icol1,icol2,input_line_number
     LOGICAL :: first_separator
     TYPE(cp_para_env_type), POINTER :: para_env
  END TYPE cp_parser_type
!!***
  !****************************************************************************

! *****************************************************************************
  INTERFACE parser_get_object
    MODULE PROCEDURE parser_get_integer,&
                     parser_get_logical,&
                     parser_get_real,&
                     parser_get_string
  END INTERFACE
!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE broadcast_input_information(parser,error)

!   Purpose: Broadcast the input information.

!   History: - Creation (02.03.2001,MK)

!   ***************************************************************************
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = "broadcast_input_information", &
      routineP = moduleN//":"//routineN

    LOGICAL                                  :: failure
    TYPE(cp_para_env_type), POINTER          :: para_env

!   ---------------------------------------------------------------------------

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(parser),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CPPrecondition(parser%ref_count>0,cp_failure_level,routineP,error,failure)
       para_env => parser%para_env
    END IF
    IF (.not.failure .AND. para_env%num_pe>1) THEN
       CALL mp_bcast(parser%input_line,para_env%source,para_env%group)
       CALL mp_bcast(parser%input_line_number,para_env%source,para_env%group)
       CALL mp_bcast(parser%icol,para_env%source,para_env%group)
       CALL cp_error_synchronize_error(error,para_env=para_env)
    END IF

  END SUBROUTINE broadcast_input_information

! *****************************************************************************

!!****f* cp_parser/parser_retain *
!!
!!   NAME
!!     parser_retain
!!
!!   FUNCTION
!!     retains the given parser
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - parser: the parser to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE parser_retain(parser,error)
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'parser_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  CPPrecondition(ASSOCIATED(parser),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPreconditionNoFail(parser%ref_count>0,cp_failure_level,routineP,error)
     parser%ref_count=parser%ref_count+1
  END IF
END SUBROUTINE parser_retain
!***************************************************************************


  SUBROUTINE parser_release(parser, error)

!   Purpose: releases the parser

!   History: - Creation (14.02.2001,MK)

!   ***************************************************************************
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "parser_release", &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(parser)) THEN
       CPPreconditionNoFail(parser%ref_count>0,cp_failure_level,routineP,error)
       parser%ref_count=parser%ref_count-1
       IF (parser%ref_count==0) THEN
          IF (parser%input_unit >= 0) THEN
             CALL close_file(unit_number=parser%input_unit)
          END IF
          CALL cp_para_env_release(parser%para_env)
          DEALLOCATE(parser,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF
    NULLIFY(parser)
  END SUBROUTINE parser_release

! *****************************************************************************

  SUBROUTINE parser_get_next_line(parser,nline,at_end,error)

!   Purpose: Read the next input line and broadcast the input information.
!            Skip (nline-1) lines and skip also all comment lines.

!   History: - Creation (22.11.1999,MK)

!   ***************************************************************************
    
    TYPE(cp_parser_type), POINTER            :: parser
    INTEGER, INTENT(IN)                      :: nline
    LOGICAL, INTENT(out), OPTIONAL           :: at_end
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = "parser_get_next_line", &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: my_at_end

!   ---------------------------------------------------------------------------

    IF (nline>0) THEN
       CALL parser_read_line(parser,nline,at_end=my_at_end,error=error)
       CALL mp_bcast(my_at_end, parser%para_env%source,&
            parser%para_env%group)
       IF (PRESENT(at_end)) THEN
          at_end=my_at_end
       ELSE
          CALL cp_assert(.not.my_at_end,cp_failure_level,cp_assertion_failed,&
               routineP, "unexpected EOF "//&
               TRIM(parser_location(parser,error=error)), error)
       END IF
       CALL broadcast_input_information(parser,error=error)
    ELSEIF (PRESENT(at_end)) THEN
       at_end=.FALSE.
    END IF

  END SUBROUTINE parser_get_next_line

! *****************************************************************************

!!****f* cp_parser/parser_location *
!!
!!   NAME
!!     parser_location
!!
!!   FUNCTION
!!     return a description of the part of the file acually parsed
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - parser: the parser
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
FUNCTION parser_location(parser,error) RESULT(res)
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    CHARACTER(len=default_path_length+&
      default_string_length)                 :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'parser_location', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  CPPrecondition(ASSOCIATED(parser),cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     CPPrecondition(parser%ref_count>0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT. failure) THEN
     res="file:'"//TRIM(parser%input_file_name)//"' line:"//&
          cp_to_string(parser%input_line_number)//" col:"//&
          cp_to_string(parser%icol)
     IF (parser%icol==-1) THEN
        res(LEN_TRIM(res):)=" (EOF)"
     ELSE IF (MAX(1,parser%icol1)<=parser%icol2) THEN
        res(LEN_TRIM(res):)=" chunk:'"//&
             parser%input_line(MAX(1,parser%icol1):parser%icol2)//"'"
     END IF
  END IF
END FUNCTION parser_location
!***************************************************************************

SUBROUTINE parser_next_token(parser,string_length,error)

  !   Purpose: Get the next input string from the input line.

  !   History: - Creation (19.02.2001,MK)

  !   ***************************************************************************
    TYPE(cp_parser_type), POINTER            :: parser
    INTEGER, INTENT(IN), OPTIONAL            :: string_length
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "parser_next_token", &
      routineP = moduleN//":"//routineN

    INTEGER                                  :: i, length
    LOGICAL                                  :: failure

!   ---------------------------------------------------------------------------

  failure=.FALSE.
  CPPrecondition(ASSOCIATED(parser),cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     CPPrecondition(parser%ref_count>0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT.failure) THEN
     IF (PRESENT(string_length)) THEN
        length = string_length
     ELSE
        length = 0
     END IF

     IF (length > 0) THEN

        !     *** Fixed input string length ***

        CALL cp_assert(parser%icol/=-1,cp_failure_level,cp_assertion_failed,&
             routineP,"Unexpectetly reached EOF "//&
             TRIM(parser_location(parser,error=error)),error,failure)
        length=MIN(LEN_TRIM(parser%input_line)-parser%icol1+1,length)
        IF (.NOT.failure) THEN
           parser%icol1 = parser%icol + 1
           parser%icol2 = parser%icol + length
           parser%icol = parser%icol2
        END IF
     ELSE

        !     *** Variable input string length (automatic search) ***

        ! check for EOF
        IF (parser%icol==-1) THEN
           parser%icol1 = 1
           parser%icol2 = -1
           RETURN
        END IF

        !     *** Search for the beginning of the next input string ***

        DO

           !       *** Increment the column counter ***

           parser%icol = parser%icol + 1

           !       *** Quick return, if the end of line is found ***

           IF ((parser%icol > LEN_TRIM(parser%input_line)).OR.&
                (parser%input_line(parser%icol:parser%icol) == parser%comment_character)) THEN
              parser%icol1 = 1
              parser%icol2 = -1
              RETURN
           END IF

           !       *** Check for input line continuation ***

           IF (parser%input_line(parser%icol:parser%icol) == &
                parser%continuation_character) THEN
              DO i=LEN(parser%input_line),1,-1
                 IF ((parser%input_line(i:i) /= blank_character).AND.&
                      (parser%input_line(i:i) /= horizontal_tab)) EXIT
              END DO
              IF (parser%icol == i.or.parser%input_line(i:i)==&
                   parser%comment_character) THEN
                 CALL parser_get_next_line(parser,1,error=error) 
                 ! does not pass at_end, it is an error to continue to EOF, change this behavoiur?
                 CYCLE
              ELSE
                 parser%icol1 = i
                 parser%icol2 = i
                 CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                      routineP, "Found non-blank tokens after the "//&
                      "line continuation character '"// &
                      parser%continuation_character//&
                      "' "//TRIM(parser_location(parser,error=error)),error,failure)
                 EXIT
              END IF
           END IF
           !       *** Ignore all white space and accept only one     ***
           !       *** separator token or a string in quotation marks ***

           IF ((parser%input_line(parser%icol:parser%icol) == blank_character).OR.&
                (parser%input_line(parser%icol:parser%icol) == horizontal_tab)) THEN
              CYCLE
           ELSE IF (INDEX(parser%separators,parser%input_line(parser%icol:parser%icol)) > 0) THEN
              IF (parser%first_separator) THEN
                 parser%first_separator = .FALSE.
                 CYCLE
              ELSE
                 parser%icol1 = parser%icol
                 parser%icol2 = parser%icol
                 CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                      routineP, "Unexpected separator token <"//&
                      parser%input_line(parser%icol:parser%icol)//"> found"//&
                      TRIM(parser_location(parser,error=error)),error,failure)
                 EXIT
              END IF
           ELSE IF (parser%input_line(parser%icol:parser%icol) == '"') THEN
              parser%first_separator = .TRUE.
              parser%icol1 = parser%icol + 1
              parser%icol2 = parser%icol + INDEX(parser%input_line(parser%icol1:),'"')
              IF (parser%icol2 == parser%icol) THEN
                 parser%icol1 = parser%icol
                 parser%icol2 = parser%icol
                 CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                      routineP, "Unmatched quotation mark found"//&
                      TRIM(parser_location(parser,error=error)),error,failure)
              ELSE
                 parser%icol = parser%icol2
                 parser%icol2 = parser%icol2 - 1
                 RETURN
              END IF
           ELSE
              parser%first_separator = .TRUE.
              parser%icol1 = parser%icol
              EXIT
           END IF

        END DO

        !     *** Search for the end of the next input string ***
        IF(.NOT.failure) THEN
           DO
              IF (parser%icol > LEN_TRIM(parser%input_line)) EXIT
              parser%icol = parser%icol + 1
              IF ((parser%icol > LEN_TRIM(parser%input_line)).OR.&
                   (parser%input_line(parser%icol:parser%icol) == blank_character).OR.&
                   (parser%input_line(parser%icol:parser%icol) == parser%comment_character).OR.&
                   (parser%input_line(parser%icol:parser%icol) == parser%continuation_character).OR.&
                   (parser%input_line(parser%icol:parser%icol) == horizontal_tab)) THEN
                 EXIT
              ELSE IF (INDEX(parser%separators,parser%input_line(parser%icol:parser%icol)) > 0) THEN
                 parser%first_separator = .FALSE.
                 EXIT
              END IF
           END DO

           parser%icol2 = parser%icol - 1

           IF (parser%input_line(parser%icol:parser%icol) == &
                parser%continuation_character) parser%icol = parser%icol2

        END IF
     END IF
  END IF

END SUBROUTINE parser_next_token

! *****************************************************************************


  FUNCTION integer_object(string,error) RESULT(contains_integer_object)

!   Purpose: Check, if the string object contains an object of type integer.

!   History: - Creation (22.11.1999,MK)

!   ***************************************************************************

    CHARACTER(LEN=*), INTENT(IN)             :: string
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
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

  SUBROUTINE parser_get_integer(parser,object,lower_to_upper,newline,skip_lines,&
       string_length, at_end,error)

    !   Purpose: Read an integer number.

    !   History: - Creation (22.11.1999,MK)

    !   ***************************************************************************

    TYPE(cp_parser_type), POINTER            :: parser
    INTEGER, INTENT(OUT)                     :: object
    LOGICAL, INTENT(IN), OPTIONAL            :: lower_to_upper, newline
    INTEGER, INTENT(IN), OPTIONAL            :: skip_lines, string_length
    LOGICAL, INTENT(out), OPTIONAL           :: at_end
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "parser_get_integer", &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nline
    LOGICAL                                  :: failure, my_at_end

!   ---------------------------------------------------------------------------

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(parser),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CPPrecondition(parser%ref_count>0,cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT.failure) THEN
       IF (PRESENT(skip_lines)) THEN
          nline = skip_lines
       ELSE
          nline = 0
       END IF

       IF (PRESENT(newline)) THEN
          IF (newline) nline = nline + 1
       END IF

       CALL parser_get_next_line(parser,nline,at_end=my_at_end,error=error)
       IF (PRESENT(at_end)) THEN
          at_end=my_at_end
          IF (my_at_end) RETURN
       ELSE IF (my_at_end) THEN
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP, "unexpected EOF "//&
               TRIM(parser_location(parser,error=error)),error,failure)
       END IF

       IF (PRESENT(string_length)) THEN
          CALL parser_next_token(parser,string_length,error=error)
       ELSE
          CALL parser_next_token(parser,error=error)
       END IF

       IF (parser%icol1 > parser%icol2) THEN
          parser%icol1 = parser%icol
          parser%icol2 = parser%icol
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP, "An integer type object was expected, "//&
               "found end of line"// TRIM(parser_location(parser,error=error)),&
               error,failure)
       END IF
    END IF
    IF (.NOT.failure) THEN
       IF (integer_object(parser%input_line(parser%icol1:parser%icol2))) THEN
          READ (UNIT=parser%input_line(parser%icol1:parser%icol2),FMT=*) object
       ELSE
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP, "An integer type object was expected, found <"//&
               parser%input_line(parser%icol1:parser%icol2)//">"//&
               TRIM(parser_location(parser,error=error)),error,failure)
       END IF
    END IF
  END SUBROUTINE parser_get_integer

! *****************************************************************************

  SUBROUTINE parser_get_logical(parser,object,newline,skip_lines,string_length,&
       at_end,error)

    !   Purpose: Read a string representing logical object.

    !   History: - Creation (01.04.2003,FM)
    !            - New version (08.07.2003,MK)

    !   ***************************************************************************

    TYPE(cp_parser_type), POINTER            :: parser
    LOGICAL, INTENT(OUT)                     :: object
    LOGICAL, INTENT(IN), OPTIONAL            :: newline
    INTEGER, INTENT(IN), OPTIONAL            :: skip_lines, string_length
    LOGICAL, INTENT(out), OPTIONAL           :: at_end
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "parser_get_logical", &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=max_line_length)           :: input_string
    INTEGER                                  :: input_string_length, nline
    LOGICAL                                  :: failure, my_at_end

!   ---------------------------------------------------------------------------

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(parser),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CPPrecondition(parser%ref_count>0,cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT.failure) THEN
       IF (PRESENT(skip_lines)) THEN
          nline = skip_lines
       ELSE
          nline = 0
       END IF

       IF (PRESENT(newline)) THEN
          IF (newline) nline = nline + 1
       END IF

       CALL parser_get_next_line(parser,nline,at_end=my_at_end,error=error)
       IF (PRESENT(at_end)) THEN
          at_end=my_at_end
          IF (my_at_end) RETURN
       ELSE IF (my_at_end) THEN
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP, "unexpected EOF "//&
               TRIM(parser_location(parser,error=error)), error,failure)
       END IF

       IF (PRESENT(string_length)) THEN
          CALL parser_next_token(parser,string_length,error=error)
       ELSE
          CALL parser_next_token(parser,error=error)
       END IF

       input_string_length = parser%icol2 - parser%icol1 + 1

       IF (input_string_length == 0) THEN
          parser%icol1 = parser%icol
          parser%icol2 = parser%icol
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP, "A string representing a logical object was expected, "//&
               "found end of line"//TRIM(parser_location(parser,error=error)),&
               error,failure)
       ELSE
          input_string=""
          input_string(:input_string_length) = parser%input_line(parser%icol1:parser%icol2)
       END IF
    END IF
    IF (.NOT.failure) THEN
       CALL uppercase(input_string)

       SELECT CASE (TRIM(input_string))
       CASE ("0","F",".F.","FALSE",".FALSE.","N","NO","OFF")
          object = .FALSE.
       CASE ("1","T",".T.","TRUE",".TRUE.","Y","YES","ON")
          object = .TRUE.
       CASE DEFAULT
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"A string representing a logical object was expected, "//&
               "found <"//TRIM(input_string)//">"//&
               TRIM(parser_location(parser,error=error)),error,failure)
       END SELECT
    END IF

  END SUBROUTINE parser_get_logical

! *****************************************************************************

  SUBROUTINE parser_read_line(parser,nline,at_end,error)

    !   Purpose: Read the next line from a logical unit "unit" (I/O node only).
    !            Skip (nline-1) lines and skip also all comment lines.

    !   History: - Creation (22.11.1999,MK)

    !   ***************************************************************************

    TYPE(cp_parser_type), POINTER            :: parser
    INTEGER, INTENT(IN)                      :: nline
    LOGICAL, INTENT(out), OPTIONAL           :: at_end
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "parser_read_line", &
      routineP = moduleN//":"//routineN

    INTEGER                                  :: icomment, iline, istat
    LOGICAL                                  :: failure

!   ---------------------------------------------------------------------------

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(parser),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CPPrecondition(parser%ref_count>0,cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT.failure) THEN
       IF (PRESENT(at_end)) at_end=.FALSE.
       IF (parser%para_env%mepos==parser%para_env%source) THEN

          IF (parser%icol==-1.and.nline==0) THEN
             CALL cp_assert(PRESENT(at_end),cp_failure_level,cp_assertion_failed,&
                  routineP,"An unexpected EOF "//TRIM(parser_location(parser,error=error)),&
                  error,failure)
             IF (PRESENT(at_end)) at_end=.TRUE.
             RETURN
          END IF

          iline = 0
          istat = 0

          DO WHILE (iline /= nline)

             iline = iline + 1

             !       *** Increment the input line counter ***

             parser%input_line_number = parser%input_line_number + 1

             !       *** Read next line from file ***

             READ (UNIT=parser%input_unit,FMT="(A)",IOSTAT=istat)&
                  parser%input_line

             !       *** Handle read errors ***

             IF (istat /= 0) THEN
                IF (istat<0) THEN ! EOF/EOR is negative other errors positive
                   CALL cp_assert(PRESENT(at_end),cp_failure_level,cp_assertion_failed,&
                        routineP,"Unexpected EOF "//TRIM(parser_location(parser,error=error)),&
                        error,failure)
                   IF (PRESENT(at_end)) at_end=.TRUE.
                   parser%icol=-1
                   parser%icol1=0
                   parser%icol2=-1
                ELSE
                   CALL cp_assert(.FALSE.,cp_failure_level,istat,routineP,&
                        "An input/output error occurred "//TRIM(parser_location(parser,error=error))//&
                        "(IOSTAT = "//cp_to_string(istat)//")",error,failure)
                END IF
                RETURN
             END IF

             !       *** Check for comments ***

             icomment = INDEX(parser%input_line,parser%comment_character)

             !       *** Ignore empty lines ***

             IF (icomment > 0) THEN
                !MK xlf 8.1.1 bug -> len_trim returns -1 for an empty string
                !         IF (LEN_TRIM(parser%input_line(:icomment-1)) == 0) iline = iline - 1
                IF (LEN_TRIM(parser%input_line(:icomment-1)) <= 0) iline = iline - 1
             ELSE
                !MK xlf 8.1.1 bug -> len_trim returns -1 for an empty string
                !         IF (LEN_TRIM(parser%input_line) == 0) iline = iline - 1
                IF (LEN_TRIM(parser%input_line) <= 0) iline = iline - 1
             END IF

          END DO

          !     *** Reset column pointer, if a new line was read ***

          IF (nline > 0) parser%icol = 0

       END IF
    END IF

  END SUBROUTINE parser_read_line

! *****************************************************************************

  SUBROUTINE parser_get_real(parser,object,lower_to_upper,newline,skip_lines,&
       string_length,multiline,at_end,error)

    !   Purpose: Read a floating point number.

    !   History: - Creation (22.11.1999,MK)

    !   ***************************************************************************

    TYPE(cp_parser_type), POINTER            :: parser
    REAL(KIND=dp), INTENT(OUT)               :: object
    LOGICAL, INTENT(IN), OPTIONAL            :: lower_to_upper, newline
    INTEGER, INTENT(IN), OPTIONAL            :: skip_lines, string_length
    LOGICAL, INTENT(IN), OPTIONAL            :: multiline
    LOGICAL, INTENT(out), OPTIONAL           :: at_end
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "parser_get_real", &
      routineP = moduleN//":"//routineN

    CHARACTER(LEN=max_line_length), &
      DIMENSION(2)                           :: string
    INTEGER                                  :: idot_first, idot_last, &
                                                islash, istat, iz, nline, nz
    LOGICAL                                  :: failure, my_at_end, &
                                                my_multiline
    REAL(KIND=dp), DIMENSION(2)              :: z

!   ---------------------------------------------------------------------------

    failure=.FALSE.
    my_multiline=.FALSE.
    IF (PRESENT(multiline)) my_multiline=multiline
    CPPrecondition(ASSOCIATED(parser),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CPPrecondition(parser%ref_count>0,cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT.failure) THEN
       IF (PRESENT(skip_lines)) THEN
          nline = skip_lines
       ELSE
          nline = 0
       END IF

       IF (PRESENT(newline)) THEN
          IF (newline) nline = nline + 1
       END IF

       CALL parser_get_next_line(parser,nline,at_end=my_at_end,error=error)
       IF (PRESENT(at_end)) THEN
          at_end=my_at_end
          IF (my_at_end) RETURN
       ELSE IF (my_at_end) THEN
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP, "unexpected EOF "//&
               TRIM(parser_location(parser,error=error)), error,failure)
       END IF

       CALL parser_next_token(parser,string_length,error=error)

       IF (parser%icol1 > parser%icol2.and.my_multiline) THEN
          CALL parser_get_next_line(parser,1,at_end=my_at_end,error=error)
          IF (PRESENT(at_end)) THEN
             at_end=my_at_end
             IF (my_at_end) RETURN
          ELSE IF (my_at_end) THEN
             CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                  routineP, "unexpected EOF "//&
                  TRIM(parser_location(parser,error=error)), error,failure)
          END IF
          CALL parser_next_token(parser,string_length,error=error)
       END IF

       IF (parser%icol1 > parser%icol2) THEN
          parser%icol1 = parser%icol
          parser%icol2 = parser%icol
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP, "A real type object was expected, found then end of the line "//&
               TRIM(parser_location(parser,error=error)),error,failure)
       END IF
    END IF
    IF (.NOT.failure) THEN
       islash = parser%icol1 + INDEX(parser%input_line(parser%icol1:parser%icol2),"/") - 1

       IF (islash > parser%icol1) THEN
          nz = 2
          string(1) = parser%input_line(parser%icol1:islash-1)
          string(2) = parser%input_line(islash+1:parser%icol2)
       ELSE
          nz = 1
          string(1) = parser%input_line(parser%icol1:parser%icol2)
          z(2) = 1.0_dp
       END IF

       DO iz=1,nz

          IF (LEN_TRIM(string(iz)) == 0) THEN
             CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                  routineP, "A real type object was expected, found end of line "//&
                  TRIM(parser_location(parser,error=error)),error,failure)
             ! try harder to find the missing nr?
             RETURN
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
          ELSE IF (integer_object(string(iz))) THEN
             READ (UNIT=string(iz),FMT=*,IOSTAT=istat) z(iz)
          ELSE
             istat = 1
          END IF
          IF (istat /= 0) THEN
             CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                  routineP, "A real type object was expected, found <"//&
                  TRIM(string(iz))//">"//TRIM(parser_location(parser,error=error)),&
                  error,failure)
             RETURN
          END IF

       END DO

       object = z(1)/z(2)
    END IF
  END SUBROUTINE parser_get_real

! *****************************************************************************

  SUBROUTINE parser_get_string(parser,object,lower_to_upper,newline,skip_lines,&
                                string_length,at_end,error)

!   Purpose: Read a string.

!   History: - Creation (22.11.1999,MK)

!   ***************************************************************************

    TYPE(cp_parser_type), POINTER            :: parser
    CHARACTER(LEN=*), INTENT(OUT)            :: object
    LOGICAL, INTENT(IN), OPTIONAL            :: lower_to_upper, newline
    INTEGER, INTENT(IN), OPTIONAL            :: skip_lines, string_length
    LOGICAL, INTENT(out), OPTIONAL           :: at_end
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "parser_get_string", &
      routineP = moduleN//":"//routineN

    INTEGER                                  :: input_string_length, nline
    LOGICAL                                  :: failure, my_at_end

!   ---------------------------------------------------------------------------

    object = ""
    failure=.FALSE.
    CPPrecondition(ASSOCIATED(parser),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CPPrecondition(parser%ref_count>0,cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT.failure) THEN
       IF (PRESENT(skip_lines)) THEN
          nline = skip_lines
       ELSE
          nline = 0
       END IF
       
       IF (PRESENT(newline)) THEN
          IF (newline) nline = nline + 1
       END IF
       
       CALL parser_get_next_line(parser,nline,at_end=my_at_end,error=error)
       IF (PRESENT(at_end)) THEN
          at_end=my_at_end
          IF (my_at_end) RETURN
       ELSE IF (my_at_end) THEN
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP, "unexpected EOF "//&
               TRIM(parser_location(parser,error=error)), error,failure)
       END IF

       CALL parser_next_token(parser,string_length,error=error)

       input_string_length = parser%icol2 - parser%icol1 + 1
       
       IF (input_string_length <= 0) THEN
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
               "A string type object was expected, found end of line "//&
               TRIM(parser_location(parser,error=error)),error,failure)
       ELSE IF (input_string_length > LEN(object)) THEN
          CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
               "The input string <"//parser%input_line(parser%icol1:parser%icol2)//&
               "> has more than "//cp_to_string(LEN(object))//&
               " characters and is therefore too long to fit in the "//&
               "specified variable"//TRIM(parser_location(parser,error=error)),&
               error)
          object=parser%input_line(parser%icol1:parser%icol1+LEN(object)-1)
       ELSE
          object(:input_string_length) = parser%input_line(parser%icol1:parser%icol2)
       END IF

!   *** Convert lowercase to uppercase, if requested ***
       
       IF (PRESENT(lower_to_upper)) THEN
          IF (lower_to_upper) CALL uppercase(object)
       END IF
    END IF
  END SUBROUTINE parser_get_string

! *****************************************************************************

  SUBROUTINE parser_search_string(parser,string,ignore_case,found,line,begin_line,error)

    !   Purpose: Search a string pattern in a file defined by its logical unit
    !            number "unit". A case sensitive search is performed, if
    !            ignore_case is .FALSE..
    !            begin_line: give back the parser at the beginning of the line
    !            matching the search

    !   History: - Creation (05.10.1999,MK)

    !   ***************************************************************************

    TYPE(cp_parser_type), POINTER            :: parser
    CHARACTER(LEN=*), INTENT(IN)             :: string
    LOGICAL, INTENT(IN)                      :: ignore_case
    LOGICAL, INTENT(OUT)                     :: found
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: line
    LOGICAL, INTENT(IN), OPTIONAL            :: begin_line
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "parser_search_string", &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=LEN(string))               :: pattern
    CHARACTER(LEN=max_line_length+1)         :: current_line
    INTEGER                                  :: ipattern
    LOGICAL                                  :: at_end, failure, begin

!   ---------------------------------------------------------------------------

    found = .FALSE.
    failure = .FALSE.
    begin = .FALSE.
    IF (PRESENT(begin_line)) begin=begin_line

    CPPrecondition(ASSOCIATED(parser),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CPPrecondition(parser%ref_count>0,cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT.failure) THEN
       IF (PRESENT(line)) line = ""

       !   *** Search for string pattern ***

       IF (parser%para_env%mepos==parser%para_env%source) THEN

          pattern = string

          IF (ignore_case) CALL uppercase(pattern)

          DO

             CALL parser_read_line(parser,1,at_end=at_end,error=error)

             !       *** Exit loop, if the end of file is reached ***

             IF (at_end) EXIT

             !       *** Check the current line for string pattern ***

             current_line = parser%input_line

             IF (ignore_case) CALL uppercase(current_line)

             ipattern = INDEX(current_line,TRIM(pattern))

             IF (ipattern > 0) THEN
                found = .TRUE.
                parser%icol = ipattern - 1
                IF (PRESENT(line)) THEN
                   CALL cp_assert(LEN(line) > LEN_TRIM(parser%input_line),&
                        cp_warning_level,cp_assertion_failed,&
                        routineP, "The returned input line has more than "//&
                        cp_to_string(LEN(line))//&
                        " characters and is therefore too long to fit in the "//&
                        "specified variable"//&
                        TRIM(parser_location(parser,error=error)),error,&
                        failure)
                END IF
                EXIT
             END IF

          END DO

          IF (found) THEN
             IF (begin) parser%icol = 0
          END IF
       END IF
      
       IF (parser%para_env%num_pe>1) THEN
          CALL mp_bcast(found,parser%para_env%source,parser%para_env%group)
          CALL broadcast_input_information(parser,error)
          CALL cp_error_synchronize_error(error,para_env=parser%para_env)
       END IF

       IF (found) THEN
          IF (PRESENT(line)) line = parser%input_line
          IF (.NOT.begin) CALL parser_next_token(parser,error=error)
       END IF
    END IF
  END SUBROUTINE parser_search_string

! *****************************************************************************

  SUBROUTINE parser_create(parser,file_name,unit_nr,para_env,end_section_label,&
       separator_chars,comment_char,continuation_char, section_char,error)

    !   Purpose: Start a parser run.

    !   History: - Creation (14.02.2001,MK)

    !   ***************************************************************************

    TYPE(cp_parser_type), POINTER            :: parser
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: file_name
    INTEGER, INTENT(in), OPTIONAL            :: unit_nr
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: end_section_label, &
                                                separator_chars
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL   :: comment_char, &
                                                continuation_char, &
                                                section_char
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "parser_create", &
      routineP = moduleN//":"//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

!   ---------------------------------------------------------------------------

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(parser),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       ALLOCATE(parser,stat=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
    END IF
    IF (.NOT.failure) THEN
       last_parser_id=last_parser_id+1
       parser%id_nr=last_parser_id
       parser%ref_count=1

       parser%input_unit=-1
       parser%input_file_name=""

       !   *** Load the default values and overwrite them, if requested ***

       parser%separators = default_separators
       IF (PRESENT(separator_chars)) parser%separators = separator_chars
       parser%comment_character = default_comment_character
       IF (PRESENT(comment_char)) parser%comment_character = comment_char
       parser%continuation_character = default_continuation_character
       IF (PRESENT(continuation_char)) parser%continuation_character = continuation_char
       parser%section_character = default_section_character
       IF (PRESENT(section_char)) parser%section_character = section_char
       parser%end_section = parser%section_character//default_end_section_label
       IF (PRESENT(end_section_label)) THEN
          parser%end_section = parser%section_character//TRIM(end_section_label)
       END IF

       CALL compress(parser%end_section) ! needed?

       ! para_env
       IF (PRESENT(para_env)) THEN
          parser%para_env => para_env
          CALL cp_para_env_retain(para_env,error=error)
       ELSE
          NULLIFY(parser%para_env)
          CALL cp_para_env_create(parser%para_env, group=MPI_COMM_SELF, source=0,&
               mepos=0, num_pe=1, owns_group=.FALSE.,error=error)
       END IF

       !   *** Get the logical output unit number for error messages ***

       IF (parser%para_env%mepos==parser%para_env%source) THEN
          IF (PRESENT(unit_nr)) THEN
             parser%input_unit=unit_nr
             IF (PRESENT(file_name)) parser%input_file_name=file_name
          ELSE
             CALL cp_assert(PRESENT(file_name),cp_failure_level,cp_assertion_failed,&
                  routineP,"at least one of filename and unit_nr must be present"//&
                  CPSourceFileRef,&
                  error,failure)
             IF (.NOT.failure) THEN
                CALL open_file(file_name=TRIM(file_name),&
                     unit_number=parser%input_unit)
                parser%input_file_name=file_name
             END IF
          END IF
       END IF

       parser%input_line_number = 0
       parser%icol=0
       parser%icol1=0
       parser%icol2=0
       parser%first_separator=.TRUE.
    END IF
  END SUBROUTINE parser_create

! *****************************************************************************

  FUNCTION parser_test_next_token(parser,newline,skip_lines,string_length,error)&
       RESULT(test_result)

    !   Purpose: Test next input object.

    !            test_result = "EOF": End of file
    !            test_result = "EOL": End of line
    !            test_result = "EOS": End of section
    !            test_result = "FLT": Floating point number
    !            test_result = "INT": Integer number
    !            test_result = "STR": String

    !   History: - Creation (23.11.1999,MK)

    !   ***************************************************************************

    TYPE(cp_parser_type), POINTER            :: parser
    LOGICAL, INTENT(IN), OPTIONAL            :: newline
    INTEGER, INTENT(IN), OPTIONAL            :: skip_lines, string_length
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    CHARACTER(LEN=3)                         :: test_result

    CHARACTER(len=*), PARAMETER :: routineN = "parser_test_next_token", &
      routineP = moduleN//":"//routineN

    CHARACTER(LEN=max_line_length)           :: old_input_line
    CHARACTER(LEN=max_line_length), &
      DIMENSION(2)                           :: string
    INTEGER :: idot_first, idot_last, iline, islash, istat, iz, nline, nz, &
      old_icol, old_icol1, old_icol2, old_input_line_number
    LOGICAL                                  :: at_end, failure
    REAL(KIND=dp)                            :: z
    TYPE(cp_error_type)                      :: suberror

!   ---------------------------------------------------------------------------

    failure=.FALSE.
    test_result = ""
    CALL cp_error_init(suberror,template_error=error)

    CPPrecondition(ASSOCIATED(parser),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CPPrecondition(parser%ref_count>0,cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT.failure) THEN
       !   *** Store current status ***

       old_input_line = parser%input_line
       old_input_line_number = parser%input_line_number
       old_icol = parser%icol
       old_icol1 = parser%icol1
       old_icol2 = parser%icol2

       IF (PRESENT(skip_lines)) THEN
          nline = skip_lines
       ELSE
          nline = 0
       END IF

       IF (PRESENT(newline)) THEN
          IF (newline) nline = nline + 1
       END IF

       CALL parser_get_next_line(parser,nline,at_end=at_end,error=suberror)
       CALL cp_error_propagate_error(suberror, fromWhere=routineP, error=error, &
            failure=failure, failure_level=cp_warning_level)

       IF (.NOT.failure) THEN
          IF (at_end) THEN
             test_result = "EOF"
             RETURN
          END IF

          CALL parser_next_token(parser,string_length,error=suberror)
          CALL cp_error_propagate_error(suberror, fromWhere=routineP, error=error, &
               failure=failure, failure_level=cp_warning_level)
       END IF
       IF (.NOT.failure) THEN
          IF (parser%icol1 > parser%icol2) THEN
             test_result = "EOL"
             RETURN
          END IF

          islash = parser%icol1 + INDEX(parser%input_line(parser%icol1:parser%icol2),"/") - 1

          IF (islash > parser%icol1) THEN
             nz = 2
             string(1) = parser%input_line(parser%icol1:islash-1)
             string(2) = parser%input_line(islash+1:parser%icol2)
          ELSE
             nz = 1
             string(1) = parser%input_line(parser%icol1:parser%icol2)
          END IF

          DO iz=1,nz

             IF (LEN_TRIM(string(iz)) == 0) THEN
                test_result = "STR"
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
                IF (string(iz) == parser%end_section) THEN
                   test_result = "EOS"
                ELSE
                   test_result = "STR"
                END IF
             END IF

          END DO
       END IF

       !   *** Reset to old status ***

       IF (parser%para_env%mepos==parser%para_env%source) THEN
          IF (parser%input_line_number>old_input_line_number) THEN
             DO iline=1,parser%input_line_number-old_input_line_number
                BACKSPACE (parser%input_unit)
             END DO
          END IF
       END IF


       parser%input_line = old_input_line
       parser%input_line_number = old_input_line_number
       parser%icol = old_icol
       parser%icol1 = old_icol1
       parser%icol2 = old_icol2
    END IF
    CALL cp_error_dealloc_ref(suberror)
  END FUNCTION parser_test_next_token

! *****************************************************************************

END MODULE cp_parser
