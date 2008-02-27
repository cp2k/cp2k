!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Utility routines to read data from files.
!>      Kept as close as possible to the old parser because
!>        1. string handling is a weak point of fortran compilers, and it is
!>           easy to write correct things that do not work
!>        2. conversion of old code
!> \par History
!>      22.11.1999 first version of the old parser (called qs_parser)
!>                 Matthias Krack
!>      06.2004 removed module variables, cp_parser_type, new module [fawzi]
!> \author fawzi
! *****************************************************************************
MODULE cp_parser_methods
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_parser_types,                 ONLY: blank_character,&
                                             cp_parser_type,&
                                             horizontal_tab,&
                                             max_line_length
  USE f77_blas
  USE inpp_methods,                    ONLY: inpp_end_include,&
                                             inpp_expand_variables,&
                                             inpp_process_directive
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE message_passing,                 ONLY: mp_bcast
  USE string_utilities,                ONLY: uppercase
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: parser_next_token, parser_test_next_token
  PUBLIC :: parser_get_object, parser_get_real, parser_get_integer,&
       parser_get_logical, parser_get_string, parser_location,&
       parser_search_string,parser_get_next_line, parser_skip_space

  ! *** Global parameters ***
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_parser_methods'

  INTERFACE parser_get_object
    MODULE PROCEDURE parser_get_integer,&
                     parser_get_logical,&
                     parser_get_real,&
                     parser_get_string
  END INTERFACE

CONTAINS

! *****************************************************************************
!> \brief   Broadcast the input information.
!> \author  MK
!> \date    02.03.2001
!> \version 1.0
! *****************************************************************************
  SUBROUTINE broadcast_input_information(parser,error)
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'broadcast_input_information', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(cp_para_env_type), POINTER          :: para_env

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
!> \brief   Read the next input line and broadcast the input information.
!>          Skip (nline-1) lines and skip also all comment lines.
!> \author  MK
!> \date    22.11.1999
!> \version 1.0
! *****************************************************************************
  SUBROUTINE parser_get_next_line(parser,nline,at_end,error)
    TYPE(cp_parser_type), POINTER            :: parser
    INTEGER, INTENT(IN)                      :: nline
    LOGICAL, INTENT(out), OPTIONAL           :: at_end
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'parser_get_next_line', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: my_at_end

    IF (nline>0) THEN
       CALL parser_read_line(parser,nline,at_end=my_at_end,error=error)
       CALL mp_bcast(my_at_end, parser%para_env%source,&
            parser%para_env%group)
       IF (PRESENT(at_end)) THEN
          at_end=my_at_end
       ELSE
          IF (my_at_end) THEN
             CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                            routineP, "unexpected EOF "//&
                            TRIM(parser_location(parser,error=error)), error)
          ENDIF
       END IF
       CALL broadcast_input_information(parser,error=error)
    ELSEIF (PRESENT(at_end)) THEN
       at_end=.FALSE.
    END IF

  END SUBROUTINE parser_get_next_line

! *****************************************************************************
!> \brief return a description of the part of the file acually parsed
!> \param parser the parser
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  FUNCTION parser_location(parser,error) RESULT(res)
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(cp_error_type), INTENT(inout)       :: error
    CHARACTER(len=&
      default_path_length+default_string_length&
      )                                      :: res

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

! *****************************************************************************
!> \brief   skips the whitespaces
!> \author  MK
!> \date    02.03.2001
!> \version 1.0
! *****************************************************************************
  SUBROUTINE parser_skip_space(parser,multiline,error)
    TYPE(cp_parser_type), POINTER            :: parser
    LOGICAL, INTENT(in), OPTIONAL            :: multiline
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'parser_skip_space', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure, my_continue

    failure=.FALSE.
    my_continue=.TRUE.
    IF (PRESENT(multiline)) my_continue=multiline
    CPPrecondition(ASSOCIATED(parser),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CPPrecondition(parser%ref_count>0,cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT.failure) THEN
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
               ANY(parser%input_line(parser%icol:parser%icol) == parser%comment_character)) THEN
             parser%icol1 = 1
             parser%icol2 = -1
             RETURN
          END IF
  
          !       *** Check for input line continuation ***
  
          !       *** Ignore all white space and accept only one     ***
          !       *** separator token or a string in quotation marks ***
  
          IF (.NOT.((parser%input_line(parser%icol:parser%icol) == blank_character).OR.&
               (parser%input_line(parser%icol:parser%icol) == horizontal_tab))) THEN
  
             IF (parser%input_line(parser%icol:parser%icol) == &
                  parser%continuation_character.and.my_continue) THEN
                DO i=LEN(parser%input_line),1,-1
                   IF ((parser%input_line(i:i) /= blank_character).AND.&
                        (parser%input_line(i:i) /= horizontal_tab)) EXIT
                END DO
                IF (parser%icol == i.OR. ANY(parser%input_line(i:i)==&
                     parser%comment_character)) THEN
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
             ELSE
                parser%icol=parser%icol-1
                parser%icol1=parser%icol
                parser%icol2=parser%icol
                EXIT
             END IF
          END IF
       END DO
    END IF
  END SUBROUTINE parser_skip_space

! *****************************************************************************
!> \brief   Get the next input string from the input line.
!> \author  MK
!> \date    19.02.2001
!> \version 1.0
! *****************************************************************************
  SUBROUTINE parser_next_token(parser,string_length,error)
    TYPE(cp_parser_type), POINTER            :: parser
    INTEGER, INTENT(IN), OPTIONAL            :: string_length
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'parser_next_token', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, len_inputline, &
                                                len_trim_inputline, length
    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(parser),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CPPrecondition(parser%ref_count>0,cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT.failure) THEN
       IF (PRESENT(string_length)) THEN
          length = MIN(string_length,max_line_length)
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
  
          ! we precompute those, they are expensive to get
          len_inputline =      LEN(parser%input_line)
          len_trim_inputline = LEN_TRIM(parser%input_line)
  
          DO
  
             !       *** Increment the column counter ***
  
             parser%icol = parser%icol + 1
  
             !       *** Quick return, if the end of line is found ***
  
             IF ((parser%icol > len_trim_inputline).OR.&
                  ANY(parser%input_line(parser%icol:parser%icol) == parser%comment_character)) THEN
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
                IF (parser%icol == i.or.ANY(parser%input_line(i:i)==&
                     parser%comment_character)) THEN
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
                IF (parser%icol > len_trim_inputline .or.parser%icol==len_inputline) EXIT
                parser%icol = parser%icol + 1
                IF ((parser%icol > len_trim_inputline).OR.&
                     (parser%input_line(parser%icol:parser%icol) == blank_character).OR.&
                     (parser%input_line(parser%icol:parser%icol) == horizontal_tab).OR.&
                     ANY(parser%input_line(parser%icol:parser%icol) == parser%comment_character).OR.&
                     (parser%input_line(parser%icol:parser%icol) == parser%continuation_character)) THEN
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
!> \brief   Check, if the string object contains an object of type integer.
!> \author  MK
!> \date    22.11.1999
!> \version 1.0
! *****************************************************************************
  FUNCTION integer_object(string,error) RESULT(contains_integer_object)

    CHARACTER(LEN=*), INTENT(IN)             :: string
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL                                  :: contains_integer_object

    INTEGER                                  :: i, length

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
!> \brief   Read an integer number.
!> \author  MK
!> \date    22.11.1999
!> \version 1.0
! *****************************************************************************
  SUBROUTINE parser_get_integer(parser,object,lower_to_upper,newline,skip_lines,&
       string_length, at_end,error)
    TYPE(cp_parser_type), POINTER            :: parser
    INTEGER, INTENT(OUT)                     :: object
    LOGICAL, INTENT(IN), OPTIONAL            :: lower_to_upper, newline
    INTEGER, INTENT(IN), OPTIONAL            :: skip_lines, string_length
    LOGICAL, INTENT(out), OPTIONAL           :: at_end
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'parser_get_integer', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nline
    LOGICAL                                  :: failure, my_at_end

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
       IF (integer_object(parser%input_line(parser%icol1:parser%icol2),error=error)) THEN
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
!> \brief   Read a string representing logical object.
!> \author  FM
!> \date    01.04.2003
!> \par History
!>      - New version (08.07.2003,MK)
!> \version 1.0
! *****************************************************************************
  SUBROUTINE parser_get_logical(parser,object,newline,skip_lines,string_length,&
       at_end,error)
    TYPE(cp_parser_type), POINTER            :: parser
    LOGICAL, INTENT(OUT)                     :: object
    LOGICAL, INTENT(IN), OPTIONAL            :: newline
    INTEGER, INTENT(IN), OPTIONAL            :: skip_lines, string_length
    LOGICAL, INTENT(out), OPTIONAL           :: at_end
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'parser_get_logical', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=max_line_length)           :: input_string
    INTEGER                                  :: input_string_length, nline
    LOGICAL                                  :: failure, my_at_end

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
!> \brief   Read the next line from a logical unit "unit" (I/O node only).
!>          Skip (nline-1) lines and skip also all comment lines.
!> \author  MK
!> \date    22.11.1999
!> \version 1.0
! *****************************************************************************
  SUBROUTINE parser_read_line(parser,nline,at_end,error)
    TYPE(cp_parser_type), POINTER            :: parser
    INTEGER, INTENT(IN)                      :: nline
    LOGICAL, INTENT(out), OPTIONAL           :: at_end
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'parser_read_line', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=1)                         :: thischar
    INTEGER                                  :: icol, iline, imark, istat
    LOGICAL                                  :: failure, non_white_found

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
             ! Try to read the next line from file
             parser%input_line_number = parser%input_line_number + 1
             READ (UNIT=parser%input_unit,FMT="(A)",IOSTAT=istat) parser%input_line

             ! Pre-processing steps:
             ! 1. Expand variables 2. Process directives and read next line.
             ! On read failure try to go back from included file to previous i/o-stream.
             IF (istat==0) THEN
                imark = INDEX(parser%input_line,"${")
                IF (imark/=0) THEN
                   CALL inpp_expand_variables(parser%inpp, parser, error)
                END IF
                imark = INDEX(parser%input_line,"@")
                IF (imark/=0) THEN
                   CALL inpp_process_directive(parser%inpp, parser, error)
                   CYCLE
                END IF
             ELSE IF (istat<0) THEN   ! handle EOF
                IF (parser%inpp%io_stack_level > 0) THEN
                   ! We were reading from an included file. Go back one level.
                   CALL inpp_end_include(parser%inpp, parser, error)
                   CYCLE
                END IF
             END IF

             ! Handle (persisting) read errors
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

             ! Pre-processing and error checking done. Ready for parsing.
             iline = iline + 1
             non_white_found=.FALSE.
             DO icol=1,LEN(parser%input_line)
                thischar=parser%input_line(icol:icol)
                IF ( (thischar.NE.blank_character).AND.&
                     (thischar.NE.horizontal_tab)) THEN
                   IF (.NOT. ANY(parser%comment_character.EQ.thischar)) non_white_found=.TRUE.
                   EXIT
                ENDIF
             ENDDO
             IF (.NOT. non_white_found)  iline = iline -1

          END DO

          !     *** Reset column pointer, if a new line was read ***

          IF (nline > 0) parser%icol = 0

       END IF
    END IF

  END SUBROUTINE parser_read_line

! *****************************************************************************
!> \brief   Read a floating point number.
!> \author  MK
!> \date    22.11.1999
!> \version 1.0
! *****************************************************************************
  SUBROUTINE parser_get_real(parser,object,lower_to_upper,newline,skip_lines,&
       string_length,multiline,at_end,error)
    TYPE(cp_parser_type), POINTER            :: parser
    REAL(KIND=dp), INTENT(OUT)               :: object
    LOGICAL, INTENT(IN), OPTIONAL            :: lower_to_upper, newline
    INTEGER, INTENT(IN), OPTIONAL            :: skip_lines, string_length
    LOGICAL, INTENT(IN), OPTIONAL            :: multiline
    LOGICAL, INTENT(out), OPTIONAL           :: at_end
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'parser_get_real', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=max_line_length), &
      DIMENSION(2)                           :: string
    INTEGER                                  :: islash, istat, iz, nline, nz
    LOGICAL                                  :: failure, my_at_end, &
                                                my_multiline
    REAL(KIND=dp), DIMENSION(2)              :: z

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

          READ(UNIT=string(iz),FMT=*,IOSTAT=istat) z(iz)

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
!> \brief   Read a string.
!> \author  MK
!> \date    22.11.1999
!> \version 1.0
! *****************************************************************************
  SUBROUTINE parser_get_string(parser,object,lower_to_upper,newline,skip_lines,&
                                string_length,at_end,error)
    TYPE(cp_parser_type), POINTER            :: parser
    CHARACTER(LEN=*), INTENT(OUT)            :: object
    LOGICAL, INTENT(IN), OPTIONAL            :: lower_to_upper, newline
    INTEGER, INTENT(IN), OPTIONAL            :: skip_lines, string_length
    LOGICAL, INTENT(out), OPTIONAL           :: at_end
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'parser_get_string', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: input_string_length, nline
    LOGICAL                                  :: failure, my_at_end

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
!> \brief   Search a string pattern in a file defined by its logical unit
!>          number "unit". A case sensitive search is performed, if
!>          ignore_case is .FALSE..
!>          begin_line: give back the parser at the beginning of the line
!>          matching the search
!> \author  MK
!> \date    05.10.1999
!> \version 1.0
! *****************************************************************************
  SUBROUTINE parser_search_string(parser,string,ignore_case,found,line,begin_line,error)
    TYPE(cp_parser_type), POINTER            :: parser
    CHARACTER(LEN=*), INTENT(IN)             :: string
    LOGICAL, INTENT(IN)                      :: ignore_case
    LOGICAL, INTENT(OUT)                     :: found
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: line
    LOGICAL, INTENT(IN), OPTIONAL            :: begin_line
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'parser_search_string', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=LEN(string))               :: pattern
    CHARACTER(LEN=max_line_length+1)         :: current_line
    INTEGER                                  :: ipattern
    LOGICAL                                  :: at_end, begin, failure

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
!> \brief   Test next input object.
!>           -  test_result : "EOF": End of file
!>           -  test_result : "EOL": End of line
!>           -  test_result : "EOS": End of section
!>           -  test_result : "FLT": Floating point number
!>           -  test_result : "INT": Integer number
!>           -  test_result : "STR": String
!> \author  MK
!> \date    23.11.1999
!> \version 1.0
! *****************************************************************************
  FUNCTION parser_test_next_token(parser,newline,skip_lines,string_length,error)&
       RESULT(test_result)
    TYPE(cp_parser_type), POINTER            :: parser
    LOGICAL, INTENT(IN), OPTIONAL            :: newline
    INTEGER, INTENT(IN), OPTIONAL            :: skip_lines, string_length
    TYPE(cp_error_type), INTENT(inout)       :: error
    CHARACTER(LEN=3)                         :: test_result

    CHARACTER(len=*), PARAMETER :: routineN = 'parser_test_next_token', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=max_line_length)           :: old_input_line
    CHARACTER(LEN=max_line_length), &
      DIMENSION(2)                           :: string
    INTEGER :: idot_first, idot_last, iline, islash, istat, iz, nline, nz, &
      old_icol, old_icol1, old_icol2, old_input_line_number
    LOGICAL                                  :: at_end, failure
    REAL(KIND=dp)                            :: z
    TYPE(cp_error_type)                      :: suberror

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
             CALL cp_error_dealloc_ref(suberror)
             RETURN
          END IF

          CALL parser_next_token(parser,string_length,error=suberror)
          CALL cp_error_propagate_error(suberror, fromWhere=routineP, error=error, &
               failure=failure, failure_level=cp_warning_level)
       END IF
       IF (.NOT.failure) THEN
          IF (parser%icol1 > parser%icol2) THEN
             test_result = "EOL"
             CALL cp_error_dealloc_ref(suberror)
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
             ELSE IF (integer_object(string(iz),error=error)) THEN
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

END MODULE cp_parser_methods
