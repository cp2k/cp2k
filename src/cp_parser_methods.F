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
  USE cp_parser_buffer_types,          ONLY: copy_buffer_type
  USE cp_parser_ilist_methods,         ONLY: ilist_reset,&
                                             ilist_setup,&
                                             ilist_update
  USE cp_parser_inpp_methods,          ONLY: inpp_end_include,&
                                             inpp_expand_variables,&
                                             inpp_process_directive
  USE cp_parser_status_types,          ONLY: create_status_type,&
                                             release_status_type
  USE cp_parser_types,                 ONLY: cp_parser_type
  USE f77_blas
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp,&
                                             max_line_length
  USE message_passing,                 ONLY: mp_bcast
  USE string_utilities,                ONLY: is_whitespace,&
                                             uppercase
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: parser_test_next_token, parser_get_object, parser_location,&
            parser_search_string,parser_get_next_line, parser_skip_space,&
            parser_read_line

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
    CPPrecondition(parser%ref_count>0,cp_failure_level,routineP,error,failure)
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
!> \brief   store the present status of the parser
!> \author  Teodoro Laino [tlaino] - University of Zurich
!> \date    08.2008
! *****************************************************************************
  SUBROUTINE parser_store_status(parser, error)
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'parser_store_status', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(parser),cp_failure_level,routineP,error,failure)
    CPPrecondition(parser%ref_count>0,cp_failure_level,routineP,error,failure)
    CALL create_status_type(parser%status, error)
    parser%status%in_use                = .TRUE.
    parser%status%old_input_line        = parser%input_line
    parser%status%old_input_line_number = parser%input_line_number
    parser%status%old_icol              = parser%icol
    parser%status%old_icol1             = parser%icol1
    parser%status%old_icol2             = parser%icol2
    ! Store buffer info
    CALL copy_buffer_type(parser%buffer, parser%status%buffer, error)
  END SUBROUTINE parser_store_status

! *****************************************************************************
!> \brief   retrieve the original status of the parser
!> \author  Teodoro Laino [tlaino] - University of Zurich
!> \date    08.2008
! *****************************************************************************
  SUBROUTINE parser_retrieve_status(parser, error)
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'parser_retrieve_status', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iline
    LOGICAL                                  :: check, failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(parser),cp_failure_level,routineP,error,failure)
    CPPrecondition(parser%ref_count>0,cp_failure_level,routineP,error,failure)

    ! We are always on the same file (no possibility of switch file at this level)
    IF (parser%para_env%mepos==parser%para_env%source) THEN
       check = ALL(parser%buffer%input_line_numbers==parser%status%buffer%input_line_numbers)
       IF (.NOT.check) THEN
          ! Parser advanced.. retrieve its original location in the file
          check = ALL(parser%buffer%input_line_numbers>=parser%status%buffer%input_line_numbers)
          CPPrecondition(check,cp_failure_level,routineP,error,failure)
          DO iline=1,parser%buffer%input_line_numbers(1)-parser%status%buffer%input_line_numbers(1)
             BACKSPACE (parser%input_unit)
          END DO
       END IF
    END IF

    parser%input_line        = parser%status%old_input_line
    parser%input_line_number = parser%status%old_input_line_number
    parser%icol              = parser%status%old_icol
    parser%icol1             = parser%status%old_icol1
    parser%icol2             = parser%status%old_icol2
    CALL copy_buffer_type(parser%status%buffer, parser%buffer, error)
    CALL release_status_type(parser%status, error)

  END SUBROUTINE parser_retrieve_status

! *****************************************************************************
!> \brief   Read the next line from a logical unit "unit" (I/O node only).
!>          Skip (nline-1) lines and skip also all comment lines.
!> \author  MK
!> \date    22.11.1999
!> \version 1.0
!> \History 08.2008 [tlaino] - Teodoro Laino UZH : updated for buffer
! *****************************************************************************
  SUBROUTINE parser_read_line(parser,nline,at_end,error)
    TYPE(cp_parser_type), POINTER            :: parser
    INTEGER, INTENT(IN)                      :: nline
    LOGICAL, INTENT(out), OPTIONAL           :: at_end
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'parser_read_line', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: icol, iline, istat
    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(parser),cp_failure_level,routineP,error,failure)
    CPPrecondition(parser%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       IF (PRESENT(at_end)) at_end=.FALSE.

       DO iline = 1, nline
          ! Try to read the next line from the buffer
          CALL parser_get_line_from_buffer(parser, istat, error)

          ! Handle (persisting) read errors
          IF (istat /= 0) THEN
             IF (istat<0) THEN ! EOF/EOR is negative other errors positive
                CALL cp_assert(PRESENT(at_end),cp_failure_level,cp_assertion_failed,&
                     routineP,"Unexpected EOF "//TRIM(parser_location(parser,error=error)),&
                     error,failure)
                IF (PRESENT(at_end)) at_end=.TRUE.
                parser%icol =-1
                parser%icol1= 0
                parser%icol2=-1
             ELSE
                CALL cp_assert(.FALSE.,cp_failure_level,istat,routineP,&
                     "An input/output error occurred "//TRIM(parser_location(parser,error=error))//&
                     "(IOSTAT = "//cp_to_string(istat)//")",error,failure)
             END IF
             RETURN
          END IF
       END DO
       !     *** Reset column pointer, if a new line was read ***
       IF (nline > 0) parser%icol = 0
    END IF

  END SUBROUTINE parser_read_line

! *****************************************************************************
!> \brief   Retrieving lines from buffer
!> \author  Teodoro Laino [tlaino] - University of Zurich
!> \date    08.2008
! *****************************************************************************
  SUBROUTINE parser_get_line_from_buffer(parser, istat, error)
    TYPE(cp_parser_type), POINTER            :: parser
    INTEGER, INTENT(OUT)                     :: istat
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'parser_get_line_from_buffer', &
      routineP = moduleN//':'//routineN

    istat = 0
    ! Check buffer
    IF (parser%buffer%present_line_number==parser%buffer%size) THEN
       ! Rebuffer input file if required
       CALL parser_read_line_low(parser,error)
    END IF
    parser%buffer%present_line_number = parser%buffer%present_line_number+1
    parser%input_line_number = parser%buffer%input_line_numbers(parser%buffer%present_line_number)
    parser%input_line        = parser%buffer%input_lines(parser%buffer%present_line_number)
    IF ( (parser%buffer%istat/=0).AND.&
         (parser%buffer%last_line_number==parser%buffer%present_line_number) ) THEN
       istat = parser%buffer%istat
    END IF
  END SUBROUTINE parser_get_line_from_buffer

! *****************************************************************************
!> \brief   Low level reading subroutine with buffering
!> \author  Teodoro Laino [tlaino] - University of Zurich
!> \date    08.2008
! *****************************************************************************
  SUBROUTINE parser_read_line_low(parser,error)
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'parser_read_line_low', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=1)                         :: thischar
    INTEGER                                  :: icol, iline, imark, istat, &
                                                last_buffered_line_number
    LOGICAL                                  :: failure, non_white_found

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(parser),cp_failure_level,routineP,error,failure)
    CPPrecondition(parser%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       IF (parser%para_env%mepos==parser%para_env%source) THEN
          iline = 0
          istat = 0
          parser%buffer%present_line_number = 0
          parser%buffer%last_line_number    = parser%buffer%size
          last_buffered_line_number         = parser%buffer%input_line_numbers(parser%buffer%size)
          DO WHILE (iline /= parser%buffer%size)
             ! Increment counters by 1
             iline                     = iline + 1
             last_buffered_line_number = last_buffered_line_number  + 1

             ! Try to read the next line from file
             parser%buffer%input_line_numbers(iline) = last_buffered_line_number
             READ(UNIT=parser%input_unit,FMT="(A)",IOSTAT=istat) parser%buffer%input_lines(iline)

             ! Pre-processing steps:
             ! 1. Expand variables 2. Process directives and read next line.
             ! On read failure try to go back from included file to previous i/o-stream.
             IF (istat==0) THEN
                imark = INDEX(parser%buffer%input_lines(iline),"${")
                IF (imark/=0) THEN
                   CALL inpp_expand_variables(parser%inpp, parser%buffer%input_lines(iline),&
                        parser%input_file_name, parser%buffer%input_line_numbers(iline),error)
                END IF
                imark = INDEX(parser%buffer%input_lines(iline),"@")
                IF (imark/=0) THEN
                   CALL inpp_process_directive(parser%inpp, parser%buffer%input_lines(iline),&
                        parser%input_file_name,      parser%buffer%input_line_numbers(iline),&
                        parser%input_unit, error)
                   ! Handle index and cycle
                   last_buffered_line_number = 0
                   iline                     = iline - 1
                   CYCLE
                END IF
             ELSE IF (istat<0) THEN   ! handle EOF
                IF (parser%inpp%io_stack_level > 0) THEN
                   ! We were reading from an included file. Go back one level.
                   CALL inpp_end_include ( parser%inpp,     parser%input_file_name,&
                        parser%buffer%input_line_numbers(iline), parser%input_unit,&
                        error)
                   ! Handle index and cycle
                   last_buffered_line_number = parser%buffer%input_line_numbers(iline)
                   iline                     = iline - 1
                   CYCLE
                END IF
             END IF

             ! Saving persisting read errors
             IF (istat /= 0) THEN
                parser%buffer%istat                      = istat
                parser%buffer%last_line_number           = iline
                parser%buffer%input_line_numbers(iline:) = 0
                parser%buffer%input_lines(iline:)        = ""
                EXIT
             END IF

             ! Pre-processing and error checking done. Ready for parsing.
             non_white_found=.FALSE.
             DO icol=1,LEN(parser%buffer%input_lines(iline))
                thischar=parser%buffer%input_lines(iline)(icol:icol)
                IF (.NOT.is_whitespace(thischar)) THEN
                   IF (.NOT. is_comment(parser,thischar)) non_white_found=.TRUE.
                   EXIT
                ENDIF
             ENDDO
             IF (.NOT. non_white_found) THEN 
                iline                      = iline -1
                last_buffered_line_number  = last_buffered_line_number - 1
             END IF
          END DO
       END IF
       ! Broadcast buffer informations
       CALL broadcast_input_information(parser,error)
    END IF
  END SUBROUTINE parser_read_line_low

! *****************************************************************************
!> \brief   Broadcast the input information.
!> \author  MK
!> \date    02.03.2001
!> \History 08.2008 [tlaino] - Teodoro Laino UZH : updated for buffer
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
    CPPrecondition(parser%ref_count>0,cp_failure_level,routineP,error,failure)
    para_env => parser%para_env
    IF (.NOT.failure .AND. para_env%num_pe>1) THEN
       CALL mp_bcast(parser%buffer%present_line_number,para_env%source,para_env%group)
       CALL mp_bcast(parser%buffer%last_line_number,para_env%source,para_env%group)
       CALL mp_bcast(parser%buffer%istat,para_env%source,para_env%group)
       CALL mp_bcast(parser%buffer%input_line_numbers,para_env%source,para_env%group)
       CALL mp_bcast(parser%buffer%input_lines,para_env%source,para_env%group)
       CALL cp_error_synchronize_error(error,para_env=para_env)
    END IF

  END SUBROUTINE broadcast_input_information

! *****************************************************************************
!> \brief returns .true. if the character passed is a comment character
!> \par History
!>      02.2008 created, AK
! *****************************************************************************
  FUNCTION is_comment(parser,testchar) RESULT(resval)
    TYPE(cp_parser_type), POINTER            :: parser
    CHARACTER(LEN=1), INTENT(IN)             :: testchar
    LOGICAL                                  :: resval

    resval=.FALSE.
    ! we're in a private function, and parser has been tested before...
    IF (ANY(parser%comment_character.EQ.testchar)) THEN
       resval=.TRUE.
    END IF
  END FUNCTION is_comment

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
       IF (PRESENT(at_end)) THEN
          at_end=my_at_end
       ELSE
          IF (my_at_end) THEN
             CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                  routineP, "unexpected EOF "//&
                  TRIM(parser_location(parser,error=error)), error)
          ENDIF
       END IF
    ELSEIF (PRESENT(at_end)) THEN
       at_end=.FALSE.
    END IF
  END SUBROUTINE parser_get_next_line

! *****************************************************************************
!> \brief   Skips the whitespaces
!> \author  MK
!> \date    02.03.2001
!> \version 1.0
! *****************************************************************************
  SUBROUTINE parser_skip_space(parser,error)
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'parser_skip_space', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(parser),cp_failure_level,routineP,error,failure)
    CPPrecondition(parser%ref_count>0,cp_failure_level,routineP,error,failure)
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
               is_comment(parser,parser%input_line(parser%icol:parser%icol))) THEN
             parser%icol1 = 1
             parser%icol2 = -1
             RETURN
          END IF
  
          !       *** Check for input line continuation ***
          !       *** Ignore all white space and accept only one     ***
          !       *** separator token or a string in quotation marks ***
          IF (.NOT.is_whitespace(parser%input_line(parser%icol:parser%icol))) THEN
             
             IF (parser%input_line(parser%icol:parser%icol) == parser%continuation_character) THEN
                DO i=LEN(parser%input_line),1,-1
                   IF (.NOT.is_whitespace(parser%input_line(i:i))) EXIT
                END DO
                IF (parser%icol == i.OR. &
                     is_comment(parser,parser%input_line(i:i))) THEN
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
!> \notes   -) this function MUST be private in this module!
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
    CPPrecondition(parser%ref_count>0,cp_failure_level,routineP,error,failure)
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
                  is_comment(parser,parser%input_line(parser%icol:parser%icol))) THEN
                parser%icol1 = 1
                parser%icol2 = -1
                RETURN
             END IF
  
             !       *** Check for input line continuation ***
             IF (parser%input_line(parser%icol:parser%icol) == &
                  parser%continuation_character) THEN
                DO i=LEN(parser%input_line),1,-1
                   IF (.NOT. is_whitespace(parser%input_line(i:i))) EXIT
                END DO
                IF (parser%icol == i.OR. &
                     is_comment(parser,parser%input_line(i:i))) THEN
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
             IF (is_whitespace(parser%input_line(parser%icol:parser%icol))) THEN
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
                     is_whitespace(parser%input_line(parser%icol:parser%icol)).OR.&
                     is_comment(parser,parser%input_line(parser%icol:parser%icol)) .OR.&
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
!> \brief   Test next input object.
!>           -  test_result : "EOL": End of line
!>           -  test_result : "EOS": End of section
!>           -  test_result : "FLT": Floating point number
!>           -  test_result : "INT": Integer number
!>           -  test_result : "STR": String
!> \author  MK
!> \date    23.11.1999
!> \History 08.2008 [tlaino] - Teodoro Laino UZH : updated for buffer
! *****************************************************************************
  FUNCTION parser_test_next_token(parser,string_length,error) RESULT(test_result)
    TYPE(cp_parser_type), POINTER            :: parser
    INTEGER, INTENT(IN), OPTIONAL            :: string_length
    TYPE(cp_error_type), INTENT(inout)       :: error
    CHARACTER(LEN=3)                         :: test_result

    CHARACTER(len=*), PARAMETER :: routineN = 'parser_test_next_token', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=max_line_length)           :: string(2)
    INTEGER                                  :: idot2_first, idot_first, &
                                                idot_last, islash, istat, iz, &
                                                nz
    LOGICAL                                  :: failure, ilist_in_use
    REAL(KIND=dp)                            :: z

    failure=.FALSE.
    test_result = ""
    CPPrecondition(ASSOCIATED(parser),cp_failure_level,routineP,error,failure)
    CPPrecondition(parser%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       !   *** Store current status ***
       CALL parser_store_status(parser, error)

       ! Handle possible list of integers
       ilist_in_use = parser%ilist%in_use.AND.(parser%ilist%ipresent<parser%ilist%iend)
       IF (ilist_in_use) THEN
          test_result = "INT"
          CALL parser_retrieve_status(parser, error)
          RETURN
       END IF

       ! Otherwise continue normally
       CALL parser_next_token(parser,string_length=string_length,error=error)

       ! End of line
       IF (parser%icol1 > parser%icol2) THEN
          test_result = "EOL"
          CALL parser_retrieve_status(parser, error)
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
          idot2_first= INDEX(string(iz),"..")
          
          IF ((idot_first>0).AND.(idot2_first==0)) THEN
             idot_last = INDEX(string(iz),".",.TRUE.)
             istat = 0
             IF (idot_first == idot_last) THEN
                READ(UNIT=string(iz),FMT=*,IOSTAT=istat) z
             ELSE
                istat = 1
             END IF
             IF (istat /= 0) THEN
                test_result = "STR"
             ELSE
                test_result = "FLT"
             END IF
          ELSE IF (integer_object(string(iz),error=error)) THEN
             IF      (nz == 1)          THEN
                test_result = "INT"
             ELSE IF (idot2_first == 0) THEN
                test_result = "FLT"
             ELSE
                test_result = "STR"
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
    CALL parser_retrieve_status(parser, error)

  END FUNCTION parser_test_next_token

! *****************************************************************************
!> \brief   Search a string pattern in a file defined by its logical unit
!>          number "unit". A case sensitive search is performed, if
!>          ignore_case is .FALSE..
!>          begin_line: give back the parser at the beginning of the line
!>          matching the search
!> \author  MK
!> \date    05.10.1999
!> \History 08.2008 [tlaino] - Teodoro Laino UZH : updated for buffer
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

    found   = .FALSE.
    failure = .FALSE.
    begin   = .FALSE.
    IF (PRESENT(begin_line)) begin=begin_line

    CPPrecondition(ASSOCIATED(parser),cp_failure_level,routineP,error,failure)
    CPPrecondition(parser%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       IF (PRESENT(line)) line = ""

       !   *** Search for string pattern ***
       pattern = string
       IF (ignore_case) CALL uppercase(pattern)
       DO
          ! This call is buffered.. so should not represent any bottleneck
          CALL parser_get_next_line(parser,1,at_end=at_end,error=error)

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

       IF (found) THEN
          IF (PRESENT(line)) line = parser%input_line
          IF (.NOT.begin) CALL parser_next_token(parser,error=error)
       END IF
    END IF
  END SUBROUTINE parser_search_string

! *****************************************************************************
!> \brief   Check, if the string object contains an object of type integer.
!> \author  MK
!> \date    22.11.1999
!> \version 1.0
!> \History Introducing the possibility to parse a range of integers INT1..INT2
!>          Teodoro Laino [tlaino] - University of Zurich - 08.2008
! *****************************************************************************
  RECURSIVE FUNCTION integer_object(string,error) RESULT(contains_integer_object)

    CHARACTER(LEN=*), INTENT(IN)             :: string
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL                                  :: contains_integer_object

    INTEGER                                  :: i, j, length

    contains_integer_object = .TRUE.
    length = LEN_TRIM(string)

    IF (length == 0) THEN
      contains_integer_object = .FALSE.
      RETURN
    END IF

    IF (INDEX(string(1:length),"..") == 0) THEN
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
    ELSE
       j = INDEX(string(1:length),"..")
       contains_integer_object = integer_object(string(  1:   j-1),error).AND.&
                                 integer_object(string(j+2:length),error)
    END IF

  END FUNCTION integer_object

! *****************************************************************************
!> \brief   Read an integer number.
!> \author  MK
!> \date    22.11.1999
!> \version 1.0
! *****************************************************************************
  SUBROUTINE parser_get_integer(parser,object,lower_to_upper,newline,skip_lines,&
       string_length, at_end, error)
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
    CPPrecondition(parser%ref_count>0,cp_failure_level,routineP,error,failure)
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

       IF (parser%ilist%in_use) THEN
          CALL ilist_update(parser%ilist, error)
       END IF
       
       IF (.NOT.parser%ilist%in_use) THEN
          CALL parser_next_token(parser,string_length=string_length,error=error)

          IF (parser%icol1 > parser%icol2) THEN
             parser%icol1 = parser%icol
             parser%icol2 = parser%icol
             CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                  routineP, "An integer type object was expected, "//&
                  "found end of line"// TRIM(parser_location(parser,error=error)),&
                  error,failure)
          END IF
          ! Checks for possible lists of integers
          
          IF (INDEX(parser%input_line(parser%icol1:parser%icol2),"..")/=0) THEN
             CALL ilist_setup(parser%ilist, parser%input_line(parser%icol1:parser%icol2), error)
          END IF
       END IF

       IF (integer_object(parser%input_line(parser%icol1:parser%icol2),error=error)) THEN
          IF (parser%ilist%in_use) THEN
             object = parser%ilist%ipresent
             CALL ilist_reset(parser%ilist, error)
          ELSE
             READ(UNIT=parser%input_line(parser%icol1:parser%icol2),FMT=*) object
          END IF
       ELSE
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP, "An integer type object was expected, found <"//&
               parser%input_line(parser%icol1:parser%icol2)//">"//&
               TRIM(parser_location(parser,error=error)),error,failure)
          STOP
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
    CPPrecondition(parser%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.parser%ilist%in_use,cp_failure_level,routineP,error,failure)
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

       CALL parser_next_token(parser,string_length=string_length,error=error)

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
!> \brief   Read a floating point number.
!> \author  MK
!> \date    22.11.1999
!> \version 1.0
! *****************************************************************************
  SUBROUTINE parser_get_real(parser,object,lower_to_upper,newline,skip_lines,&
       string_length,at_end,error)
    TYPE(cp_parser_type), POINTER            :: parser
    REAL(KIND=dp), INTENT(OUT)               :: object
    LOGICAL, INTENT(IN), OPTIONAL            :: lower_to_upper, newline
    INTEGER, INTENT(IN), OPTIONAL            :: skip_lines, string_length
    LOGICAL, INTENT(out), OPTIONAL           :: at_end
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'parser_get_real', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=max_line_length)           :: string(2)
    INTEGER                                  :: islash, istat, iz, nline, nz
    LOGICAL                                  :: failure, my_at_end
    REAL(KIND=dp)                            :: z(2)

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(parser),cp_failure_level,routineP,error,failure)
    CPPrecondition(parser%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.parser%ilist%in_use,cp_failure_level,routineP,error,failure)
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

       CALL parser_next_token(parser,string_length=string_length,error=error)

       IF (parser%icol1 > parser%icol2) THEN
          parser%icol1 = parser%icol
          parser%icol2 = parser%icol
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP, "A real type object was expected, found then end of the line "//&
               TRIM(parser_location(parser,error=error)),error,failure)
       END IF

       ! Possibility to have real numbers described in the input as division between two numbers
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
                  routineP, "A real type object was expected, found white spaces! "//&
                  TRIM(parser_location(parser,error=error)),error,failure)
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
    CPPrecondition(parser%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.parser%ilist%in_use,cp_failure_level,routineP,error,failure)
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

END MODULE cp_parser_methods
