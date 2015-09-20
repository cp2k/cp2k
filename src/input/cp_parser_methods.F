!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
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
!> \author Fawzi Mohamed, Matthias Krack
! *****************************************************************************
MODULE cp_parser_methods

  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_parser_buffer_types,          ONLY: copy_buffer_type,&
                                             finalize_sub_buffer,&
                                             initialize_sub_buffer
  USE cp_parser_ilist_methods,         ONLY: ilist_reset,&
                                             ilist_setup,&
                                             ilist_update
  USE cp_parser_inpp_methods,          ONLY: inpp_end_include,&
                                             inpp_expand_variables,&
                                             inpp_process_directive
  USE cp_parser_types,                 ONLY: cp_parser_type,&
                                             parser_reset
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp,&
                                             max_line_length
  USE message_passing,                 ONLY: mp_bcast
  USE string_utilities,                ONLY: is_whitespace,&
                                             uppercase
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: parser_test_next_token, parser_get_object, parser_location,&
            parser_search_string,parser_get_next_line, parser_skip_space,&
            parser_read_line, read_float_object

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
!> \retval res ...
!> \author fawzi
! *****************************************************************************
  FUNCTION parser_location(parser) RESULT(res)

    TYPE(cp_parser_type), POINTER            :: parser
    CHARACTER(len=&
      default_path_length+default_string_length&
      )                                      :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'parser_location', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPASSERT(ASSOCIATED(parser))
    CPASSERT(parser%ref_count>0)
    res = ", File: '"//TRIM(parser%input_file_name)//"', Line: "//&
          TRIM(ADJUSTL(cp_to_string(parser%input_line_number)))//&
          ", Column: "//TRIM(ADJUSTL(cp_to_string(parser%icol)))
    IF (parser%icol == -1) THEN
       res(LEN_TRIM(res):) = " (EOF)"
    ELSE IF (MAX(1,parser%icol1) <= parser%icol2) THEN
       res(LEN_TRIM(res):) = ", Chunk: <"//&
       parser%input_line(MAX(1,parser%icol1):parser%icol2)//">"
    END IF

  END FUNCTION parser_location

! *****************************************************************************
!> \brief   store the present status of the parser
!> \param parser ...
!> \date    08.2008
!> \author  Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE parser_store_status(parser)

    TYPE(cp_parser_type), POINTER            :: parser

    CHARACTER(LEN=*), PARAMETER :: routineN = 'parser_store_status', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPASSERT(ASSOCIATED(parser))
    CPASSERT(parser%ref_count>0)
    CPASSERT(ASSOCIATED(parser%status))
    parser%status%in_use                = .TRUE.
    parser%status%old_input_line        = parser%input_line
    parser%status%old_input_line_number = parser%input_line_number
    parser%status%old_icol              = parser%icol
    parser%status%old_icol1             = parser%icol1
    parser%status%old_icol2             = parser%icol2
    ! Store buffer info
    CALL copy_buffer_type(parser%buffer, parser%status%buffer)

  END SUBROUTINE parser_store_status

! *****************************************************************************
!> \brief   retrieve the original status of the parser
!> \param parser ...
!> \date    08.2008
!> \author  Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE parser_retrieve_status(parser)

    TYPE(cp_parser_type), POINTER            :: parser

    CHARACTER(LEN=*), PARAMETER :: routineN = 'parser_retrieve_status', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPASSERT(ASSOCIATED(parser))
    CPASSERT(parser%ref_count>0)

    ! Always store the new buffer (if it is really newly read)
    IF (parser%buffer%buffer_id/=parser%status%buffer%buffer_id) THEN
       CALL initialize_sub_buffer(parser%buffer%sub_buffer, parser%buffer)
    END IF
    parser%status%in_use     = .FALSE.
    parser%input_line        = parser%status%old_input_line
    parser%input_line_number = parser%status%old_input_line_number
    parser%icol              = parser%status%old_icol
    parser%icol1             = parser%status%old_icol1
    parser%icol2             = parser%status%old_icol2

    ! Retrieve buffer info
    CALL copy_buffer_type(parser%status%buffer, parser%buffer)

  END SUBROUTINE parser_retrieve_status

! *****************************************************************************
!> \brief   Read the next line from a logical unit "unit" (I/O node only).
!>          Skip (nline-1) lines and skip also all comment lines.
!> \param parser ...
!> \param nline ...
!> \param at_end ...
!> \date    22.11.1999
!> \author  Matthias Krack (MK)
!> \version 1.0
!> \note 08.2008 [tlaino] - Teodoro Laino UZH : updated for buffer 
! *****************************************************************************
  SUBROUTINE parser_read_line(parser,nline,at_end)

    TYPE(cp_parser_type), POINTER            :: parser
    INTEGER, INTENT(IN)                      :: nline
    LOGICAL, INTENT(out), OPTIONAL           :: at_end

    CHARACTER(LEN=*), PARAMETER :: routineN = 'parser_read_line', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iline, istat
    LOGICAL                                  :: failure

    CALL timeset(routineN, handle)

    failure = .FALSE.
    CPASSERT(ASSOCIATED(parser))
    CPASSERT(parser%ref_count>0)
    IF (PRESENT(at_end)) at_end = .FALSE.

    DO iline=1,nline
       ! Try to read the next line from the buffer
       CALL parser_get_line_from_buffer(parser, istat)

       ! Handle (persisting) read errors
       IF (istat /= 0) THEN
          IF (istat < 0) THEN ! EOF/EOR is negative other errors positive
             CALL cp_assert(PRESENT(at_end),cp_failure_level,cp_assertion_failed,&
                  routineP,"Unexpected EOF"//TRIM(parser_location(parser)),&
                  failure)
             IF (PRESENT(at_end)) at_end = .TRUE.
             parser%icol  = -1
             parser%icol1 =  0
             parser%icol2 = -1
          ELSE
             CALL cp_assert(.FALSE.,cp_failure_level,istat,routineP,&
                  "An I/O error occurred (IOSTAT = "//&
                  TRIM(ADJUSTL(cp_to_string(istat)))//")"//&
                  TRIM(parser_location(parser)),failure)
          END IF
          CALL timestop(handle)
          RETURN
       END IF
    END DO

    ! Reset column pointer, if a new line was read
    IF (nline > 0) parser%icol = 0


    CALL timestop(handle)
  END SUBROUTINE parser_read_line

! *****************************************************************************
!> \brief   Retrieving lines from buffer
!> \param parser ...
!> \param istat ...
!> \date    08.2008
!> \author  Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE parser_get_line_from_buffer(parser, istat)

    TYPE(cp_parser_type), POINTER            :: parser
    INTEGER, INTENT(OUT)                     :: istat

    CHARACTER(LEN=*), PARAMETER :: routineN = 'parser_get_line_from_buffer', &
      routineP = moduleN//':'//routineN

    istat = 0
    ! Check buffer
    IF (parser%buffer%present_line_number==parser%buffer%size) THEN
       IF(ASSOCIATED(parser%buffer%sub_buffer)) THEN
          ! If the sub_buffer is initialized let's restore its buffer
          CALL finalize_sub_buffer(parser%buffer%sub_buffer,parser%buffer)
       ELSE
          ! Rebuffer input file if required
          CALL parser_read_line_low(parser)
       END IF
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
!> \param parser ...
!> \date    08.2008
!> \author  Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE parser_read_line_low(parser)

    TYPE(cp_parser_type), POINTER            :: parser

    CHARACTER(LEN=*), PARAMETER :: routineN = 'parser_read_line_low', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iline, imark, islen, &
                                                istat, &
                                                last_buffered_line_number
    LOGICAL                                  :: failure, non_white_found, &
                                                this_line_is_white_or_comment

    CALL timeset(routineN,handle)

    failure=.FALSE.
    CPASSERT(ASSOCIATED(parser))
    CPASSERT(parser%ref_count>0)
    parser%buffer%input_lines = ""
    IF (parser%para_env%mepos==parser%para_env%source) THEN
       iline = 0
       istat = 0
       parser%buffer%buffer_id           = parser%buffer%buffer_id + 1
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
          islen=LEN_TRIM(parser%buffer%input_lines(iline))

          ! Pre-processing steps:
          ! 1. Expand variables 2. Process directives and read next line.
          ! On read failure try to go back from included file to previous i/o-stream.
          IF (istat==0) THEN
             this_line_is_white_or_comment = is_comment_line(parser, parser%buffer%input_lines(iline))
             IF (.NOT.this_line_is_white_or_comment .AND. parser%apply_preprocessing) THEN
                imark = INDEX(parser%buffer%input_lines(iline)(1:islen),"$")
                IF (imark/=0) THEN
                   CALL inpp_expand_variables(parser%inpp, parser%buffer%input_lines(iline),&
                        parser%input_file_name, parser%buffer%input_line_numbers(iline))
                   islen=LEN_TRIM(parser%buffer%input_lines(iline))
                END IF
                imark = INDEX(parser%buffer%input_lines(iline)(1:islen),"@")
                IF (imark/=0) THEN
                   CALL inpp_process_directive(parser%inpp, parser%buffer%input_lines(iline),&
                        parser%input_file_name,      parser%buffer%input_line_numbers(iline),&
                        parser%input_unit)
                   islen=LEN_TRIM(parser%buffer%input_lines(iline))
                   ! Handle index and cycle
                   last_buffered_line_number = 0
                   iline                     = iline - 1
                   CYCLE
                END IF
             END IF
          ELSE IF (istat<0) THEN   ! handle EOF
             IF (parser%inpp%io_stack_level > 0) THEN
                ! We were reading from an included file. Go back one level.
                CALL inpp_end_include ( parser%inpp,     parser%input_file_name,&
                     parser%buffer%input_line_numbers(iline), parser%input_unit)
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
          IF (.NOT.parser%parse_white_lines) THEN
             non_white_found=.NOT.this_line_is_white_or_comment
          ELSE
             non_white_found=.TRUE.
          END IF
          IF (.NOT. non_white_found) THEN
             iline                      = iline -1
             last_buffered_line_number  = last_buffered_line_number - 1
          END IF
       END DO
    END IF
    ! Broadcast buffer informations
    CALL broadcast_input_information(parser)

    CALL timestop(handle)

  END SUBROUTINE parser_read_line_low

! *****************************************************************************
!> \brief   Broadcast the input information.
!> \param parser ...
!> \date    02.03.2001
!> \author  Matthias Krack (MK)
!> \note 08.2008 [tlaino] - Teodoro Laino UZH : updated for buffer 
! *****************************************************************************
  SUBROUTINE broadcast_input_information(parser)

    TYPE(cp_parser_type), POINTER            :: parser

    CHARACTER(len=*), PARAMETER :: routineN = 'broadcast_input_information', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    TYPE(cp_para_env_type), POINTER          :: para_env

    CALL timeset(routineN,handle)

    failure=.FALSE.
    CPASSERT(ASSOCIATED(parser))
    CPASSERT(parser%ref_count>0)
    para_env => parser%para_env
    IF (.NOT.failure .AND. para_env%num_pe>1) THEN
       CALL mp_bcast(parser%buffer%buffer_id,para_env%source,para_env%group)
       CALL mp_bcast(parser%buffer%present_line_number,para_env%source,para_env%group)
       CALL mp_bcast(parser%buffer%last_line_number,para_env%source,para_env%group)
       CALL mp_bcast(parser%buffer%istat,para_env%source,para_env%group)
       CALL mp_bcast(parser%buffer%input_line_numbers,para_env%source,para_env%group)
       CALL mp_bcast(parser%buffer%input_lines,para_env%source,para_env%group)
    END IF

    CALL timestop(handle)

  END SUBROUTINE broadcast_input_information

! *****************************************************************************
!> \brief returns .true. if the line is a comment line or an empty line
!> \param parser ...
!> \param line ...
!> \retval resval ...
!> \par History
!>      03.2009 [tlaino] - Teodoro Laino
! *****************************************************************************
  FUNCTION is_comment_line(parser, line) RESULT(resval)

    TYPE(cp_parser_type), POINTER            :: parser
    CHARACTER(LEN=*), INTENT(IN)             :: line
    LOGICAL                                  :: resval

    CHARACTER(len=*), PARAMETER :: routineN = 'is_comment_line', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=1)                         :: thischar
    INTEGER                                  :: icol

    resval = .TRUE.
    DO icol=1,LEN(line)
       thischar=line(icol:icol)
       IF (.NOT.is_whitespace(thischar)) THEN
          IF (.NOT.is_comment(parser,thischar)) resval=.FALSE.
          EXIT
       END IF
    END DO

  END FUNCTION is_comment_line

! *****************************************************************************
!> \brief returns .true. if the character passed is a comment character
!> \param parser ...
!> \param testchar ...
!> \retval resval ...
!> \par History
!>      02.2008 created, AK
!> \author AK
! *****************************************************************************
  FUNCTION is_comment(parser,testchar) RESULT(resval)

    TYPE(cp_parser_type), POINTER            :: parser
    CHARACTER(LEN=1), INTENT(IN)             :: testchar
    LOGICAL                                  :: resval

    CHARACTER(len=*), PARAMETER :: routineN = 'is_comment', &
      routineP = moduleN//':'//routineN

    resval = .FALSE.
    ! We are in a private function, and parser has been tested before...
    IF (ANY(parser%comment_character.EQ.testchar)) THEN
       resval = .TRUE.
    END IF

  END FUNCTION is_comment

! *****************************************************************************
!> \brief   Read the next input line and broadcast the input information.
!>          Skip (nline-1) lines and skip also all comment lines.
!> \param parser ...
!> \param nline ...
!> \param at_end ...
!> \date    22.11.1999
!> \author  Matthias Krack (MK)
!> \version 1.0
! *****************************************************************************
  SUBROUTINE parser_get_next_line(parser,nline,at_end)

    TYPE(cp_parser_type), POINTER            :: parser
    INTEGER, INTENT(IN)                      :: nline
    LOGICAL, INTENT(out), OPTIONAL           :: at_end

    CHARACTER(len=*), PARAMETER :: routineN = 'parser_get_next_line', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: my_at_end

    IF (nline > 0) THEN
       CALL parser_read_line(parser,nline,at_end=my_at_end)
       IF (PRESENT(at_end)) THEN
          at_end = my_at_end
       ELSE
          IF (my_at_end) THEN
             CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                            "Unexpected EOF"//TRIM(parser_location(parser)))
          END IF
       END IF
    ELSE IF (PRESENT(at_end)) THEN
       at_end = .FALSE.
    END IF

  END SUBROUTINE parser_get_next_line

! *****************************************************************************
!> \brief   Skips the whitespaces
!> \param parser ...
!> \date    02.03.2001
!> \author  Matthias Krack (MK)
!> \version 1.0
! *****************************************************************************
  SUBROUTINE parser_skip_space(parser)
    TYPE(cp_parser_type), POINTER            :: parser

    CHARACTER(LEN=*), PARAMETER :: routineN = 'parser_skip_space', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPASSERT(ASSOCIATED(parser))
    CPASSERT(parser%ref_count>0)
    ! Variable input string length (automatic search)

    ! Check for EOF
    IF (parser%icol == -1) THEN
       parser%icol1 = 1
       parser%icol2 = -1
       RETURN
    END IF

    ! Search for the beginning of the next input string
    DO
       ! Increment the column counter
       parser%icol = parser%icol + 1

       ! Quick return, if the end of line is found
       IF ((parser%icol > LEN_TRIM(parser%input_line)).OR.&
            is_comment(parser,parser%input_line(parser%icol:parser%icol))) THEN
          parser%icol1 = 1
          parser%icol2 = -1
          RETURN
       END IF

       ! Check for input line continuation
       ! Ignore all white space and accept only one
       ! separator token or a string in quotation marks
       IF (.NOT.is_whitespace(parser%input_line(parser%icol:parser%icol))) THEN

          IF (parser%input_line(parser%icol:parser%icol) == parser%continuation_character) THEN
             DO i=LEN(parser%input_line),1,-1
                IF (.NOT.is_whitespace(parser%input_line(i:i))) EXIT
             END DO
             IF (parser%icol == i.OR. &
                  is_comment(parser,parser%input_line(i:i))) THEN
                CALL parser_get_next_line(parser,1)
                ! does not pass at_end, it is an error to continue to EOF, change this behavoiur?
                CYCLE
             ELSE
                parser%icol1 = i
                parser%icol2 = i
                CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                     "Found non-blank tokens after the line continuation character '"//&
                     parser%continuation_character//"' "//&
                     TRIM(parser_location(parser)),failure)
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

  END SUBROUTINE parser_skip_space

! *****************************************************************************
!> \brief   Get the next input string from the input line.
!> \param parser ...
!> \param string_length ...
!> \date    19.02.2001
!> \author  Matthias Krack (MK)
!> \version 1.0
!> \notes   -) this function MUST be private in this module!
! *****************************************************************************
  SUBROUTINE parser_next_token(parser,string_length)

    TYPE(cp_parser_type), POINTER            :: parser
    INTEGER, INTENT(IN), OPTIONAL            :: string_length

    CHARACTER(LEN=*), PARAMETER :: routineN = 'parser_next_token', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, len_inputline, &
                                                len_trim_inputline, length
    LOGICAL                                  :: failure

    failure=.FALSE.
    CPASSERT(ASSOCIATED(parser))
    CPASSERT(parser%ref_count>0)
    IF (PRESENT(string_length)) THEN
       length = MIN(string_length,max_line_length)
    ELSE
       length = 0
    END IF

    IF (length > 0) THEN
       ! Fixed input string length
       CALL cp_assert(parser%icol/=-1,cp_failure_level,cp_assertion_failed,routineP,&
            "Unexpectetly reached EOF"//TRIM(parser_location(parser)),&
            failure)

       length=MIN(LEN_TRIM(parser%input_line)-parser%icol1+1,length)
       parser%icol1 = parser%icol + 1
       parser%icol2 = parser%icol + length
       parser%icol = parser%icol2
    ELSE

       ! Variable input string length (automatic search)
       ! Check for EOF
       IF (parser%icol == -1) THEN
          parser%icol1 =  1
          parser%icol2 = -1
          RETURN
       END IF

       ! Search for the beginning of the next input string
       ! we precompute those, they are expensive to get
       len_inputline =      LEN(parser%input_line)
       len_trim_inputline = LEN_TRIM(parser%input_line)

       DO

          ! Increment the column counter
          parser%icol = parser%icol + 1

          ! Quick return, if the end of line is found
          IF ((parser%icol > len_trim_inputline).OR.&
               is_comment(parser,parser%input_line(parser%icol:parser%icol))) THEN
             parser%icol1 = 1
             parser%icol2 = -1
             RETURN
          END IF

          ! Check for input line continuation
          IF (parser%input_line(parser%icol:parser%icol) == &
               parser%continuation_character) THEN
             DO i=LEN(parser%input_line),1,-1
                IF (.NOT. is_whitespace(parser%input_line(i:i))) EXIT
             END DO
             IF (parser%icol == i.OR. &
                  is_comment(parser,parser%input_line(i:i))) THEN
                CALL parser_get_next_line(parser,1)
                ! does not pass at_end, it is an error to continue to EOF, change this behavoiur?
                CYCLE
             ELSE
                parser%icol1 = i
                parser%icol2 = i
                CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                     "Found non-blank tokens after the line continuation character '"//&
                     parser%continuation_character//"' "//TRIM(parser_location(parser)),&
                     failure)
                EXIT
             END IF
          END IF
          ! Ignore all white space and accept only one
          ! separator token or a string in quotation marks
          IF (is_whitespace(parser%input_line(parser%icol:parser%icol))) THEN
             CYCLE
          ELSE IF (INDEX(parser%separators,parser%input_line(parser%icol:parser%icol)) > 0) THEN
             IF (parser%first_separator) THEN
                parser%first_separator = .FALSE.
                CYCLE
             ELSE
                parser%icol1 = parser%icol
                parser%icol2 = parser%icol
                CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                     "Unexpected separator token <"//parser%input_line(parser%icol:parser%icol)//&
                     "> found"//TRIM(parser_location(parser)),failure)
                EXIT
             END IF
          ELSE IF (parser%input_line(parser%icol:parser%icol) == '"') THEN
             ! identify quotes and preserve them through-out the execution..
             parser%first_separator = .TRUE.
             parser%icol1 = parser%icol + 1
             parser%icol2 = parser%icol + INDEX(parser%input_line(parser%icol1:),'"')
             IF (parser%icol2 == parser%icol) THEN
                parser%icol1 = parser%icol
                parser%icol2 = parser%icol
                CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                     "Unmatched quotation mark found"//TRIM(parser_location(parser)),&
                     failure)
             ELSE
                parser%icol1 = parser%icol1 - 1
                parser%icol  = parser%icol2
                RETURN
             END IF
          ELSE
             parser%first_separator = .TRUE.
             parser%icol1 = parser%icol
             EXIT
          END IF

       END DO

       ! Search for the end of the next input string
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

  END SUBROUTINE parser_next_token

! *****************************************************************************
!> \brief   Test next input object.
!>           -  test_result : "EOL": End of line
!>           -  test_result : "EOS": End of section
!>           -  test_result : "FLT": Floating point number
!>           -  test_result : "INT": Integer number
!>           -  test_result : "STR": String
!> \param parser ...
!> \param string_length ...
!> \retval test_result ...
!> \date    23.11.1999
!> \author  Matthias Krack (MK)
!> \note - 08.2008 [tlaino] - Teodoro Laino UZH : updated for buffer 
!>          - Major rewrite to parse also (multiple) products of integer or 
!>            floating point numbers (23.11.2012,MK)
! *****************************************************************************
  FUNCTION parser_test_next_token(parser,string_length) RESULT(test_result)

    TYPE(cp_parser_type), POINTER            :: parser
    INTEGER, INTENT(IN), OPTIONAL            :: string_length
    CHARACTER(LEN=3)                         :: test_result

    CHARACTER(len=*), PARAMETER :: routineN = 'parser_test_next_token', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=max_line_length)           :: error_message, string
    INTEGER                                  :: iz, n
    LOGICAL                                  :: failure, ilist_in_use
    REAL(KIND=dp)                            :: fz

    failure = .FALSE.
    test_result = ""
    CPASSERT(ASSOCIATED(parser))
    CPASSERT(parser%ref_count>0)


    ! Store current status
    CALL parser_store_status(parser)

    ! Handle possible list of integers
    ilist_in_use = parser%ilist%in_use.AND.(parser%ilist%ipresent < parser%ilist%iend)
    IF (ilist_in_use) THEN
       test_result = "INT"
       CALL parser_retrieve_status(parser)
       RETURN
    END IF

    ! Otherwise continue normally
    CALL parser_next_token(parser,string_length=string_length)

    ! End of line
    IF (parser%icol1 > parser%icol2) THEN
       test_result = "EOL"
       CALL parser_retrieve_status(parser)
       RETURN
    END IF

    string = parser%input_line(parser%icol1:parser%icol2)
    n = LEN_TRIM(string)

    IF (n == 0) THEN
       test_result = "STR"
       CALL parser_retrieve_status(parser)
       RETURN
    END IF

    ! Check for end section string
    IF (string(1:n) == parser%end_section) THEN
       test_result = "EOS"
       CALL parser_retrieve_status(parser)
       RETURN
    END IF

    ! Check for integer object
    error_message = ""
    CALL read_integer_object(string(1:n),iz,error_message)
    IF (LEN_TRIM(error_message) == 0) THEN
       test_result = "INT"
       CALL parser_retrieve_status(parser)
       RETURN
    END IF

    ! Check for floating point object
    error_message = ""
    CALL read_float_object(string(1:n),fz,error_message)
    IF (LEN_TRIM(error_message) == 0) THEN
       test_result = "FLT"
       CALL parser_retrieve_status(parser)
       RETURN
    END IF

    test_result = "STR"
    CALL parser_retrieve_status(parser)


  END FUNCTION parser_test_next_token

! *****************************************************************************
!> \brief   Search a string pattern in a file defined by its logical unit
!>          number "unit". A case sensitive search is performed, if
!>          ignore_case is .FALSE..
!>          begin_line: give back the parser at the beginning of the line
!>          matching the search
!> \param parser ...
!> \param string ...
!> \param ignore_case ...
!> \param found ...
!> \param line ...
!> \param begin_line ...
!> \param search_from_begin_of_file ...
!> \date    05.10.1999
!> \author  MK
!> \note 08.2008 [tlaino] - Teodoro Laino UZH : updated for buffer 
! *****************************************************************************
  SUBROUTINE parser_search_string(parser,string,ignore_case,found,line,begin_line,&
             search_from_begin_of_file)

    TYPE(cp_parser_type), POINTER            :: parser
    CHARACTER(LEN=*), INTENT(IN)             :: string
    LOGICAL, INTENT(IN)                      :: ignore_case
    LOGICAL, INTENT(OUT)                     :: found
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: line
    LOGICAL, INTENT(IN), OPTIONAL            :: begin_line, &
                                                search_from_begin_of_file

    CHARACTER(LEN=*), PARAMETER :: routineN = 'parser_search_string', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=LEN(string))               :: pattern
    CHARACTER(LEN=max_line_length+1)         :: current_line
    INTEGER                                  :: ipattern
    LOGICAL                                  :: at_end, begin, do_reset, &
                                                failure

    found   = .FALSE.
    failure = .FALSE.
    begin   = .FALSE.
    do_reset= .FALSE.
    IF (PRESENT(begin_line)) begin=begin_line
    IF (PRESENT(search_from_begin_of_file)) do_reset = search_from_begin_of_file
    CPASSERT(ASSOCIATED(parser))
    CPASSERT(parser%ref_count>0)
    IF (PRESENT(line)) line = ""

    ! Search for string pattern
    pattern = string
    IF (ignore_case) CALL uppercase(pattern)
    IF (do_reset)    CALL parser_reset(parser)
    DO
       ! This call is buffered.. so should not represent any bottleneck
       CALL parser_get_next_line(parser,1,at_end=at_end)

       ! Exit loop, if the end of file is reached
       IF (at_end) EXIT

       ! Check the current line for string pattern
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
                  TRIM(ADJUSTL(cp_to_string(LEN(line))))//&
                  " characters and is therefore too long to fit in the "//&
                  "specified variable"//&
                  TRIM(parser_location(parser)),failure)
          END IF
          EXIT
       END IF

    END DO

    IF (found) THEN
       IF (begin) parser%icol = 0
    END IF

    IF (found) THEN
       IF (PRESENT(line)) line = parser%input_line
       IF (.NOT.begin) CALL parser_next_token(parser)
    END IF

  END SUBROUTINE parser_search_string

! *****************************************************************************
!> \brief   Check, if the string object contains an object of type integer.
!> \param string ...
!> \retval contains_integer_object ...
!> \date    22.11.1999
!> \author  Matthias Krack (MK)
!> \version 1.0
!> \note - Introducing the possibility to parse a range of integers INT1..INT2 
!>            Teodoro Laino [tlaino] - University of Zurich - 08.2008 
!>          - Parse also a product of integer numbers (23.11.2012,MK) 
! *****************************************************************************
  FUNCTION integer_object(string) RESULT(contains_integer_object)

    CHARACTER(LEN=*), INTENT(IN)             :: string
    LOGICAL                                  :: contains_integer_object

    INTEGER                                  :: i, idots, istar, n

    contains_integer_object = .TRUE.
    n = LEN_TRIM(string)

    IF (n == 0) THEN
      contains_integer_object = .FALSE.
      RETURN
    END IF

    idots  = INDEX(string(1:n),"..")
    istar  = INDEX(string(1:n),"*")

    IF (idots /= 0) THEN
       contains_integer_object = is_integer(string(1:idots-1)).AND.&
                                 is_integer(string(idots+2:n))
    ELSE IF (istar /= 0) THEN
       i = 1
       DO WHILE (istar /= 0)
          IF (.NOT.is_integer(string(i:i+istar-2))) THEN
             contains_integer_object = .FALSE.
             RETURN
          END IF
          i = i + istar
          istar = INDEX(string(i:n),"*")
       END DO
       contains_integer_object = is_integer(string(i:n))
    ELSE
       contains_integer_object = is_integer(string(1:n))
    END IF

  END FUNCTION integer_object

! *****************************************************************************
!> \brief ...
!> \param string ...
!> \retval check ...
! *****************************************************************************
  FUNCTION is_integer(string) RESULT(check)

    CHARACTER(LEN=*), INTENT(IN)             :: string
    LOGICAL                                  :: check

    INTEGER                                  :: i, n

    check = .TRUE.
    n = LEN_TRIM(string)

    IF (n == 0) THEN
      check = .FALSE.
      RETURN
    END IF

    IF ((INDEX("+-",string(1:1)) > 0).AND.(n == 1)) THEN
       check = .FALSE.
       RETURN
    END IF

    IF (INDEX("+-0123456789",string(1:1)) == 0) THEN
       check = .FALSE.
       RETURN
    END IF

    DO i=2,n
       IF (INDEX("0123456789",string(i:i)) == 0) THEN
          check = .FALSE.
          RETURN
       END IF
    END DO

  END FUNCTION is_integer

! *****************************************************************************
!> \brief   Read an integer number.
!> \param parser ...
!> \param object ...
!> \param newline ...
!> \param skip_lines ...
!> \param string_length ...
!> \param at_end ...
!> \date    22.11.1999
!> \author  Matthias Krack (MK)
!> \version 1.0
! *****************************************************************************
  SUBROUTINE parser_get_integer(parser,object,newline,skip_lines,&
                                string_length, at_end)

    TYPE(cp_parser_type), POINTER            :: parser
    INTEGER, INTENT(OUT)                     :: object
    LOGICAL, INTENT(IN), OPTIONAL            :: newline
    INTEGER, INTENT(IN), OPTIONAL            :: skip_lines, string_length
    LOGICAL, INTENT(out), OPTIONAL           :: at_end

    CHARACTER(LEN=*), PARAMETER :: routineN = 'parser_get_integer', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=max_line_length)           :: error_message
    INTEGER                                  :: nline
    LOGICAL                                  :: failure, my_at_end

    failure = .FALSE.
    CPASSERT(ASSOCIATED(parser))
    CPASSERT(parser%ref_count>0)
    IF (PRESENT(skip_lines)) THEN
       nline = skip_lines
    ELSE
       nline = 0
    END IF

    IF (PRESENT(newline)) THEN
       IF (newline) nline = nline + 1
    END IF

    CALL parser_get_next_line(parser,nline,at_end=my_at_end)
    IF (PRESENT(at_end)) THEN
       at_end = my_at_end
       IF (my_at_end) RETURN
    ELSE IF (my_at_end) THEN
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
            "Unexpected EOF"//TRIM(parser_location(parser)),&
            failure)
    END IF

    IF (parser%ilist%in_use) THEN
       CALL ilist_update(parser%ilist)
    ELSE
       CALL parser_next_token(parser,string_length=string_length)
       IF (parser%icol1 > parser%icol2) THEN
          parser%icol1 = parser%icol
          parser%icol2 = parser%icol
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
               "An integer type object was expected, found end of line"//&
               TRIM(parser_location(parser)),failure)
       END IF
       ! Checks for possible lists of integers
       IF (INDEX(parser%input_line(parser%icol1:parser%icol2),"..") /= 0) THEN
          CALL ilist_setup(parser%ilist, parser%input_line(parser%icol1:parser%icol2))
       END IF
    END IF

    IF (integer_object(parser%input_line(parser%icol1:parser%icol2))) THEN
       IF (parser%ilist%in_use) THEN
          object = parser%ilist%ipresent
          CALL ilist_reset(parser%ilist)
       ELSE
          CALL read_integer_object(parser%input_line(parser%icol1:parser%icol2),object,error_message)
          IF (LEN_TRIM(error_message) > 0) THEN
             CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                  TRIM(error_message)//TRIM(parser_location(parser)),&
                  failure)
          END IF
       END IF
    ELSE
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
            "An integer type object was expected, found <"//&
            parser%input_line(parser%icol1:parser%icol2)//">"//&
            TRIM(parser_location(parser)),failure)
    END IF


  END SUBROUTINE parser_get_integer

! *****************************************************************************
!> \brief   Read a string representing logical object.
!> \param parser ...
!> \param object ...
!> \param newline ...
!> \param skip_lines ...
!> \param string_length ...
!> \param at_end ...
!> \date    01.04.2003
!> \par History
!>      - New version (08.07.2003,MK)
!> \author  FM
!> \version 1.0
! *****************************************************************************
  SUBROUTINE parser_get_logical(parser,object,newline,skip_lines,&
                                string_length,at_end)

    TYPE(cp_parser_type), POINTER            :: parser
    LOGICAL, INTENT(OUT)                     :: object
    LOGICAL, INTENT(IN), OPTIONAL            :: newline
    INTEGER, INTENT(IN), OPTIONAL            :: skip_lines, string_length
    LOGICAL, INTENT(out), OPTIONAL           :: at_end

    CHARACTER(LEN=*), PARAMETER :: routineN = 'parser_get_logical', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=max_line_length)           :: input_string
    INTEGER                                  :: input_string_length, nline
    LOGICAL                                  :: failure, my_at_end

    failure=.FALSE.
    CPASSERT(ASSOCIATED(parser))
    CPASSERT(parser%ref_count>0)
    CPASSERT(.NOT.parser%ilist%in_use)
    IF (PRESENT(skip_lines)) THEN
       nline = skip_lines
    ELSE
       nline = 0
    END IF

    IF (PRESENT(newline)) THEN
       IF (newline) nline = nline + 1
    END IF

    CALL parser_get_next_line(parser,nline,at_end=my_at_end)
    IF (PRESENT(at_end)) THEN
       at_end=my_at_end
       IF (my_at_end) RETURN
    ELSE IF (my_at_end) THEN
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
            "Unexpected EOF"//TRIM(parser_location(parser)),&
            failure)
    END IF

    CALL parser_next_token(parser,string_length=string_length)

    input_string_length = parser%icol2 - parser%icol1 + 1

    IF (input_string_length == 0) THEN
       parser%icol1 = parser%icol
       parser%icol2 = parser%icol
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
            "A string representing a logical object was expected, found end of line"//&
            TRIM(parser_location(parser)),failure)
    ELSE
       input_string = ""
       input_string(:input_string_length) = parser%input_line(parser%icol1:parser%icol2)
    END IF
    CALL uppercase(input_string)

    SELECT CASE (TRIM(input_string))
    CASE ("0","F",".F.","FALSE",".FALSE.","N","NO","OFF")
       object = .FALSE.
    CASE ("1","T",".T.","TRUE",".TRUE.","Y","YES","ON")
       object = .TRUE.
    CASE DEFAULT
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
            "A string representing a logical object was expected, found <"//&
            TRIM(input_string)//">"//TRIM(parser_location(parser)),&
            failure)
    END SELECT

  END SUBROUTINE parser_get_logical

! *****************************************************************************
!> \brief   Read a floating point number.
!> \param parser ...
!> \param object ...
!> \param newline ...
!> \param skip_lines ...
!> \param string_length ...
!> \param at_end ...
!> \date    22.11.1999
!> \author  Matthias Krack (MK)
!> \version 1.0
! *****************************************************************************
  SUBROUTINE parser_get_real(parser,object,newline,skip_lines,string_length,&
                             at_end)

    TYPE(cp_parser_type), POINTER            :: parser
    REAL(KIND=dp), INTENT(OUT)               :: object
    LOGICAL, INTENT(IN), OPTIONAL            :: newline
    INTEGER, INTENT(IN), OPTIONAL            :: skip_lines, string_length
    LOGICAL, INTENT(out), OPTIONAL           :: at_end

    CHARACTER(LEN=*), PARAMETER :: routineN = 'parser_get_real', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=max_line_length)           :: error_message
    INTEGER                                  :: nline
    LOGICAL                                  :: failure, my_at_end

    failure=.FALSE.
    CPASSERT(ASSOCIATED(parser))
    CPASSERT(parser%ref_count>0)
    CPASSERT(.NOT.parser%ilist%in_use)
    

    IF (PRESENT(skip_lines)) THEN
       nline = skip_lines
    ELSE
       nline = 0
    END IF

    IF (PRESENT(newline)) THEN
       IF (newline) nline = nline + 1
    END IF

    CALL parser_get_next_line(parser,nline,at_end=my_at_end)
    IF (PRESENT(at_end)) THEN
       at_end=my_at_end
       IF (my_at_end) RETURN
    ELSE IF (my_at_end) THEN
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
            "Unexpected EOF"//TRIM(parser_location(parser)),&
            failure)
    END IF

    CALL parser_next_token(parser,string_length=string_length)

    IF (parser%icol1 > parser%icol2) THEN
       parser%icol1 = parser%icol
       parser%icol2 = parser%icol
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
            "A floating point type object was expected, found end of the line"//&
            TRIM(parser_location(parser)),failure)
    END IF

    ! Possibility to have real numbers described in the input as division between two numbers
    CALL read_float_object(parser%input_line(parser%icol1:parser%icol2),object,error_message)
    IF (LEN_TRIM(error_message) > 0) THEN
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
            TRIM(error_message)//TRIM(parser_location(parser)),&
            failure)
    END IF


  END SUBROUTINE parser_get_real

! *****************************************************************************
!> \brief   Read a string.
!> \param parser ...
!> \param object ...
!> \param lower_to_upper ...
!> \param newline ...
!> \param skip_lines ...
!> \param string_length ...
!> \param at_end ...
!> \date    22.11.1999
!> \author  Matthias Krack (MK)
!> \version 1.0
! *****************************************************************************
  SUBROUTINE parser_get_string(parser,object,lower_to_upper,newline,skip_lines,&
                                string_length,at_end)

    TYPE(cp_parser_type), POINTER            :: parser
    CHARACTER(LEN=*), INTENT(OUT)            :: object
    LOGICAL, INTENT(IN), OPTIONAL            :: lower_to_upper, newline
    INTEGER, INTENT(IN), OPTIONAL            :: skip_lines, string_length
    LOGICAL, INTENT(out), OPTIONAL           :: at_end

    CHARACTER(LEN=*), PARAMETER :: routineN = 'parser_get_string', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: input_string_length, nline
    LOGICAL                                  :: failure, my_at_end

    object = ""
    failure=.FALSE.
    CPASSERT(ASSOCIATED(parser))
    CPASSERT(parser%ref_count>0)
    CPASSERT(.NOT.parser%ilist%in_use)
    IF (PRESENT(skip_lines)) THEN
       nline = skip_lines
    ELSE
       nline = 0
    END IF

    IF (PRESENT(newline)) THEN
       IF (newline) nline = nline + 1
    END IF

    CALL parser_get_next_line(parser,nline,at_end=my_at_end)
    IF (PRESENT(at_end)) THEN
       at_end=my_at_end
       IF (my_at_end) RETURN
    ELSE IF (my_at_end) THEN
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
            "Unexpected EOF"//TRIM(parser_location(parser)),&
            failure)
    END IF

    CALL parser_next_token(parser,string_length)

    input_string_length = parser%icol2 - parser%icol1 + 1

    IF (input_string_length <= 0) THEN
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
            "A string type object was expected, found end of line"//&
            TRIM(parser_location(parser)),failure)
    ELSE IF (input_string_length > LEN(object)) THEN
       CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
            "The input string <"//parser%input_line(parser%icol1:parser%icol2)//&
            "> has more than "//cp_to_string(LEN(object))//&
            " characters and is therefore too long to fit in the "//&
            "specified variable"//TRIM(parser_location(parser)),&
            failure)
       object = parser%input_line(parser%icol1:parser%icol1+LEN(object)-1)
    ELSE
       object(:input_string_length) = parser%input_line(parser%icol1:parser%icol2)
    END IF

    ! Convert lowercase to uppercase, if requested
    IF (PRESENT(lower_to_upper)) THEN
       IF (lower_to_upper) CALL uppercase(object)
    END IF

  END SUBROUTINE parser_get_string

! *****************************************************************************
!> \brief   Returns a floating point number read from a string including
!>          fraction like z1/z2.
!> \param string ...
!> \param object ...
!> \param error_message ...
!> \date    11.01.2011 (MK)
!> \author  Matthias Krack
!> \version 1.0
!> \note - Parse also multiple products and fractions of floating point
!>            numbers (23.11.2012,MK)
! *****************************************************************************
  SUBROUTINE read_float_object(string,object,error_message)

    CHARACTER(LEN=*), INTENT(IN)             :: string
    REAL(KIND=dp), INTENT(OUT)               :: object
    CHARACTER(LEN=*), INTENT(OUT)            :: error_message

    INTEGER                                  :: i, iop, islash, istar, istat, &
                                                n
    LOGICAL                                  :: parsing_done
    REAL(KIND=dp)                            :: z

    error_message = ""

    i = 1
    iop = 0
    n = LEN_TRIM(string)

    parsing_done = .FALSE.

    DO WHILE (.NOT.parsing_done)
       i = i + iop
       islash = INDEX(string(i:n),"/")
       istar  = INDEX(string(i:n),"*")
       IF ((islash == 0).AND.(istar == 0)) THEN
          ! Last factor found: read it and then exit the loop
          iop = n - i + 2
          parsing_done = .TRUE.
       ELSE IF ((islash > 0).AND.(istar > 0)) THEN
          iop = MIN(islash,istar)
       ELSE IF (islash > 0) THEN
          iop = islash
       ELSE IF (istar > 0) THEN
          iop = istar
       END IF
       READ (UNIT=string(i:i+iop-2),FMT=*,IOSTAT=istat) z
       IF (istat /= 0) THEN
          error_message = "A floating point type object was expected, found <"//&
                          string(i:i+iop-2)//">"
          RETURN
       END IF
       IF (i == 1) THEN
          object = z
       ELSE IF (string(i-1:i-1) == "*") THEN
          object = object*z
       ELSE
          IF (z == 0.0_dp) THEN
             error_message = "Division by zero found <"//&
                             string(i:i+iop-2)//">"
             RETURN
          ELSE
             object = object/z
          END IF
       END IF
    END DO

  END SUBROUTINE read_float_object

! *****************************************************************************
!> \brief   Returns an integer number read from a string including products of
!>          integer numbers like iz1*iz2*iz3
!> \param string ...
!> \param object ...
!> \param error_message ...
!> \date    23.11.2012 (MK)
!> \author  Matthias Krack
!> \version 1.0
!> \note - Parse also (multiple) products of integer numbers (23.11.2012,MK)
! *****************************************************************************
  SUBROUTINE read_integer_object(string,object,error_message)

    CHARACTER(LEN=*), INTENT(IN)             :: string
    INTEGER, INTENT(OUT)                     :: object
    CHARACTER(LEN=*), INTENT(OUT)            :: error_message

    CHARACTER(LEN=20)                        :: fmtstr
    INTEGER                                  :: i, iop, istat, iz, n
    LOGICAL                                  :: parsing_done

    error_message = ""

    i = 1
    iop = 0
    n = LEN_TRIM(string)

    parsing_done = .FALSE.

    DO WHILE (.NOT.parsing_done)
       i = i + iop
       iop = INDEX(string(i:n),"*")
       IF (iop == 0) THEN
          ! Last factor found: read it and then exit the loop
          iop = n - i + 2
          parsing_done= .TRUE.
       END IF
       IF (iop-1>0) THEN
          ! need an explicit fmtstr here. With 'FMT=*' compilers from intel and pgi will also
          ! read float numbers as integers, without setting istat non-zero, i.e. string="0.3", istat=0, iz=0
          ! this leads to wrong CP2K results (e.g. parsing force fields).
          WRITE(fmtstr,FMT='(A,I0,A)') '(I',iop-1,')'
          READ (UNIT=string(i:i+iop-2),FMT=fmtstr,IOSTAT=istat) iz
       ELSE
          istat=1
       ENDIF
       IF (istat /= 0) THEN
          error_message = "An integer type object was expected, found <"//&
                          string(i:i+iop-2)//">"
          RETURN
       END IF
       IF (i == 1) THEN
          object = iz
       ELSE
          object = object*iz
       END IF
    END DO

  END SUBROUTINE read_integer_object

END MODULE cp_parser_methods
