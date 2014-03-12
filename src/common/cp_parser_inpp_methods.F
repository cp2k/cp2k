!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief a module to allow simple internal preprocessing in input files.
!> \par History
!>      - standalone proof-of-concept implementation (20.02.2008,AK)
!>      - integration into cp2k (22.02.2008,tlaino)
!>      - variables added (23.02.2008,AK)
!>      - @IF/@ENDIF added (25.02.2008,AK)
!>      - @PRINT and debug ifdefs added (26.02.2008,AK)
!> \author Axel Kohlmeyer [AK] - CMM/UPenn Philadelphia
!> \date 20.02.2008
! *****************************************************************************
MODULE cp_parser_inpp_methods
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_parser_inpp_types,            ONLY: inpp_type
  USE f77_blas
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length
  USE machine,                         ONLY: default_output_unit
  USE memory_utilities,                ONLY: reallocate
  USE string_utilities,                ONLY: is_whitespace,&
                                             uppercase
  USE termination,                     ONLY: stop_program
#include "../cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_parser_inpp_methods'
  LOGICAL, PARAMETER, PRIVATE          :: debug_this_module=.FALSE.
  INTEGER, PARAMETER, PRIVATE          :: max_message_length = 400

  PUBLIC  :: inpp_process_directive, inpp_end_include, inpp_expand_variables
  PRIVATE :: inpp_find_variable, inpp_list_variables

CONTAINS

! *****************************************************************************
!> \brief process internal preprocessor directives like @INCLUDE, @SET, @IF/@ENDIF
!> \par History
!>      - standalone proof-of-concept implementation (20.02.2008,AK)
!>      - integration into cp2k (22.02.2008,tlaino)
!>      - variables added (23.02.2008,AK)
!>      - @IF/@ENDIF added (25.02.2008,AK)
!> \author AK
! *****************************************************************************
  SUBROUTINE inpp_process_directive(inpp, input_line, input_file_name, input_line_number,&
       input_unit, error)
    TYPE(inpp_type), POINTER                 :: inpp
    CHARACTER(LEN=*), INTENT(INOUT)          :: input_line, input_file_name
    INTEGER, INTENT(INOUT)                   :: input_line_number, input_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'inpp_process_directive', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: cond1, cond2, filename, &
                                                mytag, value, varname
    CHARACTER(LEN=max_message_length)        :: message
    INTEGER                                  :: i, indf, indi, istat, pos1, &
                                                pos2, unit
    LOGICAL                                  :: check, failure

    failure = .TRUE.
    CPPostcondition(ASSOCIATED(inpp),cp_failure_level,routineP,error,failure)

    ! find location of directive in line and check whether it is commented out
    indi = INDEX(input_line,"@")
    pos1 = INDEX(input_line,"!")
    pos2 = INDEX(input_line,"#")
    IF(((pos1>0).AND.(pos1<indi)) .OR. ((pos2>0).AND.(pos2<indi)) ) THEN
       ! nothing to do here.
       RETURN
    ENDIF

    ! Get the start of the instruction and find "@KEYWORD" (or "@")
    indf = indi
    DO WHILE (.NOT. is_whitespace(input_line(indf:indf)))
       indf = indf + 1
    END DO
    mytag = input_line(indi:indf-1)
    CALL uppercase(mytag)

    SELECT CASE(mytag)

    CASE ("@INCLUDE")
       ! Get the filename.. allow for " or ' or nothing..
       filename = TRIM(input_line(indf:))
       IF(LEN_TRIM(filename)==0) THEN
          WRITE(UNIT=message,FMT="(3A,I6)") &
               "INPP_@INCLUDE: Incorrect @INCLUDE directive in file: ",&
               TRIM(input_file_name), "  Line:",input_line_number
          CALL stop_program(routineN,moduleN,__LINE__,TRIM(message))
       ENDIF
       indi=1
       DO WHILE (is_whitespace(filename(indi:indi)))
          indi = indi + 1
       END DO
       filename=TRIM(filename(indi:))

       ! handle quoting of the filename
       pos1 = INDEX(filename,'"')
       pos2 = INDEX(filename(pos1+1:),'"')
       IF((pos1/=0).AND.(pos2/=0)) THEN
          filename = filename(pos1+1:pos1+pos2-1)
       ELSE
          pos1 = INDEX(filename,"'")
          pos2 = INDEX(filename(pos1+1:),"'")
          IF((pos1/=0).AND.(pos2/=0)) THEN
             filename = filename(pos1+1:pos1+pos2-1)
          ELSE
             ! incorrect quotes (only one of ' or ").
             pos2 = INDEX(filename,'"')
             IF((pos1/=0).OR.(pos2/=0)) THEN
                WRITE(UNIT=message,FMT="(3A,I6)") &
                     "INPP_@INCLUDE: Incorrect quoting of include file in file: ",&
                     TRIM(input_file_name), "  Line:",input_line_number
                CALL stop_program(routineN,moduleN,__LINE__,TRIM(message))
             ENDIF
             ! nothing to do. unquoted filename.
          ENDIF
       ENDIF

       ! Let's check that files already opened won't be again opened
       DO i = 1, inpp%io_stack_level
          check = TRIM(filename)/=TRIM(inpp%io_stack_filename(i))
          CPPostcondition(check,cp_failure_level,routineP,error,failure)
       END DO

       ! this stops on error, so we can always assume success
       CALL open_file(file_name=TRIM(filename),&
                      file_status="OLD",&
                      file_form="FORMATTED",&
                      file_action="READ",&
                      unit_number=unit)

       IF (debug_this_module) THEN
          WRITE(UNIT=message,FMT="(3A,I6,2A)") "INPP_@INCLUDE: in file: ",&
               TRIM(input_file_name),"  Line:",input_line_number,&
               " Opened include file: ",TRIM(filename)
          WRITE(default_output_unit,*) TRIM(message)
       END IF

       ! make room, save status and position the parser at the beginning of new file.
       inpp%io_stack_level = inpp%io_stack_level+1
       CALL reallocate(inpp%io_stack_channel,1,inpp%io_stack_level)
       CALL reallocate(inpp%io_stack_lineno,1,inpp%io_stack_level)
       CALL reallocate(p_long=inpp%io_stack_filename,lb_new=1,ub_new=inpp%io_stack_level)

       inpp%io_stack_channel(inpp%io_stack_level)  = input_unit
       inpp%io_stack_lineno(inpp%io_stack_level)   = input_line_number
       inpp%io_stack_filename(inpp%io_stack_level) = input_file_name

       input_file_name = TRIM(filename)
       input_line_number = 0
       input_unit = unit

    CASE ("@SET")
       ! split directive into variable name and value data.
       varname = TRIM(input_line(indf:))
       IF(LEN_TRIM(varname)==0) THEN
          WRITE(UNIT=message,FMT="(3A,I6)") &
               "INPP_@SET: Incorrect @SET directive in file: ",&
               TRIM(input_file_name), "  Line:",input_line_number
          CALL stop_program(routineN,moduleN,__LINE__,TRIM(message))
       ENDIF

       indi=1
       DO WHILE (is_whitespace(varname(indi:indi)))
          indi = indi + 1
       END DO
       indf=indi
       DO WHILE (.NOT.is_whitespace(varname(indf:indf)))
          indf = indf + 1
       END DO
       value=TRIM(varname(indf:))
       varname=TRIM(varname(indi:indf-1))
       indi=1
       DO WHILE (is_whitespace(value(indi:indi)))
          indi = indi + 1
       END DO
       value=TRIM(value(indi:))

       IF(LEN_TRIM(value)==0) THEN
          WRITE(UNIT=message,FMT="(3A,I6)") &
               "INPP_@SET: Incorrect @SET directive in file: ",&
               TRIM(input_file_name), "  Line:",input_line_number
          CALL stop_program(routineN,moduleN,__LINE__,TRIM(message))
       ENDIF

       ! sort into table of variables.
       indi=inpp_find_variable(inpp,varname)
       IF(indi == 0) THEN
          ! create new variable
          inpp%num_variables=inpp%num_variables+1
          CALL reallocate(p_long=inpp%variable_name,lb_new=1,ub_new=inpp%num_variables)
          CALL reallocate(p_long=inpp%variable_value,lb_new=1,ub_new=inpp%num_variables)
          inpp%variable_name(inpp%num_variables)=varname
          inpp%variable_value(inpp%num_variables)=value
          IF (debug_this_module) THEN
             WRITE(UNIT=message,FMT="(3A,I6,4A)") "INPP_@SET: in file: ",&
                  TRIM(input_file_name),"  Line:",input_line_number,&
                  " Set new variable ", TRIM(varname), " to value: ",TRIM(value)
             WRITE(default_output_unit,*) TRIM(message)
          END IF
       ELSE
          ! reassign variable
          IF (debug_this_module) THEN
             WRITE(UNIT=message,FMT="(3A,I6,6A)") "INPP_@SET: in file: ",&
                  TRIM(input_file_name),"  Line:",input_line_number,&
                  " Change variable ", TRIM(varname), " from value: ",&
                  TRIM(inpp%variable_value(indi)), " to value: ",TRIM(value)
             WRITE(default_output_unit,*) TRIM(message)
          END IF
          inpp%variable_value(indi)=value
       ENDIF

       IF (debug_this_module) CALL inpp_list_variables(inpp,6)

    CASE ("@IF")
       ! detect IF expression.
       ! we recognize lexical equality or inequality, and presence of
       ! a string (true) vs. blank (false). in case the expression resolves
       ! to "false" we read lines here until we reach an @ENDIF or EOF.
       indi=indf
       pos1=INDEX(input_line,"==")
       pos2=INDEX(input_line,"/=")
       ! shave off leading whitespace
       DO WHILE (is_whitespace(input_line(indi:indi)))
          indi = indi + 1
          IF (indi > LEN_TRIM(input_line)) EXIT
       END DO
       check=.FALSE.
       IF(pos1>0) THEN
          cond1=input_line(indi:pos1-1)
          cond2=input_line(pos1+2:)
          check=.TRUE.
          IF((pos2>0).OR.(INDEX(cond2,"==")>0)) THEN
             WRITE(UNIT=message,FMT="(3A,I6)") &
                  "INPP_@IF: Incorrect @IF directive in file: ",&
                  TRIM(input_file_name), "  Line:",input_line_number
             CALL stop_program(routineN,moduleN,__LINE__,TRIM(message))
          ENDIF
       ELSEIF(pos2>0) THEN
          cond1=input_line(indi:pos2-1)
          cond2=input_line(pos2+2:)
          check=.FALSE.
          IF((pos1>0).OR.(INDEX(cond2,"/=")>0)) THEN
             WRITE(UNIT=message,FMT="(3A,I6)") &
                  "INPP_@IF: Incorrect @IF directive in file: ",&
                  TRIM(input_file_name), "  Line:",input_line_number
             CALL stop_program(routineN,moduleN,__LINE__,TRIM(message))
          ENDIF
       ELSE
          IF(LEN_TRIM(input_line(indi:))>0) THEN
             IF(TRIM(input_line(indi:))=='0') THEN
                cond1='XXX'
                cond2='XXX'
                check=.FALSE.
             ELSE
                cond1='XXX'
                cond2='XXX'
                check=.TRUE.
             ENDIF
          ELSE
             cond1='XXX'
             cond2='XXX'
             check=.FALSE.
          ENDIF
       ENDIF

       ! Get rid of possible parentheses
       IF (INDEX(cond1,"(") /= 0) cond1 = cond1(INDEX(cond1,"(")+1:)
       IF (INDEX(cond2,")") /= 0) cond2 = cond2(1:INDEX(cond2,")")-1)

       ! Shave off leading whitespace from cond1
       indi=1
       DO WHILE (is_whitespace(cond1(indi:indi)))
          indi = indi + 1
       END DO
       cond1=cond1(indi:)

       ! Shave off leading whitespace from cond2
       indi=1
       DO WHILE (is_whitespace(cond2(indi:indi)))
          indi = indi + 1
       END DO
       cond2=cond2(indi:)

       IF(LEN_TRIM(cond2)==0) THEN
          WRITE(UNIT=message,FMT="(3A,I6)") &
               "INPP_@IF: Incorrect @IF directive in file: ",&
               TRIM(input_file_name), "  Line:",input_line_number
          CALL stop_program(routineN,moduleN,__LINE__,TRIM(message))
       ENDIF

       IF ((TRIM(cond1)==TRIM(cond2)) .EQV. check) THEN
          IF (debug_this_module) THEN
             WRITE(UNIT=message,FMT="(3A,I6,A)") "INPP_@IF: in file: ",&
                  TRIM(input_file_name),"  Line:",input_line_number,&
                  " Conditional ("//TRIM(cond1)//","//TRIM(cond2)//&
                  ") resolves to true. Continuing parsing."
             WRITE(default_output_unit,*) TRIM(message)
          END IF
          ! resolves to true. keep on reading normally...
          RETURN
       ELSE
          IF (debug_this_module) THEN
             WRITE(UNIT=message,FMT="(3A,I6,A)") "INPP_@IF: in file: ",&
                  TRIM(input_file_name),"  Line:",input_line_number,&
                  " Conditional ("//TRIM(cond1)//","//TRIM(cond2)//&
                  ") resolves to false. Skipping Lines."
             WRITE(default_output_unit,*) TRIM(message)
          END IF
          istat=0
          DO WHILE(istat==0)
             input_line_number = input_line_number + 1
             READ (UNIT=input_unit,FMT="(A)",IOSTAT=istat) input_line
             IF (debug_this_module) THEN
                WRITE(UNIT=message,FMT="(1A,I6,2A)") "INPP_@IF: skipping line ",&
                  input_line_number, ": ",TRIM(input_line)
                WRITE(default_output_unit,*) TRIM(message)
             END IF

             indi = INDEX(input_line,"@")
             pos1 = INDEX(input_line,"!")
             pos2 = INDEX(input_line,"#")
             IF(((pos1>0).AND.(pos1<indi)) .OR. ((pos2>0).AND.(pos2<indi)) ) THEN
                ! comment. nothing to do here.
                CYCLE
             ENDIF

             ! Get the start of the instruction and find "@KEYWORD"
             indi = MAX(1,indi)
             indf = indi
             DO WHILE (input_line(indf:indf)/=" ")
                indf = indf + 1
             END DO
             CPPostcondition((indf-indi)<=default_string_length,cp_failure_level,routineP,error,failure)
             mytag =input_line(indi:indf-1)
             CALL uppercase(mytag)
             IF(INDEX(mytag,"@ENDIF")>0) THEN
                ! ok found it. go back to normal
                IF (debug_this_module) WRITE(default_output_unit,*) "INPP_@IF: found @ENDIF. End of skipping."
                RETURN
             ENDIF
          END DO
          IF(istat/=0) THEN
             WRITE(UNIT=message,FMT="(3A,I6)") &
                  "INPP_@IF: Error while looking for @ENDIF directive in file: ",&
                  TRIM(input_file_name), "  Line:",input_line_number
             CALL stop_program(routineN,moduleN,__LINE__,TRIM(message))
          ENDIF
       ENDIF

    CASE ("@ENDIF")
       IF (debug_this_module) WRITE(default_output_unit,*) "INPP_@IF: found @ENDIF in normal parsing. Ignoring it."
       ! nothing to do. just return to skip the line.
       RETURN

    CASE ("@PRINT")
       ! for debugging of variables etc.
       WRITE(UNIT=message,FMT="(3A,I6,2A)") "INPP_@PRINT: in file: ",&
            TRIM(input_file_name),"  Line:",input_line_number,&
            " Text: ",TRIM(input_line(indf:))
       WRITE(default_output_unit,*) TRIM(message)
       RETURN
       ! Do Nothing..
    END SELECT

  END SUBROUTINE inpp_process_directive


! *****************************************************************************
!> \brief Restore older file status from stack after EOF on include file.
!> \par History
!>      - standalone proof-of-concept implemenation (20.02.2008,AK)
!>      - integrated into cp2k (21.02.2008)
!> \author AK
! *****************************************************************************
  SUBROUTINE inpp_end_include(inpp, input_file_name, input_line_number, input_unit, error)
    TYPE(inpp_type), POINTER                 :: inpp
    CHARACTER(LEN=*), INTENT(INOUT)          :: input_file_name
    INTEGER, INTENT(INOUT)                   :: input_line_number, input_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'inpp_end_include', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(inpp),cp_failure_level,routineP,error,failure)
    IF(inpp%io_stack_level > 0) THEN
       CALL close_file(input_unit)
       input_unit=inpp%io_stack_channel(inpp%io_stack_level)
       input_line_number=inpp%io_stack_lineno(inpp%io_stack_level)
       input_file_name=TRIM(inpp%io_stack_filename(inpp%io_stack_level))
       inpp%io_stack_level=inpp%io_stack_level-1
       CALL reallocate(inpp%io_stack_channel,1,inpp%io_stack_level)
       CALL reallocate(inpp%io_stack_lineno,1,inpp%io_stack_level)
       CALL reallocate(p_long=inpp%io_stack_filename,lb_new=1,ub_new=inpp%io_stack_level)
    ENDIF

  END SUBROUTINE inpp_end_include

! *****************************************************************************
!> \brief expand all ${VAR} or $VAR variable entries on the input string
!> \par History
!>      - standalone proof-of-concept implemenation (22.02.2008,AK)
!>      - integrated into cp2k (23.02.2008)
!> \author AK
! *****************************************************************************
  SUBROUTINE inpp_expand_variables(inpp, input_line, input_file_name, input_line_number,&
       error)
    TYPE(inpp_type), POINTER                 :: inpp
    CHARACTER(LEN=*), INTENT(INOUT)          :: input_line, input_file_name
    INTEGER, INTENT(IN)                      :: input_line_number
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'inpp_expand_variables', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: newline, varname
    CHARACTER(LEN=max_message_length)        :: message
    INTEGER                                  :: idx, pos1, pos2
    LOGICAL                                  :: failure

    failure = .TRUE.
    CPPostcondition(ASSOCIATED(inpp),cp_failure_level,routineP,error,failure)

    ! process line until all variables named with the convention ${VAR} are expanded
    DO WHILE (INDEX(input_line,'${') > 0)
       pos1=INDEX(input_line,'${')
       pos1=pos1+2
       pos2=INDEX(input_line(pos1:),'}')
       CPPostcondition(pos2>0,cp_failure_level,routineP,error,failure)
       IF(pos2==0) THEN
          WRITE(UNIT=message,FMT="(3A,I6)") &
               "INPP_@SET: Missing '}' in file: ",&
               TRIM(input_file_name), "  Line:",input_line_number
          CALL stop_program(routineN,moduleN,__LINE__,TRIM(message))
       ENDIF

       pos2=pos1+pos2-2
       varname=input_line(pos1:pos2)
       idx=inpp_find_variable(inpp,varname)
       CPPostcondition(idx>0,cp_failure_level,routineP,error,failure)
       IF(idx==0) THEN
          WRITE(UNIT=message,FMT="(5A,I6)") &
               "INPP_@SET: Variable ${",TRIM(varname),"} not defined in file: ",&
               TRIM(input_file_name), "  Line:",input_line_number
          CALL stop_program(routineN,moduleN,__LINE__,TRIM(message))
       ENDIF

       IF(pos1>3) THEN
          newline=input_line(1:pos1-3) // TRIM(inpp%variable_value(idx)) &
               // input_line(pos2+2:)
       ELSE
          newline=TRIM(inpp%variable_value(idx)) // input_line(pos2+2:)
       ENDIF

       input_line=newline
    END DO

    ! process line until all variables named with the convention $VAR are expanded
    DO WHILE (INDEX(input_line,'$') > 0)
       pos1=INDEX(input_line,'$')
       pos1=pos1+1
       pos2=INDEX(input_line(pos1:),' ')
       CPPostcondition(pos2>=0,cp_failure_level,routineP,error,failure)
       IF(pos2==0) pos2 = LEN_TRIM(input_line(pos1:))+1

       pos2=pos1+pos2-1
       varname=input_line(pos1:pos2)
       idx=inpp_find_variable(inpp,varname)
       CPPostcondition(idx>0,cp_failure_level,routineP,error,failure)
       IF(idx==0) THEN
          WRITE(UNIT=message,FMT="(5A,I6)") &
               "INPP_@SET: Variable $",TRIM(varname)," not defined in file: ",&
               TRIM(input_file_name), "  Line:",input_line_number
          CALL stop_program(routineN,moduleN,__LINE__,TRIM(message))
       ENDIF

       IF(pos1>2) THEN
          newline=input_line(1:pos1-2) // TRIM(inpp%variable_value(idx)) &
               // input_line(pos2+1:)
       ELSE
          newline=TRIM(inpp%variable_value(idx)) // input_line(pos2+1:)
       ENDIF

       input_line=newline
    END DO
  END SUBROUTINE inpp_expand_variables

! *****************************************************************************
!> \brief return index position of a variable in dictionary. 0 if not found.
!> \par History
!>      - standalone proof-of-concept implemenation (22.02.2008,AK)
!>      - integrated into cp2k (23.02.2008)
!> \author AK
! *****************************************************************************
  FUNCTION inpp_find_variable(inpp, varname) RESULT(idx)
    TYPE(inpp_type), POINTER                 :: inpp
    CHARACTER(len=default_path_length)       :: varname
    INTEGER                                  :: idx

    INTEGER                                  :: i

    idx=0
    DO i=1,inpp%num_variables
       IF(TRIM(varname) == TRIM(inpp%variable_name(i))) THEN
          idx=i
          RETURN
       ENDIF
    END DO
    RETURN
  END FUNCTION inpp_find_variable

! *****************************************************************************
!> \brief print a list of the variable/value table
!> \par History
!>      - standalone proof-of-concept implemenation (22.02.2008,AK)
!>      - integrated into cp2k (23.02.2008)
!> \author AK
! *****************************************************************************
  SUBROUTINE inpp_list_variables(inpp,iochan)
    TYPE(inpp_type), POINTER                 :: inpp
    INTEGER, INTENT(IN)                      :: iochan

    INTEGER                                  :: i

    WRITE(iochan,'(A)') '   #   NAME                   VALUE'
    DO i=1,inpp%num_variables
       WRITE(iochan,'(I4," | ",A,T30," | ",A," |")') &
            i,TRIM(inpp%variable_name(i)),TRIM(inpp%variable_value(i))
    END DO
  END SUBROUTINE inpp_list_variables

END MODULE cp_parser_inpp_methods
