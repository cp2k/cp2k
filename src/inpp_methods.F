!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!
 
! *****************************************************************************
!> \brief a module to allow simple internal preprocessing in input files.
!> \par History
!>      - standalone proof-of-concept implemenation (20.02.2008,AK)
!> \author AK
! *****************************************************************************
MODULE inpp_methods
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_parser_types,                 ONLY: cp_parser_type
  USE inpp_types,                      ONLY: inpp_type
  USE kinds,                           ONLY: default_path_length, default_string_length
  USE memory_utilities,                ONLY: reallocate
  USE string_utilities,                ONLY: uppercase
#include "cp_common_uses.h"

  IMPLICIT NONE
  
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'inpp_methods'

  PUBLIC :: inpp_start_include, inpp_end_include
  
CONTAINS

! *****************************************************************************
!> \brief process an '@INCLUDE "file.inc"' statement. Store old status on stack.
!> \par History
!>      - standalone proof-of-concept implemenation (20.02.2008,AK)
!> \author AK
! *****************************************************************************
  SUBROUTINE inpp_start_include(inpp, parser, error)
    TYPE(inpp_type), POINTER                 :: inpp
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'inpp_start_include', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: filename
    CHARACTER(LEN=default_string_length)     :: mytag
    INTEGER                                  :: i, indf, indi, pos1, pos2, &
                                                unit
    LOGICAL                                  :: check, failure

    failure = .TRUE.
    CPPostcondition(ASSOCIATED(inpp),cp_failure_level,routineP,error,failure)
    CPPostcondition(ASSOCIATED(parser),cp_failure_level,routineP,error,failure)

    ! Get the first instruction and check for "@INCLUDE"
    indi = INDEX(parser%input_line,"@")
    indf = indi
    DO WHILE (parser%input_line(indf:indf)/=" ")
       indf = indf + 1
       CPPostcondition(indf<=LEN_TRIM(parser%input_line),cp_failure_level,routineP,error,failure)
    END DO
    CPPostcondition((indf-indi)<=default_string_length,cp_failure_level,routineP,error,failure)
    mytag =parser%input_line(indi:indf-1)
    CALL uppercase(mytag)
    SELECT CASE(mytag)
    CASE ("@INCLUDE")
       inpp%io_stack_level=inpp%io_stack_level+1
       CALL reallocate(inpp%io_stack_channel,1,inpp%io_stack_level)
       CALL reallocate(inpp%io_stack_lineno,1,inpp%io_stack_level)
       CALL reallocate(p_long=inpp%io_stack_filename,lb_new=1,ub_new=inpp%io_stack_level)

       ! save current status
       inpp%io_stack_channel(inpp%io_stack_level)=parser%input_unit
       inpp%io_stack_lineno(inpp%io_stack_level)=parser%input_line_number
       inpp%io_stack_filename(inpp%io_stack_level)=parser%input_file_name

       ! Get the filename.. allow for " or ' or nothing..
       filename = parser%input_line(indf:)
       indf = 1
       DO WHILE (filename(indf:indf)==" ")
          indf = indf + 1
          CPPostcondition(indf<=LEN_TRIM(filename),cp_failure_level,routineP,error,failure)
       END DO
       filename = filename(indf:)

       pos1=INDEX(filename,'"')
       pos2=INDEX(filename(pos1+1:),'"')
       check = (pos1==0).EQV.(pos2==0)
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
       IF (pos1/=0) filename = filename(pos1+1:pos1+pos2-1)
 
       pos1=INDEX(filename,"'")
       pos2=INDEX(filename(pos1+1:),"'")
       check = (pos1==0).EQV.(pos2==0)
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
       IF (pos1/=0) filename = filename(pos1+1:pos1+pos2-1)

       CALL open_file(file_name=TRIM(filename),&
                      file_status="OLD",&
                      file_form="FORMATTED",&
                      file_action="READ",&
                      unit_number=unit)

       ! Let's check that files already opened won't be again opened
       DO i = 1, inpp%io_stack_level
          check = TRIM(filename)/=TRIM(inpp%io_stack_filename(i))
          CPPostcondition(check,cp_failure_level,routineP,error,failure)
       END DO
       parser%input_file_name=TRIM(filename)
       parser%input_line_number=0
       parser%input_unit=unit
    CASE DEFAULT
       ! Do Nothing..
    END SELECT
 
  END SUBROUTINE inpp_start_include
  

! *****************************************************************************
!> \brief Restore older file status from stack after EOF on include file.
!> \par History
!>      - standalone proof-of-concept implemenation (20.02.2008,AK)
!> \author AK
! *****************************************************************************
  SUBROUTINE inpp_end_include(inpp, parser, error)
    TYPE(inpp_type), POINTER                 :: inpp
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'inpp_end_include', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(inpp),cp_failure_level,routineP,error,failure)
    CPPostcondition(ASSOCIATED(parser),cp_failure_level,routineP,error,failure)
    IF(inpp%io_stack_level > 0) THEN
       CALL close_file(parser%input_unit)
       parser%input_unit=inpp%io_stack_channel(inpp%io_stack_level)
       parser%input_line_number=inpp%io_stack_lineno(inpp%io_stack_level)
       parser%input_file_name=TRIM(inpp%io_stack_filename(inpp%io_stack_level))
       inpp%io_stack_level=inpp%io_stack_level-1
       CALL reallocate(inpp%io_stack_channel,1,inpp%io_stack_level)
       CALL reallocate(inpp%io_stack_lineno,1,inpp%io_stack_level)
       CALL reallocate(p_long=inpp%io_stack_filename,lb_new=1,ub_new=inpp%io_stack_level)
    ENDIF

  END SUBROUTINE inpp_end_include

END MODULE inpp_methods
