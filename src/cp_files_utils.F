!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004 CP2K developer group
!-----------------------------------------------------------------------------!

!!****** cp2k/cp_files_utils [1.0] *
!!
!!   NAME
!!     cp_files_utils
!!
!!   FUNCTION
!!     Utility routines to manipulate files
!!
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     01.2006 created [teo]
!!
!!   SOURCE
!******************************************************************************

MODULE cp_files_utils
  USE cp_files,                        ONLY: open_file, close_file
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE cp_parser,                       ONLY: parser_release,&
                                             parser_get_object,&
                                             parser_get_next_line,&
                                             parser_search_string,&
                                             parser_create,&
                                             cp_parser_type,&
                                             parser_test_next_token
  USE global_types,                    ONLY: global_environment_type
                                            
#include "cp_common_uses.h"
  
  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: moduleN = "cp_files_utils"

  PRIVATE
  PUBLIC :: move_file

CONTAINS
  
  !******************************************************************************
  !!****** cp_files_utils/move_file [1.0] *
  !!
  !!   NAME
  !!     move_file
  !!
  !!   FUNCTION
  !!
  !!
  !!   AUTHOR
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE move_file(filename, globenv, error)
    IMPLICIT NONE
    ! Arguments
    CHARACTER (LEN=*), INTENT(IN)            :: filename
    TYPE(global_environment_type), &
         POINTER                             :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error    
    ! Local Variables
    TYPE(cp_parser_type), POINTER            :: parser    
    CHARACTER(LEN=3)                         :: advance
    CHARACTER(LEN=default_string_length)     :: string
    REAL(KIND=dp)                            :: r_num
    INTEGER                                  :: i_num, file_unit

    NULLIFY(parser)
    CALL open_file(file_name=TRIM(filename)//".bak",&
                   file_status="REPLACE",&
                   file_action="WRITE",&
                   unit_number=file_unit)

    CALL parser_create(parser,filename,para_env=globenv%para_env)    
    CALL parser_get_next_line(parser,1,error=error)
    File: DO WHILE ( .TRUE. )

       IF     (parser_test_next_token(parser) == "EOL") THEN
          IF (parser_test_next_token(parser,newline=.TRUE.) == "EOF") EXIT File
          CALL parser_get_next_line(parser,1,error=error)
       ELSEIF (parser_test_next_token(parser) == "FLT") THEN
          CALL parser_get_object(parser,r_num)
          advance="NO"
          IF (parser_test_next_token(parser) == "EOL") advance="YES"
          WRITE(file_unit,'(E30.16)',ADVANCE=advance)r_num
       ELSEIF (parser_test_next_token(parser) == "INT") THEN
          CALL parser_get_object(parser,i_num)
          advance="NO"
          IF (parser_test_next_token(parser) == "EOL") advance="YES"
          WRITE(file_unit,'(I15)',ADVANCE=advance)i_num
       ELSEIF (parser_test_next_token(parser) == "STR") THEN
          CALL parser_get_object(parser,string)
          advance="NO"
          IF (parser_test_next_token(parser) == "EOL") advance="YES"
          WRITE(file_unit,'(A)',ADVANCE=advance)TRIM(string)
       END IF

    END DO File
    CALL parser_release(parser,error=error)
    CALL close_file(unit_number=file_unit)
  END SUBROUTINE move_file

END MODULE cp_files_utils
  
  
