!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2006  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/checkpoint [1.0] *
!!
!!   NAME
!!     checkpoint
!!
!!   FUNCTION
!!     Checkpoint handling.
!!
!!   AUTHOR
!!     Matthias Krack (26.09.2002)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE checkpoint_handler
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_cputime,&
                                             m_flush,&
                                             m_getpid
  USE memory_utilities,                ONLY: write_memory
  USE message_passing,                 ONLY: mp_sum,&
                                             mp_sync
  USE string_utilities,                ONLY: ascii_to_string,&
                                             compress,&
                                             string_to_ascii
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: write_checkpoint_information,&
            Pwarning

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE write_checkpoint_information(checkpoint,root_section,para_env,error)

!   Purpose: Write checkpoint information.

!   History: - Creation (29.09.2002,MK)

!   ***************************************************************************

    CHARACTER(LEN=*), INTENT(IN)             :: checkpoint
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: iw
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,root_section,"GLOBAL%PERFORMANCE/ROUTINES",&
         extension=".Log",error=error)
!   ---------------------------------------------------------------------------
!   *** Wait for ionode and empty IO buffer ***
    IF (iw>0) WRITE(iw,*)checkpoint
    IF (iw>0) CALL m_flush(iw)
    CALL mp_sync(para_env%group)

!   *** Write checkpoint informations ***

    CALL write_cputime(checkpoint,root_section,para_env,error)
    CALL write_memory(checkpoint,root_section,para_env,error)

    IF (iw>0) CALL m_flush(iw)
    CALL mp_sync(para_env%group)
    CALL cp_print_key_finished_output(iw,logger,root_section,&
         "GLOBAL%PERFORMANCE/ROUTINES",error=error)    
  END SUBROUTINE write_checkpoint_information

! *****************************************************************************

  SUBROUTINE write_cputime(checkpoint,root_section,para_env,error)

!   Purpose: Write the current CPU time of all processes in seconds.

!   History: - Creation (25.09.2002,MK)

!   ***************************************************************************

    CHARACTER(LEN=*), INTENT(IN)             :: checkpoint
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=LEN_TRIM(checkpoint)+55)   :: string2
    CHARACTER(LEN=LEN_TRIM(checkpoint)+7)    :: string1
    INTEGER                                  :: group, ipe, iw, mype, npe, &
                                                output_unit, pid
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: myconvertedstring
    LOGICAL                                  :: ionode
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    IF (BTEST(cp_print_key_should_output(logger%iter_info,root_section,&
         "GLOBAL%PERFORMANCE/CPUTIME",error=error),cp_p_file)) THEN
       iw = cp_print_key_unit_nr(logger,root_section,"GLOBAL%PERFORMANCE/CPUTIME",&
            extension=".Log",error=error)

       group = para_env%group
       ionode = para_env%ionode
       npe = para_env%num_pe
       mype = para_env%mepos
       output_unit = iw
       
       ALLOCATE(myconvertedstring(LEN(string2),npe))
       myconvertedstring=0
       
       CALL m_getpid(pid)
       WRITE (UNIT=string1,FMT="(I6,A)") mype,":"//TRIM(checkpoint)
       string1 = ADJUSTL(string1)
       WRITE (UNIT=string2,FMT="(A,I8,A,F12.3,A)")&
            "CPU time of process",pid," <"//TRIM(string1)//">: ",m_cputime(),&
            " sec"
       CALL string_to_ascii(string2,myconvertedstring(:,mype+1))
       CALL mp_sum(myconvertedstring,group)
       IF (iw>0) THEN
          DO ipe=1,npe
             CALL ascii_to_string(myconvertedstring(:,ipe),string2)
             WRITE (UNIT=output_unit,FMT="(T2,A)") TRIM(string2)
          ENDDO
          CALL m_flush(output_unit)
       ENDIF       
       DEALLOCATE(myconvertedstring)
       CALL cp_print_key_finished_output(iw,logger,root_section,&
            "GLOBAL%PERFORMANCE/CPUTIME",error=error)       
    END IF
  END SUBROUTINE write_cputime

!!*****
!******************************************************************************
  SUBROUTINE Pwarning(message, input_file, error)

!   Purpose: Write warning messages

!   History: - Creation (01.2006, Teo)

!   ***************************************************************************

    CHARACTER(LEN=*), INTENT(IN)             :: message
    TYPE(section_vals_type), POINTER         :: input_file
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=1024)                      :: my_message
    INTEGER                                  :: cyc, i, iw, j
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,input_file,"GLOBAL%WARNINGS",&
         extension=".Log",error=error)
    IF (iw>0) THEN
       my_message = message(1:MIN(LEN(message),1024))
       CALL compress(my_message)
       cyc = FLOOR(REAL(LEN_TRIM(my_message),KIND=dp) / 69.0_dp) + 1
       WRITE(iw,'(/)')
       DO i = 1, cyc
          j =0
          IF (my_message((i-1)*69+1:(i-1)*69+1)==" ") j = 1
          WRITE(iw,'(T2,A,2X,A)')"WARNING|",my_message((i-1)*69+1+j:MIN((i-1)*69+69+j,LEN_TRIM(my_message)))
       END DO
       WRITE(iw,'(/)')
    END IF
    CALL cp_print_key_finished_output(iw,logger,input_file,&
         "GLOBAL%WARNINGS",error=error)    
  END SUBROUTINE Pwarning

END MODULE checkpoint_handler
