!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
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

  USE kinds,                           ONLY: dp
  USE global_types,                    ONLY: global_environment_type
  USE machine,                         ONLY: m_cputime,&
                                             m_flush,&
                                             m_getpid
  USE memory_utilities,                ONLY: write_memory
  USE message_passing,                 ONLY: mp_sync,&
                                             mp_sum
  USE string_utilities,                ONLY: compress,&
                                             string_to_ascii,&
                                             ascii_to_string
  USE cp_output_handling,              ONLY: cp_print_key_unit_nr, &
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output, cp_p_file

#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: write_checkpoint_information,&
            Pwarning

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE write_checkpoint_information(checkpoint,globenv,error)

!   Purpose: Write checkpoint information.

!   History: - Creation (29.09.2002,MK)

!   ***************************************************************************

    CHARACTER(LEN=*), INTENT(IN)                :: checkpoint
    TYPE(global_environment_type), &
         POINTER                                :: globenv
    TYPE(cp_error_type), INTENT(inout)          :: error

    TYPE(cp_logger_type), POINTER            :: logger
    INTEGER :: iw
    
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,globenv%input_file,"GLOBAL%PERFORMANCE/ROUTINES",&
         extension=".Log",error=error)
!   ---------------------------------------------------------------------------
!   *** Wait for ionode and empty IO buffer ***
    
    IF (iw>0) WRITE(iw,*)checkpoint
    IF (globenv%ionode) CALL m_flush(globenv%scr)
    CALL mp_sync(globenv%group)

!   *** Write checkpoint informations ***

    CALL write_cputime(checkpoint,globenv,error)
    CALL write_memory(checkpoint,globenv,error)

    IF (globenv%ionode) CALL m_flush(globenv%scr)
    CALL mp_sync(globenv%group)
    CALL cp_print_key_finished_output(iw,logger,globenv%input_file,&
         "GLOBAL%PERFORMANCE/ROUTINES",error=error)    
  END SUBROUTINE write_checkpoint_information

! *****************************************************************************

  SUBROUTINE write_cputime(checkpoint,globenv,error)

!   Purpose: Write the current CPU time of all processes in seconds.

!   History: - Creation (25.09.2002,MK)

!   ***************************************************************************

    CHARACTER(LEN=*), INTENT(IN)             :: checkpoint
    TYPE(global_environment_type), &
      pointer                             :: globenv

    CHARACTER(LEN=LEN_TRIM(checkpoint)+55)   :: string2
    CHARACTER(LEN=LEN_TRIM(checkpoint)+7)    :: string1
    INTEGER                                  :: group, ipe, mype, npe, &
                                                output_unit, pid
    INTEGER, DIMENSION(:,:), ALLOCATABLE     :: myconvertedstring
    LOGICAL                                  :: ionode
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(cp_logger_type), POINTER            :: logger
    INTEGER :: iw
    
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    IF (BTEST(cp_print_key_should_output(logger%iter_info,globenv%input_file,&
         "GLOBAL%PERFORMANCE/CPUTIME",error=error),cp_p_file)) THEN
       iw = cp_print_key_unit_nr(logger,globenv%input_file,"GLOBAL%PERFORMANCE/CPUTIME",&
            extension=".Log",error=error)

       group = globenv%group
       ionode = globenv%ionode
       npe = globenv%num_pe
       mype = globenv%mepos
       output_unit = globenv%scr
       
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
       CALL cp_print_key_finished_output(iw,logger,globenv%input_file,&
            "GLOBAL%PERFORMANCE/CPUTIME",error=error)       
    END IF
  END SUBROUTINE write_cputime

!!*****
!******************************************************************************
  SUBROUTINE Pwarning(message, globenv, error)

!   Purpose: Write warning messages

!   History: - Creation (01.2006, Teo)

!   ***************************************************************************

    CHARACTER(LEN=*), INTENT(IN)                :: message
    TYPE(global_environment_type), &
         POINTER                                :: globenv
    TYPE(cp_error_type), INTENT(inout)          :: error

    TYPE(cp_logger_type), POINTER               :: logger
    INTEGER :: iw, cyc, i, j
    CHARACTER(LEN=1024)                         :: my_message

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,globenv%input_file,"GLOBAL%WARNINGS",&
         extension=".Log",error=error)
!   ---------------------------------------------------------------------------
    
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
    CALL cp_print_key_finished_output(iw,logger,globenv%input_file,&
         "GLOBAL%WARNINGS",error=error)    
  END SUBROUTINE Pwarning

END MODULE checkpoint_handler
