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

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: write_checkpoint_information

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE write_checkpoint_information(checkpoint,globenv)

!   Purpose: Write checkpoint information.

!   History: - Creation (29.09.2002,MK)

!   ***************************************************************************

    CHARACTER(LEN=*), INTENT(IN)             :: checkpoint
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

!   ---------------------------------------------------------------------------
!   *** Wait for ionode and empty IO buffer ***

    IF (globenv%ionode) CALL m_flush(globenv%scr)
    CALL mp_sync(globenv%group)

!   *** Write checkpoint informations ***

    CALL write_cputime(checkpoint,globenv)
    CALL write_memory(checkpoint,globenv)

    IF (globenv%ionode) CALL m_flush(globenv%scr)
    CALL mp_sync(globenv%group)

  END SUBROUTINE write_checkpoint_information

! *****************************************************************************

  SUBROUTINE write_cputime(checkpoint,globenv)

!   Purpose: Write the current CPU time of all processes in seconds.

!   History: - Creation (25.09.2002,MK)

!   ***************************************************************************

    CHARACTER(LEN=*), INTENT(IN)             :: checkpoint
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=LEN_TRIM(checkpoint)+55)   :: string2
    CHARACTER(LEN=LEN_TRIM(checkpoint)+7)    :: string1
    INTEGER                                  :: group, ipe, mype, npe, &
                                                output_unit, pid
    INTEGER, DIMENSION(:,:), ALLOCATABLE     :: myconvertedstring
    LOGICAL                                  :: ionode

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%print%cputime) RETURN

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
    IF (ionode) THEN
      DO ipe=1,npe
         CALL ascii_to_string(myconvertedstring(:,ipe+1),string2)
         WRITE (UNIT=output_unit,FMT="(T2,A)") TRIM(string2)
      ENDDO
      CALL m_flush(output_unit)
    ENDIF

    DEALLOCATE(myconvertedstring)

  END SUBROUTINE write_cputime

!!*****
!******************************************************************************

END MODULE checkpoint_handler
