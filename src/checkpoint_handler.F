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
  USE machine,                         ONLY: m_flush
  USE memory_utilities,                ONLY: write_memory
  USE message_passing,                 ONLY: mp_sync
  USE timings,                         ONLY: write_cputime

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

END MODULE checkpoint_handler
