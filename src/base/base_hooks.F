!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Central dispatch for basic hooks
!> \author Ole Schuett
! *****************************************************************************
MODULE base_hooks
  USE kinds,                           ONLY: default_string_length
  USE machine,                         ONLY: default_output_unit,&
                                             m_abort,&
                                             m_flush

  IMPLICIT NONE
  PRIVATE
  
  !API
  PUBLIC :: cp_abort, cp_warn, timeset, timestop
  !API
  PUBLIC :: cp_abort_hook, cp_warn_hook, timeset_hook, timestop_hook
  !API
  PUBLIC :: cp__a, cp__b, cp__w, cp__l

  ! this interface (with subroutines in it) musst to be defined right before
  ! the regular subroutines/functions - otherwise prettify.py will screw up.
  INTERFACE
    SUBROUTINE cp_abort_interface(location, message)
    CHARACTER(len=*), INTENT(in)             :: location, message

    END SUBROUTINE cp_abort_interface

    SUBROUTINE cp_warn_interface(location, message)
    CHARACTER(len=*), INTENT(in)             :: location, message

    END SUBROUTINE cp_warn_interface

    SUBROUTINE timeset_interface(routineN, handle)
    CHARACTER(LEN=*), INTENT(IN)             :: routineN
    INTEGER, INTENT(OUT)                     :: handle

    END SUBROUTINE timeset_interface

    SUBROUTINE timestop_interface(handle)
    INTEGER, INTENT(IN)                      :: handle

    END SUBROUTINE timestop_interface
  END INTERFACE

  PROCEDURE(cp_abort_interface), POINTER :: cp_abort_hook => Null()
  PROCEDURE(cp_warn_interface), POINTER :: cp_warn_hook => Null()
  PROCEDURE(timeset_interface), POINTER :: timeset_hook => Null()
  PROCEDURE(timestop_interface), POINTER :: timestop_hook => Null()

CONTAINS

! *****************************************************************************
!> \brief Terminate the program
!> \param location ...
!> \param message ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE cp_abort(location, message)
    CHARACTER(len=*), INTENT(in)             :: location, message

    IF(ASSOCIATED(cp_abort_hook)) THEN
       CALL cp_abort_hook(location, message)
    ELSE
       WRITE (default_output_unit,*) "ABORT in "//TRIM(location)//" "//TRIM(message)
       CALL m_flush(default_output_unit)
       CALL m_abort()
    ENDIF
    ! compiler hint
    STOP "Never return from here"
  END SUBROUTINE cp_abort

! *****************************************************************************
!> \brief Issue a warning
!> \param location ...
!> \param message ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE cp_warn(location, message)
    CHARACTER(len=*), INTENT(in)             :: location, message

    IF(ASSOCIATED(cp_warn_hook)) THEN
       CALL cp_warn_hook(location, message)
    ELSE
       WRITE (default_output_unit,*) "WARNING in "//TRIM(location)//" "//TRIM(message)
       CALL m_flush(default_output_unit)
    ENDIF
  END SUBROUTINE cp_warn

! *****************************************************************************
!> \brief Start timer
!> \param routineN ...
!> \param handle ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE timeset(routineN, handle)
    CHARACTER(LEN=*), INTENT(IN)             :: routineN
    INTEGER, INTENT(OUT)                     :: handle

    IF(ASSOCIATED(timeset_hook)) THEN
       CALL timeset_hook(routineN, handle)
    ELSE
       handle = -1
    ENDIF
  END SUBROUTINE timeset

! *****************************************************************************
!> \brief Stop timer
!> \param handle ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE timestop(handle)
    INTEGER, INTENT(IN)                      :: handle

    IF(ASSOCIATED(timestop_hook)) THEN
       CALL timestop_hook(handle)
    ELSE
       IF(handle/=-1) &
          CALL cp_abort(cp__l("base_hooks.F",__LINE__),"Got wrong handle")
    ENDIF
  END SUBROUTINE timestop

! *****************************************************************************
!> \brief CPASSERT handler
!> \param filename ...
!> \param lineNr ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE cp__a(filename,lineNr)
    CHARACTER(len=*), INTENT(in)             :: filename
    INTEGER, INTENT(in)                      :: lineNr

    CALL cp_abort(location=cp__l(filename,lineNr), message="CPASSERT failed")
    ! compiler hint
    STOP "Never return from here"
  END SUBROUTINE cp__a

! *****************************************************************************
!> \brief CPABORT handler
!> \param filename ...
!> \param lineNr ...
!> \param message ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE cp__b(filename,lineNr,message)
    CHARACTER(len=*), INTENT(in)             :: filename
    INTEGER, INTENT(in)                      :: lineNr
    CHARACTER(len=*), INTENT(in)             :: message

    CALL cp_abort(location=cp__l(filename,lineNr), message=message)
    ! compiler hint
    STOP "Never return from here"
  END SUBROUTINE cp__b

! *****************************************************************************
!> \brief CPWARN handler
!> \param filename ...
!> \param lineNr ...
!> \param message ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE cp__w(filename, lineNr,message)
    CHARACTER(len=*), INTENT(in)             :: filename
    INTEGER, INTENT(in)                      :: lineNr
    CHARACTER(len=*), INTENT(in)             :: message

    CALL cp_warn(location=cp__l(filename,lineNr), message=message)
  END SUBROUTINE cp__w

! *****************************************************************************
!> \brief Helper routine to assemble __LOCATION__
!> \param filename ...
!> \param lineNr ...
!> \retval location ...
!> \author Ole Schuett
! *****************************************************************************
  FUNCTION cp__l(filename, lineNr) RESULT(location)
    CHARACTER(len=*), INTENT(in)             :: filename
    INTEGER, INTENT(in)                      :: lineNr
    CHARACTER(len=default_string_length)     :: location

    CHARACTER(len=15)                        :: lineNr_str

    WRITE (lineNr_str, FMT='(I10)') lineNr
    location = TRIM(filename)//":"//TRIM(ADJUSTL(lineNr_str))

  END FUNCTION cp__l

END MODULE base_hooks
