!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2005                                                        !
!-----------------------------------------------------------------------------!
!!****** cp2k/cp_files [1.0] *
!!
!!   NAME
!!     cp_files
!!
!!   FUNCTION
!!     Utility routines to open and close files
!!
!!   AUTHOR
!!     CP2K_WORKSHOP 1.0 TEAM
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
MODULE cp_files
  USE string_utilities, ONLY : uppercase
  IMPLICIT NONE
  PRIVATE


  PUBLIC :: close_file, &
            open_file, &
            get_unit_number

!!***
  INTEGER, PARAMETER, PRIVATE :: max_message_length=400
  INTEGER, PARAMETER, PRIVATE :: max_unit_number=999
! *****************************************************************************
  CHARACTER(LEN=*), PARAMETER :: module_name = "cp_files"
  INTEGER, DIMENSION(2), PARAMETER :: reserved_unit_numbers = (/5,6/)
CONTAINS
  !!****f* cp_log_handling/close_file *
  !!
  !!   NAME
  !!     close_file
  !!
  !!   FUNCTION
  !!     closes the given unit
  !!
  !!   NOTES
  !!     -
  !!
  !!   ARGUMENTS
  !!     -
  !!
  !!   AUTHOR
  !!     MK
  !!
  !!*** *********************************************************************
  SUBROUTINE close_file(unit_number,file_status)

    INTEGER, INTENT(IN)                      :: unit_number
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: file_status

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "close_file"

    CHARACTER(LEN=6)                         :: status_string
    CHARACTER(LEN=max_message_length)        :: message
    INTEGER                                  :: istat
    LOGICAL                                  :: exists, opened

!   *** Check the specified input file name ***

    INQUIRE (UNIT=unit_number,EXIST=exists,OPENED=opened,IOSTAT=istat)

    IF (istat /= 0) THEN
      WRITE (UNIT=message,FMT="(A,I6,A,I6,A)")&
        "An error occurred inquiring the unit with the number ",unit_number,&
        " (IOSTAT = ",istat,")"
      WRITE(6,*) TRIM(message)
      STOP  ! should actually do error recovery....
    ELSE IF (.NOT.exists) THEN
      WRITE (UNIT=message,FMT="(A,I6,A,I6,A)")&
        "The specified unit number ",unit_number," does not exist"
      WRITE(6,*) TRIM(message)
      STOP  ! should actually do error recovery....
    END IF

!   *** Close the specified file ***

    IF (opened) THEN

      IF (PRESENT(file_status)) THEN
        status_string = file_status
        CALL uppercase(status_string)
      ELSE
        status_string = "KEEP"
      END IF

      CLOSE (UNIT=unit_number,IOSTAT=istat,STATUS=TRIM(status_string))

      IF (istat /= 0) THEN
        WRITE (UNIT=message,FMT="(A,I6,A,I6,A)")&
          "An error occurred closing the file with the unit number ",unit_number," (IOSTAT = ",istat,")"
        WRITE(6,*) TRIM(message)
        STOP  ! should actually do error recovery....
      END IF

    END IF

  END SUBROUTINE close_file

!!****f* cp_files/get_unit_number *
!!
!!   NAME
!!      get_unit_number
!!   FUNCTION
!!      returns the first fortran unit that is not preconnected
!!   NOTES
!!      -1 if no free unit exists
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************
  FUNCTION get_unit_number() RESULT(unit_number)
    INTEGER                                  :: unit_number

    INTEGER                                  :: istat
    LOGICAL                                  :: exists, opened

    DO unit_number=1,max_unit_number
      IF (ANY(unit_number == reserved_unit_numbers)) CYCLE
      INQUIRE (UNIT=unit_number,EXIST=exists,OPENED=opened,IOSTAT=istat)
      IF (exists.AND.(.NOT.opened).AND.(istat == 0)) RETURN
    END DO

    unit_number = -1

   END FUNCTION get_unit_number
  !!****f* cp_files/open_file *
  !!
  !!   NAME
  !!     open_file
  !!
  !!   FUNCTION
  !!     opens the requested file using a free unit number
  !!
  !!   NOTES
  !!     -
  !!
  !!   ARGUMENTS
  !!     -
  !!
  !!   AUTHOR
  !!     MK
  !!
  !!*** *********************************************************************
  SUBROUTINE open_file(file_name,file_status,file_form,file_action,&
                       file_position,unit_number)

    CHARACTER(LEN=*), INTENT(IN)             :: file_name
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: file_status, file_form, &
                                                file_action, file_position
    INTEGER, INTENT(OUT)                     :: unit_number

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "open_file"

    CHARACTER(LEN=11)                        :: action_string, form_string, &
                                                position_string, status_string
    CHARACTER(LEN=max_message_length)        :: message
    INTEGER                                  :: istat
    LOGICAL                                  :: exists, opened

!   ---------------------------------------------------------------------------

    IF (PRESENT(file_status)) THEN
      status_string = file_status
      CALL uppercase(status_string)
    ELSE
      status_string = "OLD"
    END IF

    IF (PRESENT(file_form)) THEN
      form_string = file_form
      CALL uppercase(form_string)
    ELSE
      form_string = "FORMATTED"
    END IF

    IF (PRESENT(file_action)) THEN
      action_string = file_action
      CALL uppercase(action_string)
    ELSE
      action_string = "READ"
    END IF

    IF (PRESENT(file_position)) THEN
      position_string = file_position
      CALL uppercase(position_string)
    ELSE
      position_string = "REWIND"
    END IF

!   *** Check the specified input file name ***

    INQUIRE (FILE=TRIM(file_name),EXIST=exists,OPENED=opened,IOSTAT=istat)

    IF (istat /= 0) THEN
      WRITE (UNIT=message,FMT="(A,I6,A)")&
        "An error occurred inquiring the file <"//TRIM(file_name)//&
        "> (IOSTAT = ",istat,")"
      WRITE(6,*) TRIM(message)
      STOP  ! should actually do error recovery....
    ELSE IF (status_string == "OLD") THEN
      IF (.NOT.exists) THEN
         write(6,*) "The specified file "//TRIM(file_name)//" can not be opened, it does not exist"
         STOP  ! should actually do error recovery....
      END IF
    END IF

!   *** Open the specified input file ***

    IF (opened) THEN

       INQUIRE (FILE=TRIM(file_name),NUMBER=unit_number)
       REWIND (UNIT=unit_number)
             
    ELSE

!     *** Find an unused unit number ***

      unit_number = get_unit_number()

      IF (unit_number < 0) THEN
         write(6,*) "Problems opening file, there are no free units left"
         STOP  ! should actually do error recovery....
      END IF

      OPEN (UNIT=unit_number,&
            FILE=TRIM(file_name),&
            STATUS=TRIM(status_string),&
            ACCESS="SEQUENTIAL",&
            FORM=TRIM(form_string),&
            POSITION=TRIM(position_string),&
            ACTION=TRIM(action_string),&
            IOSTAT=istat)

      IF (istat /= 0) THEN
        WRITE (UNIT=message,FMT="(A,I6,A,I6,A)")&
          "An error occurred opening the file <"//TRIM(file_name)//&
          "> with the unit number ",unit_number," (IOSTAT = ",istat,")"
        WRITE(6,*) TRIM(message)
        STOP  ! should actually do error recovery....
      END IF

    END IF

  END SUBROUTINE open_file
END MODULE cp_files

