!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
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
  USE machine,                         ONLY: default_input_unit,&
                                             default_output_unit,&
                                             m_flush

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: close_file, &
            open_file, &
            get_unit_number,&
            get_eof_stat

!!***
  LOGICAL, SAVE, PRIVATE :: did_init=.FALSE.
  INTEGER, SAVE, PRIVATE :: eof_stat
  INTEGER, PARAMETER, PRIVATE :: max_message_length=400
  INTEGER, PARAMETER, PRIVATE :: max_unit_number=999
! *****************************************************************************
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_files'
  INTEGER, DIMENSION(2), PARAMETER :: reserved_unit_numbers = (/default_input_unit,default_output_unit/)
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

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'close_file'

    CHARACTER(LEN=6)                         :: status_string
    CHARACTER(LEN=max_message_length)        :: message
    INTEGER                                  :: istat
    LOGICAL                                  :: exists, opened

!   *** Check the specified input file name ***

    INQUIRE (UNIT=unit_number,EXIST=exists,OPENED=opened,IOSTAT=istat)

    IF (istat /= 0) THEN
      WRITE (UNIT=message,FMT="(A,I6,A,I6,A)")&
        "CP2K: An error occurred inquiring the unit with the number ",unit_number,&
        " (IOSTAT = ",istat,")"
      WRITE(default_output_unit,*) TRIM(message)
      CALL m_flush(default_output_unit)
      STOP 1  ! should actually do error recovery....
    ELSE IF (.NOT.exists) THEN
      WRITE (UNIT=message,FMT="(A,I6,A,I6,A)")&
        "CP2K: The specified unit number ",unit_number," can not be closed, it does not exist."
      WRITE(default_output_unit,*) TRIM(message)
      CALL m_flush(default_output_unit)
      STOP 1  ! should actually do error recovery....
    END IF

!   *** Close the specified file ***

    IF (opened) THEN
       IF (unit_number==6) STOP "tryed to close unit 6"

      IF (PRESENT(file_status)) THEN
        status_string = file_status
      ELSE
        status_string = "KEEP"
      END IF

      CLOSE (UNIT=unit_number,IOSTAT=istat,STATUS=TRIM(status_string))

      IF (istat /= 0) THEN
        WRITE (UNIT=message,FMT="(A,I6,A,I6,A)")&
          "CP2K: An error occurred closing the file with the unit number ",unit_number," (IOSTAT = ",istat,")"
        WRITE(default_output_unit,*) TRIM(message)
      CALL m_flush(default_output_unit)
        STOP 1 ! should actually do error recovery....
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
                       file_position,file_pad,unit_number)

    CHARACTER(LEN=*), INTENT(IN)             :: file_name
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: file_status, file_form, &
                                                file_action, file_position, &
                                                file_pad
    INTEGER, INTENT(OUT)                     :: unit_number

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'open_file'

    CHARACTER(LEN=11)                        :: action_string, form_string, &
                                                pad_string, position_string, &
                                                status_string
    CHARACTER(LEN=max_message_length)        :: message
    INTEGER                                  :: istat
    LOGICAL                                  :: exists, opened

!   ---------------------------------------------------------------------------

    IF (PRESENT(file_status)) THEN
      status_string = file_status
    ELSE
      status_string = "OLD"
    END IF

    IF (PRESENT(file_form)) THEN
      form_string = file_form
    ELSE
      form_string = "FORMATTED"
    END IF

    IF (PRESENT(file_pad)) THEN
       pad_string = file_pad
       IF (form_string=="UNFORMATTED") THEN
          WRITE (UNIT=message,FMT="(A)")&
               "CP2K: The PAD= specifier is not allowed for an UNFORMATTED file!"
          WRITE(default_output_unit,*) TRIM(message)
      CALL m_flush(default_output_unit)
          STOP 1  ! should actually do error recovery....          
       END IF
    ELSE
       pad_string = "YES"
    END IF    

    IF (PRESENT(file_action)) THEN
      action_string = file_action
    ELSE
      action_string = "READ"
    END IF

    IF (PRESENT(file_position)) THEN
      position_string = file_position
    ELSE
      position_string = "REWIND"
    END IF

!   *** Check the specified input file name ***

    INQUIRE (FILE=TRIM(file_name),EXIST=exists,OPENED=opened,IOSTAT=istat)

    IF (istat /= 0) THEN
      WRITE (UNIT=message,FMT="(A,I6,A)")&
        "CP2K: An error occurred inquiring the file <"//TRIM(file_name)//&
        "> (IOSTAT = ",istat,")"
      WRITE(default_output_unit,*) TRIM(message)
      CALL m_flush(default_output_unit)
      STOP 1  ! should actually do error recovery....
    ELSE IF (status_string == "OLD") THEN
      IF (.NOT.exists) THEN
         WRITE(default_output_unit,*) "CP2K: The specified file "//TRIM(file_name)//" can not be opened, it does not exist."
      CALL m_flush(default_output_unit)
         STOP 1  ! should actually do error recovery....
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
         WRITE(default_output_unit,*) "CP2K: Problems opening file, there are no free units left"
      CALL m_flush(default_output_unit)
         STOP 1  ! should actually do error recovery....
      END IF

      IF (TRIM(form_string)== "FORMATTED") THEN
         OPEN (UNIT=unit_number,&
              FILE=TRIM(file_name),&
              STATUS=TRIM(status_string),&
              ACCESS="SEQUENTIAL",&
              FORM=TRIM(form_string),&
              POSITION=TRIM(position_string),&
              ACTION=TRIM(action_string),&
              PAD=TRIM(pad_string),&
              IOSTAT=istat)
      ELSE
         OPEN (UNIT=unit_number,&
              FILE=TRIM(file_name),&
              STATUS=TRIM(status_string),&
              ACCESS="SEQUENTIAL",&
              FORM=TRIM(form_string),&
              POSITION=TRIM(position_string),&
              ACTION=TRIM(action_string),&
              IOSTAT=istat)
      END IF

      IF (istat /= 0) THEN
        WRITE (UNIT=message,FMT="(A,I6,A,I6,A)")&
          "CP2K: An error occurred opening the file <"//TRIM(file_name)//&
          "> with the unit number ",unit_number," (IOSTAT = ",istat,")"
        WRITE(default_output_unit,*) TRIM(message)
      CALL m_flush(default_output_unit)
        STOP 1  ! should actually do error recovery....
      END IF

    END IF

  END SUBROUTINE open_file

!!****f* cp_files/get_eof_stat *
!!
!!   NAME
!!      get_eof_stat
!!   FUNCTION
!!      returns the number returned by iostat on eof
!!
!!*** **********************************************************************
  FUNCTION get_eof_stat() RESULT(res)
    INTEGER                                  :: res

    CHARACTER(len=5)                         :: rdBuf
    INTEGER                                  :: stat_n, unit

    IF (.not.did_init) THEN
       unit=get_unit_number()

       WRITE(unit,fmt=*)'z'
       REWIND(unit)

       eof_stat=0
       stat_n=0
       DO
          READ (unit,'(a1)', iostat=stat_n)rdBuf
          IF (stat_n/=0) EXIT
       END DO
       READ (unit,'(a1)', iostat=eof_stat)rdBuf
       IF (eof_stat >= 0) THEN
          STOP 'ERROR,cp_files:get_eof_stat failed identification of EOF IOSTAT'
          eof_stat=-1
       END IF
       CALL close_file(unit,'delete')
       did_init=.TRUE.
    END IF
    res=eof_stat
  END FUNCTION get_eof_stat
    
END MODULE cp_files

