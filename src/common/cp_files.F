!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Utility routines to open and close files. Tracking of preconnections.
!> \par History
!>      - Creation CP2K_WORKSHOP 1.0 TEAM
!>      - Revised (18.02.2011,MK)
!>      - Enhanced error checking (22.02.2011,MK)
!> \author Matthias Krack (MK)
! *****************************************************************************
MODULE cp_files

  
  USE kinds,                           ONLY: default_path_length
  USE machine,                         ONLY: default_input_unit,&
                                             default_output_unit
  USE termination,                     ONLY: stop_program

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: close_file,&
            init_preconnection_list,&
            open_file,&
            get_unit_number,&
            file_exists

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_files'

  INTEGER, PARAMETER :: max_preconnections = 10,&
                        max_unit_number = 999

  TYPE preconnection_type
     PRIVATE
     CHARACTER(LEN=default_path_length) :: file_name
     INTEGER                            :: unit_number
  END TYPE preconnection_type

  TYPE(preconnection_type), DIMENSION(max_preconnections) :: preconnected

CONTAINS

! *****************************************************************************
!> \brief Add an entry to the list of preconnected units
!> \param file_name ...
!> \param unit_number ...
!> \par History
!>      - Creation (22.02.2011,MK)
!> \author Matthias Krack (MK)
! *****************************************************************************
  SUBROUTINE assign_preconnection(file_name,unit_number)

    CHARACTER(LEN=*), INTENT(IN)             :: file_name
    INTEGER, INTENT(IN)                      :: unit_number

    CHARACTER(LEN=*), PARAMETER :: routineN = 'assign_preconnection', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ic, islot, nc

    IF ((unit_number < 1).OR.(unit_number > max_unit_number)) THEN
       CALL stop_program(routineN,moduleN,__LINE__,&
                         "An invalid logical unit number was specified.")
    END IF

    IF (LEN_TRIM(file_name) == 0) THEN
       CALL stop_program(routineN,moduleN,__LINE__,&
                         "No valid file name was specified")
    END IF

    nc = SIZE(preconnected)

    ! Check if a preconnection already exists
    DO ic=1,nc
       IF (TRIM(preconnected(ic)%file_name) == TRIM(file_name)) THEN
          ! Return if the entry already exists
          IF (preconnected(ic)%unit_number == unit_number) THEN
             RETURN
          ELSE
             CALL print_preconnection_list()
             CALL stop_program(routineN,moduleN,__LINE__,&
                               "Attempt to connect the already connected file <"//&
                               TRIM(ADJUSTL(file_name))//"> to another unit")
          END IF
       END IF
    END DO

    ! Seach for an unused entry
    islot = -1
    DO ic=1,nc
       IF (preconnected(ic)%unit_number == -1) THEN
          islot = ic
          EXIT
       END IF
    END DO

    IF (islot == -1) THEN
       CALL print_preconnection_list()
       CALL stop_program(routineN,moduleN,__LINE__,&
                         "No free slot found in the list of preconnected units")
    END IF

    preconnected(islot)%file_name = TRIM(ADJUSTL(file_name))
    preconnected(islot)%unit_number = unit_number

  END SUBROUTINE assign_preconnection

! *****************************************************************************
!> \brief Close an open file given by its logical unit number.
!>        Optionally, keep the file and unit preconnected.
!> \param unit_number ...
!> \param file_status ...
!> \param keep_preconnection ...
!> \author Matthias Krack (MK)
! *****************************************************************************
  SUBROUTINE close_file(unit_number,file_status,keep_preconnection)

    INTEGER, INTENT(IN)                      :: unit_number
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: file_status
    LOGICAL, INTENT(IN), OPTIONAL            :: keep_preconnection

    CHARACTER(LEN=*), PARAMETER :: routineN = 'close_file', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2*default_path_length)     :: message
    CHARACTER(LEN=6)                         :: status_string
    CHARACTER(LEN=default_path_length)       :: file_name
    INTEGER                                  :: istat
    LOGICAL                                  :: exists, is_open, &
                                                keep_file_connection

    keep_file_connection = .FALSE.
    IF (PRESENT(keep_preconnection)) keep_file_connection = keep_preconnection

    INQUIRE (UNIT=unit_number,EXIST=exists,OPENED=is_open,IOSTAT=istat)

    IF (istat /= 0) THEN
       WRITE (UNIT=message,FMT="(A,I0,A,I0,A)")&
         "An error occurred inquiring the unit with the number ",unit_number,&
         " (IOSTAT = ",istat,")"
       CALL stop_program(routineN,moduleN,__LINE__,TRIM(message))
    ELSE IF (.NOT.exists) THEN
       WRITE (UNIT=message,FMT="(A,I0,A)")&
         "The specified unit number ",unit_number,&
         " cannot be closed, because it does not exist."
       CALL stop_program(routineN,moduleN,__LINE__,TRIM(message))
    END IF

    ! Close the specified file

    IF (is_open) THEN
       ! Refuse to close any preconnected system unit
       IF (unit_number == default_input_unit) THEN
          WRITE (UNIT=message,FMT="(A,I0)")&
            "Attempt to close the default input unit number ",unit_number
          CALL stop_program(routineN,moduleN,__LINE__,TRIM(message))
       END IF
       IF (unit_number == default_output_unit) THEN
          WRITE (UNIT=message,FMT="(A,I0)")&
            "Attempt to close the default output unit number ",unit_number
          CALL stop_program(routineN,moduleN,__LINE__,TRIM(message))
       END IF
       ! Define status after closing the file
       IF (PRESENT(file_status)) THEN
          status_string = TRIM(ADJUSTL(file_status))
       ELSE
          status_string = "KEEP"
       END IF
       ! Optionally, keep this unit preconnected
       INQUIRE (UNIT=unit_number,NAME=file_name,IOSTAT=istat)
       IF (istat /= 0) THEN
          WRITE (UNIT=message,FMT="(A,I0,A,I0,A)")&
            "An error occurred inquiring the unit with the number ",unit_number,&
            " (IOSTAT = ",istat,")"
          CALL stop_program(routineN,moduleN,__LINE__,TRIM(message))
       END IF
       ! Manage preconnections
       IF (keep_file_connection) THEN
          CALL assign_preconnection(file_name,unit_number)
       ELSE
          CALL delete_preconnection(file_name,unit_number)
          CLOSE (UNIT=unit_number,IOSTAT=istat,STATUS=TRIM(status_string))
          IF (istat /= 0) THEN
             WRITE (UNIT=message,FMT="(A,I0,A,I0,A)")&
               "An error occurred closing the file with the logical unit number ",&
               unit_number," (IOSTAT = ",istat,")"
             CALL stop_program(routineN,moduleN,__LINE__,TRIM(message))
          END IF
       END IF
    END IF

  END SUBROUTINE close_file

! *****************************************************************************
!> \brief Remove an entry from the list of preconnected units
!> \param file_name ...
!> \param unit_number ...
!> \par History
!>      - Creation (22.02.2011,MK)
!> \author Matthias Krack (MK)
! *****************************************************************************
  SUBROUTINE delete_preconnection(file_name,unit_number)

    CHARACTER(LEN=*), INTENT(IN)             :: file_name
    INTEGER                                  :: unit_number

    CHARACTER(LEN=*), PARAMETER :: routineN = 'delete_preconnection', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ic, nc

    nc = SIZE(preconnected)

    ! Seach for preconnection entry and delete it when found
    DO ic=1,nc
       IF (TRIM(preconnected(ic)%file_name) == TRIM(file_name)) THEN
          IF (preconnected(ic)%unit_number == unit_number) THEN
             preconnected(ic)%file_name = ""
             preconnected(ic)%unit_number = -1
             EXIT
          ELSE
             CALL print_preconnection_list()
             CALL stop_program(routineN,moduleN,__LINE__,&
                               "Attempt to disconnect the file <"//&
                               TRIM(ADJUSTL(file_name))//&
                               "> from an unlisted unit")
          END IF
       END IF
    END DO

  END SUBROUTINE delete_preconnection

! *****************************************************************************
!> \brief Returns the first logical unit that is not preconnected
!> \param file_name ...
!> \retval unit_number ...
!> \author Matthias Krack (MK)
!> \note
!>       -1 if no free unit exists
! *****************************************************************************
  FUNCTION get_unit_number(file_name) RESULT(unit_number)

    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: file_name
    INTEGER                                  :: unit_number

    INTEGER                                  :: ic, istat, nc
    LOGICAL                                  :: exists, is_open

    IF (PRESENT(file_name)) THEN
       nc = SIZE(preconnected)
       ! Check for preconnected units
       DO ic=3,nc ! Exclude the preconnected system units (< 3)
          IF (TRIM(preconnected(ic)%file_name) == TRIM(file_name)) THEN
             unit_number = preconnected(ic)%unit_number
             RETURN
          END IF
       END DO
    END IF

    ! Get a new unit number
    DO unit_number=1,max_unit_number
       IF (ANY(unit_number == preconnected(:)%unit_number)) CYCLE
       INQUIRE (UNIT=unit_number,EXIST=exists,OPENED=is_open,IOSTAT=istat)
       IF (exists.AND.(.NOT.is_open).AND.(istat == 0)) RETURN
    END DO

    unit_number = -1

  END FUNCTION get_unit_number

! *****************************************************************************
!> \brief Allocate and initialise the list of preconnected units
!> \par History
!>      - Creation (22.02.2011,MK)
!> \author Matthias Krack (MK)
! *****************************************************************************
  SUBROUTINE init_preconnection_list()

    CHARACTER(LEN=*), PARAMETER :: routineN = 'init_preconnection_list', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ic, nc

    nc = SIZE(preconnected)

    DO ic=1,nc
       preconnected(ic)%file_name = ""
       preconnected(ic)%unit_number = -1
    END DO

    ! Define reserved unit numbers
    preconnected(1)%file_name = "stdin"
    preconnected(1)%unit_number = default_input_unit
    preconnected(2)%file_name = "stdout"
    preconnected(2)%unit_number = default_output_unit

  END SUBROUTINE init_preconnection_list

! *****************************************************************************
!> \brief Opens the requested file using a free unit number
!> \param file_name ...
!> \param file_status ...
!> \param file_form ...
!> \param file_action ...
!> \param file_position ...
!> \param file_pad ...
!> \param unit_number ...
!> \param debug ...
!> \param skip_get_unit_number ...
!> \author Matthias Krack (MK)
! *****************************************************************************
  SUBROUTINE open_file(file_name,file_status,file_form,file_action,&
                       file_position,file_pad,unit_number,debug,&
                       skip_get_unit_number)

    CHARACTER(LEN=*), INTENT(IN)             :: file_name
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: file_status, file_form, &
                                                file_action, file_position, &
                                                file_pad
    INTEGER, INTENT(INOUT)                   :: unit_number
    INTEGER, INTENT(IN), OPTIONAL            :: debug
    LOGICAL, INTENT(IN), OPTIONAL            :: skip_get_unit_number

    CHARACTER(LEN=*), PARAMETER :: routineN = 'open_file', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=11) :: access_string, action_string, current_action, &
      current_form, form_string, pad_string, position_string, status_string
    CHARACTER(LEN=2*default_path_length)     :: message
    CHARACTER(LEN=default_path_length)       :: real_file_name
    INTEGER                                  :: debug_unit, istat
    LOGICAL                                  :: exists, get_a_new_unit, &
                                                is_open

    access_string = "SEQUENTIAL"

    IF (PRESENT(file_status)) THEN
       status_string = TRIM(file_status)
    ELSE
       status_string = "OLD"
    END IF

    IF (PRESENT(file_form)) THEN
       form_string = TRIM(file_form)
    ELSE
       form_string = "FORMATTED"
    END IF

    IF (PRESENT(file_pad)) THEN
       pad_string = file_pad
       IF (form_string=="UNFORMATTED") THEN
          WRITE (UNIT=message,FMT="(A)")&
            "The PAD specifier is not allowed for an UNFORMATTED file"
          CALL stop_program(routineN,moduleN,__LINE__,TRIM(message))
       END IF
    ELSE
       pad_string = "YES"
    END IF

    IF (PRESENT(file_action)) THEN
       action_string = TRIM(file_action)
    ELSE
       action_string = "READ"
    END IF

    IF (PRESENT(file_position)) THEN
       position_string = TRIM(file_position)
    ELSE
       position_string = "REWIND"
    END IF

    IF (PRESENT(debug)) THEN
       debug_unit = debug
    ELSE
       debug_unit = 0 ! use default_output_unit for debugging
    END IF

    real_file_name = file_name
    IF (status_string == "OLD")&
       real_file_name = discover_file(file_name)

    ! Check the specified input file name
    INQUIRE (FILE=TRIM(real_file_name),EXIST=exists,OPENED=is_open,IOSTAT=istat)

    IF (istat /= 0) THEN
       WRITE (UNIT=message,FMT="(A,I0,A)")&
         "An error occurred inquiring the file <"//TRIM(real_file_name)//&
         "> (IOSTAT = ",istat,")"
       CALL stop_program(routineN,moduleN,__LINE__,TRIM(message))
    ELSE IF (status_string == "OLD") THEN
       IF (.NOT.exists) THEN
          WRITE (UNIT=message,FMT="(A)")&
            "The specified OLD file <"//TRIM(ADJUSTL(real_file_name))//&
            "> cannot be opened. It does not exist."
          CALL stop_program(routineN,moduleN,__LINE__,TRIM(message))
       END IF
    END IF

    ! Open the specified input file
    IF (is_open) THEN
       INQUIRE (FILE=TRIM(real_file_name),NUMBER=unit_number,&
                ACTION=current_action,FORM=current_form)
       IF (TRIM(position_string) == "REWIND") REWIND (UNIT=unit_number)
       IF (TRIM(status_string) == "NEW") THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "Attempt to re-open the existing OLD file <"//&
                            TRIM(real_file_name)//"> with status attribute NEW.")
       END IF
       IF (TRIM(current_form) /= TRIM(form_string)) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "Attempt to re-open the existing "//&
                            TRIM(current_form)//" file <"//TRIM(real_file_name)//&
                            "> as "//TRIM(form_string)//" file.")
       END IF
       IF (TRIM(current_action) /= TRIM(action_string)) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "Attempt to re-open the existing file <"//&
                            TRIM(real_file_name)//"> with the modified ACTION attribute "//&
                            TRIM(action_string)//". The current ACTION attribute is "//&
                            TRIM(current_action)//".")
       END IF
    ELSE
       ! Find an unused unit number
       get_a_new_unit = .TRUE.
       IF (PRESENT(skip_get_unit_number)) THEN
          IF (skip_get_unit_number) get_a_new_unit = .FALSE.
       END IF
       IF (get_a_new_unit) unit_number = get_unit_number(TRIM(real_file_name))
       IF (unit_number < 1) THEN
          WRITE (UNIT=message,FMT="(A)")&
            "Cannot open the file <"//TRIM(ADJUSTL(real_file_name))//&
            ">, because no unused logical unit number could be obtained."
          CALL stop_program(routineN,moduleN,__LINE__,TRIM(message))
       END IF
       IF (TRIM(form_string)== "FORMATTED") THEN
          OPEN (UNIT=unit_number,&
                FILE=TRIM(real_file_name),&
                STATUS=TRIM(status_string),&
                ACCESS=TRIM(access_string),&
                FORM=TRIM(form_string),&
                POSITION=TRIM(position_string),&
                ACTION=TRIM(action_string),&
                PAD=TRIM(pad_string),&
                IOSTAT=istat)
       ELSE
          OPEN (UNIT=unit_number,&
                FILE=TRIM(real_file_name),&
                STATUS=TRIM(status_string),&
                ACCESS=TRIM(access_string),&
                FORM=TRIM(form_string),&
                POSITION=TRIM(position_string),&
                ACTION=TRIM(action_string),&
                IOSTAT=istat)
       END IF
       IF (istat /= 0) THEN
          WRITE (UNIT=message,FMT="(A,I0,A,I0,A)")&
            "An error occurred opening the file <"//TRIM(ADJUSTL(real_file_name))//&
            "> with the unit number ",unit_number," (IOSTAT = ",istat,")"
          CALL stop_program(routineN,moduleN,__LINE__,TRIM(message))
       END IF
    END IF

    IF (debug_unit > 0) THEN
       INQUIRE (FILE=TRIM(real_file_name),OPENED=is_open,NUMBER=unit_number,&
                POSITION=position_string,NAME=message,ACCESS=access_string,&
                FORM=form_string,ACTION=action_string)
       WRITE (UNIT=debug_unit,FMT="(T2,A)") "BEGIN DEBUG "//TRIM(ADJUSTL(routineN))
       WRITE (UNIT=debug_unit,FMT="(T3,A,I0)") "NUMBER  : ",unit_number
       WRITE (UNIT=debug_unit,FMT="(T3,A,L1)") "OPENED  : ",is_open
       WRITE (UNIT=debug_unit,FMT="(T3,A)") "NAME    : "//TRIM(ADJUSTL(message))
       WRITE (UNIT=debug_unit,FMT="(T3,A)") "POSITION: "//TRIM(ADJUSTL(position_string))
       WRITE (UNIT=debug_unit,FMT="(T3,A)") "ACCESS  : "//TRIM(ADJUSTL(access_string))
       WRITE (UNIT=debug_unit,FMT="(T3,A)") "FORM    : "//TRIM(ADJUSTL(form_string))
       WRITE (UNIT=debug_unit,FMT="(T3,A)") "ACTION  : "//TRIM(ADJUSTL(action_string))
       WRITE (UNIT=debug_unit,FMT="(T2,A)") "END DEBUG "//TRIM(ADJUSTL(routineN))
       CALL print_preconnection_list()
    END IF

  END SUBROUTINE open_file

! *****************************************************************************
!> \brief Checks if file exists, considering also the file discovery mechanism.
!> \param file_name ...
!> \retval exist
!> \author Ole Schuett
! *****************************************************************************
  FUNCTION file_exists(file_name) RESULT(exist)
    CHARACTER(LEN=*), INTENT(IN)             :: file_name
    LOGICAL                                  :: exist

    CHARACTER(LEN=default_path_length)       :: real_file_name

    real_file_name = discover_file(file_name)
    INQUIRE(FILE=TRIM(real_file_name), exist=exist)
  END FUNCTION file_exists

! *****************************************************************************
!> \brief Checks various locations for a file name.
!> \param file_name ...
!> \retval real_file_name ...
!> \author Ole Schuett
! *****************************************************************************
  FUNCTION discover_file(file_name) RESULT(real_file_name)
    CHARACTER(LEN=*), INTENT(IN)             :: file_name
    CHARACTER(LEN=default_path_length)       :: real_file_name

    CHARACTER(LEN=default_path_length)       :: candidates(3), data_path
    INTEGER                                  :: i, stat
    LOGICAL                                  :: exists

    candidates(:) = file_name

    CALL GET_ENVIRONMENT_VARIABLE("CP2K_DATA_DIR", data_path, status=stat)
    IF(stat == 0) &
        candidates(2) = join_paths(data_path, file_name)

#if defined(__DATA_DIR)
     candidates(3) = join_paths(__DATA_DIR, file_name)
#endif

    real_file_name = file_name !fall back in case file not found

    DO i=1, SIZE(candidates)
       INQUIRE(file=TRIM(candidates(i)), exist=exists, iostat=stat)
       IF(stat==0 .AND. exists) THEN
           real_file_name = candidates(i)
           EXIT
       ENDIF
    ENDDO

  END FUNCTION discover_file

! *****************************************************************************
!> \brief Joins two file-paths, inserting '/' as needed.
!> \param path1 ...
!> \param path2 ...
!> \retval joined_path ...
!> \author Ole Schuett
! *****************************************************************************
  FUNCTION join_paths(path1, path2) RESULT(joined_path)
    CHARACTER(LEN=*), INTENT(IN)             :: path1, path2
    CHARACTER(LEN=default_path_length)       :: joined_path

    INTEGER                                  :: n

    n = LEN_TRIM(path1)
    IF(path2(1:1)=='/') THEN
        joined_path = path2
    ELSE IF(n==0 .OR. path1(n:n)=='/') THEN
        joined_path = TRIM(path1) // path2
    ELSE
        joined_path = TRIM(path1) // '/' // path2
    ENDIF
  END FUNCTION join_paths

! *****************************************************************************
!> \brief Print the list of preconnected units
!> \par History
!>      - Creation (22.02.2011,MK)
!> \author Matthias Krack (MK)
! *****************************************************************************
  SUBROUTINE print_preconnection_list()

    CHARACTER(LEN=*), PARAMETER :: routineN = 'print_preconnection_list', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ic, nc

    nc = SIZE(preconnected)

    WRITE (UNIT=default_output_unit,FMT="(A,/,A)")&
      " LIST OF PRECONNECTED LOGICAL UNITS",&
      "  Slot   Unit number   File name"
    DO ic=1,nc
       IF (preconnected(ic)%unit_number > 0) THEN
          WRITE (UNIT=default_output_unit,FMT="(I6,3X,I6,8X,A)")&
            ic,preconnected(ic)%unit_number,&
            TRIM(ADJUSTL(preconnected(ic)%file_name))
       ELSE
          WRITE (UNIT=default_output_unit,FMT="(I6,17X,A)")&
            ic,"UNUSED"
       END IF
    END DO

  END SUBROUTINE print_preconnection_list

END MODULE cp_files
