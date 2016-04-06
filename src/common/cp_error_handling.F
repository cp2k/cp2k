!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Module that contains the routines for error handling
!> \author Ole Schuett
! **************************************************************************************************
MODULE cp_error_handling
   USE base_hooks,                      ONLY: cp_abort_hook,&
                                              cp_warn_hook
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: default_output_unit,&
                                              m_flush,&
                                              m_walltime
   USE message_passing,                 ONLY: mp_abort
   USE termination,                     ONLY: print_message
   USE timings,                         ONLY: print_stack

!$ USE OMP_LIB,                        ONLY: omp_get_thread_num

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_error_handling'

   !API public routines
   PUBLIC :: cp_error_handling_setup

   !API (via pointer assignment to hook, PR67982, not meant to be called directly)
   PUBLIC :: cp_abort_handler, cp_warn_handler

   INTEGER, PUBLIC, SAVE :: warning_counter = 0

CONTAINS

! **************************************************************************************************
!> \brief Registers handlers with base_hooks.F
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE cp_error_handling_setup()
      cp_abort_hook => cp_abort_handler
      cp_warn_hook => cp_warn_handler
   END SUBROUTINE cp_error_handling_setup

! **************************************************************************************************
!> \brief Abort program with error message
!> \param location ...
!> \param message ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE cp_abort_handler(location, message)
      CHARACTER(len=*), INTENT(in)                       :: location, message

      INTEGER                                            :: unit_nr

      CALL delay_non_master() ! cleaner output if all ranks abort simultaneously

      unit_nr = cp_logger_get_default_io_unit()
      IF (unit_nr <= 0) &
         unit_nr = default_output_unit ! fall back to stdout

      CALL print_abort_message(message, location, unit_nr)
      CALL print_stack(unit_nr)
      CALL m_flush(unit_nr)

      CALL mp_abort()
   END SUBROUTINE cp_abort_handler

! **************************************************************************************************
!> \brief Signal a warning
!> \param location ...
!> \param message ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE cp_warn_handler(location, message)
      CHARACTER(len=*), INTENT(in)                       :: location, message

      INTEGER                                            :: unit_nr

!$OMP MASTER
      warning_counter = warning_counter+1
!$OMP END MASTER

      unit_nr = cp_logger_get_default_io_unit()
      IF (unit_nr > 0) THEN
         CALL print_message("WARNING in "//TRIM(location)//' :: '//TRIM(ADJUSTL(message)), unit_nr, 1, 1, 1)
         CALL m_flush(unit_nr)
      ENDIF
   END SUBROUTINE cp_warn_handler

! **************************************************************************************************
!> \brief Delay non-master ranks/threads, used by cp_abort_handler()
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE delay_non_master()
      INTEGER                                            :: unit_nr
      REAL(KIND=dp)                                      :: t1, wait_time

      wait_time = 0.0_dp

      ! we (ab)use the logger to determine the first MPI rank
      unit_nr = cp_logger_get_default_io_unit()
      IF (unit_nr <= 0) &
         wait_time = wait_time+1.0_dp ! rank-0 gets a head start of one second.

!$    IF (omp_get_thread_num() /= 0) &
!$       wait_time = wait_time+1.0_dp ! master threads gets another second

      ! sleep
      IF (wait_time > 0.0_dp) THEN
         t1 = m_walltime()
         DO
            IF (m_walltime()-t1 > wait_time .OR. t1 < 0) EXIT
         ENDDO
      ENDIF

   END SUBROUTINE delay_non_master

! **************************************************************************************************
!> \brief Prints a nicely formatted abort message box
!> \param message ...
!> \param location ...
!> \param output_unit ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE print_abort_message(message, location, output_unit)
      CHARACTER(LEN=*), INTENT(IN)                       :: message, location
      INTEGER, INTENT(IN)                                :: output_unit

      CHARACTER(LEN=*), PARAMETER :: img = "   ___   "//"  /   \  "//" [ABORT] "//"  \___/  "// &
         "    |    "//"  O/|    "//" /| |    "//" / \     "
      INTEGER, PARAMETER :: img_height = 8, img_width = 9, screen_width = 80, &
         txt_width = screen_width-img_width-5

      CHARACTER(LEN=screen_width)                        :: msg_line
      INTEGER                                            :: a, b, c, fill, i, img_start, indent, &
                                                            msg_height, msg_start

! count message lines

      a = 1; b = -1; msg_height = 0
      DO WHILE (b < LEN_TRIM(message))
         b = next_linebreak(message, a, txt_width)
         a = b+1
         msg_height = msg_height+1
      ENDDO

      ! calculate message and image starting lines
      IF (img_height > msg_height) THEN
         msg_start = (img_height-msg_height)/2+1
         img_start = 1
      ELSE
         msg_start = 1
         img_start = msg_height-img_height+2
      ENDIF

      ! print empty line
      WRITE (UNIT=output_unit, FMT="(A)") ""

      ! print opening line
      WRITE (UNIT=output_unit, FMT="(T2,A)") REPEAT("*", screen_width-1)

      ! print body
      a = 1; b = -1; c = 1
      DO i = 1, MAX(img_height-1, msg_height)
         WRITE (UNIT=output_unit, FMT="(A)", advance='no') " *"
         IF (i < img_start) THEN
            WRITE (UNIT=output_unit, FMT="(A)", advance='no') REPEAT(" ", img_width)
         ELSE
            WRITE (UNIT=output_unit, FMT="(A)", advance='no') img(c:c+img_width-1)
            c = c+img_width
         ENDIF
         IF (i < msg_start) THEN
            WRITE (UNIT=output_unit, FMT="(A)", advance='no') REPEAT(" ", txt_width+2)
         ELSE
            b = next_linebreak(message, a, txt_width)
            msg_line = message(a:b)
            a = b+1
            fill = (txt_width-LEN_TRIM(msg_line))/2+1
            indent = txt_width-LEN_TRIM(msg_line)-fill+2
            WRITE (UNIT=output_unit, FMT="(A)", advance='no') REPEAT(" ", indent)
            WRITE (UNIT=output_unit, FMT="(A)", advance='no') TRIM(msg_line)
            WRITE (UNIT=output_unit, FMT="(A)", advance='no') REPEAT(" ", fill)
         ENDIF
         WRITE (UNIT=output_unit, FMT="(A)", advance='yes') "*"
      ENDDO

      ! print location line
      WRITE (UNIT=output_unit, FMT="(A)", advance='no') " *"
      WRITE (UNIT=output_unit, FMT="(A)", advance='no') img(c:c+img_width-1)
      indent = txt_width-LEN_TRIM(location)+1
      WRITE (UNIT=output_unit, FMT="(A)", advance='no') REPEAT(" ", indent)
      WRITE (UNIT=output_unit, FMT="(A)", advance='no') TRIM(location)
      WRITE (UNIT=output_unit, FMT="(A)", advance='yes') " *"

      ! print closing line
      WRITE (UNIT=output_unit, FMT="(T2,A)") REPEAT("*", screen_width-1)

      ! print empty line
      WRITE (UNIT=output_unit, FMT="(A)") ""

   END SUBROUTINE print_abort_message

! **************************************************************************************************
!> \brief Helper routine for print_abort_message()
!> \param message ...
!> \param pos ...
!> \param rowlen ...
!> \retval ibreak ...
!> \author Ole Schuett
! **************************************************************************************************
   FUNCTION next_linebreak(message, pos, rowlen) RESULT(ibreak)
      CHARACTER(LEN=*), INTENT(IN)                       :: message
      INTEGER, INTENT(IN)                                :: pos, rowlen
      INTEGER                                            :: ibreak

      INTEGER                                            :: i, n

      n = LEN_TRIM(message)
      IF (n-pos <= rowlen) THEN
         ibreak = n ! remaining message shorter than line
      ELSE
         i = INDEX(message(pos+1:pos+1+rowlen), " ", BACK=.TRUE.)
         IF (i == 0) THEN
            ibreak = pos+rowlen-1 ! no space found, break mid-word
         ELSE
            ibreak = pos+i ! break at space closest to rowlen
         ENDIF
      ENDIF
   END FUNCTION next_linebreak

END MODULE cp_error_handling
