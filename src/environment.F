!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/environment [1.0] *
!!
!!   NAME
!!     environment
!!
!!   FUNCTION
!!     Sets up and terminates the global environment vaiables
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   USED BY
!!     cp2k
!!
!!   SOURCE
!******************************************************************************

MODULE environment
  
  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl
  USE machine, ONLY : m_hostnm, m_getcwd, m_getlog, &
       m_getuid, m_getpid, m_getarg
  USE mathconstants, ONLY : zero
  USE physcon, ONLY : init_physcon
  USE stop_program, ONLY : set_stop, stop_prg
  USE output_utilities, ONLY : set_output
  USE timesl, ONLY : walltime, cputime, datum
  USE timings, ONLY : timeprint, timeset, timestop
  USE util, ONLY : ran2
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: initialisation, trailer
  
  CHARACTER ( LEN = 32 ) :: hname = ' ' !! name of the computer
  CHARACTER ( LEN = 60 ) :: curdir = ' ' !! current directory name
  CHARACTER ( LEN = 10 ) :: user = ' ' !! user name
  INTEGER :: my_pid = 0 !! process ID
  INTEGER :: my_uid = 0 !! user ID
  INTEGER :: handle

!!*****
!******************************************************************************
  
CONTAINS

!******************************************************************************
!!****** environment/initialisation [1.0] *
!!
!!   NAME
!!     initialisation
!!
!!   FUNCTION
!!     Initializes global environment vaiables
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE initialisation ( globenv )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv ! check intent
  
! Locals
  CHARACTER ( LEN = 26 ) :: datx
  INTEGER :: iw, l, initial_random, i1, i2
  
!------------------------------------------------------------------------------
  
  globenv % scr = 6  
  
  CALL set_stop ( globenv % scr )
  CALL set_output ( globenv % scr, globenv % ionode )
  
  IF ( globenv % ionode ) THEN
     iw = globenv % scr
     
     CALL set_environment ( globenv )
     
     CALL datum ( datx )
     
     WRITE ( iw, &
          '( A, T55, A26, /, A, T49, A32, /, A, T71, A10, /, A, T71, I10 )' ) &
          '  **** **** ******  **  PROGRAM STARTED AT ', ADJUSTR ( datx ), &
          ' ***** ** ***  *** **   PROGRAM STARTED ON ', ADJUSTR ( hname ), &
          ' **    ****   ******    PROGRAM STARTED BY ', ADJUSTR ( user ), &
          ' ***** **    ** ** **   PROGRAM PROCESS ID ', my_pid
     
     l = LEN_TRIM ( curdir )
     IF ( l <= 38 ) THEN
        WRITE ( iw, '( A, T43, A38 )' ) &
             '  **** **  *******  **  PROGRAM STARTED IN ', &
             ADJUSTR ( curdir ( 1 : l ) )
     ELSE
        WRITE ( iw, '( A )' ) &
             '  **** **  *******  **  PROGRAM STARTED IN '
        WRITE ( iw, '( T21, A60 )' ) ADJUSTR ( curdir )
     END IF
     
     WRITE ( iw, '( )' )
  END IF
  
! initialize timing 
  CALL timeset ( 'CP2K', 'I', ' ', handle )
  
! initialize random number generator
  globenv % idum = -1
  initial_random = ran2 ( globenv % idum )
  
! initialize physical constants
  CALL init_physcon()

END SUBROUTINE initialisation

!!*****
!******************************************************************************
!!****** environment/trailer [1.0] *
!!
!!   NAME
!!     trailer
!!
!!   FUNCTION
!!     Writes final timings and banner
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE trailer ( globenv )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv ! check intent
  
! Locals
  CHARACTER ( LEN = 26 ) :: datx
  INTEGER :: iw, l
  
!------------------------------------------------------------------------------
  
  CALL timestop ( zero, handle )
  
  IF ( globenv % ionode ) THEN
     iw = globenv % scr
     
     CALL timeprint ( iw, globenv % print_level )
     
     WRITE ( iw, '( )' )
     
     CALL datum ( datx )
     WRITE ( iw, &
          '( A, T55, A26, /, A, T49, A32, /, A, T71, A10, /, A, T71, I10 )' ) &
          '  **** **** ******  **  PROGRAM ENDED AT   ', ADJUSTR ( datx ), &
          ' ***** ** ***  *** **   PROGRAM RAN ON     ', ADJUSTR ( hname ), &
          ' **    ****   ******    PROGRAM RAN BY     ', ADJUSTR ( user ), &
          ' ***** **    ** ** **   PROGRAM PROCESS ID ', my_pid
     
     l = LEN_TRIM ( curdir )
     IF ( l <= 38 ) THEN
        WRITE ( iw, '( A, T43, A38 )' ) &
             '  **** **  *******  **  PROGRAM RAN IN     ', &
             ADJUSTR ( curdir ( 1 : l ) )
     ELSE
        WRITE ( iw, '( A )' ) &
             '  **** **  *******  **  PROGRAM RAN IN     '
        WRITE ( iw, '( T21, A60 )' ) ADJUSTR ( curdir )
     END IF
  END IF
  
END SUBROUTINE trailer

!!*****
!******************************************************************************
!!****** environment/set_environment [1.0] *
!!
!!   NAME
!!     set_environment
!!
!!   FUNCTION
!!     Calls for UNIX environment variabls
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE set_environment ( globenv )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
  
! Locals
  INTEGER :: narg
  
  INTEGER :: iargc
  
!------------------------------------------------------------------------------
  
  CALL m_hostnm ( hname )
  CALL m_getcwd ( curdir )
  CALL m_getlog ( user )
  CALL m_getuid ( my_uid )
  CALL m_getpid ( my_pid )
  
  narg = iargc()
  
  SELECT CASE ( narg )
  CASE ( 0 )
     globenv % input_file_name = 'input'
  CASE ( 1 )
     CALL m_getarg ( 1, globenv % input_file_name )
  CASE DEFAULT
     CALL stop_prg ( 'set_environment', &
          'illegal number of command line arguments passed' )
  END SELECT
  
END SUBROUTINE set_environment

!!*****
!******************************************************************************

END MODULE environment

!******************************************************************************
