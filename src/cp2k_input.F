!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/cp2k_input [1.0] *
!!
!!   NAME
!!     cp2k_input
!!
!!   FUNCTION
!!     Master input
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE cp2k_input
  
  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl
  USE message_passing, ONLY : mp_bcast
  USE parser, ONLY : parser_init, parser_end, read_line, test_next, &
       cfield, p_error, get_real, get_int
  USE stop_program, ONLY : stop_prg
  USE string_utilities, ONLY : uppercase, xstring
  USE fft_tools, ONLY : init_fft
  
  PRIVATE
  PUBLIC :: read_cp2k_section

!!*****
!******************************************************************************
  
CONTAINS

!******************************************************************************
!!****** cp2k_input/read_cp2k_section [1.0] *
!!
!!   NAME
!!     read_cp2k_section
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     2-Dec-2000 (JGH) added default fft library
!!
!!   INPUTS
!!I---------------------------------------------------------------------------I
!!I                                                                           I
!!I SECTION: &cp2k ... &end                                                   I
!!I                                                                           I
!!I program prg_name                                                          I
!!I project name                                                              I
!!I iolevel n                                                                 I
!!I fftlib  default_fft_library                                               I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE read_cp2k_section ( globenv )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
  
! Locals
  INTEGER :: ierror, ilen, iw, source, allgrp, i1, i2
  CHARACTER ( LEN = 6 ) :: string
  CHARACTER ( LEN = 5 ) :: label
  
! ..defaults
  globenv%program_name = 'NONE'
  globenv%project_name = 'project'
  globenv%print_level = 0
  globenv%default_fft_library = 'FFTSG'
  
   iw = globenv % scr
! ..parse the input section
   label = '&CP2K'
   ierror = 0
   CALL parser_init(globenv%input_file_name,label,ierror,globenv)
   IF (ierror /= 0) THEN

      IF (globenv % ionode) &
        WRITE ( iw, '( a )' ) ' No input section &CP2K found '

   ELSE
        
      CALL read_line()
      DO WHILE ( test_next() /= 'X' )
         ilen = 6
         CALL cfield ( string, ilen )
         CALL uppercase ( string )
         
         SELECT CASE (string)
         CASE DEFAULT
            CALL p_error()
            CALL stop_prg( 'input_cntl','unknown option')
            
         CASE ( 'PROGRA' )
            ilen = 0
            CALL cfield ( globenv % program_name, ilen )
            CALL uppercase ( globenv % program_name )
            
         CASE ( 'PROJEC' )
            ilen = 0
            CALL cfield ( globenv % project_name, ilen )
            
         CASE ( 'IOLEVE' )
            globenv%print_level = get_int()

         CASE ( 'FFTLIB' )
            ilen = 0
            CALL cfield ( globenv % default_fft_library, ilen )
            CALL uppercase ( globenv % default_fft_library )

         END SELECT
         
         CALL read_line()
      END DO
      
   END IF
   CALL parser_end
! ..end of parsing the input section
     
   IF ( globenv % program_name == 'QUICKSTEP' ) globenv % program_name = 'QS'

! initialize FFT library
   CALL xstring ( globenv % default_fft_library, i1, i2 )
   CALL init_fft ( fftlib = globenv % default_fft_library ( i1:i2 ) )

     
   IF ( globenv % ionode ) THEN
! ..write some information to output
     WRITE ( iw, '( A, T41, A )' ) ' CP2K| Program name ', &
          ADJUSTR ( globenv % program_name ( 1:40 ) )
     WRITE ( iw, '( A, T41, A )' ) ' CP2K| Input file name ', &
          ADJUSTR ( globenv % input_file_name ( 1:40 ) )
     WRITE ( iw, '( A, T41, A )' ) ' CP2K| Project name ', &
          ADJUSTR ( globenv % project_name ( 1:40 ) )
     WRITE ( iw, '( A, T41, A )' ) ' CP2K| Default FFT library ', &
          ADJUSTR ( globenv % default_fft_library ( 1:40 ) )
     WRITE ( iw, '( A, T77, I4 )' ) ' CP2K| Global print level ', &
          globenv % print_level
     WRITE ( iw, '( A, T75, I6 )' ) ' CP2K| Total number of processors ', &
          globenv % num_pe
     WRITE ( iw, '( A, T75, I6 )' ) ' CP2K| This output from processor ', &
          globenv % mepos
     WRITE (iw,*)
  END IF
  
END SUBROUTINE read_cp2k_section

!!*****
!******************************************************************************

END MODULE cp2k_input

!******************************************************************************
