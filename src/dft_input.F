!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/dft_input [1.0] *
!!
!!   NAME
!!     dft_input
!!
!!   FUNCTION
!!     Read in the DFT specific input
!!     Read additional input for the DFT programs
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     JGH (21.11.2001) : add Quickstep routine
!!                        clean up and embrace method section from QS
!!
!!   SOURCE
!!
!******************************************************************************

MODULE dft_input
  
  USE dft_types, ONLY : dft_control_type, wave_control_type
  USE ewald_parameters_types, ONLY : ewald_parameters_type
  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl, YES, NO
  USE message_passing, ONLY : mp_bcast
  USE parser, ONLY : parser_init, parser_end, read_line, test_next, &
       p_error, get_next, search_label
  USE termination, ONLY : stop_program
  USE string_utilities, ONLY : uppercase, str_comp, xstring
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: read_dft_section, read_wave_section, read_qs_section
  
!!*****

CONTAINS

!!****** dft_input/read_dft_section [1.0] *
!!
!!   NAME
!!     read_dft_section
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!!
!******************************************************************************
!!>----------------------------------------------------------------------------
!!  SECTION: &dft ... &end                                                    !
!!                                                                            !
!!  functional       [ ... ]                                                  !
!!  density_cutoff   den_cut                                                  !
!!  gradient_cutoff  grad_cut                                                 !
!!                                                                            !
!!<----------------------------------------------------------------------------

SUBROUTINE read_dft_section ( dft_control, dftpar )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( dft_control_type ), INTENT ( OUT ) :: dft_control
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: dftpar
  
! Locals
  CHARACTER ( LEN = 80 ) :: string, string2
  CHARACTER ( LEN = 6 ) :: label
  INTEGER :: iw, ierror, ilen, ios, i, n, source, group, ia, ie
  
!------------------------------------------------------------------------------
  
! defaults
  dft_control % functional = "PZ"
  dft_control % density_cut = 1.0E-10_dbl
  dft_control % gradient_cut = 1.0E-8_dbl
  dft_control % nspins = 1
  dft_control % xc_mesh ( : ) = -1
  
  iw = dftpar % scr
  
! parse the input section
  label = "&DFT"
  CALL parser_init ( dftpar % input_file_name, dftpar )
  CALL search_label ( label, ierror, ignore_case=.TRUE. )
  IF ( ierror /= 0 ) THEN
     IF ( dftpar % ionode ) &
          WRITE ( iw, '( A )' ) &
          " read_dft_section| No input section &DFT found "
  ELSE
     CALL read_line
     
     DO WHILE ( test_next() /= "X" )
        
        ilen = 40
        CALL get_next ( string, ilen )
        CALL uppercase ( string )
        
        SELECT CASE ( string )
        CASE DEFAULT
           CALL p_error()
           CALL stop_program ( "read_dft_section", "unknown option" // string )
           
        CASE ( "FUNCTIONAL" )
           ilen = 8
           CALL get_next ( string2, ilen )
           CALL uppercase ( string2 )
           
           SELECT CASE ( string2 )
           CASE DEFAULT
              CALL p_error()
              CALL stop_program ( "read_dft_section", &
                   "unknown suboption to FUNCTIONAL" )
           CASE ( "LDA", "PZ" )
              dft_control % functional = "PZ"
           CASE ( "BP", "BECKE-PERDEW" )
              dft_control % functional = "BP"
           CASE ( "BLYP", "BECKE-LEE-YANG-PARR" )
              dft_control % functional = "BLYP"
           CASE ( "PBE" )
              dft_control % functional = "PBE"
           END SELECT
           
        CASE ( "DENSITY_CUTOFF" )
           CALL get_next ( dft_control % density_cut )

        CASE ( "GRADIENT_CUTOFF" )
           CALL get_next ( dft_control % gradient_cut )
           
        CASE ( "MESH" )
           CALL get_next ( dft_control % xc_mesh ( 1 ) )
           CALL get_next ( dft_control % xc_mesh ( 2 ) )
           CALL get_next ( dft_control % xc_mesh ( 3 ) )
           
        CASE ( "SPIN", "SPINS", "NSPINS" )
           CALL get_next ( dft_control % nspins )
           
        CASE ( "LSD", "SPIN POLARISATION", "SPIN POLARIZATION" )
           dft_control % nspins = 2
           
        END SELECT
        
        CALL read_line
     END DO
     
  END IF
  CALL parser_end()
! end of parsing the input section
  
  IF ( dft_control % nspins < 0 ) THEN
     CALL stop_program ( "read_dft_section", &
          "please supply a valid value for nspins" )
  END IF

! write some information to output
  IF ( dftpar % ionode .AND. dftpar % print_level >= 0 ) THEN
     WRITE ( iw, '( A, T41, A )' ) &
          " DFT| Functional", ADJUSTR ( dft_control % functional )
     IF ( dft_control % nspins == 1 ) then
        !do nothing
     ELSEIF ( dft_control % nspins == 2 ) then
        WRITE ( iw, '( A, T60, A )' ) " DFT| ","spin polarized system"
     ELSE
        WRITE ( iw, '( A, T60, A, T77, I4 )' ) " DFT| ","sets of spin states:",&
           dft_control % nspins
     END IF
     WRITE ( iw, '( A, T68, E13.5 )' ) &
          " DFT| Cutoffs : density ", dft_control % density_cut
     WRITE ( iw, '( A, T68, E13.5 )' ) &
          " DFT|           gradient", dft_control % gradient_cut
     IF ( dft_control % xc_mesh ( 1 ) == -1 ) then
        WRITE ( iw, '( A, T74, A )' ) " DFT| XC mesh (PW density)", "default"
     ELSE
        WRITE ( iw, '( A, T51, 3I10 )' ) &
             " DFT| XC mesh (PW density)", dft_control % xc_mesh
     END IF
     WRITE ( iw, '()' )
  END IF
  
END SUBROUTINE read_dft_section

!!*****
!******************************************************************************
!!****** dft_input/read_wave_section [1.0] *
!!
!!   NAME
!!     read_wave_section
!!
!!   FUNCTION
!!     Reads the section wave from the input file
!!
!!   AUTHOR
!!     APS
!!
!!   MODIFICATION HISTORY
!!     JGH (18/11/2001) : clean and expand
!!
!!   SOURCE
!!
!******************************************************************************
!!>----------------------------------------------------------------------------
!!  SECTION: &wave ... &end                                                    !
!!                                                                            !
!!  ecut             e_cutoff_wf                                              !
!!  dual             e_cutoff_dual                                            !
!!                                                                            !
!!<----------------------------------------------------------------------------

SUBROUTINE read_wave_section ( wave_control, wavepar )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: wavepar
  TYPE ( wave_control_type ), INTENT ( OUT ) :: wave_control
  
! Locals
  CHARACTER ( LEN = 20 ) :: string, string2
  CHARACTER ( LEN = 6 ) :: label
  INTEGER :: iw, ierror, ilen, ios, i, n, source, group, ia, ie
  
!------------------------------------------------------------------------------
  
! defaults
  wave_control % e_cutoff_wf = -1.0_dbl
  wave_control % e_cutoff_dual = 4.0_dbl
  
  iw = wavepar % scr
  
! parse the input section
  label = "&WAVE"
  CALL parser_init ( wavepar % input_file_name, wavepar )
  CALL search_label ( label, ierror, ignore_case=.TRUE. )
  IF ( ierror /= 0 ) THEN
     WRITE ( iw, '( A )' ) " No input section &WAVE found"
  ELSE
     CALL read_line
     DO WHILE ( test_next() /= "X" )
        ilen = 40
        CALL get_next ( string, ilen )
        CALL uppercase ( string )
        
        SELECT CASE ( string )
        CASE DEFAULT
           CALL p_error()
           CALL stop_program ( "read_wave_section", "unknown option" // string )
           
        CASE ( "ECUT", "ECUTOFF", "CUTOFF" )
           CALL get_next ( wave_control % e_cutoff_wf )
           
        CASE ( "DUAL" )
           CALL get_next ( wave_control % e_cutoff_dual )
           
        END SELECT
        
        CALL read_line
     END DO
     
  END IF
  CALL parser_end()
! end of parsing the input section
  
! write some information to output
  IF ( wavepar % ionode .AND. wavepar % print_level >= 0 ) THEN

     WRITE ( iw, '( A, T71, F10.3 )' ) " WAVE| Plane wave cutoff [a.u.]", &
           wave_control % e_cutoff_wf
           
     WRITE ( iw, '( A, T71, F10.2 )' ) " WAVE| Density cutoff (DUAL) ", &
           wave_control % e_cutoff_dual
     
     WRITE ( iw, '()' )
  END IF
  
  IF ( wave_control % e_cutoff_wf < 0.0_dbl ) THEN
     CALL stop_program ( "read_wave_section", &
          "please set a cutoff energy for the wave functions" )
  END IF
  
END SUBROUTINE read_wave_section

!!*****
!******************************************************************************
!!****** dft_input/read_qs_section [1.0] *
!!
!!   NAME
!!     read_qs_section
!!
!!   FUNCTION
!!     Reads the section QS from the input file
!!
!!   AUTHOR
!!     JGH (21.11.2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!!
!******************************************************************************
!!>----------------------------------------------------------------------------
!!  SECTION: &qs ... &end                                                     !
!!                                                                            !
!!  method           qs_method                                                !
!!                                                                            !
!!<----------------------------------------------------------------------------

SUBROUTINE read_qs_section ( qs_control, qspar )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: qspar
  TYPE ( wave_control_type ), INTENT ( OUT ) :: qs_control
  
! Locals
  CHARACTER ( LEN = 20 ) :: string, string2
  CHARACTER ( LEN = 6 ) :: label
  INTEGER :: iw, ierror, ilen, ios, i, n, source, group, ia, ie
  
!------------------------------------------------------------------------------
  
! defaults
! qs_control % method = -1
  
  iw = qspar % scr
  
! parse the input section
  label = "&QS"
  CALL parser_init ( qspar % input_file_name, qspar )
  CALL search_label ( label, ierror, ignore_case=.TRUE. )
  IF ( ierror /= 0 ) THEN
     WRITE ( iw, '( A )' ) " No input section &QS found"
  ELSE
     CALL read_line
     DO WHILE ( test_next() /= "X" )
        ilen = 40
        CALL get_next ( string, ilen )
        CALL uppercase ( string )
        
        SELECT CASE ( string )
        CASE DEFAULT
           CALL p_error()
           CALL stop_program ( "read_qs_section", "unknown option" // string )
           
        CASE ( "MESH" )
           
        END SELECT
        
        CALL read_line
     END DO
     
  END IF
  CALL parser_end()
! end of parsing the input section
  
! write some information to output
  IF ( qspar % ionode .AND. qspar % print_level >= 0 ) THEN
     
!    WRITE ( iw, '( A, T71, F10.3 )' ) " WAVE| Plane wave cutoff [a.u.]", &
!          wave_control % e_cutoff_wf
           
     WRITE ( iw, '()' )
  END IF
  
END SUBROUTINE read_qs_section

!!*****
!******************************************************************************

END MODULE dft_input

!******************************************************************************
