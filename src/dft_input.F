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
!!     Read in the DFT/WAVE specific input
!!
!!   AUTHOR
!!     JGH?
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!!
!******************************************************************************

MODULE dft_input
  
  USE dft_types, ONLY : dft_control_type, wave_control_type
  USE ewald_parameters_types, ONLY : ewald_parameters_type
  USE global_types, ONLY : global_environment_type
  USE input_types, ONLY : setup_parameters_type
  USE kinds, ONLY : dbl, YES, NO
  USE message_passing, ONLY : mp_bcast
  USE parser, ONLY : parser_init, parser_end, read_line, test_next, &
       cfield, p_error, get_real, get_int
  USE termination, ONLY : stop_program
  USE string_utilities, ONLY : uppercase, str_comp, xstring
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: read_dft_section, read_wave_section
  
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
!!     JGH?
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!!
!******************************************************************************
!!>----------------------------------------------------------------------------
!!  SECTION: &dft ... &end                                                  !
!!                                                                            !
!!  simulation       [md,debug]                                               !
!!  printlevel       globenv%print_level                                      !
!!  units            [kelvin,atomic]                                          !
!!  functional       [ ... ]                                                  !
!!  mesh             xc_mesh(3)                                               !
!!  spherical_grid   n                                                        !
!!      rad(i) ang(i) i=1..n                                                  !
!!  density_cutoff   den_cut                                                  !
!!  gradient_cutoff  grad_cut                                                 !
!!  set_file         "filename"                                               !
!!  input_file       "filename"                                               !
!!                                                                            !
!!<----------------------------------------------------------------------------
!!*****

SUBROUTINE read_dft_section ( setup, dft_control, dftpar )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( setup_parameters_type ), INTENT ( INOUT ) :: setup
  TYPE ( dft_control_type ), INTENT ( OUT ) :: dft_control
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: dftpar
  
! Locals
  CHARACTER ( LEN = 80 ) :: string, string2
  CHARACTER ( LEN = 6 ) :: label
  INTEGER :: iw, ierror, ilen, ios, i, n, source, group, ia, ie
  
!------------------------------------------------------------------------------
  
! defaults
  setup % run_type = "ENERGY"
  setup % unit_type = "ATOMIC"
  dft_control % functional = "PZ"
  dft_control % density_cut = 1.0E-10_dbl
  dft_control % gradient_cut = 1.0E-8_dbl
  dft_control % fixedspinmoment_amount = 0.0_dbl
  dft_control % broadening_tperform = NO
  dft_control % broadening_method = "NONE"
  dft_control % broadening_width = 0.0_dbl
  dft_control % nstates = -1
  dft_control % nspins = 1
  dft_control % charge = 0.0_dbl
  dft_control % nelectrons = 0
  
  iw = dftpar % scr
  
! parse the input section
  label = "&DFT"
  CALL parser_init ( dftpar % input_file_name, label, ierror, dftpar )
  IF ( ierror /= 0 ) THEN
     IF ( dftpar % ionode ) &
          WRITE ( iw, '( A )' ) &
          " read_dft_section| No input section &DFT found "
  ELSE
     CALL read_line
     
     DO WHILE ( test_next() /= "X" )
        
        ilen = 40
        CALL cfield ( string, ilen )
        CALL uppercase ( string )
        
        SELECT CASE ( string )
        CASE DEFAULT
           CALL p_error()
           CALL stop_program ( "read_dft_section", "unknown option" // string )
           
        CASE ( "SIMULATION_TYPE" )
           ilen = 40
           CALL cfield ( setup % run_type, ilen )
           CALL uppercase ( setup % run_type )
           
        CASE ( "PRINTLEVEL" )
           dftpar % print_level = get_int()
           
        CASE ( "UNITS" )
           ilen = 40
           CALL cfield ( setup % unit_type, ilen )
           CALL uppercase ( setup % unit_type )
           
        CASE ( "FUNCTIONAL" )
           ilen = 3
           CALL cfield ( string2, ilen )
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
           
!*apsi* 140201: Are these necessary?
!        CASE ( "SPHERICAL_GRID", "SPHERICAL_GRIDS" )
!           n = get_int()
!           dft_control % spherical_grids = n
!           ALLOCATE ( dft_control % acs_rad ( n ), STAT = ios )
!           IF ( ios /= 0 ) CALL stop_program ( "read_dft_section", &
!                "failed to allocate dft_control%acs_rad" )
!           ALLOCATE ( dft_control % acs_ang ( n ), STAT = ios )
!           IF ( ios /= 0 ) CALL stop_program ( "read_dft_section", &
!                "failed to allocate dft_control%acs_ang" )
!           DO i = 1, n
!              dft_control % acs_rad ( i ) = get_int()
!              dft_control % acs_ang ( i ) = get_int()
!           END DO
           
        CASE ( "DENSITY_CUTOFF" )
           dft_control % density_cut = get_real()
        CASE ( "GRADIENT_CUTOFF" )
           dft_control % gradient_cut = get_real()
           
        CASE ( "SET_FILE" )
           ilen = 132
           CALL cfield ( setup % set_file_name, ilen )
           
        CASE ( "INPUT_FILE" )
           ilen = 132
           CALL cfield ( setup % input_file_name, ilen )
           
        CASE ( "STATES", "NSTATES" )
           dft_control % nstates = get_int()
           
        CASE ( "SPIN", "SPINS", "NSPINS" )
           dft_control % nspins = get_int()
           
        CASE ( "SPINMOMENT", "FIXEDSPINMOMENT" )
           dft_control % fixedspinmoment_amount = get_real()
           
        CASE ( "CHARGE" )
           dft_control % charge = get_real()
           
        CASE ( "BROADENING" )
           dft_control % broadening_tperform = YES
           
           ilen = 40
           CALL cfield ( string2, ilen )
           CALL uppercase ( string2 )
           
           dft_control % broadening_method = string2
           dft_control % broadening_width = get_real()
           
        END SELECT
        
        CALL read_line
     END DO
     
  END IF
  CALL parser_end()
! end of parsing the input section
  
! write some information to output
  IF ( dftpar % ionode .AND. dftpar % print_level >= 0 ) THEN
     WRITE ( iw, '()' )
     WRITE ( iw, '( A, T41, A )' ) &
          " DFT| Functional", ADJUSTR ( dft_control % functional )
     
     WRITE ( iw, '( A, T68, E13.5 )' ) &
          " DFT| Cutoffs : density ", dft_control % density_cut
     WRITE ( iw, '( A, T68, E13.5 )' ) &
          " DFT|           gradient", dft_control % gradient_cut
     
!*apsi* 140201: Are these necessary?
!     IF ( dft_control % xc_mesh ( 1 ) == -1 ) then
!        WRITE ( iw, '( A, T74, A )' ) " DFT| XC mesh ", "default"
!     ELSE
!        WRITE ( iw, '( A, T51, 3I10 )' ) &
!             " DFT| XC mesh ", dft_control % xc_mesh
!     END IF
!     WRITE ( iw, '( A, T77, I4 )' ) " DFT| Number of spherical grids ", &
!          dft_control % spherical_grids
!     DO i = 1, dft_control % spherical_grids
!        WRITE ( iw, '( A, I4, T35, A, I6, T61, A, I6 )' ) &
!             " DFT| Grid #:", i, &
!             " radial pts. ", dft_control % acs_rad ( i ), &
!             " angular pts. ", dft_control % acs_ang ( i )
!     END DO
  END IF
  
!*apsi* 140201: Are these necessary?
!! update defaults
!  IF ( dft_control % spherical_grids == -1 ) THEN
!     dft_control % spherical_grids = 1
!     ALLOCATE ( dft_control % acs_rad ( 1 ), STAT = ios )
!     IF ( ios /= 0 ) CALL stop_program ( "read_dft_section", &
!          "failed to allocate dft_control % acs_rad" )
!     ALLOCATE ( dft_control % acs_ang ( 1 ), STAT = ios )
!     IF ( ios /= 0 ) CALL stop_program ( "read_dft_section", &
!          "failed to allocate dft_control%acs_ang" )
!     dft_control % acs_rad ( 1 ) = 50
!     dft_control % acs_ang ( 1 ) = 96
!  END IF
  
  IF ( dft_control % nstates < 0 ) THEN
     CALL stop_program ( "read_dft_section", &
          "please supply a value for nstates" )
  END IF
  IF ( dft_control % nspins /= 1 .AND. dft_control % nspins /= 2 ) THEN
     CALL stop_program ( "read_dft_section", &
          "please supply a value 1 or 2 for nspins" )
  END IF
  
END SUBROUTINE read_dft_section

!******************************************************************************

SUBROUTINE read_wave_section ( setup, wave_control, wavepar )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( setup_parameters_type ), INTENT ( INOUT ) :: setup
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: wavepar
  TYPE ( wave_control_type ), INTENT ( OUT ) :: wave_control
  
! Locals
  CHARACTER ( LEN = 20 ) :: string, string2
  CHARACTER ( LEN = 6 ) :: label
  INTEGER :: iw, ierror, ilen, ios, i, n, source, group, ia, ie
  
!------------------------------------------------------------------------------
  
! defaults
  wave_control % xc_mesh ( : ) = -1
  wave_control % wf_optimisation_method = "SD"
  wave_control % e_cutoff_wf = -1.0_dbl
  wave_control % e_cutoff_dual = 4.0_dbl
  
  iw = wavepar % scr
  
! parse the input section
  label = "&WAVE"
  CALL parser_init ( wavepar % input_file_name, label, ierror, wavepar )
  IF ( ierror /= 0 ) THEN
     WRITE ( iw, '( A )' ) " No input section &WAVE found"
  ELSE
     CALL read_line
     DO WHILE ( test_next() /= "X" )
        ilen = 40
        CALL cfield ( string, ilen )
        CALL uppercase ( string )
        
        SELECT CASE ( string )
        CASE DEFAULT
           CALL p_error()
           CALL stop_program ( "read_wave_section", "unknown option" // string )
           
        CASE ( "PRINTLEVEL" )
           wavepar % print_level = get_int()
           
        CASE ( "UNITS" )
           ilen = 40
           CALL cfield ( setup % unit_type, ilen )
           CALL uppercase ( setup % unit_type )
           
        CASE ( "MESH" )
           wave_control % xc_mesh ( 1 ) = get_int()
           wave_control % xc_mesh ( 2 ) = get_int()
           wave_control % xc_mesh ( 3 ) = get_int()
           
        CASE ( "SET_FILE" )
           ilen = 132
           CALL cfield ( setup % set_file_name, ilen )
           
        CASE ( "INPUT_FILE" )
           ilen = 132
           CALL cfield ( setup % input_file_name, ilen )
           
        CASE ( "ECUT", "ECUTOFF", "CUTOFF" )
           wave_control % e_cutoff_wf = get_real()
           
        CASE ( "DUAL" )
           wave_control % e_cutoff_dual = get_real()
           
        CASE ( "PCG" )
           wave_control % wf_optimisation_method = "PCG"
           
        CASE ( "DIIS" )
           wave_control % wf_optimisation_method = "DIIS"
           
        CASE ( "ELECTRONIC ITERATIONS", "ITERATIONS ELECTRONS" )
           wave_control % niter_elc = get_int()
           
        CASE ( "PRINT FREQUENCY" )
           wave_control % wf_print_frequency = get_int()
           
        CASE ( "TIMESTEP" )
           wave_control % wf_timestep = get_real()
           
        CASE ( "EMASS" )
           wave_control % wf_emass = get_real()
           
        CASE ( "SUBSPACE ROTATION" )
           wave_control % wf_subspacerotation_tperform = YES
           wave_control % wf_subspacerotation_frequency = get_int()
           wave_control % wf_subspacerotation_first = get_int()
           
        END SELECT
        
        CALL read_line
     END DO
     
  END IF
  CALL parser_end()
! end of parsing the input section
  
! write some information to output
  IF ( wavepar % ionode .AND. wavepar % print_level >= 0 ) THEN
     WRITE ( iw, '()' )
     
     IF ( wave_control % xc_mesh ( 1 ) == -1 ) then
        WRITE ( iw, '( A, T74, A )' ) " WAVE| XC mesh ", "default"
     ELSE
        WRITE ( iw, '( A, T51, 3I10 )' ) &
             " WAVE| XC mesh ", wave_control % xc_mesh
     END IF
     
     WRITE ( iw, '()' )
  END IF
  
  IF ( wave_control % e_cutoff_wf < 0.0_dbl ) THEN
     CALL stop_program ( "read_wave_section", &
          "please set a cutoff energy for the wave functions" )
  END IF
  
END SUBROUTINE read_wave_section

!******************************************************************************

END MODULE dft_input
