!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/wave [1.0] *
!!
!!   NAME
!!     wave
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH, 18-NOV-2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE wave_module

  USE band, ONLY : band_structure_type
  USE brillouin, ONLY : kpoint_type, brillouin_info, kpoint_input
  USE cntl_input, ONLY : read_cntl_section
  USE cp_control_types, ONLY : dft_control_type, read_dft_control, read_wave_control, &
                        write_dft_control, write_wave_control
  USE extended_system_types, ONLY : extended_system_type
  USE fermi, ONLY : fermi_distribution_type
  USE force_fields, ONLY : read_force_field_section
  USE global_types, ONLY : global_environment_type
  USE header, ONLY : wave_header
  USE input_types, ONLY : setup_parameters_type
  USE kinds, ONLY : dbl
  USE mathconstants, ONLY : zero, twopi
  USE md, ONLY : simulation_parameters_type
  USE molecule_input, ONLY : read_molecule_section, read_setup_section
  USE molecule_types, ONLY : molecule_type, intra_parameters_type
  USE particle_types, ONLY : particle_prop_type
  USE read_pimd, ONLY : read_pimd_section
  USE structure_types, ONLY : structure_type, interaction_type
  USE termination, ONLY : stop_program, stop_memory
  USE timings, ONLY : timeset, timestop, trace_debug

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: wave
  
!!*****
!-----------------------------------------------------------------------------!

CONTAINS

!-----------------------------------------------------------------------------!
! CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL  !
!-----------------------------------------------------------------------------!
!!****** wave_module/wave [1.0] *
!!
!!   NAME
!!     wave
!!
!!   FUNCTION
!!     Controls program flow for WAVE
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE wave ( globenv )

  IMPLICIT NONE

! Argument
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv

! Locals
  INTEGER :: handle, beads, isos
  
  TYPE ( molecule_type ), DIMENSION ( : ), POINTER :: mol_setup
  TYPE ( simulation_parameters_type ) :: simpar
  TYPE ( structure_type ), ALLOCATABLE, DIMENSION ( : ) :: struc
  TYPE ( setup_parameters_type ) :: setup
  TYPE ( kpoint_type ) :: kp
  TYPE ( fermi_distribution_type ) :: fd
  TYPE ( band_structure_type ) :: bs
  TYPE ( dft_control_type ), POINTER :: dft_control
  TYPE ( particle_prop_type ), DIMENSION ( : ), POINTER :: pstat
  TYPE ( interaction_type ) :: inter
  TYPE ( intra_parameters_type ) :: intra_param

!------------------------------------------------------------------------------

! uncomment the next line for a traced executable
! IF( globenv % ionode ) CALL trace_debug ( "START" )

  CALL timeset ( 'WAVE', 'I', ' ', handle )

  CALL wave_header ( globenv % scr )

! read control section
  CALL read_cntl_section ( setup, globenv )

! path integrals; special input section
  IF ( setup % path_integrals ) THEN
    CALL read_pimd_section ( simpar % pimd_params, globenv )
    beads = simpar % pimd_params % beads
  ELSE
    beads = 1
    simpar % pimd_params % beads = 0
  END IF

! the structure that holds all the system information
  ALLOCATE ( struc ( beads ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'allocstruc', 'struc', beads )

! read the dft and wave section inputs
  CALL read_dft_control ( dft_control, globenv )
  CALL write_dft_control ( dft_control, globenv )
  CALL read_wave_control ( dft_control % wave_control, globenv )
  CALL write_wave_control ( dft_control % wave_control, globenv )

! read the K-Point information
  CALL kpoint_input ( kp, globenv )
  CALL brillouin_info ( kp, globenv )

! read from the setup and molecule section
  CALL read_setup_section ( mol_setup, setup, globenv )
  CALL read_molecule_section ( mol_setup, setup, globenv )

! read possible force fields
  CALL read_force_field_section ( setup, mol_setup, setup%set_file_name,&
       intra_param, inter%potparm, pstat, globenv )
        
  CALL timestop ( zero, handle )
  
END SUBROUTINE wave


!-----------------------------------------------------------------------------!
! CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL  !
!-----------------------------------------------------------------------------!

!******************************************************************************

SUBROUTINE set_var ( struc, extended_type, beads )
  IMPLICIT NONE
! Arguments
  TYPE ( structure_type ), DIMENSION ( : ), INTENT ( INOUT ) :: struc
  TYPE ( extended_system_type ), INTENT ( INOUT ) :: extended_type
  INTEGER, INTENT ( IN ) :: beads
! Local
  INTEGER :: isos, nll, i, ibead

  ALLOCATE ( extended_type % nhc_part ( beads ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'allocstruc', 'nhc_part', beads )

  ALLOCATE ( extended_type % nhc_baro ( beads ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'allocstruc', 'nhc_baro', beads )

! total number of neighbor linked lists needed
  nll = 1

  DO ibead = 1, beads
    ALLOCATE ( struc ( ibead ) % ll_data ( nll ), STAT = isos )
    IF ( isos /=0 ) CALL stop_memory ( 'control', 'll_data', 0 )
    DO i = 1, nll
      struc ( ibead ) % ll_data ( i ) % natom_types = 0
      struc ( ibead ) % ll_data ( i ) % list_type = 0
      struc ( ibead ) % ll_data ( i ) % counter = 0
      struc ( ibead ) % ll_data ( i ) % last_update = 0
      struc ( ibead ) % ll_data ( i ) % num_update = 0
      struc ( ibead ) % ll_data ( i ) % print_level = 1
      NULLIFY ( struc ( ibead ) % ll_data ( i ) % r_last_update )
      NULLIFY ( struc ( ibead ) % ll_data ( i ) % rlist_cut )
      NULLIFY ( struc ( ibead ) % ll_data ( i ) % rlist_cutsq )
      NULLIFY ( struc ( ibead ) % ll_data ( i ) % pp_images )
      NULLIFY ( struc ( ibead ) % ll_data ( i ) % pp_ncell )
      NULLIFY ( struc ( ibead ) % ll_data ( i ) % pp_startlist_im )
      NULLIFY ( struc ( ibead ) % ll_data ( i ) % pp_startlist_nl )
      NULLIFY ( struc ( ibead ) % ll_data ( i ) % pp_neighbor )
      NULLIFY ( struc ( ibead ) % ll_data ( i ) % pp_startlist_cell )
      NULLIFY ( struc ( ibead ) % ll_data ( i ) % pp_cell_ll )
    END DO
  END DO

END SUBROUTINE set_var

!******************************************************************************

END MODULE wave_module

!******************************************************************************
