!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/control_module [1.0] *
!!
!!   NAME
!!     control_module
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH, 14-NOV-2000
!!
!!   MODIFICATION HISTORY
!!     CJM, 20-Feb-2001:  Now contains logic for restart options
!!     'INIT' and 'ALL'. Also deallocates new ammendments to
!!      system_type. Passes box_ref to ewald routines for
!!      consistent restarts
!!   SOURCE
!******************************************************************************

MODULE control_module

  USE atoms_input, ONLY : read_coord_vel, system_type
  USE band, ONLY : band_structure_type, init_band_structure
  USE brillouin, ONLY : kpoint_type, brillouin_info, kpoint_input
  USE cntl_input, ONLY : read_cntl_section
  USE convert_units, ONLY : convert
  USE dft_input, ONLY : read_dft_section, read_wave_section
  USE dft_types, ONLY : dft_control_type, wave_control_type
  USE dump, ONLY : dump_variables
  USE ewalds, ONLY : ewald_print, ewald_correction
  USE ewald_parameters_types, ONLY : ewald_parameters_type
  USE fermi, ONLY : fermi_distribution_type, init_fermi_dist, fermi_info
  USE fist_debug, ONLY : fist_debug_control => debug_control
  USE force_fields, ONLY : read_force_field_section, ATOMNAMESLENGTH
  USE global_types, ONLY : global_environment_type
  USE header, ONLY : fist_header, tbmd_header, wave_header
  USE initialize_extended_types, ONLY : initialize_extended_system, &
                                        assign_extended_parameters
  USE initialize_molecule_types, ONLY : initialize_molecule_type
  USE initialize_particle_types, ONLY : initialize_particle_type
  USE input_types, ONLY : setup_parameters_type
  USE integrator, ONLY : velocity_verlet, force, set_energy_parm, energy, &
       set_integrator
  USE kinds, ONLY : dbl
  USE kpoint_initialization, ONLY : initialize_kpoints
  USE linklist_control, ONLY : set_ll_parm
  USE mathconstants, ONLY : zero
  USE md, ONLY : read_md_section, simulation_parameters_type, &
       initialize_velocities, thermodynamic_type, mdio_parameters_type
  USE molecule_input, ONLY : read_molecule_section, read_setup_section, &
       charge
  USE molecule_types, ONLY : molecule_type, intra_parameters_type
  USE nose, ONLY : lnhc_parameters_type, npt_info_type
  USE pair_potential, ONLY : spline_nonbond_control
  USE particle_types, ONLY : particle_prop_type, particle_type
  USE simulation_cell, ONLY : cell_type, get_hinv
  USE termination, ONLY : stop_program, stop_memory
  USE structure_types, ONLY : structure_type, interaction_type
  USE tbmd_debug, ONLY : tbmd_debug_control => debug_control
  USE tbmd_initialize, ONLY : tbmd_init, tb_get_numel
  USE tbmd_input, ONLY : read_tb_hamiltonian, read_tb_hopping_elements
  USE timings, ONLY : timeset, timestop, trace_debug
  USE unit, ONLY : unit_convert_type, set_units
  USE util, ONLY : close_unit, get_share

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: control
  
!!*****
!-----------------------------------------------------------------------------!

CONTAINS

!-----------------------------------------------------------------------------!
! CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL  !
!-----------------------------------------------------------------------------!
!!****** control_module/control [1.0] *
!!
!!   NAME
!!     control
!!
!!   FUNCTION
!!     Controls program flow for FIST and TBMD
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     Harald Forbert (Dec-2000): enable multiple linked lists
!!
!!   USED BY
!!     cp2k
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE control ( globenv )

  IMPLICIT NONE

! Argument
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv

! Locals
  REAL ( dbl ), DIMENSION ( :, : ), ALLOCATABLE :: rcut
  INTEGER :: handle1, handle2, isos
  INTEGER :: nel
  CHARACTER ( LEN = ATOMNAMESLENGTH ), DIMENSION ( : ), POINTER :: atom_names
  CHARACTER ( LEN = 40 ) :: set_fn
  LOGICAL :: tbmd, fist, dft, wave, quickstep
  
  TYPE ( particle_prop_type ), DIMENSION ( : ), POINTER :: pstat
  TYPE ( molecule_type ), DIMENSION ( : ), POINTER :: mol_setup
  TYPE ( unit_convert_type ) :: units
  TYPE ( simulation_parameters_type ) :: simpar
  TYPE ( structure_type ) :: struc
  TYPE ( interaction_type ) :: inter
  TYPE ( lnhc_parameters_type ) :: nhc_part, nhc_baro
  TYPE ( npt_info_type ), DIMENSION ( :, : ), POINTER :: npt_info
  TYPE ( system_type ) :: ainp
  TYPE ( ewald_parameters_type ) :: ewald_param
  TYPE ( setup_parameters_type ) :: setup
  TYPE ( intra_parameters_type ) :: intra_param
  TYPE ( kpoint_type ) :: kp
  TYPE ( fermi_distribution_type ) :: fd
  TYPE ( band_structure_type ) :: bs
  TYPE ( dft_control_type ) :: dft_control
  TYPE ( mdio_parameters_type ) :: mdio
  
!------------------------------------------------------------------------------

! IF( globenv % ionode ) CALL trace_debug ( "START" )

  CALL timeset ( 'CONTROL', 'I', ' ', handle1 )
  CALL timeset ( 'CNTL_INIT', 'I', ' ', handle2 )

  struc % ll_data % natom_types = 0
  struc % ll_data % list_type = 0
  struc % ll_data % counter = 0
  struc % ll_data % last_update = 0
  struc % ll_data % num_update = 0
  struc % ll_data % print_level = 1
  NULLIFY ( struc % ll_data % r_last_update )
  NULLIFY ( struc % ll_data % rlist_cut )
  NULLIFY ( struc % ll_data % rlist_cutsq )
  NULLIFY ( struc % ll_data % pp_images )
  NULLIFY ( struc % ll_data % pp_ncell )
  NULLIFY ( struc % ll_data % pp_startlist_im )
  NULLIFY ( struc % ll_data % pp_startlist_nl )
  NULLIFY ( struc % ll_data % pp_neighbor )
  NULLIFY ( struc % ll_data % pp_startlist_cell )
  NULLIFY ( struc % ll_data % pp_cell_ll )
  
!*apsi* 130201: Waiting for Matthias + Jürg :*)
  quickstep = .FALSE.
  
  IF ( globenv % program_name == 'FIST' ) THEN
    fist = .true.
    tbmd = .false.
    wave = .false.
  ELSE IF ( globenv % program_name == 'TBMD' ) THEN
    tbmd = .true.
    fist = .false.
    wave = .false.
  ELSE IF ( globenv % program_name == 'WAVE' ) THEN
    tbmd = .false.
    fist = .false.
    wave = .true.
  ELSE
    CALL stop_program ( ' control ',' program_name not specified' )
  END IF
  
  dft = wave .OR. quickstep
  
  IF ( globenv % ionode ) THEN
    IF ( fist ) CALL fist_header ( globenv % scr )
    IF ( tbmd ) CALL tbmd_header ( globenv % scr )
    IF ( wave ) CALL wave_header ( globenv % scr )
  END IF
  
! read control section
  CALL read_cntl_section ( setup, ewald_param, globenv )

! read from the setup and molecule section of *.set

  set_fn = setup % set_file_name

  CALL read_setup_section ( mol_setup, set_fn, globenv )

  CALL read_molecule_section ( mol_setup, set_fn, globenv )

! read force_field information for classical MD
! read pair potential information for TB

  CALL read_force_field_section ( setup, mol_setup, set_fn, &
       intra_param, inter%potparm, atom_names, pstat, globenv )

!..read Hamiltonian section

  IF ( tbmd ) THEN
     CALL read_tb_hamiltonian ( setup, atom_names, inter%tbatom, globenv )
     CALL read_tb_hopping_elements ( setup, atom_names, &
          inter%tbatom, inter%tbhop, globenv )
  END IF
  
  IF ( dft ) THEN
     
     IF ( wave ) THEN
        
        CALL read_dft_section ( setup, dft_control, globenv )
        
        CALL read_wave_section ( setup, dft_control % wave_control, globenv )
        
     ELSE
        CALL stop_program ( "control", "DFT method not implemented" )
     END IF
     
  END IF
  
! read the input of the molecular dynamics section
  CALL read_md_section ( simpar, globenv, mdio )
  simpar % program = globenv % program_name

!..read atomic coordinates, velocities (optional) and the simulation box
  ainp % rtype = simpar % read_type
  CALL read_coord_vel ( ainp, setup % input_file_name, globenv )

!..initialize box, perd
  struc % box % hmat = ainp % box
  struc % box % perd = setup % perd

  IF ( simpar % read_type /= 'INIT' ) THEN
    struc % box_ref % hmat = ainp % box_ref
  ENDIF

! If run is a debug.  Make box_ref = box to work under
! all ensembles and restart options
  IF ( setup % run_type == 'DEBUG' ) THEN
    struc % box_ref % hmat = ainp % box
  ENDIF

!..K-Points
  IF ( SUM ( struc % box % perd ) == 0 ) THEN
    kp % scheme = "NULL"
  ELSE IF ( tbmd .OR. wave ) THEN
    CALL kpoint_input ( kp, globenv )
  ELSE
    kp % scheme = "NULL"
  END IF
  
! initialize working units
  CALL set_units ( setup % unit_type, units )

  CALL set_energy_parm ( units % pconv, units % econv, units % l_label, &
       units % vol_label, units % e_label, units % pres_label, &
       units % temp_label, units % angl_label )

!..allocate memory for atoms and molecules
  CALL allocmem ( ainp, mol_setup, struc, globenv )

!..initialize particle_type
  CALL initialize_particle_type ( atom_names, simpar, mol_setup, &
       ainp, pstat, struc % part )

!..convert the units ( includes part and box ONLY when read_type is 'INIT' )
  IF ( simpar % read_type == 'INIT' ) THEN
    CALL convert ( units = units, simpar = simpar, &
                 part = struc % part, pstat = pstat, box = struc % box, &
                 potparm = inter%potparm, intra_param = intra_param, &
                 ewald_param = ewald_param )
  ELSE
    CALL convert ( units = units, simpar = simpar, &
                 pstat = pstat, potparm = inter%potparm, &
                 intra_param = intra_param, ewald_param = ewald_param )
  ENDIF

!..calculate the inverse box matrix now after the unit conversion,
!  so it also has the right units and assign the reference box
  CALL get_hinv ( struc % box )
  IF ( simpar % read_type == 'INIT' ) THEN
    struc % box_ref = struc % box
  ENDIF
  CALL get_hinv ( struc % box_ref )


!..initialize molecule_type
  CALL initialize_molecule_type ( mol_setup, intra_param, struc % pnode, &
       struc % part, struc % molecule, globenv )

!..allocate lnhc_parameters_type and npt_info and get
!  number of degrees of freedom
  CALL initialize_extended_system ( struc % box, simpar, &
     struc % molecule, mol_setup, globenv, nhc_part, nhc_baro, npt_info )

! initialize velocities if needed
  IF ( simpar % read_type == 'POS' .OR. simpar % read_type == 'INIT' ) THEN
     CALL initialize_velocities ( simpar, struc % part, globenv )
  END IF

! Assign extended system variables
  IF ( simpar % read_type == 'ALL' ) THEN
     CALL assign_extended_parameters ( struc % box % deth, simpar % ensemble, &
                                      ainp, npt_info, nhc_part, nhc_baro )
  ENDIF

!...initialize splines

  inter%potparm ( :, : ) % energy_cutoff = 0.0_dbl
  inter%potparm ( :, : ) % e_cutoff_coul = 0.0_dbl
  CALL spline_nonbond_control ( inter%potparm, pstat, 5000, ewald_param )

!..set linklist control parameters
  ALLOCATE ( rcut ( setup % natom_type, setup % natom_type ), STAT = isos )
  IF ( isos /=0 ) CALL stop_memory ( 'control', 'rcut', 0 )

  rcut ( :, : ) = inter%potparm ( :, : ) % rcutsq

  CALL set_ll_parm ( struc % ll_data, globenv, simpar % verlet_skin, &
       setup % natom_type, rcut, simpar % n_cell )

  CALL set_ll_parm ( struc % ll_data, globenv, printlevel = globenv % print_level, &
       ltype = 'NONBOND' )

  DEALLOCATE ( rcut, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'control', 'rcut' )

!..deallocate arrays needed for atom input
  IF ( ASSOCIATED ( ainp % c ) ) THEN
     DEALLOCATE ( ainp % c, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'control', 'ainp%c' )
  END IF

  IF ( ASSOCIATED ( ainp % v ) ) THEN
     DEALLOCATE ( ainp % v, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'control', 'ainp%v' )
  END IF

  IF ( ASSOCIATED ( ainp % eta_part ) ) THEN
     DEALLOCATE ( ainp % eta_part, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'control', 'ainp%eta_part' )
  END IF

  IF ( ASSOCIATED ( ainp % veta_part ) ) THEN
     DEALLOCATE ( ainp % veta_part, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'control', 'ainp%veta_part' )
  END IF

  IF ( ASSOCIATED ( ainp % mass_eta_part ) ) THEN
     DEALLOCATE ( ainp % mass_eta_part, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'control', 'ainp%mass_eta_part' )
  END IF

  IF ( ASSOCIATED ( ainp % eta_baro ) ) THEN
     DEALLOCATE ( ainp % eta_baro, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'control', 'ainp%eta_baro' )
  END IF

  IF ( ASSOCIATED ( ainp % veta_baro ) ) THEN
     DEALLOCATE ( ainp % veta_baro, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'control', 'ainp%veta_baro' )
  END IF

  IF ( ASSOCIATED ( ainp % mass_eta_baro ) ) THEN
     DEALLOCATE ( ainp % mass_eta_baro, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'control', 'ainp%mass_eta_baro' )
  END IF

  IF ( ASSOCIATED ( ainp % veps ) ) THEN
     DEALLOCATE ( ainp % veps, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'control', 'ainp%veps' )
  END IF

  IF ( ASSOCIATED ( ainp % mass_eps ) ) THEN
     DEALLOCATE ( ainp % mass_eps, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'control', 'ainp%mass_eps' )
  END IF
!
!..initialize the on-site terms for TB
  IF ( tbmd ) CALL tbmd_init ( struc % part, inter%tbatom, inter%tbhop )
  
!
! Symmetry and K-points
!
  IF ( SUM ( struc % box % perd ) == 0 ) THEN
!!!!!! symmetry setup for molecules
  ELSE IF ( tbmd ) THEN
    CALL initialize_kpoints ( globenv, kp, setup % symmetry, &
         struc % box % hmat, struc % part )
    IF ( globenv%ionode .AND. globenv%print_level > 0) &
       CALL brillouin_info ( kp, globenv%scr )
  END IF
  
! nota bene: we use the defaults for the elctron temperature and
!             spin polarisation
  IF ( tbmd ) then
    nel = tb_get_numel ( inter % tbatom, charge )
    CALL init_fermi_dist ( fd, nel )
    IF ( globenv % ionode ) &
         CALL fermi_info ( fd, globenv % print_level, globenv % scr )
    CALL init_band_structure ( bs, fd, kp )
  END IF
  
  CALL timestop ( zero, handle2 )
  
  CALL control_work ( globenv, units, simpar, struc, inter, nhc_part, &
       nhc_baro, npt_info, ewald_param, dft_control, setup, mdio, tbmd, fist, &
       kp )
  
! deallocate memory for atoms and molecules
  CALL deallocmem ( struc )
  
  CALL timestop ( zero, handle1 )
  
END SUBROUTINE control

!!*****
!-----------------------------------------------------------------------------!
! CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL  !
!-----------------------------------------------------------------------------!

!******************************************************************************

SUBROUTINE control_work ( globenv, units, simpar, struc, inter, nhc_part, &
     nhc_baro, npt_info, ewald_param, dft_control, setup, mdio, tbmd, fist, &
     kp )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
  TYPE ( unit_convert_type ) :: units
  TYPE ( simulation_parameters_type ) :: simpar
  TYPE ( structure_type ) :: struc
  TYPE ( interaction_type ) :: inter
  TYPE ( lnhc_parameters_type ) :: nhc_part, nhc_baro
  TYPE ( npt_info_type ), DIMENSION ( :, : ), POINTER :: npt_info
  TYPE ( ewald_parameters_type ) :: ewald_param
  TYPE ( setup_parameters_type ) :: setup
  TYPE ( mdio_parameters_type ) :: mdio
  TYPE ( dft_control_type ) :: dft_control
  LOGICAL, INTENT ( IN ) :: tbmd, fist
  TYPE ( kpoint_type ), INTENT ( IN ) :: kp
  
! Locals
  INTEGER :: handle
  
!------------------------------------------------------------------------------
  
  CALL timeset ( 'CNTL_WORK', 'I', ' ', handle )
  
  SELECT CASE ( setup % run_type )
  CASE DEFAULT
     CALL stop_program ( "control_work", "no suitable run_type" &
                         //  setup % run_type )
     
  CASE ( "DEBUG" )
     
! debug the forces
     CALL control_debug_work ( globenv, simpar, struc, inter, ewald_param, &
          mdio, tbmd, fist )
     
  CASE ( "ENERGY" )
     
     CALL control_energy_work ( globenv, units, simpar, struc, inter, &
          nhc_part, nhc_baro, npt_info, ewald_param, mdio, dft_control, kp )
     
  CASE ( "MD" )
     
     CALL control_md_work ( globenv, units, simpar, struc, inter, nhc_part, &
          nhc_baro, npt_info, ewald_param, mdio, dft_control, kp )
     
  END SELECT
  
  CALL timestop ( zero, handle )
  
END SUBROUTINE control_work

!******************************************************************************

SUBROUTINE control_debug_work ( globenv, simpar, struc, inter, ewald_param, &
     mdio, tbmd, fist )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
  TYPE ( simulation_parameters_type ) :: simpar
  TYPE ( structure_type ) :: struc
  TYPE ( interaction_type ) :: inter
  TYPE ( ewald_parameters_type ) :: ewald_param
  TYPE ( mdio_parameters_type ) :: mdio
  LOGICAL :: tbmd, fist
  
! Locals
  TYPE ( thermodynamic_type ) :: thermo
  
!------------------------------------------------------------------------------
  
! initialize integrator
  CALL set_integrator ( globenv, mdio )
  
  IF ( fist ) then
     CALL fist_debug_control ( globenv, ewald_param, struc%part, &
          struc%pnode, struc%molecule, struc%box, struc%box_ref, thermo, &
          inter%potparm, simpar % ensemble, struc % ll_data )
  ELSE IF ( tbmd ) THEN
     CALL tbmd_debug_control ( globenv, ewald_param, struc%part, &
          struc%pnode, struc%molecule, struc%box, struc%box_ref, thermo, &
          inter%potparm, simpar % ensemble, struc % ll_data )
  ELSE
     CALL stop_program ( "control_debug_work", "no method to debug found" )
  END IF
  
END SUBROUTINE control_debug_work

!******************************************************************************

SUBROUTINE control_energy_work ( globenv, units, simpar, struc, inter, &
     nhc_part, nhc_baro, npt_info, ewald_param, mdio, dft_control, kp )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
  TYPE ( unit_convert_type ) :: units
  TYPE ( simulation_parameters_type ) :: simpar
  TYPE ( structure_type ) :: struc
  TYPE ( interaction_type ) :: inter
  TYPE ( lnhc_parameters_type ) :: nhc_part, nhc_baro
  TYPE ( npt_info_type ), DIMENSION ( :, : ), POINTER :: npt_info
  TYPE ( ewald_parameters_type ) :: ewald_param
  TYPE ( mdio_parameters_type ) :: mdio
  TYPE ( dft_control_type ) :: dft_control
  TYPE ( kpoint_type ), INTENT ( IN ) :: kp
  
! Locals
  REAL ( dbl ) :: cons
  INTEGER :: itimes
  
  TYPE ( thermodynamic_type ) :: thermo
  
!------------------------------------------------------------------------------
  
! initialize integrator
  CALL set_integrator ( globenv, mdio )
  
! MD
  itimes = 0
  
  CALL force ( struc, inter, thermo, simpar, ewald_param, globenv, &
       dft_control, kp )
  IF ( globenv % ionode .AND. ewald_param % ewald_type /= 'NONE' ) &
       CALL ewald_print ( globenv % scr, thermo, struc % box, &
       units % e_label )
  
  CALL energy ( itimes, cons, simpar, struc, thermo, nhc_part, &
       nhc_baro, npt_info )
  
  CALL velocity_verlet ( itimes, cons, simpar, inter, thermo, &
       struc, ewald_param, nhc_part, nhc_baro, npt_info, dft_control, kp )
  
  CALL dump_variables ( struc, nhc_part, nhc_baro, npt_info, &
                        mdio % dump_file_name, globenv )
  
  IF ( globenv % ionode ) CALL close_unit ( 10, 99 )
  
END SUBROUTINE control_energy_work

!******************************************************************************

SUBROUTINE control_md_work ( globenv, units, simpar, struc, inter, nhc_part, &
     nhc_baro, npt_info, ewald_param, mdio, dft_control, kp )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
  TYPE ( unit_convert_type ) :: units
  TYPE ( simulation_parameters_type ) :: simpar
  TYPE ( structure_type ) :: struc
  TYPE ( interaction_type ) :: inter
  TYPE ( lnhc_parameters_type ) :: nhc_part, nhc_baro
  TYPE ( npt_info_type ), DIMENSION ( :, : ), POINTER :: npt_info
  TYPE ( ewald_parameters_type ) :: ewald_param
  TYPE ( mdio_parameters_type ) :: mdio
  TYPE ( dft_control_type ) :: dft_control
  TYPE ( kpoint_type ), INTENT ( IN ) :: kp
  
! Locals
  REAL ( dbl ) :: cons
  INTEGER :: itimes
  
  TYPE ( thermodynamic_type ) :: thermo
  
!------------------------------------------------------------------------------
  
! initialize integrator
  CALL set_integrator ( globenv, mdio )
  
! MD
  itimes = 0
  
  CALL force ( struc, inter, thermo, simpar, ewald_param, globenv, &
       dft_control, kp )
  IF ( globenv % ionode .AND. ewald_param % ewald_type /= 'NONE' ) &
       CALL ewald_print ( globenv % scr, thermo, struc % box, &
       units % e_label )
  
  CALL energy ( itimes, cons, simpar, struc, thermo, nhc_part, &
       nhc_baro, npt_info )
  
  DO itimes = 1, simpar % nsteps
     CALL velocity_verlet ( itimes, cons, simpar, inter, thermo, &
          struc, ewald_param, nhc_part, nhc_baro, npt_info, dft_control, kp )
     
     IF ( MOD ( itimes, mdio % idump ) == 0 ) &
          CALL dump_variables ( struc, nhc_part, nhc_baro, npt_info, &
                                mdio % dump_file_name, globenv )
  END DO
  
  CALL dump_variables ( struc, nhc_part, nhc_baro, npt_info, &
                        mdio % dump_file_name, globenv )
  
  IF ( globenv % ionode ) CALL close_unit ( 10, 99 )
  
END SUBROUTINE control_md_work

!******************************************************************************

#if defined ( __JUST_A_SECURITY_COPY )
SUBROUTINE control_work ( globenv, units, simpar, struc, inter, nhc_part, &
     nhc_baro, npt_info, ewald_param, setup )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
  TYPE ( unit_convert_type ) :: units
  TYPE ( simulation_parameters_type ) :: simpar
  TYPE ( structure_type ) :: struc
  TYPE ( interaction_type ) :: inter
  TYPE ( lnhc_parameters_type ) :: nhc_part, nhc_baro
  TYPE ( npt_info_type ), DIMENSION ( :, : ), POINTER :: npt_info
  TYPE ( ewald_parameters_type ) :: ewald_param
  TYPE ( setup_parameters_type ) :: setup
  LOGICAL :: tbmd, fist
  
! Locals
  REAL ( dbl ) :: cons, ecut
  INTEGER :: handle, itimes
  
  TYPE ( thermodynamic_type ) :: thermo
  
!------------------------------------------------------------------------------
  
  CALL timeset ( 'CNTL_WORK', 'I', ' ', handle )
  
  IF ( setup % run_type == 'DEBUG' ) THEN

! debug the forces
    
! initialize integrator
    CALL set_integrator ( globenv, mdio )
    
    IF ( fist ) then
      CALL fist_debug_control( globenv, ewald_param, struc%part, struc%pnode, &
                 struc%molecule, struc%box, struc%box_ref, thermo,  &
                 inter%potparm, simpar % ensemble, struc % ll_data )
    ELSE IF ( tbmd ) THEN
      CALL tbmd_debug_control( globenv, ewald_param, struc%part, struc%pnode, &
                 struc%molecule, struc%box, struc%box_ref, thermo,  &
                 inter%potparm, simpar % ensemble, struc % ll_data )
    ENDIF

  ELSE

!..initialize integrator
    CALL set_integrator ( globenv, mdio )

!..MD
    itimes = 0

    CALL force ( struc, inter, thermo, simpar, ewald_param, globenv )
    IF ( globenv % ionode .AND. ewald_param % ewald_type /= 'NONE' ) &
         CALL ewald_print ( globenv % scr, thermo, struc % box, &
         units % e_label )

    CALL energy ( itimes, cons, simpar, struc, thermo, nhc_part, &
                    nhc_baro, npt_info )
    
    DO itimes = 1, simpar % nsteps
       CALL velocity_verlet ( itimes, cons, simpar, inter, thermo, &
            struc, ewald_param, nhc_part, nhc_baro, npt_info )

       IF ( MOD ( itimes, mdio % idump ) == 0 ) &
            CALL dump_variables ( struc, nhc_part, nhc_baro, &
                           npt_info, mdio % dump_file_name, globenv )
    END DO

    CALL dump_variables ( struc, nhc_part, nhc_baro, npt_info, &
                          mdio % dump_file_name, globenv )

    IF ( globenv % ionode ) CALL close_unit ( 10, 99 )
  END IF

!..deallocate memory for atoms and molecules
  CALL deallocmem ( struc )
  
  CALL timestop ( zero, handle )
  
END SUBROUTINE control_work
#endif

!******************************************************************************

SUBROUTINE allocmem ( ainp, mol_setup, struc, globenv )

  IMPLICIT NONE

! Arguments
  TYPE ( molecule_type ), DIMENSION ( : ), INTENT ( IN ) :: mol_setup
  TYPE ( structure_type ), INTENT ( INOUT ) :: struc
  TYPE ( system_type ), INTENT ( IN ) :: ainp
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv

! Locals
  INTEGER :: iw, natoms, nnodes, nmol, nmoltype, ios, iat, i, nsh

!------------------------------------------------------------------------------

  struc % name = globenv % program_name // ' MOLECULAR SYSTEM'
  
  IF ( globenv % num_pe == 1 ) THEN
     natoms = SIZE ( ainp % c ( 1, : ) )
     ALLOCATE ( struc % part ( 1:natoms ), STAT = ios )
     IF ( ios /= 0 ) CALL stop_memory ( 'control', 'part', natoms )
     ALLOCATE ( struc % pnode ( 1:natoms ), STAT = ios )
     IF ( ios /= 0 ) CALL stop_memory ( 'control', 'pnode', natoms )
     nmol = SUM ( mol_setup ( : ) % num_mol )

     ALLOCATE ( struc % molecule ( 1:nmol ), STAT = ios )
     IF ( ios /= 0 ) CALL stop_memory ( 'control', 'molecule', nmol )

     IF ( globenv % ionode .AND. globenv % print_level > 3 ) THEN
        iw = globenv % scr

        WRITE ( iw, '( A )' )
        WRITE ( iw, '( A, T71, I10 )' ) &
             ' CONTROL| Number of allocated particles ', natoms
        WRITE ( iw,'( A, T71, I10 )' ) &
             ' CONTROL| Number of allocated particle nodes ', natoms
        WRITE ( iw, '( A, T71, I10 )' ) &
             ' CONTROL| Number of allocated molecules ', nmol
        WRITE ( iw, '( A )' )
     END IF
  ELSE

!..replicated data
     natoms = SIZE ( ainp % c ( 1, : ) )
     ALLOCATE ( struc % part ( 1:natoms ), STAT = ios )
     IF ( ios /= 0 ) CALL stop_memory ( 'control', 'part', natoms )
     nmoltype = SIZE ( mol_setup )
     nmol = 0
     nnodes = 0
     DO i = 1, nmoltype
        nsh = get_share ( mol_setup ( i ) % num_mol, &
             globenv % num_pe, globenv % mepos )
        nmol = nmol + nsh
        nnodes = nnodes + nsh * mol_setup ( i ) % molpar % natom
     END DO

     ALLOCATE ( struc % molecule ( 1:nmol ), STAT = ios )
     IF ( ios /= 0 ) CALL stop_memory ( 'control', 'molecule' , nmol )
     ALLOCATE ( struc % pnode ( 1:nnodes ), STAT = ios )
     IF ( ios /= 0 ) CALL stop_memory ( 'control', 'pnode', nnodes )

     IF ( globenv % ionode .AND. globenv % print_level > 3 ) THEN
        iw = globenv % scr
        WRITE ( iw, '( A )' )
        WRITE ( iw, '( A, T71, I10 )' ) &
             ' CONTROL| Number of allocated particles ', natoms
        WRITE ( iw, '( A, T71, I10 )' ) &
             ' CONTROL| Number of allocated particle nodes ', nnodes
        WRITE ( iw, '( A, I5, T71, I10 )' ) &
             ' CONTROL| Number of allocated molecules on processor ', &
             globenv % mepos, nmol
        WRITE ( iw, '( A )' )
     END IF

  END IF

END SUBROUTINE allocmem

!******************************************************************************

SUBROUTINE deallocmem ( struc )
  IMPLICIT NONE

! Arguments
  TYPE ( structure_type ), INTENT ( INOUT ) :: struc

! Locals
  INTEGER :: ios

!------------------------------------------------------------------------------

  DEALLOCATE ( struc % part, STAT = ios )
  IF ( ios /= 0 ) CALL stop_memory ( 'control', 'part' )

  DEALLOCATE ( struc % pnode, STAT = ios )
  IF ( ios /= 0 ) CALL stop_memory ( 'control', 'pnode' )

  DEALLOCATE ( struc % molecule, STAT = ios )
  IF ( ios /= 0 ) CALL stop_memory ( 'control', 'molecule' )

END SUBROUTINE deallocmem

!******************************************************************************

END MODULE control_module
