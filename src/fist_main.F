!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/fist_main [1.0] *
!!
!!   NAME
!!     fist_main
!!
!!   FUNCTION
!!     perform classical molecular dynamics and path integral simulations
!!   AUTHOR
!!     CJM-Sept-01-02
!!
!!   MODIFICATION HISTORY
!!     gt SEPT-23-2002: part is allocated/deallocated/initialized in
!!                      read_coord_vel
!!   SOURCE
!******************************************************************************

MODULE fist_main
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE atoms_input,                     ONLY: read_coord_vel,&
                                             set_particle_properties,&
                                             system_type
  USE cntl_input,                      ONLY: read_cntl_section
  USE convert_units,                   ONLY: convert
  USE dump,                            ONLY: dump_variables
  USE ewald_parameters_types,          ONLY: ewald_parameters_type
  USE ewalds,                          ONLY: ewald_print
  USE extended_system_types,           ONLY: extended_system_type
  USE fist_debug,                      ONLY: fist_debug_control => debug_control
  USE fist_environment_types,          ONLY: fist_environment_type,&
                                             fist_thermodynamic_type,&
                                             init_fist_env,&
                                             replica_environment_type,&
                                             set_fist_env
  USE force_control,                   ONLY: force,&
                                             set_force_env
  USE force_fields,                    ONLY: read_force_field_section
  USE global_types,                    ONLY: global_environment_type
  USE header,                          ONLY: fist_header
  USE initialize_extended_types,       ONLY: assign_extended_parameters,&
                                             initialize_nhc_baro,&
                                             initialize_nhc_part,&
                                             initialize_npt_type
  USE initialize_molecule_types,       ONLY: initialize_molecule_type
  USE initialize_pimd_types,           ONLY: initialize_pimd
  USE input_types,                     ONLY: setup_parameters_type
  USE integrator,                      ONLY: set_integrator
  USE kinds,                           ONLY: dbl
  USE linklist_types,                  ONLY: initialize_linklist_data,&
                                             linklist_internal_data_type,&
                                             set_linklist_internal_data
  USE mathconstants,                   ONLY: twopi,&
                                             zero
  USE md,                              ONLY: initialize_velocities,&
                                             mdio_parameters_type,&
                                             read_md_section,&
                                             simulation_parameters_type,&
                                             virial_type
  USE md_environment_types,            ONLY: init_md_env,&
                                             md_environment_type,&
                                             set_md_env,&
                                             zero_virial
  USE md_fist_energies,                ONLY: md_energy_fist,&
                                             md_energy_fist_pimd,&
                                             set_fist_energies,&
                                             set_fist_energies_param
  USE molecule_input,                  ONLY: read_molecule_section,&
                                             read_setup_section,&
                                             read_topology_section
  USE molecule_types,                  ONLY: intra_parameters_type,&
                                             molecule_type,&
                                             topology_type
  USE pair_potential,                  ONLY: potentialparm_type,&
                                             spline_nonbond_control
  USE particle_types,                  ONLY: particle_prop_type,&
                                             particle_type
  USE physcon,                         ONLY: boltzmann,&
                                             h_planck
  USE read_pimd,                       ONLY: read_pimd_section
  USE simulation_cell,                 ONLY: get_hinv,&
                                             init_cell
  USE string_utilities,                ONLY: integer_to_string,&
                                             xstring
  USE structure_types,                 ONLY: init_structure_type,&
                                             set_structure_type,&
                                             structure_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop,&
                                             trace_debug
  USE transformations,                 ONLY: fr2fu,&
                                             ud2v,&
                                             v2ud
  USE unit,                            ONLY: set_units,&
                                             unit_convert_type
  USE util,                            ONLY: close_unit,&
                                             get_share
  USE velocity_verlet_control,         ONLY: velocity_verlet
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: fist
  
!!*****
!-----------------------------------------------------------------------------!

CONTAINS

!-----------------------------------------------------------------------------!
! FIST FIST FIST FIST FIST FIST FIST FIST FIST FIST FIST FIST FIST FIST FIST  !
!-----------------------------------------------------------------------------!
!!****** fist_main/fist [1.0] *
!!
!!   NAME
!!     fist
!!
!!   FUNCTION
!!     Controls program flow for classical MD and path-integrals
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!
!!   USED BY
!!     cp2k
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE fist ( globenv )

  IMPLICIT NONE

! Argument
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv

! Locals
  REAL ( dbl ), DIMENSION ( :, : ), ALLOCATABLE :: rcut
  REAL ( dbl ), DIMENSION ( :, : ), ALLOCATABLE :: rcut_cgf
  INTEGER :: handle1, handle2, isos
  INTEGER :: ibead, beads, ia, ib
  INTEGER :: n, nll, i
  INTEGER :: maxder

  CHARACTER ( LEN = 40 ) :: set_fn, project_name
  CHARACTER ( LEN = 5 ) :: tag
  LOGICAL, POINTER :: pimd
  
  TYPE ( extended_system_type ) :: extended_type
  TYPE ( ewald_parameters_type ), POINTER :: ewald_param
  TYPE ( fist_environment_type ) :: fist_env
  TYPE ( intra_parameters_type ) :: intra_param
  TYPE ( linklist_internal_data_type ), POINTER :: ll_data
  TYPE ( md_environment_type ) :: md_env
  TYPE ( mdio_parameters_type ) :: mdio
  TYPE ( molecule_type ), DIMENSION ( : ), POINTER :: mol_setup
  TYPE ( particle_prop_type ), DIMENSION ( : ), POINTER :: pstat
  TYPE ( potentialparm_type ), DIMENSION ( :, : ), POINTER :: potparm
  TYPE ( replica_environment_type ), DIMENSION ( : ), POINTER :: rep_env
  TYPE ( simulation_parameters_type ), POINTER :: simpar
  TYPE ( setup_parameters_type ) :: setup
  TYPE ( system_type ) :: ainp
  TYPE ( fist_thermodynamic_type ), POINTER :: thermo
  TYPE ( topology_type ) :: topo
  TYPE ( unit_convert_type ) :: units
  TYPE ( virial_type ), POINTER :: virial
  
!------------------------------------------------------------------------------

! IF( globenv % ionode ) CALL trace_debug ( "START" )

  CALL timeset ( 'FIST_MAIN', 'I', ' ', handle1 )
  CALL timeset ( 'FIST_INIT', 'I', ' ', handle2 )
  
! initialize fist_env
  CALL init_fist_env ( fist_env )

! allocating or nullifying local pointers
  ALLOCATE ( simpar, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'allocstruc', 'simpar', 0 )
  ALLOCATE ( pimd, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'allocstruc', 'pimd', 0 )
  ALLOCATE ( ewald_param, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'allocstruc', 'ewald_param', 0 )
  ALLOCATE ( thermo, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'allocstruc', 'thermo', 0 )
  ALLOCATE ( virial, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'allocstruc', 'virial', 0 )
  NULLIFY ( potparm )
  NULLIFY ( rep_env )
! initialize the virial
  CALL zero_virial ( virial )
! point from fist_thermo % virial to virial
  thermo % virial => virial

! read control section
  CALL read_cntl_section ( setup, globenv, ewald_param )

  IF ( setup % path_integrals ) THEN
    pimd = .TRUE.
    CALL read_pimd_section ( simpar % pimd_params, globenv )
    beads = simpar % pimd_params % beads
  ELSE
    pimd = .FALSE.
    beads = 1
    simpar % pimd_params % beads = 0
  END IF


! Allocate particle thermostat for general use with path-integrals
  ALLOCATE ( extended_type % nhc_part ( beads ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'allocstruc', 'nhc_part', beads )

! Allocate barostat thermostat for general use with path-integrals
  ALLOCATE ( extended_type % nhc_baro ( beads ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'allocstruc', 'nhc_baro', beads )

! Allocate replica_environment_type for general use with path integrals
  ALLOCATE ( rep_env ( beads ), STAT=isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'allocstruc', 'rep_env', beads )

! Allocating and setting up the INTERNAL_DATA_TYPE and CELL_TYPE for  
! general use with path integrals
  DO ibead = 1, beads
    ALLOCATE ( rep_env ( ibead ) % box, STAT = isos )
    IF ( isos /=0 ) CALL stop_memory ( 'allocstruc', 'box', 1 )
    ALLOCATE ( rep_env ( ibead ) % box_ref, STAT = isos )
    IF ( isos /=0 ) CALL stop_memory ( 'allocstruc', 'box_ref', 1 )
    ALLOCATE ( rep_env ( ibead ) % ll_data ( 1 ), STAT = isos )
    IF ( isos /=0 ) CALL stop_memory ( 'allocstruc', 'll_data', 1 )
    ll_data => rep_env ( ibead ) % ll_data ( 1 ) 
    CALL initialize_linklist_data ( ll_data )
  END DO

  IF ( globenv % ionode ) THEN
    CALL fist_header ( globenv % scr )
  END IF
  

! read from the &SETUP and &MOLECULE section of *.set

  set_fn = setup % set_file_name

  CALL read_setup_section ( mol_setup, setup, globenv )

  CALL read_molecule_section ( mol_setup, setup, globenv )

  CALL read_topology_section ( topo, setup, globenv )

! read force_field information for classical MD in *.set

  CALL read_force_field_section ( setup, set_fn, intra_param, &
                                  potparm, globenv )

! read the input of the molecular dynamics section

  CALL read_md_section ( simpar, globenv, mdio )
  simpar % program = globenv % program_name

! setup the particle properties
  
  CALL set_particle_properties ( setup, mol_setup, pstat )

!..read atomic coordinates, velocities (optional) and the simulation box
  ainp % rtype = simpar % read_type

! initialize working units
  CALL set_units ( setup % unit_type, units )

!..1.0e-15 because convert_unit expects time in [fs]^-1 not [s]^-1
  simpar % pimd_params % wp = 1.0e-15_dbl * SQRT ( 1.0_dbl * beads ) * twopi * &
                              simpar % temp_ext * boltzmann / h_planck
  simpar % pimd_params % beta = 1.0_dbl / simpar % temp_ext

  CALL convert ( units = units, simpar = simpar, &
                 pstat = pstat, potparm = potparm,  &
                 intra_param = intra_param,  &
                 ewald_param = ewald_param )

  CALL xstring ( setup % coord_file_name, ia, ib )

! Filling up the replica environment one bead at a time
  DO ibead = 1, beads

    IF ( pimd ) THEN
      CALL integer_to_string( ibead, tag )
      project_name = setup % coord_file_name ( ia : ib ) // '_' // ADJUSTL ( tag )
    ELSE
      project_name = setup % coord_file_name ( ia : ib )
    END IF

    CALL read_coord_vel ( rep_env (ibead) % part, rep_env ( ibead ) % atomic_kind_set,&
                          ainp, project_name, mol_setup, setup % atom_names, &
                          pstat, globenv )

!..initialize box, perd
    IF (ainp%n > 0) THEN
      CALL init_cell ( rep_env ( ibead ) % box, ainp % box, &
                       setup % simulation_cell % perd )
    ELSE 
      CALL init_cell ( rep_env ( ibead ) % box, setup % simulation_cell % hmat, &
                       setup % simulation_cell % perd )
    END IF

! dbg TEMPORARY FIX
    rep_env ( ibead ) % box_ref % hmat = setup % reference_cell % hmat
    rep_env ( ibead ) % box % unit_of_length_name = "ANGSTROM"
    rep_env ( ibead ) % box % unit_of_length = 1.0_dbl
    rep_env ( ibead ) % box % scaled_coordinates = .FALSE.
    rep_env ( ibead ) % box % subcells = simpar % subcells
! dbg TEMPORARY FIX

! If run is a debug.  Make box_ref = box to work under
! all ensembles and restart options
    IF ( setup % run_type == 'DEBUG' ) THEN
      IF (ainp%n > 0) THEN
         rep_env ( ibead ) % box_ref % hmat = ainp % box
      ELSE 
         rep_env ( ibead ) % box_ref % hmat = setup % simulation_cell % hmat
      END IF
    ENDIF

!..allocate memory for atoms and molecules
    CALL allocmem ( ainp, mol_setup, rep_env ( ibead ), globenv )

!..convert the units ( includes part and box ONLY when read_type is 'INIT' )
    IF ( simpar % read_type == 'INIT' ) &
    CALL convert ( units = units, part = rep_env ( ibead ) % part, &
                   box = rep_env ( ibead ) % box )
!..calculate the inverse box matrix now after the unit conversion,
!  so it also has the right units and assign the reference box
    CALL get_hinv ( rep_env ( ibead ) % box )
    IF ( simpar % read_type == 'INIT' ) THEN
       rep_env ( ibead ) % box_ref = rep_env ( ibead ) % box
    ENDIF
    CALL get_hinv ( rep_env ( ibead ) % box_ref )

!..initialize molecule_type
    CALL initialize_molecule_type ( mol_setup, intra_param, &
                                    rep_env ( ibead ) % pnode, &
                                    rep_env ( ibead ) % part, &
                                    rep_env ( ibead ) % molecule, &
                                    globenv )

!..allocate lnhc_parameters_type for particles and get
!  number of degrees of freedom and initialize if necessary
    CALL initialize_nhc_part ( rep_env ( ibead ) % box, simpar, &
      rep_env ( ibead ) % molecule, mol_setup, globenv, &
      extended_type % nhc_part ( ibead ) )

!..allocate lnhc_parameters_type for barostat and initialize if
!  necessary
    CALL initialize_nhc_baro ( simpar, globenv, extended_type % nhc_baro ( ibead ) )

!..allocate npt_info_type and initialize if necessary
    CALL initialize_npt_type ( simpar, globenv, extended_type % npt_info, &
                               rep_env ( ibead ) % box )


! initialize velocities if needed
    IF ( simpar % read_type == 'POS' .OR. simpar % read_type == 'INIT' ) &
    CALL initialize_velocities ( simpar, rep_env ( ibead ) % part, globenv )

! Assign extended system variables
    IF ( simpar % read_type == 'ALL' ) &
    CALL assign_extended_parameters ( rep_env ( ibead ) % box % deth,    &
                                      simpar % ensemble, ainp,           &
                                      extended_type % npt_info,           &
                                      extended_type % nhc_part ( ibead ), &
                                      extended_type % nhc_baro ( ibead ) )

    CALL release_mem ( ainp )

  END DO

  CALL set_fist_energies_param ( units % pconv, units % econv, units % l_label, &
       units % vol_label, units % e_label, units % pres_label, &
       units % temp_label, units % angl_label )


!...initialize splines

  potparm ( :, : ) % energy_cutoff = 0.0_dbl
  potparm ( :, : ) % e_cutoff_coul = 0.0_dbl
  CALL spline_nonbond_control ( potparm, pstat, 2500, ewald_param )

!..set linklist control parameters
  ALLOCATE ( rcut ( setup % natom_type, setup % natom_type ), STAT = isos )
  IF ( isos /=0 ) CALL stop_memory ( 'fist_main', 'rcut', 0 )

  rcut ( :, : ) = potparm ( :, : ) % rcutsq

  DO ibead = 1, beads 
    CALL set_linklist_internal_data ( rep_env ( ibead ) % ll_data ( 1 ), globenv,      &
                                      simpar % verlet_skin, setup % natom_type, rcut,  &
                                      simpar % subcells )

    CALL set_linklist_internal_data ( rep_env ( ibead ) % ll_data ( 1 ), globenv,  &
                                      printlevel = globenv % print_level )
               
  END DO

  DEALLOCATE ( rcut, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'fist_main', 'rcut' )

!
  CALL timestop ( zero, handle2 )

! set the fist_env
  CALL set_fist_env ( fist_env, pimd, rep_env, ewald_param, potparm, thermo )

! begin the MD process
  CALL control_work ( globenv, units, simpar, fist_env, md_env, extended_type, &
                      setup, mdio )

! deallocate memory for atoms and molecules
  DO ibead = 1, beads
    CALL deallocmem ( rep_env ( ibead ) )
  END DO

END SUBROUTINE fist


!-----------------------------------------------------------------------------!
! FIST FIST FIST FIST FIST FIST FIST FIST FIST FIST FIST FIST FIST FIST FIST  !
!-----------------------------------------------------------------------------!

SUBROUTINE release_mem ( ainp )
  IMPLICIT NONE
! Arguments
  TYPE ( system_type ) :: ainp
! Local
  INTEGER :: isos

!..deallocate arrays needed for atom input

  IF ( ASSOCIATED ( ainp % c ) ) THEN
     DEALLOCATE ( ainp % c, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'fist_main', 'ainp%c' )
  END IF

  IF ( ASSOCIATED ( ainp % v ) ) THEN
     DEALLOCATE ( ainp % v, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'fist_main', 'ainp%v' )
  END IF

  IF ( ASSOCIATED ( ainp % eta_part ) ) THEN
     DEALLOCATE ( ainp % eta_part, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'fist_main', 'ainp%eta_part' )
  END IF

  IF ( ASSOCIATED ( ainp % veta_part ) ) THEN
     DEALLOCATE ( ainp % veta_part, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'fist_main', 'ainp%veta_part' )
  END IF

  IF ( ASSOCIATED ( ainp % mass_eta_part ) ) THEN
     DEALLOCATE ( ainp % mass_eta_part, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'fist_main', 'ainp%mass_eta_part' )
  END IF

  IF ( ASSOCIATED ( ainp % eta_baro ) ) THEN
     DEALLOCATE ( ainp % eta_baro, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'fist_main', 'ainp%eta_baro' )
  END IF

  IF ( ASSOCIATED ( ainp % veta_baro ) ) THEN
     DEALLOCATE ( ainp % veta_baro, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'fist_main', 'ainp%veta_baro' )
  END IF

  IF ( ASSOCIATED ( ainp % mass_eta_baro ) ) THEN
     DEALLOCATE ( ainp % mass_eta_baro, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'fist_main', 'ainp%mass_eta_baro' )
  END IF

  IF ( ASSOCIATED ( ainp % veps ) ) THEN
     DEALLOCATE ( ainp % veps, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'fist_main', 'ainp%veps' )
  END IF

  IF ( ASSOCIATED ( ainp % mass_eps ) ) THEN
     DEALLOCATE ( ainp % mass_eps, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'fist_main', 'ainp%mass_eps' )
  END IF 
END SUBROUTINE release_mem
!******************************************************************************

SUBROUTINE control_work ( globenv, units, simpar, fist_env, md_env, &
                          extended_type, setup, mdio )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
  TYPE ( unit_convert_type ) :: units
  TYPE ( simulation_parameters_type ), POINTER :: simpar
  TYPE ( fist_environment_type ), INTENT ( INOUT ) :: fist_env 
  TYPE ( md_environment_type ), INTENT ( INOUT ) :: md_env
  TYPE ( extended_system_type ) :: extended_type
  TYPE ( setup_parameters_type ) :: setup
  TYPE ( mdio_parameters_type ) :: mdio
  
! Locals
  INTEGER :: handle, ibead, beads, isos
  TYPE ( structure_type ), POINTER :: struc ( : )
  
!------------------------------------------------------------------------------
  
  CALL timeset ( 'CNTL_WORK', 'I', ' ', handle )

! allocate structure_type
  beads = SIZE ( fist_env % rep_env )
  ALLOCATE ( struc ( beads ), STAT = isos )  
  IF ( isos /= 0 ) CALL stop_memory ( 'allocstruc', 'struc', beads )

! set up structure_type
  DO ibead = 1, beads
    CALL init_structure_type ( struc ( ibead ) )
    CALL set_structure_type ( struc ( ibead ), rep_env = fist_env % rep_env ( ibead ) )
  ENDDO

! set the force environment
  CALL set_force_env ( fist_env_type = fist_env )

! set the integrator environment
  CALL init_md_env ( md_env )
  CALL set_md_env ( md_env, simpar=simpar, struc=struc,           &
                            virial=fist_env % thermo % virial,    &
                            nhc_part=extended_type % nhc_part,    &
                            nhc_baro=extended_type % nhc_baro,    &
                            npt=extended_type % npt_info,         &
                            cell=fist_env % rep_env ( 1 ) % box )
  ALLOCATE ( md_env % itimes, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'fist_main', 'itimes', 0 )
  ALLOCATE ( md_env % constant, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'fist_main', 'constant', 0 )

  SELECT CASE ( setup % run_type )
  CASE DEFAULT
     CALL stop_program ( "control_work", "no suitable run_type" &
                         //  setup % run_type )
     
  CASE ( "DEBUG" )
     
! debug the forces
     CALL control_debug_work ( globenv, fist_env, mdio )
     
  CASE ( "MD" )
     
     CALL control_md_work ( globenv, fist_env, md_env, units, mdio )
     
  END SELECT
  
  DEALLOCATE ( struc, STAT = isos )  
  IF ( isos /= 0 ) CALL stop_memory ( 'deallocstruc', 'struc' )
  CALL timestop ( zero, handle )
  
END SUBROUTINE control_work

!******************************************************************************

SUBROUTINE control_debug_work ( globenv, fist_env, mdio )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
  TYPE ( fist_environment_type ), INTENT ( INOUT ) :: fist_env
  TYPE ( mdio_parameters_type ) :: mdio
  
!------------------------------------------------------------------------------
  
! initialize integrator
  CALL set_integrator ( globenv )
! debug routine  
  CALL fist_debug_control ( globenv,  fist_env % rep_env ( 1 ), &
                            fist_env % ewald_param,             &
                            fist_env % potparm, fist_env % thermo )
  
END SUBROUTINE control_debug_work

!******************************************************************************

SUBROUTINE control_md_work ( globenv, fist_env, md_env, units, mdio )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
  TYPE ( fist_environment_type ), INTENT ( INOUT ) :: fist_env
  TYPE ( md_environment_type ), INTENT ( INOUT ) :: md_env
  TYPE ( unit_convert_type ) :: units
  TYPE ( mdio_parameters_type ) :: mdio
  
! Locals
  INTEGER, POINTER :: itimes
  LOGICAL, PARAMETER :: box_change = .FALSE.
  
!------------------------------------------------------------------------------

! initialize the environment for md_energies
  CALL set_fist_energies ( globenv, mdio )

! initialize integrator
  CALL set_integrator ( globenv )

  IF ( fist_env % pimd ) THEN
!..initialize pimd_data (transformation related variables)
    CALL initialize_pimd ( md_env % struc , md_env % simpar )

    CALL v2ud ( md_env % struc, md_env % simpar )
  END IF
  
! MD
  itimes => md_env % itimes 
  itimes = 0
  
  CALL force ( globenv, box_change )

  IF ( globenv % ionode .AND. fist_env % ewald_param % ewald_type /= 'NONE' ) &
       CALL ewald_print ( globenv % scr, md_env % cell,                 &
                          fist_env % thermo % e_gspace,                &
                          fist_env % thermo % e_self,                  &
                          fist_env % thermo % e_neut,                  &
                          fist_env % thermo % e_bonded,                &
                          units % e_label )
  
  IF ( fist_env % pimd ) CALL fr2fu ( md_env % struc, md_env % simpar )


! Initial call to md energy
  IF ( fist_env % pimd ) THEN
    CALL md_energy_fist_pimd ( md_env % itimes, md_env % constant,       &
                               fist_env % rep_env ( 1 ) % box,           & 
                               fist_env % rep_env ( 1 ) % ll_data ( 1 ), &
                               md_env % nhc_part, md_env % nhc_baro,     &
                               md_env % npt, md_env % simpar,            &
                               md_env % struc, fist_env % thermo )
  ELSE 
    CALL md_energy_fist ( md_env % itimes, md_env % constant,       &
                          fist_env % rep_env ( 1 ) % box,           & 
                          fist_env % rep_env ( 1 ) % ll_data ( 1 ), &
                          md_env % nhc_part ( 1 ), md_env % nhc_baro ( 1 ),  &
                          md_env % npt, md_env % simpar,            &
                          md_env % struc ( 1 ), fist_env % thermo )
  ENDIF

! Main MD loop
  DO itimes = 1, md_env % simpar % nsteps

! Call the integrator
     CALL velocity_verlet ( md_env )

! Call the energy routines
     IF ( fist_env % pimd ) THEN
       CALL md_energy_fist_pimd ( md_env % itimes, md_env % constant,       &
                                  fist_env % rep_env ( 1 ) % box,           & 
                                  fist_env % rep_env ( 1 ) % ll_data ( 1 ), &
                                  md_env % nhc_part, md_env % nhc_baro,     &
                                  md_env % npt, md_env % simpar,            &
                                  md_env % struc, fist_env % thermo )
     ELSE 
       CALL md_energy_fist ( md_env % itimes, md_env % constant,       &
                             fist_env % rep_env ( 1 ) % box,           & 
                             fist_env % rep_env ( 1 ) % ll_data ( 1 ), &
                             md_env % nhc_part ( 1 ), md_env % nhc_baro ( 1 ),  &
                             md_env % npt, md_env % simpar,            &
                             md_env % struc ( 1 ), fist_env % thermo )
     ENDIF
     
! Test for the dump cycle
     IF ( MOD ( itimes, mdio % idump ) == 0 ) THEN
          IF ( fist_env % pimd ) CALL ud2v ( md_env % struc, md_env % simpar )
          CALL dump_variables ( md_env, fist_env % rep_env, mdio % dump_file_name, &
                                globenv, fist_env % pimd )
          IF ( fist_env % pimd ) CALL v2ud ( md_env % struc, md_env % simpar )
     END IF

  END DO
  
  IF ( fist_env % pimd ) CALL ud2v ( md_env % struc, md_env % simpar )

  CALL dump_variables ( md_env, fist_env % rep_env, mdio % dump_file_name, globenv,  &
                        fist_env % pimd )
  
  IF ( fist_env % pimd ) CALL v2ud ( md_env % struc, md_env % simpar )

  IF ( globenv % ionode ) CALL close_unit ( 10, 99 )
  
END SUBROUTINE control_md_work

!******************************************************************************

SUBROUTINE allocmem ( ainp, mol_setup, rep_env, globenv )

  IMPLICIT NONE

! Arguments
  TYPE ( molecule_type ), DIMENSION ( : ), INTENT ( IN ) :: mol_setup
  TYPE ( replica_environment_type ), INTENT ( INOUT ) :: rep_env
  TYPE ( system_type ), INTENT ( IN ) :: ainp
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv

! Locals
  INTEGER :: iw, natoms, nnodes, nmol, nmoltype, ios, iat, i, nsh

!------------------------------------------------------------------------------

  IF ( globenv % num_pe == 1 ) THEN
     IF ( ainp % n > 0 ) THEN
        natoms = SIZE ( ainp % c, 2 )
     ELSE
        natoms = SIZE ( rep_env % part )
     END IF

     ALLOCATE ( rep_env % pnode ( natoms ), STAT = ios )
     IF ( ios /= 0 ) CALL stop_memory ( 'fist_main', 'pnode', natoms )

     nmol = SUM ( mol_setup ( : ) % num_mol )
     ALLOCATE ( rep_env % molecule ( nmol ), STAT = ios )
     IF ( ios /= 0 ) CALL stop_memory ( 'fist_main', 'molecule', nmol )

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
     IF ( ainp % n > 0 ) THEN
        natoms = SIZE ( ainp % c, 2  )
     ELSE
        natoms = SIZE ( rep_env % part )
     END IF
     nmoltype = SIZE ( mol_setup )
     nmol = 0
     nnodes = 0
     DO i = 1, nmoltype
        nsh = get_share ( mol_setup ( i ) % num_mol, &
             globenv % num_pe, globenv % mepos )
        nmol = nmol + nsh
        nnodes = nnodes + nsh * mol_setup ( i ) % molpar % natom
     END DO

     ALLOCATE ( rep_env % molecule ( nmol ), STAT = ios )
     IF ( ios /= 0 ) CALL stop_memory ( 'fist_main', 'molecule' , nmol )
     ALLOCATE ( rep_env % pnode ( nnodes ), STAT = ios )
     IF ( ios /= 0 ) CALL stop_memory ( 'fist_main', 'pnode', nnodes )

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

SUBROUTINE deallocmem ( rep_env )
  IMPLICIT NONE

! Arguments
  TYPE ( replica_environment_type ), INTENT ( INOUT ) :: rep_env

! Locals
  INTEGER :: ios

!------------------------------------------------------------------------------

  DEALLOCATE ( rep_env % pnode, STAT = ios )
  IF ( ios /= 0 ) CALL stop_memory ( 'fist_main', 'pnode' )

  DEALLOCATE ( rep_env % molecule, STAT = ios )
  IF ( ios /= 0 ) CALL stop_memory ( 'fist_main', 'molecule' )

END SUBROUTINE deallocmem

!******************************************************************************

END MODULE fist_main
