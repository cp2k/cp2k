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
  USE atoms_input,                     ONLY: read_atoms_input,&
                                             read_md_restart,&
                                             read_pimd_atoms_input
  USE cntl_input,                      ONLY: read_cntl_section
  USE distribution_1d_types,           ONLY: distribution_1d_build_gdist,&
                                             distribution_1d_destroy_gdist,&
                                             distribution_1d_p_type,&
                                             distribution_1d_type
  USE dump,                            ONLY: dump_variables
  USE ewald_parameters_types,          ONLY: ewald_parameters_type, &
                                             read_ewald_section
  USE ewalds,                          ONLY: ewald_print
  USE extended_system_dynamics,        ONLY: lnhc_parameters_type,&
                                             npt_info_type
  USE extended_system_types,           ONLY: extended_system_type
  USE fist_debug,                      ONLY: fist_debug_control => debug_control
  USE fist_environment_types,          ONLY: fist_environment_type,&
                                             fist_thermodynamic_type,&
                                             init_fist_env,&
                                             replica_env_init,&
                                             replica_environment_type,&
                                             set_fist_env
  USE force_env_types,                   ONLY: force_env_calculate_force,&
                                             force_env_create,&
                                             force_env_release,&
                                             force_env_type
  USE force_fields,                    ONLY: read_force_field_section
  USE global_types,                    ONLY: global_environment_type
  USE header,                          ONLY: fist_header
  USE initialize_extended_types,       ONLY: initialize_nhc_baro,&
                                             initialize_nhc_forces,&
                                             initialize_nhc_part,&
                                             initialize_npt_type
  USE initialize_molecule_types,       ONLY: initialize_molecule_type,&
                                             mol_local_particles_create
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
  USE md_environment_types,            ONLY: destroy_md_env,&
                                             init_md_env,&
                                             md_environment_type,&
                                             set_md_env,&
                                             zero_virial
  USE md_fist_energies,                ONLY: md_energy_fist,&
                                             set_fist_energies
  USE molecule_kind_types,             ONLY: allocate_molecule_kind_set
  USE molecule_types,                  ONLY: intra_parameters_type,&
                                             molecule_type
  USE pair_potential_types,            ONLY: pair_potential_type
  USE pair_potential,                  ONLY: spline_nonbond_control
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: boltzmann,&
                                             h_planck
  USE read_pimd,                       ONLY: read_pimd_section
  USE setup_input,                     ONLY: read_setup_section
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_hinv,&
                                             read_cell,&
                                             write_cell
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
  USE topology,                        ONLY: topology_control
  USE topology_input,                  ONLY: read_topology_section
  USE topology_types,                  ONLY: topology_parameters_type
  USE transformations,                 ONLY: fr2fu,&
                                             ud2v,&
                                             v2ud
  USE util,                            ONLY: close_unit,&
                                             get_share
  USE velocity_verlet_control,         ONLY: velocity_verlet
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: fist
  
!!***
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
!!   SYNOPSIS
!!     Subroutine fist(globenv)
!!       Implicit None
!!       Type(global_environment_type), Intent (INOUT):: globenv
!!     End Subroutine fist
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
!!*** *************************************************************************

SUBROUTINE fist ( globenv )

  IMPLICIT NONE

! Argument
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv

! Locals
  REAL ( dbl ), DIMENSION ( :, : ), ALLOCATABLE :: rcut
  INTEGER :: handle2, isos
  INTEGER :: ibead, beads, ia, ib
  INTEGER :: i
  INTEGER :: maxder

  CHARACTER ( LEN = 40 ) :: set_fn, project_name
  CHARACTER ( LEN = 5 ) :: tag
  LOGICAL, POINTER :: pimd
  
  TYPE ( atomic_kind_type ), DIMENSION ( : ), POINTER :: atomic_kind_set
  TYPE ( particle_type ), POINTER :: particle_set ( : )
  TYPE ( cell_type ), POINTER :: cell, cell_ref

  TYPE ( extended_system_type ) :: extended_type
  TYPE ( ewald_parameters_type ), POINTER :: ewald_param
  TYPE ( fist_environment_type ), POINTER :: fist_env
  TYPE ( intra_parameters_type ) :: intra_param
  TYPE ( linklist_internal_data_type ), POINTER :: ll_data
  TYPE ( lnhc_parameters_type ), POINTER :: nhc_part, nhc_baro
  TYPE ( md_environment_type ) :: md_env
  TYPE ( mdio_parameters_type ) :: mdio
  TYPE ( molecule_type ), DIMENSION ( : ), POINTER :: mol_setup
  TYPE ( npt_info_type ), DIMENSION ( :, : ), POINTER :: npt_info
  TYPE ( pair_potential_type ), DIMENSION ( :, : ), POINTER :: potparm
  TYPE ( replica_environment_type ), DIMENSION ( : ), POINTER :: rep_env
  TYPE ( simulation_parameters_type ), POINTER :: simpar
  TYPE ( setup_parameters_type ) :: setup
  TYPE ( topology_parameters_type ) :: topology
  TYPE ( fist_thermodynamic_type ), POINTER :: thermo
  TYPE ( virial_type ), POINTER :: virial
  
!------------------------------------------------------------------------------

! IF( globenv % ionode ) CALL trace_debug ( "START" )

  CALL timeset ( 'FIST_INIT', 'I', ' ', handle2 )
  
! initialize fist_env
  ALLOCATE(fist_env,stat=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'allocstruc', 'fist_env', 0 )
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
  NULLIFY ( particle_set )
  NULLIFY ( atomic_kind_set )
  NULLIFY ( cell )
  NULLIFY ( cell_ref )
  NULLIFY ( ll_data )
  NULLIFY ( nhc_part )
  NULLIFY ( nhc_baro )
  NULLIFY ( npt_info )
  NULLIFY ( mol_setup )
! initialize the virial
  CALL zero_virial ( virial )
! point from fist_thermo % virial to virial
  thermo % virial => virial

! read control section
  CALL read_cntl_section ( setup, globenv )
! read ewald_section
  CALL read_ewald_section ( ewald_param, globenv )

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
    CALL replica_env_init( rep_env ( ibead ))
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
  
!*** Read the input section with the cell parameters ***

  CALL read_cell ( cell, globenv )

!*** Print the cell parameters ***

  CALL write_cell ( cell, globenv )

  topology%cell => cell

! read from the &SETUP and &MOLECULE section of *.set

  set_fn = setup % set_file_name

  CALL read_setup_section ( setup, globenv )

! topology_control
! CALL topology_control (atomic_kind_set,particle_set,mol_setup,&
!                        topology,setup,globenv)
  
! read force_field information for classical MD in *.set

  CALL read_force_field_section ( topology, set_fn, intra_param, &
                                  potparm, globenv )

! read the input of the molecular dynamics section

  CALL read_md_section ( simpar, globenv, mdio )
  simpar % program = globenv % program_name


! define project name string
  CALL xstring ( topology % coord_file_name, ia, ib )
  project_name = topology % coord_file_name ( ia : ib )

!!!!!! Needs to be integrated into cp2k_units.F
!..1.0e-15 because convert_unit expects time in [fs]^-1 not [s]^-1
!  simpar % pimd_params % wp = 1.0e-15_dbl * SQRT ( 1.0_dbl * beads ) * twopi * &
!                              simpar % temp_ext * boltzmann / h_planck
!  simpar % pimd_params % beta = 1.0_dbl / simpar % temp_ext
!!!!!! Needs to be integrated into cp2k_units.F

! Filling up the replica environment one bead at a time
  DO ibead = 1, beads

    IF ( pimd ) THEN
      CALL integer_to_string( ibead, tag )
      project_name = topology % coord_file_name ( ia : ib ) // '_' // ADJUSTL ( tag )
    END IF

!..initialize particle_set
    rep_env ( ibead ) % part => particle_set

!..initialize atomic_kind_set
    rep_env ( ibead ) % atomic_kind_set => atomic_kind_set

!..initialize cell parameters
    cell_ref => cell
    rep_env ( ibead ) % box => cell
    rep_env ( ibead ) % box_ref => cell_ref

! If run is a debug.  Make box_ref = box to work under
! all ensembles and restart options
    IF ( globenv % run_type == 'DEBUG' ) THEN
      rep_env ( ibead ) % box_ref => cell
    ENDIF

!..allocate memory for atoms and molecules

    CALL allocmem ( mol_setup, rep_env ( ibead ), globenv )

!..calculate the inverse box matrix now after the unit conversion,
!  so it also has the right units and assign the reference box
    CALL get_hinv ( rep_env ( ibead ) % box )
    IF ( simpar % read_type == 'INIT' ) THEN
       rep_env ( ibead ) % box_ref = rep_env ( ibead ) % box
    ENDIF
    CALL get_hinv ( rep_env ( ibead ) % box_ref )

!..initialize molecule_type
    CALL initialize_molecule_type ( mol_setup, rep_env ( ibead ) % pnode, &
                                    rep_env ( ibead ) % part, &
                                    rep_env ( ibead ) % molecule, &
                                    globenv, intra_param )

    CALL mol_local_particles_create ( &
      local_particles=rep_env(ibead)%local_particles,&
      atomic_kind_set=rep_env(ibead)%atomic_kind_set,&
      local_molecules=rep_env(ibead)%molecule,&
      particles=rep_env(ibead)%part,&
      para_env=globenv%para_env)
!MK   *** Distribute molecules and atoms using the new data structures ***
!MK
!MK    CALL allocate_molecule_kind_set(rep_env(ibead)%molecule_kind_set,&
!MK                                    SIZE(mol_setup))
!MK
!MK    CALL init_molecule_kind_set(rep_env(ibead)%molecule_kind_set,mol_setup)
!MK
!MK    CALL distribute_molecules(&
!MK      particle_kind_set=rep_env(ibead)%atomic_kind_set,&
!MK      particle_set=rep_env(ibead)%part,&
!MK      local_particles=rep_env(ibead)%local_particles,&
!MK      molecule_kind_set=rep_env(ibead)%molecule_kind_set,&
!MK      molecule_set=rep_env(ibead)%molecule_set,&
!MK      local_molecules=rep_env(ibead)%local_molecules,&
!MK      globenv=globenv)
!MK
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

! initialize velocities and read restart files if needed
    npt_info => extended_type % npt_info 
    nhc_part => extended_type % nhc_part ( ibead )
    nhc_baro => extended_type % nhc_baro ( ibead ) 
    SELECT CASE ( simpar % read_type )
    CASE ( 'INIT' )
      IF ( pimd )  CALL read_pimd_atoms_input ( simpar % read_type, &
                   globenv, rep_env ( ibead ) % part, rep_env ( ibead ) % box, &
                   rep_env ( ibead ) % box_ref, nhc_part, nhc_baro, npt_info )
      CALL initialize_velocities ( simpar, rep_env ( ibead ) % part, globenv )
    CASE ( 'POS' )
      IF ( pimd ) THEN
        CALL read_pimd_atoms_input ( simpar % read_type, globenv, &
                                    rep_env ( ibead ) %part,      &
                                    rep_env ( ibead ) % box,      &
                                    rep_env ( ibead ) % box_ref,  &
                                    nhc_part, nhc_baro, npt_info )
      ELSE
        CALL read_md_restart ( simpar % read_type, globenv, &
                               rep_env ( ibead ) %part,      &
                               rep_env ( ibead ) % box,      &
                               rep_env ( ibead ) % box_ref,  &
                               nhc_part, nhc_baro, npt_info )
      ENDIF
      CALL initialize_velocities ( simpar, rep_env ( ibead ) % part, globenv )
    CASE DEFAULT
      IF ( pimd ) THEN
        CALL read_pimd_atoms_input ( simpar % read_type, globenv, &
                                    rep_env ( ibead ) %part,      &
                                    rep_env ( ibead ) % box,      &
                                    rep_env ( ibead ) % box_ref,  &
                                    nhc_part, nhc_baro, npt_info )
      ELSE
        CALL read_md_restart ( simpar % read_type, globenv, &
                               rep_env ( ibead ) %part,      &
                               rep_env ( ibead ) % box,      &
                               rep_env ( ibead ) % box_ref,  &
                               nhc_part, nhc_baro, npt_info )
      ENDIF
    END SELECT
    CALL initialize_nhc_forces ( nhc_part, nhc_baro )
  END DO

!...initialize splines

  potparm ( :, : ) % energy_cutoff = 0.0_dbl
  potparm ( :, : ) % e_cutoff_coul = 0.0_dbl
  CALL spline_nonbond_control ( potparm, rep_env ( 1 ) % atomic_kind_set, 2500, ewald_param )

!..set linklist control parameters
  ALLOCATE ( rcut ( topology % natom_type, topology % natom_type ), STAT = isos )
  IF ( isos /=0 ) CALL stop_memory ( 'fist_main', 'rcut', 0 )

  rcut ( :, : ) = potparm ( :, : ) % rcutsq

  DO ibead = 1, beads 
    CALL set_linklist_internal_data ( rep_env ( ibead ) % ll_data ( 1 ), globenv,      &
                                      simpar % verlet_skin, topology % natom_type, rcut,  &
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
  CALL control_work ( globenv, simpar, fist_env, md_env, extended_type, &
                      setup, mdio )

! deallocate memory for atoms and molecules
  DO ibead = 1, beads
    CALL deallocmem ( rep_env ( ibead ) )
  END DO

END SUBROUTINE fist


!-----------------------------------------------------------------------------!
! FIST FIST FIST FIST FIST FIST FIST FIST FIST FIST FIST FIST FIST FIST FIST  !
!-----------------------------------------------------------------------------!

SUBROUTINE control_work ( globenv, simpar, fist_env, md_env, &
                          extended_type, setup, mdio )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
  TYPE ( simulation_parameters_type ), POINTER :: simpar
  TYPE ( fist_environment_type ), POINTER :: fist_env 
  TYPE ( md_environment_type ), INTENT(inout) :: md_env
  TYPE ( extended_system_type ) :: extended_type
  TYPE ( setup_parameters_type ) :: setup
  TYPE ( mdio_parameters_type ) :: mdio
  TYPE ( force_env_type ), POINTER :: force_env
  
! Locals
  INTEGER :: handle, ibead, beads, isos
  TYPE ( structure_type ), POINTER :: struc ( : )
  
!------------------------------------------------------------------------------
  
  CALL timeset ( 'CNTL_WORK', 'I', ' ', handle )
  NULLIFY ( force_env )

! allocate structure_type
  beads = SIZE ( fist_env % rep_env )
  ALLOCATE ( struc ( beads ), STAT = isos )  
  IF ( isos /= 0 ) CALL stop_memory ( 'allocstruc', 'struc', beads )

! set up structure_type
  DO ibead = 1, beads
    CALL init_structure_type ( struc ( ibead ) )
    CALL set_structure_type ( struc ( ibead ), rep_env = fist_env % rep_env ( ibead ) )
  ENDDO

! create the force environment
  CALL force_env_create ( force_env, fist_env = fist_env, &
       globenv=globenv)

! set the integrator environment
  CALL init_md_env ( md_env, globenv%para_env )
  CALL set_md_env ( md_env, simpar=simpar, struc=struc,           &
                            virial=fist_env % thermo % virial,    &
                            nhc_part=extended_type % nhc_part,    &
                            nhc_baro=extended_type % nhc_baro,    &
                            npt=extended_type % npt_info,         &
                            cell=fist_env % rep_env ( 1 ) % box,  &
                            force_env=force_env)
  ALLOCATE ( md_env % itimes, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'fist_main', 'itimes', 0 )
  ALLOCATE ( md_env % constant, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'fist_main', 'constant', 0 )
  CALL force_env_release(force_env)

  SELECT CASE ( globenv % run_type )
  CASE DEFAULT
     CALL stop_program ( "control_work", "no suitable run_type" &
                         //  globenv % run_type )
     
  CASE ( "DEBUG" )
     
! debug the forces
     CALL control_debug_work ( globenv, fist_env, mdio )
     
  CASE ( "MOLECULAR DYNAMICS" )
     
     CALL control_md_work ( globenv, fist_env, md_env, mdio )
     
  END SELECT

  
  DEALLOCATE ( struc, STAT = isos )  
  IF ( isos /= 0 ) CALL stop_memory ( 'deallocstruc', 'struc' )

  CALL destroy_md_env ( md_env )

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

SUBROUTINE control_md_work ( globenv, fist_env, md_env, mdio )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
  TYPE ( fist_environment_type ), POINTER :: fist_env
  TYPE ( md_environment_type ), INTENT ( INOUT ) :: md_env
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
  
  
  CALL force_env_calculate_force ( md_env % force_env, box_change=box_change )

  IF ( globenv % ionode .AND. fist_env % ewald_param % ewald_type /= 'NONE' ) &
       CALL ewald_print ( globenv % scr, md_env % cell,                &
                          fist_env % thermo % e_gspace,                &
                          fist_env % thermo % e_self,                  &
                          fist_env % thermo % e_neut,                  &
                          fist_env % thermo % e_bonded )
  
  IF ( fist_env % pimd ) CALL fr2fu ( md_env % struc, md_env % simpar )


! Initial call to md energy
  CALL md_energy_fist ( md_env , fist_env % rep_env ( 1 ) % box,  &
                               fist_env % rep_env ( 1 ) % ll_data ( 1 ), &
                               fist_env % thermo, fist_env % pimd )

! Main MD loop
  DO itimes = 1, md_env % simpar % nsteps

! Call the integrator
     CALL velocity_verlet ( md_env )

! Call the energy routines
     CALL md_energy_fist ( md_env , fist_env % rep_env ( 1 ) % box,  &
                               fist_env % rep_env ( 1 ) % ll_data ( 1 ), &
                               fist_env % thermo, fist_env % pimd )
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

SUBROUTINE allocmem ( mol_setup, rep_env, globenv )

  IMPLICIT NONE

! Arguments
  TYPE ( molecule_type ), DIMENSION ( : ), INTENT ( IN ) :: mol_setup
  TYPE ( replica_environment_type ), INTENT ( INOUT ) :: rep_env
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv

! Locals
  INTEGER :: iw, natoms, nnodes, nmol, nmoltype, ios, iat, i, nsh

!------------------------------------------------------------------------------

  IF ( globenv % num_pe == 1 ) THEN
     natoms = SIZE ( rep_env % part )
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
