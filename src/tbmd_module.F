!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/tbmd_module [1.0] *
!!
!!   NAME
!!     tbmd_module
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH, 22-NOV-2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE tbmd_module

  USE atoms_input,                     ONLY: read_coord_vel,&
                                             system_type
  USE band,                            ONLY: band_structure_type,&
                                             init_band_structure
  USE brillouin,                       ONLY: brillouin_info,&
                                             kpoint_input,&
                                             kpoint_type
  USE cntl_input,                      ONLY: read_cntl_section
  USE convert_units,                   ONLY: convert
  USE cp_control_types,                ONLY: dft_control_type
  USE dump,                            ONLY: dump_variables
  USE energy_input,                    ONLY: energy_min_type,&
                                             read_energy_section
  USE energy_optimize,                 ONLY: optimize_ao
  USE energy_output,                   ONLY: print_convergence,&
                                             print_energy_info
  USE ewald_parameters_types,          ONLY: ewald_parameters_type
  USE ewalds,                          ONLY: ewald_print
  USE extended_system_types,           ONLY: extended_system_type
  USE fermi,                           ONLY: fermi_distribution_type,&
                                             fermi_info,&
                                             init_fermi_dist
  USE force_env_types,                 ONLY: force
  USE force_fields,                    ONLY: read_force_field_section
  USE global_types,                    ONLY: global_environment_type
  USE header,                          ONLY: tbmd_header
  USE initialize_extended_types,       ONLY: assign_extended_parameters,&
                                             initialize_nhc_baro,&
                                             initialize_nhc_part,&
                                             initialize_npt_type
  USE initialize_molecule_types,       ONLY: initialize_molecule_type
  USE initialize_particle_types,       ONLY: initialize_particle_type
  USE initialize_pimd_types,           ONLY: initialize_pimd
  USE input_types,                     ONLY: setup_parameters_type
  USE integrator,                      ONLY: energy,&
                                             set_energy_parm,&
                                             set_integrator,&
                                             velocity_verlet
  USE kinds,                           ONLY: dbl
  USE kpoint_initialization,           ONLY: initialize_kpoints
  USE linklist_control,                ONLY: set_ll_parm
  USE mathconstants,                   ONLY: twopi,&
                                             zero
  USE md,                              ONLY: initialize_velocities,&
                                             mdio_parameters_type,&
                                             read_md_section,&
                                             simulation_parameters_type,&
                                             thermodynamic_type
  USE molecule_input,                  ONLY: read_molecule_section,&
                                             read_setup_section
  USE molecule_types,                  ONLY: intra_parameters_type,&
                                             molecule_type
  USE pair_potential,                  ONLY: spline_nonbond_control
  USE particle_types,                  ONLY: particle_prop_type
  USE physcon,                         ONLY: boltzmann,&
                                             h_planck
  USE read_pimd,                       ONLY: read_pimd_section
  USE simulation_cell,                 ONLY: get_hinv,&
                                             init_cell
  USE string_utilities,                ONLY: integer_to_string,&
                                             xstring
  USE structure_types,                 ONLY: interaction_type,&
                                             structure_type
  USE tbmd_debug,                      ONLY: tbmd_debug_control => debug_control
  USE tbmd_initialize,                 ONLY: tb_get_numel,&
                                             tbmd_init
  USE tbmd_input,                      ONLY: read_tb_hamiltonian,&
                                             read_tb_hopping_elements
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

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: tbmd
  
!!*****
!-----------------------------------------------------------------------------!

CONTAINS

!-----------------------------------------------------------------------------!
! CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL  !
!-----------------------------------------------------------------------------!
!!****** tbmd_module/tbmd [1.0] *
!!
!!   NAME
!!     tbmd
!!
!!   FUNCTION
!!     Controls program flow for TBMD
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE tbmd ( globenv )


! Argument
    TYPE(global_environment_type), &
      INTENT(INOUT)                          :: globenv

    CHARACTER(LEN=40)                        :: project_name, set_fn
    CHARACTER(LEN=5)                         :: tag
    INTEGER                                  :: beads, handle1, handle2, i, &
                                                ia, ib, ibead, isos, nel, nll
    LOGICAL                                  :: pimd
    REAL(dbl), ALLOCATABLE, DIMENSION(:, :)  :: rcut
    TYPE(band_structure_type)                :: bs
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(energy_min_type)                    :: energy_min
    TYPE(ewald_parameters_type)              :: ewald_param
    TYPE(extended_system_type)               :: extended_type
    TYPE(fermi_distribution_type)            :: fd
    TYPE(interaction_type)                   :: inter
    TYPE(intra_parameters_type)              :: intra_param
    TYPE(kpoint_type)                        :: kp
    TYPE(mdio_parameters_type)               :: mdio
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: mol_setup
    TYPE(particle_prop_type), DIMENSION(:), &
      POINTER                                :: pstat
    TYPE(setup_parameters_type)              :: setup
    TYPE(simulation_parameters_type)         :: simpar
    TYPE(structure_type), ALLOCATABLE, &
      DIMENSION(:)                           :: struc
    TYPE(system_type)                        :: ainp
    TYPE(unit_convert_type)                  :: units

!------------------------------------------------------------------------------
! IF( globenv % ionode ) CALL trace_debug ( "START" )

  CALL timeset ( 'CONTROL', 'I', ' ', handle1 )
  CALL timeset ( 'CNTL_INIT', 'I', ' ', handle2 )

! read control section
  CALL read_cntl_section ( setup, globenv, ewald_param )

  IF ( setup % path_integrals ) THEN
    pimd = .true.
    CALL read_pimd_section ( simpar % pimd_params, globenv )
    beads = simpar % pimd_params % beads
    CALL stop_program ( ' control ','PIMD not implemented with program' )
  ELSE
    pimd = .false.
    beads = 1
    simpar % pimd_params % beads = 0
  END IF

  ALLOCATE ( struc ( beads ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'allocstruc', 'struc', beads )

  ALLOCATE ( extended_type % nhc_part ( beads ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'allocstruc', 'nhc_part', beads )

  ALLOCATE ( extended_type % nhc_baro ( beads ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'allocstruc', 'nhc_baro', beads )

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

  CALL tbmd_header ( globenv % scr )
  

! read from the setup and molecule section of *.set

  set_fn = setup % set_file_name

  CALL read_setup_section ( mol_setup, setup, globenv )

  CALL read_molecule_section ( mol_setup, setup, globenv )

! read force_field information for classical MD
! read pair potential information for TB

  CALL read_force_field_section ( setup, mol_setup, set_fn, &
       intra_param, inter%potparm, pstat, globenv )

!..read Hamiltonian section

  CALL read_tb_hamiltonian ( setup, setup % atom_names, inter%tbatom, globenv )
  CALL read_tb_hopping_elements ( setup, setup % atom_names, &
       inter%tbatom, inter%tbhop, globenv )
  
 IF ( setup % run_type == 'ENERGY') THEN 
    CALL read_energy_section ( energy_min, globenv )
 END IF
      
! read the input of the molecular dynamics section
  CALL read_md_section ( simpar, globenv, mdio )
  simpar % program = globenv % program_name

!..read atomic coordinates, velocities (optional) and the simulation box
  ainp % rtype = simpar % read_type

! initialize working units
  CALL set_units ( setup % unit_type, units )

!..1.0e-15 because convert_unit expects time in [fs]^-1 not [s]^-1
  simpar % pimd_params % wp = 1.0e-15_dbl * SQRT ( 1.0_dbl * beads ) * twopi * &
        simpar % temp_ext * boltzmann / h_planck
  simpar % pimd_params % beta = 1.0_dbl / simpar % temp_ext
  CALL convert ( units = units, simpar = simpar, &
                pstat = pstat, potparm = inter % potparm,  &
                intra_param = intra_param, ewald_param = ewald_param )

  CALL xstring ( setup % coord_file_name, ia, ib )
  DO ibead = 1, beads

    IF ( pimd ) THEN
      CALL integer_to_string( ibead, tag )
      project_name = setup % coord_file_name ( ia : ib ) // '_' // ADJUSTL ( tag )
    ELSE
      project_name = setup % coord_file_name ( ia : ib )
    END IF

    CALL read_coord_vel ( ainp, project_name, globenv )

!..initialize box, perd
    CALL init_cell ( struc ( ibead ) % box, ainp % box, setup % simulation_cell % perd )
    struc ( ibead ) % box_ref % hmat = setup % reference_cell % hmat

! If run is a debug.  Make box_ref = box to work under
! all ensembles and restart options
    IF ( setup % run_type == 'DEBUG' ) THEN
      struc ( ibead ) % box_ref % hmat = ainp % box
    ENDIF

!..allocate memory for atoms and molecules
    CALL allocmem ( ainp, mol_setup, struc ( ibead ), globenv )

!..initialize particle_type
    CALL initialize_particle_type ( setup % atom_names, simpar, mol_setup, &
       ainp, pstat, struc ( ibead ) % part )

!..convert the units ( includes part and box ONLY when read_type is 'INIT' )
    IF ( simpar % read_type == 'INIT' ) &
    CALL convert ( units = units, part = struc  ( ibead ) % part, &
                   box = struc ( ibead ) % box )
!..calculate the inverse box matrix now after the unit conversion,
!  so it also has the right units and assign the reference box
    CALL get_hinv ( struc ( ibead ) % box )
    IF ( simpar % read_type == 'INIT' ) THEN
      struc ( ibead ) % box_ref = struc ( ibead ) % box
    ENDIF
    CALL get_hinv ( struc ( ibead ) % box_ref )

!..initialize molecule_type
    CALL initialize_molecule_type ( mol_setup, intra_param, struc ( ibead )% pnode, &
       struc ( ibead ) % part, struc ( ibead ) % molecule, globenv )

!..allocate lnhc_parameters_type for particles and get
!  number of degrees of freedom and initialize if necessary
    CALL initialize_nhc_part ( struc ( ibead ) % box, simpar, &
     struc ( ibead ) % molecule, mol_setup, globenv, &
     extended_type % nhc_part ( ibead ) )

!..allocate lnhc_parameters_type for barostat and initialize if
!  necessary
    CALL initialize_nhc_baro ( simpar, globenv, extended_type % nhc_baro ( ibead ) )

!..allocate npt_info_type and initialize if necessary
    CALL initialize_npt_type ( simpar, globenv, extended_type % npt_info, &
                             struc ( ibead ) % box )


! initialize velocities if needed
    IF ( simpar % read_type == 'POS' .OR. simpar % read_type == 'INIT' ) &
    CALL initialize_velocities ( simpar, struc ( ibead ) % part, globenv )

! Assign extended system variables
    IF ( simpar % read_type == 'ALL' ) &
    CALL assign_extended_parameters ( struc ( ibead ) % box % deth, simpar % ensemble, &
                 ainp, extended_type % npt_info, extended_type % nhc_part ( ibead ), &
                 extended_type % nhc_baro ( ibead ) )

    CALL release_mem ( ainp )

  END DO

!..K-Points
  IF ( SUM ( struc ( 1 ) % box % perd ) == 0 ) THEN
    kp % scheme = "NULL"
  ELSE 
    CALL kpoint_input ( kp, globenv )
  END IF
  
  CALL set_energy_parm ( units % pconv, units % econv, units % l_label, &
       units % vol_label, units % e_label, units % pres_label, &
       units % temp_label, units % angl_label )

  IF ( pimd ) THEN
!..initialize pimd_data (transformation related variables)
    CALL initialize_pimd ( struc , simpar )

    CALL v2ud ( struc, simpar )
  END IF

  DO i = 1, beads
    NULLIFY ( struc ( i ) % coef_pos )
    NULLIFY ( struc ( i ) % coef_vel )
    NULLIFY ( struc ( i ) % coef_force )
  ENDDO

!...initialize splines

  inter%potparm ( :, : ) % energy_cutoff = 0.0_dbl
  inter%potparm ( :, : ) % e_cutoff_coul = 0.0_dbl
  CALL spline_nonbond_control ( inter%potparm, pstat, 2500, ewald_param )

!..set linklist control parameters
  ALLOCATE ( rcut ( setup % natom_type, setup % natom_type ), STAT = isos )
  IF ( isos /=0 ) CALL stop_memory ( 'control', 'rcut', 0 )

  rcut ( :, : ) = inter%potparm ( :, : ) % rcutsq

  DO ibead = 1, beads 
    CALL set_ll_parm ( struc ( ibead ) % ll_data ( 1 ), globenv, simpar % verlet_skin, &
       setup % natom_type, rcut, simpar % n_cell )

    CALL set_ll_parm ( struc ( ibead ) % ll_data ( 1 ), globenv,  &
       printlevel = globenv % print_level, ltype = 'NONBOND' )
  END DO

  DEALLOCATE ( rcut, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'control', 'rcut' )

!
!..initialize the on-site terms for TB
  CALL tbmd_init ( struc ( 1 ) % part, inter%tbatom, inter%tbhop )
  
!
! Symmetry and K-points
!
  IF ( SUM ( struc ( 1 ) % box % perd ) == 0 ) THEN
!!!!!! symmetry setup for molecules
  ELSE
    CALL initialize_kpoints ( globenv, kp, setup % symmetry, &
         struc ( 1 ) % box % hmat, struc ( 1 ) % part )
    CALL brillouin_info ( kp, globenv )
  END IF
  
! nota bene: we use the defaults for the elctron temperature and
!             spin polarisation and assume the charge to be 0
  nel = tb_get_numel ( inter % tbatom, 0._dbl )
  CALL init_fermi_dist ( fd, REAL(nel,dbl) )
  IF ( globenv % ionode ) &
       CALL fermi_info ( fd, globenv % print_level, globenv % scr )
  CALL init_band_structure ( bs, fd, kp )
  
  CALL timestop ( zero, handle2 )

  CALL control_work ( globenv, units, simpar, struc, inter, extended_type, &
                  ewald_param, dft_control, setup, mdio, .true., .false., .false., &
                  kp, energy_min )

! deallocate memory for atoms and molecules
  DO ibead = 1, beads
    CALL deallocmem ( struc ( ibead ) )
  END DO
  DEALLOCATE ( struc, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'control', 'struc' )
  CALL timestop ( zero, handle1 )
  
END SUBROUTINE tbmd


!-----------------------------------------------------------------------------!
! CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL  !
!-----------------------------------------------------------------------------!

!******************************************************************************
SUBROUTINE release_mem ( ainp )
    TYPE(system_type)                        :: ainp

    INTEGER                                  :: isos

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
END SUBROUTINE release_mem
!******************************************************************************

SUBROUTINE control_work ( globenv, units, simpar, struc, inter, &
     extended_type, ewald_param, dft_control, setup, mdio, tbmd, fist, &
     pol, kp, energy_min )
  
  
    TYPE(global_environment_type), &
      INTENT(INOUT)                          :: globenv
    TYPE(unit_convert_type)                  :: units
    TYPE(simulation_parameters_type)         :: simpar
    TYPE(structure_type), DIMENSION(:)       :: struc
    TYPE(interaction_type)                   :: inter
    TYPE(extended_system_type)               :: extended_type
    TYPE(ewald_parameters_type)              :: ewald_param
    TYPE(dft_control_type)                   :: dft_control
    TYPE(setup_parameters_type)              :: setup
    TYPE(mdio_parameters_type)               :: mdio
    LOGICAL, INTENT(IN)                      :: tbmd, fist, pol
    TYPE(kpoint_type), INTENT(IN)            :: kp
    TYPE(energy_min_type), INTENT(INOUT)     :: energy_min

    INTEGER                                  :: handle

!------------------------------------------------------------------------------

  CALL timeset ( 'CNTL_WORK', 'I', ' ', handle )
  
  SELECT CASE ( setup % run_type )
  CASE DEFAULT
     CALL stop_program ( "control_work", "no suitable run_type" &
                         //  setup % run_type )
     
  CASE ( "DEBUG" )
     
! debug the forces
     CALL control_debug_work ( globenv, simpar, struc ( 1 ), inter, ewald_param, &
          mdio, tbmd, fist, pol )
     
  CASE ( "ENERGY" )
     
     CALL control_energy_work ( globenv, units, simpar, struc (1), inter, &
          extended_type, ewald_param, mdio, pol, dft_control, kp, energy_min )
     
  CASE ( "MD" )
     
     CALL control_md_work ( globenv, units, simpar, struc, inter, &
          extended_type, ewald_param, mdio, dft_control, kp )
     
  END SELECT
  
  CALL timestop ( zero, handle )
  
END SUBROUTINE control_work

!******************************************************************************

SUBROUTINE control_debug_work ( globenv, simpar, struc, inter, ewald_param, &
     mdio, tbmd, fist, pol )
  
  
    TYPE(global_environment_type), &
      INTENT(INOUT)                          :: globenv
    TYPE(simulation_parameters_type)         :: simpar
    TYPE(structure_type)                     :: struc
    TYPE(interaction_type)                   :: inter
    TYPE(ewald_parameters_type)              :: ewald_param
    TYPE(mdio_parameters_type)               :: mdio
    LOGICAL                                  :: tbmd, fist, pol

    TYPE(thermodynamic_type)                 :: thermo

!------------------------------------------------------------------------------
! initialize integrator

  CALL set_integrator ( globenv, mdio )
  
  CALL tbmd_debug_control ( globenv, ewald_param, struc%part, &
          struc%pnode, struc%molecule, struc%box, struc%box_ref, thermo, &
          inter%potparm, struc % ll_data ( 1 ) )
  
END SUBROUTINE control_debug_work

!******************************************************************************

SUBROUTINE control_energy_work ( globenv, units, simpar, struc, inter, &
     extended_type, ewald_param, mdio, pol, dft_control, kp, energy_min )
  
  
    TYPE(global_environment_type), &
      INTENT(INOUT)                          :: globenv
    TYPE(unit_convert_type)                  :: units
    TYPE(simulation_parameters_type)         :: simpar
    TYPE(structure_type)                     :: struc
    TYPE(interaction_type)                   :: inter
    TYPE(extended_system_type)               :: extended_type
    TYPE(ewald_parameters_type)              :: ewald_param
    TYPE(mdio_parameters_type)               :: mdio
    LOGICAL, INTENT(IN)                      :: pol
    TYPE(dft_control_type)                   :: dft_control
    TYPE(kpoint_type), INTENT(IN)            :: kp
    TYPE(energy_min_type), INTENT(INOUT)     :: energy_min

    LOGICAL, PARAMETER                       :: box_change = .FALSE.

    INTEGER                                  :: itimes
    LOGICAL                                  :: coef_flag, conv_flag
    TYPE(thermodynamic_type)                 :: thermo

!------------------------------------------------------------------------------

  coef_flag = .FALSE.
  conv_flag = .FALSE.

  coef_flag = ( ASSOCIATED ( struc % coef_pos ) .AND.  &
                ASSOCIATED ( struc % coef_vel ) .AND.  &
                ASSOCIATED ( struc % coef_force ) ) 

  itimes = 0

  IF ( pol) THEN 

    CALL force ( struc, inter, thermo, simpar, ewald_param, box_change, &
       globenv, dft_control, kp )
  ELSE 
    CALL stop_program ( "control_energy_work", "nothing in place" )
  END IF

  energy_min % min_info % e = thermo%pot

  DO itimes = 1, energy_min % iter   

     thermo%pot = 0._dbl
     CALL optimize_ao ( struc, inter, thermo, simpar, ewald_param, box_change, &
                      globenv, dft_control, kp, energy_min, itimes)

     CALL print_energy_info( itimes, energy_min % min_info, energy_min % mintype, globenv )

!   convergence acheived
    IF (energy_min % min_info % tol <= energy_min % max_tol) THEN
       conv_flag = .TRUE.
    EXIT
    END IF

  END DO

  CALL print_convergence ( conv_flag, globenv )
  
END SUBROUTINE control_energy_work

!******************************************************************************

SUBROUTINE control_md_work ( globenv, units, simpar, struc, inter, &
                 extended_type, ewald_param, mdio, dft_control, kp )
  
  
    TYPE(global_environment_type), &
      INTENT(INOUT)                          :: globenv
    TYPE(unit_convert_type)                  :: units
    TYPE(simulation_parameters_type)         :: simpar
    TYPE(structure_type), DIMENSION(:)       :: struc
    TYPE(interaction_type)                   :: inter
    TYPE(extended_system_type)               :: extended_type
    TYPE(ewald_parameters_type)              :: ewald_param
    TYPE(mdio_parameters_type)               :: mdio
    TYPE(dft_control_type)                   :: dft_control
    TYPE(kpoint_type), INTENT(IN)            :: kp

    LOGICAL, PARAMETER                       :: box_change = .FALSE.

    INTEGER                                  :: itimes
    LOGICAL                                  :: coef_flag, pimd
    REAL(dbl)                                :: cons
    TYPE(thermodynamic_type)                 :: thermo

!------------------------------------------------------------------------------
! Path integrals?

  pimd = simpar % pimd_params % beads /= 0

  coef_flag = ( ASSOCIATED ( struc ( 1 ) % coef_pos ) .AND.  &
                ASSOCIATED ( struc ( 1 ) % coef_vel ) .AND.  &
                ASSOCIATED ( struc ( 1 ) % coef_force ) ) 

! initialize integrator
  CALL set_integrator ( globenv, mdio )
  
! MD
  itimes = 0
  
  IF ( pimd ) THEN 
    CALL force ( struc, inter, thermo, simpar, ewald_param, box_change, &
       globenv, dft_control, kp )
  ELSE
    CALL force ( struc ( 1 ), inter, thermo, simpar, ewald_param, box_change, &
       globenv, dft_control, kp )
  END IF

  IF ( globenv % ionode .AND. ewald_param % ewald_type /= 'NONE' ) &
       CALL ewald_print ( globenv % scr, thermo, struc ( 1 ) % box, &
       units % e_label )
  
  IF ( pimd ) CALL fr2fu ( struc, simpar )

  IF ( coef_flag ) THEN
    CALL energy ( itimes, cons, simpar, struc ( 1 ), thermo, extended_type % nhc_part ( 1 ), &
       extended_type % nhc_baro ( 1 ), extended_type % nhc_pol, extended_type % npt_info )

  ELSEIF ( pimd ) THEN
    CALL energy ( itimes, cons, simpar, struc, thermo, extended_type % nhc_part, &
       extended_type % nhc_baro, extended_type % npt_info )
  ELSE
    CALL energy ( itimes, cons, simpar, struc ( 1 ), thermo, extended_type % nhc_part ( 1 ), &
       extended_type % nhc_baro ( 1 ), extended_type % npt_info )
  ENDIF

  DO itimes = 1, simpar % nsteps
     IF ( coef_flag ) THEN
       CALL velocity_verlet ( itimes, cons, simpar, inter, thermo, &
          struc ( 1 ), ewald_param, extended_type % nhc_part ( 1 ), &
          extended_type % nhc_baro ( 1 ), extended_type % nhc_pol, &
          extended_type % npt_info, dft_control, kp )
     ELSE
       CALL velocity_verlet ( itimes, cons, simpar, inter, thermo, &
          struc, ewald_param, extended_type % nhc_part, extended_type % nhc_baro, &
          extended_type % npt_info, dft_control, kp )
     ENDIF
     
     IF ( MOD ( itimes, mdio % idump ) == 0 ) THEN
          IF ( pimd ) CALL ud2v ( struc, simpar )
          CALL dump_variables ( struc, extended_type % nhc_part, &
               extended_type % nhc_baro, extended_type % npt_info, & 
               mdio % dump_file_name, globenv, pimd )
          IF ( pimd ) CALL v2ud ( struc, simpar )
     END IF
  END DO
  
  IF ( pimd ) CALL ud2v ( struc, simpar )

  CALL dump_variables ( struc, extended_type % nhc_part, &
       extended_type % nhc_baro, extended_type % npt_info, & 
       mdio % dump_file_name, globenv, pimd )
  
  IF(pimd ) CALL v2ud ( struc, simpar )

  IF ( globenv % ionode ) CALL close_unit ( 10, 99 )
  
END SUBROUTINE control_md_work

!******************************************************************************

SUBROUTINE allocmem ( ainp, mol_setup, struc, globenv )


    TYPE(system_type), INTENT(IN)            :: ainp
    TYPE(molecule_type), DIMENSION(:), &
      INTENT(IN)                             :: mol_setup
    TYPE(structure_type), INTENT(INOUT)      :: struc
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    INTEGER                                  :: i, ios, iw, natoms, nmol, &
                                                nmoltype, nnodes, nsh

!------------------------------------------------------------------------------

  struc % name = globenv % program_name // ' MOLECULAR SYSTEM'
  
  IF ( globenv % num_pe == 1 ) THEN
     natoms = SIZE ( ainp % c ,2 )
     ALLOCATE ( struc % part ( natoms ), STAT = ios )
     IF ( ios /= 0 ) CALL stop_memory ( 'control', 'part', natoms )
     ALLOCATE ( struc % pnode ( natoms ), STAT = ios )
     IF ( ios /= 0 ) CALL stop_memory ( 'control', 'pnode', natoms )
     nmol = SUM ( mol_setup ( : ) % num_mol )

     ALLOCATE ( struc % molecule ( nmol ), STAT = ios )
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
     natoms = SIZE ( ainp % c ,2 )
     ALLOCATE ( struc % part ( natoms ), STAT = ios )
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

     ALLOCATE ( struc % molecule ( nmol ), STAT = ios )
     IF ( ios /= 0 ) CALL stop_memory ( 'control', 'molecule' , nmol )
     ALLOCATE ( struc % pnode ( nnodes ), STAT = ios )
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

    TYPE(structure_type), INTENT(INOUT)      :: struc

    INTEGER                                  :: ios

!------------------------------------------------------------------------------

  DEALLOCATE ( struc % part, STAT = ios )
  IF ( ios /= 0 ) CALL stop_memory ( 'control', 'part' )

  DEALLOCATE ( struc % pnode, STAT = ios )
  IF ( ios /= 0 ) CALL stop_memory ( 'control', 'pnode' )

  DEALLOCATE ( struc % molecule, STAT = ios )
  IF ( ios /= 0 ) CALL stop_memory ( 'control', 'molecule' )

END SUBROUTINE deallocmem

!******************************************************************************

END MODULE tbmd_module
