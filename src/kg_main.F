!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/kg_main [1.0] *
!!
!!   NAME
!!     kg_main
!!
!!   FUNCTION
!!     perform kim-gordon moecular dynamics 
!!   AUTHOR
!!     gt SEPT-23-2002
!!   SOURCE
!******************************************************************************

MODULE kg_main 

  USE atomic_kind_types, ONLY: atomic_kind_type,&
                                             get_atomic_kind_set,&
                                             init_atomic_kind_set,&
                                             read_atomic_kind_set,&
                                             write_atomic_kind_set,&
                                             write_gto_basis_sets
  USE atoms_input, ONLY : read_coord_vel, system_type
  USE cntl_input, ONLY : read_cntl_section
  USE convert_units, ONLY : convert
  USE dft_types, ONLY : dft_control_type, read_dft_control, read_kg_control, &
                        write_dft_control, write_kg_control
  USE dump, ONLY : dump_variables
  USE ewalds, ONLY : ewald_print
  USE ewald_parameters_types, ONLY : ewald_parameters_type 
  USE extended_system_types, ONLY : extended_system_type
  USE force_control, ONLY : set_force_env, force
  USE force_fields, ONLY : read_force_field_section
  USE global_types, ONLY : global_environment_type
  USE initialize_extended_types, ONLY : initialize_nhc_part,  &
                                        initialize_nhc_baro,  &
                                        initialize_npt_type,  &
                                        assign_extended_parameters
  USE initialize_molecule_types, ONLY : initialize_molecule_type
  USE initialize_pimd_types, ONLY : initialize_pimd
  USE input_types, ONLY : setup_parameters_type
  USE integrator, ONLY : set_integrator
  USE kinds, ONLY : dbl
  USE kg_environment_types, ONLY : kg_environment_type,  &
                                     kg_thermodynamic_type, &
                                     set_kg_env, init_kg_env,               &
                                     kg_replica_environment_type
!gt not yet written!
!  USE kg_debug, ONLY : kg_debug_control => debug_control
!gt 
  USE linklist_control, ONLY : set_ll_parm
  USE mathconstants, ONLY : zero, twopi
  USE md, ONLY : read_md_section, simulation_parameters_type, &
       initialize_velocities, mdio_parameters_type, virial_type
!gt not yet written!
!  USE md_fist_energies, ONLY : set_fist_energies, set_fist_energies_param, &
!                          md_energy_fist, md_energy_fist_pimd
!gt
  USE md_environment_types, ONLY :  set_md_env, init_md_env, & 
                                    md_environment_type, zero_virial
  USE molecule_input, ONLY : read_molecule_section, read_setup_section, &
       read_topology_section
  USE molecule_types, ONLY : molecule_type, intra_parameters_type,&
       topology_type
  USE orbital_transformation_matrices, ONLY: init_spherical_harmonics
  USE pair_potential, ONLY : potentialparm_type, spline_nonbond_control
  USE particle_types, ONLY : particle_prop_type, particle_type, init_particle_set
  USE physcon, ONLY : boltzmann, h_planck
  USE qs_interactions,                 ONLY: init_interaction_radii
  USE read_pimd, ONLY : read_pimd_section
  USE simulation_cell, ONLY :  cell_type, get_hinv, init_cell
  USE string_utilities, ONLY : xstring, integer_to_string
  USE structure_types, ONLY : structure_type, init_structure_type,  &
                              set_structure_type
  USE termination, ONLY : stop_program, stop_memory
  USE transformations, ONLY : fr2fu, v2ud, ud2v
  USE timings, ONLY : timeset, timestop, trace_debug
  USE unit, ONLY : unit_convert_type, set_units
  USE util, ONLY : close_unit, get_share
  USE velocity_verlet_control, ONLY : velocity_verlet

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: kg
  
!!*****
!-----------------------------------------------------------------------------!

CONTAINS

!-----------------------------------------------------------------------------!
! KG KG KG                                                                    !
!-----------------------------------------------------------------------------!
!!****** kg_main/kg [1.0] *
!!
!!   NAME
!!     kg
!!
!!   FUNCTION
!!     Controls program flow for Kim-Gordon type MD
!!
!!   AUTHOR
!!     GT-CJM
!!
!!   MODIFICATION HISTORY
!!
!!   USED BY
!!     cp2k
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE kg ( globenv )

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
  INTEGER :: ikind,istat,lmax,maxl,maxlppl,maxlppnl,nkind,output_unit
  LOGICAL :: ionode

  CHARACTER ( LEN = 40 ) :: set_fn, project_name
  CHARACTER ( LEN = 5 ) :: tag
  LOGICAL, POINTER :: pimd
  
  TYPE ( atomic_kind_type ), DIMENSION(:), POINTER:: atomic_kind_set
  TYPE ( dft_control_type ), POINTER :: dft_control
  TYPE ( extended_system_type ) :: extended_type
  TYPE ( ewald_parameters_type ), POINTER :: ewald_param
  TYPE ( kg_environment_type ) :: kg_env
  TYPE ( intra_parameters_type ) :: intra_param
  TYPE ( md_environment_type ) :: md_env
  TYPE ( mdio_parameters_type ) :: mdio
  TYPE ( molecule_type ), DIMENSION ( : ), POINTER :: mol_setup
  TYPE ( particle_prop_type ), DIMENSION ( : ), POINTER :: pstat
  TYPE ( particle_type ), DIMENSION(:), POINTER    :: particle_set
  TYPE ( potentialparm_type ), DIMENSION ( :, : ), POINTER :: potparm
  TYPE ( kg_replica_environment_type ), DIMENSION ( : ), POINTER :: rep_env
  TYPE ( simulation_parameters_type ), POINTER :: simpar
  TYPE ( setup_parameters_type ) :: setup
  TYPE ( system_type ) :: ainp
  TYPE ( kg_thermodynamic_type ), POINTER :: thermo
  TYPE ( topology_type ) :: topo
  TYPE ( unit_convert_type ) :: units
  TYPE ( virial_type ), POINTER :: virial

  
!------------------------------------------------------------------------------

! IF( globenv % ionode ) CALL trace_debug ( "START" )

  CALL timeset ( 'KG_MAIN', 'I', ' ', handle1 )
  CALL timeset ( 'KG_INIT', 'I', ' ', handle2 )
  
! initialize kg_env
  CALL init_kg_env ( kg_env )

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
  NULLIFY ( atomic_kind_set )
  NULLIFY ( dft_control )
  NULLIFY ( particle_set )
  NULLIFY ( potparm )
  NULLIFY ( rep_env )
! initialize the virial
  CALL zero_virial ( virial )
! point from kg_thermo % virial to virial
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

  CALL read_dft_control ( dft_control, globenv )

  CALL write_dft_control ( dft_control, globenv )

  CALL read_kg_control ( dft_control % kg_control, globenv )

  CALL write_kg_control ( dft_control % kg_control, globenv )

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
    rep_env ( ibead ) % ll_data ( 1 ) % natom_types = 0
    rep_env ( ibead ) % ll_data ( 1 ) % list_type = 0
    rep_env ( ibead ) % ll_data ( 1 ) % counter = 0
    rep_env ( ibead ) % ll_data ( 1 ) % last_update = 0
    rep_env ( ibead ) % ll_data ( 1 ) % num_update = 0
    rep_env ( ibead ) % ll_data ( 1 ) % print_level = 1
    NULLIFY ( rep_env ( ibead ) % ll_data ( 1 ) % r_last_update )
    NULLIFY ( rep_env ( ibead ) % ll_data ( 1 ) % rlist_cut )
    NULLIFY ( rep_env ( ibead ) % ll_data ( 1 ) % rlist_cutsq )
    NULLIFY ( rep_env ( ibead ) % ll_data ( 1 ) % pp_images )
    NULLIFY ( rep_env ( ibead ) % ll_data ( 1 ) % pp_ncell )
    NULLIFY ( rep_env ( ibead ) % ll_data ( 1 ) % pp_startlist_im )
    NULLIFY ( rep_env ( ibead ) % ll_data ( 1 ) % pp_startlist_nl )
    NULLIFY ( rep_env ( ibead ) % ll_data ( 1 ) % pp_neighbor )
    NULLIFY ( rep_env ( ibead ) % ll_data ( 1 ) % pp_startlist_cell )
    NULLIFY ( rep_env ( ibead ) % ll_data ( 1 ) % pp_cell_ll )
    NULLIFY ( rep_env ( ibead ) % atomic_kind_set )
  END DO

!  IF ( globenv % ionode ) THEN
!    CALL fist_header ( globenv % scr )
!  END IF
  

! read from the &SETUP and &MOLECULE section of *.set

  set_fn = setup % set_file_name

  CALL read_setup_section ( mol_setup, setup, globenv )

  CALL read_molecule_section ( mol_setup, setup, globenv )

  CALL read_topology_section ( topo, setup, globenv )

! read force_field information for classical MD in *.set

  CALL read_force_field_section ( setup, mol_setup, set_fn, &
                        intra_param, potparm, pstat, globenv )

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
                pstat = pstat, potparm = potparm,  &
                intra_param = intra_param, ewald_param = ewald_param )

  CALL xstring ( setup % coord_file_name, ia, ib )

! Filling up the replica environment one bead at a time
  DO ibead = 1, beads
   

    IF ( pimd ) THEN
      CALL integer_to_string( ibead, tag )
      project_name = setup % coord_file_name ( ia : ib ) // '_' // ADJUSTL ( tag )
    ELSE
      project_name = setup % coord_file_name ( ia : ib )
    END IF

    CALL read_coord_vel ( particle_set, atomic_kind_set,&
                           ainp, project_name, mol_setup, &
                           setup%atom_names, pstat, globenv )
! Read the atomic kind set
    CALL read_atomic_kind_set(atomic_kind_set,globenv)

! Write the atomic kind set
    CALL write_gto_basis_sets(atomic_kind_set,globenv)

! Initialize the spherical harmonics
! the orbital transformation matrices
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxl=maxl,&
                             maxlppl=maxlppl,&
                             maxlppnl=maxlppnl)

    lmax = MAX(maxl,maxlppl,maxlppnl)

    CALL init_spherical_harmonics(lmax,globenv)

! Initialise the atomic kind set
    CALL init_atomic_kind_set(atomic_kind_set,globenv)

! Print the atomic kind set
    CALL write_atomic_kind_set(atomic_kind_set,globenv)
! Initialize the atomic data sets completely ***
    CALL init_particle_set(particle_set)

!..initialize box, perd
!  set defaults for working units (a.u)
    rep_env(ibead)%box%unit_of_length_name = "BOHR"
    rep_env(ibead)%box%unit_of_length = 1.0_dbl
    rep_env(ibead)%box%scaled_coordinates = .FALSE.

    IF (ainp%n > 0) THEN
      CALL init_cell ( rep_env (ibead) % box, ainp % box, &
                       setup % simulation_cell % perd )
    ELSE 
      CALL init_cell ( rep_env (ibead) % box, setup % simulation_cell % hmat, &
                       setup % simulation_cell % perd )
    END IF
    rep_env ( ibead ) % box_ref % hmat = setup % reference_cell % hmat

! If run is a debug.  Make box_ref = box to work under
! all ensembles and restart options
    IF ( setup % run_type == 'DEBUG' ) THEN
      IF (ainp%n > 0) THEN
         rep_env ( ibead ) % box_ref % hmat = ainp % box
      ELSE 
         rep_env ( ibead ) % box_ref % hmat = setup % simulation_cell % hmat
      END IF
    ENDIF

    CALL init_interaction_radii(dft_control%kg_control, rep_env (ibead) % box, &
                                atomic_kind_set,globenv)

!..allocate memory for atoms and molecules
    CALL allocmem ( ainp, mol_setup, rep_env ( ibead ), globenv )

!..convert the units ( includes part and box ONLY when read_type is 'INIT' )
    IF ( simpar % read_type == 'INIT' ) &
    CALL convert ( units = units, part = particle_set, &
                   box = rep_env (ibead) % box )
!..calculate the inverse box matrix now after the unit conversion,
!  so it also has the right units and assign the reference box
    CALL get_hinv (  rep_env (ibead) % box )
    IF ( simpar % read_type == 'INIT' ) THEN
       rep_env ( ibead ) % box_ref = rep_env (ibead) % box
    ENDIF
    CALL get_hinv ( rep_env ( ibead ) % box_ref )

    rep_env (ibead) % atomic_kind_set => atomic_kind_set
    rep_env (ibead) % part=> particle_set

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
                               rep_env (ibead) % box)


! initialize velocities if needed
    IF ( simpar % read_type == 'POS' .OR. simpar % read_type == 'INIT' ) &
    CALL initialize_velocities ( simpar, particle_set, globenv )

! Assign extended system variables
    IF ( simpar % read_type == 'ALL' ) &
    CALL assign_extended_parameters ( rep_env ( ibead ) % box % deth,    &
                                      simpar % ensemble, ainp,           &
                                      extended_type % npt_info,           &
                                      extended_type % nhc_part ( ibead ), &
                                      extended_type % nhc_baro ( ibead ) )

    CALL release_mem ( ainp )

  END DO

!gt not yet written !!
!  CALL set_kg_energies_param ( units % pconv, units % econv, units % l_label, &
!       units % vol_label, units % e_label, units % pres_label, &
!       units % temp_label, units % angl_label )
!gt


!...initialize splines

!
  CALL timestop ( zero, handle2 )

! set the fist_env
  CALL set_kg_env ( kg_env, pimd, dft_control, ewald_param,&
                      rep_env, potparm, thermo )

! begin the MD process
!  gt not yet written
!  CALL control_work ( globenv, units, simpar, kg_env, md_env, extended_type, &
!                      setup, mdio )

! deallocate memory for atoms and molecules
  DO ibead = 1, beads
    CALL deallocmem ( rep_env ( ibead ) )
  END DO

END SUBROUTINE kg

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
     IF ( isos /= 0 ) CALL stop_memory ( 'kg_main', 'ainp%c' )
  END IF

  IF ( ASSOCIATED ( ainp % v ) ) THEN
     DEALLOCATE ( ainp % v, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'kg_main', 'ainp%v' )
  END IF

  IF ( ASSOCIATED ( ainp % eta_part ) ) THEN
     DEALLOCATE ( ainp % eta_part, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'kg_main', 'ainp%eta_part' )
  END IF

  IF ( ASSOCIATED ( ainp % veta_part ) ) THEN
     DEALLOCATE ( ainp % veta_part, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'kg_main', 'ainp%veta_part' )
  END IF

  IF ( ASSOCIATED ( ainp % mass_eta_part ) ) THEN
     DEALLOCATE ( ainp % mass_eta_part, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'kg_main', 'ainp%mass_eta_part' )
  END IF

  IF ( ASSOCIATED ( ainp % eta_baro ) ) THEN
     DEALLOCATE ( ainp % eta_baro, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'kg_main', 'ainp%eta_baro' )
  END IF

  IF ( ASSOCIATED ( ainp % veta_baro ) ) THEN
     DEALLOCATE ( ainp % veta_baro, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'kg_main', 'ainp%veta_baro' )
  END IF

  IF ( ASSOCIATED ( ainp % mass_eta_baro ) ) THEN
     DEALLOCATE ( ainp % mass_eta_baro, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'kg_main', 'ainp%mass_eta_baro' )
  END IF

  IF ( ASSOCIATED ( ainp % veps ) ) THEN
     DEALLOCATE ( ainp % veps, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'kg_main', 'ainp%veps' )
  END IF

  IF ( ASSOCIATED ( ainp % mass_eps ) ) THEN
     DEALLOCATE ( ainp % mass_eps, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'kg_main', 'ainp%mass_eps' )
  END IF 
END SUBROUTINE release_mem
!******************************************************************************

SUBROUTINE allocmem ( ainp, mol_setup, rep_env, globenv )

  IMPLICIT NONE

! Arguments
  TYPE ( molecule_type ), DIMENSION ( : ), INTENT ( IN ) :: mol_setup
  TYPE ( kg_replica_environment_type ), INTENT ( INOUT ) :: rep_env
  TYPE ( system_type ), INTENT ( IN ) :: ainp
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv

! Locals
  INTEGER :: iw, natoms, nnodes, nmol, nmoltype, ios, iat, i, nsh

!------------------------------------------------------------------------------

  IF ( globenv % num_pe == 1 ) THEN
     IF ( ainp % n > 0 ) THEN
        natoms = SIZE ( ainp % c ( 1, : ) )
     ELSE
        natoms = SIZE ( rep_env % part ( : ) )
     END IF

     ALLOCATE ( rep_env % pnode ( natoms ), STAT = ios )
     IF ( ios /= 0 ) CALL stop_memory ( 'kg_main', 'pnode', natoms )

     nmol = SUM ( mol_setup ( : ) % num_mol )
     ALLOCATE ( rep_env % molecule ( nmol ), STAT = ios )
     IF ( ios /= 0 ) CALL stop_memory ( 'kg_main', 'molecule', nmol )

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
        natoms = SIZE ( ainp % c ( 1, : ) )
     ELSE
        natoms = SIZE ( rep_env % part ( : ) )
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
     IF ( ios /= 0 ) CALL stop_memory ( 'kg_main', 'molecule' , nmol )
     ALLOCATE ( rep_env % pnode ( nnodes ), STAT = ios )
     IF ( ios /= 0 ) CALL stop_memory ( 'kg_main', 'pnode', nnodes )

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
  TYPE ( kg_replica_environment_type ), INTENT ( INOUT ) :: rep_env

! Locals
  INTEGER :: ios

!------------------------------------------------------------------------------

  DEALLOCATE ( rep_env % pnode, STAT = ios )
  IF ( ios /= 0 ) CALL stop_memory ( 'kg_main', 'pnode' )

  DEALLOCATE ( rep_env % molecule, STAT = ios )
  IF ( ios /= 0 ) CALL stop_memory ( 'kg_main', 'molecule' )

END SUBROUTINE deallocmem

!******************************************************************************

END MODULE kg_main
!******************************************************************************

!******************************************************************************
