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

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set,&
                                             init_atomic_kind_set,&
                                             read_atomic_kind_set,&
                                             write_atomic_kind_set,&
                                             write_gto_basis_sets
  USE atoms_input,                     ONLY: read_atoms_input,&
                                             read_md_restart
  USE cntl_input,                      ONLY: read_cntl_section
  USE coefficient_types,               ONLY: coeff_type
  USE cp2k_units,                      ONLY: init_cp2k_units
  USE dft_types,                       ONLY: dft_control_type,&
                                             read_dft_control,&
                                             read_kg_control,&
                                             write_dft_control,&
                                             write_kg_control
  USE dump,                            ONLY: dump_variables
  USE extended_system_dynamics,        ONLY: lnhc_parameters_type,&
                                             npt_info_type
  USE extended_system_types,           ONLY: extended_system_type
  USE force_env_types,                   ONLY: force_env_calculate_force,&
                                             force_env_create,&
                                             force_env_release,&
                                             force_env_type
  USE force_fields,                    ONLY: read_force_field_section
  USE global_types,                    ONLY: global_environment_type
  USE initialize_extended_types,       ONLY: initialize_nhc_baro,&
                                             initialize_nhc_coef,&
                                             initialize_nhc_forces,&
                                             initialize_nhc_part,&
                                             initialize_npt_type
  USE initialize_molecule_types,       ONLY: initialize_molecule_type
  USE input_types,                     ONLY: setup_parameters_type
  USE integrator,                      ONLY: set_integrator
  USE kg_debug,                        ONLY: kg_debug_control
  USE kg_energy_optimize,              ONLY: optimize_ao
  USE kg_energy_utils,                 ONLY: energy_min_type,&
                                             min_info_type,&
                                             print_convergence,&
                                             print_energy_info,&
                                             read_energy_section
  USE kg_environment_types,            ONLY: init_kg_env,&
                                             kg_environment_type,&
                                             kg_thermodynamic_type,&
                                             set_kg_env
  USE kg_force,                        ONLY: kg_force_release
  USE kg_pol_coefs,                    ONLY: allocate_coefs,&
                                             dump_coefs,&
                                             initialize_coefs
  USE kinds,                           ONLY: dbl
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
  USE md_kg_energies,                  ONLY: md_energy_kg,&
                                             set_kg_energies
  USE molecule_input,                  ONLY: read_molecule_section,&
                                             read_setup_section,&
                                             read_topology_section
  USE molecule_types,                  ONLY: intra_parameters_type,&
                                             molecule_structure_type,&
                                             molecule_type,&
                                             particle_node_type,&
                                             topology_type
  USE orbital_transformation_matrices, ONLY: init_spherical_harmonics
  USE pair_potential_types,            ONLY: pair_potential_type
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: boltzmann,&
                                             h_planck
  USE qs_interactions,                 ONLY: init_interaction_radii
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_hinv,&
                                             init_cell,&
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
  USE util,                            ONLY: close_unit,&
                                             get_share
  USE velocity_verlet_control,         ONLY: velocity_verlet
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
!!   SYNOPSIS
!!     Subroutine kg(globenv)
!!       Implicit None
!!       Type(global_environment_type), Intent (INOUT):: globenv
!!     End Subroutine kg
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
!!*** *************************************************************************

SUBROUTINE kg ( globenv )

  IMPLICIT NONE

! Argument
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv

! Locals
  REAL ( dbl ), DIMENSION ( :, : ), ALLOCATABLE :: rcut
  INTEGER :: handle1, handle2, isos
  INTEGER :: ia, ib
  INTEGER :: n, nll, i
  INTEGER :: maxder
  INTEGER :: ncoef
  INTEGER :: ikind,istat,maxl,maxlgto,maxlppl,maxlppnl,nkind,output_unit
  LOGICAL :: ionode, kgpol

  CHARACTER ( LEN = 40 ) :: set_fn, project_name
  CHARACTER ( LEN = 5 ) :: tag
  
  TYPE ( atomic_kind_type ), DIMENSION(:), POINTER:: atomic_kind_set
  TYPE ( cell_type ), POINTER ::  box, box_ref
  TYPE ( coeff_type ), DIMENSION (:), POINTER :: coef_pos
  TYPE ( coeff_type ), DIMENSION (:), POINTER :: coef_vel
  TYPE ( coeff_type ), DIMENSION (:), POINTER :: coef_force
  TYPE ( dft_control_type ), POINTER :: dft_control
  TYPE ( extended_system_type ) :: extended_type
  TYPE ( energy_min_type ) :: energy_min
  TYPE ( intra_parameters_type ) :: intra_param
  TYPE ( kg_environment_type ), POINTER :: kg_env
  TYPE ( lnhc_parameters_type ), POINTER :: nhc_part, nhc_baro
  TYPE ( md_environment_type ) :: md_env
  TYPE ( mdio_parameters_type ) :: mdio
  TYPE ( molecule_type ), DIMENSION ( : ), POINTER :: mol_setup
  TYPE ( molecule_structure_type ), DIMENSION ( : ), POINTER :: molecule
  TYPE ( npt_info_type ), DIMENSION ( :, : ), POINTER :: npt_info
  TYPE ( particle_type ), DIMENSION(:), POINTER    :: particle_set
  TYPE ( particle_node_type ), DIMENSION(:), POINTER    :: pnode
  TYPE ( pair_potential_type ), DIMENSION ( :, : ), POINTER :: potparm
  TYPE ( simulation_parameters_type ), POINTER :: simpar
  TYPE ( setup_parameters_type ) :: setup
  TYPE ( kg_thermodynamic_type ), POINTER :: thermo
  TYPE ( topology_type ) :: topo
  TYPE ( virial_type ), POINTER :: virial

!------------------------------------------------------------------------------

! IF( globenv % ionode ) CALL trace_debug ( "START" )

  CALL timeset ( 'KG_INIT', 'I', ' ', handle2 )
  
  CALL init_cp2k_units ()

  NULLIFY ( kg_env )
  ALLOCATE ( kg_env, STAT = isos )
  IF ( isos /=0 ) CALL stop_memory ( 'allocstruc', 'kg_env', 0 )
  
! initialize kg_env
  CALL init_kg_env ( kg_env )

! allocating or nullifying local pointers
! potparm
  NULLIFY ( potparm )
! simpar
  ALLOCATE ( simpar, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'allocstruc', 'simpar', 0 )
! thermodynamics
  ALLOCATE ( thermo, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'allocstruc', 'thermo', 0 )
! virial
  ALLOCATE ( virial, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'allocstruc', 'virial', 0 )
! nullify the rest of the pointers
  NULLIFY ( atomic_kind_set )
  NULLIFY ( coef_pos, coef_vel, coef_force )
  NULLIFY ( dft_control )
  NULLIFY ( particle_set, pnode, molecule )
  NULLIFY ( nhc_part )
  NULLIFY ( nhc_baro )
  NULLIFY ( npt_info )
  NULLIFY ( mol_setup )
  NULLIFY ( box )
  NULLIFY ( box_ref )

! initialize the virial
  CALL zero_virial ( virial )
! point from kg_thermo % virial to virial
  thermo % virial => virial

  kgpol = .FALSE.

! read control section
  CALL read_cntl_section ( setup, globenv )

  CALL read_dft_control ( dft_control, globenv )

  CALL write_dft_control ( dft_control, globenv )

  CALL read_kg_control ( dft_control % kg_control, globenv )

  CALL write_kg_control ( dft_control % kg_control, globenv )
   
  kgpol = dft_control % kg_control % polarization

! Allocate particle thermostat for general use with path-integrals
  ALLOCATE ( extended_type % nhc_part ( 1 ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'allocstruc', 'nhc_part', 1 )

! Allocate barostat thermostat for general use with path-integrals
  ALLOCATE ( extended_type % nhc_baro ( 1 ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'allocstruc', 'nhc_baro', 1 )

! read from the &SETUP and &MOLECULE section of *.set

  set_fn = setup % set_file_name

  CALL read_setup_section ( mol_setup, setup, globenv )

  CALL read_molecule_section ( mol_setup, setup, globenv )

  CALL read_topology_section ( topo, setup, globenv )
   
! read force_field information for classical MD in *.set

  CALL read_force_field_section ( setup, set_fn, intra_param,  &
                                  potparm, globenv )

! read the input of the molecular dynamics section
  CALL read_md_section ( simpar, globenv, mdio )
  simpar % program = globenv % program_name

  CALL xstring ( setup % coord_file_name, ia, ib )

  project_name = setup % coord_file_name ( ia : ib )

!*** Read the input section with the cell parameters ***

  CALL read_cell ( box, globenv )

!*** Print the cell parameters ***

  CALL write_cell ( box, globenv )

!*** Read &COORD section ***
  CALL read_atoms_input ( box % unit_of_length_name, particle_set, &
                          atomic_kind_set, mol_setup, setup, globenv )


! Read the atomic kind set
  CALL read_atomic_kind_set(atomic_kind_set,globenv)

! Write the atomic kind set
  CALL write_gto_basis_sets(atomic_kind_set,globenv)

! Initialize the spherical harmonics
! the orbital transformation matrices
  CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxlgto=maxlgto,&
                             maxlppl=maxlppl,&
                             maxlppnl=maxlppnl)

  maxl = MAX(maxlgto,maxlppl,maxlppnl)

  CALL init_spherical_harmonics(maxl,globenv)

! Initialise the atomic kind set
  CALL init_atomic_kind_set(atomic_kind_set,globenv)

! Print the atomic kind set
  CALL write_atomic_kind_set(atomic_kind_set,globenv)

!..initialize cell parameters
    box_ref => box

  CALL init_interaction_radii ( dft_control%kg_control, box, &
                                atomic_kind_set, globenv)

!..allocate memory for pnodes and molecules
  CALL allocmem ( mol_setup, particle_set, pnode, molecule, globenv )

! If polarization is also present, allocate and initialize coefs pointer  

  IF ( kgpol ) THEN

     ALLOCATE ( coef_pos ( 1 ) , STAT = isos )
     IF ( isos /=0 ) CALL stop_memory ( 'kg', 'coef_pos', 1 )
     ALLOCATE ( coef_vel ( 1 ) , STAT = isos )
     IF ( isos /=0 ) CALL stop_memory ( 'kg', 'coef_vel', 1 )
     ALLOCATE ( coef_force ( 1 ) , STAT = isos )
     IF ( isos /=0 ) CALL stop_memory ( 'kg', 'coef_force', 1 )

!    Allocate coef thermostats
     ALLOCATE ( extended_type % nhc_coef, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'allocstruc', 'nhc_coef', 0 )


     CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                              ncgf_aux=ncoef)

     CALL allocate_coefs (ncoef, coef_pos ( 1 ), coef_vel ( 1 ), &
                              coef_force ( 1 ) )

     CALL initialize_coefs (simpar, coef_pos ( 1 ), coef_vel ( 1 ),&
                              coef_force ( 1 ), globenv )

     CALL initialize_nhc_coef ( simpar,  coef_vel ( 1 ),  &
                                globenv, extended_type % nhc_coef )

     IF ( globenv % run_type == 'WAVEFUNCTION OPTIMIZATION') THEN
       CALL read_energy_section ( energy_min, globenv )
     END IF

  END IF

!..initialize molecule_type
  CALL initialize_molecule_type ( mol_setup, pnode,  particle_set,  &
                                  molecule, globenv, intra_param) 

!..allocate lnhc_parameters_type for particles and get
!  number of degrees of freedom and initialize if necessary
  CALL initialize_nhc_part ( box, simpar, molecule, mol_setup, &
                               globenv, extended_type % nhc_part (1) ) 

!..allocate lnhc_parameters_type for barostat and initialize if
!  necessary
  CALL initialize_nhc_baro ( simpar, globenv, extended_type % nhc_baro (1) )

!..allocate npt_info_type and initialize if necessary
  CALL initialize_npt_type ( simpar, globenv, extended_type % npt_info, box )

! initialize velocities if needed
  npt_info => extended_type % npt_info 
  nhc_part => extended_type % nhc_part ( 1 )
  nhc_baro => extended_type % nhc_baro ( 1 ) 
  SELECT CASE ( simpar % read_type )
    CASE ( 'INIT' )
      CALL initialize_velocities ( simpar, particle_set, globenv )
    CASE ( 'POS' )
      CALL read_md_restart ( simpar % read_type, globenv, particle_set,    &
                             box, box_ref, nhc_part, nhc_baro, npt_info )
      CALL initialize_velocities ( simpar, particle_set, globenv )
    CASE DEFAULT
      CALL read_md_restart ( simpar % read_type, globenv, particle_set,    &
                             box, box_ref, nhc_part, nhc_baro, npt_info )
  END SELECT
!   CALL get_hinv ( box )
!   CALL get_hinv ( box_ref )

  CALL initialize_nhc_forces ( nhc_part, nhc_baro )

  CALL timestop ( zero, handle2 )

! set the fist_env
  IF ( kgpol ) THEN 
    CALL set_kg_env ( kg_env=kg_env, atomic_kind_set=atomic_kind_set, box=box, &
                        box_ref=box_ref, coef_pos=coef_pos, coef_vel=coef_vel, &
                        coef_force=coef_force, dft_control=dft_control, &
                        molecule=molecule, part=particle_set, &
                        pnode=pnode, thermo=thermo )  
  ELSE
    CALL set_kg_env ( kg_env=kg_env, atomic_kind_set=atomic_kind_set, box=box, &
                        box_ref=box_ref, dft_control=dft_control, &
                        molecule=molecule, part=particle_set, &
                        pnode=pnode, thermo=thermo )  
  END IF
    
! begin the MD process
  CALL control_work ( globenv, simpar, kg_env, md_env, extended_type, &
                      setup, mdio, energy_min, kgpol )

! deallocate memory for atoms and molecules
  CALL deallocmem ( pnode, molecule )
! gt: to be done:
! dealocation of the coef pointer if kgpol

END SUBROUTINE kg

!-----------------------------------------------------------------------------!
! KG KG KG KG KG KG KG KG KG KG KG KG KG KG KG  !
!-----------------------------------------------------------------------------!
!******************************************************************************

SUBROUTINE control_work ( globenv, simpar, kg_env, md_env, &
                          extended_type, setup, mdio, energy_min, kgpol )

  IMPLICIT NONE

! Arguments
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
  TYPE ( simulation_parameters_type ), POINTER :: simpar
  TYPE ( kg_environment_type ), POINTER :: kg_env
  TYPE ( md_environment_type ), INTENT ( INOUT ) :: md_env
  TYPE ( extended_system_type ) :: extended_type
  TYPE ( setup_parameters_type ) :: setup
  TYPE ( mdio_parameters_type ) :: mdio
  TYPE ( energy_min_type ), INTENT (inout) :: energy_min
  LOGICAL, INTENT (in) :: kgpol

! Locals
  INTEGER :: handle, isos
  TYPE ( structure_type ), DIMENSION (:), POINTER :: struc 
  TYPE ( force_env_type ), POINTER :: force_env

!------------------------------------------------------------------------------

  CALL timeset ( 'CNTL_WORK', 'I', ' ', handle )

  NULLIFY(force_env)
! allocate structure_type
  ALLOCATE ( struc (1) , STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'allocstruc', 'struc', 1 )

! set up structure_type
  CALL init_structure_type ( struc (1) )
  CALL set_structure_type ( struc=struc(1), kg_env=kg_env )

! create the force environment
  CALL force_env_create ( force_env, kg_env=kg_env, &
       globenv= globenv )

! set the integrator environment
  CALL init_md_env ( md_env , globenv % para_env )

  IF (kgpol) THEN
    CALL set_md_env ( md_env, simpar=simpar, struc=struc,           &
                              virial=kg_env % thermo % virial,      &
                              nhc_part=extended_type % nhc_part,    &
                              nhc_coef=extended_type % nhc_coef,    &
                              nhc_baro=extended_type % nhc_baro,    &
                              npt=extended_type % npt_info,         &
                              cell = kg_env % box,                  &
                              force_env=force_env)
  ELSE
    CALL set_md_env ( md_env, simpar=simpar, struc=struc,           &
                              virial=kg_env % thermo % virial,      &
                              nhc_part=extended_type % nhc_part,    &
                              nhc_baro=extended_type % nhc_baro,    &
                              npt=extended_type % npt_info,         &
                              cell = kg_env % box,                  &
                              force_env=force_env)
  END IF
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
     CALL control_debug_work ( globenv, kg_env, mdio )

  CASE ( "WAVEFUNCTION OPTIMIZATION" )

     CALL control_energy_work ( globenv, kg_env, md_env, energy_min, mdio )

  CASE ( "MOLECULAR DYNAMICS" )

     CALL control_md_work ( globenv, kg_env, md_env, mdio )

  END SELECT

  CALL kg_force_release (kg_env, globenv)

  DEALLOCATE ( struc, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'deallocstruc', 'struc' )
  CALL destroy_md_env(md_env)
  CALL timestop ( zero, handle )

END SUBROUTINE control_work

!******************************************************************************

 SUBROUTINE control_debug_work ( globenv, kg_env, mdio )

  IMPLICIT NONE

! Arguments
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
  TYPE ( kg_environment_type ), POINTER :: kg_env
  TYPE ( mdio_parameters_type ) :: mdio

!------------------------------------------------------------------------------

! initialize integrator
  CALL set_integrator ( globenv )

! debug routine
  CALL kg_debug_control ( globenv, kg_env )

 END SUBROUTINE control_debug_work

!******************************************************************************

 SUBROUTINE control_energy_work ( globenv, kg_env, md_env, &
                                  energy_min, mdio )

  IMPLICIT NONE

! Arguments
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
  TYPE ( kg_environment_type ), POINTER :: kg_env
  TYPE ( md_environment_type ), INTENT ( INOUT ) :: md_env
  TYPE ( energy_min_type ), INTENT ( INOUT ) :: energy_min
  TYPE ( mdio_parameters_type ) :: mdio

! Locals
  INTEGER :: iw, itimes
  LOGICAL, PARAMETER :: box_change = .FALSE.
  LOGICAL :: conv_flag


!------------------------------------------------------------------------------

! initialize integrator
  CALL set_integrator ( globenv )

  iw = globenv % scr
  itimes = 0

  CALL force_env_calculate_force ( md_env%force_env, box_change=box_change )

  energy_min % min_info %de = 0.0_dbl
  energy_min % min_info % e = kg_env % thermo%pot
  energy_min % min_info % e_hartree = kg_env % thermo % e_hartree
  energy_min % min_info % e_xc = kg_env %  thermo % e_xc
  energy_min % min_info % e_ppl = kg_env %  thermo % e_ppl
  energy_min % min_info % e_ppl2 = kg_env %  thermo % e_ppl2


  CALL print_energy_info ( itimes, energy_min % min_info,  &
                           energy_min % mintype, globenv )

  DO itimes = 1, energy_min % iter

    kg_env% thermo%pot = 0._dbl
    CALL optimize_ao (  kg_env, globenv, energy_min, itimes )

    CALL print_energy_info ( itimes, energy_min % min_info, energy_min % mintype, globenv )

!    write  coefficents to restart file
    IF ( MOD ( itimes, mdio % idump ) == 0 ) THEN
       CALL dump_coefs (kg_env % coef_pos (1) , kg_env % coef_vel (1), &
                          kg_env % coef_force (1))
    END IF

!   convergence acheived
    IF (energy_min % min_info % tol <= energy_min % max_tol) THEN
      conv_flag = .TRUE.
      EXIT
    END IF

  END DO

  CALL dump_coefs (kg_env % coef_pos (1) , kg_env % coef_vel (1), &
                           kg_env % coef_force (1))
  CALL print_convergence ( conv_flag, globenv )

 END SUBROUTINE control_energy_work


!******************************************************************************

 SUBROUTINE control_md_work ( globenv, kg_env, md_env, mdio )

  IMPLICIT NONE

! Arguments
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
  TYPE ( kg_environment_type ), POINTER :: kg_env
  TYPE ( md_environment_type ), INTENT ( INOUT ) :: md_env
  TYPE ( mdio_parameters_type ) :: mdio

! Locals
  INTEGER, POINTER :: itimes
  LOGICAL, PARAMETER :: box_change = .FALSE.

!------------------------------------------------------------------------------

! initialize the environment for md_energies
  CALL set_kg_energies ( globenv, mdio )

! initialize integrator
  CALL set_integrator ( globenv )

! MD
  itimes => md_env % itimes
  itimes = 0

  CALL force_env_calculate_force ( md_env%force_env, box_change=box_change )

! Initial call to md energy
  CALL md_energy_kg ( md_env , kg_env % thermo ) 

! Main MD loop
  DO itimes = 1, md_env % simpar % nsteps

! Call the integrator
    CALL velocity_verlet ( md_env )

! Call the energy routines
    CALL md_energy_kg ( md_env , kg_env % thermo ) 

! Test for the dump cycle
    IF ( MOD ( itimes, mdio % idump ) == 0 ) THEN
       CALL dump_variables ( md_env, kg_env, mdio % dump_file_name, globenv )
    END IF

  END DO
  CALL dump_variables ( md_env, kg_env, mdio % dump_file_name, globenv )
  IF ( globenv % ionode ) CALL close_unit ( 10, 99 )
 
END SUBROUTINE control_md_work

!******************************************************************************

SUBROUTINE allocmem ( mol_setup, part, pnode, molecule, globenv )

  IMPLICIT NONE

! Arguments
  TYPE ( molecule_type ), DIMENSION ( : ), INTENT ( IN ) :: mol_setup
  TYPE (molecule_structure_type), DIMENSION (:), POINTER :: molecule
  TYPE (particle_type), DIMENSION(:), INTENT (IN) :: part
  TYPE (particle_node_type), DIMENSION (:), POINTER :: pnode
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv

! Locals
  INTEGER :: iw, natoms, nnodes, nmol, nmoltype, ios, iat, i, nsh

!------------------------------------------------------------------------------

  IF ( globenv % num_pe == 1 ) THEN

     natoms = SIZE ( part  )

     ALLOCATE ( pnode ( natoms ), STAT = ios )
     IF ( ios /= 0 ) CALL stop_memory ( 'kg_main', 'pnode', natoms )

     nmol = SUM ( mol_setup ( : ) % num_mol )
     ALLOCATE ( molecule ( nmol ), STAT = ios )
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

     natoms = SIZE ( part  )
     nmoltype = SIZE ( mol_setup )
     nmol = 0
     nnodes = 0
     DO i = 1, nmoltype
        nsh = get_share ( mol_setup ( i ) % num_mol, &
             globenv % num_pe, globenv % mepos )
        nmol = nmol + nsh
        nnodes = nnodes + nsh * mol_setup ( i ) % molpar % natom
     END DO

     ALLOCATE ( molecule ( nmol ), STAT = ios )
     IF ( ios /= 0 ) CALL stop_memory ( 'kg_main', 'molecule' , nmol )
     ALLOCATE ( pnode ( nnodes ), STAT = ios )
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

SUBROUTINE deallocmem ( pnode, molecule )
  IMPLICIT NONE

! Arguments
  TYPE (molecule_structure_type),DIMENSION(:), POINTER :: molecule
  TYPE (particle_node_type), DIMENSION(:), POINTER :: pnode
! Locals
  INTEGER :: ios

!------------------------------------------------------------------------------

  DEALLOCATE ( pnode, STAT = ios )
  IF ( ios /= 0 ) CALL stop_memory ( 'kg_main', 'pnode' )

  DEALLOCATE ( molecule, STAT = ios )
  IF ( ios /= 0 ) CALL stop_memory ( 'kg_main', 'molecule' )

END SUBROUTINE deallocmem

!******************************************************************************

END MODULE kg_main
!******************************************************************************

!******************************************************************************
