!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

MODULE control_module
  
  USE atoms_input, ONLY : read_coord_vel, system_type
  USE convert_units, ONLY : convert
  USE global_types, ONLY : global_environment_type
  USE dump, ONLY : dump_variables
  USE ewalds, ONLY : ewald_print, ewald_correction
  USE ewald_parameters_types, ONLY : ewald_parameters_type
  USE header, ONLY : fist_header, tbmd_header
  USE kinds, ONLY : dbl
  USE band, ONLY : band_structure_type, init_band_structure
  USE brillouin, ONLY : kpoint_type, brillouin_info, kpoint_input
  USE kpoint_initialization, ONLY : initialize_kpoints
  USE fermi, ONLY : fermi_distribution_type, init_fermi_dist, fermi_info
  USE fist_debug, ONLY : fist_debug_control => debug_control
  USE tbmd_debug, ONLY : tbmd_debug_control => debug_control
  USE tbmd_initialize, ONLY : tbmd_init, tb_get_numel
  USE tbmd_global, ONLY : tbatom, tbhop
  USE tbmd_input, ONLY : read_tb_hamiltonian, read_tb_hopping_elements
  USE cntl_input, ONLY : read_cntl_section
  USE force_fields, ONLY : read_force_field_section, ATOMNAMESLENGTH
  USE initialize_extended_types, ONLY : initialize_extended_type
  USE initialize_molecule_types, ONLY : initialize_molecule_type
  USE initialize_particle_types, ONLY : initialize_particle_type
  USE integrator, ONLY : velocity_verlet, force, set_energy_parm, energy, &
       set_integrator
  USE input_types, ONLY : setup_parameters_type
  USE linklist_control, ONLY : set_ll_parm
  USE mathconstants, ONLY : zero
  USE md, ONLY : read_md_section, simulation_parameters_type, &
       initialize_velocities, thermodynamic_type, mdio_parameters_type
  USE molecule_input, ONLY : read_molecule_section, read_setup_section, &
       charge
  USE molecule_types, ONLY : molecule_type, intra_parameters_type
  USE nose, ONLY : extended_parameters_type
  USE pair_potential, ONLY : spline_nonbond_control
  USE particle_types, ONLY : particle_prop_type, particle_type
  USE simulation_cell, ONLY : cell_type, get_hinv
  USE stop_program, ONLY : stop_prg, stop_memory
  USE structure_types, ONLY : structure_type, interaction_type
  USE timings, ONLY : timeset, timestop, trace_debug
  USE unit, ONLY : unit_convert_type, set_units
  USE util, ONLY : close_unit, get_share
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: control
  
  TYPE ( mdio_parameters_type ) :: mdio
  
CONTAINS

!-----------------------------------------------------------------------------!
! CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL  !
!-----------------------------------------------------------------------------!

SUBROUTINE control ( globenv )
  
  IMPLICIT NONE

  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
  
! Locals
  REAL ( dbl ) :: cons, ecut, qi, qj
  REAL ( dbl ), DIMENSION ( :, : ), ALLOCATABLE :: rcut
  INTEGER :: handle1, handle2, handle3 , itimes, isos, para_comm_cart
  INTEGER :: iat, jat, nel
  CHARACTER ( LEN = ATOMNAMESLENGTH ), DIMENSION ( : ), POINTER :: atom_names
  CHARACTER ( LEN = 20 ) :: set_fn
  LOGICAL :: tbmd, fist
  
  TYPE ( particle_prop_type ), DIMENSION ( : ), POINTER :: pstat
  TYPE ( molecule_type ), DIMENSION ( : ), POINTER :: mol_setup
  TYPE ( unit_convert_type ) :: units
  TYPE ( simulation_parameters_type ) :: simpar
  TYPE ( structure_type ) :: struc
  TYPE ( interaction_type ) :: inter
  TYPE ( extended_parameters_type ) :: nhcp
  TYPE ( thermodynamic_type ) :: thermo
  TYPE ( system_type ) :: ainp
  TYPE ( ewald_parameters_type ) :: ewald_param
  TYPE ( setup_parameters_type ) :: setup
  TYPE ( intra_parameters_type ) :: intra_param
  TYPE ( kpoint_type ) :: kp
  TYPE ( fermi_distribution_type ) :: fd
  TYPE ( band_structure_type ) :: bs
  
!------------------------------------------------------------------------------

! IF( globenv % ionode ) CALL trace_debug ( "START" )
  
  CALL timeset ( 'CONTROL', 'I', ' ', handle1 )
  CALL timeset ( 'CNTL_INIT', 'I', ' ', handle2 )
  
  IF ( globenv % program_name == 'FIST' ) THEN
    fist = .true.
    tbmd = .false.
  ELSE IF ( globenv % program_name == 'TBMD' ) THEN
    tbmd = .true.
    fist = .false.
  ELSE
    call stop_prg ( ' control ',' program_name not specified' )
  ENDIF

  IF ( globenv % ionode ) THEN
    IF ( fist ) CALL fist_header ( globenv % scr )
    IF ( tbmd ) CALL tbmd_header ( globenv % scr )
  END IF
  
!..read cntl section
  
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

!..read the input of the molecular dynamics section
  
  CALL read_md_section ( simpar, globenv, mdio )
  simpar % program = globenv % program_name
  
!..read atomic coordinates, velocities (optional) and the simulation box
  ainp % rtype = simpar % read_type
  CALL read_coord_vel ( ainp, setup % input_file_name, globenv )
  
!..initialize box, perd
  struc % box % hmat = ainp % box
  struc % box % perd = setup % perd
  
!..K-Points
  IF ( SUM ( struc % box % perd ) == 0 ) THEN
    kp%scheme = "NULL"
  ELSE IF ( tbmd ) THEN
    CALL kpoint_input ( kp, globenv )
  ELSE
    kp%scheme = "NULL"
  END IF

!..initialize working units
  CALL set_units ( setup % unit_type, units )
  
  CALL set_energy_parm ( units % pconv, units % econv, units % l_label, &
       units % vol_label, units % e_label, units % pres_label, &
       units % temp_label, units % angl_label )
  
!..allocate memory for atoms and molecules
  CALL allocmem ( ainp, mol_setup, struc, globenv )
  
!..initialize particle_type
  CALL initialize_particle_type ( atom_names, simpar, mol_setup, &
       ainp, pstat, struc % part )
  
!..convert the units
  CALL convert ( units = units, simpar = simpar, &
                 part = struc % part, pstat = pstat, box = struc % box, &
                 potparm = inter%potparm, intra_param = intra_param, &
                 ewald_param = ewald_param )

!..calculate the inverse box matrix now after the unit conversion, 
!  so it also has the right units
  CALL get_hinv ( struc % box )
  
!..initialize molecule_type
  CALL initialize_molecule_type ( mol_setup, intra_param, struc % pnode, &
       struc % part, struc % molecule, globenv )
  
!..initialize extended_parameters_type and get number of degrees of freedom
  CALL initialize_extended_type ( struc % box, simpar, &
       struc % molecule, mol_setup, nhcp, globenv )
  
! initialize velocities if needed
  IF ( simpar % read_type == 'POS' ) THEN
     CALL initialize_velocities ( simpar, struc % part, globenv )
  END IF
  
!...initialize splines
  
  inter%potparm ( :, : ) % energy_cutoff = 0.0_dbl
  inter%potparm ( :, : ) % e_cutoff_coul = 0.0_dbl
  CALL spline_nonbond_control ( inter%potparm, pstat, 5000, ewald_param )
  
!..set linklist control parameters
  ALLOCATE ( rcut ( setup % natom_type, setup % natom_type ), STAT = isos )
  IF ( isos /=0 ) CALL stop_memory ( 'control', 'rcut', 0 )
  
  rcut ( :, : ) = inter%potparm ( :, : ) % rcutsq
  
  CALL set_ll_parm ( globenv, simpar % verlet_skin, &
       setup % natom_type, rcut, simpar % n_cell )
  
  CALL set_ll_parm ( globenv, printlevel = globenv % print_level, &
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
!
!..initialize the on-site terms for TB
  IF ( tbmd ) CALL tbmd_init ( struc % part, tbatom, tbhop )
!
!..Symmetry and K-points
!
  IF ( SUM ( struc % box % perd ) == 0 ) THEN
!!!!!! symmetry setup for molecules
  ELSE IF ( tbmd ) THEN
    CALL initialize_kpoints ( globenv, kp, setup % symmetry, &
         struc % box % hmat, struc % part )
    IF ( globenv%ionode .AND. globenv%print_level > 0) &
       CALL brillouin_info ( kp, globenv%scr )
  END IF
!..nota bene: we use the defaults for the elctron temperature and
!             spin polarisation
  IF ( tbmd ) then
    nel = tb_get_numel ( tbatom, charge )
    CALL init_fermi_dist(fd,nel)
    IF (globenv%ionode) CALL fermi_info(fd,globenv%print_level,globenv%scr)
    CALL init_band_structure(bs,fd,kp)
  ENDIF
  
  CALL timestop ( zero, handle2 )
  CALL timeset ( 'CNTL_WORK', 'I', ' ', handle3 )
  IF ( setup % run_type == 'DEBUG' ) THEN
!..debug the forces

!..initialize integrator
    CALL set_integrator ( globenv, mdio )

    IF ( fist ) then
      CALL fist_debug_control( globenv, ewald_param, struc%part, struc%pnode, &
                 struc%molecule, struc%box, thermo, inter%potparm, simpar % ensemble )
    ELSE IF ( tbmd ) THEN
      CALL tbmd_debug_control( globenv, ewald_param, struc%part, struc%pnode, &
                 struc%molecule, struc%box, thermo, inter%potparm, simpar % ensemble )
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
    CALL energy ( itimes, cons, simpar, struc, thermo, nhcp )
    
    DO itimes = 1, simpar % nsteps
       CALL velocity_verlet ( itimes, cons, simpar, inter, thermo, &
            struc, ewald_param, nhcp )
       
       IF ( MOD ( itimes, mdio % idump ) == 0 ) &
            CALL dump_variables ( struc, mdio % dump_file_name, globenv )
    END DO
    
    CALL dump_variables ( struc, mdio % dump_file_name, globenv )
    
    IF ( globenv % ionode ) CALL close_unit ( 10, 99 )
  END IF
  
!..deallocate memory for atoms and molecules
  CALL deallocmem ( struc )
  
  CALL timestop ( zero, handle3 )
  CALL timestop ( zero, handle1 )
  
END SUBROUTINE control

!-----------------------------------------------------------------------------!
! CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL  !
!-----------------------------------------------------------------------------!

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
