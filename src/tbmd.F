!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

MODULE tbmd
  
  USE kinds, ONLY : dbl
  USE tbmd_global, ONLY : tbmdpar
  USE tbmd_input, ONLY : read_tbmd_section, read_tb_hamiltonian, &
                read_tb_hopping_elements
  USE tbmd_debug, ONLY : debug_control
  USE mathconstants, ONLY : zero
  USE timings, ONLY : timeset, timestop
  USE header, ONLY : tbmd_header
  USE ewalds, ONLY : ewald_print, ewald_correction
  USE ewald_parameters_types, ONLY : ewald_parameters_type
  USE input_types, ONLY : setup_parameters_type
  USE molecule_input, ONLY : read_molecule_section, read_setup_section
  USE molecule_types, ONLY : molecule_type, intra_parameters_type
  USE particle_types, ONLY : particle_prop_type, particle_type
  USE pair_potential, ONLY : get_potrange, potentialparm_type, &
       potential_f, spline_nonbond_control
  USE force_fields, ONLY : read_force_field_section
  USE atoms_input, ONLY : read_coord_vel, system_type
  USE md, ONLY : read_md_section, simulation_parameters_type, &
       initialize_velocities, thermodynamic_type, mdio_parameters_type
  USE unit, ONLY : unit_convert_type, set_units
  USE structure_types, ONLY : structure_type
  USE convert_units, ONLY : convert
  USE simulation_cell, ONLY : cell_type, get_hinv
  USE stop_program, ONLY : stop_memory
  USE integrator, ONLY : velocity_verlet, force, set_energy_parm, energy
  USE initialize_particle_types, ONLY : initialize_particle_type
  USE initialize_extended_types, ONLY : initialize_extended_type
  USE initialize_molecule_types, ONLY : initialize_molecule_type
  USE nose, ONLY : extended_parameters_type
  USE util, ONLY : close_unit, get_share
  USE linklist_control, ONLY : set_ll_parm
  USE dump, ONLY : dump_variables
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: tbmd_main

  TYPE ( mdio_parameters_type ) :: mdio
  
CONTAINS

!-----------------------------------------------------------------------------!
! TBMD TBMD TBMD TBMD TBMD TBMD TBMD TBMD TBMD TBMD TBMD TBMD TBMD TBMD TBMD  !
!-----------------------------------------------------------------------------!

SUBROUTINE tbmd_main
  
  IMPLICIT NONE
  
! Locals
  REAL (dbl) :: cons
  REAL (dbl), DIMENSION(:,:), ALLOCATABLE :: rcut
  INTEGER :: handle, itimes, isos
  INTEGER :: iat, jat
  CHARACTER ( LEN = 20 ), POINTER :: atom_names ( : )
  CHARACTER ( LEN = 20 ) :: set_fn
! 
  TYPE (potentialparm_type), POINTER :: potparm ( :, : )
  TYPE (particle_prop_type), POINTER :: pstat ( : )
  TYPE (molecule_type), POINTER :: mol_setup ( : )
  TYPE (unit_convert_type) :: units
  TYPE (simulation_parameters_type) :: simpar
  TYPE (structure_type) :: struc
  TYPE (extended_parameters_type) :: nhcp
  TYPE (thermodynamic_type) :: thermo
  TYPE (system_type) :: ainp
  TYPE (ewald_parameters_type) :: ewald_param
  TYPE (setup_parameters_type) :: setup
  TYPE (intra_parameters_type) :: intra_param      
  
!------------------------------------------------------------------------------
! 
  CALL timeset ( 'TBMD', 'I', ' ', handle )
! 
  IF ( tbmdpar % ionode ) CALL tbmd_header ( tbmdpar % scr )
! 
!..read tbmd section
! 
  CALL read_tbmd_section ( setup, ewald_param, tbmdpar )
! 
! read from the setup and molecule section of *.set
! 
  set_fn = setup % set_file_name
  
  CALL read_setup_section(mol_setup,set_fn,tbmdpar)
  
  CALL read_molecule_section(mol_setup,set_fn,tbmdpar)
! 
! read pair potential information 
! 
  CALL read_force_field_section ( setup, mol_setup, set_fn, &
       intra_param, potparm, atom_names, pstat, tbmdpar )
!
!..read Hamiltonian section
!
  CALL read_tb_hamiltonian ( setup, atom_names, tbmdpar )
  CALL read_tb_hopping_elements ( setup, atom_names, tbmdpar )
! 
!..read the input of the molecular dynamics section
! 
  CALL read_md_section(simpar,tbmdpar,mdio)
  simpar % program = 'TBMD'
! 
!..read atomic coordinates, velocities (optional) and the simulation box
  ainp % rtype = simpar % read_type
  CALL read_coord_vel ( ainp, setup % input_file_name, tbmdpar )
! 
!..initialize box, perd
  struc % box % hmat = ainp % box
  struc % box % perd = setup % perd
  
  CALL get_hinv ( struc % box )
! 
!..initialize working units
  CALL set_units ( setup % unit_type, units )
  
  CALL set_energy_parm ( units % pconv, units % econv, units % l_label, &
       units % vol_label, units % e_label, units % pres_label, &
       units % temp_label, units % angl_label )
! 
!..allocate memory for atoms and molecules
  CALL allocmem ( ainp, mol_setup, struc )
!
!..initialize particle_type
  CALL initialize_particle_type ( atom_names, simpar, mol_setup, &
       ainp, pstat, struc % part )
! 
!..convert the units
  CALL convert ( units = units, simpar = simpar, &
                 part = struc % part, pstat = pstat, box = struc % box, &
                 potparm = potparm, intra_param = intra_param, &
                 ewald_param = ewald_param )
!       
!..initialize molecule_type
  CALL initialize_molecule_type ( mol_setup, intra_param, struc % pnode, &
       struc % part, struc % molecule, tbmdpar )
! 
!..initialize extended_parameters_type and get number of degrees of freedom
  CALL initialize_extended_type ( struc % box, simpar, &
       struc % molecule, mol_setup, nhcp, tbmdpar )
! 
! initialize velocities if needed
  IF ( simpar % read_type == 'POS' ) THEN
     CALL initialize_velocities ( simpar, struc % part, tbmdpar )
  END IF
! 
!...initialize splines
! 
  potparm ( :, : ) % energy_cutoff = 0.0_dbl
  potparm ( :, : ) % e_cutoff_coul = 0.0_dbl
  CALL spline_nonbond_control ( potparm, pstat, 5000, ewald_param )
! 
!..set linklist control parameters
  ALLOCATE ( rcut ( setup % natom_type, setup % natom_type ), STAT = isos )
  IF ( isos /=0 ) CALL stop_memory ( 'tbmd', 'rcut' , 0)
  
  rcut ( :, : ) = potparm ( :, : ) % rcutsq
  
  CALL set_ll_parm ( tbmdpar, simpar % verlet_skin, &
       setup % natom_type, rcut, simpar % n_cell )

  CALL set_ll_parm ( tbmdpar, printlevel = tbmdpar % print_level, &
       ltype = 'NONBOND' )
  
  DEALLOCATE ( rcut, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'tbmd', 'rcut' )
! 
!..deallocate arrays needed for atom input
  IF ( ASSOCIATED ( ainp % c ) ) THEN
     DEALLOCATE ( ainp % c, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'tbmd', 'ainp%c' )
  END IF
  
  IF ( ASSOCIATED ( ainp % v ) ) THEN
     DEALLOCATE (ainp % v, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'tbmd', 'ainp%v' )
  END IF
! 
  IF ( setup % run_type == 'DEBUG' ) THEN

!..debug the forces
    CALL debug_control( tbmdpar, ewald_param, struc%part, struc%pnode, &
               struc%molecule, struc%box, thermo, potparm, simpar % ensemble )
     
  ELSE
!    
!..MD
     itimes = 0
     CALL force ( struc, potparm, thermo, simpar, ewald_param )
     IF ( tbmdpar % ionode .AND. ewald_param % ewald_type /= 'NONE' ) &
          CALL ewald_print( tbmdpar % scr, thermo, struc % box, &
          units % e_label )
     CALL energy ( itimes, cons, simpar, struc, thermo, nhcp )

     DO itimes = 1, simpar%nsteps
        CALL velocity_verlet ( itimes, cons, simpar, potparm,thermo, &
             struc, ewald_param, nhcp )
        IF ( MOD ( itimes, mdio % idump ) == 0 ) &
             CALL dump_variables ( struc, mdio % dump_file_name, tbmdpar )
     END DO

     CALL dump_variables ( struc, mdio % dump_file_name, tbmdpar )
     IF (tbmdpar % ionode) CALL close_unit(10,99)

  END IF
  
!..deallocate memory for atoms and molecules
  CALL deallocmem ( struc )
  
  CALL timestop ( zero, handle )
  
END SUBROUTINE tbmd_main

!-----------------------------------------------------------------------------!
! TBMD TBMD TBMD TBMD TBMD TBMD TBMD TBMD TBMD TBMD TBMD TBMD TBMD TBMD TBMD  !
!-----------------------------------------------------------------------------!

 SUBROUTINE allocmem ( ainp, mol_setup, struc )
   IMPLICIT NONE
  
! Arguments
  TYPE(molecule_type), DIMENSION(:), INTENT(IN) :: mol_setup
  TYPE(structure_type), INTENT(INOUT) :: struc
  TYPE(system_type), INTENT(IN) :: ainp
  
! Locals
  INTEGER :: natoms, nnodes, nmol, nmoltype, ios, iat, i, scr
  
!------------------------------------------------------------------------------
  
  struc % name = 'TBMD MOLECULAR SYSTEM'
  IF ( tbmdpar % num_pe == 1 ) THEN
     
     natoms = SIZE ( ainp % c ( 1, : ) )
     ALLOCATE ( struc % part ( 1 : natoms ), STAT = ios )
     IF ( ios /= 0 ) CALL stop_memory ( 'tbmd', 'part', natoms )
     ALLOCATE ( struc % pnode ( 1 :natoms ), STAT = ios )
     IF ( ios /= 0 ) CALL stop_memory ( 'tbmd', 'pnode' , natoms )
     nmol = SUM ( mol_setup ( : ) % num_mol )
     
     ALLOCATE ( struc % molecule ( 1 : nmol ), STAT = ios )
     IF ( ios /= 0 ) CALL stop_memory ( 'tbmd', 'molecule' , nmol )
     
     IF ( tbmdpar % ionode .AND. tbmdpar % print_level > 3 ) THEN
        scr = tbmdpar % scr
        WRITE ( scr, '( A )' )
        WRITE (scr,'(A,T71,I10)') ' TBMD| Number of allocated particles ', &
             natoms
        WRITE (scr,'(A,T71,I10)') &
             ' TBMD| Number of allocated particle nodes ', natoms
        WRITE ( scr, '( A, T71, I10 )' ) &
             ' TBMD| Number of allocated molecules ', nmol
        WRITE ( scr, '( A )' )
     END IF
  ELSE
!    
!..replicated data
     
     natoms = size(ainp%c(1,:))
     ALLOCATE (struc%part(1:natoms),STAT=ios)
     IF ( ios /= 0 ) CALL stop_memory ( 'tbmd', 'part' , natoms )
     nmoltype = SIZE ( mol_setup )
     nmol = 0
     nnodes = 0
     DO i = 1, nmoltype
        nmol = nmol + get_share ( mol_setup ( i ) % num_mol, &
             tbmdpar % num_pe, tbmdpar % mepos )
        nnodes = nnodes + nmol*mol_setup(i) %molpar%natom
     END DO
     
     ALLOCATE (struc%molecule(1:nmol),STAT=ios)
     IF (ios/=0) CALL stop_memory('tbmd','molecule',nmol)
     ALLOCATE (struc%pnode(1:nnodes),STAT=ios)
     IF (ios/=0) CALL stop_memory('tbmd','pnode',nnodes)
     
     IF ( tbmdpar%ionode .AND. tbmdpar%print_level>3) THEN
        WRITE ( scr, '( A )' )
        WRITE (scr,'(A,T71,I10)') ' TBMD| Number of allocated particles ' &
             , natoms
        WRITE (scr,'(A,T71,I10)') &
             ' TBMD| Number of allocated particle nodes ', nnodes
        WRITE (scr,'(A,I5,T71,I10)') &
             ' TBMD| Number of allocated molecules on processor ', &
             tbmdpar % mepos, nmol
        WRITE ( scr, '( A )' )
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
  
  DEALLOCATE ( struc % part, STAT = ios )
  IF ( ios /= 0 ) CALL stop_memory ( 'tbmd','part' )
  DEALLOCATE ( struc % pnode, STAT = ios )
  IF ( ios /= 0 ) CALL stop_memory ( 'tbmd','pnode' )
  DEALLOCATE ( struc % molecule, STAT = ios )
  IF ( ios /= 0 ) CALL stop_memory ( 'tbmd','molecule' )
  
 END SUBROUTINE deallocmem

!******************************************************************************

END MODULE tbmd
