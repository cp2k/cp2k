!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/fit_module [1.0] *
!!
!!   NAME
!!     fit_module
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     gloria, 29-MAY-2002
!!
!!   MODIFICATION HISTORY
!!   SOURCE
!******************************************************************************

MODULE fit_module

  USE atoms_input,                     ONLY: read_coord_vel,&
                                             system_type
  USE cntl_input,                      ONLY: read_cntl_section
  USE convert_units,                   ONLY: convert
!  USE dft_input, ONLY : read_dft_section, read_wave_section

  USE empirical_parameters,            ONLY: read_empirical_parameters
  USE energy_input,                    ONLY: energy_min_type,&
                                             read_energy_section
  USE energy_optimize,                 ONLY: optimize_ao,&
                                             optimize_pw
  USE energy_output,                   ONLY: print_convergence,&
                                             print_energy_info
  USE ewald_parameters_types,          ONLY: ewald_parameters_type
  USE extended_system_types,           ONLY: extended_system_type
  USE fit_force_control,               ONLY: fit_force
  USE fit_input,                       ONLY: read_fit_section
  USE fit_pair_potential,              ONLY: fit_potparm
  USE fit_parameters_types,            ONLY: fit_info_type,&
                                             fit_parameter_type
  USE fit_utilities,                   ONLY: force_ref
  USE force_fields,                    ONLY: read_force_field_section
  USE global_types,                    ONLY: global_environment_type
  USE initialize_molecule_types,       ONLY: initialize_molecule_type
  USE initialize_particle_types,       ONLY: initialize_particle_type
  USE initialize_pimd_types,           ONLY: initialize_pimd
  USE input_types,                     ONLY: setup_parameters_type
  USE integrator,                      ONLY: energy,&
                                             set_energy_parm
  USE kinds,                           ONLY: dp
  USE linklist_control,                ONLY: set_ll_parm
  USE mathconstants,                   ONLY: twopi,&
                                             zero
  USE md,                              ONLY: mdio_parameters_type,&
                                             read_md_section,&
                                             simulation_parameters_type,&
                                             thermodynamic_type
  USE method_specifications,           ONLY: get_maxder
  USE molecule_input,                  ONLY: read_molecule_section,&
                                             read_setup_section,&
                                             read_topology_section
  USE molecule_types,                  ONLY: intra_parameters_type,&
                                             molecule_type,&
                                             topology_type
  USE pair_potential,                  ONLY: spline_nonbond_control
  USE particle_types,                  ONLY: particle_prop_type
  USE physcon,                         ONLY: boltzmann,&
                                             h_planck
  USE pol_coefs,                       ONLY: allocate_pol_coefs,&
                                             dump_pol_coefs,&
                                             initialize_pol_coefs
  USE pol_setup,                       ONLY: get_rcutsq_cgf,&
                                             read_pol_basis_sets
  USE pw_grid_types,                   ONLY: HALFSPACE
  USE pw_grids,                        ONLY: pw_grid_setup
  USE read_pimd,                       ONLY: read_pimd_section
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_hinv,&
                                             init_cell
  USE string_utilities,                ONLY: integer_to_string,&
                                             xstring
  USE structure_types,                 ONLY: interaction_type,&
                                             structure_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop,&
                                             trace_debug
  USE transformations,                 ONLY: v2ud
  USE unit,                            ONLY: set_units,&
                                             unit_convert_type
  USE util,                            ONLY: get_share

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: fit
  
!!*****
!-----------------------------------------------------------------------------!

CONTAINS

!-----------------------------------------------------------------------------!
! CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL CNTL  !
!-----------------------------------------------------------------------------!
!!****** fit_module/fit [1.0] *
!!
!!   NAME
!!     fit
!!
!!   FUNCTION
!!     Controls program flow for POL              
!!
!!
!!   USED BY
!!     cp2k
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE fit ( globenv )


! Argument
    TYPE(global_environment_type), &
      INTENT(INOUT)                          :: globenv

    CHARACTER(LEN=40)                        :: project_name, set_fn
    CHARACTER(LEN=5)                         :: tag
    INTEGER                                  :: beads, handle1, handle2, i, &
                                                ia, ib, ibead, isos, maxder, &
                                                n, nbasis_types, ncoef, nll, &
                                                nproj
    LOGICAL                                  :: pimd, pol
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)  :: rcut, rcut_cgf
    TYPE(cell_type), DIMENSION(:), POINTER   :: box_ref
    TYPE(energy_min_type)                    :: energy_min
    TYPE(ewald_parameters_type)              :: ewald_param
    TYPE(extended_system_type)               :: extended_type
    TYPE(fit_info_type), DIMENSION(:, :), &
      POINTER                                :: finfo
    TYPE(fit_parameter_type)                 :: fit_param
    TYPE(interaction_type)                   :: inter
    TYPE(intra_parameters_type)              :: intra_param
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
    TYPE(topology_type)                      :: topo
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
    CALL stop_program ( ' fit ','PIMD not implemented with program' )
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

  CALL read_fit_section ( fit_param, globenv )

  IF (fit_param%fit_type == 'FORCES_POL') THEN 
    pol = .TRUE.
    nll = 2
  ELSE
    pol = .FALSE.
    nll = 1
  END IF

  DO ibead = 1, beads
    ALLOCATE ( struc ( ibead ) % ll_data ( nll ), STAT = isos )
    IF ( isos /=0 ) CALL stop_memory ( 'fit', 'll_data', 0 )
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

! read from the setup and molecule section of *.set

  set_fn = setup % set_file_name

  CALL read_setup_section ( mol_setup, setup, globenv )

  CALL read_molecule_section ( mol_setup, setup, globenv )

  CALL read_topology_section ( topo, setup, globenv )

! read force_field information for classical MD
! read pair potential information for TB

  CALL read_force_field_section ( setup, mol_setup, set_fn, &
       intra_param, inter%potparm, pstat, globenv )

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
  simpar % pimd_params % wp = 1.0e-15_dp * SQRT ( 1.0_dp * beads ) * twopi * &
        simpar % temp_ext * boltzmann / h_planck
  simpar % pimd_params % beta = 1.0_dp / simpar % temp_ext
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
    CALL init_cell ( struc ( ibead ) % box, ainp % box, &
       setup % simulation_cell % perd )
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

  END DO

    CALL release_mem ( ainp )

  CALL set_energy_parm ( units % pconv, units % econv, units % l_label, &
       units % vol_label, units % e_label, units % pres_label, &
       units % temp_label, units % angl_label )

  IF ( pimd ) THEN
!..initialize pimd_data (transformation related variables)
    CALL initialize_pimd ( struc , simpar )

    CALL v2ud ( struc, simpar )
  END IF

  
  IF (pol) THEN
    maxder = 1    
    CALL get_maxder(maxder)

! reads basis set and set up coefficents arrays 

     NULLIFY ( struc ( 1 ) % coef_pos )
     NULLIFY ( struc ( 1 ) % coef_vel )
     NULLIFY ( struc ( 1 ) % coef_force )

     ALLOCATE ( struc ( 1 ) % coef_pos ( 1 ) , STAT = isos )
     IF ( isos /=0 ) CALL stop_memory ( 'fit', 'coef_pos', 1 )

     ALLOCATE ( struc ( 1 ) % coef_vel ( 1 ) , STAT = isos )
     IF ( isos /=0 ) CALL stop_memory ( 'fit', 'coef_vel', 1 )

     ALLOCATE ( struc ( 1 ) % coef_force ( 1 ) , STAT = isos )
     IF ( isos /=0 ) CALL stop_memory ( 'fit', 'coef_force', 1 )

     IF ( simpar % pol_basis == "PW" ) THEN
       ALLOCATE ( struc ( 1 ) % projector ( 1 ) )
       IF ( isos /=0 ) CALL stop_memory ( 'control', 'projector', 1 )

       ALLOCATE ( struc ( 1 ) % pw_grid ( 1 ) )
       IF ( isos /=0 ) CALL stop_memory ( 'control', 'pw_grid', 1 )
     ENDIF

     CALL read_pol_basis_sets( struc ( 1 ) % rho0_basis_info,  &
          struc ( 1 ) % drho_basis_info, mol_setup, struc ( 1 ) % part, &
          setup % atom_names, n, globenv )

     nbasis_types  = size ( struc ( 1 ) % drho_basis_info )

     CALL read_empirical_parameters ( struc (1) % drho_basis_info, set_fn, &
                setup % atom_names, inter%empparm, globenv )

     IF ( simpar % pol_basis == "AO" ) THEN
       ncoef = n
       CALL allocate_pol_coefs ( ncoef, struc ( 1 ) % coef_pos ( 1 ),  &
          struc ( 1 ) % coef_vel ( 1 ), struc ( 1 ) % coef_force ( 1 ) )

       CALL initialize_pol_coefs ( simpar, struc ( 1 ) % drho_basis_info, &
          struc ( 1 ) % part, struc ( 1 ) % coef_pos ( 1 ),  &
          struc ( 1 ) % coef_vel ( 1 ), struc ( 1 ) % coef_force ( 1 ), globenv )

     ELSEIF ( simpar % pol_basis == "PW" ) THEN

       struc ( 1 ) % pw_grid ( 1 ) % grid_span = HALFSPACE
       struc ( 1 ) % pw_grid ( 1 ) % bounds = 0
       CALL pw_grid_setup ( struc ( 1 ) % box, struc ( 1 ) % pw_grid ( 1 ),  &
                            simpar % pol_ecut )
       nproj = n
       CALL allocate_pol_coefs ( nproj, struc ( 1 ) % coef_pos ( 1 ),  &
          struc ( 1 ) % coef_vel ( 1 ), struc ( 1 ) % coef_force ( 1 ), &
          struc ( 1 ) % projector ( 1 ),  struc ( 1 ) % pw_grid ( 1 ) )

       CALL initialize_pol_coefs ( simpar, struc ( 1 ) % drho_basis_info, &
          struc ( 1 ) % part, struc ( 1 ) % coef_pos ( 1 ),  &
          struc ( 1 ) % coef_vel ( 1 ), struc ( 1 ) % coef_force ( 1 ), globenv, &
          struc ( 1 ) % projector ( 1 ) )
     ENDIF

  ENDIF

!...initialize splines

  inter%potparm ( :, : ) % energy_cutoff = 0.0_dp
  inter%potparm ( :, : ) % e_cutoff_coul = 0.0_dp
  CALL spline_nonbond_control ( inter%potparm, pstat, 2500, ewald_param )

!..set linklist control parameters
  ALLOCATE ( rcut ( setup % natom_type, setup % natom_type ), STAT = isos )
  IF ( isos /=0 ) CALL stop_memory ( 'fit', 'rcut', 0 )

  rcut ( :, : ) = inter%potparm ( :, : ) % rcutsq
     
  

  DO ibead = 1, beads 
! construct the nonbonded  lists for the particles  
    CALL set_ll_parm ( struc ( ibead ) % ll_data ( 1 ), globenv, simpar % verlet_skin, &
       setup % natom_type, rcut, simpar % n_cell )

    CALL set_ll_parm ( struc ( ibead ) % ll_data ( 1 ), globenv,  &
       printlevel = globenv % print_level, ltype = 'NONBOND' )

  END DO
! construct the lists for the overlaps
  IF (POL) then
    ALLOCATE (rcut_cgf ( nbasis_types, nbasis_types ), STAT = isos)
    IF ( isos /=0 ) CALL stop_memory ( 'fit', 'rcut_cgf', 0 )

    CALL get_rcutsq_cgf ( struc ( 1 ) % drho_basis_info, rcut_cgf )
    CALL set_ll_parm ( struc ( 1 ) % ll_data ( 2 ), globenv, simpar % verlet_skin, &
         nbasis_types, rcut_cgf, simpar % n_cell )

    CALL set_ll_parm ( struc ( 1 ) % ll_data ( 2 ), globenv,  &
          printlevel = globenv % print_level, ltype = 'OVERLAP' )
  END IF

  DEALLOCATE ( rcut, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'fit', 'rcut' )

  IF (POL) then
     DEALLOCATE ( rcut_cgf, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'fit', 'rcut_cgf' )
  END IF

  CALL timestop ( zero, handle2 )

! Get reference force

  CALL force_ref(struc(1), inter % potparm,simpar,fit_param,finfo,box_ref)

  DO i =1, size(box_ref)
     CALL init_cell ( box_ref ( i ), box_ref (i) % hmat, setup % simulation_cell % perd )
  END DO

  CALL fit_work ( globenv, units, simpar, struc, inter, &
                  ewald_param, setup, mdio, fit_param,  box_ref, finfo, pol, &
                  pstat, energy_min )

! deallocate memory for atoms and molecules
  DO ibead = 1, beads
    CALL deallocmem ( struc ( ibead ) )
  END DO
  DEALLOCATE ( struc, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'fit', 'struc' )
  CALL timestop ( zero, handle1 )
  
END SUBROUTINE fit

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
     IF ( isos /= 0 ) CALL stop_memory ( 'fit', 'ainp%c' )
  END IF

  IF ( ASSOCIATED ( ainp % v ) ) THEN
     DEALLOCATE ( ainp % v, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'fit', 'ainp%v' )
  END IF

  IF ( ASSOCIATED ( ainp % eta_part ) ) THEN
     DEALLOCATE ( ainp % eta_part, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'fit', 'ainp%eta_part' )
  END IF

  IF ( ASSOCIATED ( ainp % veta_part ) ) THEN
     DEALLOCATE ( ainp % veta_part, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'fit', 'ainp%veta_part' )
  END IF

  IF ( ASSOCIATED ( ainp % mass_eta_part ) ) THEN
     DEALLOCATE ( ainp % mass_eta_part, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'fit', 'ainp%mass_eta_part' )
  END IF

  IF ( ASSOCIATED ( ainp % eta_baro ) ) THEN
     DEALLOCATE ( ainp % eta_baro, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'fit', 'ainp%eta_baro' )
  END IF

  IF ( ASSOCIATED ( ainp % veta_baro ) ) THEN
     DEALLOCATE ( ainp % veta_baro, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'fit', 'ainp%veta_baro' )
  END IF

  IF ( ASSOCIATED ( ainp % mass_eta_baro ) ) THEN
     DEALLOCATE ( ainp % mass_eta_baro, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'fit', 'ainp%mass_eta_baro' )
  END IF

  IF ( ASSOCIATED ( ainp % veps ) ) THEN
     DEALLOCATE ( ainp % veps, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'fit', 'ainp%veps' )
  END IF

  IF ( ASSOCIATED ( ainp % mass_eps ) ) THEN
     DEALLOCATE ( ainp % mass_eps, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'fit', 'ainp%mass_eps' )
  END IF 
END SUBROUTINE release_mem
!******************************************************************************

SUBROUTINE fit_work ( globenv, units, simpar, struc, inter, &
     ewald_param,  setup, mdio,   &
     fit_param, box_ref, finfo, pol, pstat, energy_min )
  
  
    TYPE(global_environment_type), &
      INTENT(INOUT)                          :: globenv
    TYPE(unit_convert_type)                  :: units
    TYPE(simulation_parameters_type)         :: simpar
    TYPE(structure_type), DIMENSION(:)       :: struc
    TYPE(interaction_type)                   :: inter
    TYPE(ewald_parameters_type)              :: ewald_param
    TYPE(setup_parameters_type)              :: setup
    TYPE(mdio_parameters_type)               :: mdio
    TYPE(fit_parameter_type), INTENT(INOUT)  :: fit_param
    TYPE(cell_type), DIMENSION(:), POINTER   :: box_ref
    TYPE(fit_info_type), DIMENSION(:, :), &
      POINTER                                :: finfo
    LOGICAL, INTENT(in)                      :: POL
    TYPE(particle_prop_type), DIMENSION(:), &
      POINTER                                :: pstat
    TYPE(energy_min_type), INTENT(INOUT)     :: energy_min

    INTEGER                                  :: handle, ii, ipart, itime, jj, &
                                                npart
    LOGICAL                                  :: box_change
    REAL(KIND=dp)                                :: boxtol, deltabox

!------------------------------------------------------------------------------

  CALL timeset ( 'FIT_WORK', 'I', ' ', handle )
  box_change = .FALSE.
  deltabox = zero
  boxtol = 1.E-4_dp
  
  SELECT CASE ( fit_param%fit_type )
  CASE DEFAULT
     CALL stop_program ( "fit_work", "no suitable run_type" &
                         //  fit_param % fit_type )
     
  CASE ( "FORCES_POL" )

     npart = SIZE ( struc (1) % part )
     DO itime = 1, fit_param%nsteps  ! loop over the sampled forces
        box_change = .FALSE.
        DO ipart = 1, npart
           struc (1) %pnode(ipart) %p%r(:) = finfo(itime,ipart) %r(:)
        END DO

!check if the box changes           
outer: DO ii = 1,3
inner:    DO JJ = 1,3 
             deltabox = box_ref(itime)%hmat(ii,jj) - struc (1) % box % hmat (ii,jj)
             IF (ABS(deltabox) >= boxtol ) THEN 
                box_change = .TRUE.
                EXIT outer
             END IF
          END DO inner
       END DO outer
!pass the new box, if the box changes
       IF ( box_change ) struc (1) %box=box_ref(itime)

       CALL fit_energy_work ( globenv, units, simpar, struc (1), inter, &
            ewald_param, mdio, fit_param, box_change, energy_min )
       DO ipart = 1, npart
          finfo(itime,ipart) %f_ref(:) = finfo(itime,ipart) %f_ref(:) - &
                    struc (1) %pnode(ipart) %p%f(:)
          write(55,'(i7,6f18.10)') itime, &
                 finfo(itime,ipart)%r(:),finfo(itime,ipart)%f_ref(:)
       END DO
       DO ii=1,3
          WRITE (55,'(i7,3f12.6)') itime, (box_ref(itime)%hmat(ii,jj),jj=1,3)
       ENDDO
    END DO

  CASE ( "AMOEBA_RI" )
     
     CALL fit_potparm ( struc ( 1 ), inter,  simpar, fit_param, pstat, &
        ewald_param, box_ref, finfo, setup%atom_names, globenv)
  END SELECT
  
  CALL timestop ( zero, handle )
  
END SUBROUTINE fit_work

!******************************************************************************


SUBROUTINE fit_energy_work ( globenv, units, simpar, struc, inter, &
      ewald_param, mdio, fit_param, box_change, energy_min )
  
  
    TYPE(global_environment_type), &
      INTENT(INOUT)                          :: globenv
    TYPE(unit_convert_type)                  :: units
    TYPE(simulation_parameters_type)         :: simpar
    TYPE(structure_type)                     :: struc
    TYPE(interaction_type)                   :: inter
    TYPE(ewald_parameters_type)              :: ewald_param
    TYPE(mdio_parameters_type)               :: mdio
    TYPE(fit_parameter_type)                 :: fit_param
    LOGICAL, INTENT(in)                      :: box_change
    TYPE(energy_min_type), INTENT(INOUT)     :: energy_min

    INTEGER                                  :: ipart, itimes, iw, npart
    LOGICAL                                  :: coef_flag, conv_flag
    TYPE(thermodynamic_type)                 :: thermo

!------------------------------------------------------------------------------

  iw = globenv% scr 
  npart = SIZE ( struc %pnode)
  coef_flag = .FALSE.
  conv_flag = .FALSE.

  coef_flag = ( ASSOCIATED ( struc % coef_pos ) .AND.  &
                ASSOCIATED ( struc % coef_vel ) .AND.  &
                ASSOCIATED ( struc % coef_force ) ) 

  itimes = 0


  CALL fit_force ( struc, inter, thermo, simpar, ewald_param, box_change, &
       fit_param, globenv)

  energy_min % min_info % e = thermo%pot

  DO itimes = 1, energy_min % iter   

     thermo%pot = 0.0_dp
     IF ( simpar % pol_basis == "AO" ) THEN
        CALL optimize_ao ( struc, inter, thermo, simpar, ewald_param, box_change, &
                        globenv, fit_param, energy_min, itimes)
     ELSEIF ( simpar % pol_basis == "PW" ) THEN
        CALL optimize_pw ( struc, inter, thermo, simpar, ewald_param, box_change, &
                        globenv, fit_param, energy_min, itimes)
     END IF

     CALL print_energy_info ( itimes, energy_min % min_info, energy_min % mintype, globenv )

!    write  coefficents to restart file
     IF (( MOD ( itimes, mdio % idump ) == 0 ) .AND. (coef_flag)) THEN
        CALL dump_pol_coefs (struc % coef_pos (1) , struc % coef_vel (1), &
                           struc % coef_force (1))
     END IF 

!   convergence acheived
    IF (energy_min % min_info % tol <= energy_min % max_tol) THEN
       conv_flag = .TRUE.
    WRITE (iw,*)
    WRITE (iw,'(A,A)') ' **************************************', &
      '*****************************************'
    WRITE (iw,'(A)') ' FORCES'                                
    DO ipart = 1, npart
       write(iw,'(i7,6f18.10)')  ipart, &
                struc %pnode(ipart) %p%f(:)                    
    END DO
    EXIT
    END IF

  END DO

  CALL dump_pol_coefs (struc % coef_pos (1) , struc % coef_vel (1), &
                           struc % coef_force (1))
  CALL print_convergence ( conv_flag, globenv )
  
END SUBROUTINE fit_energy_work

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
     IF ( ios /= 0 ) CALL stop_memory ( 'fit', 'part', natoms )
     ALLOCATE ( struc % pnode ( natoms ), STAT = ios )
     IF ( ios /= 0 ) CALL stop_memory ( 'fit', 'pnode', natoms )
     nmol = SUM ( mol_setup ( : ) % num_mol )

     ALLOCATE ( struc % molecule ( nmol ), STAT = ios )
     IF ( ios /= 0 ) CALL stop_memory ( 'fit', 'molecule', nmol )

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
     IF ( ios /= 0 ) CALL stop_memory ( 'fit', 'part', natoms )
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
     IF ( ios /= 0 ) CALL stop_memory ( 'fit', 'molecule' , nmol )
     ALLOCATE ( struc % pnode ( nnodes ), STAT = ios )
     IF ( ios /= 0 ) CALL stop_memory ( 'fit', 'pnode', nnodes )

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
  IF ( ios /= 0 ) CALL stop_memory ( 'fit', 'part' )

  DEALLOCATE ( struc % pnode, STAT = ios )
  IF ( ios /= 0 ) CALL stop_memory ( 'fit', 'pnode' )

  DEALLOCATE ( struc % molecule, STAT = ios )
  IF ( ios /= 0 ) CALL stop_memory ( 'fit', 'molecule' )

END SUBROUTINE deallocmem

!******************************************************************************

END MODULE fit_module
!******************************************************************************
