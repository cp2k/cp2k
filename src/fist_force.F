!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C ) 2000  CP2K developers group                                !
!-----------------------------------------------------------------------------!
!!****** cp2k/fist_force [1.0] 
!!
!!   NAME
!!     fist_force
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM & JGH
!!
!!   MODIFICATION HISTORY
!!     Torsions added (DG) 05-Dec-2000
!!     Variable names changed (DG) 05-Dec-2000
!!     CJM SEPT-12-2002: int_env is now passed
!!     CJM NOV-30-2003: only uses fist_env
!!
!!   SOURCE
!******************************************************************************

MODULE fist_force
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_assertion_failed
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_logger_type,&
                                             cp_warning_level
  USE cp_output_handling, ONLY: cp_should_output
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE ewalds, ONLY : ewald_evaluate, ewald_print
  USE ewald_environment_types, ONLY : ewald_environment_type, ewald_env_get
  USE ewald_pw_methods, ONLY        : ewald_pw_grid_change
  USE ewald_pw_types, ONLY : ewald_pw_type
  USE fist_energy_types, ONLY : fist_energy_type
  USE fist_environment_types, ONLY : fist_environment_type, get_fist_env
  USE fist_intra_force, ONLY : force_intra_control
  USE fist_neighbor_list_control, ONLY : list_control
  USE fist_nonbond_env_types, ONLY : fist_nonbond_env_type
  USE fist_nonbond_force, ONLY : force_nonbond, density_nonbond, &
                                 bonded_correct_gaussian, energy_manybody
  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY: dp, default_string_length
  USE mathconstants, ONLY : pi, zero
  USE message_passing, ONLY : mp_sum
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_types, ONLY : particle_type
  USE pme, ONLY: pme_evaluate
  USE simulation_cell, ONLY : cell_type, get_hinv
  USE spme, ONLY : spme_evaluate
  USE termination, ONLY : stop_program, stop_memory
  USE timings, ONLY : timeset, timestop
  USE virial_methods, ONLY : virial_evaluate
  USE virial_types, ONLY : virial_type
  USE input_section_types, ONLY: section_type, section_create, &
       section_add_keyword, section_add_subsection, section_release,&
       section_describe_html, section_describe_index_html, section_vals_type,&
       section_vals_get_subs_vals, section_vals_get, section_vals_val_get
  PRIVATE
  PUBLIC :: fist_force_control, debug_variables_type

  TYPE debug_variables_type
     REAL (KIND=dp) :: pot_nonbond, pot_g, pot_bond, pot_bend, pot_torsion
     REAL (KIND=dp) :: pot_imptors, pot_onef, pot_urey_bradley
     REAL (KIND=dp), DIMENSION ( :, : ), POINTER :: &
          f_nonbond, f_g, f_bond, f_bend, f_torsion, f_imptors, f_onef, f_ub
     REAL (KIND=dp), DIMENSION ( 3, 3 ) :: pv_nonbond, pv_g, pv_bond, &
          pv_bend, pv_torsion, pv_imptors, pv_onef, pv_ub
  END TYPE debug_variables_type

!******************************************************************************

 CONTAINS

!******************************************************************************
!!****** fist_force/fist_force_control [1.1] *
!!
!!   NAME
!!    fist_force_control 
!!
!!   FUNCTION
!!     Calculates the total potential energy, total force, and the
!!     total pressure tensor from the potentials
!!
!!   AUTHOR
!!     CJM & JGH
!!
!!   MODIFICATION HISTORY
!!     Harald Forbert (Dec-2000): Changes for multiple linked lists
!!     cjm, 20-Feb-2001: box_ref used to initialize ewald.  Now
!!                       have consistent restarts with npt and ewald
!!     JGH (15-Mar-2001) : box_change replaces ensemble parameter
!!                         Call ewald_setup if box has changed
!!                         Consistent setup for PME and SPME
!!
!!   SOURCE
!!******************************************************************************

SUBROUTINE fist_force_control ( fist_env, virial, globenv, box_change, debug ,error)

  IMPLICIT NONE

! Arguments
  TYPE ( fist_environment_type ), POINTER :: fist_env
  TYPE ( virial_type ), POINTER :: virial
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv
  LOGICAL, INTENT ( IN ) :: box_change
  TYPE ( debug_variables_type ), INTENT ( INOUT ), OPTIONAL :: debug
  TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                               :: error
! Locals
  INTEGER :: id, i, ii, natoms, nnodes, handle, isos, n, nlocal_particles
  INTEGER :: ikind, iparticle_kind, nparticle_kind, iparticle_local, nparticle_local
  INTEGER :: node, iw, j
  REAL (KIND=dp) :: pot_nonbond, pot_bond, pot_bend, pot_urey_bradley,  &
                  pot_imptors, pot_torsion, pot_onef,pot_manybody,vg_coulomb
  REAL (KIND=dp), DIMENSION ( 3,3 ) :: pv_nonbond, pv_bond, pv_bend,  &
                 pv_torsion, pv_urey_bradley, pv_imptors, pv_onef
  REAL (KIND=dp), DIMENSION ( :,: ), ALLOCATABLE :: fg_coulomb
  REAL (KIND=dp), DIMENSION ( :,: ), ALLOCATABLE :: f_total
  REAL (KIND=dp), DIMENSION ( 3,3 ) :: pv_g, ident
  REAL (KIND=dp), DIMENSION ( 3,3 ) :: pv_bc
  CHARACTER ( len = default_string_length ) :: ewald_type
  TYPE ( fist_energy_type ), POINTER :: thermo
  TYPE ( particle_type ), DIMENSION(:), POINTER :: particle_set 
  TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                              local_particles
  TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
  TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
  TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
  TYPE ( fist_nonbond_env_type ), POINTER    :: fist_nonbond_env
  TYPE ( cell_type ), POINTER                :: cell
  TYPE ( ewald_pw_type ),   POINTER          :: ewald_pw
  TYPE ( ewald_environment_type ),   POINTER :: ewald_env
  REAL (KIND=dp), DIMENSION ( :,: ), ALLOCATABLE :: f_nonbond
  LOGICAL,SAVE :: first_time=.true.
  LOGICAL                                  :: PM1, PM2, PM3, PM4
  TYPE(section_vals_type), POINTER         :: fist_section, print_section

  CALL timeset ( 'FORCE','I',' ',handle )
  ! Print Debug Info
  ! Da riaggiustare dopo con il nuovo input...
  PM1 = .FALSE.
  PM2 = .FALSE.
  PM3 = .FALSE.
  PM4 = .FALSE.
  iw = globenv%scr
  IF((globenv%ionode).AND.(globenv%print%level>0)) PM1= .TRUE.
  IF((globenv%ionode).AND.(globenv%print%level>1)) PM2= .TRUE.
  IF((globenv%ionode).AND.(globenv%print%level>2)) PM3= .TRUE.
  IF((globenv%ionode).AND.(globenv%print%level>3)) PM4= .TRUE.
  NULLIFY(fist_section)
  fist_section => section_vals_get_subs_vals(globenv%input_file,"FORCE_EVAL%MM",error=error)

! associcating local pointers
  CALL get_fist_env ( fist_env, ewald_pw = ewald_pw, ewald_env = ewald_env, &
                      local_particles = local_particles, particle_set = particle_set, &
                      atomic_kind_set = atomic_kind_set, molecule_set = molecule_set, &
                      local_molecules = local_molecules, thermo = thermo, &
                      molecule_kind_set = molecule_kind_set,  &
                      fist_nonbond_env = fist_nonbond_env, cell = cell )

  CALL ewald_env_get ( ewald_env, ewald_type = ewald_type )

  isos = 0
  natoms = SIZE ( particle_set )
  nlocal_particles = 0
  nparticle_kind = SIZE ( atomic_kind_set )
  DO ikind = 1, nparticle_kind
     nlocal_particles = nlocal_particles + local_particles%n_el(ikind) 
  ENDDO 

  ALLOCATE ( f_nonbond ( 3, natoms ), STAT = isos )
  IF ( isos /= 0 ) &
      CALL stop_memory ( 'fist_force', 'f_nonbond', 3 * natoms )

! reinitialize the gspace for the new box
  IF ( ( box_change ) .OR. ( first_time ) ) THEN
    CALL ewald_pw_grid_change ( ewald_pw, ewald_env )
  ENDIF

! first check with list_control to update neighbor lists

  CALL list_control ( atomic_kind_set, particle_set, local_particles, &
                      cell, fist_nonbond_env, globenv ) 

!
! initial force, energy and pressure tensor arrays
!
  DO i = 1, natoms
    particle_set ( i ) % f ( 1 ) = 0.0_dp
    particle_set ( i ) % f ( 2 ) = 0.0_dp
    particle_set ( i ) % f ( 3 ) = 0.0_dp
    particle_set ( i ) % rho = 0.0_dp
    particle_set ( i ) % f_embed = 0.0_dp
  ENDDO

  pv_bc = 0.0_dp
  pv_bond = 0.0_dp
  pv_bend = 0.0_dp
  pv_torsion = 0.0_dp
  pv_imptors = 0.0_dp
  pv_onef = 0.0_dp
  pv_urey_bradley = 0.0_dp
  pot_nonbond = 0.0_dp
  pot_manybody = 0.0_dp
  pot_bond = 0.0_dp
  pot_bend = 0.0_dp
  pot_torsion = 0.0_dp
  pot_imptors = 0.0_dp
  pot_onef = 0.0_dp
  pot_urey_bradley = 0.0_dp
  vg_coulomb = 0.0_dp
  thermo % pot = 0.0_dp
  virial % pv_virial = 0.0_dp
!
! get real-space non-bonded forces:
!
  IF (PM2) THEN
     WRITE(iw,'(A)')" FIST:: FORCES IN INPUT..."
     WRITE(iw,'(3f15.9)')((particle_set ( i ) % f ( j ),j=1,3), i = 1, SIZE(particle_set))
  END IF

!  compute density for EAM 
  CALL density_nonbond ( fist_nonbond_env, particle_set, cell, globenv )

! compute embedding function and manybody energy 
  CALL energy_manybody ( fist_nonbond_env, atomic_kind_set, &
                         local_particles, particle_set, pot_manybody, globenv )

! nonbond contribution + manybody forces
  CALL force_nonbond ( fist_nonbond_env, particle_set, cell, &
                       pot_nonbond, f_nonbond, pv_nonbond  )

  IF (PM2) THEN
     WRITE(iw,'(A)')" FIST:: NONBOND + R-SPACE ELECTROSTATIC FORCES ..."
     WRITE(iw,'(3f15.9)')f_nonbond
  END IF

!
! get g-space non-bonded forces:
!
  IF ( ewald_type /= 'NONE' ) THEN

     ! compute g-space part of the ewald sum
     SELECT CASE ( ewald_type )

     CASE ( "EWALD" )
        ! Parallelisation over atoms --> allocate local atoms
        ALLOCATE ( fg_coulomb ( 3, nlocal_particles ), STAT=isos )
        IF ( isos /= 0 ) &
        CALL stop_memory ( 'fist_force', 'fg_coulomb', 3 * nlocal_particles )
        fg_coulomb = 0.0_dp
        CALL ewald_evaluate ( ewald_env, ewald_pw, atomic_kind_set, particle_set, &
                              local_particles, fg_coulomb, vg_coulomb, pv_g )
        
     CASE ( "PME" )
        ! Parallelisation over grids --> allocate all atoms
        ALLOCATE ( fg_coulomb ( 3, natoms ), STAT=isos )
        IF ( isos /= 0 ) &
        CALL stop_memory ( 'fist_force', 'fg_coulomb', 3 * natoms )
        fg_coulomb = 0.0_dp
        CALL pme_evaluate ( ewald_env, ewald_pw, particle_set, &
                            vg_coulomb, fg_coulomb, pv_g )
        CALL mp_sum ( fg_coulomb, globenv % group ) 
     CASE ( "SPME" )
        ! Parallelisation over grids --> allocate all atoms
        ALLOCATE ( fg_coulomb ( 3, natoms ), STAT=isos )
        IF ( isos /= 0 ) &
        CALL stop_memory ( 'fist_force', 'fg_coulomb', 3 * natoms )
        fg_coulomb = 0.0_dp
        CALL spme_evaluate ( ewald_env, ewald_pw, particle_set, &
                             fg_coulomb, vg_coulomb, pv_g)
        CALL mp_sum ( fg_coulomb, globenv % group ) 
        
     CASE DEFAULT
        CALL stop_program ( "fist_force", "illegal value of ewald_type: "//&
                            ewald_type )
     END SELECT

     CALL bonded_correct_gaussian ( atomic_kind_set, local_particles, &
                                    particle_set, ewald_env, &
                                    thermo % e_bonded, pv_bc )
     IF ( ( first_time ) .AND. ( globenv % ionode ) ) &
     CALL ewald_print ( globenv % scr, cell, vg_coulomb, &
                       thermo % e_self, thermo % e_neut, thermo % e_bonded )

  END IF

  IF (PM2) THEN
     IF (ALLOCATED(fg_coulomb)) THEN
        WRITE(iw,'(A)')" FIST:: NONBONDED ELECTROSTATIC FORCES IN G-SPACE..."
        WRITE(iw,'(3f15.9)')((fg_coulomb(j,i) ,j=1,3),i = 1, SIZE(particle_set))
     END IF
  END IF

!
! get intramolecular forces
!
  IF ( PRESENT ( debug ) ) THEN
    CALL force_intra_control ( molecule_set, molecule_kind_set, local_molecules, &
                   particle_set, pot_bond, pot_bend, pot_urey_bradley, &
                   pot_torsion, pot_imptors, pot_onef, pv_bond, pv_bend, &
                   pv_urey_bradley, pv_torsion, pv_imptors, pv_onef,    &
                   debug % f_bond, debug % f_bend, debug % f_torsion,  &
                   debug % f_ub, debug % f_imptors, debug % f_onef )
  ELSE
    CALL force_intra_control ( molecule_set, molecule_kind_set, local_molecules, &
                   particle_set, pot_bond, pot_bend, pot_urey_bradley, &
                   pot_torsion, pot_imptors, pot_onef, pv_bond, pv_bend, &
                   pv_urey_bradley, pv_torsion, pv_imptors, pv_onef )
  ENDIF

  IF (PM2) THEN
     WRITE(iw,'(A)')" FIST:: CORRECTED BONDED ELECTROSTATIC FORCES + INTERNAL FORCES..."
     WRITE(iw,'(3f15.9)')((particle_set ( i ) % f ( j ),j=1,3), i = 1, SIZE(particle_set))
  END IF

!
! add up all the potential energies
!
  thermo % pot = pot_nonbond + pot_bond + pot_bend + pot_torsion + pot_onef + &
                 pot_imptors + pot_urey_bradley + pot_manybody
  CALL mp_sum ( thermo % pot, globenv % group )
!
! add g-space contributions if needed
  IF ( ewald_type /= 'NONE' ) THEN
! e_self, e_neut, and ebonded are already summed over all processors
! vg_coulomb is not calculated in parallel
     thermo % e_gspace = vg_coulomb
     thermo % pot = thermo % pot + thermo % e_self + thermo % e_neut/cell % deth
     thermo % pot = thermo % pot + vg_coulomb + thermo % e_bonded
  END IF
!
! add up all the forces
!
! nonbonded forces might be claculated for atoms not on this node
! ewald forces are strictly local -> sum only over pnode
! We first sum the forces in f_nonbond, this allows for a more efficient
! global sum in the parallel code and in the end copy them back to part
  ALLOCATE ( f_total ( 3,natoms ), STAT=isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'fist_force', 'f_total', 3 * natoms )
  f_total = 0.0_dp
  DO i = 1, natoms
     f_total ( 1, i ) = particle_set ( i ) % f ( 1 ) + f_nonbond ( 1, i )
     f_total ( 2, i ) = particle_set ( i ) % f ( 2 ) + f_nonbond ( 2, i )
     f_total ( 3, i ) = particle_set ( i ) % f ( 3 ) + f_nonbond ( 3, i )
  END DO

  IF (PM2) THEN
     WRITE(iw,'(A)')" FIST:: (1) INTERNAL + ELECTROSTATIC BONDED + NONBONDED"
     WRITE(iw,'(3f15.9)')((f_total ( j, i ),j=1,3), i = 1, natoms)
  END IF  

! Adding in the reciprocal forces: EWALD is a special
! case because of distrubted data
  SELECT CASE ( ewald_type )
  CASE ( "EWALD" )
     node = 0
     DO iparticle_kind = 1, nparticle_kind
        nparticle_local = local_particles % n_el ( iparticle_kind )
        DO iparticle_local = 1, nparticle_local
          ii = local_particles%list(iparticle_kind)%array(iparticle_local)
          node = node + 1
          f_total ( 1, ii ) = f_total ( 1, ii ) + fg_coulomb ( 1, node )
          f_total ( 2, ii ) = f_total ( 2, ii ) + fg_coulomb ( 2, node )
          f_total ( 3, ii ) = f_total ( 3, ii ) + fg_coulomb ( 3, node )
          IF ( PRESENT ( debug )  ) THEN
            debug % f_g ( 1, ii ) =  debug % f_g ( 1, ii ) + fg_coulomb ( 1, node )
            debug % f_g ( 2, ii ) =  debug % f_g ( 2, ii ) + fg_coulomb ( 2, node )
            debug % f_g ( 3, ii ) =  debug % f_g ( 3, ii ) + fg_coulomb ( 3, node )
          ENDIF
        END DO
     END DO
  END SELECT

  IF (PM2) THEN
     WRITE(iw,'(A)')" FIST:: (2) TOTAL FORCES (1) + ELECTROSTATIC FORCES"
     WRITE(iw,'(3f15.9)')((f_total ( j, i ),j=1,3), i = 1, natoms)
  END IF  

!
! add up all the pressure tensors
!
  IF ( ewald_type == 'NONE' ) THEN
     virial % pv_virial = pv_nonbond + pv_bond + pv_bend + &
                          pv_torsion + pv_imptors + pv_urey_bradley
     CALL mp_sum ( virial % pv_virial, globenv % group )
  ELSE
     ident = 0.0_dp
     DO i = 1, 3
        ident ( i, i ) = 1.0_dp
     END DO

     virial % pv_virial = pv_nonbond + pv_bond + pv_bend + &
                          pv_torsion + pv_imptors + pv_urey_bradley +pv_bc

     CALL mp_sum ( virial % pv_virial, globenv % group )

     virial % pv_virial = virial % pv_virial + &
                          ident * thermo % e_neut / cell % deth
     virial % pv_virial = virial % pv_virial + pv_g
  END IF

  IF ( first_time ) CALL virial_evaluate ( atomic_kind_set, particle_set, &
                                           local_particles, virial, globenv % group )

! Sum total forces
  CALL mp_sum ( f_total, globenv % group )

! Assign to particle_set
  SELECT CASE ( ewald_type )
  CASE ( "SPME", "PME" )
    DO i = 1, natoms
      particle_set ( i ) % f ( 1 ) = f_total ( 1, i ) + fg_coulomb ( 1, i )
      particle_set ( i ) % f ( 2 ) = f_total ( 2, i ) + fg_coulomb ( 2, i )
      particle_set ( i ) % f ( 3 ) = f_total ( 3, i ) + fg_coulomb ( 3, i )
    END DO
  CASE ( "EWALD", "NONE" )
    DO i = 1, natoms
      particle_set ( i ) % f ( 1 ) = f_total ( 1, i ) 
      particle_set ( i ) % f ( 2 ) = f_total ( 2, i ) 
      particle_set ( i ) % f ( 3 ) = f_total ( 3, i ) 
    END DO
  CASE DEFAULT
    CALL stop_program ( "fist_force", "illegal value of ewald_type: "//&
                         ewald_type )
  END SELECT
  ! Check if we have to delete some gradients.. i.e. some atom is considered fixed...
  IF (ANY(particle_set%is_fixed)) THEN
     DO i = 1, natoms
        IF (particle_set(i)%is_fixed) particle_set(i)%f = 0.0_dp
     END DO
  END IF
  IF (PM2) THEN
     WRITE(iw,'(A)')" FIST:: (3) TOTAL FORCES - THE END..."
     WRITE(iw,'(3f15.9)')((particle_set ( i ) % f ( j ),j=1,3), i = 1, natoms)
     WRITE(iw,'(A,f15.9)')"Energy after FIST calculation.. exiting now ::",thermo%pot
  END IF  

!
! if we are doing debugging, check if variables are present and assign
!

  IF ( PRESENT ( debug )  ) THEN
     CALL mp_sum ( pot_nonbond, globenv % group )
     debug % pot_nonbond = pot_nonbond
     CALL mp_sum ( pot_bond, globenv % group )
     debug % pot_bond = pot_bond
     CALL mp_sum ( pot_bend, globenv % group )
     debug % pot_bend = pot_bend
     CALL mp_sum ( pot_torsion, globenv % group )
     debug % pot_torsion = pot_torsion
     CALL mp_sum ( pot_imptors, globenv % group )
     debug % pot_imptors = pot_imptors
     CALL mp_sum ( pot_onef, globenv % group )
     debug % pot_onef = pot_onef
     CALL mp_sum ( pot_urey_bradley, globenv % group )
     debug % pot_urey_bradley = pot_urey_bradley
     CALL mp_sum ( f_nonbond, globenv % group )
     debug % f_nonbond = f_nonbond
     CALL mp_sum ( pv_nonbond, globenv % group )
     debug % pv_nonbond = pv_nonbond
     CALL mp_sum ( pv_bond, globenv % group )
     debug % pv_bond = pv_bond
     CALL mp_sum ( pv_bend, globenv % group )
     debug % pv_bend = pv_bend
     CALL mp_sum ( pv_torsion, globenv % group )
     debug % pv_torsion = pv_torsion
     CALL mp_sum ( pv_imptors, globenv % group )
     debug % pv_imptors = pv_imptors
     CALL mp_sum ( pv_onef, globenv % group )
     debug % pv_onef = pv_onef
     CALL mp_sum ( pv_urey_bradley, globenv % group )
     debug % pv_ub = pv_urey_bradley
     SELECT CASE ( ewald_type )
     CASE ( 'PME', 'SPME' )
        debug % pot_g = vg_coulomb
        debug % pv_g = pv_g
        debug % f_g = fg_coulomb
     CASE ( 'EWALD' )
        debug % pot_g = vg_coulomb
        debug % pv_g = pv_g
        ! debug % fg is computed above because of
        ! distributed particle data
     CASE default
        debug % pot_g = 0.0_dp
        debug % f_g = 0.0_dp
        debug % pv_g = 0.0_dp
     END SELECT
  END IF

  ! print properties if requested
  print_section => section_vals_get_subs_vals(fist_section,"PRINT",error=error)
  CALL print_fist(print_section, atomic_kind_set, particle_set, cell, error)

! deallocating all local variables
  IF ( ALLOCATED ( fg_coulomb ) ) THEN
    DEALLOCATE ( fg_coulomb, STAT = isos )
    IF ( isos /= 0 ) CALL stop_memory ( 'fist_force', 'fg_coulomb' )
  END IF
  IF ( ALLOCATED ( f_total ) ) THEN
    DEALLOCATE ( f_total, STAT = isos )
    IF ( isos /= 0 ) CALL stop_memory ( 'fist_force', 'f_total' )
  END IF
  DEALLOCATE ( f_nonbond, STAT = isos )
  IF ( isos /= 0 ) &
      CALL stop_memory ( 'fist_force', 'f_nonbond', 3 * natoms )
  first_time = .false.
  CALL timestop ( zero, handle )

END SUBROUTINE fist_force_control

!******************************************************************************
!!****** fist_force/fist_force_control [1.1] *
!!
!!   NAME
!!    print_fist
!!
!!   FUNCTION
!!     Print properties number according the requests in input file
!!
!!   AUTHOR
!!     Teo
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!!******************************************************************************
SUBROUTINE print_fist(print_section, atomic_kind_set, particle_set, cell, error)
  IMPLICIT NONE
  ! Arguments
  TYPE(section_vals_type), POINTER         :: print_section
  TYPE ( particle_type ), DIMENSION(:), POINTER :: particle_set 
  TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set  
  TYPE ( cell_type ), POINTER                :: cell
  TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                               :: error
  ! Local Variables
  TYPE(cp_logger_type), POINTER :: logger
  TYPE(section_vals_type), POINTER :: print_key

  logger => cp_error_get_logger(error)
  print_key => section_vals_get_subs_vals(print_section,"dipole",error=error)  
  IF (cp_should_output(logger%iter_info,print_key,error=error).AND.&
       logger%para_env%mepos==logger%para_env%source) THEN 
!     CALL fist_dipole(atomic_kind_set, particle_set, cell, error)
  END IF

END SUBROUTINE print_fist
  
!******************************************************************************
END MODULE fist_force
!******************************************************************************
