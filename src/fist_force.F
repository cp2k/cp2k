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
  USE cp_output_handling, ONLY: cp_print_key_should_output, cp_p_file
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
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
  USE timings, ONLY : timeset, timestop
  USE virial_methods, ONLY : virial_evaluate
  USE virial_types, ONLY : virial_type
  USE physcon, ONLY: debye
  USE input_section_types, ONLY: section_type, section_create, &
       section_add_keyword, section_add_subsection, section_release,&
       section_describe_html, section_describe_index_html, section_vals_type,&
       section_vals_get_subs_vals, section_vals_get, section_vals_val_get
  USE input_constants,                 ONLY: do_ewald_none,&
                                             do_ewald_ewald,&
                                             do_ewald_pme,&
                                             do_ewald_spme
  USE cp_output_handling,              ONLY: cp_print_key_unit_nr, &
                                             cp_print_key_finished_output
  USE checkpoint_handler,              ONLY: write_checkpoint_information
#include "cp_common_uses.h"
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: fist_force_control, debug_variables_type
  CHARACTER(len=*), PARAMETER :: moduleN="fist_force"

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
  TYPE ( global_environment_type ), pointer :: globenv
  LOGICAL, INTENT ( IN ) :: box_change
  TYPE ( debug_variables_type ), INTENT ( INOUT ), OPTIONAL :: debug
  TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                               :: error
! Locals
    CHARACTER(len=*), PARAMETER :: routineN="fist_force_control",&
         routineP=moduleN//":"//routineN  
  INTEGER :: id, i, ii, natoms, nnodes, handle, stat, n, nlocal_particles
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
  INTEGER                           :: ewald_type
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
  TYPE(section_vals_type), POINTER         :: fist_section, print_section
  LOGICAL :: failure
  TYPE(cp_logger_type), POINTER            :: logger

  failure = .FALSE.
  NULLIFY(logger)
  logger => cp_error_get_logger(error)
  CALL write_checkpoint_information("Entering "//routineN,globenv,error)
  iw = cp_print_key_unit_nr(logger,globenv%input_file,"FORCE_EVAL%MM%PRINT%DERIVATIVES",&
       extension=".mmLog",error=error)
  CALL timeset(routineN,'I','',handle)

  NULLIFY(fist_section,thermo,particle_set,local_molecules,local_particles,atomic_kind_set,&
          molecule_kind_set,molecule_set,fist_nonbond_env,cell,ewald_pw,ewald_env,fist_section,print_section)
  fist_section => section_vals_get_subs_vals(globenv%input_file,"FORCE_EVAL%MM",error=error)

! associcating local pointers
  CALL get_fist_env ( fist_env, ewald_pw = ewald_pw, ewald_env = ewald_env, &
                      local_particles = local_particles, particle_set = particle_set, &
                      atomic_kind_set = atomic_kind_set, molecule_set = molecule_set, &
                      local_molecules = local_molecules, thermo = thermo, &
                      molecule_kind_set = molecule_kind_set,  &
                      fist_nonbond_env = fist_nonbond_env, cell = cell )

  CALL ewald_env_get ( ewald_env, ewald_type = ewald_type )

  natoms = SIZE ( particle_set )
  nlocal_particles = 0
  nparticle_kind = SIZE ( atomic_kind_set )
  DO ikind = 1, nparticle_kind
     nlocal_particles = nlocal_particles + local_particles%n_el(ikind) 
  ENDDO 

  ALLOCATE ( f_nonbond ( 3, natoms ), STAT = stat )
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
! reinitialize the gspace for the new box
  IF ( ( box_change ) .OR. ( first_time ) ) THEN
    CALL ewald_pw_grid_change ( ewald_pw, ewald_env )
  ENDIF

! first check with list_control to update neighbor lists

  CALL list_control ( atomic_kind_set, particle_set, local_particles, &
                      cell, fist_nonbond_env, globenv, error) 

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
  IF (iw>0) THEN
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

  IF (iw>0) THEN
     WRITE(iw,'(A)')" FIST:: NONBOND + R-SPACE ELECTROSTATIC FORCES ..."
     WRITE(iw,'(3f15.9)')f_nonbond
  END IF

!
! get g-space non-bonded forces:
!
  IF ( ewald_type /= do_ewald_none ) THEN

     ! compute g-space part of the ewald sum
     SELECT CASE ( ewald_type )

     CASE ( do_ewald_ewald )
        ! Parallelisation over atoms --> allocate local atoms
        ALLOCATE ( fg_coulomb ( 3, nlocal_particles ), STAT=stat )
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        fg_coulomb = 0.0_dp
        CALL ewald_evaluate ( ewald_env, ewald_pw, atomic_kind_set, particle_set, &
                              local_particles, fg_coulomb, vg_coulomb, pv_g )
        
     CASE ( do_ewald_pme )
        ! Parallelisation over grids --> allocate all atoms
        ALLOCATE ( fg_coulomb ( 3, natoms ), STAT=stat )
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        fg_coulomb = 0.0_dp
        CALL pme_evaluate ( ewald_env, ewald_pw, particle_set, &
                            vg_coulomb, fg_coulomb, pv_g ,error)
        CALL mp_sum ( fg_coulomb, globenv % group ) 
     CASE ( do_ewald_spme )
        ! Parallelisation over grids --> allocate all atoms
        ALLOCATE ( fg_coulomb ( 3, natoms ), STAT=stat )
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        fg_coulomb = 0.0_dp
        CALL spme_evaluate ( ewald_env, ewald_pw, particle_set, &
                             fg_coulomb, vg_coulomb, pv_g)
        CALL mp_sum ( fg_coulomb, globenv % group ) 
        
     END SELECT

     CALL bonded_correct_gaussian ( atomic_kind_set, local_particles, &
                                    particle_set, ewald_env, &
                                    thermo % e_bonded, pv_bc )
     IF ( ( first_time ) .AND. ( globenv % ionode ) ) &
     CALL ewald_print ( globenv % scr, cell, vg_coulomb, &
                       thermo % e_self, thermo % e_neut, thermo % e_bonded )

  END IF

  IF (iw>0) THEN
     IF (ALLOCATED(fg_coulomb)) THEN
        WRITE(iw,'(A)')" FIST:: NONBONDED ELECTROSTATIC FORCES IN G-SPACE..."
        WRITE(iw,'(3f15.9)')((fg_coulomb(j,i) ,j=1,3),i = 1, SIZE(fg_coulomb,2))
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

  IF (iw>0) THEN
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
  IF ( ewald_type /= do_ewald_none ) THEN
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
  ALLOCATE ( f_total ( 3,natoms ), STAT=stat )
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  f_total = 0.0_dp
  DO i = 1, natoms
     f_total ( 1, i ) = particle_set ( i ) % f ( 1 ) + f_nonbond ( 1, i )
     f_total ( 2, i ) = particle_set ( i ) % f ( 2 ) + f_nonbond ( 2, i )
     f_total ( 3, i ) = particle_set ( i ) % f ( 3 ) + f_nonbond ( 3, i )
  END DO

  IF (iw>0) THEN
     WRITE(iw,'(A)')" FIST:: (1) INTERNAL + ELECTROSTATIC BONDED + NONBONDED"
     WRITE(iw,'(3f15.9)')((f_total ( j, i ),j=1,3), i = 1, natoms)
  END IF  

! Adding in the reciprocal forces: EWALD is a special
! case because of distrubted data
  SELECT CASE ( ewald_type )
  CASE ( do_ewald_ewald )
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


  IF (iw>0) THEN
     WRITE(iw,'(A)')" FIST:: (2) TOTAL FORCES (1) + ELECTROSTATIC FORCES"
     WRITE(iw,'(3f15.9)')((f_total ( j, i ),j=1,3), i = 1, natoms)
  END IF  

!
! add up all the pressure tensors
!
  IF ( ewald_type == do_ewald_none ) THEN
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


! Sum total forces
  CALL mp_sum ( f_total, globenv % group )

! Assign to particle_set
  SELECT CASE ( ewald_type )
  CASE ( do_ewald_spme, do_ewald_pme )
    DO i = 1, natoms
      particle_set ( i ) % f ( 1 ) = f_total ( 1, i ) + fg_coulomb ( 1, i )
      particle_set ( i ) % f ( 2 ) = f_total ( 2, i ) + fg_coulomb ( 2, i )
      particle_set ( i ) % f ( 3 ) = f_total ( 3, i ) + fg_coulomb ( 3, i )
    END DO
  CASE ( do_ewald_ewald, do_ewald_none )
    DO i = 1, natoms
      particle_set ( i ) % f ( 1 ) = f_total ( 1, i ) 
      particle_set ( i ) % f ( 2 ) = f_total ( 2, i ) 
      particle_set ( i ) % f ( 3 ) = f_total ( 3, i ) 
    END DO
  END SELECT
  IF (iw>0) THEN
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
     CASE ( do_ewald_spme, do_ewald_pme )
        debug % pot_g = vg_coulomb
        debug % pv_g = pv_g
        debug % f_g = fg_coulomb
     CASE ( do_ewald_ewald )
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
    DEALLOCATE ( fg_coulomb, STAT = stat )
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF ( ALLOCATED ( f_total ) ) THEN
    DEALLOCATE ( f_total, STAT = stat )
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  DEALLOCATE ( f_nonbond, STAT = stat )
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  first_time = .false.
  CALL timestop(0.0_dp,handle)
  CALL write_checkpoint_information("Leaving "//routineN,globenv,error)
  CALL cp_print_key_finished_output(iw,logger,globenv%input_file,&
       "FORCE_EVAL%MM%PRINT%DERIVATIVES",error=error)

END SUBROUTINE fist_force_control

!******************************************************************************
!!****** fist_force/print_fist [1.1] *
!!
!!   NAME
!!    print_fist
!!
!!   FUNCTION
!!     Print properties number according the requests in input file
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     [01.2006] created
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
  INTEGER :: unit_nr

  NULLIFY(logger,print_key)
  logger => cp_error_get_logger(error)
  print_key => section_vals_get_subs_vals(print_section,"dipole",error=error)  
  IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),&
       cp_p_file)) THEN
     unit_nr=cp_print_key_unit_nr(logger,print_section,"dipole",&
          extension=".data",middle_name="MM_DIPOLE",log_filename=.FALSE.,error=error)
     CALL fist_dipole(atomic_kind_set, particle_set, cell, unit_nr, error)
     CALL cp_print_key_finished_output(unit_nr,logger,print_key,error=error)     
  END IF

END SUBROUTINE print_fist

!******************************************************************************
!!****** fist_force/fist_dipole [1.1] *
!!
!!   NAME
!!    fist_dipole
!!
!!   FUNCTION
!!     Evaluates the Dipole of a classical charge distribution (point-like)
!!     using the berry phase formalism
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     [01.2006] created
!!
!!   SOURCE
!!******************************************************************************
SUBROUTINE fist_dipole(atomic_kind_set, particle_set, cell, unit_nr, error)
  IMPLICIT NONE
  ! Arguments
  TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
  TYPE ( particle_type ), DIMENSION(:), POINTER :: particle_set 
  TYPE (cell_type), POINTER :: cell
  INTEGER, INTENT(IN) :: unit_nr
  TYPE(cp_error_type), INTENT(INOUT), &
       OPTIONAL                               :: error
  ! Local Variables
  INTEGER                          :: i, j, iparticle_kind
  REAL(KIND=dp)                    :: dipole(3), theta, gvec(3), q
  COMPLEX(KIND=dp), DIMENSION(3)   :: ggamma, zeta
  INTEGER, DIMENSION(:), POINTER   :: atom_list
  TYPE(atomic_kind_type), POINTER          :: atomic_kind

  NULLIFY(atomic_kind)
  DO i = 1,3
     gvec(i) = 2.0_dp*pi/cell%hmat(i,i)
  ENDDO

  ggamma = CMPLX(0.0_dp,0.0_dp,KIND=dp)
  DO iparticle_kind = 1, SIZE(atomic_kind_set)     
     atomic_kind => atomic_kind_set(iparticle_kind)
     CALL get_atomic_kind(atomic_kind=atomic_kind, qeff=q,&
          atom_list=atom_list)
     
     DO i = 1, SIZE(atom_list)
        DO j = 1, 3
           theta = particle_set(i)%r(j)*gvec(j)
           zeta(j) = q * CMPLX(COS(theta),SIN(theta),KIND=dp)
        END DO
        ggamma = ggamma + zeta
     ENDDO
  END DO
  
  dipole = 0.0_dp
  DO j = 1, 3
     IF (REAL(ggamma(j),KIND=dp) /= 0.0_dp) &
     dipole(j) = 1.0_dp/gvec(j) * ATAN( REAL(AIMAG(ggamma(j)),KIND=dp)/REAL(ggamma(j),KIND=dp) )  
  END DO
  
  IF (unit_nr>0) THEN
     WRITE(unit_nr,'(1X,A,3F12.6)')"DIPOLE BERRY PHASE ( A.U.)|",dipole
     WRITE(unit_nr,'(1X,A,3F12.6)')"DIPOLE BERRY PHASE (Debye)|",dipole*debye
  END IF
  
END SUBROUTINE fist_dipole

!******************************************************************************
END MODULE fist_force
!******************************************************************************
