!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C ) 2000  CP2K developers group                                !
!-----------------------------------------------------------------------------!
!!****** cp2k/fist_force [1.0] *
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
  USE fist_nonbond_force, ONLY : force_nonbond, bonded_correct_gaussian
  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl, default_string_length
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
  USE virial_types, ONLY : virial_type

  PRIVATE
  PUBLIC :: fist_force_control, debug_variables_type

  TYPE debug_variables_type
     REAL ( dbl ) :: pot_nonbond, pot_g, pot_bond, pot_bend, pot_torsion
     REAL ( dbl ) :: pot_imptors, pot_onef, pot_urey_bradley
     REAL ( dbl ), DIMENSION ( :, : ), POINTER :: &
          f_nonbond, f_g, f_bond, f_bend, f_torsion, f_imptors, f_onef, f_ub
     REAL ( dbl ), DIMENSION ( 3, 3 ) :: pv_nonbond, pv_g, pv_bond, &
          pv_bend, pv_torsion, pv_imptors, pv_onef, pv_ub
  END TYPE debug_variables_type

  REAL ( dbl ), DIMENSION ( :,: ), ALLOCATABLE :: f_nonbond

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

SUBROUTINE fist_force_control ( fist_env, virial, globenv, box_change, debug )

  IMPLICIT NONE

! Arguments
  TYPE ( fist_environment_type ), POINTER :: fist_env
  TYPE ( virial_type ), POINTER :: virial
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv
  LOGICAL, INTENT ( IN ) :: box_change
  TYPE ( debug_variables_type ), INTENT ( INOUT ), OPTIONAL :: debug

! Locals
  INTEGER :: id, i, ii, natoms, nnodes, handle, isos, n, nlocal_particles
  INTEGER :: ikind, iparticle_kind, nparticle_kind, nparticle_local, iparticle_local
  INTEGER :: node
  REAL ( dbl ) :: pot_nonbond, pot_bond, pot_bend, pot_urey_bradley,  &
                  pot_imptors, pot_torsion, pot_onef, vg_coulomb
  REAL ( dbl ), DIMENSION ( 3,3 ) :: pv_nonbond, pv_bond, pv_bend,  &
                 pv_torsion, pv_urey_bradley, pv_imptors, pv_onef
  REAL ( dbl ), DIMENSION ( :,: ), ALLOCATABLE :: fg_coulomb
  REAL ( dbl ), DIMENSION ( :,: ), ALLOCATABLE :: f_total
  REAL ( dbl ), DIMENSION ( 3,3 ) :: pv_g, ident
  REAL ( dbl ), DIMENSION ( 3,3 ) :: pv_bc
  CHARACTER ( len = default_string_length ) :: ewald_type
  TYPE ( fist_energy_type ), POINTER :: thermo
  TYPE ( particle_type ), POINTER :: particle_set ( : )
  TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
  TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
  TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
  TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
  TYPE ( fist_nonbond_env_type ), POINTER    :: fist_nonbond_env
  TYPE ( cell_type ), POINTER          :: cell
  TYPE ( ewald_pw_type ),   POINTER          :: ewald_pw
  TYPE ( ewald_environment_type ),   POINTER :: ewald_env

  LOGICAL :: first_time


  CALL timeset ( 'FORCE','I',' ',handle )

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

  first_time = .NOT. ALLOCATED ( f_nonbond )

  IF ( .NOT. ALLOCATED ( f_nonbond )  ) THEN
    ALLOCATE ( f_nonbond ( 3, natoms ), STAT = isos )
    IF ( isos /= 0 ) &
       CALL stop_memory ( 'fist_force', 'f_nonbond', 3 * natoms )
  ELSE IF ( SIZE ( f_nonbond, 2  ) < natoms ) THEN
    DEALLOCATE ( f_nonbond, STAT = isos )
    IF ( isos /= 0 ) CALL stop_memory ( 'fist_force', 'f_nonbond' )
    ALLOCATE ( f_nonbond ( 3,natoms ), STAT = isos )
    IF ( isos /= 0 ) &
    CALL stop_memory ( 'fist_force', 'f_nonbond', 3 * natoms )
! change in system, do a full reinitialisation
    first_time = .TRUE.
  END IF


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
     particle_set ( i ) % f ( 1 ) = 0.0_dbl
     particle_set ( i ) % f ( 2 ) = 0.0_dbl
     particle_set ( i ) % f ( 3 ) = 0.0_dbl
  END DO
  pv_bc = 0._dbl
  pv_bond = 0.0_dbl
  pv_bend = 0.0_dbl
  pv_urey_bradley = 0.0_dbl
  pv_torsion = 0.0_dbl
  pv_imptors = 0.0_dbl
  pv_onef = 0.0_dbl
  pot_nonbond = 0._dbl
  pot_bond = 0._dbl
  pot_bend = 0._dbl
  pot_urey_bradley = 0._dbl
  pot_torsion = 0._dbl
  pot_imptors = 0._dbl
  pot_onef = 0._dbl
  vg_coulomb = 0._dbl
  thermo % pot = 0.0_dbl
  virial % pv_virial = 0.0_dbl
!
! get real-space non-bonded forces:
!
  f_nonbond = zero
  CALL force_nonbond ( fist_nonbond_env, particle_set, cell, &
                       pot_nonbond, f_nonbond, pv_nonbond  )

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
        CALL ewald_evaluate ( ewald_env, ewald_pw, atomic_kind_set, particle_set, &
                              local_particles, fg_coulomb, vg_coulomb, pv_g )
        
     CASE ( "PME" )
        ! Parallelisation over grids --> allocate all atoms
        ALLOCATE ( fg_coulomb ( 3, natoms ), STAT=isos )
        IF ( isos /= 0 ) &
        CALL stop_memory ( 'fist_force', 'fg_coulomb', 3 * natoms )
        fg_coulomb = 0._dbl
        CALL pme_evaluate ( ewald_env, ewald_pw, particle_set, &
                            vg_coulomb, fg_coulomb, pv_g )
        CALL mp_sum ( fg_coulomb, globenv % group ) 
     CASE ( "SPME" )
        ! Parallelisation over grids --> allocate all atoms
        ALLOCATE ( fg_coulomb ( 3, natoms ), STAT=isos )
        IF ( isos /= 0 ) &
        CALL stop_memory ( 'fist_force', 'fg_coulomb', 3 * natoms )
        fg_coulomb = 0._dbl
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

!
! get intramolecular forces
!
  IF ( PRESENT ( debug )  ) THEN
     CALL force_intra_control ( molecule_set, molecule_kind_set, local_molecules, &
                      particle_set, pot_bond, pot_bend, pot_urey_bradley, &
                      pot_torsion, pot_imptors, pot_onef, pv_bond, pv_bend, &
                      pv_urey_bradley, pv_torsion, pv_imptors, pv_onef,    &
                      debug % f_bond, debug % f_bend, debug % f_torsion,  &
                      debug % f_ub, debug % f_imptors, debug % f_onef )
  ELSE
     CALL force_intra_control ( molecule_set, molecule_kind_set, local_molecules, &
                               particle_set, pot_bond, pot_bend, pot_urey_bradley, &
                               pot_torsion, pot_imptors, pot_onef, pv_bond, pv_bend,     &
                               pv_urey_bradley, pv_torsion, pv_imptors,pv_onef )
  END IF

!
! add up all the potential energies
!
  IF ( ewald_type == 'NONE' ) THEN
     thermo % pot = pot_nonbond + pot_bond + pot_bend + pot_torsion + &
                    pot_imptors + pot_urey_bradley
     CALL mp_sum ( thermo % pot, globenv % group )
  ELSE
     thermo % pot = pot_nonbond + pot_bond + pot_bend + pot_torsion + &
                    pot_imptors + pot_urey_bradley
     thermo % e_gspace = vg_coulomb
     CALL mp_sum ( thermo % pot,globenv % group )
! e_self, e_neut, and ebonded are already summed over all processors
! vg_coulomb is not calculated in parallel
     thermo % pot = thermo % pot + thermo % e_self + thermo % e_neut/cell % deth
     thermo % pot = thermo % pot + vg_coulomb + thermo % e_bonded
  END IF

! add up all the forces
! nonbonded forces might be claculated for atoms not on this node
! ewald forces are strictly local -> sum only over pnode
! We first sum the forces in f_nonbond, this allows for a more efficient
! global sum in the parallel code and in the end copy them back to part
  ALLOCATE ( f_total ( 3,natoms ), STAT=isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'fist_force', 'f_total', 3 * natoms )
  f_total = 0._dbl
  DO i = 1, natoms
     f_total ( 1, i ) = particle_set ( i ) % f ( 1 ) + f_nonbond ( 1, i )
     f_total ( 2, i ) = particle_set ( i ) % f ( 2 ) + f_nonbond ( 2, i )
     f_total ( 3, i ) = particle_set ( i ) % f ( 3 ) + f_nonbond ( 3, i )
  END DO
  IF ( ewald_type /= 'NONE' ) THEN
     SELECT CASE ( ewald_type )
     CASE ( "EWALD" )
        node = 0
        DO iparticle_kind = 1, nparticle_kind
           nparticle_local = local_particles % n_el ( iparticle_kind )
           DO iparticle_local = 1, nparticle_local
             node = node + 1
             ii = local_particles%list(iparticle_kind)%array(iparticle_local)
             f_total ( 1, ii ) = f_total ( 1, ii ) + fg_coulomb ( 1, node )
             f_total ( 2, ii ) = f_total ( 2, ii ) + fg_coulomb ( 2, node )
             f_total ( 3, ii ) = f_total ( 3, ii ) + fg_coulomb ( 3, node )
           END DO
        END DO
     CASE ( "PME" )
!        f_total ( :, : ) = f_total ( :, : ) + fg_coulomb ( :, : )
     CASE ( "SPME" )
!        f_total ( :, : ) = f_total ( :, : ) + fg_coulomb ( :, : )
     CASE DEFAULT
        CALL stop_program ( "fist_force", "illegal value of ewald_type: "//&
                             ewald_type )
     END SELECT
  END IF

! add up all the pressure tensors
  IF ( ewald_type == 'NONE' ) THEN
     virial % pv_virial = pv_nonbond + pv_bond + pv_bend + &
                          pv_torsion + pv_imptors + pv_urey_bradley
     CALL mp_sum ( virial % pv_virial, globenv % group )
  ELSE
     ident = 0.0_dbl
     DO i = 1, 3
        ident ( i, i ) = 1.0_dbl
     END DO

     virial % pv_virial = pv_nonbond + pv_bond + pv_bend + &
                          pv_torsion + pv_imptors + pv_urey_bradley +pv_bc
     CALL mp_sum ( virial % pv_virial, globenv % group )

     virial % pv_virial = virial % pv_virial + &
                                   ident * thermo % e_neut / cell % deth
     virial % pv_virial = virial % pv_virial + pv_g
  END IF
!
! if we are doing debugging, check if variables are present and assign
!
  IF ( PRESENT ( debug )  ) THEN
     debug % pot_bond = pot_bond
     debug % pot_bend = pot_bend
     debug % pot_torsion = pot_torsion
     debug % pot_imptors = pot_imptors
     debug % pot_nonbond = pot_nonbond
     debug % f_nonbond = f_nonbond
     debug % pv_nonbond = pv_nonbond
     debug % pv_bond = pv_bond
     debug % pv_bend = pv_bend
     debug % pv_torsion = pv_torsion
     debug % pv_imptors = pv_imptors
     IF ( ewald_type /= 'NONE' ) THEN
        debug % pot_g = vg_coulomb
        debug % f_g = fg_coulomb
        debug % pv_g = pv_g
     ELSE
        debug % pot_g = 0.0_dbl
        debug % f_g = 0.0_dbl
        debug % pv_g = 0.0_dbl
     END IF
  END IF

  CALL mp_sum ( f_total, globenv % group )

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

! deallocating all local variables
  IF ( ALLOCATED ( fg_coulomb ) ) THEN
    DEALLOCATE ( fg_coulomb, STAT = isos )
    IF ( isos /= 0 ) CALL stop_memory ( 'fist_force', 'fg_coulomb' )
  END IF
  IF ( ALLOCATED ( f_total ) ) THEN
    DEALLOCATE ( f_total, STAT = isos )
    IF ( isos /= 0 ) CALL stop_memory ( 'fist_force', 'f_total' )
  END IF

  CALL timestop ( zero, handle )

END SUBROUTINE fist_force_control

!******************************************************************************
END MODULE fist_force
!******************************************************************************
