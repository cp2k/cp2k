!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C ) 2000  CP2K developers group                                !
!-----------------------------------------------------------------------------!
MODULE fist_force
  
  USE dg_types, ONLY : dg_type
  USE ewald_parameters_types, ONLY : ewald_parameters_type
  USE ewalds, ONLY : ewald_gaussian, ewald_initialize
  USE fist_intra_force, ONLY : force_intra_control
  USE fist_nonbond_force, ONLY : force_nonbond, bonded_correct_gaussian
  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl
  USE linklist_control, ONLY : list_control
  USE mathconstants, ONLY : pi, zero
  USE md, ONLY : thermodynamic_type
  USE molecule_types, ONLY : molecule_structure_type, particle_node_type, &
       linklist_neighbor, linklist_exclusion, linklist_bonds, linklist_bends
  USE message_passing, ONLY : mp_sum
  USE pair_potential, ONLY : potentialparm_type
  USE particle_types, ONLY : particle_type
  USE pme, ONLY: pme_evaluate
  USE pw_grid_types, ONLY : pw_grid_type
  USE pw_grids, ONLY : pw_grid_change
  USE simulation_cell, ONLY : cell_type, get_hinv
  USE stop_program, ONLY : stop_prg, stop_memory
  USE timings, ONLY : timeset, timestop
  
  PRIVATE
  PUBLIC :: force_control, debug_variables_type
  
  TYPE debug_variables_type
     REAL ( dbl ) :: pot_nonbond, pot_g, pot_bond, pot_bend
     REAL ( dbl ), DIMENSION ( :, : ), POINTER :: &
          f_nonbond, f_g, f_bond, f_bend
     REAL ( dbl ), DIMENSION ( 3, 3 ) :: pv_nonbond, pv_g, pv_bond, pv_bend
  END TYPE debug_variables_type
  
CONTAINS

!******************************************************************************

!! Calculates the total potential energy, total force, and the
!! total pressure tensor from the potentials

SUBROUTINE force_control ( molecule, pnode, part, box, thermo, &
     potparm, ewald_param, ensemble, fc_global, debug )

  IMPLICIT NONE

! Arguments
  TYPE ( particle_type ), DIMENSION ( : ), INTENT ( INOUT ) :: part
  TYPE ( particle_node_type ), DIMENSION ( : ), INTENT ( INOUT ) :: pnode
  TYPE ( molecule_structure_type ), DIMENSION ( : ), INTENT ( IN ) :: molecule
  TYPE ( cell_type ), INTENT ( INOUT ) :: box
  TYPE ( thermodynamic_type ), INTENT ( INOUT ) :: thermo
  TYPE ( potentialparm_type ), DIMENSION ( :,: ), INTENT ( IN ) :: potparm
  TYPE ( ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param
  CHARACTER ( LEN = * ), INTENT ( IN ) :: ensemble
  TYPE ( global_environment_type ), INTENT ( IN ) :: fc_global
  TYPE ( debug_variables_type ), INTENT ( OUT ), OPTIONAL :: debug
  
! Locals
  INTEGER :: id, i, natoms, nnodes, handle, isos
  REAL ( dbl ) :: pot_nonbond, pot_bond, pot_bend, vg_coulomb
  REAL ( dbl ), DIMENSION ( :,: ), ALLOCATABLE, SAVE :: f_nonbond
  REAL ( dbl ), DIMENSION ( 3,3 ) :: pv_nonbond, pv_bond, pv_bend
  REAL ( dbl ), DIMENSION ( :,: ), ALLOCATABLE :: fg_coulomb
  REAL ( dbl ), DIMENSION ( 3,3 ) :: pv_g, ident
  REAL ( dbl ), DIMENSION ( 3,3 ) :: pv_bc
  TYPE ( pw_grid_type ), SAVE :: grid_s, grid_b, grid_ewald
  TYPE ( dg_type ), SAVE :: dg
  LOGICAL :: first_time
  
!------------------------------------------------------------------------------
  
  CALL timeset ( 'FORCE','I',' ',handle )
  
  nnodes = SIZE ( pnode )
  natoms = SIZE ( part )
  isos = 0
  first_time = .NOT. ALLOCATED ( f_nonbond )
  IF ( .NOT. ALLOCATED ( f_nonbond )  ) &
       ALLOCATE ( f_nonbond ( 3,natoms ), STAT = isos )
  IF ( isos /= 0 ) &
       CALL stop_memory ( 'force_control', 'f_nonbond', 3 * natoms )

! initialize ewalds
  IF ( first_time ) THEN
     SELECT CASE ( ewald_param % ewald_type )
     CASE ( 'EWALD_GAUSS' )
        CALL ewald_initialize ( dg, part, pnode, fc_global % group, &
             ewald_param, box, thermo, fc_global % scr, &
             ewald_grid = grid_ewald )
     CASE ( 'PME_GAUSS' )
        CALL ewald_initialize ( dg, part, pnode, fc_global % group, &
             ewald_param, box, thermo, fc_global % scr, &
             pme_small_grid = grid_s, pme_big_grid = grid_b )
     END SELECT
  END IF

! reinitialize the gspace for the new box
  SELECT CASE ( ensemble ( 1:3 )  )
  CASE ( 'NPT' )
     CALL get_hinv ( box )
     
     SELECT CASE ( ewald_param % ewald_type )
     CASE ( 'EWALD_GAUSS' )
        CALL pw_grid_change ( box, grid_s )
     CASE ( 'PME_GAUSS' )
        CALL stop_prg ( 'force_control',  &
             'constant pressure with PME not implemented' )
     END SELECT
  END SELECT

!
! first check with list_control to update neighbor lists
!
  CALL list_control ( pnode, part, box )
!
! initial force, energy and pressure tensor arrays
!
  DO i = 1, natoms
     part ( i ) % f ( 1 ) = 0.0_dbl
     part ( i ) % f ( 2 ) = 0.0_dbl
     part ( i ) % f ( 3 ) = 0.0_dbl
  END DO
  pv_bond = 0.0_dbl
  pv_bend = 0.0_dbl
  thermo % pot = 0.0_dbl
  thermo % pv = 0.0_dbl
!
! get real-space non-bonded forces:
!
  f_nonbond = zero
  CALL force_nonbond ( ewald_param,pnode,box,potparm, &
       pot_nonbond,f_nonbond,pv_nonbond )
!
! get g-space non-bonded forces:
!
  IF ( ewald_param % ewald_type /= 'NONE' ) THEN
     IF ( .NOT. ALLOCATED ( fg_coulomb )  ) &
          ALLOCATE ( fg_coulomb ( 3,nnodes ), STAT=isos )
     IF ( isos /= 0 ) &
          CALL stop_memory ( 'force_control', 'fg_coulomb', 3 * nnodes )
     
! compute g-space part of the ewald sum
     SELECT CASE ( ewald_param % ewald_type )
        
     CASE ( "EWALD_GAUSS" )
        CALL ewald_gaussian ( dg, ewald_param, fg_coulomb, vg_coulomb, &
             pv_g, pnode, box )
        CALL bonded_correct_gaussian ( ewald_param, molecule, &
             thermo % e_bonded, pv_bc )
        
     CASE ( "PME_GAUSS" )
        CALL pme_evaluate ( dg, fg_coulomb, vg_coulomb, pv_g, box, &
             grid_s, grid_b, ewald_param )
        CALL bonded_correct_gaussian ( ewald_param, molecule, &
             thermo % e_bonded, pv_bc )
        
     CASE DEFAULT
        CALL stop_prg ( "force_control", "illegal value of ewald_type:", &
             ewald_param % ewald_type )
     END SELECT
  END IF
  
!
! get intramolecular forces
!
  IF ( PRESENT ( debug )  ) THEN
     CALL force_intra_control ( molecule, pot_bond, pot_bend, &
          pv_bond, pv_bend, debug % f_bond,debug % f_bend )
  ELSE
     CALL force_intra_control ( molecule, pot_bond, pot_bend, &
          pv_bond, pv_bend )
  END IF
  
!
! add up all the potential energies
!
  IF ( ewald_param % ewald_type == 'NONE' ) THEN
     thermo % pot = pot_nonbond + pot_bond + pot_bend
#if defined ( __parallel )
     CALL mp_sum ( thermo % pot,fc_global % group )
#endif
  ELSE
     thermo % pot = pot_nonbond + pot_bond + pot_bend + thermo % e_bonded
     thermo % gspace = vg_coulomb
#if defined ( __parallel )
     CALL mp_sum ( thermo % pot,fc_global % group )
     CALL mp_sum ( thermo % e_bonded,fc_global % group )
#endif
! e_self and e_neut are already summed over all processors
! vg_coulomb is not calculated in parallel
     thermo % pot = thermo % pot + thermo % e_self + thermo % e_neut/box % deth
     thermo % pot = thermo % pot + vg_coulomb
  END IF

! add up all the forces
! nonbonded forces might be claculated for atoms not on this node
! ewald forces are strictly local -> sum only over pnode
  DO i = 1, natoms
     part ( i ) % f ( 1 ) = part ( i ) % f ( 1 ) + f_nonbond ( 1, i )
     part ( i ) % f ( 2 ) = part ( i ) % f ( 2 ) + f_nonbond ( 2, i )
     part ( i ) % f ( 3 ) = part ( i ) % f ( 3 ) + f_nonbond ( 3, i )
  END DO
  IF ( ewald_param % ewald_type /= 'NONE' ) THEN
     DO i = 1, nnodes
        pnode ( i ) % p % f ( 1 ) = pnode ( i ) % p % f ( 1 ) &
             + fg_coulomb ( 1, i )
        pnode ( i ) % p % f ( 2 ) = pnode ( i ) % p % f ( 2 ) &
             + fg_coulomb ( 2, i )
        pnode ( i ) % p % f ( 3 ) = pnode ( i ) % p % f ( 3 ) &
             + fg_coulomb ( 3, i )
     END DO
  END IF

! add up all the pressure tensors
  IF ( ewald_param % ewald_type == 'NONE' ) THEN
     thermo % pv = pv_nonbond + pv_bond + pv_bend
#if defined ( __parallel )
     CALL mp_sum ( thermo % pv, fc_global % group )
#endif
  ELSE
     ident = 0.0_dbl
     DO i = 1, 3
        ident ( i, i ) = 1.0_dbl
     END DO
     
     thermo % pv = pv_nonbond + pv_bond + pv_bend + pv_bc
#if defined ( __parallel )
     CALL mp_sum ( thermo % pv,fc_global % group )
#endif
     
     thermo % pv = thermo % pv + ident * thermo % e_neut / box % deth
     thermo % pv = thermo % pv + pv_g
  END IF

!
! if we are doing debugging, check if variables are present and assign
!
  IF ( PRESENT ( debug )  ) THEN
     debug % pot_bond = pot_bond
     debug % pot_bend = pot_bend
     debug % pot_nonbond = pot_nonbond
     debug % f_nonbond = f_nonbond
     debug % pv_nonbond = pv_nonbond
     debug % pv_bond = pv_bond
     debug % pv_bend = pv_bend
     IF ( ewald_param % ewald_type /= 'NONE' ) THEN
        debug % pot_g = vg_coulomb
        debug % f_g = fg_coulomb
        debug % pv_g = pv_g
     ELSE
        debug % pot_g = 0.0_dbl
        debug % f_g = 0.0_dbl
        debug % pv_g = 0.0_dbl
     END IF
  END IF
  
! deallocating all local variables
  isos = 0
  IF ( ALLOCATED ( fg_coulomb )  ) DEALLOCATE ( fg_coulomb, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'force_control', 'fg_coulomb' )
  
#if defined ( __parallel )
  DO i = 1, natoms
     CALL mp_sum ( part ( i ) % f, fc_global % group )
  END DO
#endif
  
  CALL timestop ( zero, handle )
  
END SUBROUTINE force_control

!******************************************************************************

END MODULE fist_force
