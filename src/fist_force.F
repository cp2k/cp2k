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
!!
!!   SOURCE
!******************************************************************************

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
       linklist_neighbor, linklist_exclusion, linklist_bonds, &
       linklist_bends, linklist_torsions
  USE message_passing, ONLY : mp_sum
  USE pair_potential, ONLY : potentialparm_type
  USE particle_types, ONLY : particle_type
  USE pme, ONLY: pme_evaluate
  USE pw_grid_types, ONLY : pw_grid_type
  USE pw_grids, ONLY : pw_grid_change
  USE simulation_cell, ONLY : cell_type, get_hinv
  USE stop_program, ONLY : stop_prg, stop_memory
  USE timings, ONLY : timeset, timestop
  USE linklist_types, ONLY : linklist_internal_data_type
  
  PRIVATE
  PUBLIC :: force_control, debug_variables_type
  
  TYPE debug_variables_type
     REAL ( dbl ) :: pot_nonbond, pot_g, pot_bond, pot_bend, pot_torsion
     REAL ( dbl ), DIMENSION ( :, : ), POINTER :: &
          f_nonbond, f_g, f_bond, f_bend, f_torsion
     REAL ( dbl ), DIMENSION ( 3, 3 ) :: pv_nonbond, pv_g, pv_bond, &
          pv_bend, pv_torsion
  END TYPE debug_variables_type
  
!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** fist_force/force_control [1.0] *
!!
!!   NAME
!!     force_control
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
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE force_control ( molecule, pnode, part, box, thermo, &
     potparm, ewald_param, ensemble, fc_global, lldata,  debug)

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
  TYPE ( linklist_internal_data_type), INTENT (INOUT) :: lldata
  TYPE ( debug_variables_type ), INTENT ( OUT ), OPTIONAL :: debug
  
! Locals
  INTEGER :: id, i, ii, natoms, nnodes, handle, isos
  REAL ( dbl ) :: pot_nonbond, pot_bond, pot_bend, pot_torsion, vg_coulomb
  REAL ( dbl ), DIMENSION ( :,: ), ALLOCATABLE, SAVE :: f_nonbond
  REAL ( dbl ), DIMENSION ( 3,3 ) :: pv_nonbond, pv_bond, pv_bend, pv_torsion
  REAL ( dbl ), DIMENSION ( :,: ), ALLOCATABLE :: fg_coulomb
  REAL ( dbl ), DIMENSION ( :,: ), ALLOCATABLE :: f_total
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
  IF ( .NOT. ALLOCATED ( f_nonbond )  ) THEN
    ALLOCATE ( f_nonbond ( 3,natoms ), STAT = isos )
    IF ( isos /= 0 ) &
       CALL stop_memory ( 'force_control', 'f_nonbond', 3 * natoms )
  ELSE IF ( SIZE ( f_nonbond ( 1, : ) ) < natoms ) THEN
    DEALLOCATE ( f_nonbond, STAT = isos )
    IF ( isos /= 0 ) CALL stop_memory ( 'force_control', 'f_nonbond' )
    ALLOCATE ( f_nonbond ( 3,natoms ), STAT = isos )
    IF ( isos /= 0 ) &
       CALL stop_memory ( 'force_control', 'f_nonbond', 3 * natoms )
  END IF

! initialize ewalds
  IF ( first_time ) THEN
     SELECT CASE ( ewald_param % ewald_type )
     CASE ( 'EWALD' )
        CALL ewald_initialize ( dg, part, pnode, fc_global, &
             ewald_param, box, thermo, ewald_grid = grid_ewald )
     CASE ( 'PME' )
        CALL ewald_initialize ( dg, part, pnode, fc_global, &
             ewald_param, box, thermo, &
             pme_small_grid = grid_s, pme_big_grid = grid_b )
     END SELECT
  END IF

! reinitialize the gspace for the new box
  SELECT CASE ( ensemble ( 1:3 )  )
  CASE ( 'NPT' )
     CALL get_hinv ( box )
     
     SELECT CASE ( ewald_param % ewald_type )
     CASE ( 'EWALD' )
        CALL pw_grid_change ( box, grid_ewald )
     CASE ( 'PME' )
        CALL stop_prg ( 'force_control',  &
             'constant pressure with PME not implemented' )
     END SELECT
  END SELECT

!
! first check with list_control to update neighbor lists
!
  CALL list_control ( lldata, pnode, part, box )
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
  pv_torsion = 0.0_dbl
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
        
     CASE ( "EWALD" )
        CALL ewald_gaussian ( dg, ewald_param, fg_coulomb, vg_coulomb, &
             pv_g, pnode, box )
        CALL bonded_correct_gaussian ( ewald_param, molecule, &
             thermo % e_bonded, pv_bc )
        
     CASE ( "PME" )
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
     CALL force_intra_control ( molecule, pot_bond, pot_bend, pot_torsion, &
          pv_bond, pv_bend, pv_torsion, &
          debug % f_bond, debug % f_bend, debug % f_torsion )
  ELSE
     CALL force_intra_control ( molecule, pot_bond, pot_bend, pot_torsion, &
          pv_bond, pv_bend, pv_torsion )
  END IF
  
!
! add up all the potential energies
!
  IF ( ewald_param % ewald_type == 'NONE' ) THEN
     thermo % pot = pot_nonbond + pot_bond + pot_bend + pot_torsion
#if defined ( __parallel )
     CALL mp_sum ( thermo % pot,fc_global % group )
#endif
  ELSE
     thermo % pot = pot_nonbond + pot_bond + pot_bend + pot_torsion + thermo % e_bonded
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
! We first sum the forces in f_nonbond, this allows for a more efficient
! global sum in the parallel code and in the end copy them back to part
  isos = 0
  IF ( .NOT. ALLOCATED ( f_total ) ) ALLOCATE ( f_total ( 3,natoms ), STAT=isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'force_control', 'f_total', 3 * natoms )
  DO i = 1, natoms
     f_total ( 1, i ) = part ( i ) % f ( 1 ) + f_nonbond ( 1, i )
     f_total ( 2, i ) = part ( i ) % f ( 2 ) + f_nonbond ( 2, i )
     f_total ( 3, i ) = part ( i ) % f ( 3 ) + f_nonbond ( 3, i )
  END DO
  IF ( ewald_param % ewald_type /= 'NONE' ) THEN
     DO i = 1, nnodes
        ii = pnode ( i ) % p % iatom
        f_total ( 1, ii ) = f_total ( 1, ii ) + fg_coulomb ( 1, i )
        f_total ( 2, ii ) = f_total ( 2, ii ) + fg_coulomb ( 2, i )
        f_total ( 3, ii ) = f_total ( 3, ii ) + fg_coulomb ( 3, i )
     END DO
  END IF

! add up all the pressure tensors
  IF ( ewald_param % ewald_type == 'NONE' ) THEN
     thermo % pv = pv_nonbond + pv_bond + pv_bend + pv_torsion
#if defined ( __parallel )
     CALL mp_sum ( thermo % pv, fc_global % group )
#endif
  ELSE
     ident = 0.0_dbl
     DO i = 1, 3
        ident ( i, i ) = 1.0_dbl
     END DO
     
     thermo % pv = pv_nonbond + pv_bond + pv_bend + pv_torsion + pv_bc
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
     debug % pot_torsion = pot_torsion
     debug % pot_nonbond = pot_nonbond
     debug % f_nonbond = f_nonbond
     debug % pv_nonbond = pv_nonbond
     debug % pv_bond = pv_bond
     debug % pv_bend = pv_bend
     debug % pv_torsion = pv_torsion
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
  
#if defined ( __parallel )
  CALL mp_sum ( f_total, fc_global % group )
#endif

  DO i = 1, natoms
    part ( i ) % f ( 1 ) = f_total ( 1, i )
    part ( i ) % f ( 2 ) = f_total ( 2, i )
    part ( i ) % f ( 3 ) = f_total ( 3, i )
  END DO
  
! deallocating all local variables
  isos = 0
  IF ( ALLOCATED ( fg_coulomb )  ) DEALLOCATE ( fg_coulomb, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'force_control', 'fg_coulomb' )
  isos = 0
  IF ( ALLOCATED ( f_total )  ) DEALLOCATE ( f_total, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'force_control', 'f_total' )

  CALL timestop ( zero, handle )
  
END SUBROUTINE force_control

!!*****
!******************************************************************************

END MODULE fist_force

!******************************************************************************
