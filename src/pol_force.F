!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C ) 2000  CP2K developers group                                !
!-----------------------------------------------------------------------------!
!!****** cp2k/fist_force [1.0] *
!!
!!   NAME
!!     pol_force
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     gloria   
!!
!!   SOURCE
!******************************************************************************

MODULE pol_force

  USE atomic_kinds, ONLY : kind_info_type
  USE dg_types, ONLY : dg_type
  USE coefficient_types, ONLY : coeff_type
  USE empirical_parameters, ONLY : empirical_parameter_type
  USE ewald_parameters_types, ONLY : ewald_parameters_type
  USE fist_intra_force, ONLY : force_intra_control
  USE fist_nonbond_force, ONLY : force_nonbond, bonded_correct_gaussian
  USE global_types, ONLY : global_environment_type
  USE greens_fn, ONLY : pw_green_fn_init
  USE kinds, ONLY : dbl
  USE linklist_control, ONLY : list_control
  USE mathconstants, ONLY : pi, zero
  USE md, ONLY : thermodynamic_type
  USE molecule_types, ONLY : molecule_structure_type, particle_node_type, &
       linklist_neighbor, linklist_exclusion, linklist_bonds, &
       linklist_bends, linklist_torsions
  USE message_passing, ONLY : mp_sum
  USE pair_potential, ONLY : potentialparm_type
  USE pol_electrostatics, ONLY : electrostatics
  USE pol_overlap, ONLY : force_overlap
  USE particle_types, ONLY : particle_type
  USE pol_grids, ONLY : pol_grids_allocate, pol_grids_count, &
       pol_grids_map, pol_grids_initialize, assign_coefs_to_grids  
  USE pme, ONLY: pme_evaluate
  USE pw_grid_types, ONLY : pw_grid_type
  USE pw_grids, ONLY : pw_grid_change
  USE simulation_cell, ONLY : cell_type, get_hinv
  USE spme, ONLY : spme_evaluate
  USE termination, ONLY : stop_program, stop_memory
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
!!****** pol_force/force_control [1.0] *
!!
!!   NAME
!!     force_control
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE force_control ( molecule, pnode, part, box, box_ref, &
     drho_basis_info, rho0_basis_info, coef_pos, coef_vel, coef_force, &
     thermo, potparm, empparm, ewald_param, ensemble, fc_global, lldata,  debug)

  IMPLICIT NONE

! Arguments
  TYPE ( particle_type ), DIMENSION ( : ), INTENT ( INOUT ) :: part
  TYPE ( particle_node_type ), DIMENSION ( : ), INTENT ( INOUT ) :: pnode
  TYPE ( molecule_structure_type ), DIMENSION ( : ), INTENT ( IN ) :: molecule
  TYPE ( cell_type ), INTENT ( INOUT ) :: box
  TYPE ( cell_type ), INTENT ( IN ) :: box_ref
  TYPE ( kind_info_type ), POINTER, DIMENSION ( : ) :: drho_basis_info
  TYPE ( kind_info_type ), POINTER, DIMENSION ( : ) :: rho0_basis_info
  TYPE ( coeff_type ), INTENT (inout) :: coef_pos
  TYPE ( coeff_type ), INTENT (inout) :: coef_vel
  TYPE ( coeff_type ), INTENT (inout) :: coef_force
  TYPE ( thermodynamic_type ), INTENT ( INOUT ) :: thermo
  TYPE ( potentialparm_type ), DIMENSION ( :,: ), INTENT ( IN ) :: potparm
  TYPE ( empirical_parameter_type ), DIMENSION ( : ), INTENT ( IN ) :: empparm
  TYPE ( ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param
  CHARACTER ( LEN = * ), INTENT ( IN ) :: ensemble
  TYPE ( global_environment_type ), INTENT ( IN ) :: fc_global
  TYPE ( linklist_internal_data_type), DIMENSION (:), INTENT (INOUT) :: lldata
  TYPE ( debug_variables_type ), INTENT ( OUT ), OPTIONAL :: debug

! Locals
  INTEGER :: id, i, ii, natoms, nnodes, handle, isos, ncoefs
  INTEGER :: ng, igrid, ng_coef, ng_atom, nt_coef, nt_atom, iref
  REAL ( dbl ) :: pot_nonbond                                                
  REAL ( dbl ), DIMENSION ( :,: ), ALLOCATABLE, SAVE :: f_nonbond
  REAL ( dbl ), DIMENSION ( :,: ), ALLOCATABLE :: fc_part
  REAL ( dbl ), DIMENSION ( :,: ), ALLOCATABLE :: fo_part
  REAL ( dbl ), DIMENSION (:), ALLOCATABLE :: fc_coef
  REAL ( dbl ), DIMENSION (:), ALLOCATABLE :: fo_coef
  REAL ( dbl ), DIMENSION ( :,: ), ALLOCATABLE :: f_total
  REAL ( dbl ), DIMENSION ( 3,3 ) :: pv_nonbond
  TYPE ( pw_grid_type ),POINTER, DIMENSION (:), SAVE :: big_grid, small_grid
  TYPE ( dg_type ), POINTER, DIMENSION (:), SAVE :: dg_coef
  TYPE ( dg_type ), POINTER, DIMENSION (:), SAVE :: dg_atom
  LOGICAL :: first_time

!------------------------------------------------------------------------------

  CALL timeset ( 'FORCE','I',' ',handle )

  nnodes = SIZE ( pnode )
  natoms = SIZE ( part )
  ncoefs = size ( coef_pos % ao % cr )
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

  IF ( first_time ) THEN

     CALL pol_grids_count ( drho_basis_info, rho0_basis_info, ng_coef, ng_atom, &
                             nt_coef, nt_atom )

     ng = ng_coef + ng_atom
     
     CALL pol_grids_allocate (dg_coef, dg_atom, small_grid, big_grid, &
              nt_coef, nt_atom, ng) 
     
     igrid = 0
     
     CALL pol_grids_map (dg_coef, drho_basis_info, igrid)
     
     CALL pol_grids_map (dg_atom, rho0_basis_info, igrid)

     CALL assign_coefs_to_grids (coef_pos % ao, part, drho_basis_info, nt_coef)

     CALL pol_grids_initialize ( dg_coef, dg_atom, part, pnode,  coef_pos % ao, &
            drho_basis_info, rho0_basis_info, fc_global,  box, thermo, &
            small_grid, big_grid, iref  )    

     CALL pw_green_fn_init ( box % green, big_grid ( iref ) )
	   
  END IF

!
! first check with list_control to update neighbor lists
!
  CALL list_control ( lldata (1), pnode, part, box )
!
! check to update the overlap neighbor list
!
  CALL list_control ( lldata (2), pnode, part, box, drho_basis_info )
 
!
! initial force, energy and pressure tensor arrays
!
  DO i = 1, natoms
     part ( i ) % f ( 1 ) = 0.0_dbl
     part ( i ) % f ( 2 ) = 0.0_dbl
     part ( i ) % f ( 3 ) = 0.0_dbl
  END DO
  thermo % pot = 0.0_dbl
  thermo % pv = 0.0_dbl


 IF ( .NOT. ALLOCATED ( fc_part )  ) &
        ALLOCATE ( fc_part ( 3,nnodes ), STAT=isos )
 IF ( isos /= 0 ) &
     CALL stop_memory ( 'force_control', 'fc_part', 3 * nnodes )

 IF ( .NOT. ALLOCATED ( fc_coef )  ) &
        ALLOCATE ( fc_coef ( ncoefs ), STAT=isos )
 IF ( isos /= 0 ) &
     CALL stop_memory ( 'force_control', 'fc_coef',  ncoefs )

 IF ( .NOT. ALLOCATED ( fo_part )  ) &
         ALLOCATE ( fo_part ( 3,nnodes ), STAT=isos )
 IF ( isos /= 0 ) &
      CALL stop_memory ( 'force_control', 'fc_part', 3 * nnodes )
             
 IF ( .NOT. ALLOCATED ( fo_coef )  ) &
         ALLOCATE ( fo_coef ( ncoefs ), STAT=isos )
 IF ( isos /= 0 ) &
    CALL stop_memory ( 'force_control', 'fo_coef',  ncoefs )

! get overlap forces

     CALL force_overlap( coef_pos%ao, drho_basis_info, part, pnode, box, thermo, & 
                       empparm, fo_coef, fo_part, fc_global)


! get electrostatics forces   

   CALL electrostatics (dg_atom, dg_coef, thermo, part, coef_pos, box,  &
               small_grid, big_grid, ewald_param, iref, fc_coef, fc_part )
 
!
! get real-space non-bonded forces:
!
  f_nonbond = zero
  CALL force_nonbond ( ewald_param,pnode,box,potparm, &
       pot_nonbond,f_nonbond,pv_nonbond )

!
! add up all the potential energies
!

! add up all the forces
! nonbonded forces might be calculated for atoms not on this node
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
     DO i = 1, nnodes
        ii = pnode ( i ) % p % iatom
        f_total ( 1, ii ) = f_total ( 1, ii ) + fc_part ( 1, i )
        f_total ( 2, ii ) = f_total ( 2, ii ) + fc_part ( 2, i )
        f_total ( 3, ii ) = f_total ( 3, ii ) + fc_part ( 3, i )
     END DO


! deallocating all local variables
  isos = 0
  IF ( ALLOCATED ( fc_part )  ) DEALLOCATE ( fc_part, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'force_control', 'fc_part' )
  isos = 0
  IF ( ALLOCATED ( f_total )  ) DEALLOCATE ( f_total, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'force_control', 'f_total' )

  CALL timestop ( zero, handle )

END SUBROUTINE force_control

!!*****
!******************************************************************************

END MODULE pol_force

!******************************************************************************
