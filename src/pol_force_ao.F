!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C ) 2000  CP2K developers group                                !
!-----------------------------------------------------------------------------!
!!****** cp2k/pol_force_ao [1.0] *
!!
!!   NAME
!!     pol_force_ao
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     gloria   
!!
!!   SOURCE
!******************************************************************************

MODULE pol_force_ao

  USE atomic_kinds, ONLY : kind_info_type
  USE dg_types, ONLY : dg_type
  USE coefficient_types, ONLY : coeff_type
  USE empirical_parameters, ONLY : empirical_parameter_type
  USE ewald_parameters_types, ONLY : ewald_parameters_type
  USE ewalds, ONLY : ewald_evaluate, ewald_setup
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
  USE particle_types, ONLY : particle_type
  USE pme, ONLY : pme_evaluate
  USE pol_electrostatics_ao, ONLY : electrostatics
  USE pol_overlap_ao, ONLY : force_overlap
  USE pol_thk_ao, ONLY : force_thk
  USE pol_grids, ONLY : pol_grids_allocate, pol_grids_count, &
       pol_grids_map, pol_grids_initialize, assign_coefs_to_grids, &
       pol_grids_ref_grid
  USE pw_grid_types, ONLY : pw_grid_type
  USE pw_grids, ONLY : pw_grid_change
  USE simulation_cell, ONLY : cell_type, get_hinv
  USE spme, ONLY : spme_evaluate
  USE termination, ONLY : stop_program, stop_memory
  USE timings, ONLY : timeset, timestop
  USE linklist_types, ONLY : linklist_internal_data_type

  PRIVATE
  PUBLIC :: force_env_types, debug_variables_type, get_pol_grids

  TYPE debug_variables_type
     REAL ( dbl ) :: pot_electrostatics, pot_overlap
     REAL ( dbl ), DIMENSION ( :, : ), POINTER :: f_ovlp_part, f_elec_part
     REAL ( dbl ), DIMENSION ( : ), POINTER :: f_ovlp_coef, f_elec_coef
     REAL ( dbl ), DIMENSION ( 3, 3 ) :: pv_electrostatic, pv_overlap
  END TYPE debug_variables_type

! globals to module
  TYPE ( pw_grid_type ),POINTER, DIMENSION ( : ) :: big_grid, small_grid
  TYPE ( pw_grid_type ),POINTER  :: small_grid_pme, big_grid_pme
  TYPE ( pw_grid_type ),POINTER  :: grid_ewald
  TYPE ( dg_type ), POINTER, DIMENSION ( : ) :: dg_coef
  TYPE ( dg_type ), POINTER, DIMENSION ( : ) :: dg_atom
  TYPE ( dg_type ), POINTER :: dg_ewald
  REAL ( dbl ), DIMENSION ( :, : ), ALLOCATABLE :: f_nonbond
  INTEGER :: iref

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** pol_force_ao/force_env_types [1.0] *
!!
!!   NAME
!!     force_env_types
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE force_env_types ( molecule, pnode, part, box, box_ref, &
     drho_basis_info, rho0_basis_info, coef_pos, coef_vel, coef_force, &
     thermo, potparm, empparm, ewald_param, box_change, fc_global, lldata,&
     debug)

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
  LOGICAL, INTENT ( IN ) :: box_change
  TYPE ( global_environment_type ), INTENT ( IN ) :: fc_global
  TYPE ( linklist_internal_data_type), DIMENSION (:), INTENT (INOUT) :: lldata
  TYPE ( debug_variables_type ), INTENT ( OUT ), OPTIONAL :: debug

! Locals
  INTEGER :: id, i, ii, natoms, nnodes, handle, isos, ncoefs
  INTEGER :: ng, igrid, ng_coef, ng_atom, nt_coef, nt_atom, npts_max ( 3 )
  REAL ( dbl ) :: pot_nonbond, vg_coulomb                                                
  REAL ( dbl ) :: pot_overlap, pot_electrostatics 
  REAL ( dbl ) :: pot_thk                                                
  REAL ( dbl ) :: e0, e1, e2
  REAL ( dbl ) :: e_ovl_ii, e_ovl_ij
  REAL ( dbl ), DIMENSION ( :,: ), ALLOCATABLE :: fe_part
  REAL ( dbl ), DIMENSION ( :,: ), ALLOCATABLE :: fo_part
  REAL ( dbl ), DIMENSION (:), ALLOCATABLE :: fe_coef
  REAL ( dbl ), DIMENSION (:), ALLOCATABLE :: fo_coef
  REAL ( dbl ), DIMENSION (:), ALLOCATABLE :: fthk
  REAL ( dbl ), DIMENSION ( :,: ), ALLOCATABLE :: f_total
  REAL ( dbl ), DIMENSION ( :, : ), ALLOCATABLE :: fg_coulomb
  REAL ( dbl ), DIMENSION ( 3,3 ) :: pv_nonbond, pv_g, pv_bc, ident
  LOGICAL :: first_time
  REAL ( dbl ) :: max_cut, cutoff_radius_max

!------------------------------------------------------------------------------

  CALL timeset ( 'FORCE','I',' ',handle )

  nnodes = SIZE ( pnode )
  natoms = SIZE ( part )
  ncoefs = size ( coef_pos % ao % cr )
  isos = 0
  pot_overlap = 0._dbl
  pot_electrostatics = 0._dbl
  e0 = 0._dbl
  e1 = 0._dbl
  e2 = 0._dbl
  e_ovl_ii = 0._dbl
  e_ovl_ij = 0._dbl
  pot_nonbond = 0._dbl
  pot_thk = 0._dbl
  vg_coulomb = 0._dbl

  first_time = .NOT. ALLOCATED ( f_nonbond )

  IF ( .NOT. ALLOCATED ( f_nonbond )  ) THEN
    ALLOCATE ( f_nonbond ( 3,natoms ), STAT = isos )
    IF ( isos /= 0 ) &
       CALL stop_memory ( 'pol_force_control_ao', 'f_nonbond', 3 * natoms )
  ELSE IF ( SIZE ( f_nonbond ,2 ) < natoms ) THEN
    DEALLOCATE ( f_nonbond, STAT = isos )
    IF ( isos /= 0 ) CALL stop_memory ( 'pol_force_control_ao', 'f_nonbond' )
    ALLOCATE ( f_nonbond ( 3, natoms ), STAT = isos )
    IF ( isos /= 0 ) &
       CALL stop_memory ( 'pol_force_control_ao', 'f_nonbond', 3 * natoms )
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

     CALL pol_grids_ref_grid ( dg_coef, dg_atom, drho_basis_info,  &
           rho0_basis_info, box, small_grid, big_grid, iref, npts_max, &
           cutoff_radius_max, max_cut  )    

     SELECT CASE ( ewald_param % ewald_type )
     CASE ( "NONE" )
       CALL pol_grids_initialize ( ewald_param, dg_coef, dg_atom, part, pnode, &
                thermo, coef_pos % ao, drho_basis_info, rho0_basis_info,       & 
                fc_global, box, small_grid, big_grid, iref, npts_max,          &
                cutoff_radius_max )
     CASE ( "EWALD" )
       CALL pol_grids_initialize ( ewald_param, dg_coef, dg_atom, part, pnode, &
                thermo, coef_pos % ao, drho_basis_info, rho0_basis_info,       & 
                fc_global, box, small_grid, big_grid, iref, npts_max,          &
                cutoff_radius_max,                                             & 
                ewald_grid = grid_ewald,                                       &
                dg_ewald = dg_ewald )                                         
     CASE ( "PME" )
       CALL pol_grids_initialize ( ewald_param, dg_coef, dg_atom, part, pnode, &
                thermo, coef_pos % ao, drho_basis_info, rho0_basis_info,       & 
                fc_global, box, small_grid, big_grid, iref, npts_max,          &
                cutoff_radius_max,                                             & 
                pme_big_grid = big_grid_pme,                                   &
                dg_ewald = dg_ewald )                                        
     CASE ( "SPME" )
        CALL stop_program ( "pol_force_control", "SPME cannot be used with Pol" )
     END SELECT

     SELECT CASE ( ewald_param % ewald_type )
     CASE ( 'EWALD' )
        NULLIFY ( dg_ewald )
        ALLOCATE ( dg_ewald, STAT=isos )
        IF ( isos /= 0 ) CALL stop_memory ( 'pol_force_control_pw', 'dg_ewald' )
        CALL ewald_setup ( grid_ewald, ewald_param, dg_ewald )
     CASE ( 'PME' )
        IF ( big_grid ( iref ) % cutoff < big_grid_pme % cutoff ) THEN
          WRITE ( *, * ) "REF CUTOFF = ",  big_grid ( iref ) % cutoff
          WRITE ( *, * ) "PME CUTOFF = ", big_grid_pme % cutoff
          CALL stop_program ( 'pol_force_control_ao', 'REF cutoff too small' )
        ENDIF
! Note Bene: dg_ewald is ALLOCATED on pol_grids_initialize. 
        CALL ewald_setup ( small_grid_pme, ewald_param, dg_ewald )
     CASE ( 'SPME' )
        CALL stop_program ( "pol_force_control_ao", "SPME cannot be used with Pol" )
     END SELECT

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
  coef_force % ao % cr (:) = 0.0_dbl
  thermo % pot = 0.0_dbl
  thermo % e_class = 0.0_dbl
  thermo % e_electro = 0.0_dbl
  thermo % e_overlap = 0.0_dbl
  thermo % e_thk = 0.0_dbl
  thermo % pv = 0.0_dbl

  IF ( .NOT. ALLOCATED ( fe_part )  ) &
         ALLOCATE ( fe_part ( 3, nnodes ), STAT=isos )
  IF ( isos /= 0 ) &
      CALL stop_memory ( 'pol_force_control_ao', 'fe_part', 3 * nnodes )

  IF ( .NOT. ALLOCATED ( fo_part )  ) &
         ALLOCATE ( fo_part ( 3, nnodes ), STAT=isos )
  IF ( isos /= 0 ) &
    CALL stop_memory ( 'pol_force_control_ao', 'fo_part',  3 * nnodes )

  IF ( .NOT. ALLOCATED ( fe_coef )  ) &
         ALLOCATE ( fe_coef ( ncoefs ), STAT=isos )
  IF ( isos /= 0 ) &
    CALL stop_memory ( 'pol_force_control_ao', 'fe_coef',  ncoefs )

  IF ( .NOT. ALLOCATED ( fo_coef )  ) &
         ALLOCATE ( fo_coef ( ncoefs ), STAT=isos )
  IF ( isos /= 0 ) &
    CALL stop_memory ( 'pol_force_control_ao', 'fo_coef',  ncoefs )

  IF ( .NOT. ALLOCATED ( fthk )  ) &
         ALLOCATE ( fthk ( ncoefs ), STAT=isos )
  IF ( isos /= 0 ) &
    CALL stop_memory ( 'pol_force_control_ao', 'fthk',  ncoefs )

! get overlap forces

   CALL force_overlap( coef_pos % ao, drho_basis_info, part, pnode, box, & 
                       empparm, pot_overlap, fo_coef, fo_part, e_ovl_ii, e_ovl_ij )

! get first order hk contribution

   CALL force_thk( coef_pos % ao, empparm, fthk, pot_thk)

! get electrostatics forces   

   CALL electrostatics (dg_atom, dg_coef, part, coef_pos, box,  &
                        small_grid, big_grid, ewald_param, iref, &
                        pot_electrostatics, e0, e1, e2, fe_coef, fe_part )
 
!
! get real-space non-bonded forces:
!
  f_nonbond = zero
  IF ( fc_global % program_name == 'POL' ) &
  CALL force_nonbond ( ewald_param,pnode,box,potparm, &
                        pot_nonbond,f_nonbond,pv_nonbond )

! compute g-space part of the ewald sum
  SELECT CASE ( ewald_param % ewald_type )

  CASE ( "NONE" )
     ! nothing to be done
  CASE ( "EWALD" )
     ! Parallelisation over atoms --> allocate only nnodes
     ALLOCATE ( fg_coulomb ( 3, nnodes ), STAT=isos )
     IF ( isos /= 0 ) &
       CALL stop_memory ( 'pol_force_control_ao', 'fg_coulomb', 3 * nnodes )
     CALL ewald_evaluate ( dg_ewald, ewald_param, fg_coulomb, vg_coulomb, &
                           pv_g, pnode, box )
     CALL bonded_correct_gaussian ( ewald_param, molecule, &
          thermo % e_bonded, pv_bc )
        
  CASE ( "PME" )
     ! Parallelisation over grids --> allocate all atoms
     ALLOCATE ( fg_coulomb ( 3, natoms ), STAT=isos )
     IF ( isos /= 0 ) &
       CALL stop_memory ( 'pol_force_control_ao', 'fg_coulomb', 3 * nnodes )
     CALL pme_evaluate ( dg_ewald, part, box, small_grid_pme, big_grid_pme, &
                         ewald_param, vg_coulomb, fg_coulomb, pv_g )
     CALL bonded_correct_gaussian ( ewald_param, molecule, &
          thermo % e_bonded, pv_bc )
  CASE ( "SPME" )
     CALL stop_program ( "force_env_types", "SPME cannot be used with Pol" )
  CASE DEFAULT
     CALL stop_program ( "force_env_types", "illegal value of ewald_type: "//&
                         ewald_param % ewald_type )
  END SELECT


!
! add up all the potential energies
!

! add up all the forces:
   DO i = 1, ncoefs
     coef_force % ao % cr ( i )  = fo_coef ( i ) + fthk ( i )  + fe_coef ( i )  
   END DO

! Nonbonded forces might be calculated for atoms not on this node
! electrostatic forces are strictly local -> sum only over pnode
! We first sum the forces in f_nonbond, this allows for a more efficient
! global sum in the parallel code and in the end copy them back to part

  isos = 0
  IF ( .NOT. ALLOCATED ( f_total ) ) ALLOCATE ( f_total ( 3,natoms ), STAT=isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'pol_force_control_ao', 'f_total', 3 * natoms )

  DO i = 1, natoms
     f_total ( 1, i ) = part ( i ) % f ( 1 ) + f_nonbond ( 1, i )
     f_total ( 2, i ) = part ( i ) % f ( 2 ) + f_nonbond ( 2, i )
     f_total ( 3, i ) = part ( i ) % f ( 3 ) + f_nonbond ( 3, i )
  END DO

  SELECT CASE ( ewald_param % ewald_type )
  CASE ( "NONE" )
     ! nothing to be done
  CASE ( "EWALD" )
     DO i = 1, nnodes
        ii = pnode ( i ) % p % iatom
        f_total ( 1, ii ) = f_total ( 1, ii ) + fg_coulomb ( 1, i )
        f_total ( 2, ii ) = f_total ( 2, ii ) + fg_coulomb ( 2, i )
        f_total ( 3, ii ) = f_total ( 3, ii ) + fg_coulomb ( 3, i )
     END DO
  CASE ( "PME" )
     f_total ( :, : ) = f_total ( :, : ) + fg_coulomb ( :, : )
  CASE ( "SPME" )
     f_total ( :, : ) = f_total ( :, : ) + fg_coulomb ( :, : )
  CASE DEFAULT
     CALL stop_program ( "pol_force_ao", "illegal value of ewald_type: "//&
                         ewald_param % ewald_type )
  END SELECT

  DO i = 1, nnodes
     ii = pnode ( i ) % p % iatom
     f_total ( 1, ii ) = f_total ( 1, ii ) + fe_part ( 1, i ) + fo_part ( 1, i )
     f_total ( 2, ii ) = f_total ( 2, ii ) + fe_part ( 2, i ) + fo_part ( 2, i )
     f_total ( 3, ii ) = f_total ( 3, ii ) + fe_part ( 3, i ) + fo_part ( 3, i )
  END DO

  IF ( PRESENT ( debug )  ) THEN
     debug % pot_electrostatics = pot_electrostatics
     debug % pot_overlap  = pot_overlap
     debug % f_ovlp_part = fo_part
     debug % f_ovlp_coef = fo_coef
     debug % f_elec_part = fe_part
     debug % f_elec_coef = fe_coef
  END IF

  thermo % e_electro = pot_electrostatics
  thermo % e_electro_0 = e0
  thermo % e_electro_1 = e1
  thermo % e_electro_2 = e2
  thermo % e_overlap = pot_overlap
  thermo % e_overlap_ii = e_ovl_ii
  thermo % e_overlap_ij = e_ovl_ij
  thermo % e_thk = pot_thk
  thermo % e_pol = pot_electrostatics + pot_overlap + pot_thk - e0

  SELECT CASE ( ewald_param % ewald_type )
  CASE ( "NONE" )
     thermo % e_class = pot_nonbond 
     thermo % pot = pot_nonbond + pot_electrostatics + pot_overlap + pot_thk
#if defined ( __parallel )
     CALL mp_sum ( thermo % pot, fc_global % group )
#endif
  CASE DEFAULT
     thermo % pot = pot_nonbond + pot_electrostatics + pot_overlap + pot_thk  &
                  + thermo % e_bonded
     thermo % gspace = vg_coulomb
#if defined ( __parallel )
     CALL mp_sum ( thermo % pot, fc_global % group )
     CALL mp_sum ( thermo % e_bonded, fc_global % group )
#endif
! e_self and e_neut are already summed over all processors
! vg_coulomb is not calculated in parallel
     thermo % pot = thermo % pot + thermo % e_self +  &
                    thermo % e_neut/box % deth + vg_coulomb
     thermo % e_class = pot_nonbond + vg_coulomb + thermo % e_self +  &
                        thermo % e_neut / box % deth 
  END SELECT

#if defined ( __parallel )
  CALL mp_sum ( f_total, fc_global % group )
#endif

  DO i = 1, natoms
    part ( i ) % f ( 1 ) = f_total ( 1, i )
    part ( i ) % f ( 2 ) = f_total ( 2, i )
    part ( i ) % f ( 3 ) = f_total ( 3, i )
  END DO


! deallocating all local variables
  IF ( ALLOCATED ( f_total )  ) DEALLOCATE ( f_total, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'pol_force_control_ao', 'f_total' )
  IF ( ALLOCATED ( fg_coulomb ) ) DEALLOCATE ( fg_coulomb, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'pol_force_control_ao', 'fg_coulomb' )

  CALL timestop ( zero, handle )

END SUBROUTINE force_env_types

!******************************************************************************
!!****** pol_force_ao/get_pol_grids [1.0] *
!!
!!   NAME
!!     get_pol_grids
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE get_pol_grids ( pw_big, pw_small, dg1, dg2, ii )
  IMPLICIT NONE
  TYPE ( pw_grid_type ), POINTER, DIMENSION ( : ), OPTIONAL :: pw_big
  TYPE ( pw_grid_type ), POINTER, DIMENSION ( : ), OPTIONAL :: pw_small
  TYPE ( dg_type ), POINTER, DIMENSION ( : ), OPTIONAL :: dg1
  TYPE ( dg_type ), POINTER, DIMENSION ( : ), OPTIONAL :: dg2
  INTEGER, INTENT ( OUT ), OPTIONAL :: ii
  
  IF ( PRESENT ( pw_big ) ) pw_big  => big_grid 
  IF ( PRESENT ( pw_small ) ) pw_small => small_grid 
  IF ( PRESENT ( dg1 ) ) dg1 => dg_atom 
  IF ( PRESENT ( dg2 ) ) dg2 => dg_coef 
  IF ( PRESENT ( ii ) ) ii = iref

END SUBROUTINE get_pol_grids
!!*****
!******************************************************************************

END MODULE pol_force_ao

!******************************************************************************
