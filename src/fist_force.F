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
!!
!!   SOURCE
!******************************************************************************

MODULE fist_force

!  USE dg_types, ONLY : dg_type
!  USE dgs, ONLY : dg_grid_change
!  USE ewald_parameters_types, ONLY : ewald_parameters_type
!  USE ewalds, ONLY : ewald_evaluate, ewald_initialize, ewald_setup
!  USE fist_environment_types, ONLY : replica_environment_type, &
!                                     fist_thermodynamic_type
!  USE fist_intra_force, ONLY : force_intra_control
!  USE fist_nonbond_force, ONLY : force_nonbond, bonded_correct_gaussian
!  USE global_types, ONLY : global_environment_type
!  USE greens_fn, ONLY : pw_green_fn_init
!  USE kinds, ONLY : dbl
!  USE linklist_control, ONLY : list_control
!  USE mathconstants, ONLY : pi, zero
!  USE molecule_types, ONLY : molecule_structure_type, particle_node_type
!  USE message_passing, ONLY : mp_sum
!  USE pair_potential_types, ONLY : pair_potential_type
!  USE particle_types, ONLY : particle_type
!  USE pme, ONLY: pme_evaluate
!  USE pw_grid_types, ONLY : pw_grid_type
!  USE pw_grids, ONLY : pw_grid_change
!  USE simulation_cell, ONLY : cell_type, get_hinv
!  USE spme, ONLY : spme_evaluate, spme_coeff_calculate
!  USE termination, ONLY : stop_program, stop_memory
!  USE timings, ONLY : timeset, timestop
!
!  PRIVATE
!  PUBLIC :: fist_force_control, debug_variables_type, get_ewald_grids
!
!  TYPE debug_variables_type
!     REAL ( dbl ) :: pot_nonbond, pot_g, pot_bond, pot_bend, pot_torsion
!     REAL ( dbl ), DIMENSION ( :, : ), POINTER :: &
!          f_nonbond, f_g, f_bond, f_bend, f_torsion
!     REAL ( dbl ), DIMENSION ( 3, 3 ) :: pv_nonbond, pv_g, pv_bond, &
!          pv_bend, pv_torsion
!  END TYPE debug_variables_type
!
!  REAL ( dbl ), DIMENSION ( :,: ), ALLOCATABLE :: f_nonbond
!  TYPE ( pw_grid_type ), TARGET :: grid_s, grid_b, grid_ewald, grid_spme
!  TYPE ( dg_type ), TARGET :: dg
!
!!******************************************************************************
!
! CONTAINS
!
!!******************************************************************************
!!!****** fist_force/force_env_types [1.1] *
!!!
!!!   NAME
!!!     force_env_types
!!!
!!!   FUNCTION
!!!     Calculates the total potential energy, total force, and the
!!!     total pressure tensor from the potentials
!!!
!!!   AUTHOR
!!!     CJM & JGH
!!!
!!!   MODIFICATION HISTORY
!!!     Harald Forbert (Dec-2000): Changes for multiple linked lists
!!!     cjm, 20-Feb-2001: box_ref used to initialize ewald.  Now
!!!                       have consistent restarts with npt and ewald
!!!     JGH (15-Mar-2001) : box_change replaces ensemble parameter
!!!                         Call ewald_setup if box has changed
!!!                         Consistent setup for PME and SPME
!!!
!!!   SOURCE
!!******************************************************************************
!
!SUBROUTINE fist_force_control ( rep_env, ewald_param, potparm, thermo, &
!                           box_change, globenv, debug )
!
!  IMPLICIT NONE
!
!! Arguments
!  TYPE ( replica_environment_type ), INTENT ( INOUT ):: rep_env
!  TYPE ( ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param
!  TYPE ( pair_potential_type ), INTENT ( IN )  :: potparm ( :, : )
!  TYPE ( fist_thermodynamic_type ), INTENT ( INOUT ) :: thermo
!  LOGICAL, INTENT ( IN ) :: box_change
!  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv
!  TYPE ( debug_variables_type ), INTENT ( OUT ), OPTIONAL :: debug
!
!! Locals
!  TYPE ( cell_type ), POINTER :: box
!  TYPE ( cell_type ), POINTER :: box_ref
!  TYPE ( molecule_structure_type ), POINTER :: molecule ( : )
!  TYPE ( particle_type ), POINTER :: part ( : )
!  TYPE ( particle_node_type ), POINTER :: pnode ( : )
!  INTEGER :: id, i, ii, natoms, nnodes, handle, isos, n
!  REAL ( dbl ) :: pot_nonbond, pot_bond, pot_bend, pot_torsion, vg_coulomb
!  REAL ( dbl ), DIMENSION ( 3,3 ) :: pv_nonbond, pv_bond, pv_bend, pv_torsion
!  REAL ( dbl ), DIMENSION ( :,: ), ALLOCATABLE :: fg_coulomb
!  REAL ( dbl ), DIMENSION ( :,: ), ALLOCATABLE :: f_total
!  REAL ( dbl ), DIMENSION ( 3,3 ) :: pv_g, ident
!  REAL ( dbl ), DIMENSION ( 3,3 ) :: pv_bc
!  LOGICAL :: first_time
!
!! associcating local pointers
!  molecule => rep_env % molecule
!  part => rep_env % part
!  pnode => rep_env % pnode
!  box => rep_env % box
!  box_ref => rep_env % box_ref
!
!!------------------------------------------------------------------------------
!
!  CALL timeset ( 'FORCE','I',' ',handle )
!
!  nnodes = SIZE ( pnode )
!  natoms = SIZE ( part )
!  isos = 0
!
!  first_time = .NOT. ALLOCATED ( f_nonbond )
!
!  IF ( .NOT. ALLOCATED ( f_nonbond )  ) THEN
!    ALLOCATE ( f_nonbond ( 3,natoms ), STAT = isos )
!    IF ( isos /= 0 ) &
!       CALL stop_memory ( 'force_env_types', 'f_nonbond', 3 * natoms )
!  ELSE IF ( SIZE ( f_nonbond, 2  ) < natoms ) THEN
!    DEALLOCATE ( f_nonbond, STAT = isos )
!    IF ( isos /= 0 ) CALL stop_memory ( 'force_env_types', 'f_nonbond' )
!    ALLOCATE ( f_nonbond ( 3,natoms ), STAT = isos )
!    IF ( isos /= 0 ) &
!       CALL stop_memory ( 'force_env_types', 'f_nonbond', 3 * natoms )
!    ! change in system, do a full reinitialisation
!    first_time = .TRUE.
!  END IF
!
!! initialize ewalds with the REFERENCE box ( for consistent NPT restarts )
!  IF ( first_time ) THEN
!
!     SELECT CASE ( ewald_param % ewald_type )
!     CASE ( 'EWALD' )
!        CALL ewald_initialize ( dg, part, pnode, globenv, &
!             ewald_param, box_ref, thermo % e_self, thermo % e_neut,  &
!             ewald_grid = grid_ewald )
!     CASE ( 'PME' )
!        CALL ewald_initialize ( dg, part, pnode, globenv, &
!             ewald_param, box_ref, thermo % e_self, thermo % e_neut, &
!             pme_small_grid = grid_s, pme_big_grid = grid_b )
!     CASE ( 'SPME' )
!        CALL ewald_initialize ( dg, part, pnode, globenv, &
!             ewald_param, box_ref, thermo % e_self, thermo % e_neut, &
!             spme_grid = grid_spme )
!        n = ewald_param%o_spline
!        ALLOCATE ( box % green % p3m_coeff ( -(n-1):n-1, 0:n-1 ), STAT = isos )
!        IF ( isos /= 0 ) CALL stop_memory ( 'spme','coeff', 2*n*n )
!        CALL spme_coeff_calculate ( n,  box % green % p3m_coeff )
!     END SELECT
!
!  END IF
!
!! reinitialize the gspace for the new box
!  IF ( box_change .OR. first_time ) THEN
!
!     CALL get_hinv ( box )
!
!     SELECT CASE ( ewald_param % ewald_type )
!     CASE ( 'EWALD' )
!        CALL pw_grid_change ( box, grid_ewald )
!        CALL ewald_setup ( grid_ewald, ewald_param, dg )
!     CASE ( 'PME' )
!        CALL pw_grid_change ( box, grid_b )
!        CALL pw_green_fn_init ( box % green, grid_b )
!        CALL dg_grid_change ( box, grid_b, grid_s )
!        CALL ewald_setup ( grid_s, ewald_param, dg )
!     CASE ( 'SPME' )
!        CALL pw_grid_change ( box, grid_spme )
!        CALL pw_green_fn_init ( box % green, grid_spme, p3m=ewald_param%o_spline, &
!             alpha = ewald_param%alpha )
!     END SELECT
!
!  END IF
!
!!
!! first check with list_control to update neighbor lists
!!
!  CALL list_control ( rep_env )
!
!!
!! initial force, energy and pressure tensor arrays
!!
!  DO i = 1, natoms
!     part ( i ) % f ( 1 ) = 0.0_dbl
!     part ( i ) % f ( 2 ) = 0.0_dbl
!     part ( i ) % f ( 3 ) = 0.0_dbl
!  END DO
!  pv_bond = 0.0_dbl
!  pv_bend = 0.0_dbl
!  pv_torsion = 0.0_dbl
!  thermo % pot = 0.0_dbl
!  thermo % virial % pv_virial = 0.0_dbl
!!
!! get real-space non-bonded forces:
!!
!  f_nonbond = zero
!  CALL force_nonbond ( ewald_param, part, pnode, box, potparm, &
!                       pot_nonbond, f_nonbond, pv_nonbond, &
!                       rep_env%ll_data(1)%nonbonded, &
!                       rep_env%ll_data(1)%r_last_update )
!!
!! get g-space non-bonded forces:
!!
!  IF ( ewald_param % ewald_type /= 'NONE' ) THEN
!
!     ! compute g-space part of the ewald sum
!     SELECT CASE ( ewald_param % ewald_type )
!
!     CASE ( "EWALD" )
!        ! Parallelisation over atoms --> allocate only nnodes
!        ALLOCATE ( fg_coulomb ( 3, nnodes ), STAT=isos )
!        IF ( isos /= 0 ) &
!          CALL stop_memory ( 'force_env_types', 'fg_coulomb', 3 * nnodes )
!        CALL ewald_evaluate ( dg, ewald_param, fg_coulomb, vg_coulomb, &
!             pv_g, pnode, box )
!        
!     CASE ( "PME" )
!        ! Parallelisation over grids --> allocate all atoms
!        ALLOCATE ( fg_coulomb ( 3, natoms ), STAT=isos )
!        IF ( isos /= 0 ) &
!          CALL stop_memory ( 'force_env_types', 'fg_coulomb', 3 * nnodes )
!        CALL pme_evaluate ( dg, part, box, grid_s, grid_b, ewald_param, &
!                            vg_coulomb, fg_coulomb, pv_g )
!        
!     CASE ( "SPME" )
!        ! Parallelisation over grids --> allocate all atoms
!        ALLOCATE ( fg_coulomb ( 3, natoms ), STAT=isos )
!        IF ( isos /= 0 ) &
!          CALL stop_memory ( 'force_env_types', 'fg_coulomb', 3 * nnodes )
!        CALL spme_evaluate ( dg, part, fg_coulomb, vg_coulomb, pv_g, box, &
!             grid_spme, ewald_param )
!        
!     CASE DEFAULT
!        CALL stop_program ( "force_env_types", "illegal value of ewald_type: "//&
!                            ewald_param % ewald_type )
!     END SELECT
!
!     CALL bonded_correct_gaussian ( ewald_param, molecule, &
!                                    thermo % e_bonded, pv_bc )
!
!  END IF
!
!!
!! get intramolecular forces
!!
!  IF ( PRESENT ( debug )  ) THEN
!     CALL force_intra_control ( molecule, pot_bond, pot_bend, pot_torsion, &
!          pv_bond, pv_bend, pv_torsion, &
!          debug % f_bond, debug % f_bend, debug % f_torsion )
!  ELSE
!     CALL force_intra_control ( molecule, pot_bond, pot_bend, pot_torsion, &
!          pv_bond, pv_bend, pv_torsion )
!  END IF
!
!!
!! add up all the potential energies
!!
!  IF ( ewald_param % ewald_type == 'NONE' ) THEN
!     thermo % pot = pot_nonbond + pot_bond + pot_bend + pot_torsion
!#if defined ( __parallel )
!     CALL mp_sum ( thermo % pot,globenv % group )
!#endif
!  ELSE
!     thermo % pot = pot_nonbond + pot_bond + pot_bend + pot_torsion + thermo % e_bonded
!     thermo % e_gspace = vg_coulomb
!#if defined ( __parallel )
!     CALL mp_sum ( thermo % pot,globenv % group )
!     CALL mp_sum ( thermo % e_bonded,globenv % group )
!#endif
!! e_self and e_neut are already summed over all processors
!! vg_coulomb is not calculated in parallel
!     thermo % pot = thermo % pot + thermo % e_self + thermo % e_neut/box % deth
!     thermo % pot = thermo % pot + vg_coulomb
!  END IF
!
!! add up all the forces
!! nonbonded forces might be claculated for atoms not on this node
!! ewald forces are strictly local -> sum only over pnode
!! We first sum the forces in f_nonbond, this allows for a more efficient
!! global sum in the parallel code and in the end copy them back to part
!  ALLOCATE ( f_total ( 3,natoms ), STAT=isos )
!  IF ( isos /= 0 ) CALL stop_memory ( 'force_env_types', 'f_total', 3 * natoms )
!  DO i = 1, natoms
!     f_total ( 1, i ) = part ( i ) % f ( 1 ) + f_nonbond ( 1, i )
!     f_total ( 2, i ) = part ( i ) % f ( 2 ) + f_nonbond ( 2, i )
!     f_total ( 3, i ) = part ( i ) % f ( 3 ) + f_nonbond ( 3, i )
!  END DO
!  IF ( ewald_param % ewald_type /= 'NONE' ) THEN
!     SELECT CASE ( ewald_param % ewald_type )
!     CASE ( "EWALD" )
!        DO i = 1, nnodes
!           ii = pnode ( i ) % p % iatom
!           f_total ( 1, ii ) = f_total ( 1, ii ) + fg_coulomb ( 1, i )
!           f_total ( 2, ii ) = f_total ( 2, ii ) + fg_coulomb ( 2, i )
!           f_total ( 3, ii ) = f_total ( 3, ii ) + fg_coulomb ( 3, i )
!        END DO
!     CASE ( "PME" )
!        f_total ( :, : ) = f_total ( :, : ) + fg_coulomb ( :, : )
!     CASE ( "SPME" )
!        f_total ( :, : ) = f_total ( :, : ) + fg_coulomb ( :, : )
!     CASE DEFAULT
!        CALL stop_program ( "force_env_types", "illegal value of ewald_type: "//&
!                            ewald_param % ewald_type )
!     END SELECT
!  END IF
!
!! add up all the pressure tensors
!  IF ( ewald_param % ewald_type == 'NONE' ) THEN
!     thermo % virial % pv_virial = pv_nonbond + pv_bond + pv_bend + pv_torsion
!#if defined ( __parallel )
!     CALL mp_sum ( thermo % virial % pv_virial, globenv % group )
!#endif
!  ELSE
!     ident = 0.0_dbl
!     DO i = 1, 3
!        ident ( i, i ) = 1.0_dbl
!     END DO
!
!     thermo % virial % pv_virial = pv_nonbond + pv_bond + pv_bend + pv_torsion + pv_bc
!#if defined ( __parallel )
!     CALL mp_sum ( thermo % virial % pv_virial, globenv % group )
!#endif
!
!     thermo % virial % pv_virial = thermo % virial % pv_virial + &
!                                      ident * thermo % e_neut / box % deth
!     thermo % virial % pv_virial = thermo % virial % pv_virial + pv_g
!  END IF
!
!!
!! if we are doing debugging, check if variables are present and assign
!!
!  IF ( PRESENT ( debug )  ) THEN
!     debug % pot_bond = pot_bond
!     debug % pot_bend = pot_bend
!     debug % pot_torsion = pot_torsion
!     debug % pot_nonbond = pot_nonbond
!     debug % f_nonbond = f_nonbond
!     debug % pv_nonbond = pv_nonbond
!     debug % pv_bond = pv_bond
!     debug % pv_bend = pv_bend
!     debug % pv_torsion = pv_torsion
!     IF ( ewald_param % ewald_type /= 'NONE' ) THEN
!        debug % pot_g = vg_coulomb
!        debug % f_g = fg_coulomb
!        debug % pv_g = pv_g
!     ELSE
!        debug % pot_g = 0.0_dbl
!        debug % f_g = 0.0_dbl
!        debug % pv_g = 0.0_dbl
!     END IF
!  END IF
!
!#if defined ( __parallel )
!  CALL mp_sum ( f_total, globenv % group )
!#endif
!
!  DO i = 1, natoms
!    part ( i ) % f ( 1 ) = f_total ( 1, i )
!    part ( i ) % f ( 2 ) = f_total ( 2, i )
!    part ( i ) % f ( 3 ) = f_total ( 3, i )
!  END DO
!
!! deallocating all local variables
!  IF ( ALLOCATED ( fg_coulomb ) ) THEN
!    DEALLOCATE ( fg_coulomb, STAT = isos )
!    IF ( isos /= 0 ) CALL stop_memory ( 'force_env_types', 'fg_coulomb' )
!  END IF
!  IF ( ALLOCATED ( f_total ) ) THEN
!    DEALLOCATE ( f_total, STAT = isos )
!    IF ( isos /= 0 ) CALL stop_memory ( 'force_env_types', 'f_total' )
!  END IF
!
!  CALL timestop ( zero, handle )
!
!END SUBROUTINE fist_force_control
!
!!******************************************************************************

!SUBROUTINE get_ewald_grids ( pw_big, pw_small, pw_ewald, pw_spme, dblgrd )
!  IMPLICIT NONE
!  TYPE ( pw_grid_type ), POINTER, OPTIONAL :: pw_big, pw_small, pw_ewald, pw_spme 
!  TYPE ( dg_type ), POINTER, OPTIONAL :: dblgrd

!  IF ( PRESENT ( pw_big ) ) pw_big => grid_b
!  IF ( PRESENT ( pw_small ) ) pw_small => grid_s
!  IF ( PRESENT ( pw_ewald ) ) pw_ewald => grid_ewald
!  IF ( PRESENT ( pw_spme ) ) pw_spme => grid_spme
!  IF ( PRESENT ( dblgrd ) ) dblgrd => dg
!
!END SUBROUTINE get_ewald_grids

!******************************************************************************

END MODULE fist_force
!
!******************************************************************************
