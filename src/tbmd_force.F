!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C ) 2000  CP2K developers group                                !
!-----------------------------------------------------------------------------!
!!****** cp2k/tbmd_force [1.0] *
!!
!!   NAME
!!     tbmd_force
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE tbmd_force

  USE dg_types,                        ONLY: dg_type
  USE dgs,                             ONLY: dg_grid_change
  USE ewald_parameters_types,          ONLY: ewald_parameters_type
  USE ewalds,                          ONLY: ewald_evaluate,&
                                             ewald_initialize,&
                                             ewald_setup
  USE fist_nonbond_force,              ONLY: force_nonbond
  USE global_types,                    ONLY: global_environment_type
  USE greens_fn,                       ONLY: pw_green_fn_init
  USE kinds,                           ONLY: dbl
  USE linklist_control,                ONLY: list_control
  USE linklist_types,                  ONLY: linklist_internal_data_type
  USE mathconstants,                   ONLY: zero
  USE md,                              ONLY: thermodynamic_type
  USE message_passing,                 ONLY: mp_sum
  USE molecule_types,                  ONLY: molecule_structure_type,&
                                             particle_node_type
  USE pair_potential,                  ONLY: potentialparm_type
  USE particle_types,                  ONLY: particle_type
  USE pme,                             ONLY: pme_evaluate
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_grids,                        ONLY: pw_grid_change
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_hinv
  USE spme,                            ONLY: spme_evaluate
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  PRIVATE
  PUBLIC :: force_env_types, debug_variables_type

  TYPE debug_variables_type
     REAL ( dbl ) :: pot_nonbond, pot_g
     REAL ( dbl ), DIMENSION ( :, : ), POINTER :: f_nonbond, f_g
     REAL ( dbl ), DIMENSION ( 3, 3 ) :: pv_nonbond, pv_g
  END TYPE debug_variables_type

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** tbmd_force/force_env_types [1.0] *
!!
!!   NAME
!!     force_env_types
!!
!!   FUNCTION
!!     Calculates the total potential energy, total force, and the
!!     total pressure tensor from the potentials
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     Harald Forbert (Dec-2000): Changes for multiple linked lists
!!     CJM, FEB-20-2001: box_ref used to initialize ewalds
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE force_env_types ( molecule, pnode, part, box, box_ref, thermo, &
     potparm, ewald_param, box_change, fc_global, lldata, debug )


    TYPE(molecule_structure_type), &
      DIMENSION(:), INTENT(IN)               :: molecule
    TYPE(particle_node_type), DIMENSION(:), &
      INTENT(INOUT)                          :: pnode
    TYPE(particle_type), DIMENSION(:), &
      INTENT(INOUT)                          :: part
    TYPE(cell_type), INTENT(INOUT)           :: box, box_ref
    TYPE(thermodynamic_type), INTENT(INOUT)  :: thermo
    TYPE(potentialparm_type), &
      DIMENSION(:, :), INTENT(IN)            :: potparm
    TYPE(ewald_parameters_type), &
      INTENT(INOUT)                          :: ewald_param
    LOGICAL, INTENT(IN)                      :: box_change
    TYPE(global_environment_type), &
      INTENT(IN)                             :: fc_global
    TYPE(linklist_internal_data_type), &
      INTENT(INOUT)                          :: lldata
    TYPE(debug_variables_type), &
      INTENT(OUT), OPTIONAL                  :: debug

    INTEGER                                  :: handle, i, isos, natoms, &
                                                nnodes
    LOGICAL                                  :: first_time
    REAL(dbl)                                :: pot_nonbond, vg_coulomb
    REAL(dbl), ALLOCATABLE, DIMENSION(:, :)  :: fg_coulomb
    REAL(dbl), ALLOCATABLE, &
      DIMENSION(:, :), SAVE                  :: f_nonbond
    REAL(dbl), DIMENSION(3, 3)               :: ident, pv_bc, pv_g, pv_nonbond
    TYPE(dg_type), SAVE                      :: dg
    TYPE(pw_grid_type), SAVE                 :: grid_b, grid_ewald, grid_s, &
                                                grid_spme

!------------------------------------------------------------------------------

  CALL timeset ( 'FORCE','I',' ',handle )

  nnodes = SIZE ( pnode )
  natoms = SIZE ( part )
  isos = 0

  first_time = .NOT. ALLOCATED ( f_nonbond )

  IF ( .NOT. ALLOCATED ( f_nonbond )  ) &
       ALLOCATE ( f_nonbond ( 3,natoms ), STAT = isos )
  IF ( isos /= 0 ) &
       CALL stop_memory ( 'force_env_types', 'f_nonbond', 3 * natoms )

! initialize ewalds with reference box ( consistent NPT restarts )
  IF ( first_time ) THEN

     SELECT CASE ( ewald_param % ewald_type )
     CASE ( 'EWALD' )
        CALL ewald_initialize ( dg, part, pnode, fc_global, &
             ewald_param, box_ref, thermo, ewald_grid = grid_ewald )
     CASE ( 'PME' )
        CALL ewald_initialize ( dg, part, pnode, fc_global, &
             ewald_param, box_ref, thermo, &
             pme_small_grid = grid_s, pme_big_grid = grid_b )
     CASE ( 'SPME' )
        CALL ewald_initialize ( dg, part, pnode, fc_global, &
             ewald_param, box_ref, thermo, &
             spme_grid = grid_spme )
     END SELECT

  END IF

! reinitialize the gspace for the new box
  IF ( box_change ) THEN

     CALL get_hinv ( box )

     SELECT CASE ( ewald_param % ewald_type )
     CASE ( 'EWALD' )
        CALL get_hinv ( box )
        CALL pw_grid_change ( box, grid_ewald )
        CALL ewald_setup ( grid_ewald, ewald_param, dg )
     CASE ( 'PME' )
        CALL get_hinv ( box )
        CALL pw_grid_change ( box, grid_b )
        CALL pw_green_fn_init ( box % green, grid_b )
        CALL dg_grid_change ( box, grid_b, grid_s )
        CALL ewald_setup ( grid_s, ewald_param, dg )
     CASE ( 'SPME' )
        CALL get_hinv ( box )
        CALL pw_grid_change ( box, grid_spme )
        CALL pw_green_fn_init ( box % green, grid_spme, p3m=ewald_param%o_spline, &
             alpha = ewald_param%alpha )
        CALL ewald_setup ( grid_spme, ewald_param, dg )
     END SELECT

  END IF

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
          CALL stop_memory ( 'force_env_types', 'fg_coulomb', 3 * nnodes )

! compute g-space part of the ewald sum
     SELECT CASE ( ewald_param % ewald_type )

     CASE ( "EWALD" )
        CALL ewald_evaluate ( dg, ewald_param, fg_coulomb, vg_coulomb, &
             pv_g, pnode, box )

     CASE ( "PME" )
        CALL pme_evaluate ( dg, part, box, grid_s, grid_b,  &
                           ewald_param, vg_coulomb, fg_coulomb, pv_g )
        
     CASE ( "SPME" )
        CALL spme_evaluate ( dg, part, fg_coulomb, vg_coulomb, pv_g, box, &
             grid_spme, ewald_param )
        
     CASE DEFAULT
        CALL stop_program ( "force_env_types", "illegal value of ewald_type: "//&
             ewald_param % ewald_type )
     END SELECT
  END IF
!
! add up all the potential energies
!
  IF ( ewald_param % ewald_type == 'NONE' ) THEN
     thermo % pot = pot_nonbond
#if defined ( __parallel )
     CALL mp_sum ( thermo % pot,fc_global % group )
#endif
  ELSE
     thermo % pot = pot_nonbond + thermo % e_bonded
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
     thermo % pv = pv_nonbond
#if defined ( __parallel )
     CALL mp_sum ( thermo % pv, fc_global % group )
#endif
  ELSE
     ident = 0.0_dbl
     DO i = 1, 3
        ident ( i, i ) = 1.0_dbl
     END DO

     thermo % pv = pv_nonbond + pv_bc
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
     debug % pot_nonbond = pot_nonbond
     debug % f_nonbond = f_nonbond
     debug % pv_nonbond = pv_nonbond
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
  IF ( isos /= 0 ) CALL stop_memory ( 'force_env_types', 'fg_coulomb' )

#if defined ( __parallel )
  DO i = 1, natoms
     CALL mp_sum ( part ( i ) % f, fc_global % group )
  END DO
#endif

  CALL timestop ( zero, handle )

END SUBROUTINE force_env_types

!!*****
!******************************************************************************

END MODULE tbmd_force

!******************************************************************************
