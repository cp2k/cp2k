!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C ) 2000  CP2K developers group                                !
!-----------------------------------------------------------------------------!
!!****** cp2k/pol_force [1.0] *
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

  USE atomic_kinds,                    ONLY: kind_info_type
  USE coefficient_types,               ONLY: coeff_type
  USE dg_types,                        ONLY: dg_type
  USE empirical_parameters,            ONLY: empirical_parameter_type
  USE ewald_parameters_types,          ONLY: ewald_parameters_type
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
  USE pol_electrostatics,              ONLY: electrostatics
  USE pol_grids,                       ONLY: assign_coefs_to_grids,&
                                             pol_grids_allocate,&
                                             pol_grids_count,&
                                             pol_grids_initialize,&
                                             pol_grids_map
  USE pol_overlap,                     ONLY: force_overlap
  USE pol_thk,                         ONLY: force_thk
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE simulation_cell,                 ONLY: cell_type
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop

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
  TYPE ( dg_type ), POINTER, DIMENSION ( : ) :: dg_coef
  TYPE ( dg_type ), POINTER, DIMENSION ( : ) :: dg_atom
  REAL ( dbl ), DIMENSION ( :, : ), ALLOCATABLE :: f_nonbond
  INTEGER :: iref

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** pol_force/force_env_types [1.0] *
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


    TYPE(molecule_structure_type), &
      DIMENSION(:), INTENT(IN)               :: molecule
    TYPE(particle_node_type), DIMENSION(:), &
      INTENT(INOUT)                          :: pnode
    TYPE(particle_type), DIMENSION(:), &
      INTENT(INOUT)                          :: part
    TYPE(cell_type), INTENT(INOUT)           :: box
    TYPE(cell_type), INTENT(IN)              :: box_ref
    TYPE(kind_info_type), DIMENSION(:), &
      POINTER                                :: drho_basis_info, &
                                                rho0_basis_info
    TYPE(coeff_type), INTENT(inout)          :: coef_pos, coef_vel, coef_force
    TYPE(thermodynamic_type), INTENT(INOUT)  :: thermo
    TYPE(potentialparm_type), &
      DIMENSION(:, :), INTENT(IN)            :: potparm
    TYPE(empirical_parameter_type), &
      DIMENSION(:), INTENT(IN)               :: empparm
    TYPE(ewald_parameters_type), &
      INTENT(INOUT)                          :: ewald_param
    LOGICAL, INTENT(IN)                      :: box_change
    TYPE(global_environment_type), &
      INTENT(IN)                             :: fc_global
    TYPE(linklist_internal_data_type), &
      DIMENSION(:), INTENT(INOUT)            :: lldata
    TYPE(debug_variables_type), &
      INTENT(OUT), OPTIONAL                  :: debug

    INTEGER                                  :: handle, i, igrid, ii, isos, &
                                                natoms, ncoefs, ng, ng_atom, &
                                                ng_coef, nnodes, nt_atom, &
                                                nt_coef
    LOGICAL                                  :: first_time
    REAL(dbl)                                :: pot_electrostatics, &
                                                pot_nonbond, pot_overlap, &
                                                pot_thk
    REAL(dbl), ALLOCATABLE, DIMENSION(:)     :: fe_coef, fo_coef, fthk
    REAL(dbl), ALLOCATABLE, DIMENSION(:, :)  :: f_total, fe_part, fo_part
    REAL(dbl), DIMENSION(3, 3)               :: pv_nonbond

!------------------------------------------------------------------------------

  CALL timeset ( 'FORCE','I',' ',handle )

  nnodes = SIZE ( pnode )
  natoms = SIZE ( part )
  ncoefs = size ( coef_pos % ao % cr )
  isos = 0
  pot_overlap = 0._dbl
  pot_electrostatics = 0._dbl
  pot_nonbond = 0._dbl
  pot_thk = 0

  first_time = .NOT. ALLOCATED ( f_nonbond )

  IF ( .NOT. ALLOCATED ( f_nonbond )  ) THEN
    ALLOCATE ( f_nonbond ( 3,natoms ), STAT = isos )
    IF ( isos /= 0 ) &
       CALL stop_memory ( 'force_env_types', 'f_nonbond', 3 * natoms )
  ELSE IF ( SIZE ( f_nonbond ,2 ) < natoms ) THEN
    DEALLOCATE ( f_nonbond, STAT = isos )
    IF ( isos /= 0 ) CALL stop_memory ( 'force_env_types', 'f_nonbond' )
    ALLOCATE ( f_nonbond ( 3, natoms ), STAT = isos )
    IF ( isos /= 0 ) &
       CALL stop_memory ( 'force_env_types', 'f_nonbond', 3 * natoms )
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
      CALL stop_memory ( 'force_env_types', 'fe_part', 3 * nnodes )

  IF ( .NOT. ALLOCATED ( fo_part )  ) &
         ALLOCATE ( fo_part ( 3, nnodes ), STAT=isos )
  IF ( isos /= 0 ) &
    CALL stop_memory ( 'force_env_types', 'fo_part',  3 * nnodes )

  IF ( .NOT. ALLOCATED ( fe_coef )  ) &
         ALLOCATE ( fe_coef ( ncoefs ), STAT=isos )
  IF ( isos /= 0 ) &
    CALL stop_memory ( 'force_env_types', 'fe_coef',  ncoefs )

  IF ( .NOT. ALLOCATED ( fo_coef )  ) &
         ALLOCATE ( fo_coef ( ncoefs ), STAT=isos )
  IF ( isos /= 0 ) &
    CALL stop_memory ( 'force_env_types', 'fo_coef',  ncoefs )

  IF ( .NOT. ALLOCATED ( fthk )  ) &
         ALLOCATE ( fthk ( ncoefs ), STAT=isos )
  IF ( isos /= 0 ) &
    CALL stop_memory ( 'force_env_types', 'fthk',  ncoefs )

! get overlap forces

   CALL force_overlap( coef_pos % ao, drho_basis_info, part, pnode, box, & 
                        empparm, pot_overlap, fo_coef, fo_part )

! get first order hk contribution

   CALL force_thk( coef_pos % ao, empparm, fthk, pot_thk)

! get electrostatics forces   

   CALL electrostatics (dg_atom, dg_coef, part, coef_pos, box,  &
                        small_grid, big_grid, ewald_param, iref, &
                        pot_electrostatics, fe_coef, fe_part )
 
!
! get real-space non-bonded forces:
!
  f_nonbond = zero
  CALL force_nonbond ( ewald_param,pnode,box,potparm, &
       pot_nonbond,f_nonbond,pv_nonbond )

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
  IF ( isos /= 0 ) CALL stop_memory ( 'force_env_types', 'f_total', 3 * natoms )

  DO i = 1, natoms
     f_total ( 1, i ) = part ( i ) % f ( 1 ) + f_nonbond ( 1, i )
     f_total ( 2, i ) = part ( i ) % f ( 2 ) + f_nonbond ( 2, i )
     f_total ( 3, i ) = part ( i ) % f ( 3 ) + f_nonbond ( 3, i )
  END DO

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

  thermo % e_class = pot_nonbond
  thermo % e_electro = pot_electrostatics
  thermo % e_overlap = pot_overlap
  thermo % e_thk = pot_thk
  thermo % pot = pot_nonbond + pot_electrostatics + pot_overlap + pot_thk

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
  IF ( isos /= 0 ) CALL stop_memory ( 'force_env_types', 'f_total' )
  IF ( ALLOCATED ( fe_part )  ) DEALLOCATE ( fe_part, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'force_env_types', 'fe_part' )
  IF ( ALLOCATED ( fo_part )  ) DEALLOCATE ( fo_part, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'force_env_types', 'fo_part' )
  IF ( ALLOCATED ( fe_coef )  ) DEALLOCATE ( fe_coef, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'force_env_types', 'fe_coef' )
  IF ( ALLOCATED ( fo_coef )  ) DEALLOCATE ( fo_coef, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'force_env_types', 'fo_coef' )
  IF ( ALLOCATED ( fthk )  ) DEALLOCATE ( fthk, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'force_env_types', 'fthk' )

  CALL timestop ( zero, handle )

END SUBROUTINE force_env_types

!******************************************************************************
!!****** pol_force/get_pol_grids [1.0] *
!!
!!   NAME
!!     get_pol_grids
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE get_pol_grids ( pw_big, pw_small, dg1, dg2, ii )
    TYPE(pw_grid_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: pw_big, pw_small
    TYPE(dg_type), DIMENSION(:), OPTIONAL, &
      POINTER                                :: dg1, dg2
    INTEGER, INTENT(OUT), OPTIONAL           :: ii

  IF ( PRESENT ( pw_big ) ) pw_big  => big_grid 
  IF ( PRESENT ( pw_small ) ) pw_small => small_grid 
  IF ( PRESENT ( dg1 ) ) dg1 => dg_atom 
  IF ( PRESENT ( dg2 ) ) dg2 => dg_coef 
  IF ( PRESENT ( ii ) ) ii = iref

END SUBROUTINE get_pol_grids
!!*****
!******************************************************************************

END MODULE pol_force

!******************************************************************************
