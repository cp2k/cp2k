!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/pme [1.0] *
!!
!!   NAME
!!     pme
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE pme

  USE coefficients, ONLY : integral
  USE coefficient_types, ONLY : coeff_type, coeff_allocate, &
       coeff_deallocate, coeff_zero, PW_COMPLEXDATA3D, PW_REALDATA3D, &
       PW_REALSPACE
  USE dgs, ONLY : dg_sum_patch, dg_sum_patch_force, &
       dg_get_1x_patch, dg_get_2x_patch, dg_get_strucfac, dg_get_rho_tot
  USE dg_rho0s, ONLY : dg_rho0_setup
  USE dg_types, ONLY : dg_type
  USE ewald_parameters_types, ONLY : ewald_parameters_type
  USE kinds, ONLY: dbl, sgl, YES, NO
  USE mathconstants, ONLY : fourpi
  USE molecule_types, ONLY : particle_node_type
  USE particle_lists, ONLY : particle_list_type, particle_list_select
  USE pw_grid_types, ONLY : pw_grid_type
  USE pws, ONLY : pw_poisson_solver
  USE simulation_cell, ONLY : cell_type, get_hinv
  USE stop_program, ONLY : stop_prg, stop_memory
  USE structure_factors, ONLY : structure_factor_allocate, &
       structure_factor_deallocate
  USE structure_factor_types, ONLY : structure_factor_type
  USE timings, ONLY: timeset, timestop

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: pme_evaluate, pme_setup

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************

SUBROUTINE pme_setup ( pnode, grid_s, ewald_param, dg )

  IMPLICIT NONE

! Arguments
  TYPE ( particle_node_type ), DIMENSION ( : ), INTENT ( IN ) :: pnode
  TYPE ( ewald_parameters_type ), INTENT ( IN ) :: ewald_param
  TYPE ( pw_grid_type ), INTENT ( IN ) :: grid_s
  TYPE ( dg_type ), INTENT ( OUT ) :: dg

! Locals
  INTEGER :: isos
  
!------------------------------------------------------------------------------

! This is the first (and only) double grid
  dg % grid_index = 1

  ALLOCATE ( dg % dg_rho0 % zet (1) , STAT = isos)
  IF ( isos /= 0 ) CALL stop_memory ('pme_setup', 'zet', 0 )

! No contracted Gaussians are used here
  NULLIFY ( dg % dg_rho0 % gcc )

  dg % dg_rho0 % type = ewald_param % ewald_type
  dg % dg_rho0 % zet ( 1 ) = ewald_param % alpha

  CALL particle_list_select ( pnode, dg % plist, dg % nparts )

  CALL dg_rho0_setup ( dg % dg_rho0, grid_s )

END SUBROUTINE pme_setup

!******************************************************************************

SUBROUTINE pme_evaluate ( dg, fg_coulomb, vg_coulomb, pv_g, box, &
     grid_s, grid_b, ewald_param )

  IMPLICIT NONE

! Arguments
  TYPE ( dg_type ), INTENT ( IN ), TARGET :: dg
  TYPE ( cell_type ), INTENT ( IN ) :: box
  TYPE ( pw_grid_type ), INTENT ( IN ), TARGET :: grid_s
  TYPE ( pw_grid_type ), INTENT ( IN ), TARGET :: grid_b
  TYPE ( ewald_parameters_type ), INTENT ( IN ) :: ewald_param
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( OUT ) :: fg_coulomb
  REAL ( dbl ), INTENT ( OUT ) :: vg_coulomb
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( OUT ) :: pv_g

! Locals
  LOGICAL :: doublepack
  INTEGER :: handle, part1, part2, i, allocstatus, p_idx
  INTEGER :: lb_pbc ( 3 ), ub_pbc ( 3 ), lb_dvhartree ( 3 ), lb_rhob ( 3 )
  INTEGER, DIMENSION ( : ), POINTER :: centre1, centre2
  INTEGER :: handle2, handle3
  REAL ( dbl ) :: r1 ( 3 ), r2 ( 3 ), q1, q2, force ( 3 ), inv_fourpi_eps0
  COMPLEX ( dbl ), DIMENSION ( : ), POINTER :: ex1, ey1, ez1, ex2, ey2, ez2
  TYPE ( particle_list_type ), POINTER :: plist
  TYPE ( structure_factor_type ), TARGET :: exp_igr
  TYPE ( coeff_type ), POINTER :: rho0
  TYPE ( coeff_type ) :: rhos1_r, rhos2_r, rhob_r, phi_r, dphi_r ( 3 ), rhos_r
!*apsi*  TYPE ( coeff_type ) :: rho_r
  
! For saved FFTs
  LOGICAL, SAVE :: saved_ffts_singleprecision = YES
  INTEGER :: nsaved_ffts_maxmemory = 512 !*apsi in megabytes
  INTEGER :: nsaved_ffts_max, nsaved_ffts, highest_fft_saved, nsize
  INTEGER :: size_of_rhos_r_save
  INTEGER, SAVE :: previous_size_of_rhos_r_save = -1
  
!*apsi Change to f9x intrinsic?
  INTEGER :: bytes_per_real_sgl = 4
  INTEGER :: bytes_per_real_dbl = 8
  INTEGER :: bytes_per_real
  
  INTEGER, DIMENSION ( : ), ALLOCATABLE :: pbc1_b, pbc2_b, pbc3_b
  INTEGER, DIMENSION ( : ), POINTER :: npts_b, npts_s
  INTEGER, DIMENSION ( :, : ), POINTER :: bounds_s, bounds_b
  REAL ( dbl ), DIMENSION ( :, :, : ), POINTER :: rho_b
  COMPLEX ( dbl ), DIMENSION ( :, :, : ), POINTER :: rho_s
  
! To store the values of the allocated pointers
  REAL ( dbl ), DIMENSION ( :, :, : ), POINTER :: rhos1_ptr, rhos2_ptr
  
  REAL ( sgl ), DIMENSION ( :, :, : ), POINTER :: rhos_r_sgl
  REAL ( sgl ), DIMENSION ( :, :, :, : ), ALLOCATABLE, TARGET :: &
       rhos_r_save_sgl
  REAL ( dbl ), DIMENSION ( :, :, :, : ), ALLOCATABLE, TARGET :: &
       rhos_r_save_dbl
  
  REAL ( dbl ), DIMENSION ( :, :, :, : ), ALLOCATABLE, TARGET :: &
       dvhartree
  
!------------------------------------------------------------------------------
  
  CALL timeset ( 'PME', 'I', ' ', handle )
  
  npts_b => grid_b % npts
  npts_s => grid_s % npts
  bounds_b => grid_b % bounds
  bounds_s => grid_s % bounds
  
! Determine the number of small FFTs to be stored
  IF ( saved_ffts_singleprecision ) THEN
     bytes_per_real = bytes_per_real_sgl
  ELSE
     bytes_per_real = bytes_per_real_dbl
  END IF
  
  nsize = grid_s % ngpts
  nsaved_ffts_max = ( nsaved_ffts_maxmemory * 1024 ** 2 ) &
       / ( nsize * bytes_per_real )
  
  nsaved_ffts_max = MAX ( 0, nsaved_ffts_max )
  
  nsaved_ffts = MIN ( dg % nparts, nsaved_ffts_max )
  
! Let us make it even so it is easier to deal with in the double packing
  IF ( dg % nparts > 1 .AND. nsaved_ffts < dg % nparts ) &
       nsaved_ffts = ( nsaved_ffts / 2 ) * 2
  
  ALLOCATE ( dvhartree ( 3, &
       bounds_b ( 1, 1 ):bounds_b ( 2, 1 ), &
       bounds_b ( 1, 2 ):bounds_b ( 2, 2 ), &
       bounds_b ( 1, 3 ):bounds_b ( 2, 3 ) &
        ), STAT = allocstatus )
  IF ( allocstatus /= 0 ) CALL stop_memory ( 'pme_evaluate', 'dvhartree', 0 )
  
  IF ( saved_ffts_singleprecision ) THEN
     ALLOCATE ( rhos_r_save_sgl ( bounds_s ( 1, 1 ):bounds_s ( 2, 1 ), &
          bounds_s ( 1, 2 ):bounds_s ( 2, 2 ), &
          bounds_s ( 1, 3 ):bounds_s ( 2, 3 ), &
          nsaved_ffts ), STAT = allocstatus )
  ELSE
     ALLOCATE ( rhos_r_save_dbl ( bounds_s ( 1, 1 ):bounds_s ( 2, 1 ), &
          bounds_s ( 1, 2 ):bounds_s ( 2, 2 ), &
          bounds_s ( 1, 3 ):bounds_s ( 2, 3 ), &
          nsaved_ffts ), STAT = allocstatus )
  END IF
  IF ( allocstatus /= 0 ) CALL stop_memory ( 'pme_evaluate', 'rhos_r_save', 0 )
  
  IF ( saved_ffts_singleprecision ) THEN
     size_of_rhos_r_save = SIZE ( rhos_r_save_sgl ) * bytes_per_real
  ELSE
     size_of_rhos_r_save = SIZE ( rhos_r_save_dbl ) * bytes_per_real
  END IF
! 
! this code has to be changed in order to be compliant with the output
! unit set and with parallel computing
! IF ( size_of_rhos_r_save /= previous_size_of_rhos_r_save ) THEN
!    IF ( size_of_rhos_r_save / 1024**2 > 4 ) THEN
!       WRITE ( 6, '( /, A, I4, A, / )' ) " DG| SIZE(rho) = ", &
!            size_of_rhos_r_save / 1024**2, " MB"
!    ELSE
!       WRITE ( 6, '( /, A, I4, A, / )' ) " DG| SIZE(rho) = ", &
!            size_of_rhos_r_save / 1024, " kB"
!    END IF
!    previous_size_of_rhos_r_save = size_of_rhos_r_save
! END IF
  
  lb_pbc ( : ) = bounds_b ( 1, : ) + bounds_s ( 1, : )
  ub_pbc ( : ) = bounds_b ( 2, : ) + bounds_s ( 2, : ) + 2
  ALLOCATE ( pbc1_b ( lb_pbc ( 1 ) : ub_pbc ( 1 ) ), STAT = allocstatus )
  IF ( allocstatus /= 0 ) CALL stop_memory ( 'pme_evaluate', 'pbc1_b', 0 )
  ALLOCATE ( pbc2_b ( lb_pbc ( 2 ) : ub_pbc ( 2 ) ), STAT = allocstatus )
  IF ( allocstatus /= 0 ) CALL stop_memory ( 'pme_evaluate', 'pbc2_b', 0 )
  ALLOCATE ( pbc3_b ( lb_pbc ( 3 ) : ub_pbc ( 3 ) ), STAT = allocstatus )
  IF ( allocstatus /= 0 ) CALL stop_memory ( 'pme_evaluate', 'pbc3_b', 0 )
  
  DO i = lb_pbc ( 1 ), ub_pbc ( 1 )
     IF ( i < bounds_b ( 1, 1 ) ) THEN
        pbc1_b ( i ) = i + npts_b ( 1 )
     ELSE IF ( i > bounds_b ( 2, 1 ) ) THEN
        pbc1_b ( i ) = i - npts_b ( 1 )
     ELSE
        pbc1_b ( i ) = i
     END IF
  END DO
  DO i = lb_pbc ( 2 ), ub_pbc ( 2 )
     IF ( i < bounds_b ( 1, 2 ) ) THEN
        pbc2_b ( i ) = i + npts_b ( 2 )
     ELSE IF ( i > bounds_b ( 2, 2 ) ) THEN
        pbc2_b ( i ) = i - npts_b ( 2 )
     ELSE
        pbc2_b ( i ) = i
     END IF
  END DO
  DO i = lb_pbc ( 3 ), ub_pbc ( 3 )
     IF ( i < bounds_b ( 1, 3 ) ) THEN
        pbc3_b ( i ) = i + npts_b ( 3 )
     ELSE IF ( i > bounds_b ( 2, 3 ) ) THEN
        pbc3_b ( i ) = i - npts_b ( 3 )
     ELSE
        pbc3_b ( i ) = i
     END IF
  END DO
  
!!------------- DENSITY CALCULATION ----------------
  
! initializing TOTAL rho_r
  inv_fourpi_eps0 = 1.0_dbl / ( fourpi * ewald_param % eps0 )
  
!*apsi  CALL coeff_allocate ( rho_r, grid_b, use_data = PW_COMPLEXDATA3D )
!*apsi  CALL coeff_zero ( rho_r )
  
  rho0 => dg % dg_rho0 % density
  
! initialize rho_r on big grid
  CALL coeff_allocate ( rhob_r, grid_b, use_data = PW_REALDATA3D )
  rhob_r % pw % in_space = PW_REALSPACE
  CALL coeff_allocate ( rhos_r, grid_s, use_data = PW_COMPLEXDATA3D )
  CALL coeff_allocate ( rhos1_r, grid_s, use_data = PW_REALDATA3D )
  CALL coeff_allocate ( rhos2_r, grid_s, use_data = PW_REALDATA3D )
  
  rho_b => rhob_r % pw % cr3d
  rho_s => rhos_r % pw % cc3d
  rhos1_ptr => rhos1_r % pw % cr3d
  rhos2_ptr => rhos2_r % pw % cr3d
  
  lb_rhob = LBOUND ( rhob_r % pw % cr3d )
  
  CALL coeff_zero ( rhob_r )
  
! Allocate the space for the structure factor
  CALL structure_factor_allocate ( grid_s % bounds, dg % nparts,  &
       exp_igr, allocate_centre = .TRUE. )
  
  IF (.NOT. ASSOCIATED ( dg % plist ) ) &
       CALL stop_prg ( 'pme_evaluate', 'no charged particles')
  
  CALL timeset ( 'PME_DENSITY', 'I', ' ', handle2 )
  
  p_idx = 1
  plist => dg % plist
  DO
     doublepack = ASSOCIATED ( plist % next )
     
     r1 ( : ) = plist % r ( : )
     q1 = plist % charge
     part1 = plist % particle_index
     centre1 => exp_igr % centre ( :, part1 )
     ex1 => exp_igr % ex ( :, part1 )
     ey1 => exp_igr % ey ( :, part1 )
     ez1 => exp_igr % ez ( :, part1 )
     IF ( doublepack ) THEN
        plist => plist % next
        r2 ( : ) = plist % r ( : )
        q2 = plist % charge
        part2 = plist % particle_index
        centre2 => exp_igr % centre ( :, part2 )
        ex2 => exp_igr % ex ( :, part2 )
        ey2 => exp_igr % ey ( :, part2 )
        ez2 => exp_igr % ez ( :, part2 )
     END IF
     
!    CALL timeset ( 'PME_GETSFAC', 'I', ' ', handle3 )
     CALL dg_get_strucfac ( box, r1, grid_s % npts, grid_b % npts, &
          centre1, exp_igr % lb, ex1, ey1, ez1 )
     IF ( doublepack ) THEN
        CALL dg_get_strucfac ( box, r2, grid_s % npts, grid_b % npts, &
             centre2, exp_igr % lb, ex2, ey2, ez2 )
     END IF
!    CALL timestop ( 0.0_dbl, handle3 )
     
!    CALL timeset ( 'PME_GETPATCH', 'I', ' ', handle3 )
     IF ( doublepack ) THEN
        CALL dg_get_2x_patch ( rho0, rhos_r, q1, q2, &
             exp_igr % lb, ex1, ey1, ez1, ex2, ey2, ez2 )
     ELSE
        CALL dg_get_1x_patch ( rho0, rhos_r, q1, &
             exp_igr % lb, ex1, ey1, ez1 )
     END IF
!    CALL timestop ( 0.0_dbl, handle3 )
     
     IF ( p_idx <= nsaved_ffts ) THEN
        IF ( saved_ffts_singleprecision ) THEN
           rhos_r_save_sgl ( :, :, :, p_idx ) = REAL ( rho_s )
           IF ( doublepack ) THEN
              rhos_r_save_sgl ( :, :, :, p_idx + 1 ) = AIMAG ( rho_s )
           END IF
        ELSE
           rhos_r_save_dbl ( :, :, :, p_idx ) = REAL ( rho_s )
           IF ( doublepack ) THEN
              rhos_r_save_dbl ( :, :, :, p_idx + 1 ) = AIMAG ( rho_s )
           END IF
        END IF
     END IF
     
     p_idx = p_idx + 2
     
     rhos1_ptr = REAL ( rho_s )
     IF ( doublepack ) THEN
        rhos2_ptr = AIMAG ( rho_s )
     END IF
     
!    CALL timeset ( 'PME_SUMPATCH1', 'I', ' ', handle3 )
     CALL dg_sum_patch ( rhos1_ptr, rho_b, bounds_s, centre1, &
          lb_pbc, pbc1_b, pbc2_b, pbc3_b, lb_rhob )
     IF ( doublepack ) THEN
        CALL dg_sum_patch ( rhos2_ptr, rho_b, bounds_s, centre2, &
             lb_pbc, pbc1_b, pbc2_b, pbc3_b, lb_rhob )
     END IF
!    CALL timestop ( 0.0_dbl, handle3 )
     
     IF ( .NOT. ASSOCIATED ( plist % next ) ) THEN
        EXIT
     ELSE
        plist => plist % next
     END IF
  END DO
  
! Obtain total density on rho_r; not necessary in PME because rho_r and rhob_r
!    can be the same...
!*apsi  CALL dg_get_rho_tot ( rho_r, rhob_r )
  
  CALL timestop ( 0.0_dbl, handle2 )
  
!*apsi  CALL coeff_deallocate ( rhob_r )
  
! initialize the forces
  
!!------------- ELECTROSTATIC CALCULATION -----------
  
  DO i = 1, 3
     CALL coeff_allocate ( dphi_r ( i ), grid_b, use_data = PW_REALDATA3D )
  END DO
  CALL coeff_allocate ( phi_r, grid_b, use_data = PW_REALDATA3D )
  
  lb_dvhartree = LBOUND ( dphi_r ( 1 ) % pw % cr3d )
  
  CALL timeset ( 'PME_POISSON', 'I', ' ', handle2 )
  
!*apsi  CALL pw_poisson_solver ( rho_r, phi_r, dphi_r )
  CALL pw_poisson_solver ( rhob_r, phi_r, dphi_r )
  
!  dvhartree ( :, :, :, 1 ) = dphi_r ( 1 ) % pw % cr3d
!  dvhartree ( :, :, :, 2 ) = dphi_r ( 2 ) % pw % cr3d
!  dvhartree ( :, :, :, 3 ) = dphi_r ( 3 ) % pw % cr3d
  dvhartree ( 1, :, :, : ) = dphi_r ( 1 ) % pw % cr3d
  dvhartree ( 2, :, :, : ) = dphi_r ( 2 ) % pw % cr3d
  dvhartree ( 3, :, :, : ) = dphi_r ( 3 ) % pw % cr3d
  
!*apsi  vg_coulomb = 0.5_dbl * integral ( rho_r, phi_r ) * inv_fourpi_eps0 &
!*apsi       / grid_b % vol
  vg_coulomb = 0.5_dbl * integral ( rhob_r, phi_r ) * inv_fourpi_eps0
  
  CALL timestop ( 0.0_dbl, handle2 )
  
!*apsi  CALL coeff_deallocate ( rho_r )
  CALL coeff_deallocate ( rhob_r )
  CALL coeff_deallocate ( phi_r )
  
  fg_coulomb = 0.0_dbl
  
  CALL timeset ( 'PME_FORCES', 'I', ' ', handle2 )
  
  p_idx = 1
  plist => dg % plist
  DO
     doublepack = ASSOCIATED ( plist % next )
     
     q1 = plist % charge
     part1 = plist % particle_index
     centre1 => exp_igr % centre ( :, part1 )
     ex1 => exp_igr % ex ( :, part1 )
     ey1 => exp_igr % ey ( :, part1 )
     ez1 => exp_igr % ez ( :, part1 )
     IF ( doublepack ) THEN
        plist => plist % next
        q2 = plist % charge
        part2 = plist % particle_index
        centre2 => exp_igr % centre ( :, part2 )
        ex2 => exp_igr % ex ( :, part2 )
        ey2 => exp_igr % ey ( :, part2 )
        ez2 => exp_igr % ez ( :, part2 )
     END IF
     
! Either fetch the densities from memory ...
     IF ( p_idx <= nsaved_ffts ) THEN
        IF ( saved_ffts_singleprecision ) THEN
           rhos1_ptr = rhos_r_save_sgl ( :, :, :, p_idx )
           IF ( doublepack ) THEN
              rhos2_ptr = rhos_r_save_sgl ( :, :, :, p_idx + 1 )
           END IF
        ELSE
           rhos1_ptr = rhos_r_save_dbl ( :, :, :, p_idx )
           IF ( doublepack ) THEN
              rhos2_ptr = rhos_r_save_dbl ( :, :, :, p_idx + 1 )
           END IF
        END IF
        
     ELSE
! ... they have to be calculated
        IF ( doublepack ) THEN
           CALL dg_get_2x_patch ( rho0, rhos_r, q1, q2, &
                exp_igr % lb, ex1, ey1, ez1, ex2, ey2, ez2 )
        ELSE
           CALL dg_get_1x_patch ( rho0, rhos_r, q1, &
                exp_igr % lb, ex1, ey1, ez1 )
        END IF
        
        rhos1_ptr = REAL ( rho_s )
        IF ( doublepack ) THEN
           rhos2_ptr = AIMAG ( rho_s )
        END IF
     END IF
     
     p_idx = p_idx + 2
     
!    CALL timeset ( 'PME_SUMPATCH2', 'I', ' ', handle3 )
     CALL dg_sum_patch_force ( rhos1_ptr, dvhartree, force, &
          bounds_s, centre1, &
          lb_pbc, pbc1_b, pbc2_b, pbc3_b, lb_dvhartree )
     fg_coulomb ( :, part1 ) = fg_coulomb ( :, part1 ) &
          - force ( : ) * rhos1_r % pw % pw_grid % dvol * inv_fourpi_eps0
     
     IF ( doublepack ) THEN
        CALL dg_sum_patch_force ( rhos2_ptr, dvhartree, force, &
             bounds_s, centre2, &
             lb_pbc, pbc1_b, pbc2_b, pbc3_b, lb_dvhartree )
        fg_coulomb ( :, part2 ) = fg_coulomb ( :, part2 ) &
             - force ( : ) * rhos2_r % pw % pw_grid % dvol * inv_fourpi_eps0
     END IF
!    CALL timestop ( 0.0_dbl, handle3 )
     
     IF ( .NOT. ASSOCIATED ( plist % next ) ) THEN
        EXIT
     ELSE
        plist => plist % next
     END IF
  END DO
  
  CALL timestop ( 0.0_dbl, handle2 )
  
  DO i = 1, 3
     CALL coeff_deallocate ( dphi_r ( i ) )
  END DO
  
  CALL coeff_deallocate ( rhos2_r )
  CALL coeff_deallocate ( rhos1_r )
  CALL coeff_deallocate ( rhos_r )
  
  CALL structure_factor_deallocate ( exp_igr )
  
  DEALLOCATE ( dvhartree, STAT = allocstatus )
  IF ( allocstatus /= 0 ) CALL stop_memory ( 'pme_evaluate', 'dvhartree' )
  
  DEALLOCATE ( pbc1_b, STAT = allocstatus )
  IF ( allocstatus /= 0 ) CALL stop_memory ( 'pme_evaluate', 'pbc1_b' )
  DEALLOCATE ( pbc2_b, STAT = allocstatus )
  IF ( allocstatus /= 0 ) CALL stop_memory ( 'pme_evaluate', 'pbc2_b' )
  DEALLOCATE ( pbc3_b, STAT = allocstatus )
  IF ( allocstatus /= 0 ) CALL stop_memory ( 'pme_evaluate', 'pbc3_b' )
  
  IF ( saved_ffts_singleprecision ) THEN
     DEALLOCATE ( rhos_r_save_sgl, STAT = allocstatus )
  ELSE
     DEALLOCATE ( rhos_r_save_dbl, STAT = allocstatus )
  END IF
  IF ( allocstatus /= 0 ) CALL stop_memory ( 'pme_evaluate', 'rhos_r_save' )

! The contribution from the electro-statistics on the pressure tensor has
!    not yet been implemented
  pv_g = -1.0E+30_dbl

  CALL timestop ( 0.0_dbl, handle )

END SUBROUTINE pme_evaluate

!******************************************************************************

END MODULE pme
