!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/dft_structure_control [1.0] *
!!
!!   NAME
!!     dft_structure_control
!!
!!   FUNCTION
!!     control the calculation of the electronic structure via DFT
!!
!!   AUTHOR
!!     CJM & JGH
!!
!!   MODIFICATION HISTORY
!!     NONE
!!
!!   SOURCE
!******************************************************************************

MODULE dft_structure_control
  
  USE brillouin,                       ONLY: kpoint_type
  USE cp_control_types,                ONLY: dft_control_type
  USE densities,                       ONLY: density_update
  USE density_types,                   ONLY: density_allocate,&
                                             density_init,&
                                             density_type
  USE ewald_parameters_types,          ONLY: ewald_parameters_type
  USE global_types,                    ONLY: global_environment_type
  USE greens_fn,                       ONLY: pw_green_fn_init
  USE kinds,                           ONLY: NO,&
                                             YES,&
                                             dbl
  USE linklist_types,                  ONLY: linklist_internal_data_type
  USE md,                              ONLY: thermodynamic_type
  USE molecule_types,                  ONLY: molecule_structure_type,&
                                             particle_node_type
  USE pair_potential_types,            ONLY: pair_potential_type
  USE particle_types,                  ONLY: particle_type
  USE potential_types,                 ONLY: potential_allocate,&
                                             potential_init,&
                                             potential_type
  USE potentials,                      ONLY: potential_update
  USE pw_grid_types,                   ONLY: FULLSPACE,&
                                             pw_grid_type
  USE pw_grids,                        ONLY: pw_grid_construct,&
                                             pw_grid_setup
  USE pw_types,                        ONLY: PW_REALDATA3D => REALDATA3D
  USE pws,                             ONLY: init_pw_poisson_solver
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_hinv
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE wf_types,                        ONLY: wf_allocate,&
                                             wf_type
  USE wfs,                             ONLY: wf_randomise,&
                                             wf_update,&
                                             wf_write_restart

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: dft_force_control
  
CONTAINS

!******************************************************************************

SUBROUTINE dft_force_control ( molecule, pnode, part, box, thermo, &
     potparm, ewald_param, box_change, lldata, dft_control, kp, &
     dftpar, wavepar )
  
  
    TYPE(molecule_structure_type), &
      DIMENSION(:), INTENT(IN)               :: molecule
    TYPE(particle_node_type), DIMENSION(:), &
      INTENT(INOUT)                          :: pnode
    TYPE(particle_type), DIMENSION(:), &
      INTENT(INOUT)                          :: part
    TYPE(cell_type), INTENT(INOUT)           :: box
    TYPE(thermodynamic_type), INTENT(INOUT)  :: thermo
    TYPE(pair_potential_type), &
      DIMENSION(:, :), INTENT(IN)            :: potparm
    TYPE(ewald_parameters_type), &
      INTENT(INOUT)                          :: ewald_param
    LOGICAL, INTENT(IN)                      :: box_change
    TYPE(linklist_internal_data_type), &
      INTENT(INOUT)                          :: lldata
    TYPE(dft_control_type), INTENT(IN)       :: dft_control
    TYPE(kpoint_type), INTENT(IN)            :: kp
    TYPE(global_environment_type), &
      INTENT(IN)                             :: dftpar, wavepar

    REAL(dbl), ALLOCATABLE, DIMENSION(:, :)  :: fion, tau
    TYPE(density_type), DIMENSION(:), &
      POINTER                                :: n
    TYPE(potential_type), DIMENSION(:), &
      POINTER                                :: vloc
    TYPE(wf_type), DIMENSION(:), POINTER     :: psi

!TMPTMPTMP
!------------------------------------------------------------------------------

  CALL dft_initialise ( box, psi, n, vloc, dft_control, kp, wavepar )
  
  ALLOCATE ( tau ( 3, 1 ), fion ( 3, 1 ) ) !TMPTMPTMP
  
  CALL init_pw_poisson_solver ( box % green )

  DO
     CALL dft_scf_loop ( psi, n, vloc, tau, fion, dft_control, kp, &
          dftpar, wavepar )
     ! up-date the positions etc...
     
     EXIT
  END DO
  
  CALL dft_finalise ( dft_control )
  
END SUBROUTINE dft_force_control

!******************************************************************************

SUBROUTINE dft_initialise ( cell, psi, n, vloc, dft_control, kp, wavepar )
  
  
    TYPE(cell_type)                          :: cell
    TYPE(wf_type), DIMENSION(:), POINTER     :: psi
    TYPE(density_type), DIMENSION(:), &
      POINTER                                :: n
    TYPE(potential_type), DIMENSION(:), &
      POINTER                                :: vloc
    TYPE(dft_control_type)                   :: dft_control
    TYPE(kpoint_type), INTENT(IN)            :: kp
    TYPE(global_environment_type), &
      INTENT(IN)                             :: wavepar

    INTEGER                                  :: densityset, kpt, &
                                                ndensitysets, nkpts, nspins, &
                                                nstates, nwfsets, wfset
    INTEGER, DIMENSION(:, :), POINTER        :: bounds
    REAL(dbl)                                :: gcut_density, gcut_wf
    TYPE(density_type), POINTER              :: n_n
    TYPE(potential_type), POINTER            :: vloc_n
    TYPE(pw_grid_type), DIMENSION(:), &
      POINTER                                :: pw_grid
    TYPE(wf_type), POINTER                   :: psi_n

! TMPTMPTMP
!------------------------------------------------------------------------------

  CALL get_hinv ( cell )
  
  nkpts = kp % nkpt
  nstates = dft_control % nstates
  nspins = dft_control % nspins
  
!TMPTMPTMP
  nwfsets = 1
  ndensitysets = nwfsets
  
  write ( unit = *, fmt = '()' )
  write ( unit = *, fmt = * ) "E(cut), dual = ", &
       dft_control % wave_control % e_cutoff_wf, &
       dft_control % wave_control % e_cutoff_dual
  write ( unit = *, fmt = * ) "#st, #k, #s ", dft_control % nstates, &
       nkpts, dft_control % nspins
  do kpt = 1, kp % nkpt
     write ( unit = *, fmt = '( a, i4, 4f10.6 )' ) &
          "k #", kpt, kp % xk ( :, kpt ), kp % weight ( kpt )
  end do
  write ( unit = *, fmt = '()' )
!TMPTMPTMP
  
  gcut_wf = SQRT ( dft_control % wave_control % e_cutoff_wf )
  gcut_density = SQRT ( dft_control % wave_control % e_cutoff_wf &
       * dft_control % wave_control % e_cutoff_dual )
  
!! Currently we need different grids only for wave functions and densities;
!!    potential uses the set of the density
  ALLOCATE ( pw_grid ( nwfsets + ndensitysets ) )
  
  ALLOCATE ( psi ( nwfsets ) )
  ALLOCATE ( n ( ndensitysets ) )
  ALLOCATE ( vloc ( ndensitysets ) )
  
  DO densityset = 1, ndensitysets
     
     CALL pw_grid_construct ( pw_grid ( densityset ) )
!     grid ( densityset ) % grid_span = HALFSPACE
     pw_grid ( densityset ) % grid_span = FULLSPACE
!     pw_grid ( densityset ) % bounds ( :, : ) = mybounds ( :, : )
     pw_grid ( densityset ) % bounds ( :, : ) = 0
     
     pw_grid ( densityset ) % npts ( : ) = &
          pw_grid ( densityset ) % bounds ( 2, : ) &
          - pw_grid ( densityset ) % bounds ( 1, : ) + 1
     pw_grid ( densityset ) % ndimension = 3
     pw_grid ( densityset ) % dr ( : ) = &
          SQRT ( SUM ( cell % hmat ( :, : ), 1 ) ) &
          / REAL ( pw_grid ( densityset ) % npts ( : ), dbl )
     pw_grid ( densityset ) % dr ( : ) = &
          SQRT ( SUM ( cell % hmat ( :, : ), 1 ) ) &
          / REAL ( pw_grid ( densityset ) % npts ( : ), dbl )
     
     CALL pw_grid_setup ( cell, pw_grid ( densityset ), cutoff = gcut_density )
     
     n_n => n ( densityset )
     CALL density_allocate ( n_n, nspins, npart = 1 )
     CALL density_init ( n_n, part = 1, pw_grid = pw_grid ( densityset ), &
          use_data = PW_REALDATA3D )
     
     vloc_n => vloc ( densityset )
     CALL potential_allocate ( vloc_n, nspins, npart = 1 )
     CALL potential_init ( vloc_n, part = 1, &
          pw_grid = pw_grid ( densityset ), use_data = PW_REALDATA3D )
  END DO
  
  DO wfset = ndensitysets + 1, ndensitysets + nwfsets
     
     CALL pw_grid_construct ( pw_grid ( wfset ) )
     ! pw_grid ( wfset ) % grid_span = HALFSPACE
     pw_grid ( wfset ) % grid_span = FULLSPACE
     pw_grid ( wfset ) % bounds ( 1, : ) = pw_grid ( 1 ) % bounds ( 1, : )
     pw_grid ( wfset ) % bounds ( 2, : ) = pw_grid ( 1 ) % bounds ( 2, : )
     pw_grid ( wfset ) % npts ( : ) = pw_grid ( wfset ) % bounds ( 2, : ) &
          - pw_grid ( wfset ) % bounds ( 1, : ) + 1
     pw_grid ( wfset ) % ndimension = 3
     pw_grid ( wfset ) % dr ( : ) = SQRT ( SUM ( cell % hmat ( :, : ), 1 ) ) &
          / REAL ( pw_grid ( wfset ) % npts ( : ), dbl )
     pw_grid ( wfset ) % dr ( : ) = SQRT ( SUM ( cell % hmat ( :, : ), 1 ) ) &
          / REAL ( pw_grid ( wfset ) % npts ( : ), dbl )
     
     CALL pw_grid_setup ( cell, pw_grid ( wfset ), cutoff = gcut_wf )
     ! CALL pw_grid_setup ( cell, pw_grid ( wfset ) )
     
     psi_n => psi ( wfset - ndensitysets )
     CALL wf_allocate ( nstates, nkpts, nspins, psi_n, pw_grid ( wfset ) )
     
     CALL wf_randomise ( psi_n, wavepar )
     
  END DO

! the Greens function
! choose the set that is used to solve the Poisson equation in G space
  densityset = 1
  CALL pw_green_fn_init ( cell % green, pw_grid ( densityset ) )
  
END SUBROUTINE dft_initialise

!******************************************************************************

SUBROUTINE dft_finalise ( dft_control )
  
    TYPE(dft_control_type)                   :: dft_control

!!!!! if the pw_grid is deallocated -> do not forget to get rid of the Greens function

END SUBROUTINE dft_finalise

!******************************************************************************

SUBROUTINE dft_scf_loop ( psi, n, vloc, tau, fion, dft_control, kp, &
     dftpar, wavepar )
  
  
    TYPE(wf_type), DIMENSION(:), &
      INTENT(INOUT)                          :: psi
    TYPE(density_type), DIMENSION(:), &
      INTENT(INOUT)                          :: n
    TYPE(potential_type), DIMENSION(:), &
      INTENT(INOUT), TARGET                  :: vloc
    REAL(dbl), DIMENSION(:, :), INTENT(IN)   :: tau
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: fion
    TYPE(dft_control_type), INTENT(IN)       :: dft_control
    TYPE(kpoint_type), INTENT(IN)            :: kp
    TYPE(global_environment_type), &
      INTENT(IN)                             :: dftpar, wavepar

    INTEGER                                  :: iter_elc, niter_elc
    INTEGER, POINTER                         :: bds( :, : )
    LOGICAL                                  :: tlastloop, tprint, &
                                                tultimateloop

!TMPTMPTMP  TYPE ( potential_type ), DIMENSION ( : ), INTENT ( INOUT ) :: vloc
!! Last loop of (adiabatic) electronic structure,
!!    maybe/probably we shall still return
!! Final loop, no return here anymore
!  LOGICAL, SAVE :: first_time = .true.
!------------------------------------------------------------------------------

  IF ( dftpar % print_level > 2 ) &
       WRITE ( dftpar % scr, * ) "Entering scf loop..."
  
!  IF ( first_time .EQV. YES ) THEN
!     CALL initialise_electronics()
!  END IF
  
  !TMPTMPTMP?
  tprint = .TRUE.
  IF ( dftpar % print_level > 5 ) &
       WRITE ( dftpar % scr, * ) "Going for density & potential update..."
  CALL density_update ( psi, n, dft_control, tprint, dftpar )
  CALL potential_update ( n, vloc, dft_control, tprint )
  
  tlastloop = NO
  tultimateloop = NO
  
  iter_elc = 0
  
  main_loop: DO
     
     iter_elc = iter_elc + 1
     
     IF ( iter_elc >= dft_control % wave_control % niter_elc &
          .OR. tultimateloop ) tlastloop = YES
     
     tprint = MOD ( iter_elc, &
          dft_control % wave_control % wf_print_frequency ) == 0
     
     !TMPTMPTMP >>>
     bds => vloc ( 1 ) % pot_part ( 1, 1 ) % pw % pw_grid % bounds
     !  vloc ( 1 ) % vloc_s ( 1 ) % cr3d = 0._dbl
     !  vloc ( 1 ) % vloc_s ( SIZE ( vloc ( 1 ) % vloc_s ) ) % cr3d = 0._dbl
     vloc(1) % pot_part(1,1) % pw % cr3d ( bds(1,1)+2, bds(1,2)+1, bds(1,3)+0 ) = &
          vloc(1)%pot_part(1,1)%pw% cr3d ( bds(1,1)+2, bds(1,2)+1, bds(1,3)+0 ) &
          - 15._dbl
     !TMPTMPTMP <<<
     
     IF ( dftpar % print_level > 5 ) &
          WRITE ( dftpar % scr, * ) "Going for wf update..."
     CALL wf_update ( vloc, psi, dft_control, dftpar, tprint )
     
     IF ( dftpar % print_level > 5 ) &
          WRITE ( dftpar % scr, * ) "Going for density update..."
     CALL density_update ( psi, n, dft_control, tprint, dftpar )
     
     IF ( dftpar % print_level > 5 ) &
          WRITE ( dftpar % scr, * ) "Going for potential update..."
     CALL potential_update ( n, vloc, dft_control, tprint )
     
     IF ( tprint .OR. tlastloop ) THEN
        OPEN ( UNIT = 99, FILE = "wfs", FORM = "UNFORMATTED" )
        CALL wf_write_restart ( psi, io_unit = 99 )
        CLOSE ( UNIT = 99 )
     END IF
     
     IF ( tlastloop ) THEN
        EXIT
     END IF
     
  END DO main_loop
  
END SUBROUTINE dft_scf_loop

!******************************************************************************

END MODULE dft_structure_control

!******************************************************************************
