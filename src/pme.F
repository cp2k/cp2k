!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      JGH (21-Mar-2001) : Complete rewrite
!> \author CJM and APSI
! *****************************************************************************
MODULE pme

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE atprop_types,                    ONLY: atprop_type
  USE bibliography,                    ONLY: cite_reference,&
                                             darden1993
  USE cell_types,                      ONLY: cell_type
  USE dg_rho0_types,                   ONLY: dg_rho0_type
  USE dg_types,                        ONLY: dg_get,&
                                             dg_type
  USE dgs,                             ONLY: dg_get_patch,&
                                             dg_get_strucfac,&
                                             dg_sum_patch,&
                                             dg_sum_patch_force_1d,&
                                             dg_sum_patch_force_3d
  USE ewald_environment_types,         ONLY: ewald_env_get,&
                                             ewald_environment_type
  USE ewald_pw_types,                  ONLY: ewald_pw_get,&
                                             ewald_pw_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: fourpi
  USE particle_types,                  ONLY: particle_type
  USE pme_tools,                       ONLY: get_center,&
                                             set_list
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_methods,                      ONLY: pw_copy,&
                                             pw_derive,&
                                             pw_integral_a2b,&
                                             pw_transfer
  USE pw_poisson_methods,              ONLY: pw_poisson_solve
  USE pw_poisson_types,                ONLY: pw_poisson_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type,&
                                             pw_type
  USE realspace_grid_types,            ONLY: &
       pw2rs, realspace_grid_desc_type, realspace_grid_p_type, &
       realspace_grid_type, rs2pw, rs_grid_create, rs_grid_release, &
       rs_grid_set_box, rs_grid_zero, rs_pw_transfer
  USE shell_potential_types,           ONLY: shell_kind_type
  USE structure_factor_types,          ONLY: structure_factor_type
  USE structure_factors,               ONLY: structure_factor_allocate,&
                                             structure_factor_deallocate,&
                                             structure_factor_init
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: pme_evaluate
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pme'

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param ewald_env ...
!> \param ewald_pw ...
!> \param box ...
!> \param particle_set ...
!> \param vg_coulomb ...
!> \param fg_coulomb ...
!> \param pv_g ...
!> \param shell_particle_set ...
!> \param core_particle_set ...
!> \param fgshell_coulomb ...
!> \param fgcore_coulomb ...
!> \param use_virial ...
!> \param charges ...
!> \param atprop ...
!> \par History
!>      JGH (15-Mar-2001) : New electrostatic calculation and pressure tensor
!>      JGH (21-Mar-2001) : Complete rewrite
!>      JGH (21-Mar-2001) : Introduced real space density type for future
!>                          parallelisation
!> \author CJM and APSI
! *****************************************************************************
  SUBROUTINE pme_evaluate ( ewald_env, ewald_pw, box, particle_set, vg_coulomb, &
       fg_coulomb, pv_g, shell_particle_set, core_particle_set,&
       fgshell_coulomb, fgcore_coulomb, use_virial, charges, atprop)
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(cell_type), POINTER                 :: box
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), INTENT(OUT)               :: vg_coulomb
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT), OPTIONAL                  :: fg_coulomb, pv_g
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: shell_particle_set, &
                                                core_particle_set
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT), OPTIONAL                  :: fgshell_coulomb, &
                                                fgcore_coulomb
    LOGICAL, INTENT(IN)                      :: use_virial
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: charges
    TYPE(atprop_type), POINTER               :: atprop

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pme_evaluate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: group, handle, i, ierr, ig, &
                                                ipart, j, nd(3), npart, &
                                                nshell, p1, p2
    LOGICAL                                  :: failure, is1_core, is2_core
    REAL(KIND=dp)                            :: alpha, dvols, fat1, ffa, ffb
    REAL(KIND=dp), DIMENSION(3)              :: fat
    REAL(KIND=dp), DIMENSION(3, 3)           :: f_stress, h_stress
    TYPE(dg_rho0_type), POINTER              :: dg_rho0
    TYPE(dg_type), POINTER                   :: dg
    TYPE(pw_grid_type), POINTER              :: grid_b, grid_s
    TYPE(pw_p_type)                          :: rhos1, rhos2
    TYPE(pw_p_type), DIMENSION(3)            :: dphi_g
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_pool_type), POINTER              :: pw_big_pool, pw_small_pool
    TYPE(pw_type), POINTER                   :: phi_g, phi_r, rhob_g, rhob_r
    TYPE(realspace_grid_desc_type), POINTER  :: rs_desc
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: drpot
    TYPE(realspace_grid_type), POINTER       :: rden, rpot
    TYPE(structure_factor_type)              :: exp_igr

    CALL timeset ( routineN, handle )
    NULLIFY( poisson_env, rden, drpot)
    failure = .FALSE.
    CALL cite_reference(Darden1993)
    CALL ewald_env_get (ewald_env, alpha=alpha, group = group)
    CALL ewald_pw_get (ewald_pw, pw_big_pool=pw_big_pool, &
         pw_small_pool=pw_small_pool, rs_desc=rs_desc,&
         poisson_env=poisson_env, dg = dg )

    grid_b => pw_big_pool%pw_grid
    grid_s => pw_small_pool%pw_grid

    CALL dg_get ( dg, dg_rho0 = dg_rho0 )

    npart = SIZE ( particle_set )

    CALL structure_factor_init(exp_igr)

    IF(PRESENT(shell_particle_set)) THEN
       CPPostcondition(ASSOCIATED(shell_particle_set),cp_failure_level,routineP,failure)
       CPPostcondition(ASSOCIATED(core_particle_set),cp_failure_level,routineP,failure)
       nshell=SIZE(shell_particle_set)
       CALL structure_factor_allocate ( grid_s%bounds, npart, exp_igr, &
            allocate_centre = .TRUE., allocate_shell_e=.TRUE., allocate_shell_centre = .TRUE., nshell=nshell )

    ELSE
       CALL structure_factor_allocate ( grid_s%bounds, npart, exp_igr, &
            allocate_centre = .TRUE. )
    END IF

    CALL pw_pool_create_pw ( pw_small_pool, rhos1%pw, &
         use_data = REALDATA3D)
    CALL pw_pool_create_pw ( pw_small_pool, rhos2%pw, &
         use_data = REALDATA3D)

    CALL rs_grid_create(rden, rs_desc)
    CALL rs_grid_set_box ( grid_b, rs=rden)
    CALL rs_grid_zero ( rden )

    IF ( rden%desc%parallel .AND. rden%desc%distributed ) THEN
       CALL get_center ( particle_set, box, exp_igr%centre, grid_b%npts, 1 )
    END IF
    IF ( PRESENT(shell_particle_set) .AND. rden%desc%parallel .AND. rden%desc%distributed ) THEN
       CALL get_center ( shell_particle_set, box, exp_igr%shell_centre, grid_b%npts, 1 )
       CALL get_center ( core_particle_set, box, exp_igr%core_centre, grid_b%npts, 1 )
    END IF

    !-------------- DENSITY CALCULATION ----------------

    ipart = 0
    DO

       CALL set_list ( particle_set, npart, exp_igr%centre, p1, rden, ipart, exp_igr %core_centre)
       CALL set_list ( particle_set, npart, exp_igr%centre, p2, rden, ipart, exp_igr %core_centre)
       IF ( p1 == 0 .AND. p2 == 0 ) EXIT

       is1_core  = (particle_set(p1)%shell_index/=0)
       IF(p2 /= 0) THEN
          is2_core  = (particle_set(p2)%shell_index/=0)
       ELSE
          is2_core  = .FALSE.
       END IF

       ! calculate function on small boxes (we use double packing in FFT)
       IF(is1_core .OR. is2_core) THEN
          CALL get_patch ( dg, particle_set, exp_igr, box, p1, p2, grid_b, grid_s, &
               rhos1, rhos2 , is1_core=is1_core, is2_core=is2_core, &
               core_particle_set=core_particle_set, charges=charges)

          ! add boxes to real space grid (big box)
          IF(is1_core) THEN
             CALL dg_sum_patch ( rden, rhos1, exp_igr%core_centre ( :, particle_set(p1)%shell_index ) )
          ELSE
             CALL dg_sum_patch ( rden, rhos1, exp_igr%centre ( :, p1 ) )
          END IF
          IF ( p2 /= 0 .AND. is2_core) THEN
             CALL dg_sum_patch ( rden, rhos2, exp_igr%core_centre ( :, particle_set(p2)%shell_index ) )
          ELSE IF(p2 /= 0) THEN
             CALL dg_sum_patch ( rden, rhos2, exp_igr%centre ( :, p2 ) )
          END IF
       ELSE
          CALL get_patch ( dg, particle_set, exp_igr, box, p1, p2, grid_b, grid_s, &
               rhos1, rhos2 , charges=charges)
          ! add boxes to real space grid (big box)
          CALL dg_sum_patch ( rden, rhos1, exp_igr%centre ( :, p1 ) )
          IF ( p2 /= 0 ) CALL dg_sum_patch ( rden, rhos2, exp_igr%centre ( :, p2 ) )
       END IF

    END DO
    IF(PRESENT(shell_particle_set)) THEN
       ipart = 0
       DO
          CALL set_list ( shell_particle_set, nshell, exp_igr%shell_centre, p1, rpot, ipart )
          CALL set_list ( shell_particle_set, nshell, exp_igr%shell_centre, p2, rpot, ipart )
          IF ( p1 == 0 .AND. p2 == 0 ) EXIT
          ! calculate function on small boxes (we use double packing in FFT)
          CALL get_patch ( dg, shell_particle_set, exp_igr, box, p1, p2, grid_b, grid_s, &
               rhos1, rhos2, is1_shell=.TRUE., is2_shell=.TRUE. , charges=charges)
          ! add boxes to real space grid (big box)
          CALL dg_sum_patch ( rpot, rhos1, exp_igr%shell_centre ( :, p1 ) )
          IF ( p2 /= 0 ) CALL dg_sum_patch ( rpot, rhos2, exp_igr%shell_centre ( :, p2 ) )
       END DO
    END IF

    CALL pw_pool_create_pw ( pw_big_pool, rhob_r, use_data = REALDATA3D, in_space = REALSPACE)
    CALL rs_pw_transfer ( rden, rhob_r, rs2pw)

    !-------------- ELECTROSTATIC CALCULATION -----------

    ! allocate intermediate arrays
    DO i = 1, 3
       NULLIFY(dphi_g ( i )%pw)
       CALL pw_pool_create_pw ( pw_big_pool, dphi_g ( i )%pw, &
            use_data = COMPLEXDATA1D, &
            in_space = RECIPROCALSPACE)
    END DO
    CALL pw_pool_create_pw ( pw_big_pool, phi_r, &
         use_data = REALDATA3D, &
         in_space = REALSPACE)

    CALL pw_poisson_solve ( poisson_env, rhob_r, vg_coulomb, phi_r, dphi_g, h_stress)

    ! atomic energies
    IF (atprop%energy .OR. atprop%stress) THEN
       dvols = rhos1%pw%pw_grid%dvol
       CALL rs_grid_create (rpot, rs_desc)
       CALL rs_pw_transfer (rpot, phi_r, pw2rs)
       ipart = 0
       DO
          CALL set_list ( particle_set, npart, exp_igr%centre, p1, rden, ipart, exp_igr %core_centre)
          CALL set_list ( particle_set, npart, exp_igr%centre, p2, rden, ipart, exp_igr %core_centre)
          IF ( p1 == 0 .AND. p2 == 0 ) EXIT
          ! integrate box and potential
          CALL get_patch ( dg, particle_set, exp_igr, box, p1, p2, grid_b, grid_s, &
               rhos1, rhos2 , charges=charges)
          ! add boxes to real space grid (big box)
          CALL dg_sum_patch_force_1d ( rpot, rhos1, exp_igr%centre(:,p1), fat1 )
          IF (atprop%energy) THEN
             atprop%atener(p1) = atprop%atener(p1) + 0.5_dp*fat1*dvols
          END IF
          IF (atprop%stress) THEN
             atprop%atstress(1,1,p1) = atprop%atstress(1,1,p1) + 0.5_dp*fat1*dvols
             atprop%atstress(2,2,p1) = atprop%atstress(2,2,p1) + 0.5_dp*fat1*dvols
             atprop%atstress(3,3,p1) = atprop%atstress(3,3,p1) + 0.5_dp*fat1*dvols
          END IF
          IF ( p2 /= 0 ) THEN
             CALL dg_sum_patch_force_1d ( rpot, rhos2, exp_igr%centre(:,p2), fat1 )
             IF (atprop%energy) THEN
                atprop%atener(p2) = atprop%atener(p2) + 0.5_dp*fat1*dvols
             END IF
             IF (atprop%stress) THEN
                atprop%atstress(1,1,p2) = atprop%atstress(1,1,p2) + 0.5_dp*fat1*dvols
                atprop%atstress(2,2,p2) = atprop%atstress(2,2,p2) + 0.5_dp*fat1*dvols
                atprop%atstress(3,3,p2) = atprop%atstress(3,3,p2) + 0.5_dp*fat1*dvols
             END IF
          END IF
       END DO
       IF (atprop%stress) THEN
          CALL pw_pool_create_pw ( pw_big_pool, phi_g,  &
                                   use_data = COMPLEXDATA1D,&
                                   in_space = RECIPROCALSPACE)
          CALL pw_pool_create_pw ( pw_big_pool, rhob_g,  &
                                   use_data = COMPLEXDATA1D,&
                                   in_space = RECIPROCALSPACE)
          ffa = ( 0.5_dp / dg_rho0%zet ( 1 ) ) ** 2
          ffb = 1.0_dp / fourpi
          DO i = 1, 3
             DO ig=grid_b%first_gne0,grid_b%ngpts_cut_local
                phi_g%cc(ig) = ffb * dphi_g(i)%pw%cc(ig) * (ffa * grid_b%gsq(ig) + 1.0_dp)
                phi_g%cc(ig) = phi_g%cc(ig) * poisson_env%green_fft%influence_fn%cc(ig)
             END DO
             IF(grid_b%have_g0) phi_g%cc(1) = 0.0_dp
             DO j = 1, i
                CALL pw_copy ( phi_g, rhob_g)
                nd = 0
                nd ( j ) = 1
                CALL pw_derive ( rhob_g, nd)
                CALL pw_transfer ( rhob_g, rhob_r)
                CALL rs_pw_transfer ( rpot, rhob_r, pw2rs)

                ipart = 0
                DO
                   CALL set_list ( particle_set, npart, exp_igr%centre, p1, rden, ipart, exp_igr %core_centre)
                   CALL set_list ( particle_set, npart, exp_igr%centre, p2, rden, ipart, exp_igr %core_centre)
                   IF ( p1 == 0 .AND. p2 == 0 ) EXIT
                   ! integrate box and potential
                   CALL get_patch ( dg, particle_set, exp_igr, box, p1, p2, grid_b, grid_s, &
                        rhos1, rhos2 , charges=charges)
                   ! add boxes to real space grid (big box)
                   CALL dg_sum_patch_force_1d ( rpot, rhos1, exp_igr%centre(:,p1), fat1 )
                   atprop%atstress(i,j,p1) = atprop%atstress(i,j,p1) + fat1*dvols
                   IF (i/=j) atprop%atstress(j,i,p1) = atprop%atstress(j,i,p1) + fat1*dvols
                   IF ( p2 /= 0 ) THEN
                      CALL dg_sum_patch_force_1d ( rpot, rhos2, exp_igr%centre(:,p2), fat1 )
                      atprop%atstress(i,j,p2) = atprop%atstress(i,j,p2) + fat1*dvols
                      IF (i/=j) atprop%atstress(j,i,p2) = atprop%atstress(j,i,p2) + fat1*dvols
                   END IF
                END DO
             END DO
          END DO
          CALL pw_pool_give_back_pw ( pw_big_pool, phi_g)
          CALL pw_pool_give_back_pw ( pw_big_pool, rhob_g)
       END IF
       CALL rs_grid_release ( rpot)
    END IF

    CALL pw_pool_give_back_pw ( pw_big_pool, phi_r)

    !---------- END OF ELECTROSTATIC CALCULATION --------

    !------------- STRESS TENSOR CALCULATION ------------

    IF ((use_virial).AND.(PRESENT(pv_g))) THEN
       DO i = 1, 3
          DO j = i, 3
             f_stress ( i, j ) = pw_integral_a2b ( dphi_g ( i )%pw, dphi_g ( j )%pw)
             f_stress ( j, i ) = f_stress ( i, j )
          END DO
       END DO
       ffa = ( 1.0_dp / fourpi ) * ( 0.5_dp / dg_rho0%zet ( 1 ) ) ** 2
       f_stress = -ffa * f_stress
       pv_g = h_stress + f_stress
    END IF

    !--------END OF STRESS TENSOR CALCULATION -----------

    ALLOCATE ( drpot(1:3), STAT=ierr )
    CPPostcondition(ierr==0,cp_failure_level,routineP,failure)
    DO i = 1, 3
       CALL rs_grid_create(drpot(i)%rs_grid, rs_desc)
       CALL rs_grid_set_box ( grid_b, rs=drpot(i)%rs_grid)
       CALL pw_transfer ( dphi_g ( i )%pw, rhob_r)
       CALL pw_pool_give_back_pw ( pw_big_pool, dphi_g ( i )%pw)
       CALL rs_pw_transfer ( drpot ( i )%rs_grid, rhob_r, pw2rs)
    END DO

    CALL pw_pool_give_back_pw ( pw_big_pool, rhob_r)

    !----------------- FORCE CALCULATION ----------------

    ! initialize the forces
    IF ( PRESENT ( fg_coulomb ) ) THEN
       fg_coulomb = 0.0_dp
       dvols = rhos1%pw%pw_grid%dvol

       ipart = 0
       DO

          CALL set_list ( particle_set, npart, exp_igr%centre, p1, rden, ipart, exp_igr%core_centre )
          CALL set_list ( particle_set, npart, exp_igr%centre, p2, rden, ipart, exp_igr%core_centre )
          IF ( p1 == 0 .AND. p2 == 0 ) EXIT

          is1_core  = (particle_set(p1)%shell_index/=0)
          IF(p2 /= 0) THEN
             is2_core  = (particle_set(p2)%shell_index/=0)
          ELSE
             is2_core  = .FALSE.
          END IF

          ! calculate function on small boxes (we use double packing in FFT)

          CALL get_patch_again ( dg, particle_set, exp_igr, p1, p2, rhos1, rhos2, &
               is1_core = is1_core, is2_core=is2_core , charges=charges)

          ! sum boxes on real space grids (big box)
          IF(is1_core ) THEN
             CALL dg_sum_patch_force_3d ( drpot, rhos1, &
                  exp_igr%core_centre ( :, particle_set(p1)%shell_index ), fat )
             fgcore_coulomb ( 1, particle_set(p1)%shell_index ) = &
                  fgcore_coulomb ( 1, particle_set(p1)%shell_index ) - fat ( 1 ) * dvols
             fgcore_coulomb ( 2,particle_set( p1)%shell_index ) = &
                  fgcore_coulomb ( 2, particle_set(p1)%shell_index ) - fat ( 2 ) * dvols
             fgcore_coulomb ( 3, particle_set(p1)%shell_index ) = &
                  fgcore_coulomb ( 3, particle_set(p1)%shell_index ) - fat ( 3 ) * dvols
          ELSE
             CALL dg_sum_patch_force_3d ( drpot, rhos1, exp_igr%centre ( :, p1 ), fat )
             fg_coulomb ( 1, p1 ) = fg_coulomb ( 1, p1 ) - fat ( 1 ) * dvols
             fg_coulomb ( 2, p1 ) = fg_coulomb ( 2, p1 ) - fat ( 2 ) * dvols
             fg_coulomb ( 3, p1 ) = fg_coulomb ( 3, p1 ) - fat ( 3 ) * dvols
          END IF
          IF ( p2 /= 0 .AND. is2_core) THEN
             CALL dg_sum_patch_force_3d ( drpot, rhos1, &
                  exp_igr%core_centre ( :, particle_set(p2)%shell_index ), fat )
             fgcore_coulomb ( 1, particle_set(p2)%shell_index ) = &
                  fgcore_coulomb ( 1, particle_set(p2)%shell_index ) - fat ( 1 ) * dvols
             fgcore_coulomb ( 2,particle_set( p2)%shell_index ) = &
                  fgcore_coulomb ( 2, particle_set(p2)%shell_index ) - fat ( 2 ) * dvols
             fgcore_coulomb ( 3, particle_set(p2)%shell_index ) = &
                  fgcore_coulomb ( 3, particle_set(p2)%shell_index ) - fat ( 3 ) * dvols
          ELSEIF ( p2 /= 0 ) THEN
             CALL dg_sum_patch_force_3d ( drpot, rhos2, exp_igr%centre ( :, p2 ), fat )
             fg_coulomb ( 1, p2 ) = fg_coulomb ( 1, p2 ) - fat ( 1 ) * dvols
             fg_coulomb ( 2, p2 ) = fg_coulomb ( 2, p2 ) - fat ( 2 ) * dvols
             fg_coulomb ( 3, p2 ) = fg_coulomb ( 3, p2 ) - fat ( 3 ) * dvols
          END IF

       END DO
    ENDIF
    IF(PRESENT(fgshell_coulomb)) THEN
       fgshell_coulomb = 0.0_dp
       dvols = rhos1%pw%pw_grid%dvol

       ipart = 0
       DO
          CALL set_list ( shell_particle_set, nshell, exp_igr%shell_centre, p1, rden, ipart )
          CALL set_list ( shell_particle_set, nshell, exp_igr%shell_centre, p2, rden, ipart )
          IF ( p1 == 0 .AND. p2 == 0 ) EXIT

          ! calculate function on small boxes (we use double packing in FFT)
          CALL get_patch_again ( dg, shell_particle_set, exp_igr, p1, p2, rhos1, rhos2, &
               is1_shell=.TRUE. , is2_shell = .TRUE., charges=charges)

          ! sum boxes on real space grids (big box)
          CALL dg_sum_patch_force_3d ( drpot, rhos1, exp_igr%shell_centre ( :, p1 ), fat )
          fgshell_coulomb ( 1, p1 ) = fgshell_coulomb ( 1, p1 ) - fat ( 1 ) * dvols
          fgshell_coulomb ( 2, p1 ) = fgshell_coulomb ( 2, p1 ) - fat ( 2 ) * dvols
          fgshell_coulomb ( 3, p1 ) = fgshell_coulomb ( 3, p1 ) - fat ( 3 ) * dvols
          IF ( p2 /= 0 ) THEN
             CALL dg_sum_patch_force_3d ( drpot, rhos2, exp_igr%shell_centre ( :, p2 ), fat )
             fgshell_coulomb ( 1, p2 ) = fgshell_coulomb ( 1, p2 ) - fat ( 1 ) * dvols
             fgshell_coulomb ( 2, p2 ) = fgshell_coulomb ( 2, p2 ) - fat ( 2 ) * dvols
             fgshell_coulomb ( 3, p2 ) = fgshell_coulomb ( 3, p2 ) - fat ( 3 ) * dvols
          END IF
       END DO

    END IF
    !--------------END OF FORCE CALCULATION -------------

    !------------------CLEANING UP ----------------------

    CALL rs_grid_release(rden)
    IF (ASSOCIATED(drpot)) THEN
       DO i = 1, 3
          CALL rs_grid_release(drpot(i)%rs_grid)
       END DO
       DEALLOCATE ( drpot, STAT = ierr )
       CPPostcondition(ierr==0,cp_failure_level,routineP,failure)
    END IF

    CALL pw_pool_give_back_pw ( pw_small_pool, rhos1%pw)
    CALL pw_pool_give_back_pw ( pw_small_pool, rhos2%pw)
    CALL structure_factor_deallocate ( exp_igr )

    CALL timestop(handle)

  END SUBROUTINE pme_evaluate

! *****************************************************************************
!> \brief Calculates local density in a small box
!> \param dg ...
!> \param particle_set ...
!> \param exp_igr ...
!> \param box ...
!> \param p1 ...
!> \param p2 ...
!> \param grid_b ...
!> \param grid_s ...
!> \param rhos1 ...
!> \param rhos2 ...
!> \param is1_core ...
!> \param is2_core ...
!> \param is1_shell ...
!> \param is2_shell ...
!> \param core_particle_set ...
!> \param charges ...
!> \par History
!>      JGH (23-Mar-2001) : Switch to integer from particle list pointers
!> \author JGH (21-Mar-2001)
! *****************************************************************************
  SUBROUTINE get_patch ( dg, particle_set, exp_igr, box, p1, p2, &
       grid_b, grid_s, rhos1, rhos2 , is1_core, is2_core, is1_shell,&
       is2_shell, core_particle_set, charges)

    TYPE(dg_type), POINTER                   :: dg
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(structure_factor_type)              :: exp_igr
    TYPE(cell_type), INTENT(IN)              :: box
    INTEGER, INTENT(IN)                      :: p1, p2
    TYPE(pw_grid_type), INTENT(IN)           :: grid_b, grid_s
    TYPE(pw_p_type)                          :: rhos1, rhos2
    LOGICAL, OPTIONAL                        :: is1_core, is2_core, &
                                                is1_shell, is2_shell
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: core_particle_set
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: charges

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_patch', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: ex1, ex2, ey1, ey2, ez1, ez2
    INTEGER, DIMENSION(:), POINTER           :: center1, center2
    LOGICAL                                  :: failure, my_is1_core, &
                                                my_is1_shell, my_is2_core, &
                                                my_is2_shell, use_charge_array
    REAL(KIND=dp)                            :: q1, q2
    REAL(KIND=dp), DIMENSION(3)              :: r1, r2
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(dg_rho0_type), POINTER              :: dg_rho0
    TYPE(pw_p_type), POINTER                 :: rho0
    TYPE(shell_kind_type), POINTER           :: shell

    failure = .FALSE.
    NULLIFY(shell)
    use_charge_array=PRESENT(charges)
    IF (use_charge_array) use_charge_array=ASSOCIATED(charges)
    my_is1_core = .FALSE.
    my_is2_core = .FALSE.
    IF(PRESENT(is1_core)) my_is1_core = is1_core
    IF(PRESENT(is2_core)) my_is2_core = is2_core
    IF(my_is1_core .OR. my_is2_core) THEN
       CPPostcondition(PRESENT(core_particle_set),cp_failure_level,routineP,failure)
    END IF
    my_is1_shell = .FALSE.
    my_is2_shell = .FALSE.
    IF(PRESENT(is1_shell)) my_is1_shell = is1_shell
    IF(PRESENT(is2_shell)) my_is2_shell = is2_shell
    IF(my_is1_core .AND. my_is1_shell) THEN
       CALL stop_program(routineN,moduleN,__LINE__,&
                         "Shell-model: cannot be core and shell simultaneously")
    END IF

    CALL dg_get ( dg, dg_rho0 = dg_rho0 )
    rho0 => dg_rho0%density

    IF(my_is1_core) THEN
       r1 = core_particle_set (particle_set(p1)%shell_index )%r
    ELSE
       r1 = particle_set ( p1 )%r
    END IF
    atomic_kind => particle_set (p1)%atomic_kind
    IF(my_is1_core) THEN
       CALL get_atomic_kind ( atomic_kind=atomic_kind, shell=shell)
       q1 = shell%charge_core
    ELSE IF (my_is1_shell) THEN
       CALL get_atomic_kind ( atomic_kind=atomic_kind, shell=shell)
       q1 = shell%charge_shell
    ELSE
       CALL get_atomic_kind ( atomic_kind=atomic_kind, qeff=q1)
    END IF
    IF (use_charge_array) q1=charges(p1)

    IF(my_is1_shell) THEN
       center1 => exp_igr%shell_centre ( :, p1 )
       ex1 => exp_igr%shell_ex ( :, p1 )
       ey1 => exp_igr%shell_ey ( :, p1 )
       ez1 => exp_igr%shell_ez ( :, p1 )
    ELSEIF(my_is1_core) THEN
       center1 => exp_igr%core_centre ( :, particle_set(p1)%shell_index )
       ex1 => exp_igr%core_ex ( :, particle_set(p1)%shell_index )
       ey1 => exp_igr%core_ey ( :, particle_set(p1)%shell_index )
       ez1 => exp_igr%core_ez ( :, particle_set(p1)%shell_index )
    ELSE
       center1 => exp_igr%centre ( :, p1 )
       ex1 => exp_igr%ex ( :, p1 )
       ey1 => exp_igr%ey ( :, p1 )
       ez1 => exp_igr%ez ( :, p1 )
    END IF

    CALL dg_get_strucfac ( box%hmat, r1, grid_s%npts, grid_b%npts, center1, &
         exp_igr%lb, ex1, ey1, ez1 )

    IF ( p2 /= 0 ) THEN
       IF(my_is2_core) THEN
          r2 = core_particle_set ( particle_set(p2)%shell_index )%r
       ELSE
          r2 = particle_set ( p2 )%r
       END IF
       atomic_kind => particle_set (p2)%atomic_kind
       IF(my_is2_core) THEN
          CALL get_atomic_kind ( atomic_kind=atomic_kind, shell=shell)
          q2 = shell%charge_core
       ELSE IF (my_is2_shell) THEN
          CALL get_atomic_kind ( atomic_kind=atomic_kind, shell=shell)
          q2 = shell%charge_shell
       ELSE
          CALL get_atomic_kind ( atomic_kind=atomic_kind, qeff=q2)
       END IF
       IF (use_charge_array) q2=charges(p2)

       IF(my_is2_shell) THEN
          center2 => exp_igr%shell_centre ( :, p2 )
          ex2 => exp_igr%shell_ex ( :, p2 )
          ey2 => exp_igr%shell_ey ( :, p2 )
          ez2 => exp_igr%shell_ez ( :, p2 )
       ELSEIF(my_is2_core) THEN
          center2 => exp_igr%core_centre ( :, particle_set(p2)%shell_index )
          ex2 => exp_igr%core_ex ( :, particle_set(p2)%shell_index )
          ey2 => exp_igr%core_ey ( :, particle_set(p2)%shell_index )
          ez2 => exp_igr%core_ez ( :, particle_set(p2)%shell_index )
       ELSE
          center2 => exp_igr%centre ( :, p2 )
          ex2 => exp_igr%ex ( :, p2 )
          ey2 => exp_igr%ey ( :, p2 )
          ez2 => exp_igr%ez ( :, p2 )
       END IF
       CALL dg_get_strucfac ( box%hmat, r2, grid_s%npts, grid_b%npts, center2, &
            exp_igr%lb, ex2, ey2, ez2 )
    END IF

    IF ( p2 == 0 ) THEN
       CALL dg_get_patch ( rho0, rhos1, q1, ex1, ey1, ez1 )
    ELSE
       CALL dg_get_patch ( rho0, rhos1, rhos2, q1, q2, ex1, ey1, ez1, ex2, ey2, ez2 )
    END IF

  END SUBROUTINE get_patch

! *****************************************************************************
!> \brief ...
!> \param dg ...
!> \param particle_set ...
!> \param exp_igr ...
!> \param p1 ...
!> \param p2 ...
!> \param rhos1 ...
!> \param rhos2 ...
!> \param is1_core ...
!> \param is2_core ...
!> \param is1_shell ...
!> \param is2_shell ...
!> \param charges ...
! *****************************************************************************
  SUBROUTINE get_patch_again ( dg, particle_set, exp_igr, p1, p2, rhos1, rhos2, is1_core, &
       is2_core, is1_shell, is2_shell, charges )

    TYPE(dg_type), POINTER                   :: dg
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(structure_factor_type)              :: exp_igr
    INTEGER, INTENT(IN)                      :: p1, p2
    TYPE(pw_p_type)                          :: rhos1, rhos2
    LOGICAL, OPTIONAL                        :: is1_core, is2_core, &
                                                is1_shell, is2_shell
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: charges

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_patch_again', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: ex1, ex2, ey1, ey2, ez1, ez2
    LOGICAL                                  :: my_is1_core, my_is1_shell, &
                                                my_is2_core, my_is2_shell, &
                                                use_charge_array
    REAL(KIND=dp)                            :: q1, q2
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(dg_rho0_type), POINTER              :: dg_rho0
    TYPE(pw_p_type), POINTER                 :: rho0
    TYPE(shell_kind_type), POINTER           :: shell

    NULLIFY(shell)
    use_charge_array=PRESENT(charges)
    IF (use_charge_array) use_charge_array=ASSOCIATED(charges)
    my_is1_core = .FALSE.
    my_is2_core = .FALSE.
    IF(PRESENT(is1_core)) my_is1_core = is1_core
    IF(PRESENT(is2_core)) my_is2_core = is2_core
    my_is1_shell = .FALSE.
    my_is2_shell = .FALSE.
    IF(PRESENT(is1_shell)) my_is1_shell = is1_shell
    IF(PRESENT(is2_shell)) my_is2_shell = is2_shell

    CALL dg_get ( dg, dg_rho0 = dg_rho0 )
    rho0 => dg_rho0%density

    atomic_kind => particle_set (p1)%atomic_kind
    IF(my_is1_core) THEN
       CALL get_atomic_kind ( atomic_kind=atomic_kind, shell=shell)
       q1 = shell%charge_core
    ELSE IF (my_is1_shell) THEN
       CALL get_atomic_kind ( atomic_kind=atomic_kind, shell=shell)
       q1 = shell%charge_shell
    ELSE
       CALL get_atomic_kind ( atomic_kind=atomic_kind, qeff=q1)
    END IF
    IF (use_charge_array) q1=charges(p1)
    IF (my_is1_core) THEN
       ex1 => exp_igr%core_ex ( :, particle_set (p1)%shell_index )
       ey1 => exp_igr%core_ey ( :, particle_set (p1)%shell_index )
       ez1 => exp_igr%core_ez ( :, particle_set (p1)%shell_index )
    ELSEIF(my_is1_shell)THEN
       ex1 => exp_igr%shell_ex ( :, p1 )
       ey1 => exp_igr%shell_ey ( :, p1 )
       ez1 => exp_igr%shell_ez ( :, p1 )
    ELSE
       ex1 => exp_igr%ex ( :, p1 )
       ey1 => exp_igr%ey ( :, p1 )
       ez1 => exp_igr%ez ( :, p1 )
    END IF

    IF ( p2 /= 0 ) THEN
       atomic_kind => particle_set (p2)%atomic_kind
       IF(my_is2_core) THEN
          CALL get_atomic_kind ( atomic_kind=atomic_kind, shell=shell)
          q2 = shell%charge_core
       ELSE IF (my_is2_shell) THEN
          CALL get_atomic_kind ( atomic_kind=atomic_kind, shell=shell)
          q2 = shell%charge_shell
       ELSE
          CALL get_atomic_kind ( atomic_kind=atomic_kind, qeff=q2)
       END IF
       IF (use_charge_array) q2=charges(p2)
       IF (my_is2_core) THEN
          ex2 => exp_igr%core_ex ( :, particle_set (p2)%shell_index )
          ey2 => exp_igr%core_ey ( :, particle_set (p2)%shell_index )
          ez2 => exp_igr%core_ez ( :, particle_set (p2)%shell_index )
       ELSEIF (my_is2_shell) THEN
          ex2 => exp_igr%shell_ex ( :, p2 )
          ey2 => exp_igr%shell_ey ( :, p2 )
          ez2 => exp_igr%shell_ez ( :, p2 )
       ELSE
          ex2 => exp_igr%ex ( :, p2 )
          ey2 => exp_igr%ey ( :, p2 )
          ez2 => exp_igr%ez ( :, p2 )
       END IF
    END IF

    IF ( p2 == 0 ) THEN
       CALL dg_get_patch ( rho0, rhos1, q1, ex1, ey1, ez1 )
    ELSE
       CALL dg_get_patch ( rho0, rhos1, rhos2, q1, q2, &
            ex1, ey1, ez1, ex2, ey2, ez2 )
    END IF

  END SUBROUTINE get_patch_again

END MODULE pme

