!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****m* cp2k/xc_pot_saop *
!!
!!   NAME
!!     xc_pot_saop
!!
!!   FUNCTION
!!     Calculate the saop potential
!!
!!****

MODULE xc_pot_saop
  
  USE coefficient_types,               ONLY: coeff_deallocate,&
                                             coeff_init,&
                                             coeff_type
  USE cp_control_types,                ONLY: dft_control_type,&
                                             xc_control_type
  USE cp_control_utils,                ONLY: xc_get_functional
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_fm_types,                     ONLY: cp_fm_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_subsystem_types,              ONLY: cp_subsystem_p_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_give_back_coeff,&
                                             pw_pool_init_coeff,&
                                             pw_pool_p_type,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type,&
                                             pw_release,&
                                             pw_zero
  USE qs_collocate_density,            ONLY: calculate_wavefunction
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_integrate_potential,          ONLY: integrate_v_rspace
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type

  USE qs_rho_types,                    ONLY: qs_rho_type
  USE simulation_cell,                 ONLY: cell_type
  USE sparse_matrix_types,             ONLY: real_matrix_p_type
  USE xc,                              ONLY: xc_vxc_pw_create
  USE xc_rho_set_types,                ONLY: xc_rho_cflags_type,&
                                             xc_rho_cflags_setall,&
                                             xc_rho_set_type,&
                                             xc_rho_set_create,&
                                             xc_rho_set_release,&
                                             xc_rho_set_update
  USE xc_xbecke88,                     ONLY: xb88_lda_info,&
                                             xb88_lsd_info
  IMPLICIT NONE

  PUBLIC :: add_saop_pot

  CHARACTER(LEN=*), PARAMETER :: moduleN="xc_pot_saop"

  ! should be eliminated
  REAL(KIND=dp), PARAMETER :: alpha=1.19_dp, beta=0.01_dp, K_rho=0.42

CONTAINS

  !*****************************************************************

  SUBROUTINE add_saop_pot(ks_matrix, qs_env, error)

    ! ARGUMENTS
    TYPE(real_matrix_p_type), &
         DIMENSION(:), POINTER               :: ks_matrix
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), &
         INTENT(INOUT), OPTIONAL             :: error

    ! PARAMETERS
    CHARACTER(LEN=*), PARAMETER              :: routineN="add_saop_pot", &
                                                routineP=moduleN//"/"//routineN

    ! LOCALS
    INTEGER                                  :: ispin, nspins, n_gridpoints
    INTEGER                                  :: output_unit, stat
    INTEGER                                  :: homo, nmo, orb
    INTEGER                                  :: i,j,k
    INTEGER, DIMENSION(2, 3)                 :: bo

    LOGICAL                                  :: failure, ionode, lsd

    REAL(KIND=dp)                            :: xc_energy, efac
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mo_eigenvalues
    REAL(KIND = dp), DIMENSION(:,:,:), POINTER :: rho_to_divide

    TYPE(cell_type), POINTER                 :: cell
    TYPE(coeff_type), DIMENSION(:), POINTER  :: vxc_SAOP
    TYPE(coeff_type)                         :: orbital, orbital_g
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsystem_p_type), DIMENSION(1)  :: subsys
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: molecular_orbitals
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: vxc_GLLB, vxc_LB, vxc_tmp, vxc_tau
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_g, rho_r, tau
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_rho_type), POINTER               :: rho_struct
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_rho_cflags_type)                 :: needs

    TYPE(xc_control_type), POINTER           :: xc_control
    TYPE(xc_control_type)                    :: my_xc_control
    TYPE(dft_control_type), POINTER          :: dft_control

    failure = .FALSE.

    NULLIFY(pw_env, dft_control, xc_control, auxbas_pw_pool)
    NULLIFY(rho_g, rho_r, tau, rho_struct, rho_set)
    NULLIFY(vxc_GLLB, vxc_LB, vxc_tmp, vxc_SAOP, vxc_tau)
    NULLIFY(mo_eigenvalues, rho_to_divide)
    logger => cp_error_get_logger(error)
    ionode = (logger%para_env%mepos==logger%para_env%source)
    output_unit = cp_logger_get_default_unit_nr(logger)

    if (ionode) write (output_unit, *) "in "//routineP

    CALL get_qs_env(qs_env=qs_env,&
                    rho=rho_struct,&
                    pw_env=pw_env,&
                    subsys=subsys(1)%subsys,&
                    cell=cell,&
                    dft_control=dft_control,&
                    mos=molecular_orbitals)
    xc_control=>dft_control%xc_control
    CALL pw_env_get(pw_env,auxbas_pw_pool=auxbas_pw_pool,error=error)

    nspins = dft_control%nspins
    lsd = (nspins>1)

    ALLOCATE(rho_r(nspins), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO ispin=1,nspins
       rho_r(ispin)%pw => rho_struct%rho_r(ispin)%pw
    END DO

    my_xc_control%functionals              = xc_control%functionals
    my_xc_control%gradient_functionals     = xc_control%gradient_functionals
    my_xc_control%crossterms               = xc_control%crossterms

    !----------------------------------------!
    ! Construct the LB94 potential in vxc_LB !
    !----------------------------------------!
    xc_control%functionals = ""
    CALL xc_get_functional ( 1, "SLATER", xc_control%functionals, &
                                          xc_control%gradient_functionals, &
                                          xc_control%crossterms )
    DO i=1, 3
       IF (xc_control%functionals(i) == "") xc_control%functionals(i) = "NONE"
    END DO
    CALL xc_vxc_pw_create(vxc_LB, vxc_tau, xc_energy, rho_r, rho_g, tau,&
                          dft_control, cell, auxbas_pw_pool, subsys, &
                          .FALSE., error)
    n_gridpoints = SIZE(vxc_LB(1)%pw%cr3d,1) * &
                   SIZE(vxc_LB(1)%pw%cr3d,2) * &
                   SIZE(vxc_LB(1)%pw%cr3d,3)
    DO ispin=1, nspins
       CALL dscal(n_gridpoints, alpha, vxc_LB(ispin)%pw%cr3d, 1)
    END DO
    xc_control%functionals = ""
    CALL xc_get_functional ( 2, "VWN", xc_control%functionals, &
                                       xc_control%gradient_functionals, &
                                       xc_control%crossterms )
    DO i=1, 3
       IF (xc_control%functionals(i) == "") xc_control%functionals(i) = "NONE"
    END DO
    CALL xc_vxc_pw_create(vxc_tmp, vxc_tau, xc_energy, rho_r, rho_g, tau,&
                          dft_control, cell, auxbas_pw_pool, subsys, &
                          .FALSE., error)
    DO ispin=1, nspins
       CALL daxpy(n_gridpoints, 1.0_dp, vxc_tmp(ispin)%pw%cr3d, 1, &
                                        vxc_LB(ispin)%pw%cr3d, 1)
    END DO

    !--------------------------------------!
    ! create the density needed to calc lb !
    !--------------------------------------!
    bo = rho_r(1)%pw%pw_grid%bounds_local
    CALL xc_rho_set_create(rho_set, bo, &
                           xc_control%density_cut, &
                           xc_control%gradient_cut, &
                           xc_control%tau_cut, &
                           error=error)
    CALL xc_rho_cflags_setall(needs,.FALSE.,error)
    IF (lsd) THEN
       CALL xb88_lsd_info(needs=needs)
    ELSE
       CALL xb88_lda_info(needs=needs)
    END IF
    CALL xc_rho_set_update(rho_set, rho_r, rho_g, tau, needs, &
                           dft_control, cell, auxbas_pw_pool, subsys, error)
 !!   CALL add_lb_pot(vxc_LB, rho_set, lsd, error)

    !---------------------------------------------!
    ! Construct 2 times PBE potential in vxc_GLLB !
    !---------------------------------------------!
    xc_control%functionals = ""
    CALL xc_get_functional ( 0, "PBE", xc_control%functionals, &
                                       xc_control%gradient_functionals, &
                                       xc_control%crossterms )
    DO i=1, 3
       IF (xc_control%functionals(i) == "") xc_control%functionals(i) = "NONE"
    END DO

    CALL xc_vxc_pw_create(vxc_GLLB, vxc_tau, xc_energy, rho_r, rho_g, tau,&
                          dft_control, cell, auxbas_pw_pool, subsys, &
                          .FALSE., error)
    DO ispin=1, nspins
       CALL dscal(n_gridpoints, 2.0_dp, vxc_GLLB(ispin)%pw%cr3d, 1)
    END DO


    CALL pw_pool_init_coeff(auxbas_pw_pool,orbital,&
                            use_data = REALDATA3D,&
                            in_space = REALSPACE, error=error)
    CALL pw_pool_init_coeff(auxbas_pw_pool,orbital_g,&
                            use_data = COMPLEXDATA1D,&
                            in_space = RECIPROCALSPACE, error=error)

    DO ispin=1, nspins
       CALL get_mo_set(molecular_orbitals(ispin)%mo_set,&
                       mo_coeff=mo_coeff, &
                       eigenvalues=mo_eigenvalues,&
                       homo=homo,nmo=nmo)
       CALL pw_zero(vxc_tmp(ispin)%pw)
       DO orb=1, homo-1
          CALL pw_zero(orbital%pw)
          CALL pw_zero(orbital_g%pw)
          CALL calculate_wavefunction(mo_coeff,orb,orbital,orbital_g,qs_env,error)
          efac = K_rho*SQRT(mo_eigenvalues(homo)-mo_eigenvalues(orb))
          DO k=bo(1,3), bo(2,3)
             DO j=bo(1,2), bo(2,2)
                DO i=bo(1,1), bo(2,1)
                   vxc_tmp(ispin)%pw%cr3d(i,j,k) = vxc_tmp(ispin)%pw%cr3d(i,j,k) + &
                        efac * orbital%pw%cr3d(i,j,k)*orbital%pw%cr3d(i,j,k)
                END DO
             END DO
          END DO
       END DO
       IF (.NOT.lsd) THEN
          rho_to_divide => rho_set%rho
       ELSE
          IF (ispin==1) THEN
             rho_to_divide => rho_set%rhoa
          ELSE
             rho_to_divide => rho_set%rhob
          END IF
       END IF
       DO k=bo(1,3), bo(2,3)
          DO j=bo(1,2), bo(2,2)
             DO i=bo(1,1), bo(2,1)
                vxc_tmp(ispin)%pw%cr3d(i,j,k) = vxc_tmp(ispin)%pw%cr3d(i,j,k) / &
                                                rho_to_divide(i,j,k)
             END DO
          END DO
       END DO
       NULLIFY(rho_to_divide)
       CALL daxpy(n_gridpoints, 1.0_dp, vxc_tmp(ispin)%pw%cr3d, 1, &
                                        vxc_GLLB(ispin)%pw%cr3d, 1)
    END DO
    
       
    !---------------!
    ! Assemble SAOP !
    !---------------!
    ALLOCATE(vxc_SAOP(nspins), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO ispin=1, nspins
       CALL pw_pool_init_coeff(auxbas_pw_pool, vxc_SAOP(ispin), &
                               use_data=REALDATA3D, in_space=REALSPACE, &
                               error=error)
       CALL pw_zero(vxc_SAOP(ispin)%pw)
    END DO
    DO ispin=1, nspins
!!       DO orb=1, homo
!!          CALL pw_zero(vxc_tmp(ispin)%pw)
!!          efac = EXP(-2.0_dp*(mo_eigenvalues(homo)-mo_eigenvalues(orb))**2)
!!
!!          CALL daxpy(n_gridpoints, 1.0_dp-efac, vxc_GLLB(ispin)%pw%cr3d, 1, &
!!                                                vxc_tmp(ispin)%pw%cr3d, 1)
!!          CALL daxpy(n_gridpoints, efac, vxc_LB(ispin)%pw%cr3d, 1, &
!!                                         vxc_tmp(ispin)%pw%cr3d, 1)
!!          
!!          CALL pw_zero(orbital%pw)
!!          CALL pw_zero(orbital_g%pw)
!!          CALL calculate_wavefunction(mo_coeff,orb,orbital,orbital_g,qs_env,error)
!!          DO k=bo(1,3), bo(2,3)
!!             DO j=bo(1,2), bo(2,2)
!!                DO i=bo(1,1), bo(2,1)
!!                   vxc_SAOP(ispin)%pw%cr3d(i,j,k) = vxc_SAOP(ispin)%pw%cr3d(i,j,k) + &
!!                        orbital%pw%cr3d(i,j,k)*orbital%pw%cr3d(i,j,k)*&
!!                        vxc_tmp(ispin)%pw%cr3d(i,j,k)
!!                END DO
!!             END DO
!!          END DO
!!       END DO
!!       IF (.NOT.lsd) THEN
!!          rho_to_divide => rho_set%rho
!!       ELSE
!!          IF (ispin==1) THEN
!!             rho_to_divide => rho_set%rhoa
!!          ELSE
!!             rho_to_divide => rho_set%rhob
!!          END IF
!!       END IF
!!       DO k=bo(1,3), bo(2,3)
!!          DO j=bo(1,2), bo(2,2)
!!             DO i=bo(1,1), bo(2,1)
!!                vxc_SAOP(ispin)%pw%cr3d(i,j,k) = vxc_SAOP(ispin)%pw%cr3d(i,j,k) / &
!!                                                rho_to_divide(i,j,k)
!!             END DO
!!          END DO
!!       END DO
!!       NULLIFY(rho_to_divide)
       DO k=bo(1,3), bo(2,3)
          DO j=bo(1,2), bo(2,2)
             DO i=bo(1,1), bo(2,1)
                vxc_SAOP(ispin)%pw%cr3d(i,j,k) = vxc_LB(ispin)%pw%cr3d(i,j,k)
             END DO
          END DO
       END DO
       NULLIFY(rho_to_divide)
    END DO
       
    CALL xc_rho_set_release(rho_set, auxbas_pw_pool, error)
    CALL pw_pool_give_back_coeff(auxbas_pw_pool,orbital,&
         error=error)
    CALL pw_pool_give_back_coeff(auxbas_pw_pool,orbital_g,&
         error=error)
    DO ispin=1, nspins
       CALL pw_release(vxc_GLLB(ispin)%pw)
       CALL pw_release(vxc_LB(ispin)%pw)
       CALL pw_release(vxc_tmp(ispin)%pw)
    END DO
    DEALLOCATE(vxc_GLLB, vxc_LB, vxc_tmp, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    

    !--------------------!
    ! Do the integration !
    !--------------------!
    DO ispin=1, nspins
       CALL dscal(n_gridpoints, vxc_SAOP(ispin)%pw%pw_grid%dvol, &
                  vxc_SAOP(ispin)%pw%cr3d, 1)

       CALL integrate_v_rspace(vxc_SAOP(ispin), rho_struct%rho_ao(ispin), &
                               ks_matrix(ispin), qs_env, &
                               calculate_forces=.FALSE.,&
                               gapw=.FALSE.,error=error)

       CALL pw_pool_give_back_coeff(auxbas_pw_pool, vxc_SAOP(ispin), error=error)
    END DO
    
    DEALLOCATE(rho_r, vxc_SAOP, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    !-------------------------!
    ! restore original values !
    !-------------------------!
    xc_control%functionals = my_xc_control%functionals
    xc_control%gradient_functionals = my_xc_control%gradient_functionals
    xc_control%crossterms = my_xc_control%crossterms

  END SUBROUTINE add_saop_pot

  !*****************************************************************

  SUBROUTINE add_lb_pot(pot, rho_set, lsd, error)

    TYPE(pw_p_type), DIMENSION(:), POINTER   :: pot
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    LOGICAL, INTENT(IN)                      :: lsd
    TYPE(cp_error_type), &
         INTENT(INOUT), OPTIONAL             :: error

    INTEGER                                  :: i,j,k
    INTEGER, DIMENSION(2,3)                  :: bo
    
    REAL(KIND=dp)                            :: n_13, x, x2

    IF (.NOT.lsd) THEN
       bo(1,:) = LBOUND(rho_set%rho)
       bo(2,:) = UBOUND(rho_set%rho)
    ELSE
       bo(1,:) = LBOUND(rho_set%rhoa)
       bo(2,:) = UBOUND(rho_set%rhoa)
    END IF

    DO k=bo(1,3), bo(2,3)
       DO j=bo(1,2), bo(2,2)
          DO i=bo(1,1), bo(2,1)
             IF (.NOT.lsd) THEN
                n_13 = rho_set%rho_1_3(i,j,k)
                x    = rho_set%norm_drho(i,j,k)/(rho_set%rho(i,j,k)*n_13)
                x2   = x*x
                pot(1)%pw%cr3d(i,j,k) = pot(1)%pw%cr3d(i,j,k) - &
                     beta*x2*n_13/(1.0_dp+3.0_dp*beta*x*LOG(SQRT(x2+1.0_dp)+x))
             ELSE
                n_13 = rho_set%rhoa_1_3(i,j,k)
                x    = rho_set%norm_drhoa(i,j,k)/(rho_set%rhoa(i,j,k)*n_13)
                x2   = x*x
                pot(1)%pw%cr3d(i,j,k) = pot(1)%pw%cr3d(i,j,k) - &
                     beta*x2*n_13/(1.0_dp+3.0_dp*beta*x*LOG(SQRT(x2+1.0_dp)+x))
                n_13 = rho_set%rhob_1_3(i,j,k)
                x    = rho_set%norm_drhob(i,j,k)/(rho_set%rhob(i,j,k)*n_13)
                x2   = x*x
                pot(2)%pw%cr3d(i,j,k) = pot(2)%pw%cr3d(i,j,k) - &
                     beta*x2*n_13/(1.0_dp+3.0_dp*beta*x*LOG(SQRT(x2+1.0_dp)+x))
             END IF
          END DO
       END DO
    END DO

  END SUBROUTINE add_lb_pot

  !*****************************************************************

END MODULE xc_pot_saop
