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
  
  USE coefficient_types,               ONLY: coeff_type
  USE cp_control_types,                ONLY: dft_control_type,&
                                             xc_control_type
  USE cp_control_utils,                ONLY: xc_get_functional
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_get_submatrix,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_set_submatrix,&
                                             cp_fm_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_sm_fm_interactions,           ONLY: cp_sm_plus_fm_fm_t
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_coeff,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_init_coeff,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type,&
                                             pw_release,&
                                             pw_zero
  USE qs_collocate_density,            ONLY: calculate_rho_elec,&
                                             calculate_wavefunction
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_integrate_potential,          ONLY: integrate_v_rspace
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE simulation_cell,                 ONLY: cell_type
  USE sparse_matrix_types,             ONLY: deallocate_matrix,&
                                             real_matrix_p_type,&
                                             replicate_matrix,&
                                             set_matrix
  USE xc,                              ONLY: xc_vxc_pw_create
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                             xc_dset_create,&
                                             xc_dset_get_derivative,&
                                             xc_dset_release
  USE xc_derivative_types,             ONLY: xc_derivative_get,&
                                             xc_derivative_type
  USE xc_derivatives,                  ONLY: xc_functional_eval
  USE xc_rho_set_types,                ONLY: xc_rho_cflags_setall,&
                                             xc_rho_cflags_type,&
                                             xc_rho_set_create,&
                                             xc_rho_set_release,&
                                             xc_rho_set_type,&
                                             xc_rho_set_update
  USE xc_xbecke88,                     ONLY: xb88_lda_info,&
                                             xb88_lsd_info

  IMPLICIT NONE

  PUBLIC :: add_saop_pot

  CHARACTER(LEN=*), PARAMETER :: moduleN="xc_pot_saop"

  ! should be eliminated
  REAL(KIND=dp), PARAMETER :: alpha=1.19_dp, beta=0.01_dp, K_rho=0.42_dp
  REAL(KIND=dp), PARAMETER :: kappa = 0.804_dp, mu = 0.21951_dp, &
                              beta_ec  = 0.066725_dp, gamma_saop = 0.031091_dp


CONTAINS

  !*****************************************************************

  SUBROUTINE add_saop_pot(ks_matrix, qs_env, error)

    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "add_saop_pot", &
      routineP = moduleN//"/"//routineN

    INTEGER                                  :: homo, i, ispin, j, k, nspins, &
                                                orb, output_unit, stat
    INTEGER, DIMENSION(2)                    :: ncol, nrow
    INTEGER, DIMENSION(2, 3)                 :: bo
    LOGICAL                                  :: failure, ionode, lsd
    REAL(KIND=dp)                            :: efac, tot_rho_psi, we_GLLB, &
                                                we_LB, xc_energy
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: coeff_col
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mo_eigenvalues
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: e_uniform
    TYPE(cell_type), POINTER                 :: cell
    TYPE(coeff_type)                         :: orbital, orbital_g
    TYPE(coeff_type), DIMENSION(:), POINTER  :: vxc_SAOP
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: single_mo_coeff
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: molecular_orbitals
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_g, rho_r, tau, vxc_GLLB, &
                                                vxc_LB, vxc_tau, vxc_tmp
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_rho_type), POINTER               :: rho_struct
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: orb_dens
    TYPE(xc_control_type)                    :: xc_control
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    TYPE(xc_derivative_type), POINTER        :: deriv
    TYPE(xc_rho_cflags_type)                 :: needs
    TYPE(xc_rho_set_type), POINTER           :: rho_set

    failure = .FALSE.

    NULLIFY(pw_env, dft_control, auxbas_pw_pool)
    NULLIFY(rho_g, rho_r, tau, rho_struct, rho_set, e_uniform)
    NULLIFY(vxc_GLLB, vxc_LB, vxc_tmp, vxc_SAOP, vxc_tau)
    NULLIFY(mo_eigenvalues, deriv_set, deriv)
    NULLIFY(orb_dens)
    logger => cp_error_get_logger(error)
    ionode = (logger%para_env%mepos==logger%para_env%source)
    output_unit = cp_logger_get_default_unit_nr(logger)

    IF (ionode) WRITE (output_unit, *) " in "//routineP

    CALL get_qs_env(qs_env=qs_env,&
                    rho=rho_struct,&
                    pw_env=pw_env,&
                    cell=cell,&
                    dft_control=dft_control,&
                    mos=molecular_orbitals)
    xc_control = dft_control%xc_control
    CALL pw_env_get(pw_env,auxbas_pw_pool=auxbas_pw_pool,error=error)

    nspins = dft_control%nspins
    lsd = (nspins>1)

    ALLOCATE(rho_r(nspins), orb_dens(nspins), single_mo_coeff(nspins), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO ispin=1,nspins
       rho_r(ispin)%pw => rho_struct%rho_r(ispin)%pw
       NULLIFY(orb_dens(ispin)%matrix)
       CALL replicate_matrix(rho_struct%rho_ao(ispin)%matrix, &
                             orb_dens(ispin)%matrix, &
                             "orbital density", .TRUE.)
    END DO

    !----------------------------------------!
    ! Construct the LB94 potential in vxc_LB !
    !----------------------------------------!
    dft_control%xc_control%functionals = ""
    CALL xc_get_functional ( 1, "SLATER", dft_control%xc_control%functionals, &
                                          dft_control%xc_control%gradient_functionals, &
                                          dft_control%xc_control%crossterms )
    DO i=1, 3
       IF (dft_control%xc_control%functionals(i) == "") &
            dft_control%xc_control%functionals(i) = "NONE"
    END DO
    CALL xc_vxc_pw_create(vxc_tmp, vxc_tau, xc_energy, rho_r, rho_g, tau,&
                          dft_control, cell, auxbas_pw_pool,  &
                          error)
    dft_control%xc_control%functionals = ""
    CALL xc_get_functional ( 2, "PZ", dft_control%xc_control%functionals, &
                                      dft_control%xc_control%gradient_functionals, &
                                      dft_control%xc_control%crossterms )
    DO i=1, 3
       IF (dft_control%xc_control%functionals(i) == "") &
            dft_control%xc_control%functionals(i) = "NONE"
    END DO
    CALL xc_vxc_pw_create(vxc_LB, vxc_tau, xc_energy, rho_r, rho_g, tau,&
                          dft_control, cell, auxbas_pw_pool, &
                          error)
    DO ispin=1, nspins
       vxc_LB(ispin)%pw%cr3d = vxc_LB(ispin)%pw%cr3d + alpha*vxc_tmp(ispin)%pw%cr3d
    END DO

    !--------------------------------------!
    ! create the density needed to calc lb !
    !--------------------------------------!
    bo = rho_r(1)%pw%pw_grid%bounds_local
    CALL xc_rho_set_create(rho_set, bo, &
                           dft_control%xc_control%density_cut, &
                           dft_control%xc_control%gradient_cut, &
                           dft_control%xc_control%tau_cut, &
                           error=error)
    CALL xc_rho_cflags_setall(needs,.FALSE.,error)
    IF (lsd) THEN
       CALL xb88_lsd_info(needs=needs)
    ELSE
       CALL xb88_lda_info(needs=needs)
    END IF
    CALL xc_rho_set_update(rho_set, rho_r, rho_g, tau, needs, &
                           dft_control%xc_control, cell, auxbas_pw_pool, error)
    CALL add_lb_pot(vxc_LB, rho_set, lsd, xc_control%density_cut, xc_control%gradient_cut, error)

    !-----------------------------------------------------------------------------------!
    ! Construct 2 times PBE one particle density from the PZ correlation energy density !
    !-----------------------------------------------------------------------------------!
    CALL xc_dset_create(deriv_set, local_bounds=bo, error=error)
    DO i=1,SIZE(dft_control%xc_control%functionals)
       CALL xc_functional_eval(dft_control%xc_control%functionals(i), &
                               lsd=lsd,&
                               rho_set=rho_set, &
                               deriv_set=deriv_set,&
                               deriv_order=0, &
                               error=error)
    END DO
    deriv => xc_dset_get_derivative(deriv_set, "", error=error)
    CALL xc_derivative_get(deriv,deriv_data=e_uniform,error=error)

    ALLOCATE(vxc_GLLB(nspins))
    DO ispin=1, nspins
       CALL pw_pool_create_pw(auxbas_pw_pool, vxc_GLLB(ispin)%pw,&
                              use_data = REALDATA3D,&
                              in_space = REALSPACE, error=error)
    END DO

    DO ispin=1, nspins
       DO k=bo(1,3), bo(2,3)
          DO j=bo(1,2), bo(2,2)
             DO i=bo(1,1), bo(2,1)
                IF (nspins == 1) THEN
                      e_uniform(i,j,k) = e_uniform(i,j,k) / rho_set%rho(i,j,k)
                      vxc_GLLB(ispin)%pw%cr3d(i,j,k) = 2.0_dp * &
                           calc_ecpbe_r(rho_set%rho(i,j,k), rho_set%norm_drho(i,j,k), &
                           e_uniform(i,j,k), xc_control%density_cut, xc_control%gradient_cut) + &
                           2.0_dp * &
                           calc_expbe_r(rho_set%rho(i,j,k), rho_set%norm_drho(i,j,k), &
                                        xc_control%density_cut, xc_control%gradient_cut)
                ELSE                
                      e_uniform(i,j,k) = e_uniform(i,j,k) / rho_set%rho(i,j,k)
                      vxc_GLLB(ispin)%pw%cr3d(i,j,k) = 2.0_dp * &
                           calc_ecpbe_u(rho_set%rhoa(i,j,k), rho_set%rhob(i,j,k), rho_set%norm_drho(i,j,k), &
                           e_uniform(i,j,k), xc_control%density_cut, xc_control%gradient_cut) + &
                           2.0_dp * &
                           calc_expbe_r(rho_set%rho(i,j,k), rho_set%norm_drho(i,j,k), &
                                        xc_control%density_cut, xc_control%gradient_cut)
                END IF
             END DO
          END DO
       END DO
    END DO

    CALL xc_dset_release(deriv_set, error=error)

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
                       homo=homo)
       CALL cp_fm_create(single_mo_coeff(ispin)%matrix, &
                         mo_coeff%matrix_struct, &
                         "orbital density matrix", error)

       CALL cp_fm_get_info(single_mo_coeff(ispin)%matrix, &
            nrow_global=nrow(ispin), ncol_global=ncol(ispin))
       ALLOCATE(coeff_col(nrow(ispin),1))

       CALL pw_zero(vxc_tmp(ispin)%pw)

       DO orb=1, homo-1

          efac = K_rho*SQRT(mo_eigenvalues(homo)-mo_eigenvalues(orb))
          IF (.NOT.lsd) efac = 2.0_dp * efac

          CALL cp_fm_set_all(single_mo_coeff(ispin)%matrix, 0.0_dp)
          CALL cp_fm_get_submatrix(mo_coeff, coeff_col, &
                                   1, orb, nrow(ispin), 1, error=error)
          CALL cp_fm_set_submatrix(single_mo_coeff(ispin)%matrix, coeff_col, &
                                   1, orb, error=error)
          CALL set_matrix(orb_dens(ispin)%matrix, value=0.0_dp)
          CALL cp_sm_plus_fm_fm_t(orb_dens(ispin)%matrix,&
                                  matrix_v=single_mo_coeff(ispin)%matrix,&
                                  matrix_g=single_mo_coeff(ispin)%matrix,&
                                  ncol=ncol(ispin), &
                                  alpha=1.0_dp)
          CALL pw_zero(orbital%pw)
          CALL pw_zero(orbital_g%pw)
          CALL calculate_rho_elec(orb_dens(ispin)%matrix, &
                                  orbital, orbital_g, &
                                  tot_rho_psi, qs_env, error=error)

          vxc_tmp(ispin)%pw%cr3d = vxc_tmp(ispin)%pw%cr3d + &
               efac * orbital%pw%cr3d

       END DO
       DEALLOCATE(coeff_col)

       DO k=bo(1,3), bo(2,3)
          DO j=bo(1,2), bo(2,2)
             DO i=bo(1,1), bo(2,1)
                IF (rho_r(ispin)%pw%cr3d(i,j,k) > xc_control%density_cut) THEN
                   vxc_tmp(ispin)%pw%cr3d(i,j,k) = vxc_tmp(ispin)%pw%cr3d(i,j,k) / &
                                                   rho_r(ispin)%pw%cr3d(i,j,k)
                ELSE
                   vxc_tmp(ispin)%pw%cr3d(i,j,k) = 0.0_dp
                END IF
             END DO
          END DO
       END DO

       vxc_GLLB(ispin)%pw%cr3d = vxc_GLLB(ispin)%pw%cr3d + vxc_tmp(ispin)%pw%cr3d

    END DO
    
    !---------------!
    ! Assemble SAOP !
    !---------------!
    ALLOCATE(vxc_SAOP(nspins), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO ispin=1, nspins
       
       CALL get_mo_set(molecular_orbitals(ispin)%mo_set,&
                       mo_coeff=mo_coeff, &
                       eigenvalues=mo_eigenvalues,&
                       homo=homo)
       CALL pw_pool_init_coeff(auxbas_pw_pool, vxc_SAOP(ispin), &
                               use_data=REALDATA3D, in_space=REALSPACE, &
                               error=error)
       CALL pw_zero(vxc_SAOP(ispin)%pw)

       ALLOCATE(coeff_col(nrow(ispin),1))

       DO orb=1, homo
          
          we_LB   = EXP(-2.0_dp*(mo_eigenvalues(homo)-mo_eigenvalues(orb))**2)
          we_GLLB = 1.0_dp - we_LB
          IF (.NOT.lsd) THEN
             we_LB   = 2.0_dp * we_LB
             we_GLLB = 2.0_dp * we_GLLB
          END IF
          
          vxc_tmp(ispin)%pw%cr3d = we_LB*vxc_LB(ispin)%pw%cr3d + &
                                   we_GLLB*vxc_GLLB(ispin)%pw%cr3d

          CALL cp_fm_set_all(single_mo_coeff(ispin)%matrix, 0.0_dp)
          CALL cp_fm_get_submatrix(mo_coeff, coeff_col, &
                                   1, orb, nrow(ispin), 1, error=error)
          CALL cp_fm_set_submatrix(single_mo_coeff(ispin)%matrix, coeff_col, &
                                   1, orb, error=error)
          CALL set_matrix(orb_dens(ispin)%matrix, value=0.0_dp)
          CALL cp_sm_plus_fm_fm_t(orb_dens(ispin)%matrix,&
                                  matrix_v=single_mo_coeff(ispin)%matrix,&
                                  matrix_g=single_mo_coeff(ispin)%matrix,&
                                  ncol=ncol(ispin), &
                                  alpha=1.0_dp)
          CALL pw_zero(orbital%pw)
          CALL pw_zero(orbital_g%pw)
          CALL calculate_rho_elec(orb_dens(ispin)%matrix, &
                                  orbital, orbital_g, &
                                  tot_rho_psi, qs_env, error=error)

!TC          CALL calculate_wavefunction(mo_coeff,orb,orbital2,orbital2_g,qs_env,error)
!TC          write (*,*) orb, maxval(abs(orbital2%pw%cr3d-orbital%pw%cr3d))
          
          vxc_SAOP(ispin)%pw%cr3d = vxc_SAOP(ispin)%pw%cr3d + &
               orbital%pw%cr3d * vxc_tmp(ispin)%pw%cr3d

       END DO

       CALL cp_fm_release(single_mo_coeff(ispin)%matrix, error=error)
       CALL deallocate_matrix(orb_dens(ispin)%matrix)

       DEALLOCATE(coeff_col)

       DO k=bo(1,3), bo(2,3)
          DO j=bo(1,2), bo(2,2)
             DO i=bo(1,1), bo(2,1)
                IF (rho_r(ispin)%pw%cr3d(i,j,k) > xc_control%density_cut) THEN
                   vxc_SAOP(ispin)%pw%cr3d(i,j,k) = vxc_SAOP(ispin)%pw%cr3d(i,j,k) / &
                                                    rho_r(ispin)%pw%cr3d(i,j,k)
                ELSE
                   vxc_SAOP(ispin)%pw%cr3d(i,j,k) = 0.0_dp
                END IF
             END DO
          END DO
       END DO

    END DO
       
    CALL xc_rho_set_release(rho_set, auxbas_pw_pool, error)
    CALL pw_pool_give_back_coeff(auxbas_pw_pool, orbital, error=error)
    CALL pw_pool_give_back_coeff(auxbas_pw_pool, orbital_g, error=error)
    DO ispin=1, nspins
       CALL pw_pool_give_back_pw(auxbas_pw_pool, vxc_GLLB(ispin)%pw, error=error)
       CALL pw_release(vxc_LB(ispin)%pw)
       CALL pw_release(vxc_tmp(ispin)%pw)
    END DO
    DEALLOCATE(vxc_GLLB, vxc_LB, vxc_tmp, orb_dens, single_mo_coeff, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    !--------------------!
    ! Do the integration !
    !--------------------!
    DO ispin=1, nspins

       vxc_SAOP(ispin)%pw%cr3d = vxc_SAOP(ispin)%pw%cr3d * vxc_SAOP(ispin)%pw%pw_grid%dvol

       CALL integrate_v_rspace(vxc_SAOP(ispin), rho_struct%rho_ao(ispin), &
                               ks_matrix(ispin), qs_env, &
                               calculate_forces=.FALSE.,&
                               gapw=.FALSE.,error=error)

       CALL pw_pool_give_back_coeff(auxbas_pw_pool, vxc_SAOP(ispin), error=error)
    END DO
    
    DEALLOCATE(rho_r, vxc_SAOP, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    !----------------------------------!
    ! Restore the values of xc_control !
    !----------------------------------!
    dft_control%xc_control = xc_control

  END SUBROUTINE add_saop_pot

  !*****************************************************************

  SUBROUTINE add_lb_pot(pot, rho_set, lsd, rho_cut, grad_cut, error)

    TYPE(pw_p_type), DIMENSION(:), POINTER   :: pot
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    LOGICAL, INTENT(IN)                      :: lsd
    REAL(KIND=dp), INTENT(IN)                :: rho_cut, grad_cut
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    REAL(KIND=dp), PARAMETER                 :: ob3 = 1.0_dp/3.0_dp

    INTEGER                                  :: i, j, k
    INTEGER, DIMENSION(2, 3)                 :: bo
    REAL(KIND=dp)                            :: n, n_13, x, x2

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
                IF (rho_set%rho(i,j,k) > rho_cut) THEN
                   n    = rho_set%rho(i,j,k)/2.0_dp
                   n_13 = n**ob3
                   x    = (rho_set%norm_drho(i,j,k)/2.0_dp)/(n*n_13)
                   x2   = x*x
                   pot(1)%pw%cr3d(i,j,k) = pot(1)%pw%cr3d(i,j,k) - &
                        beta*x2*n_13/(1.0_dp+3.0_dp*beta*x*LOG(x + SQRT(x2+1.0_dp)))
                END IF
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

  FUNCTION calc_ecpbe_u(ra, rb, ngr, ec_unif, rc, ngrc) RESULT(res)

    REAL(kind=dp), INTENT(in)                :: ra, rb, ngr, ec_unif, rc, ngrc
    REAL(kind=dp)                            :: res

    REAL(kind=dp), parameter                 :: ob3 = 1.0_dp/3.0_dp, &
                                                tb3 = 2.0_dp/3.0_dp

    REAL(kind=dp)                            :: A, At2, H, kf, kl, ks, phi, &
                                                phi3, r, t2, zeta

    r = ra + rb
    H = 0.0_dp
    IF (r > rc .AND. ngr > ngrc) THEN
       zeta = (ra-rb)/r
       phi  = ((1.0_dp+zeta)**tb3 + (1.0_dp-zeta)**tb3)/2.0_dp
       phi3 = phi*phi*phi
       kf   = (3.0_dp*r*pi*pi)**ob3
       ks   = SQRT(4.0_dp*kf/pi)
       t2   = (ngr/(2.0_dp*phi*ks*r))**2
       A    = beta_ec/gamma_saop/(EXP(-ec_unif/(gamma_saop*phi3))-1.0_dp)
       At2  = A*t2
       kl   = (1.0_dp + At2)/(1.0_dp + At2 + At2*At2)
       H    = gamma_saop * LOG(1.0_dp + beta_ec/gamma_saop*t2*kl)
    END IF
    res = ec_unif + H

  END FUNCTION calc_ecpbe_u

  !*****************************************************************

  FUNCTION calc_ecpbe_r(r, ngr, ec_unif, rc, ngrc) RESULT(res)

    REAL(kind=dp), INTENT(in)                :: r, ngr, ec_unif, rc, ngrc
    REAL(kind=dp)                            :: res

    REAL(kind=dp), parameter                 :: ob3 = 1.0_dp/3.0_dp

    REAL(kind=dp)                            :: A, At2, H, kf, kl, ks, t2

    H = 0.0_dp
    IF (r > rc .AND. ngr > ngrc) THEN
       kf   = (3.0_dp*r*pi*pi)**(1.0_dp/3.0_dp)
       ks   = SQRT(4.0_dp*kf/pi)
       t2   = (ngr/(2.0_dp*ks*r))**2
       A    = beta_ec/gamma_saop/(EXP(-ec_unif/gamma_saop)-1.0_dp)
       At2  = A*t2
       kl   = (1.0_dp + At2)/(1.0_dp + At2 + At2*At2)
       H    = gamma_saop * LOG(1.0_dp + beta_ec/gamma_saop*t2*kl)
    END IF
    res  = ec_unif + H
    
  END FUNCTION calc_ecpbe_r

  !*****************************************************************

  FUNCTION calc_expbe_u(ra, rb, ngr, rc, ngrc) RESULT(res)

    REAL(kind=dp), INTENT(in)                :: ra, rb, ngr, rc, ngrc
    REAL(kind=dp)                            :: res

    REAL(kind=dp)                            :: r

    r = ra + rb
    res = calc_expbe_r(r, ngr, rc, ngrc)
    
  END FUNCTION calc_expbe_u

  !*****************************************************************

  FUNCTION calc_expbe_r(r, ngr, rc, ngrc) RESULT(res)

    REAL(kind=dp), INTENT(in)                :: r, ngr, rc, ngrc
    REAL(kind=dp)                            :: res

    REAL(kind=dp)                            :: ex_unif, fx, kf, s

    IF (r > rc) THEN
       kf      = (3.0_dp*r*pi*pi)**(1.0_dp/3.0_dp)
       ex_unif = -3.0_dp*kf/(4.0_dp*pi)
       fx = 1.0_dp
       IF (ngr > ngrc) THEN
          s     = ngr/(2.0_dp*kf*r)
          fx    = fx + kappa - kappa/(1.0_dp + mu*s*s/kappa)
       END IF
       res = ex_unif*fx
    ELSE
       res = 0.0_dp
    END IF
    
  END FUNCTION calc_expbe_r

  !*****************************************************************

END MODULE xc_pot_saop
