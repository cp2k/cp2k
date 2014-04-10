!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculate the saop potential
! *****************************************************************************
MODULE xc_pot_saop

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cp_array_r_utils,                ONLY: cp_1d_r_p_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_copy,&
                                             cp_dbcsr_init,&
                                             cp_dbcsr_set
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_allocate_matrix_set,&
                                             cp_dbcsr_deallocate_matrix,&
                                             cp_dbcsr_deallocate_matrix_set,&
                                             cp_dbcsr_plus_fm_fm_t
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_p_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_get_submatrix,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_set_submatrix,&
                                             cp_fm_type
  USE f77_blas
  USE input_constants,                 ONLY: do_method_gapw,&
                                             oe_gllb,&
                                             oe_lb,&
                                             oe_saop,&
                                             xc_funct_no_shortcut
  USE input_section_types,             ONLY: &
       section_vals_create, section_vals_duplicate, &
       section_vals_get_subs_vals, section_vals_release, section_vals_retain, &
       section_vals_set_subs_vals, section_vals_type, section_vals_val_get, &
       section_vals_val_set
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_methods,                      ONLY: pw_zero
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type,&
                                             pw_release
  USE qs_collocate_density,            ONLY: calculate_rho_elec
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_gapw_densities,               ONLY: prepare_gapw_den
  USE qs_grid_atom,                    ONLY: grid_atom_type
  USE qs_harmonics_atom,               ONLY: harmonics_atom_type
  USE qs_integrate_potential,          ONLY: integrate_v_rspace
  USE qs_ks_atom,                      ONLY: update_ks_atom
  USE qs_local_rho_types,              ONLY: local_rho_set_create,&
                                             local_rho_set_release,&
                                             local_rho_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_rho_atom_methods,             ONLY: allocate_rho_atom_internals,&
                                             calculate_rho_atom_coeff
  USE qs_rho_atom_types,               ONLY: get_rho_atom,&
                                             rho_atom_coeff,&
                                             rho_atom_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE qs_vxc_atom,                     ONLY: calc_rho_angular,&
                                             gaVxcgb_noGC
  USE util,                            ONLY: get_limit
  USE xc,                              ONLY: xc_vxc_pw_create
  USE xc_atom,                         ONLY: fill_rho_set,&
                                             vxc_of_r_new,&
                                             xc_rho_set_atom_update
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                             xc_dset_create,&
                                             xc_dset_get_derivative,&
                                             xc_dset_release,&
                                             xc_dset_zero_all
  USE xc_derivative_types,             ONLY: xc_derivative_get,&
                                             xc_derivative_type
  USE xc_derivatives,                  ONLY: xc_functionals_eval
  USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_setall,&
                                             xc_rho_cflags_type
  USE xc_rho_set_types,                ONLY: xc_rho_set_create,&
                                             xc_rho_set_release,&
                                             xc_rho_set_type,&
                                             xc_rho_set_update
  USE xc_xbecke88,                     ONLY: xb88_lda_info,&
                                             xb88_lsd_info
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: add_saop_pot

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_pot_saop'

  ! should be eliminated
  REAL(KIND=dp), PARAMETER :: alpha=1.19_dp, beta=0.01_dp, K_rho=0.42_dp
  REAL(KIND=dp), PARAMETER :: kappa = 0.804_dp, mu = 0.21951_dp, &
                              beta_ec  = 0.066725_dp, gamma_saop = 0.031091_dp

CONTAINS

! *****************************************************************************
  SUBROUTINE add_saop_pot(ks_matrix, qs_env, oe_corr, error)

    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: oe_corr
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'add_saop_pot', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: dft_method_id, homo, i, &
                                                ispin, j, k, nspins, orb, &
                                                stat, xc_deriv_method_id, &
                                                xc_rho_smooth_id
    INTEGER, DIMENSION(2)                    :: ncol, nrow
    INTEGER, DIMENSION(2, 3)                 :: bo
    LOGICAL                                  :: failure, gapw, ionode, lsd
    REAL(KIND=dp)                            :: density_cut, efac, &
                                                gradient_cut, tau_cut, &
                                                tot_rho_psi, we_GLLB, we_LB, &
                                                xc_energy
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: coeff_col
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mo_eigenvalues
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: dummy, e_uniform
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: orbital_density_matrix
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: single_mo_coeff
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: molecular_orbitals
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: orbital, orbital_g
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_g, rho_r, tau, vxc_GLLB, &
                                                vxc_LB, vxc_SAOP, vxc_tau, &
                                                vxc_tmp
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_rho_type), POINTER               :: rho_struct
    TYPE(section_vals_type), POINTER         :: input, xc_fun_section_orig, &
                                                xc_fun_section_tmp, &
                                                xc_section_orig, &
                                                xc_section_tmp
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    TYPE(xc_derivative_type), POINTER        :: deriv
    TYPE(xc_rho_cflags_type)                 :: needs
    TYPE(xc_rho_set_type), POINTER           :: rho_set

    failure = .FALSE.

    NULLIFY(pw_env, auxbas_pw_pool, input)
    NULLIFY(rho_g, rho_r, tau, rho_struct, rho_set, e_uniform)
    NULLIFY(vxc_GLLB, vxc_LB, vxc_tmp, vxc_SAOP, vxc_tau)
    NULLIFY(mo_eigenvalues, deriv_set, deriv)
    NULLIFY(orbital_density_matrix, dummy,xc_section_tmp,xc_fun_section_tmp)
    logger => cp_error_get_logger(error)
    ionode = (logger%para_env%mepos==logger%para_env%source)

    IF (ionode) WRITE (*, *) " in "//routineP

    CALL get_qs_env(qs_env=qs_env,&
                    rho=rho_struct,&
                    pw_env=pw_env,&
                    input=input,&
                    mos=molecular_orbitals,error=error)
    CALL section_vals_val_get(input,"DFT%QS%METHOD",i_val=dft_method_id,error=error)
    gapw = (dft_method_id==do_method_gapw)

    xc_section_orig => section_vals_get_subs_vals(input,"DFT%XC",error=error)
    CALL section_vals_retain(xc_section_orig,error=error)
    CALL section_vals_duplicate(xc_section_orig,xc_section_tmp,error=error)

    CALL section_vals_val_get(xc_section_orig,"DENSITY_CUTOFF",&
         r_val=density_cut,error=error)
    CALL section_vals_val_get(xc_section_orig,"GRADIENT_CUTOFF",&
         r_val=gradient_cut,error=error)
    CALL section_vals_val_get(xc_section_orig,"TAU_CUTOFF",&
         r_val=tau_cut,error=error)

    CALL pw_env_get(pw_env,auxbas_pw_pool=auxbas_pw_pool,error=error)

    CALL section_vals_val_get(input,"DFT%LSD",l_val=lsd,error=error)
    IF (lsd) THEN
       nspins=2
    ELSE
       nspins=1
    END IF

    ALLOCATE(rho_r(nspins), single_mo_coeff(nspins), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_allocate_matrix_set(orbital_density_matrix,nspins,error=error)
    DO ispin=1,nspins
       rho_r(ispin)%pw => rho_struct%rho_r(ispin)%pw
       ALLOCATE(orbital_density_matrix(ispin)%matrix)
       CALL cp_dbcsr_init(orbital_density_matrix(ispin)%matrix, error=error)
       CALL cp_dbcsr_copy(orbital_density_matrix(ispin)%matrix,&
            rho_struct%rho_ao(ispin)%matrix,"orbital density",error=error)
    END DO
    bo = rho_r(1)%pw%pw_grid%bounds_local

    !---------------------------!
    ! create the density needed !
    !---------------------------!
    CALL xc_rho_set_create(rho_set, bo, &
                           density_cut, &
                           gradient_cut, &
                           tau_cut, &
                           error=error)
    CALL xc_rho_cflags_setall(needs,.FALSE.,error)
    IF (lsd) THEN
       CALL xb88_lsd_info(needs=needs,error=error)
       needs%norm_drho = .TRUE.
    ELSE
       CALL xb88_lda_info(needs=needs,error=error)
    END IF
    CALL section_vals_val_get(xc_section_orig,"XC_GRID%XC_DERIV",&
         i_val=xc_deriv_method_id,error=error)
    CALL section_vals_val_get(xc_section_orig,"XC_GRID%XC_SMOOTH_RHO",&
         i_val=xc_rho_smooth_id,error=error)
    CALL xc_rho_set_update(rho_set, rho_r, rho_g, tau, needs, &
                           xc_deriv_method_id,&
                           xc_rho_smooth_id,&
                           auxbas_pw_pool, error)

    !----------------------------------------!
    ! Construct the LB94 potential in vxc_LB !
    !----------------------------------------!
    xc_fun_section_orig => section_vals_get_subs_vals(xc_section_orig,&
         "XC_FUNCTIONAL",error=error)
    CALL section_vals_create(xc_fun_section_tmp,xc_fun_section_orig%section,&
         error=error)
    CALL section_vals_val_set(xc_fun_section_tmp,"_SECTION_PARAMETERS_",&
         i_val=xc_funct_no_shortcut,error=error)
    CALL section_vals_val_set(xc_fun_section_tmp,"XALPHA%_SECTION_PARAMETERS_",&
         l_val=.TRUE.,error=error)
    CALL section_vals_set_subs_vals(xc_section_tmp,"XC_FUNCTIONAL",&
         xc_fun_section_tmp,error=error)

    CALL xc_vxc_pw_create(vxc_tmp, vxc_tau, xc_energy, rho_r, rho_g, tau,&
                          xc_section_tmp, auxbas_pw_pool,  &
                          error)

    CALL section_vals_val_set(xc_fun_section_tmp,"XALPHA%_SECTION_PARAMETERS_",&
         l_val=.FALSE.,error=error)
    CALL section_vals_val_set(xc_fun_section_tmp,"PZ81%_SECTION_PARAMETERS_",&
         l_val=.TRUE.,error=error)

    CALL xc_vxc_pw_create(vxc_LB, vxc_tau, xc_energy, rho_r, rho_g, tau,&
                          xc_section_tmp, auxbas_pw_pool, &
                          error)

    DO ispin=1, nspins
       vxc_LB(ispin)%pw%cr3d = vxc_LB(ispin)%pw%cr3d + alpha*vxc_tmp(ispin)%pw%cr3d
    END DO

    DO ispin=1, nspins
       dummy => vxc_tmp(ispin)%pw%cr3d
       CALL add_lb_pot(dummy, rho_set, lsd, ispin, error)
       vxc_LB(ispin)%pw%cr3d = vxc_LB(ispin)%pw%cr3d - vxc_tmp(ispin)%pw%cr3d
    END DO
    NULLIFY(dummy)

    !-----------------------------------------------------------------------------------!
    ! Construct 2 times PBE one particle density from the PZ correlation energy density !
    !-----------------------------------------------------------------------------------!
    CALL xc_dset_create(deriv_set, local_bounds=bo, error=error)
    CALL xc_functionals_eval(xc_fun_section_tmp, &
         lsd=lsd,&
         rho_set=rho_set, &
         deriv_set=deriv_set,&
         deriv_order=0, &
         error=error)

    deriv => xc_dset_get_derivative(deriv_set, "", error=error)
    CALL xc_derivative_get(deriv,deriv_data=e_uniform,error=error)

    ALLOCATE(vxc_GLLB(nspins))
    DO ispin=1, nspins
       CALL pw_pool_create_pw(auxbas_pw_pool, vxc_GLLB(ispin)%pw,&
                              use_data = REALDATA3D,&
                              in_space = REALSPACE, error=error)
    END DO

    DO ispin=1, nspins
       dummy => vxc_GLLB(ispin)%pw%cr3d
       CALL calc_2excpbe(dummy, rho_set, e_uniform, lsd)
    END DO
    NULLIFY(dummy)

    CALL xc_dset_release(deriv_set, error=error)

    CALL pw_pool_create_pw(auxbas_pw_pool,orbital%pw,&
                            use_data = REALDATA3D,&
                            in_space = REALSPACE, error=error)
    CALL pw_pool_create_pw(auxbas_pw_pool,orbital_g%pw,&
                            use_data = COMPLEXDATA1D,&
                            in_space = RECIPROCALSPACE, error=error)

    DO ispin=1, nspins

       CALL get_mo_set(molecular_orbitals(ispin)%mo_set,&
                       mo_coeff=mo_coeff, &
                       eigenvalues=mo_eigenvalues,&
                       homo=homo)
       CALL cp_fm_create(single_mo_coeff(ispin)%matrix, &
                         mo_coeff%matrix_struct, &
                         "orbital density matrix", error=error)

       CALL cp_fm_get_info(single_mo_coeff(ispin)%matrix, &
            nrow_global=nrow(ispin), ncol_global=ncol(ispin),error=error)
       ALLOCATE(coeff_col(nrow(ispin),1))

       CALL pw_zero(vxc_tmp(ispin)%pw,error=error)

       DO orb=1, homo-1

          efac = K_rho*SQRT(mo_eigenvalues(homo)-mo_eigenvalues(orb))
          IF (.NOT.lsd) efac = 2.0_dp * efac

          CALL cp_fm_set_all(single_mo_coeff(ispin)%matrix, 0.0_dp,error=error)
          CALL cp_fm_get_submatrix(mo_coeff, coeff_col, &
                                   1, orb, nrow(ispin), 1, error=error)
          CALL cp_fm_set_submatrix(single_mo_coeff(ispin)%matrix, coeff_col, &
                                   1, orb, error=error)
          CALL cp_dbcsr_set(orbital_density_matrix(ispin)%matrix,0.0_dp,error=error)
          CALL cp_dbcsr_plus_fm_fm_t(orbital_density_matrix(ispin)%matrix,&
                                  matrix_v=single_mo_coeff(ispin)%matrix,&
                                  matrix_g=single_mo_coeff(ispin)%matrix,&
                                  ncol=ncol(ispin), &
                                  alpha=1.0_dp,error=error)
          CALL pw_zero(orbital%pw,error=error)
          CALL pw_zero(orbital_g%pw,error=error)
          CALL calculate_rho_elec(orbital_density_matrix(ispin)%matrix, &
                                  orbital, orbital_g, &
                                  tot_rho_psi, qs_env, error=error)

          vxc_tmp(ispin)%pw%cr3d = vxc_tmp(ispin)%pw%cr3d + &
               efac * orbital%pw%cr3d

       END DO
       DEALLOCATE(coeff_col)

       DO k=bo(1,3), bo(2,3)
          DO j=bo(1,2), bo(2,2)
             DO i=bo(1,1), bo(2,1)
                IF (rho_r(ispin)%pw%cr3d(i,j,k) > density_cut) THEN
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
       CALL pw_pool_create_pw(auxbas_pw_pool, vxc_SAOP(ispin)%pw, &
                               use_data=REALDATA3D, in_space=REALSPACE, &
                               error=error)
       CALL pw_zero(vxc_SAOP(ispin)%pw,error=error)

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

          CALL cp_fm_set_all(single_mo_coeff(ispin)%matrix, 0.0_dp,error=error)
          CALL cp_fm_get_submatrix(mo_coeff, coeff_col, &
                                   1, orb, nrow(ispin), 1, error=error)
          CALL cp_fm_set_submatrix(single_mo_coeff(ispin)%matrix, coeff_col, &
                                   1, orb, error=error)
          CALL cp_dbcsr_set(orbital_density_matrix(ispin)%matrix,0.0_dp,error=error)
          CALL cp_dbcsr_plus_fm_fm_t(orbital_density_matrix(ispin)%matrix,&
                                  matrix_v=single_mo_coeff(ispin)%matrix,&
                                  matrix_g=single_mo_coeff(ispin)%matrix,&
                                  ncol=ncol(ispin), &
                                  alpha=1.0_dp,error=error)
          CALL pw_zero(orbital%pw,error=error)
          CALL pw_zero(orbital_g%pw,error=error)
          CALL calculate_rho_elec(orbital_density_matrix(ispin)%matrix, &
                                  orbital, orbital_g, &
                                  tot_rho_psi, qs_env, error=error)

!TC          CALL calculate_wavefunction(mo_coeff,orb,orbital2,orbital2_g,qs_env,error)
!TC          write (*,*) orb, maxval(abs(orbital2%pw%cr3d-orbital%pw%cr3d))

          vxc_SAOP(ispin)%pw%cr3d = vxc_SAOP(ispin)%pw%cr3d + &
               orbital%pw%cr3d * vxc_tmp(ispin)%pw%cr3d

       END DO

       CALL cp_fm_release(single_mo_coeff(ispin)%matrix, error=error)
       CALL cp_dbcsr_deallocate_matrix(orbital_density_matrix(ispin)%matrix,error=error)

       DEALLOCATE(coeff_col)

       DO k=bo(1,3), bo(2,3)
          DO j=bo(1,2), bo(2,2)
             DO i=bo(1,1), bo(2,1)
                IF (rho_r(ispin)%pw%cr3d(i,j,k) > density_cut) THEN
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
    CALL pw_pool_give_back_pw(auxbas_pw_pool, orbital%pw, error=error)
    CALL pw_pool_give_back_pw(auxbas_pw_pool, orbital_g%pw, error=error)

    !--------------------!
    ! Do the integration !
    !--------------------!
    DO ispin=1, nspins

       IF (oe_corr == oe_lb) THEN
          vxc_SAOP(ispin)%pw%cr3d = vxc_LB(ispin)%pw%cr3d
       ELSE IF (oe_corr == oe_gllb) THEN
          vxc_SAOP(ispin)%pw%cr3d = vxc_GLLB(ispin)%pw%cr3d
       END IF
       vxc_SAOP(ispin)%pw%cr3d = vxc_SAOP(ispin)%pw%cr3d * vxc_SAOP(ispin)%pw%pw_grid%dvol

       CALL integrate_v_rspace(vxc_SAOP(ispin), rho_struct%rho_ao(ispin), &
                               ks_matrix(ispin), qs_env, &
                               calculate_forces=.FALSE.,&
                               gapw=gapw, error=error)

    END DO

    DO ispin=1, nspins
       CALL pw_pool_give_back_pw(auxbas_pw_pool, vxc_SAOP(ispin)%pw, error=error)
       CALL pw_pool_give_back_pw(auxbas_pw_pool, vxc_GLLB(ispin)%pw, error=error)
       CALL pw_release(vxc_LB(ispin)%pw,error=error)
       CALL pw_release(vxc_tmp(ispin)%pw,error=error)
    END DO
    DEALLOCATE(vxc_GLLB, vxc_LB, vxc_tmp, orbital_density_matrix, single_mo_coeff, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE(rho_r, vxc_SAOP, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL section_vals_release(xc_fun_section_tmp,error=error)
    CALL section_vals_release(xc_section_tmp,error=error)
    CALL section_vals_release(xc_section_orig,error=error)

    !-----------------------!
    ! Call the GAPW routine !
    !-----------------------!
    IF (gapw) THEN
       CALL gapw_add_atomic_saop_pot(qs_env, oe_corr, error)
    END IF

  END SUBROUTINE add_saop_pot

! *****************************************************************************
  SUBROUTINE gapw_add_atomic_saop_pot(qs_env, oe_corr, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: oe_corr
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'gapw_add_atomic_saop_pot', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ia, iat, iatom, ikind, ir, &
                                                ispin, mepos, na, natom, nr, &
                                                nspins, num_pe, orb, stat
    INTEGER, DIMENSION(2)                    :: bo, homo, ncol, nrow
    INTEGER, DIMENSION(2, 3)                 :: bounds
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: failure, lsd, paw_atom
    REAL(dp), DIMENSION(:, :), POINTER       :: tau
    REAL(KIND=dp)                            :: density_cut, efac, exc, &
                                                gradient_cut, tau_cut, &
                                                we_GLLB, we_LB
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: coeff_col
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rho_h, rho_s, weight
    REAL(KIND=dp), DIMENSION(:, :, :), POINTER :: dummy, e_uniform, vtau, &
      vxc_GLLB_h, vxc_GLLB_s, vxc_LB_h, vxc_LB_s, vxc_SAOP_h, vxc_SAOP_s, &
      vxc_tmp_h, vxc_tmp_s
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: drho_h, drho_s, vxg
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(cp_1d_r_p_type), DIMENSION(:), &
      POINTER                                :: mo_eigenvalues
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: orbital_density_matrix
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mo_coeff, single_mo_coeff
    TYPE(grid_atom_type), POINTER            :: atomic_grid
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    TYPE(local_rho_type), POINTER            :: local_rho_set
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: molecular_orbitals
    TYPE(qs_rho_type), POINTER               :: rho_structure
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                                :: dr_h, dr_s, r_h, r_s
    TYPE(rho_atom_coeff), DIMENSION(:, :), &
      POINTER                                :: r_h_d, r_s_d
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom_set
    TYPE(rho_atom_type), POINTER             :: rho_atom
    TYPE(section_vals_type), POINTER         :: input, xc_fun_section_orig, &
                                                xc_fun_section_tmp, &
                                                xc_section_orig, &
                                                xc_section_tmp
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    TYPE(xc_derivative_type), POINTER        :: deriv
    TYPE(xc_rho_cflags_type)                 :: needs, needs_orbs
    TYPE(xc_rho_set_type), POINTER           :: orb_rho_set_h, orb_rho_set_s, &
                                                rho_set_h, rho_set_s

    failure = .FALSE.

    NULLIFY(weight, rho_h, rho_s, vxc_LB_h, vxc_LB_s, vxc_GLLB_h, vxc_GLLB_s, &
            vxc_tmp_h, vxc_tmp_s, vtau, dummy, e_uniform, drho_h, drho_s, vxg, atom_list, &
            atom_kind, atomic_kind_set,  deriv_set, deriv, atomic_grid, &
            harmonics, molecular_orbitals, rho_structure, r_h, r_s, dr_h, dr_s, &
            r_h_d, r_s_d, rho_atom_set, rho_atom, rho_set_h, rho_set_s, &
            mo_coeff, single_mo_coeff, mo_eigenvalues, local_rho_set, &
            orbital_density_matrix, orb_rho_set_h, orb_rho_set_s, vxc_SAOP_h, vxc_SAOP_s)

    ! tau is needed for fill_rho_set, but should never be used
    NULLIFY(tau)

    CALL get_qs_env(qs_env, input=input, &
                            rho=rho_structure, &
                            mos=molecular_orbitals, &
                            atomic_kind_set=atomic_kind_set, &
                            rho_atom_set=rho_atom_set,error=error)

    xc_section_orig => section_vals_get_subs_vals(input,"DFT%XC",error=error)
    CALL section_vals_retain(xc_section_orig,error=error)
    CALL section_vals_duplicate(xc_section_orig,xc_section_tmp,error=error)

    CALL section_vals_val_get(input,"DFT%LSD",l_val=lsd,error=error)
    IF (lsd) THEN
       nspins=2
    ELSE
       nspins=1
    END IF

    CALL section_vals_val_get(xc_section_orig,"DENSITY_CUTOFF",&
         r_val=density_cut,error=error)
    CALL section_vals_val_get(xc_section_orig,"GRADIENT_CUTOFF",&
         r_val=gradient_cut,error=error)
    CALL section_vals_val_get(xc_section_orig,"TAU_CUTOFF",&
         r_val=tau_cut,error=error)

    CALL calculate_rho_atom_coeff(qs_env, rho_structure%rho_ao,error=error)
    CALL prepare_gapw_den(qs_env, error=error)

    ALLOCATE(mo_coeff(nspins), single_mo_coeff(nspins), &
             mo_eigenvalues(nspins),&
             stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL cp_dbcsr_allocate_matrix_set(orbital_density_matrix,nspins,error=error)

    DO ispin=1, nspins
       CALL get_mo_set(molecular_orbitals(ispin)%mo_set, &
                       mo_coeff=mo_coeff(ispin)%matrix, &
                       eigenvalues=mo_eigenvalues(ispin)%array, &
                       homo=homo(ispin))
       CALL cp_fm_create(single_mo_coeff(ispin)%matrix, &
                         mo_coeff(ispin)%matrix%matrix_struct, &
                         "orbital density matrix", error=error)
       CALL cp_fm_get_info(single_mo_coeff(ispin)%matrix, &
                           nrow_global=nrow(ispin), ncol_global=ncol(ispin),error=error)
       ALLOCATE(orbital_density_matrix(ispin)%matrix)
       CALL cp_dbcsr_init(orbital_density_matrix(ispin)%matrix, error=error)
       CALL cp_dbcsr_copy(orbital_density_matrix(ispin)%matrix,&
            rho_structure%rho_ao(ispin)%matrix, &
            "orbital density", error=error)
    END DO
    CALL local_rho_set_create(local_rho_set, error=error)
    CALL allocate_rho_atom_internals(qs_env, local_rho_set%rho_atom_set,error=error)

    DO ikind=1, SIZE(atomic_kind_set)

       atom_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind=atom_kind, atom_list=atom_list,&
                            natom=natom, paw_atom=paw_atom,&
                            harmonics=harmonics, grid_atom=atomic_grid)

       IF(.NOT. paw_atom) CYCLE

       nr = atomic_grid%nr
       na = atomic_grid%ng_sphere
       bounds(1:2,1:3) = 1
       bounds(2,1)     = na
       bounds(2,2)     = nr

       CALL xc_dset_create(deriv_set, local_bounds=bounds, error=error)

       CALL xc_rho_set_create(rho_set_h, bounds, density_cut,&
                              gradient_cut, tau_cut,&
                              error=error)
       CALL xc_rho_set_create(rho_set_s, bounds, density_cut,&
                              gradient_cut, tau_cut,&
                              error=error)
       CALL xc_rho_set_create(orb_rho_set_h, bounds, density_cut,&
                              gradient_cut, tau_cut,&
                              error=error)
       CALL xc_rho_set_create(orb_rho_set_s, bounds, density_cut,&
                              gradient_cut, tau_cut,&
                              error=error)

       CALL xc_rho_cflags_setall(needs, .FALSE., error=error)
       IF (lsd) THEN
          CALL xb88_lsd_info(needs=needs,error=error)
          needs%norm_drho = .TRUE.
       ELSE
          CALL xb88_lda_info(needs=needs,error=error)
       END IF
       CALL xc_rho_set_atom_update(rho_set_h, needs, nspins, bounds)
       CALL xc_rho_set_atom_update(rho_set_s, needs, nspins, bounds)
       CALL xc_rho_cflags_setall(needs_orbs, .FALSE., error=error)
       needs_orbs%rho = .TRUE.
       IF (lsd) needs_orbs%rho_spin = .TRUE.
       CALL xc_rho_set_atom_update(orb_rho_set_h, needs, nspins, bounds)
       CALL xc_rho_set_atom_update(orb_rho_set_s, needs, nspins, bounds)

       ALLOCATE(rho_h(1:na,1:nspins), rho_s(1:na,1:nspins), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(weight(1:na,1:nr), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(vxc_LB_h(1:na,1:nr,1:nspins), vxc_LB_s(1:na,1:nr,1:nspins), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(vxc_GLLB_h(1:na,1:nr,1:nspins), vxc_GLLB_s(1:na,1:nr,1:nspins), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(vxc_tmp_h(1:na,1:nr,1:nspins), vxc_tmp_s(1:na,1:nr,1:nspins), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(vxc_SAOP_h(1:na,1:nr,1:nspins), vxc_SAOP_s(1:na,1:nr,1:nspins), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(drho_h(1:4,1:na,1:nr,1:nspins), drho_s(1:4,1:na,1:nr,1:nspins), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

!      Distribute the atoms of this kind
       num_pe = qs_env%para_env%num_pe
       mepos  = qs_env%para_env%mepos
       bo     = get_limit( natom, num_pe, mepos )

       DO iat = 1,natom !bo(1),bo(2)
          iatom = atom_list(iat)

          rho_atom => rho_atom_set(iatom)
          NULLIFY(r_h,r_s,dr_h,dr_s,r_h_d,r_s_d)
          CALL get_rho_atom(rho_atom=rho_atom,rho_rad_h=r_h,&
                            rho_rad_s=r_s,drho_rad_h=dr_h,&
                            drho_rad_s=dr_s,rho_rad_h_d=r_h_d,&
                            rho_rad_s_d=r_s_d)
          drho_h = 0.0_dp
          drho_s = 0.0_dp

          DO ir=1, nr
             CALL calc_rho_angular(atomic_grid, harmonics, nspins, .TRUE., &
                                   ir, r_h, r_s, rho_h, rho_s, &
                                   dr_h, dr_s, r_h_d, r_s_d, drho_h, drho_s, error)
             CALL fill_rho_set(rho_set_h,lsd,nspins,needs,rho_h,drho_h,tau,na,ir,error)
             CALL fill_rho_set(rho_set_s,lsd,nspins,needs,rho_s,drho_s,tau,na,ir,error)
             DO ia=1, na
                weight(ia,ir) = atomic_grid%wr(ir)*atomic_grid%wa(ia)
             END DO
          END DO

          !-----------------------------!
          ! 1. Slater exchange for LB94 !
          !-----------------------------!
    xc_fun_section_orig => section_vals_get_subs_vals(xc_section_orig,&
         "XC_FUNCTIONAL",error=error)
    CALL section_vals_create(xc_fun_section_tmp,xc_fun_section_orig%section,&
         error=error)
    CALL section_vals_val_set(xc_fun_section_tmp,"_SECTION_PARAMETERS_",&
         i_val=xc_funct_no_shortcut,error=error)
    CALL section_vals_val_set(xc_fun_section_tmp,"XALPHA%_SECTION_PARAMETERS_",&
         l_val=.TRUE.,error=error)
    CALL section_vals_set_subs_vals(xc_section_tmp,"XC_FUNCTIONAL",&
         xc_fun_section_tmp,error=error)

          !---------------------!
          ! Both: hard and soft !
          !---------------------!
          CALL xc_dset_zero_all(deriv_set, error)
          CALL vxc_of_r_new(xc_fun_section_tmp, rho_set_h, deriv_set, 1, needs,&
                            weight, lsd, na, nr, exc, vxc_tmp_h, vxg, vtau, error=error)
          CALL xc_dset_zero_all(deriv_set, error)
          CALL vxc_of_r_new(xc_fun_section_tmp, rho_set_s, deriv_set, 1, needs,&
                            weight, lsd, na, nr, exc, vxc_tmp_s, vxg, vtau, error=error)

          !-------------------------------------------!
          ! 2. PZ correlation for LB94 and ec_uniform !
          !-------------------------------------------!
    CALL section_vals_val_set(xc_fun_section_tmp,"XALPHA%_SECTION_PARAMETERS_",&
         l_val=.FALSE.,error=error)
    CALL section_vals_val_set(xc_fun_section_tmp,"PZ81%_SECTION_PARAMETERS_",&
         l_val=.TRUE.,error=error)

          !------!
          ! Hard !
          !------!
          CALL xc_dset_zero_all(deriv_set, error)
          CALL vxc_of_r_new(xc_fun_section_tmp, rho_set_h, deriv_set, 1, needs,&
                            weight, lsd, na, nr, exc, vxc_LB_h, vxg, vtau, error=error)
          vxc_LB_h = vxc_LB_h + alpha*vxc_tmp_h
          DO ispin=1, nspins
             dummy => vxc_tmp_h(:,:,ispin:ispin)
             CALL add_lb_pot(dummy, rho_set_h, lsd, ispin, error)
             vxc_LB_h(:,:,ispin) = vxc_LB_h(:,:,ispin) - weight(:,:)*vxc_tmp_h(:,:,ispin)
          END DO
          NULLIFY(dummy)

          vxc_GLLB_h = 0.0_dp
          deriv => xc_dset_get_derivative(deriv_set, "", error=error)
          CPPostcondition(ASSOCIATED(deriv),cp_failure_level,routineP,error,failure)
          CALL xc_derivative_get(deriv,deriv_data=e_uniform,error=error)
          DO ispin=1, nspins
             dummy => vxc_GLLB_h(:,:,ispin:ispin)
             CALL calc_2excpbe(dummy, rho_set_h, e_uniform, lsd)
             vxc_GLLB_h(:,:,ispin) = vxc_GLLB_h(:,:,ispin)*weight(:,:)
          END DO
          NULLIFY(deriv, dummy, e_uniform)

          !------!
          ! Soft !
          !------!
          CALL xc_dset_zero_all(deriv_set, error)
          CALL vxc_of_r_new(xc_fun_section_tmp, rho_set_s, deriv_set, 1, needs,&
                            weight, lsd, na, nr, exc, vxc_LB_s, vxg, vtau, error=error)

          vxc_LB_s = vxc_LB_s + alpha*vxc_tmp_s
          DO ispin=1, nspins
             dummy => vxc_tmp_s(:,:,ispin:ispin)
             CALL add_lb_pot(dummy, rho_set_s, lsd, ispin, error)
             vxc_LB_s(:,:,ispin) = vxc_LB_s(:,:,ispin) - weight(:,:)*vxc_tmp_s(:,:,ispin)
          END DO
          NULLIFY(dummy)

          vxc_GLLB_s = 0.0_dp
          deriv => xc_dset_get_derivative(deriv_set, "", error=error)
          CPPostcondition(ASSOCIATED(deriv),cp_failure_level,routineP,error,failure)
          CALL xc_derivative_get(deriv,deriv_data=e_uniform,error=error)
          DO ispin=1, nspins
             dummy => vxc_GLLB_s(:,:,ispin:ispin)
             CALL calc_2excpbe(dummy, rho_set_s, e_uniform, lsd)
             vxc_GLLB_s(:,:,ispin) = vxc_GLLB_s(:,:,ispin)*weight(:,:)
          END DO
          NULLIFY(deriv, dummy, e_uniform)

          !------------------!
          ! Now the orbitals !
          !------------------!
          vxc_tmp_h = 0.0_dp ; vxc_tmp_s = 0.0_dp

          DO ispin=1, nspins

             DO orb=1, homo(ispin)-1

                ALLOCATE(coeff_col(nrow(ispin),1), stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

                efac = K_rho * SQRT(mo_eigenvalues(ispin)%array(homo(ispin)) - &
                                    mo_eigenvalues(ispin)%array(orb))
                IF (.not.lsd) efac = 2.0_dp * efac

                CALL cp_fm_set_all(single_mo_coeff(ispin)%matrix, 0.0_dp,error=error)
                CALL cp_fm_get_submatrix(mo_coeff(ispin)%matrix, coeff_col, &
                                         1, orb, nrow(ispin), 1, error=error)
                CALL cp_fm_set_submatrix(single_mo_coeff(ispin)%matrix, coeff_col, &
                                         1, orb, error=error)
                CALL cp_dbcsr_set(orbital_density_matrix(ispin)%matrix,0.0_dp,error=error)
                CALL cp_dbcsr_plus_fm_fm_t(orbital_density_matrix(ispin)%matrix,&
                                        matrix_v=single_mo_coeff(ispin)%matrix,&
                                        matrix_g=single_mo_coeff(ispin)%matrix,&
                                        ncol=ncol(ispin), &
                                        alpha=1.0_dp,error=error)

                DEALLOCATE(coeff_col, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

                ! This calculates the CPC and density on the grids for every atom even though
                ! we need it only for iatom at the moment. It seems that to circumvent this,
                ! the routines must be adapted to calculate just iatom
                CALL calculate_rho_atom_coeff(qs_env, orbital_density_matrix, &
                                              local_rho_set%rho_atom_set,error=error)
                CALL prepare_gapw_den(qs_env,local_rho_set,.FALSE.,error)

                rho_atom => local_rho_set%rho_atom_set(iatom)
                NULLIFY(r_h,r_s,dr_h,dr_s,r_h_d,r_s_d)
                CALL get_rho_atom(rho_atom=rho_atom, &
                                  rho_rad_h=r_h, rho_rad_s=r_s)
                DO ir=1, nr
                   CALL calc_rho_angular(atomic_grid, harmonics, nspins, .FALSE., &
                                         ir, r_h, r_s, rho_h, rho_s, &
                                         dr_h, dr_s, r_h_d, r_s_d, drho_h, drho_s, error)
                   CALL fill_rho_set(orb_rho_set_h,lsd,nspins,needs_orbs,rho_h,drho_h,tau,na,ir,error=error)
                   CALL fill_rho_set(orb_rho_set_s,lsd,nspins,needs_orbs,rho_s,drho_s,tau,na,ir,error=error)
                END DO

                IF (lsd) THEN
                   IF (ispin == 1) THEN
                      vxc_tmp_h(:,:,1) = vxc_tmp_h(:,:,1) + efac*orb_rho_set_h%rhoa(:,:,1)
                      vxc_tmp_s(:,:,1) = vxc_tmp_s(:,:,1) + efac*orb_rho_set_s%rhoa(:,:,1)
                   ELSE
                      vxc_tmp_h(:,:,2) = vxc_tmp_h(:,:,2) + efac*orb_rho_set_h%rhob(:,:,1)
                      vxc_tmp_s(:,:,2) = vxc_tmp_s(:,:,2) + efac*orb_rho_set_s%rhob(:,:,1)
                   END IF
                ELSE
                   vxc_tmp_h(:,:,1) = vxc_tmp_h(:,:,1) + efac*orb_rho_set_h%rho(:,:,1)
                   vxc_tmp_s(:,:,1) = vxc_tmp_s(:,:,1) + efac*orb_rho_set_s%rho(:,:,1)
                END IF

             END DO ! orb

          END DO ! ispin

          IF (lsd) THEN
             DO ir=1, nr
                DO ia=1, na
                   IF (rho_set_h%rhoa(ia,ir,1) > rho_set_h%rho_cutoff) &
                        vxc_GLLB_h(ia,ir,1) = vxc_GLLB_h(ia,ir,1) + &
                                              weight(ia,ir)*vxc_tmp_h(ia,ir,1)/rho_set_h%rhoa(ia,ir,1)
                   IF (rho_set_h%rhob(ia,ir,1) > rho_set_h%rho_cutoff) &
                        vxc_GLLB_h(ia,ir,2) = vxc_GLLB_h(ia,ir,2) + &
                                              weight(ia,ir)*vxc_tmp_h(ia,ir,2)/rho_set_h%rhob(ia,ir,1)
                   IF (rho_set_s%rhoa(ia,ir,1) > rho_set_s%rho_cutoff) &
                        vxc_GLLB_s(ia,ir,1) = vxc_GLLB_s(ia,ir,1) + &
                                              weight(ia,ir)*vxc_tmp_s(ia,ir,1)/rho_set_s%rhoa(ia,ir,1)
                   IF (rho_set_s%rhob(ia,ir,1) > rho_set_s%rho_cutoff) &
                        vxc_GLLB_s(ia,ir,2) = vxc_GLLB_s(ia,ir,2) + &
                                              weight(ia,ir)*vxc_tmp_s(ia,ir,2)/rho_set_s%rhob(ia,ir,1)
                END DO
             END DO
          ELSE
             DO ir=1, nr
                DO ia=1, na
                   IF (rho_set_h%rho(ia,ir,1) > rho_set_h%rho_cutoff) &
                        vxc_GLLB_h(ia,ir,1) = vxc_GLLB_h(ia,ir,1) + &
                                              weight(ia,ir)*vxc_tmp_h(ia,ir,1)/rho_set_h%rho(ia,ir,1)
                   IF (rho_set_s%rho(ia,ir,1) > rho_set_s%rho_cutoff) &
                        vxc_GLLB_s(ia,ir,1) = vxc_GLLB_s(ia,ir,1) + &
                                              weight(ia,ir)*vxc_tmp_s(ia,ir,1)/rho_set_s%rho(ia,ir,1)
                END DO
             END DO
          END IF

          vxc_SAOP_h = 0.0_dp ; vxc_SAOP_s = 0.0_dp

          DO ispin=1, nspins

             DO orb=1, homo(ispin)

                ALLOCATE(coeff_col(nrow(ispin),1), stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

                we_LB   = EXP(-2.0_dp*(mo_eigenvalues(ispin)%array(homo(ispin)) - &
                                       mo_eigenvalues(ispin)%array(orb))**2)
                we_GLLB = 1.0_dp - we_LB
                IF (.not.lsd) THEN
                   we_LB   = 2.0_dp * we_LB
                   we_GLLB = 2.0_dp * we_GLLB
                END IF

                vxc_tmp_h(:,:,ispin) = we_LB*vxc_LB_h(:,:,ispin) + &
                                       we_GLLB*vxc_GLLB_h(:,:,ispin)
                vxc_tmp_s(:,:,ispin) = we_LB*vxc_LB_s(:,:,ispin) + &
                                       we_GLLB*vxc_GLLB_s(:,:,ispin)

                CALL cp_fm_set_all(single_mo_coeff(ispin)%matrix, 0.0_dp,error=error)
                CALL cp_fm_get_submatrix(mo_coeff(ispin)%matrix, coeff_col, &
                                         1, orb, nrow(ispin), 1, error=error)
                CALL cp_fm_set_submatrix(single_mo_coeff(ispin)%matrix, coeff_col, &
                                         1, orb, error=error)
                CALL cp_dbcsr_set(orbital_density_matrix(ispin)%matrix,0.0_dp,error=error)
                CALL cp_dbcsr_plus_fm_fm_t(orbital_density_matrix(ispin)%matrix,&
                                        matrix_v=single_mo_coeff(ispin)%matrix,&
                                        matrix_g=single_mo_coeff(ispin)%matrix,&
                                        ncol=ncol(ispin), &
                                        alpha=1.0_dp,error=error)

                DEALLOCATE(coeff_col, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

                ! This calculates the CPC and density on the grids for every atom even though
                ! we need it only for iatom at the moment. It seems that to circumvent this,
                ! the routines must be adapted to calculate just iatom
                CALL calculate_rho_atom_coeff(qs_env, orbital_density_matrix, &
                                              local_rho_set%rho_atom_set,error=error)
                CALL prepare_gapw_den(qs_env,local_rho_set,.FALSE.,error)

                rho_atom => local_rho_set%rho_atom_set(iatom)
                NULLIFY(r_h,r_s,dr_h,dr_s,r_h_d,r_s_d)
                CALL get_rho_atom(rho_atom=rho_atom, &
                                  rho_rad_h=r_h, rho_rad_s=r_s)
                DO ir=1, nr
                   CALL calc_rho_angular(atomic_grid, harmonics, nspins, .FALSE., &
                                         ir, r_h, r_s, rho_h, rho_s, &
                                         dr_h, dr_s, r_h_d, r_s_d, drho_h, drho_s, error)
                   CALL fill_rho_set(orb_rho_set_h,lsd,nspins,needs_orbs,rho_h,drho_h,tau,na,ir,error=error)
                   CALL fill_rho_set(orb_rho_set_s,lsd,nspins,needs_orbs,rho_s,drho_s,tau,na,ir,error=error)
                END DO

                IF (lsd) THEN
                   IF (ispin == 1) THEN
                      vxc_SAOP_h(:,:,1) = vxc_SAOP_h(:,:,1) + vxc_tmp_h(:,:,1)*orb_rho_set_h%rhoa(:,:,1)
                      vxc_SAOP_s(:,:,1) = vxc_SAOP_s(:,:,1) + vxc_tmp_s(:,:,1)*orb_rho_set_s%rhoa(:,:,1)
                   ELSE
                      vxc_SAOP_h(:,:,2) = vxc_SAOP_h(:,:,2) + vxc_tmp_h(:,:,2)*orb_rho_set_h%rhob(:,:,1)
                      vxc_SAOP_s(:,:,2) = vxc_SAOP_s(:,:,2) + vxc_tmp_s(:,:,2)*orb_rho_set_s%rhob(:,:,1)
                   END IF
                ELSE
                   vxc_SAOP_h(:,:,1) = vxc_SAOP_h(:,:,1) + vxc_tmp_h(:,:,1)*orb_rho_set_h%rho(:,:,1)
                   vxc_SAOP_s(:,:,1) = vxc_SAOP_s(:,:,1) + vxc_tmp_s(:,:,1)*orb_rho_set_s%rho(:,:,1)
                END IF

             END DO ! orb

          END DO ! ispin

          IF (lsd) THEN
             DO ir=1, nr
                DO ia=1, na
                   IF (rho_set_h%rhoa(ia,ir,1) > rho_set_h%rho_cutoff) THEN
                      vxc_SAOP_h(ia,ir,1) = vxc_SAOP_h(ia,ir,1) / rho_set_h%rhoa(ia,ir,1)
                   ELSE
                      vxc_SAOP_h(ia,ir,1) = 0.0_dp
                   END IF
                   IF (rho_set_h%rhob(ia,ir,1) > rho_set_h%rho_cutoff) THEN
                      vxc_SAOP_h(ia,ir,2) = vxc_SAOP_h(ia,ir,2) / rho_set_h%rhob(ia,ir,1)
                   ELSE
                      vxc_SAOP_h(ia,ir,2) = 0.0_dp
                   END IF
                   IF (rho_set_s%rhoa(ia,ir,1) > rho_set_s%rho_cutoff) THEN
                      vxc_SAOP_s(ia,ir,1) = vxc_SAOP_s(ia,ir,1) / rho_set_s%rhoa(ia,ir,1)
                   ELSE
                      vxc_SAOP_s(ia,ir,1) = 0.0_dp
                   END IF
                   IF (rho_set_s%rhob(ia,ir,1) > rho_set_s%rho_cutoff) THEN
                      vxc_SAOP_s(ia,ir,2) = vxc_SAOP_s(ia,ir,2) / rho_set_s%rhob(ia,ir,1)
                   ELSE
                      vxc_SAOP_s(ia,ir,2) = 0.0_dp
                   END IF
                END DO
             END DO
          ELSE
             DO ir=1, nr
                DO ia=1, na
                   IF (rho_set_h%rho(ia,ir,1) > rho_set_h%rho_cutoff) THEN
                      vxc_SAOP_h(ia,ir,1) = vxc_SAOP_h(ia,ir,1) / rho_set_h%rho(ia,ir,1)
                   ELSE
                      vxc_SAOP_h(ia,ir,1) = 0.0_dp
                   END IF
                   IF (rho_set_s%rho(ia,ir,1) > rho_set_s%rho_cutoff) THEN
                      vxc_SAOP_s(ia,ir,1) = vxc_SAOP_s(ia,ir,1) / rho_set_s%rho(ia,ir,1)
                   ELSE
                      vxc_SAOP_s(ia,ir,1) = 0.0_dp
                   END IF
                END DO
             END DO
          END IF

          rho_atom => rho_atom_set(iatom)
          atom_kind => atomic_kind_set(ikind)
          SELECT CASE (oe_corr)
          CASE(oe_lb)
             CALL gaVxcgb_noGC(vxc_LB_h,vxc_LB_s,atom_kind,rho_atom,nspins,error=error)
          CASE (oe_gllb)
             CALL gaVxcgb_noGC(vxc_GLLB_h,vxc_GLLB_s,atom_kind,rho_atom,nspins,error=error)
          CASE (oe_saop)
             CALL gaVxcgb_noGC(vxc_SAOP_h,vxc_SAOP_s,atom_kind,rho_atom,nspins,error=error)
          CASE default
             CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
          END SELECT

       END DO

       DEALLOCATE(rho_h, rho_s, weight, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(vxc_LB_h, vxc_LB_s, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(vxc_GLLB_h, vxc_GLLB_s, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(vxc_tmp_h, vxc_tmp_s, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(vxc_SAOP_h, vxc_SAOP_s, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(drho_h, drho_s, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       CALL xc_dset_release(deriv_set, error=error)
       CALL xc_rho_set_release(rho_set_h, error=error)
       CALL xc_rho_set_release(rho_set_s, error=error)
       CALL xc_rho_set_release(orb_rho_set_h, error=error)
       CALL xc_rho_set_release(orb_rho_set_s, error=error)

    END DO

    CALL update_ks_atom(qs_env, qs_env%matrix_ks,rho_structure%rho_ao,.FALSE.,error=error)

    !---------!
    ! Cleanup !
    !---------!
    CALL section_vals_release(xc_fun_section_tmp,error=error)
    CALL section_vals_release(xc_section_tmp,error=error)
    CALL section_vals_release(xc_section_orig,error=error)

    CALL local_rho_set_release(local_rho_set,error)
    DO ispin=1, nspins
       CALL cp_fm_release(single_mo_coeff(ispin)%matrix, error=error)
    END DO
    DEALLOCATE(mo_coeff, single_mo_coeff, mo_eigenvalues, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_deallocate_matrix_set(orbital_density_matrix,error=error)

  END SUBROUTINE gapw_add_atomic_saop_pot

! *****************************************************************************
  SUBROUTINE add_lb_pot(pot, rho_set, lsd, spin, error)

    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: pot
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    LOGICAL, INTENT(IN)                      :: lsd
    INTEGER, INTENT(IN)                      :: spin
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    REAL(KIND=dp), PARAMETER                 :: ob3 = 1.0_dp/3.0_dp

    INTEGER                                  :: i, j, k
    INTEGER, DIMENSION(2, 3)                 :: bo
    REAL(KIND=dp)                            :: n, n_13, x, x2

    bo = rho_set%local_bounds

    DO k=bo(1,3), bo(2,3)
       DO j=bo(1,2), bo(2,2)
          DO i=bo(1,1), bo(2,1)
             IF (.NOT.lsd) THEN
                IF (rho_set%rho(i,j,k) > rho_set%rho_cutoff) THEN
                   n    = rho_set%rho(i,j,k)/2.0_dp
                   n_13 = n**ob3
                   x    = (rho_set%norm_drho(i,j,k)/2.0_dp)/(n*n_13)
                   x2   = x*x
                   pot(i,j,k) = beta*x2*n_13/(1.0_dp+3.0_dp*beta*x*LOG(x + SQRT(x2+1.0_dp)))
                END IF
             ELSE
                IF (spin == 1) THEN
                   IF (rho_set%rhoa(i,j,k) > rho_set%rho_cutoff) THEN
                      n_13 = rho_set%rhoa_1_3(i,j,k)
                      x    = rho_set%norm_drhoa(i,j,k)/(rho_set%rhoa(i,j,k)*n_13)
                      x2   = x*x
                      pot(i,j,k) = beta*x2*n_13/(1.0_dp+3.0_dp*beta*x*LOG(SQRT(x2+1.0_dp)+x))
                   END IF
                ELSE IF (spin == 2) THEN
                   IF (rho_set%rhob(i,j,k) > rho_set%rho_cutoff) THEN
                      n_13 = rho_set%rhob_1_3(i,j,k)
                      x    = rho_set%norm_drhob(i,j,k)/(rho_set%rhob(i,j,k)*n_13)
                      x2   = x*x
                      pot(i,j,k) = beta*x2*n_13/(1.0_dp+3.0_dp*beta*x*LOG(SQRT(x2+1.0_dp)+x))
                   END IF
                END IF
             END IF
          END DO
       END DO
    END DO

  END SUBROUTINE add_lb_pot

! *****************************************************************************
  SUBROUTINE calc_2excpbe(pot, rho_set, e_uniform, lsd)

    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: pot
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: e_uniform
    LOGICAL, INTENT(IN)                      :: lsd

    INTEGER                                  :: i, j, k
    INTEGER, DIMENSION(2, 3)                 :: bo
    REAL(KIND=dp)                            :: e_unif, rho

    bo = rho_set%local_bounds

    DO k=bo(1,3), bo(2,3)
       DO j=bo(1,2), bo(2,2)
          DO i=bo(1,1), bo(2,1)
             IF (.NOT.lsd) THEN
                IF (rho_set%rho(i,j,k) > rho_set%rho_cutoff) THEN
                   e_unif = e_uniform(i,j,k)/rho_set%rho(i,j,k)
                ELSE
                   e_unif = 0.0_dp
                END IF
                pot(i,j,k) = &
                     2.0_dp * &
                     calc_ecpbe_r(rho_set%rho(i,j,k), rho_set%norm_drho(i,j,k), &
                                  e_unif, rho_set%rho_cutoff, rho_set%drho_cutoff) + &
                     2.0_dp * &
                     calc_expbe_r(rho_set%rho(i,j,k), rho_set%norm_drho(i,j,k), &
                                  rho_set%rho_cutoff, rho_set%drho_cutoff)
             ELSE
                rho = rho_set%rhoa(i,j,k) + rho_set%rhob(i,j,k)
                IF (rho > rho_set%rho_cutoff) THEN
                   e_unif = e_uniform(i,j,k)/rho
                ELSE
                   e_unif = 0.0_dp
                END IF
                pot(i,j,k) = &
                     2.0_dp * &
                     calc_ecpbe_u(rho_set%rhoa(i,j,k), rho_set%rhob(i,j,k), rho_set%norm_drho(i,j,k), &
                                  e_unif, &
                                  rho_set%rho_cutoff, rho_set%drho_cutoff) + &
                     2.0_dp * &
                     calc_expbe_u(rho_set%rhoa(i,j,k), rho_set%rhob(i,j,k), rho_set%norm_drho(i,j,k), &
                                  rho_set%rho_cutoff, rho_set%drho_cutoff)
              END IF
          END DO
       END DO
    END DO

  END SUBROUTINE calc_2excpbe

! *****************************************************************************
  FUNCTION calc_ecpbe_u(ra, rb, ngr, ec_unif, rc, ngrc) RESULT(res)

    REAL(kind=dp), INTENT(in)                :: ra, rb, ngr, ec_unif, rc, ngrc
    REAL(kind=dp)                            :: res

    REAL(kind=dp), PARAMETER                 :: ob3 = 1.0_dp/3.0_dp, &
                                                tb3 = 2.0_dp/3.0_dp

    REAL(kind=dp)                            :: A, At2, H, kf, kl, ks, phi, &
                                                phi3, r, t2, zeta

    r = ra + rb
    H = 0.0_dp
    IF (r > rc .AND. ngr > ngrc) THEN
       zeta = (ra-rb)/r
       IF (zeta > 1.0_dp) zeta = 1.0_dp ! machine precision problem
       IF (zeta < -1.0_dp) zeta = -1.0_dp ! machine precision problem
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

! *****************************************************************************
  FUNCTION calc_ecpbe_r(r, ngr, ec_unif, rc, ngrc) RESULT(res)

    REAL(kind=dp), INTENT(in)                :: r, ngr, ec_unif, rc, ngrc
    REAL(kind=dp)                            :: res

    REAL(kind=dp), PARAMETER                 :: ob3 = 1.0_dp/3.0_dp

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

! *****************************************************************************
  FUNCTION calc_expbe_u(ra, rb, ngr, rc, ngrc) RESULT(res)

    REAL(kind=dp), INTENT(in)                :: ra, rb, ngr, rc, ngrc
    REAL(kind=dp)                            :: res

    REAL(kind=dp)                            :: r

    r = ra + rb
    res = calc_expbe_r(r, ngr, rc, ngrc)

  END FUNCTION calc_expbe_u

! *****************************************************************************
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

END MODULE xc_pot_saop
