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
  
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE coefficient_types,               ONLY: coeff_type
  USE cp_array_r_utils,                ONLY: cp_1d_r_p_type
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
  USE input_constants,                 ONLY: oe_lb, oe_gllb, oe_saop
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
  USE simulation_cell,                 ONLY: cell_type
  USE sparse_matrix_types,             ONLY: deallocate_matrix,&
                                             real_matrix_p_type,&
                                             replicate_matrix,&
                                             set_matrix
  USE util,                            ONLY: get_limit
  USE xc,                              ONLY: xc_vxc_pw_create
  USE xc_atom,                         ONLY: vxc_of_r_new,&
                                             xc_rho_set_atom_update,&
                                             fill_rho_set
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                             xc_dset_create,&
                                             xc_dset_get_derivative,&
                                             xc_dset_release,&
                                             xc_dset_zero_all
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

  SUBROUTINE add_saop_pot(ks_matrix, qs_env, oe_corr, error)

    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: oe_corr
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "add_saop_pot", &
      routineP = moduleN//"/"//routineN

    INTEGER                                  :: homo, i, ispin, j, k, nspins, &
                                                orb, stat
    INTEGER, DIMENSION(2)                    :: ncol, nrow
    INTEGER, DIMENSION(2, 3)                 :: bo
    LOGICAL                                  :: failure, ionode, lsd, gapw
    REAL(KIND=dp)                            :: efac, tot_rho_psi, we_GLLB, &
                                                we_LB, xc_energy
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: coeff_col
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mo_eigenvalues
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: e_uniform, dummy
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
      POINTER                                :: orbital_density_matrix
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
    NULLIFY(orbital_density_matrix, dummy)
    logger => cp_error_get_logger(error)
    ionode = (logger%para_env%mepos==logger%para_env%source)

    IF (ionode) WRITE (*, *) " in "//routineP

    CALL get_qs_env(qs_env=qs_env,&
                    rho=rho_struct,&
                    pw_env=pw_env,&
                    cell=cell,&
                    dft_control=dft_control,&
                    mos=molecular_orbitals)
    xc_control = dft_control%xc_control
    gapw = dft_control%qs_control%gapw
    CALL pw_env_get(pw_env,auxbas_pw_pool=auxbas_pw_pool,error=error)

    nspins = dft_control%nspins
    lsd = (nspins>1)
    
    ALLOCATE(rho_r(nspins), orbital_density_matrix(nspins), single_mo_coeff(nspins), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO ispin=1,nspins
       rho_r(ispin)%pw => rho_struct%rho_r(ispin)%pw
       NULLIFY(orbital_density_matrix(ispin)%matrix)
       CALL replicate_matrix(rho_struct%rho_ao(ispin)%matrix, &
                             orbital_density_matrix(ispin)%matrix, &
                             "orbital density", .TRUE.)
    END DO
    bo = rho_r(1)%pw%pw_grid%bounds_local

    !---------------------------!
    ! create the density needed !
    !---------------------------!
    CALL xc_rho_set_create(rho_set, bo, &
                           dft_control%xc_control%density_cut, &
                           dft_control%xc_control%gradient_cut, &
                           dft_control%xc_control%tau_cut, &
                           error=error)
    CALL xc_rho_cflags_setall(needs,.FALSE.,error)
    IF (lsd) THEN
       CALL xb88_lsd_info(needs=needs)
       needs%norm_drho = .TRUE.
    ELSE
       CALL xb88_lda_info(needs=needs)
    END IF
    CALL xc_rho_set_update(rho_set, rho_r, rho_g, tau, needs, &
                           dft_control%xc_control, cell, auxbas_pw_pool, error)

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
       dummy => vxc_GLLB(ispin)%pw%cr3d
       CALL calc_2excpbe(dummy, rho_set, e_uniform, lsd)
    END DO
    NULLIFY(dummy)

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
          CALL set_matrix(orbital_density_matrix(ispin)%matrix, value=0.0_dp)
          CALL cp_sm_plus_fm_fm_t(orbital_density_matrix(ispin)%matrix,&
                                  matrix_v=single_mo_coeff(ispin)%matrix,&
                                  matrix_g=single_mo_coeff(ispin)%matrix,&
                                  ncol=ncol(ispin), &
                                  alpha=1.0_dp)
          CALL pw_zero(orbital%pw)
          CALL pw_zero(orbital_g%pw)
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
          CALL set_matrix(orbital_density_matrix(ispin)%matrix, value=0.0_dp)
          CALL cp_sm_plus_fm_fm_t(orbital_density_matrix(ispin)%matrix,&
                                  matrix_v=single_mo_coeff(ispin)%matrix,&
                                  matrix_g=single_mo_coeff(ispin)%matrix,&
                                  ncol=ncol(ispin), &
                                  alpha=1.0_dp)
          CALL pw_zero(orbital%pw)
          CALL pw_zero(orbital_g%pw)
          CALL calculate_rho_elec(orbital_density_matrix(ispin)%matrix, &
                                  orbital, orbital_g, &
                                  tot_rho_psi, qs_env, error=error)

!TC          CALL calculate_wavefunction(mo_coeff,orb,orbital2,orbital2_g,qs_env,error)
!TC          write (*,*) orb, maxval(abs(orbital2%pw%cr3d-orbital%pw%cr3d))
          
          vxc_SAOP(ispin)%pw%cr3d = vxc_SAOP(ispin)%pw%cr3d + &
               orbital%pw%cr3d * vxc_tmp(ispin)%pw%cr3d

       END DO

       CALL cp_fm_release(single_mo_coeff(ispin)%matrix, error=error)
       CALL deallocate_matrix(orbital_density_matrix(ispin)%matrix)

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

    !--------------------!
    ! Do the integration !
    !--------------------!
    DO ispin=1, nspins

       if (oe_corr == oe_lb) then
          vxc_SAOP(ispin)%pw%cr3d = vxc_LB(ispin)%pw%cr3d
       else if (oe_corr == oe_gllb) then
          vxc_SAOP(ispin)%pw%cr3d = vxc_GLLB(ispin)%pw%cr3d
       end if
       vxc_SAOP(ispin)%pw%cr3d = vxc_SAOP(ispin)%pw%cr3d * vxc_SAOP(ispin)%pw%pw_grid%dvol

       CALL integrate_v_rspace(vxc_SAOP(ispin), rho_struct%rho_ao(ispin), &
                               ks_matrix(ispin), qs_env, &
                               calculate_forces=.FALSE.,&
                               gapw=gapw, error=error)

    END DO

    DO ispin=1, nspins
       CALL pw_pool_give_back_coeff(auxbas_pw_pool, vxc_SAOP(ispin), error=error)
       CALL pw_pool_give_back_pw(auxbas_pw_pool, vxc_GLLB(ispin)%pw, error=error)
       CALL pw_release(vxc_LB(ispin)%pw)
       CALL pw_release(vxc_tmp(ispin)%pw)
    END DO
    DEALLOCATE(vxc_GLLB, vxc_LB, vxc_tmp, orbital_density_matrix, single_mo_coeff, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    
    DEALLOCATE(rho_r, vxc_SAOP, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    !----------------------------------!
    ! Restore the values of xc_control !
    !----------------------------------!
    dft_control%xc_control = xc_control

    !-----------------------!
    ! Call the GAPW routine !
    !-----------------------!
    IF (gapw) THEN
       CALL gapw_add_atomic_saop_pot(ks_matrix, qs_env, oe_corr, error)
    END IF

  END SUBROUTINE add_saop_pot

  !*****************************************************************

  SUBROUTINE gapw_add_atomic_saop_pot(ks_matrix, qs_env, oe_corr, error)

    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: oe_corr
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'gapw_add_atomic_saop_pot', &
                                                routineP = moduleN//'/'//routineN

    INTEGER                                  :: ikind, natom, ispin, nspins
    INTEGER                                  :: ia, na, ir, nr, num_pe, mepos, i
    INTEGER                                  :: ifunc, stat, iat, iatom
    INTEGER                                  :: orb
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    INTEGER, DIMENSION(2, 3)                 :: bounds
    INTEGER, DIMENSION(2)                    :: bo, ncol, nrow, homo

    LOGICAL                                  :: paw_atom, lsd, failure
    
    REAL(KIND=dp)                            :: exc, efac, we_LB, we_GLLB
    TYPE(cp_1d_r_p_type), DIMENSION(:), &
         POINTER                             :: mo_eigenvalues
    REAL(KIND=dp), DIMENSION(:, :), &
         ALLOCATABLE                         :: coeff_col
    REAL(KIND=dp), DIMENSION(:,:), &
         POINTER                             :: weight
    REAL(KIND=dp), DIMENSION(:,:), &
         POINTER                             :: rho_h, rho_s
    REAL(KIND=dp), DIMENSION(:,:,:), &
         POINTER                             :: vxc_LB_h, vxc_LB_s, &
                                                vxc_GLLB_h, vxc_GLLB_s, &
                                                vxc_tmp_h, vxc_tmp_s, &
                                                vxc_SAOP_h, vxc_SAOP_s, &
                                                dummy, &
                                                e_uniform
    REAL(KIND=dp), DIMENSION(:,:,:,:), &
         POINTER                             :: drho_h, drho_s, vxg

    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(atomic_kind_type), DIMENSION(:),&
         POINTER                             :: atomic_kind_set
    TYPE(cp_fm_p_type), DIMENSION(:), &
         POINTER                             :: mo_coeff, single_mo_coeff
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(grid_atom_type), POINTER            :: atomic_grid
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    TYPE(local_rho_type), POINTER            :: local_rho_set
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: molecular_orbitals
    TYPE(qs_rho_type), POINTER               :: rho_structure
    TYPE(real_matrix_p_type), DIMENSION(:), &
         POINTER                             :: orbital_density_matrix
    TYPE(rho_atom_coeff),DIMENSION(:),& 
         POINTER                             :: r_h, r_s, dr_h, dr_s 
    TYPE(rho_atom_coeff),DIMENSION(:,:),&
         POINTER                             :: r_h_d, r_s_d
    TYPE(rho_atom_type), DIMENSION(:), &
         POINTER                             :: rho_atom_set
    TYPE(rho_atom_type),  POINTER            :: rho_atom
    TYPE(xc_control_type)                    :: xc_control_sav
    TYPE(xc_control_type), POINTER           :: xc_control
    TYPE(xc_derivative_set_type), &
         POINTER                             :: deriv_set
    TYPE(xc_derivative_type), POINTER        :: deriv
    TYPE(xc_rho_cflags_type)                 :: needs, needs_orbs
    TYPE(xc_rho_set_type), POINTER           :: rho_set_h, rho_set_s, orb_rho_set_h, orb_rho_set_s
    
    failure = .false.
    
    nullify(weight, rho_h, rho_s, vxc_LB_h, vxc_LB_s, vxc_GLLB_h, vxc_GLLB_s, &
            vxc_tmp_h, vxc_tmp_s, dummy, e_uniform, drho_h, drho_s, vxg, atom_list, &
            atom_kind, atomic_kind_set, dft_control, deriv_set, deriv, atomic_grid, &
            harmonics, molecular_orbitals, rho_structure, r_h, r_s, dr_h, dr_s, &
            r_h_d, r_s_d, rho_atom_set, rho_atom, rho_set_h, rho_set_s, xc_control, &
            mo_coeff, single_mo_coeff, mo_eigenvalues, local_rho_set, &
            orbital_density_matrix, orb_rho_set_h, orb_rho_set_s, vxc_SAOP_h, vxc_SAOP_s)

    call get_qs_env(qs_env, dft_control=dft_control, &
                            rho=rho_structure, &
                            mos=molecular_orbitals, &
                            atomic_kind_set=atomic_kind_set, & 
                            rho_atom_set=rho_atom_set)

    xc_control => dft_control%xc_control
    xc_control_sav = xc_control
    nspins     = dft_control%nspins
    lsd        = (nspins > 1)

    call calculate_rho_atom_coeff(qs_env, rho_structure%rho_ao)
    call prepare_gapw_den(qs_env, error=error)

    allocate(mo_coeff(nspins), single_mo_coeff(nspins), &
             mo_eigenvalues(nspins),&
             orbital_density_matrix(nspins), &
             stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    do ispin=1, nspins
       call get_mo_set(molecular_orbitals(ispin)%mo_set, &
                       mo_coeff=mo_coeff(ispin)%matrix, &
                       eigenvalues=mo_eigenvalues(ispin)%array, &
                       homo=homo(ispin))
       call cp_fm_create(single_mo_coeff(ispin)%matrix, &
                         mo_coeff(ispin)%matrix%matrix_struct, &
                         "orbital density matrix", error)
       call cp_fm_get_info(single_mo_coeff(ispin)%matrix, &
                           nrow_global=nrow(ispin), ncol_global=ncol(ispin))
       nullify(orbital_density_matrix(ispin)%matrix)
       call replicate_matrix(rho_structure%rho_ao(ispin)%matrix, &
                             orbital_density_matrix(ispin)%matrix, &
                             "orbital density", .TRUE.)
    end do
    call local_rho_set_create(local_rho_set, error=error)
    call allocate_rho_atom_internals(qs_env, local_rho_set%rho_atom_set)

    do ikind=1, size(atomic_kind_set)

       atom_kind => atomic_kind_set(ikind)
       call get_atomic_kind(atomic_kind=atom_kind, atom_list=atom_list,&
                            natom=natom, paw_atom=paw_atom,&
                            harmonics=harmonics, grid_atom=atomic_grid)

       if(.not. paw_atom) cycle

       nr = atomic_grid%nr
       na = atomic_grid%ng_sphere
       bounds(1:2,1:3) = 1
       bounds(2,1)     = na
       bounds(2,2)     = nr

       call xc_dset_create(deriv_set, local_bounds=bounds, error=error)
       
       call xc_rho_set_create(rho_set_h, bounds, xc_control%density_cut,&
                              xc_control%gradient_cut, xc_control%tau_cut,&
                              error=error)
       call xc_rho_set_create(rho_set_s, bounds, xc_control%density_cut,&
                              xc_control%gradient_cut, xc_control%tau_cut,&
                              error=error)
       call xc_rho_set_create(orb_rho_set_h, bounds, xc_control%density_cut,&
                              xc_control%gradient_cut, xc_control%tau_cut,&
                              error=error)
       call xc_rho_set_create(orb_rho_set_s, bounds, xc_control%density_cut,&
                              xc_control%gradient_cut, xc_control%tau_cut,&
                              error=error)

       call xc_rho_cflags_setall(needs, .false., error=error)
       if (lsd) then
          call xb88_lsd_info(needs=needs)
          needs%norm_drho = .TRUE.
       else
          call xb88_lda_info(needs=needs)
       end if
       call xc_rho_set_atom_update(rho_set_h, needs, nspins, bounds)
       call xc_rho_set_atom_update(rho_set_s, needs, nspins, bounds)
       call xc_rho_cflags_setall(needs_orbs, .false., error=error)
       needs_orbs%rho = .true.
       if (lsd) needs_orbs%rho_spin = .true.
       call xc_rho_set_atom_update(orb_rho_set_h, needs, nspins, bounds)
       call xc_rho_set_atom_update(orb_rho_set_s, needs, nspins, bounds)
       
       allocate(rho_h(1:na,1:nspins), rho_s(1:na,1:nspins), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       allocate(weight(1:na,1:nr), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       allocate(vxc_LB_h(1:na,1:nr,1:nspins), vxc_LB_s(1:na,1:nr,1:nspins), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       allocate(vxc_GLLB_h(1:na,1:nr,1:nspins), vxc_GLLB_s(1:na,1:nr,1:nspins), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       allocate(vxc_tmp_h(1:na,1:nr,1:nspins), vxc_tmp_s(1:na,1:nr,1:nspins), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       allocate(vxc_SAOP_h(1:na,1:nr,1:nspins), vxc_SAOP_s(1:na,1:nr,1:nspins), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       allocate(drho_h(1:4,1:na,1:nr,1:nspins), drho_s(1:4,1:na,1:nr,1:nspins), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       

!      Distribute the atoms of this kind
       num_pe = qs_env%para_env%num_pe
       mepos  = qs_env%para_env%mepos
       bo     = get_limit( natom, num_pe, mepos )

       do iat = 1,natom !bo(1),bo(2) 
          iatom = atom_list(iat)
          
          rho_atom => rho_atom_set(iatom)
          nullify(r_h,r_s,dr_h,dr_s,r_h_d,r_s_d)
          call get_rho_atom(rho_atom=rho_atom,rho_rad_h=r_h,&
                            rho_rad_s=r_s,drho_rad_h=dr_h,&
                            drho_rad_s=dr_s,rho_rad_h_d=r_h_d,&
                            rho_rad_s_d=r_s_d)
          drho_h = 0.0_dp
          drho_s = 0.0_dp

          do ir=1, nr
             call calc_rho_angular(atomic_grid, harmonics, nspins, .true., &
                                   ir, r_h, r_s, rho_h, rho_s, &
                                   dr_h, dr_s, r_h_d, r_s_d, drho_h, drho_s, error)
             call fill_rho_set(rho_set_h,lsd,nspins,needs,rho_h,drho_h,na,ir,error)
             call fill_rho_set(rho_set_s,lsd,nspins,needs,rho_s,drho_s,na,ir,error)
             do ia=1, na
                weight(ia,ir) = atomic_grid%wr(ir)*atomic_grid%wa(ia)
             end do
          end do

          !-----------------------------!
          ! 1. Slater exchange for LB94 !
          !-----------------------------!
          xc_control%functionals = ""
          call xc_get_functional (1, "SLATER", xc_control%functionals, &
                                               xc_control%gradient_functionals, &
                                               xc_control%crossterms )
          do i=1, 3
             if (xc_control%functionals(i) == "") xc_control%functionals(i) = "NONE"
          end do

          !---------------------!
          ! Both: hard and soft !
          !---------------------!
          call xc_dset_zero_all(deriv_set, error)
          call vxc_of_r_new(xc_control, rho_set_h, deriv_set, 1, needs,&
                            weight, lsd, na, nr, exc, vxc_tmp_h, vxg, error=error)
          call xc_dset_zero_all(deriv_set, error)
          call vxc_of_r_new(xc_control, rho_set_s, deriv_set, 1, needs,&
                            weight, lsd, na, nr, exc, vxc_tmp_s, vxg, error=error)


          !-------------------------------------------!
          ! 2. PZ correlation for LB94 and ec_uniform !
          !-------------------------------------------!
          xc_control%functionals = ""
          call xc_get_functional ( 2, "PZ", xc_control%functionals, &
                                            xc_control%gradient_functionals, &
                                            xc_control%crossterms )
          do i=1, 3
             if (xc_control%functionals(i) == "") xc_control%functionals(i) = "NONE"
          end do

          !------!
          ! Hard !
          !------!
          call xc_dset_zero_all(deriv_set, error)
          call vxc_of_r_new(xc_control, rho_set_h, deriv_set, 1, needs,&
                            weight, lsd, na, nr, exc, vxc_LB_h, vxg, error=error)
          vxc_LB_h = vxc_LB_h + alpha*vxc_tmp_h
          do ispin=1, nspins
             dummy => vxc_tmp_h(:,:,ispin:ispin)
             call add_lb_pot(dummy, rho_set_h, lsd, ispin, error)
             vxc_LB_h(:,:,ispin) = vxc_LB_h(:,:,ispin) - weight(:,:)*vxc_tmp_h(:,:,ispin)
          end do
          nullify(dummy)
             
          vxc_GLLB_h = 0.0_dp
          deriv => xc_dset_get_derivative(deriv_set, "", error=error)
          CPPostcondition(associated(deriv),cp_failure_level,routineP,error,failure)
          call xc_derivative_get(deriv,deriv_data=e_uniform,error=error)
          do ispin=1, nspins
             dummy => vxc_GLLB_h(:,:,ispin:ispin)
             call calc_2excpbe(dummy, rho_set_h, e_uniform, lsd)
             vxc_GLLB_h(:,:,ispin) = vxc_GLLB_h(:,:,ispin)*weight(:,:)
          end do
          nullify(deriv, dummy, e_uniform)

          !------!
          ! Soft !
          !------!
          call xc_dset_zero_all(deriv_set, error)
          call vxc_of_r_new(xc_control, rho_set_s, deriv_set, 1, needs,&
                            weight, lsd, na, nr, exc, vxc_LB_s, vxg, error=error)

          vxc_LB_s = vxc_LB_s + alpha*vxc_tmp_s
          do ispin=1, nspins
             dummy => vxc_tmp_s(:,:,ispin:ispin)
             call add_lb_pot(dummy, rho_set_s, lsd, ispin, error)
             vxc_LB_s(:,:,ispin) = vxc_LB_s(:,:,ispin) - weight(:,:)*vxc_tmp_s(:,:,ispin)
          end do
          nullify(dummy)

          vxc_GLLB_s = 0.0_dp
          deriv => xc_dset_get_derivative(deriv_set, "", error=error)
          CPPostcondition(associated(deriv),cp_failure_level,routineP,error,failure)
          call xc_derivative_get(deriv,deriv_data=e_uniform,error=error)
          do ispin=1, nspins
             dummy => vxc_GLLB_s(:,:,ispin:ispin)
             call calc_2excpbe(dummy, rho_set_s, e_uniform, lsd)
             vxc_GLLB_s(:,:,ispin) = vxc_GLLB_s(:,:,ispin)*weight(:,:)
          end do
          nullify(deriv, dummy, e_uniform)


          !------------------!
          ! Now the orbitals !
          !------------------!
          vxc_tmp_h = 0.0_dp ; vxc_tmp_s = 0.0_dp

          do ispin=1, nspins

             do orb=1, homo(ispin)-1
                
                allocate(coeff_col(nrow(ispin),1), stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                
                efac = K_rho * sqrt(mo_eigenvalues(ispin)%array(homo(ispin)) - &
                                    mo_eigenvalues(ispin)%array(orb))
                if (.not.lsd) efac = 2.0_dp * efac
                
                call cp_fm_set_all(single_mo_coeff(ispin)%matrix, 0.0_dp)
                call cp_fm_get_submatrix(mo_coeff(ispin)%matrix, coeff_col, &
                                         1, orb, nrow(ispin), 1, error=error)
                call cp_fm_set_submatrix(single_mo_coeff(ispin)%matrix, coeff_col, &
                                         1, orb, error=error)
                call set_matrix(orbital_density_matrix(ispin)%matrix, value=0.0_dp)
                call cp_sm_plus_fm_fm_t(orbital_density_matrix(ispin)%matrix,&
                                        matrix_v=single_mo_coeff(ispin)%matrix,&
                                        matrix_g=single_mo_coeff(ispin)%matrix,&
                                        ncol=ncol(ispin), &
                                        alpha=1.0_dp)

                deallocate(coeff_col, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)


                ! This calculates the CPC and density on the grids for every atom even though
                ! we need it only for iatom at the moment. It seems that to circumvent this,
                ! the routines must be adapted to calculate just iatom
                call calculate_rho_atom_coeff(qs_env, orbital_density_matrix, &
                                              local_rho_set%rho_atom_set)
                call prepare_gapw_den(qs_env,local_rho_set,.FALSE.,error)
                   
                rho_atom => local_rho_set%rho_atom_set(iatom)
                nullify(r_h,r_s,dr_h,dr_s,r_h_d,r_s_d)
                call get_rho_atom(rho_atom=rho_atom, &
                                  rho_rad_h=r_h, rho_rad_s=r_s)
                do ir=1, nr
                   call calc_rho_angular(atomic_grid, harmonics, nspins, .false., &
                                         ir, r_h, r_s, rho_h, rho_s, &
                                         dr_h, dr_s, r_h_d, r_s_d, drho_h, drho_s, error)
                   call fill_rho_set(orb_rho_set_h,lsd,nspins,needs_orbs,rho_h,drho_h,na,ir,error=error)
                   call fill_rho_set(orb_rho_set_s,lsd,nspins,needs_orbs,rho_s,drho_s,na,ir,error=error)
                end do
   
                if (lsd) then
                   if (ispin == 1) then
                      vxc_tmp_h(:,:,1) = vxc_tmp_h(:,:,1) + efac*orb_rho_set_h%rhoa(:,:,1)
                      vxc_tmp_s(:,:,1) = vxc_tmp_s(:,:,1) + efac*orb_rho_set_s%rhoa(:,:,1)
                   else
                      vxc_tmp_h(:,:,2) = vxc_tmp_h(:,:,2) + efac*orb_rho_set_h%rhob(:,:,1)
                      vxc_tmp_s(:,:,2) = vxc_tmp_s(:,:,2) + efac*orb_rho_set_s%rhob(:,:,1)
                   end if
                else                                         
                   vxc_tmp_h(:,:,1) = vxc_tmp_h(:,:,1) + efac*orb_rho_set_h%rho(:,:,1)
                   vxc_tmp_s(:,:,1) = vxc_tmp_s(:,:,1) + efac*orb_rho_set_s%rho(:,:,1)
                end if

             end do ! orb

          end do ! ispin

          if (lsd) then
             do ir=1, nr
                do ia=1, na
                   if (rho_set_h%rhoa(ia,ir,1) > rho_set_h%rho_cutoff) &
                        vxc_GLLB_h(ia,ir,1) = vxc_GLLB_h(ia,ir,1) + &
                                              weight(ia,ir)*vxc_tmp_h(ia,ir,1)/rho_set_h%rhoa(ia,ir,1)
                   if (rho_set_h%rhob(ia,ir,1) > rho_set_h%rho_cutoff) &
                        vxc_GLLB_h(ia,ir,2) = vxc_GLLB_h(ia,ir,2) + &
                                              weight(ia,ir)*vxc_tmp_h(ia,ir,2)/rho_set_h%rhob(ia,ir,1)
                   if (rho_set_s%rhoa(ia,ir,1) > rho_set_s%rho_cutoff) &
                        vxc_GLLB_s(ia,ir,1) = vxc_GLLB_s(ia,ir,1) + &
                                              weight(ia,ir)*vxc_tmp_s(ia,ir,1)/rho_set_s%rhoa(ia,ir,1)
                   if (rho_set_s%rhob(ia,ir,1) > rho_set_s%rho_cutoff) &
                        vxc_GLLB_s(ia,ir,2) = vxc_GLLB_s(ia,ir,2) + &
                                              weight(ia,ir)*vxc_tmp_s(ia,ir,2)/rho_set_s%rhob(ia,ir,1)
                end do
             end do
          else
             do ir=1, nr
                do ia=1, na
                   if (rho_set_h%rho(ia,ir,1) > rho_set_h%rho_cutoff) &
                        vxc_GLLB_h(ia,ir,1) = vxc_GLLB_h(ia,ir,1) + &
                                              weight(ia,ir)*vxc_tmp_h(ia,ir,1)/rho_set_h%rho(ia,ir,1)
                   if (rho_set_s%rho(ia,ir,1) > rho_set_s%rho_cutoff) &
                        vxc_GLLB_s(ia,ir,1) = vxc_GLLB_s(ia,ir,1) + &
                                              weight(ia,ir)*vxc_tmp_s(ia,ir,1)/rho_set_s%rho(ia,ir,1)
                end do
             end do
          end if

          vxc_SAOP_h = 0.0_dp ; vxc_SAOP_s = 0.0_dp

          do ispin=1, nspins

             do orb=1, homo(ispin)
                
                allocate(coeff_col(nrow(ispin),1), stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

                we_LB   = exp(-2.0_dp*(mo_eigenvalues(ispin)%array(homo(ispin)) - &
                                       mo_eigenvalues(ispin)%array(orb))**2)
                we_GLLB = 1.0_dp - we_LB
                if (.not.lsd) then
                   we_LB   = 2.0_dp * we_LB
                   we_GLLB = 2.0_dp * we_GLLB
                end if
                
                vxc_tmp_h(:,:,ispin) = we_LB*vxc_LB_h(:,:,ispin) + &
                                       we_GLLB*vxc_GLLB_h(:,:,ispin)
                vxc_tmp_s(:,:,ispin) = we_LB*vxc_LB_s(:,:,ispin) + &
                                       we_GLLB*vxc_GLLB_s(:,:,ispin)

                call cp_fm_set_all(single_mo_coeff(ispin)%matrix, 0.0_dp)
                call cp_fm_get_submatrix(mo_coeff(ispin)%matrix, coeff_col, &
                                         1, orb, nrow(ispin), 1, error=error)
                call cp_fm_set_submatrix(single_mo_coeff(ispin)%matrix, coeff_col, &
                                         1, orb, error=error)
                call set_matrix(orbital_density_matrix(ispin)%matrix, value=0.0_dp)
                call cp_sm_plus_fm_fm_t(orbital_density_matrix(ispin)%matrix,&
                                        matrix_v=single_mo_coeff(ispin)%matrix,&
                                        matrix_g=single_mo_coeff(ispin)%matrix,&
                                        ncol=ncol(ispin), &
                                        alpha=1.0_dp)

                deallocate(coeff_col, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

                ! This calculates the CPC and density on the grids for every atom even though
                ! we need it only for iatom at the moment. It seems that to circumvent this,
                ! the routines must be adapted to calculate just iatom
                call calculate_rho_atom_coeff(qs_env, orbital_density_matrix, &
                                              local_rho_set%rho_atom_set)
                call prepare_gapw_den(qs_env,local_rho_set,.FALSE.,error)
                   
                rho_atom => local_rho_set%rho_atom_set(iatom)
                nullify(r_h,r_s,dr_h,dr_s,r_h_d,r_s_d)
                call get_rho_atom(rho_atom=rho_atom, &
                                  rho_rad_h=r_h, rho_rad_s=r_s)
                do ir=1, nr
                   call calc_rho_angular(atomic_grid, harmonics, nspins, .false., &
                                         ir, r_h, r_s, rho_h, rho_s, &
                                         dr_h, dr_s, r_h_d, r_s_d, drho_h, drho_s, error)
                   call fill_rho_set(orb_rho_set_h,lsd,nspins,needs_orbs,rho_h,drho_h,na,ir,error=error)
                   call fill_rho_set(orb_rho_set_s,lsd,nspins,needs_orbs,rho_s,drho_s,na,ir,error=error)
                end do
   
                if (lsd) then
                   if (ispin == 1) then
                      vxc_SAOP_h(:,:,1) = vxc_SAOP_h(:,:,1) + vxc_tmp_h(:,:,1)*orb_rho_set_h%rhoa(:,:,1)
                      vxc_SAOP_s(:,:,1) = vxc_SAOP_s(:,:,1) + vxc_tmp_s(:,:,1)*orb_rho_set_s%rhoa(:,:,1)
                   else
                      vxc_SAOP_h(:,:,2) = vxc_SAOP_h(:,:,2) + vxc_tmp_h(:,:,2)*orb_rho_set_h%rhob(:,:,1)
                      vxc_SAOP_s(:,:,2) = vxc_SAOP_s(:,:,2) + vxc_tmp_s(:,:,2)*orb_rho_set_s%rhob(:,:,1)
                   end if
                else
                   vxc_SAOP_h(:,:,1) = vxc_SAOP_h(:,:,1) + vxc_tmp_h(:,:,1)*orb_rho_set_h%rho(:,:,1)
                   vxc_SAOP_s(:,:,1) = vxc_SAOP_s(:,:,1) + vxc_tmp_s(:,:,1)*orb_rho_set_s%rho(:,:,1)
                end if
   
             end do ! orb
                
          end do ! ispin

          if (lsd) then
             do ir=1, nr
                do ia=1, na
                   if (rho_set_h%rhoa(ia,ir,1) > rho_set_h%rho_cutoff) then
                      vxc_SAOP_h(ia,ir,1) = vxc_SAOP_h(ia,ir,1) / rho_set_h%rhoa(ia,ir,1)
                   else
                      vxc_SAOP_h(ia,ir,1) = 0.0_dp
                   end if
                   if (rho_set_h%rhob(ia,ir,1) > rho_set_h%rho_cutoff) then
                      vxc_SAOP_h(ia,ir,2) = vxc_SAOP_h(ia,ir,2) / rho_set_h%rhob(ia,ir,1)
                   else
                      vxc_SAOP_h(ia,ir,2) = 0.0_dp
                   end if
                   if (rho_set_s%rhoa(ia,ir,1) > rho_set_s%rho_cutoff) then
                      vxc_SAOP_s(ia,ir,1) = vxc_SAOP_s(ia,ir,1) / rho_set_s%rhoa(ia,ir,1)
                   else
                      vxc_SAOP_s(ia,ir,1) = 0.0_dp
                   end if
                   if (rho_set_s%rhob(ia,ir,1) > rho_set_s%rho_cutoff) then
                      vxc_SAOP_s(ia,ir,2) = vxc_SAOP_s(ia,ir,2) / rho_set_s%rhob(ia,ir,1)
                   else
                      vxc_SAOP_s(ia,ir,2) = 0.0_dp
                   end if
                end do
             end do
          else
             do ir=1, nr
                do ia=1, na
                   if (rho_set_h%rho(ia,ir,1) > rho_set_h%rho_cutoff) then
                      vxc_SAOP_h(ia,ir,1) = vxc_SAOP_h(ia,ir,1) / rho_set_h%rho(ia,ir,1)
                   else
                      vxc_SAOP_h(ia,ir,1) = 0.0_dp
                   end if
                   if (rho_set_s%rho(ia,ir,1) > rho_set_s%rho_cutoff) then
                      vxc_SAOP_s(ia,ir,1) = vxc_SAOP_s(ia,ir,1) / rho_set_s%rho(ia,ir,1)
                   else
                      vxc_SAOP_s(ia,ir,1) = 0.0_dp
                   end if
                end do
             end do
          end if

          rho_atom => rho_atom_set(iatom)
          atom_kind => atomic_kind_set(ikind)
          select case (oe_corr)
          case(oe_lb)
             call gaVxcgb_noGC(vxc_LB_h,vxc_LB_s,atom_kind,rho_atom,nspins)
          case (oe_gllb)             
             call gaVxcgb_noGC(vxc_GLLB_h,vxc_GLLB_s,atom_kind,rho_atom,nspins)
          case (oe_saop)          
             call gaVxcgb_noGC(vxc_SAOP_h,vxc_SAOP_s,atom_kind,rho_atom,nspins)
          case default
             CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
          end select

       end do

       deallocate(rho_h, rho_s, weight, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       deallocate(vxc_LB_h, vxc_LB_s, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       deallocate(vxc_GLLB_h, vxc_GLLB_s, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       deallocate(vxc_tmp_h, vxc_tmp_s, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       deallocate(vxc_SAOP_h, vxc_SAOP_s, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       deallocate(drho_h, drho_s, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
            
       call xc_dset_release(deriv_set, error=error)
       call xc_rho_set_release(rho_set_h, error=error)
       call xc_rho_set_release(rho_set_s, error=error)
       call xc_rho_set_release(orb_rho_set_h, error=error)
       call xc_rho_set_release(orb_rho_set_s, error=error)

    end do


    ! finally add the contributions to the KS-matrix
!    do ispin = 1,nspins
       call update_ks_atom(qs_env, qs_env%matrix_ks, &
                           rho_structure%rho_ao,&
                           .false., error=error)
!    end do
    
    !---------!
    ! Cleanup !
    !---------!
    xc_control = xc_control_sav
    call local_rho_set_release(local_rho_set,error)
    do ispin=1, nspins
       call cp_fm_release(single_mo_coeff(ispin)%matrix, error=error)
       call deallocate_matrix(orbital_density_matrix(ispin)%matrix)
    end do
    deallocate(mo_coeff, single_mo_coeff, mo_eigenvalues, &
               orbital_density_matrix, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE gapw_add_atomic_saop_pot

  !*****************************************************************

  SUBROUTINE add_lb_pot(pot, rho_set, lsd, spin, error)

    REAL(KIND=dp), DIMENSION(:,:,:), POINTER   :: pot
    TYPE(xc_rho_set_type), POINTER             :: rho_set
    LOGICAL, INTENT(IN)                        :: lsd
    INTEGER, INTENT(IN)                        :: spin
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                                 :: error

    REAL(KIND=dp), PARAMETER                   :: ob3 = 1.0_dp/3.0_dp

    INTEGER                                    :: i, j, k
    INTEGER, DIMENSION(2, 3)                   :: bo
    REAL(KIND=dp)                              :: n, n_13, x, x2

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

  !*****************************************************************

  SUBROUTINE calc_2excpbe(pot, rho_set, e_uniform, lsd)

    REAL(KIND=dp), DIMENSION(:,:,:), POINTER :: pot
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    REAL(KIND=dp), DIMENSION(:,:,:), POINTER :: e_uniform
    LOGICAL, INTENT(IN)                      :: lsd

    INTEGER                                  :: i, j, k
    INTEGER, DIMENSION(2,3)                  :: bo
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
