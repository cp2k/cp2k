!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/scp_debug [1.0] *
!!
!!   NAME
!!     scp_debug
!!
!!   FUNCTION
!!     scp_debug
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!   SOURCE
!******************************************************************************

MODULE scp_debug
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cell_types,                      ONLY: cell_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_walltime
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: REALDATA3D,&
                                             REALSPACE,&
                                             pw_integral_ab,&
                                             pw_p_type,&
                                             pw_release,&
                                             pw_scale,&
                                             pw_zero
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_ks_scp_methods,               ONLY: qs_ks_scp_did_change,&
                                             qs_ks_scp_update
  USE qs_ks_scp_types,                 ONLY: qs_ks_scp_env_type
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE scp_coeff_types,                 ONLY: aux_coeff_set_type,&
                                             aux_coeff_type,&
                                             get_aux_coeff
  USE scp_density_methods,             ONLY: integrate_rhoscp_vrspace,&
                                             update_rhoscp
  USE scp_energy_types,                ONLY: scp_energy_type
  USE scp_environment_types,           ONLY: get_scp_env,&
                                             scp_environment_type
  USE scp_hartree_1center,             ONLY: integrate_a_vhscp_b,&
                                             integrate_vhscp_gcore,&
                                             integrate_vhscp_gscp
  USE scp_rho_types,                   ONLY: scp_rho_type
  USE sparse_matrix_types,             ONLY: real_matrix_p_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'scp_debug'
  PUBLIC :: scp_qs_debug

!******************************************************************************

CONTAINS

!!****f** scp_debug/scp_qs_debug [1.1] *
!!
!!   NAME
!!    scp_qs_debug
!!
!!   FUNCTION
!!     Debugs the forces on the SCP coefficients using the qs_environment
!!   AUTHOR
!!    CJM
!!
!!   MODIFICATION HISTORY
!!
!!*** *************************************************************************
SUBROUTINE scp_qs_debug ( qs_env,  error )

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE scp_qs_debug (MODULE scp_debug)", &
      routineN = 'scp_qs_debug', routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, icoef_atom, ikind, &
                                                n_els, ncoef_atom, nkind, stat
    LOGICAL                                  :: failure, ionode
    REAL(dp)                                 :: em, ep, f
    REAL(dp), DIMENSION(:, :), POINTER       :: coeff
    REAL(KIND=dp)                            :: t1
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(aux_coeff_set_type), POINTER        :: aux_coeff_set
    TYPE(aux_coeff_type), POINTER            :: local_coeffs
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type), POINTER                 :: vks
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_ks_scp_env_type), POINTER        :: ks_scp_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks
    TYPE(scp_energy_type), POINTER           :: scp_energy
    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(scp_rho_type), POINTER              :: rho_scp

  CALL timeset ( 'scp_qs_debug','I',' ',handle )

  failure = .FALSE.
  NULLIFY(atomic_kind_set,aux_coeff_set,scp_env, matrix_ks, rho)
  NULLIFY(para_env,logger, coeff, local_coeffs, local_particles)
  NULLIFY(vks, rho_scp)
  para_env=>qs_env%para_env
  logger => cp_error_get_logger(error)
  ionode=para_env%ionode
  stat = 0
  t1 = m_walltime()

  CALL get_qs_env ( qs_env = qs_env, scp_env = scp_env, &
                    atomic_kind_set=atomic_kind_set,    &
                    local_particles=local_particles, &
                    particle_set = particle_set,  &
                    matrix_ks = matrix_ks, &
                    ks_env = ks_env, &
                    rho=rho, &
                    dft_control = dft_control, &
                    cell = cell, &
                    pw_env = pw_env, &
                    error = error )
  CALL pw_env_get ( pw_env, auxbas_pw_pool = auxbas_pw_pool, error = error )

! Getting the coeffs and the density structures
  CALL get_scp_env ( scp_env=scp_env,&
                     aux_coeff_set=aux_coeff_set,&
                     ks_scp_env = ks_scp_env, &
                     energy=scp_energy, error=error)
  WRITE ( *, * ) '**********************************IN DEBUG**********************************'
!-------------
  CALL integrate_rhoscp_vrspace ( scp_env, qs_env%ks_env%v_hartree_rspace, &
                                  pw_env, atomic_kind_set, &
                                  particle_set, local_particles, cell, &
                                  dft_control%qs_control%eps_gvg_rspace, &
                                  calculate_forces = .FALSE., just_energy = .FALSE., &
                                  debug=.TRUE., error=error )

! Subtracting off the V_HSCP from the potential
  ALLOCATE ( vks, stat = stat )
  CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
  CALL pw_pool_create_pw(auxbas_pw_pool,vks%pw,&
       use_data=REALDATA3D,in_space=REALSPACE,error=error)

  CALL pw_zero(vks%pw) 
  vks%pw%cr3d(:,:,:) = ks_env%v_hartree_rspace%pw%cr3d(:,:,:) - &
                       ks_scp_env%v_scp_rspace%pw%cr3d(:,:,:)* &
                       ks_env%v_hartree_rspace%pw%pw_grid%dvol 
  CALL pw_scale(vks%pw, 1._dp/ks_env%v_hartree_rspace%pw%pw_grid%dvol)


  nkind = SIZE ( aux_coeff_set % coeffs_of_kind )
  DO ikind = 1, nkind
    local_coeffs => aux_coeff_set % coeffs_of_kind ( ikind ) % coeffs
    IF ( .NOT. ASSOCIATED ( local_coeffs ) ) CYCLE
    CALL get_aux_coeff ( coeffs = local_coeffs, c = coeff, n_els = n_els, &
                         ncoef_atom = ncoef_atom, error = error )
    DO i = 1, n_els
      DO icoef_atom = 1, ncoef_atom
        coeff ( i, icoef_atom ) = coeff ( i, icoef_atom ) + 1.e-5_dp
        CALL update_rhoscp ( qs_env, error )
        CALL qs_ks_scp_did_change ( qs_env, full_reset = .TRUE., error = error )
        CALL qs_ks_scp_update ( qs_env, just_energy = .FALSE., error = error )
        CALL get_scp_env ( scp_env = scp_env, rho_scp = rho_scp,  &
                            energy=scp_energy, error = error )
        ep = 0._dp
        ep = ep + scp_energy%e_scp_hartree
        ep = ep + pw_integral_ab ( rho_scp % rhop_r % pw, vks % pw )

        coeff ( i, icoef_atom ) = coeff ( i, icoef_atom ) - 2.0_dp * 1.e-5_dp

        CALL update_rhoscp ( qs_env, error )
        CALL qs_ks_scp_did_change ( qs_env, full_reset = .TRUE., error = error )
        CALL qs_ks_scp_update ( qs_env, just_energy = .FALSE., error = error )
        CALL get_scp_env ( scp_env = scp_env, rho_scp = rho_scp,  &
                           energy=scp_energy, error = error )
        em = 0._dp
        em = em + scp_energy%e_scp_hartree
        em = em + pw_integral_ab ( rho_scp % rhop_r % pw, vks % pw )
        f = ( em - ep )/2._dp/1.e-5_dp
        coeff ( i, icoef_atom ) = coeff ( i, icoef_atom ) + 1.e-5_dp
        WRITE ( *, * ) 'FDEBUG_E[rho_ks + rho_scp]', ikind, i, icoef_atom, f
      END DO ! icoef_atom
    END DO ! i
  END DO ! ikind
  CALL pw_release (vks%pw,error=error)
  DEALLOCATE ( vks, stat = stat )
  CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

  WRITE ( *, * ) '****************************************************************************'
!-------------
  CALL integrate_a_vhscp_b ( qs_env, matrix_ks, rho%rho_ao, &
                             just_energy = .FALSE., error=error, debug = .TRUE.)

  CALL get_scp_env ( scp_env = scp_env, aux_coeff_set = aux_coeff_set, error = error )
  nkind = SIZE ( aux_coeff_set % coeffs_of_kind )
  DO ikind = 1, nkind
    local_coeffs => aux_coeff_set % coeffs_of_kind ( ikind ) % coeffs
    IF ( .NOT. ASSOCIATED ( local_coeffs ) ) CYCLE
    CALL get_aux_coeff ( coeffs = local_coeffs, c = coeff, n_els = n_els, &
                         ncoef_atom = ncoef_atom, error = error )
    DO i = 1, n_els
      DO icoef_atom = 1, ncoef_atom
        coeff ( i, icoef_atom ) = coeff ( i, icoef_atom ) + 1.e-5_dp
        CALL integrate_a_vhscp_b ( qs_env, matrix_ks, rho%rho_ao, &
                                   just_energy = .TRUE., error=error, debug = .FALSE. )
        CALL get_scp_env ( scp_env = scp_env, energy = scp_energy , error = error )
        ep = scp_energy % e_scp_ks_self
        coeff ( i, icoef_atom ) = coeff ( i, icoef_atom ) - 2.0_dp * 1.e-5_dp

        CALL integrate_a_vhscp_b ( qs_env, matrix_ks, rho%rho_ao, &
                                  just_energy =  .TRUE., error=error, debug = .FALSE.)
        CALL get_scp_env ( scp_env = scp_env, energy = scp_energy , error = error )
        em = scp_energy % e_scp_ks_self
        f = ( em - ep )/2._dp/1.e-5_dp
        coeff ( i, icoef_atom ) = coeff ( i, icoef_atom ) + 1.e-5_dp
        WRITE ( *, * ) 'FDEBUG_INTEGRATE_A_VHSCP_B', ikind, i, icoef_atom, f
      END DO ! icoef_atom
    END DO ! i
  END DO ! ikind
!-------------
  WRITE ( *, * ) '****************************************************************************'

  CALL integrate_vhscp_gscp ( scp_env, atomic_kind_set, local_particles, &
                              just_energy = .FALSE., error=error, debug=.TRUE. )
  nkind = SIZE ( atomic_kind_set )
  DO ikind = 1, nkind
    local_coeffs => aux_coeff_set % coeffs_of_kind ( ikind ) % coeffs
    IF ( .NOT. ASSOCIATED ( local_coeffs ) ) CYCLE
    CALL get_aux_coeff ( coeffs = local_coeffs, c = coeff, n_els = n_els, &
                         ncoef_atom = ncoef_atom, error = error )
    DO i = 1, n_els
      DO icoef_atom = 1, ncoef_atom
        coeff ( i, icoef_atom ) = coeff ( i, icoef_atom ) + 1.e-5_dp

        CALL integrate_vhscp_gscp ( scp_env, atomic_kind_set, local_particles, &
                                    just_energy = .TRUE., error=error )
        CALL get_scp_env ( scp_env = scp_env, energy = scp_energy , error = error )
        ep = scp_energy % e_scp_self
        coeff ( i, icoef_atom ) = coeff ( i, icoef_atom ) - 2.0_dp * 1.e-5_dp


        CALL integrate_vhscp_gscp ( scp_env, atomic_kind_set, local_particles, &
                                    just_energy = .TRUE., error = error )
        CALL get_scp_env ( scp_env = scp_env, energy = scp_energy , error = error )
        em = scp_energy % e_scp_self
        f = ( em - ep )/2._dp/1.e-5_dp
        WRITE ( *, * ) 'FDEBUG_INTEGRATE_VHSCP_GSCP', ikind, i, icoef_atom, f
        coeff ( i, icoef_atom ) = coeff ( i, icoef_atom ) + 1.e-5_dp
      END DO ! icoef_atom
    END DO ! i
  END DO ! ikind

!-------------
  WRITE ( *, * ) '****************************************************************************'
  CALL integrate_vhscp_gcore ( scp_env, atomic_kind_set, local_particles, &
                               just_energy = .FALSE., error = error, debug = .TRUE. )
  nkind = SIZE ( atomic_kind_set )
  DO ikind = 1, nkind
    local_coeffs => aux_coeff_set % coeffs_of_kind ( ikind ) % coeffs
    IF ( .NOT. ASSOCIATED ( local_coeffs ) ) CYCLE
    CALL get_aux_coeff ( coeffs = local_coeffs, c = coeff, n_els = n_els, &
                         ncoef_atom = ncoef_atom, error = error )
    DO i = 1, n_els
      DO icoef_atom = 1, ncoef_atom
        coeff ( i, icoef_atom ) = coeff ( i, icoef_atom ) + 1.e-5_dp

        CALL integrate_vhscp_gcore ( scp_env, atomic_kind_set, local_particles, &
                                    just_energy = .TRUE., error = error )
        CALL get_scp_env ( scp_env = scp_env, energy = scp_energy , error = error )
        ep = scp_energy % e_scp_core
        coeff ( i, icoef_atom ) = coeff ( i, icoef_atom ) - 2.0_dp * 1.e-5_dp


        CALL integrate_vhscp_gcore ( scp_env, atomic_kind_set, local_particles, &
                                    just_energy = .TRUE., error = error )
        CALL get_scp_env ( scp_env = scp_env, energy = scp_energy , error = error )
        em = scp_energy % e_scp_core
        f = ( em - ep )/2._dp/1.e-5_dp
        WRITE ( *, * ) 'FDEBUG_INTEGRATE_VHSCP_GCORE', ikind, i, icoef_atom, f
        coeff ( i, icoef_atom ) = coeff ( i, icoef_atom ) + 1.e-5_dp
      END DO ! icoef_atom
    END DO ! i
  END DO ! ikind

  WRITE ( *, * ) '**********************************END DEBUG*********************************'
  CALL timestop(0.0_dp,handle)

END SUBROUTINE scp_qs_debug

!******************************************************************************
END MODULE scp_debug
