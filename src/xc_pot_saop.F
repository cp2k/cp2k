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
  USE pw_pool_types,                   ONLY: pw_pool_type
  USE pw_types,                        ONLY: pw_p_type,&
                                             pw_release
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_integrate_potential,          ONLY: integrate_v_rspace
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE simulation_cell,                 ONLY: cell_type
  USE sparse_matrix_types,             ONLY: real_matrix_p_type
  USE xc,                              ONLY: xc_vxc_pw_create

  IMPLICIT NONE

  PUBLIC :: add_saop_pot

  CHARACTER(LEN=*), PARAMETER :: moduleN="xc_pot_saop"

  ! should be eliminated
  REAL(KIND=dp), PARAMETER :: alpha=1.19_dp, beta=0.01_dp

CONTAINS

  !*****************************************************************

  !!
  !! Adds the saop potential to the ks_matrix
  !!
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

    LOGICAL                                  :: failure, ionode

    REAL(KIND=dp)                            :: xc_energy

    TYPE(cell_type), POINTER                 :: cell
    TYPE(coeff_type), DIMENSION(:), POINTER  :: vxc_SAOP
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsystem_p_type), DIMENSION(1)  :: subsys
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: vxc_GLLB, vxc_LB, vxc_LB_tmp, vxc_tau
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_g, rho_r, tau
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_rho_type), POINTER               :: rho_struct

    TYPE(xc_control_type), POINTER           :: xc_control
    TYPE(xc_control_type)                    :: my_xc_control
    TYPE(dft_control_type), POINTER          :: dft_control

    failure = .FALSE.

    NULLIFY(pw_env, dft_control, xc_control, auxbas_pw_pool)
    NULLIFY(rho_g, rho_r, tau, rho_struct)
    NULLIFY(vxc_GLLB, vxc_LB, vxc_LB_tmp, vxc_SAOP, vxc_tau)
    logger => cp_error_get_logger(error)
    ionode = (logger%para_env%mepos==logger%para_env%source)
    output_unit = cp_logger_get_default_unit_nr(logger)

    if (ionode) write (output_unit, *) "in "//routineP

    CALL get_qs_env(qs_env=qs_env,&
                    rho=rho_struct,&
                    pw_env=pw_env,&
                    subsys=subsys(1)%subsys,&
                    cell=cell,&
                    dft_control=dft_control)
    xc_control=>dft_control%xc_control
    CALL pw_env_get(pw_env,auxbas_pw_pool=auxbas_pw_pool,error=error)

    nspins = dft_control%nspins

    ALLOCATE(rho_r(nspins), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO ispin=1,nspins
       rho_r(ispin)%pw => rho_struct%rho_r(ispin)%pw
    END DO

    my_xc_control%gradient_functionals     = xc_control%gradient_functionals
    my_xc_control%crossterms               = xc_control%crossterms


    !---------------------------------------------!
    ! Construct 2 times PBE potential in vxc_GLLB !
    !---------------------------------------------!
    xc_control%functionals = ""
    CALL xc_get_functional ( 0, "PBE", xc_control%functionals, &
                                       xc_control%gradient_functionals, &
                                       xc_control%crossterms )
    WHERE (xc_control%functionals == "") xc_control%functionals = "NONE"

    CALL xc_vxc_pw_create(vxc_GLLB, vxc_tau, xc_energy, rho_r, rho_g, tau,&
                          dft_control, cell, auxbas_pw_pool, subsys, &
                          .FALSE., error)
    n_gridpoints = SIZE(vxc_GLLB(1)%pw%cr3d,1) * &
                   SIZE(vxc_GLLB(1)%pw%cr3d,2) * &
                   SIZE(vxc_GLLB(1)%pw%cr3d,3)
    DO ispin=1, nspins
       CALL dscal(n_gridpoints, 2.0_dp, vxc_GLLB(ispin)%pw%cr3d, 1)
    END DO


    !----------------------------------------!
    ! Construct the LB94 potential in vxc_LB !
    !----------------------------------------!
    xc_control%functionals = ""
    CALL xc_get_functional ( 1, "SLATER", xc_control%functionals, &
                                          xc_control%gradient_functionals, &
                                          xc_control%crossterms )
    WHERE (xc_control%functionals == "") xc_control%functionals = "NONE"
    CALL xc_vxc_pw_create(vxc_LB, vxc_tau, xc_energy, rho_r, rho_g, tau,&
                          dft_control, cell, auxbas_pw_pool, subsys, &
                          .FALSE., error)
    DO ispin=1, nspins
       CALL dscal(n_gridpoints, alpha, vxc_LB(ispin)%pw%cr3d, 1)
    END DO
    xc_control%functionals = ""
    CALL xc_get_functional ( 2, "VWN", xc_control%functionals, &
                                       xc_control%gradient_functionals, &
                                       xc_control%crossterms )
    WHERE (xc_control%functionals == "") xc_control%functionals = "NONE"
    CALL xc_vxc_pw_create(vxc_LB_tmp, vxc_tau, xc_energy, rho_r, rho_g, tau,&
                          dft_control, cell, auxbas_pw_pool, subsys, &
                          .FALSE., error)
    DO ispin=1, nspins
       CALL daxpy(n_gridpoints, 1.0_dp, vxc_LB_tmp(ispin)%pw%cr3d, 1, &
                                        vxc_LB(ispin)%pw%cr3d, 1)
    END DO

    !-------------------------!
    ! restore original values !
    !-------------------------!
    xc_control%functionals = my_xc_control%functionals
    xc_control%gradient_functionals = my_xc_control%gradient_functionals
    xc_control%crossterms = my_xc_control%crossterms

    !---------------!
    ! Assemble SAOP !
    !---------------!
    ALLOCATE(vxc_SAOP(nspins), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO ispin=1, nspins
       ! just add them for testing purposes
       CALL daxpy(n_gridpoints, 1.0_dp, vxc_LB(ispin)%pw%cr3d, 1, &
                                        vxc_GLLB(ispin)%pw%cr3d, 1)
       CALL coeff_init(vxc_SAOP(ispin),pw=vxc_GLLB(ispin)%pw,&
            error=error)
    END DO
    
    DO ispin=1, nspins
       CALL pw_release(vxc_LB(ispin)%pw)
       CALL pw_release(vxc_LB_tmp(ispin)%pw)
    END DO
    
    DEALLOCATE(vxc_GLLB, vxc_LB, vxc_LB_tmp, stat=stat)
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

       CALL coeff_deallocate(vxc_SAOP(ispin))
    END DO
    
    DEALLOCATE(rho_r, vxc_SAOP, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE add_saop_pot

  !*****************************************************************

  SUBROUTINE calc_lb_pot(pot, qs_env, glob_env, error)

    INTEGER, INTENT(OUT)                     :: pot
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), &
         INTENT(IN), TARGET                  :: glob_env
    TYPE(cp_error_type), &
         INTENT(INOUT), OPTIONAL             :: error

    INTEGER                                  :: i
    TYPE(xc_control_type), POINTER           :: xc_control
    TYPE(dft_control_type), POINTER          :: dft_control

    TYPE(xc_control_type)                    :: my_xc_control

    
    NULLIFY(dft_control, xc_control)

    CALL get_qs_env(qs_env, dft_control=dft_control)
    xc_control => dft_control%xc_control

    !--------------------------------!
    ! save original values for later !
    !--------------------------------!
    my_xc_control%gradient_functionals     = xc_control%gradient_functionals
    my_xc_control%crossterms               = xc_control%crossterms
    xc_control%functionals = ""
    CALL xc_get_functional ( 0, "PBE", xc_control%functionals, &
                                       xc_control%gradient_functionals, &
                                       xc_control%crossterms )
    DO i=1, 3
       IF (xc_control%functionals(i) == "") THEN 
          xc_control%functionals(i) = "NONE"
       END IF
    END DO
    
    ! CALL PBEX AND MULTIPLY BY ALPHA

    ! CALL PBEC AND ADD

    !-------------------------!
    ! restore original values !
    !-------------------------!
    xc_control%functionals = my_xc_control%functionals
    xc_control%gradient_functionals = my_xc_control%gradient_functionals
    xc_control%crossterms = my_xc_control%crossterms
    
    ! ADD REST

  END SUBROUTINE calc_lb_pot

  !*****************************************************************

END MODULE xc_pot_saop
