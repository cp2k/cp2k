!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Provides an interface to the velocity-verlet based integrator
!>      routines for all ensembles
!> \author CJM (11-SEPT-2002)
! *****************************************************************************
MODULE velocity_verlet_control

  
  USE force_env_types,                 ONLY: force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: &
       isokin_ensemble, langevin_ensemble, npe_f_ensemble, npe_i_ensemble, &
       nph_uniaxial_damped_ensemble, nph_uniaxial_ensemble, npt_f_ensemble, &
       npt_i_ensemble, nve_ensemble, nvt_adiabatic_ensemble, nvt_ensemble, &
       reftraj_ensemble
  USE integrator,                      ONLY: &
       isokin, langevin, nph_uniaxial, nph_uniaxial_damped, npt_f, npt_i, &
       nve, nve_respa, nvt, nvt_adiabatic, reftraj
  USE md_environment_types,            ONLY: get_md_env,&
                                             md_environment_type
  USE simpar_types,                    ONLY: simpar_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'velocity_verlet_control'
  PUBLIC :: velocity_verlet

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param md_env ...
!> \param globenv ...
!> \par History
!>      none
!> \author CJM
! *****************************************************************************
  SUBROUTINE velocity_verlet ( md_env, globenv)

    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(global_environment_type), POINTER   :: globenv

    CHARACTER(LEN=*), PARAMETER :: routineN = 'velocity_verlet', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(simpar_type), POINTER               :: simpar

    CALL timeset (routineN, handle )

    ! Get force environment
    CALL get_md_env ( md_env, force_env=force_env, simpar=simpar)

    ! RESPA implemented only for NVE
    IF(simpar%do_respa .AND. nve_ensemble.NE.simpar % ensemble) THEN
       CALL stop_program(routineN,moduleN,__LINE__,&
                         "RESPA integrator not implemented for this ensemble")
    END IF

    ! Choice of the ensemble
    SELECT CASE (simpar%ensemble)
    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,&
                         "Integrator not implemented")
    CASE (nve_ensemble)
       IF(simpar%do_respa)THEN
          CALL nve_respa(md_env)
       ELSE
          CALL nve (md_env, globenv)
       END IF
    CASE (nvt_ensemble)
       CALL nvt (md_env, globenv)
    CASE (nvt_adiabatic_ensemble)
       CALL nvt_adiabatic (md_env, globenv)
    CASE (isokin_ensemble)
       CALL isokin (md_env)
    CASE (npt_i_ensemble)
       CALL npt_i (md_env, globenv)
    CASE (npt_f_ensemble)
       CALL npt_f (md_env, globenv)
    CASE (nph_uniaxial_ensemble)
       CALL nph_uniaxial (md_env)
    CASE (nph_uniaxial_damped_ensemble)
       CALL nph_uniaxial_damped (md_env)
    CASE (reftraj_ensemble)
       CALL reftraj (md_env)
    CASE (langevin_ensemble)
       CALL langevin(md_env)
    CASE (npe_f_ensemble)
       CALL npt_f (md_env, globenv)
    CASE (npe_i_ensemble)
       CALL npt_i (md_env, globenv)
    END SELECT

    CALL timestop(handle)

  END SUBROUTINE velocity_verlet

END MODULE velocity_verlet_control
