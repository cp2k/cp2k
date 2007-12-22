!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/velocity_verlet_control.F [1.0] *
!!
!!   NAME
!!     velocity_verlet_control.F
!!
!!   FUNCTION
!!     Provides an interface to the velocity-verlet based integrator
!!     routines for all ensembles
!!
!!   AUTHOR
!!     CJM (11-SEPT-2002)
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!******************************************************************************

MODULE velocity_verlet_control

  USE f77_blas
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type,&
                                             use_fist_force
  USE input_constants,                 ONLY: &
       isokin_ensemble, langevin_ensemble, nph_uniaxial_damped_ensemble, &
       nph_uniaxial_ensemble, npt_f_ensemble, npt_i_ensemble, nve_ensemble, &
       nvt_ensemble, reftraj_ensemble
  USE integrator,                      ONLY: &
       isokin, langevin, nph_uniaxial, nph_uniaxial_damped, npt_f, npt_i, &
       nve, nve_respa, nvt, reftraj
  USE mathconstants,                   ONLY: zero
  USE md,                              ONLY: simulation_parameters_type
  USE md_environment_types,            ONLY: get_md_env,&
                                             md_environment_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'velocity_verlet_control'
  PUBLIC :: velocity_verlet

!!***
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** velocity_verlet_control/velocity_verlet [1.0] *
!!
!!   NAME
!!     velocity_verlet
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** *************************************************************************

SUBROUTINE velocity_verlet ( md_env, error )


    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'velocity_verlet', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, in_use
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar

!------------------------------------------------------------------------------

  CALL timeset ( 'velocity_verlet', 'I', ' ', handle )

! get force environment
  CALL get_md_env ( md_env, force_env=force_env, simpar=simpar, error=error)
! get in_use
  CALL force_env_get ( force_env, in_use=in_use ,error=error)

  IF(simpar%do_respa .AND. nve_ensemble.NE.simpar % ensemble)THEN
        CALL stop_program ( 'velocity_verlet','respa integrator not implemented for this ensemble')
  END IF
  SELECT CASE ( in_use )
  CASE ( use_fist_force )
     SELECT CASE ( simpar % ensemble  )
     CASE DEFAULT
        CALL stop_program ( 'velocity_verlet','integrator not implemented')
     CASE ( nve_ensemble )
        IF(simpar%do_respa)THEN
           CALL nve_respa(md_env,error)
        ELSE
           CALL nve ( md_env, error )
        END IF
     CASE ( nvt_ensemble )
        CALL nvt ( md_env, error  )
     CASE ( isokin_ensemble )
        CALL isokin ( md_env, error  )
     CASE ( npt_i_ensemble )
        CALL npt_i ( md_env, error  )
     CASE ( npt_f_ensemble)
        CALL npt_f ( md_env, error  )
     CASE ( nph_uniaxial_ensemble)
        CALL nph_uniaxial ( md_env, error  )
     CASE ( nph_uniaxial_damped_ensemble )
        CALL nph_uniaxial_damped ( md_env, error  )
     CASE ( reftraj_ensemble)
        CALL reftraj ( md_env, error  )
     CASE (langevin_ensemble)
        CALL langevin(md_env, error )
     END SELECT
  CASE default
    SELECT CASE (simpar%ensemble)
     CASE ( nve_ensemble )
        IF(simpar%do_respa)THEN
           CALL nve_respa(md_env,error)
        ELSE
           CALL nve ( md_env, error )
        END IF
    CASE (nvt_ensemble)
      CALL nvt(md_env, error )
    CASE (isokin_ensemble)
      CALL isokin ( md_env, error  )
    CASE (reftraj_ensemble)
      CALL reftraj ( md_env, error  )
    CASE (langevin_ensemble)
      CALL langevin(md_env, error)
    CASE ( nph_uniaxial_ensemble)
       CALL nph_uniaxial ( md_env,error )
    CASE ( nph_uniaxial_damped_ensemble )
       CALL nph_uniaxial_damped ( md_env,error )
    CASE ( npt_i_ensemble )
       CALL npt_i ( md_env, error  )
    CASE ( npt_f_ensemble )
       CALL npt_f ( md_env, error  )
    CASE DEFAULT
      CALL stop_program("velocity_verlet",&
                        "integrator not implemented for a generic force_env")
    END SELECT

  END SELECT

  CALL timestop ( zero, handle )

END SUBROUTINE velocity_verlet

!******************************************************************************
END MODULE velocity_verlet_control
