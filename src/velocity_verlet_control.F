!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
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

  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_p_type,&
                                             use_fist_force,&
                                             use_kg_force,&
                                             use_qs_force
  USE integrator,                      ONLY: nve,&
                                             nvt
  USE mathconstants,                   ONLY: zero
  USE md,                              ONLY: simulation_parameters_type
  USE md_environment_types,            ONLY: get_md_env,&
                                             md_environment_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE
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

SUBROUTINE velocity_verlet ( md_env )
  

    TYPE(md_environment_type), POINTER       :: md_env

    INTEGER                                  :: handle, in_use
    TYPE(force_env_p_type), POINTER          :: force_env_p ( : )
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar

!------------------------------------------------------------------------------

  CALL timeset ( 'velocity_verlet', 'I', ' ', handle )

! get force environment
  CALL get_md_env ( md_env, force_env_p=force_env_p, simpar=simpar )
! get in_use
  CALL force_env_get ( force_env=force_env_p ( 1 ) % force_env, in_use=in_use )

  SELECT CASE ( in_use )
  CASE ( use_fist_force )
     SELECT CASE ( TRIM ( simpar % ensemble ) )
     CASE DEFAULT
        CALL stop_program ( 'velocity_verlet','integrator not implemented')
     CASE ( 'NVE')
        CALL nve ( md_env )
     CASE ( 'NVT')
        CALL nvt ( md_env )
     CASE ( 'NPT_I')
       ! CALL npt_i ( md_env )
     CASE ( 'NPT_F')
       ! CALL npt_f ( md_env )
     END SELECT
 CASE default
    SELECT CASE (TRIM(simpar%ensemble))
    CASE ("NVE")
      CALL nve(md_env)
    CASE ("NVT")
      CALL nvt(md_env)
    CASE DEFAULT
      CALL stop_program("velocity_verlet",&
                        "integrator not implemented for a generic force_env")
    END SELECT
    
  END SELECT

  CALL timestop ( zero, handle )

END SUBROUTINE velocity_verlet

!******************************************************************************
END MODULE velocity_verlet_control
