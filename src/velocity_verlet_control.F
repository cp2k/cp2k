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

  USE coefficient_types,               ONLY: ATOMICORBITALS,&
                                             PLANEWAVES
  USE force_env_types,                   ONLY: force_env_get,&
                                             use_fist_force,&
                                             use_kg_force,&
                                             use_pol_force,&
                                             use_qs_force, &
                                             force_env_type
  USE integrator,                      ONLY: nve,&
                                             nvt
  USE kinds,                           ONLY: dbl
  USE mathconstants,                   ONLY: zero
  USE md_environment_types,            ONLY: md_environment_type,get_md_env
  USE md,                              ONLY: simulation_parameters_type
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
!!   SYNOPSIS
!!     Subroutine velocity_verlet(md_env)
!!       Implicit None
!!       Type(md_environment_type), Intent (INOUT):: md_env
!!     End Subroutine velocity_verlet
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
  
  IMPLICIT NONE

! Arguments
  TYPE ( md_environment_type ), POINTER :: md_env
  
! locals
  INTEGER :: handle, in_use
  TYPE ( simulation_parameters_type ), POINTER :: simpar
  TYPE ( force_env_type ), POINTER :: force_env

!------------------------------------------------------------------------------
  CALL timeset ( 'velocity_verlet', 'I', ' ', handle )

! get force environment
  CALL get_md_env ( md_env, force_env=force_env, simpar=simpar )
! get in_use
  CALL force_env_get ( force_env, in_use=in_use )

  SELECT CASE ( in_use )
  CASE ( use_fist_force )
    IF ( simpar % pimd_params % beads > 0 ) THEN
       SELECT CASE ( simpar % ensemble )
       CASE DEFAULT
          CALL stop_program ( 'velocity','PIMD not possible' ) 
       CASE ( 'NVT' )
    !      CALL nvt_pimd ( md_env )
       END SELECT
    ELSE
       SELECT CASE ( simpar % ensemble )
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
    END IF
  CASE ( use_qs_force )
    SELECT CASE (TRIM(simpar%ensemble))
    CASE ("NVE")
      CALL nve(md_env)
    CASE ("NVT")
      CALL nvt(md_env)
    CASE DEFAULT
      CALL stop_program("velocity_verlet",&
                        "integrator not implemented for QUICKSTEP")
    END SELECT
  CASE ( use_kg_force )
    SELECT CASE ( simpar % ensemble )
    CASE DEFAULT
       CALL stop_program ( 'velocity_verlet','integrator not implemented for KG')
    CASE ( 'NVE')
       CALL nve ( md_env )
    CASE ( 'NVT')
       CALL nvt ( md_env )
    END SELECT
  END SELECT

  CALL timestop ( zero, handle )

END SUBROUTINE velocity_verlet

!******************************************************************************
END MODULE velocity_verlet_control
