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
    
  USE coefficient_types, ONLY : PLANEWAVES, ATOMICORBITALS
  USE force_control, ONLY : get_force_env, FIST, QUICKSTEP, POL, KIM_GORDON
  USE kinds, ONLY : dbl
  USE integrator, ONLY : nve, nvt, npt_i, npt_f, nve_ao, nvt_ao, nve_pw, nvt_pw, &
                         nvt_pimd
  USE mathconstants, ONLY :  zero
  USE md_environment_types, ONLY : md_environment_type
  USE termination, ONLY : stop_program
  USE timings, ONLY : timeset, timestop

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
  
  IMPLICIT NONE

! Arguments
  TYPE ( md_environment_type ), INTENT ( INOUT ) :: md_env
  
! locals
  INTEGER :: handle, in_use

!------------------------------------------------------------------------------
  CALL timeset ( 'VELOCITY VERLET', 'I', ' ', handle )

! get force environment
  CALL get_force_env ( in_use )

  SELECT CASE ( in_use )
  CASE ( FIST )
    IF ( md_env % simpar % pimd_params % beads > 0 ) THEN
       SELECT CASE ( md_env % simpar % ensemble )
       CASE DEFAULT
          CALL stop_program ( 'velocity','PIMD not possible' ) 
       CASE ( 'NVT' )
          CALL nvt_pimd ( md_env )
       END SELECT
    ELSE
       SELECT CASE ( md_env % simpar % ensemble )
       CASE DEFAULT
          CALL stop_program ( 'velocity_verlet','integrator not implemented')
       CASE ( 'NVE')
          CALL nve ( md_env )
       CASE ( 'NVT')
          CALL nvt ( md_env )
       CASE ( 'NPT_I')
          CALL npt_i ( md_env )
       CASE ( 'NPT_F')
          CALL npt_f ( md_env )
       END SELECT
    END IF
  CASE ( POL )
    SELECT CASE ( md_env % simpar % ensemble_coef )
    CASE DEFAULT
       CALL stop_program ( 'velocity_verlet','integrator not implemented for POL')
    CASE ( 'NVE')
       SELECT CASE ( md_env % struc ( 1 ) % coef_pos ( 1 ) % in_use )
       CASE ( ATOMICORBITALS )
         CALL nve_ao ( md_env )
       CASE ( PLANEWAVES )
         CALL nve_pw ( md_env )
       END SELECT
    CASE ( 'NVT')
       SELECT CASE ( md_env % struc ( 1 ) % coef_pos ( 1 ) % in_use )
       CASE ( ATOMICORBITALS )
         CALL nvt_ao ( md_env )
       CASE ( PLANEWAVES )
         CALL nvt_pw ( md_env )
       END SELECT
    END SELECT
  CASE ( QUICKSTEP )
    SELECT CASE (md_env%simpar%ensemble)
    CASE ("NVE")
      CALL nve(md_env)
    CASE ("NVT")
      CALL nvt(md_env)
    CASE DEFAULT
      CALL stop_program("velocity_verlet",&
                        "integrator not implemented for QUICKSTEP")
    END SELECT
  CASE ( KIM_GORDON )
    IF (ASSOCIATED (  md_env % struc ( 1 ) % coef_pos  ) ) THEN 
      SELECT CASE ( md_env % simpar % ensemble )
      CASE DEFAULT
         CALL stop_program ( 'velocity_verlet','integrator not implemented for KG')
      CASE ( 'NVE')
         SELECT CASE ( md_env % struc ( 1 ) % coef_pos ( 1 ) % in_use )
         CASE ( ATOMICORBITALS )
           CALL nve_ao ( md_env )
         CASE ( PLANEWAVES )
           CALL nve_pw ( md_env )
         END SELECT
      CASE ( 'NVT')
         SELECT CASE ( md_env % struc ( 1 ) % coef_pos ( 1 ) % in_use )
         CASE ( ATOMICORBITALS )
           CALL nvt_ao ( md_env )
         CASE ( PLANEWAVES )
           CALL nvt_pw ( md_env )
         END SELECT
      END SELECT
    ELSE
      SELECT CASE ( md_env % simpar % ensemble )
      CASE DEFAULT
         CALL stop_program ( 'velocity_verlet','integrator not implemented for KG')
      CASE ( 'NVE')
         CALL nve ( md_env )
      CASE ( 'NVT')
         CALL nvt ( md_env )
      END SELECT
    END IF
  END SELECT

  CALL timestop ( zero, handle )

END SUBROUTINE velocity_verlet

!******************************************************************************
END MODULE velocity_verlet_control
