!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/force_control [1.0] *
!!
!!   NAME
!!     force_control
!!
!!   FUNCTION
!!     Interface routine for the force calculations
!!
!!   AUTHOR
!!     CJM & JGH
!!
!!   MODIFICATION HISTORY
!!     cjm, FEB-20-2001: pass variable box_ref
!!     cjm, SEPT-12-2002: major reorganization
!!
!!   SOURCE
!******************************************************************************

MODULE force_control
  USE fist_environment_types,          ONLY: fist_environment_type
  USE fist_force,                      ONLY: fist_force_control => force_control
  USE global_types,                    ONLY: global_environment_type
  USE kg_environment_types,            ONLY: kg_environment_type
  USE kg_force,                        ONLY: kg_force_control
  USE kinds,                           ONLY: dbl
  USE qs_environment_types,            ONLY: qs_environment_type
  USE qs_force,                        ONLY: qs_forces_update
  USE termination,                     ONLY: stop_program
  IMPLICIT NONE

  TYPE ( fist_environment_type ), POINTER :: fist_env 
  TYPE ( kg_environment_type ), POINTER :: kg_env
  TYPE ( qs_environment_type ), POINTER :: qs_env

  INTEGER, PARAMETER :: FIST=501, QUICKSTEP=502, POL=503, KIM_GORDON=504 
  INTEGER :: in_use

  PRIVATE
  PUBLIC :: force, set_force_env, get_force_env, FIST, QUICKSTEP, POL, KIM_GORDON
  
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** force_control/force [1.0] *
!!
!!   NAME
!!     force
!!
!!   SYNOPSIS
!!     Subroutine force(globenv, box_change)
!!       Implicit None
!!       Type(global_environment_type), Intent (IN):: globenv
!!       Logical, Intent (IN):: box_change
!!     End Subroutine force
!!
!!   FUNCTION
!!     Interface routine for force calculations
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE force ( globenv, box_change )
  
  IMPLICIT NONE

! Arguments
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv
  LOGICAL, INTENT ( IN ) :: box_change

! local
  REAL ( dbl ) :: e_pot
  INTEGER :: beads, i
  
!------------------------------------------------------------------------------

  SELECT CASE ( in_use )
    CASE ( FIST )
      SELECT CASE ( fist_env % pimd )
      CASE ( .TRUE. )
        beads = SIZE ( fist_env % rep_env )
        e_pot = 0._dbl
        DO i = 1, beads
          CALL fist_force_control ( fist_env % rep_env ( i ), &
                                    fist_env % ewald_param, &
                                    fist_env % potparm, fist_env % thermo, &
                                    box_change, globenv )
          e_pot = e_pot + fist_env % thermo % pot
        ENDDO
        fist_env % thermo % pot = e_pot
      CASE ( .FALSE. )
          CALL fist_force_control ( fist_env % rep_env ( 1 ), &
                                    fist_env % ewald_param, &
                                    fist_env % potparm, fist_env % thermo, &
                                    box_change, globenv )
      END SELECT 
    CASE ( QUICKSTEP )
      CALL qs_forces_update(qs_env,globenv)
    CASE ( POL )
    CASE ( KIM_GORDON )
      CALL kg_force_control ( kg_env, globenv )
  END SELECT

END SUBROUTINE force

!******************************************************************************
!!****** force_control/set_force_env [1.0] *
!!
!!   NAME
!!     set_force_env
!!
!!   SYNOPSIS
!!     Subroutine set_force_env(fist_env_type, kg_env_type, qs_env_type)
!!       Implicit None
!!       Type(fist_environment_type), Target, Optional:: fist_env_type
!!       Type(kg_environment_type), Target, Optional:: kg_env_type
!!       Type(qs_environment_type), Target, Optional:: qs_env_type
!!     End Subroutine set_force_env
!!
!!   FUNCTION
!!     sets the globals for the force routine
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE set_force_env ( fist_env_type, kg_env_type, qs_env_type )
  IMPLICIT NONE
  TYPE ( fist_environment_type ), TARGET, OPTIONAL :: fist_env_type
  TYPE ( kg_environment_type ), TARGET, OPTIONAL :: kg_env_type
  TYPE ( qs_environment_type ), TARGET, OPTIONAL :: qs_env_type

! Classical MD
  IF ( PRESENT ( fist_env_type ) ) THEN
     fist_env => fist_env_type
     in_use = FIST
  ENDIF

! KIm-Grdon MD
  IF ( PRESENT ( kg_env_type ) ) THEN
     kg_env => kg_env_type
     in_use = KIM_GORDON
  ENDIF

! Quickstep
  IF ( PRESENT ( qs_env_type ) ) THEN
     qs_env => qs_env_type
     in_use = QUICKSTEP
  ENDIF

END SUBROUTINE set_force_env 
!******************************************************************************
SUBROUTINE get_force_env ( inuse )
  IMPLICIT NONE

  INTEGER, INTENT ( OUT ) :: inuse
  inuse = in_use

END SUBROUTINE get_force_env 
!******************************************************************************
END MODULE force_control
!******************************************************************************
