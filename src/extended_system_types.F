!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!***** cp2k/extended_system_types [1.0] *
!!
!!   NAME
!!     extended_system_types
!!
!!   FUNCTION
!!   Lumps all possible extended system variables into one
!!   type for easy access and passing
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     Teodoro Laino - 09.2007 - University of Zurich
!!       Cleaned the typo.. no need to have an extended
!!       type. Thermostat and Barostat type have been created
!!
!!   SOURCE
!******************************************************************************

MODULE extended_system_types
  USE kinds,                           ONLY: dp
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: lnhc_dealloc,&
       lnhc_init,&
       lnhc_parameters_type,&
       npt_info_type,&
       nhc_info_type, &
       particle_kinetic_type,&
       particle_scale_type,&
       extended_energy_type,&
       nhc_set_type


  TYPE extended_energy_type
     REAL(KIND=dp)   ::  energy_g3x3_restraint,&
                         energy_fixd_restraint,&
                         energy_g4x6_restraint,&
                         energy_colv_restraint
  END TYPE extended_energy_type

  TYPE npt_info_type
     REAL (KIND=dp) :: eps
     REAL (KIND=dp) :: v
     REAL (KIND=dp) :: f
     REAL (KIND=dp) :: mass
  END TYPE npt_info_type

  TYPE nhc_info_type
     REAL (KIND=dp) :: eta
     REAL (KIND=dp) :: v
     REAL (KIND=dp) :: f
     REAL (KIND=dp) :: nkt
     REAL (KIND=dp) :: mass
  END TYPE nhc_info_type

  TYPE particle_kinetic_type
     REAL (KIND=dp), POINTER :: point
  END TYPE particle_kinetic_type

  TYPE particle_scale_type
     REAL (KIND=dp), POINTER :: point
  END TYPE particle_scale_type

  TYPE lnhc_parameters_type
     CHARACTER ( LEN = 10 )                  :: dis_type
     INTEGER                                 :: nyosh, nc, nhc_len, num_nhc
     REAL (KIND=dp), POINTER                 :: dt_yosh(:)
     TYPE ( nhc_info_type ), POINTER         :: nvt (:,:)
     TYPE ( particle_kinetic_type ), POINTER :: p_kin (:,:)
     TYPE ( particle_scale_type ), POINTER   :: p_scale (:,:)
     REAL (KIND=dp), POINTER                 :: v_scale (:)
     REAL (KIND=dp), POINTER                 :: s_kin (:)
     INTEGER, POINTER                        :: INDEX(:)
  END TYPE lnhc_parameters_type

  TYPE nhc_set_type
     INTEGER :: id_nr, ref_count, in_use
     TYPE ( lnhc_parameters_type ), POINTER :: nhc_part  ( : )
     TYPE ( lnhc_parameters_type ), POINTER :: nhc_baro ( : )
     TYPE ( lnhc_parameters_type ), POINTER :: nhc_coef (:)
     TYPE ( lnhc_parameters_type ), POINTER :: nhc_shell (:)
  END TYPE nhc_set_type

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'extended_system_types'

CONTAINS

! *****************************************************************************
!!****f* extended_system_types/lnhc_init [1.0] *
!!
!!   NAME
!!    lnhc_init 
!!
!!   FUNCTION
!!     Initialize type for Nose-Hoover thermostat
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     
!!
!!   MODIFICATION HISTORY
!!     
!!
!!*** **********************************************************************
  SUBROUTINE lnhc_init(lnhc, error)
    TYPE(lnhc_parameters_type), INTENT(OUT)  :: lnhc
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lnhc_init', &
      routineP = moduleN//':'//routineN

    NULLIFY(lnhc%dt_yosh)
    NULLIFY(lnhc%nvt)
    NULLIFY(lnhc%p_kin)
    NULLIFY(lnhc%p_scale)
    NULLIFY(lnhc%v_scale)
    NULLIFY(lnhc%s_kin)
    NULLIFY(lnhc%index)
    lnhc%dis_type=""
    lnhc%nyosh=0 
    lnhc%nc=0 
    lnhc%nhc_len=0 
    lnhc%num_nhc=0
  END SUBROUTINE lnhc_init


! *****************************************************************************
!!****f* extended_system_types/lnhc_dealloc [1.0] *
!!
!!   NAME
!!     lnhc_dealloc
!!
!!   FUNCTION
!!     Deallocate type for Nose-Hoover thermostat
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     
!!
!!   MODIFICATION HISTORY
!!     
!!
!!*** **********************************************************************
  SUBROUTINE lnhc_dealloc ( lnhc, error )
    TYPE(lnhc_parameters_type), &
      DIMENSION(:), POINTER                  :: lnhc
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lnhc_dealloc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure = .FALSE.

    IF (ASSOCIATED(lnhc)) THEN
       DO i = 1, SIZE(lnhc)
          IF ( ASSOCIATED(lnhc(i)%dt_yosh)) THEN
             DEALLOCATE (lnhc(i)%dt_yosh, STAT=stat)
             CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED (lnhc(i)%nvt)) THEN
             DEALLOCATE (lnhc(i)%nvt, STAT=stat)
             CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
          ENDIF
          IF (ASSOCIATED (lnhc(i)%p_kin)) THEN
             DEALLOCATE (lnhc(i)%p_kin, STAT=stat)
             CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED (lnhc(i)%p_scale)) THEN
             DEALLOCATE (lnhc(i)%p_scale, STAT=stat)
             CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED (lnhc(i)%v_scale)) THEN
             DEALLOCATE (lnhc(i)%v_scale, STAT=stat)
             CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED (lnhc(i)%s_kin)) THEN
             DEALLOCATE (lnhc(i)%s_kin, STAT=stat)
             CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED (lnhc(i)%index)) THEN
             DEALLOCATE (lnhc(i)%index, STAT=stat)
             CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
          END IF
       ENDDO
       DEALLOCATE (lnhc, STAT=stat)
       CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
    ENDIF

  END SUBROUTINE lnhc_dealloc

END MODULE extended_system_types
!******************************************************************************
