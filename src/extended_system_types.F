!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
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
!!
!!   SOURCE
!******************************************************************************

MODULE extended_system_types
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE kinds,                           ONLY: dbl

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: extended_system_type,& ! extended_system_type_destroy, &
            lnhc_parameters_type, npt_info_type, nhc_info_type, &
            particle_kinetic_type, particle_scale_type

  TYPE npt_info_type
     REAL ( dbl ) :: eps
     REAL ( dbl ) :: v
     REAL ( dbl ) :: f
     REAL ( dbl ) :: mass
  END TYPE npt_info_type

  TYPE nhc_info_type
     REAL ( dbl ) :: eta
     REAL ( dbl ) :: v
     REAL ( dbl ) :: f
     REAL ( dbl ) :: nkt
     REAL ( dbl ) :: mass
  END TYPE nhc_info_type

  TYPE particle_kinetic_type
     REAL ( dbl ), POINTER :: point
  END TYPE particle_kinetic_type

  TYPE particle_scale_type
     REAL ( dbl ), POINTER :: point
  END TYPE particle_scale_type

  TYPE lnhc_parameters_type
     CHARACTER ( LEN = 10 ) :: dis_type
     INTEGER :: nyosh, nc, nhc_len, num_nhc
     REAL ( dbl ), POINTER :: dt_yosh ( : )
     TYPE ( nhc_info_type ), POINTER :: nvt ( :, : )
     TYPE ( particle_kinetic_type ), POINTER :: p_kin ( :, : )
     TYPE ( particle_scale_type ), POINTER :: p_scale ( :, : )
     REAL ( dbl ), POINTER :: v_scale ( : )
     REAL ( dbl ), POINTER :: s_kin ( : )
     INTEGER, POINTER :: index ( : )
  END TYPE lnhc_parameters_type

  TYPE extended_system_type
     !PRIVATE
     INTEGER :: id_nr, ref_count, in_use
     TYPE ( lnhc_parameters_type ), POINTER :: nhc_part  ( : )
     TYPE ( lnhc_parameters_type ), POINTER :: nhc_baro ( : )
     TYPE ( lnhc_parameters_type ), POINTER :: nhc_coef (:)
     TYPE ( npt_info_type ), POINTER :: npt_info ( :, : )
  END TYPE extended_system_type

! *** Global parameters ***

  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='extended_systems_types'

!  CONTAINS

!!*****
!******************************************************************************
!!****** extended_system_types/extended_system_type_destroy [1.0] *
!!
!!   NAME
!!     extended_system_type_destroy
!!
!!   FUNCTION
!!     Kills lnhc_paramters_type and npt_info_type
!!
!!   AUTHOR
!!
!!   CREATION DATE
!!     14-NOV-2000
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
!SUBROUTINE extended_system_type_destroy ( lnhc, npt_info )
!  IMPLICIT NONE
!  TYPE ( lnhc_parameters_type ), POINTER, OPTIONAL :: lnhc
!  TYPE ( npt_info_type ), POINTER, OPTIONAL :: npt_info

!  INTEGER :: isos

!  IF ( PRESENT ( lnhc ) ) THEN
!    IF ( ASSOCIATED ( lnhc % dt_yosh ) ) THEN
!      DEALLOCATE ( lnhc % dt_yosh , STAT = isos )
!      IF (isos/=0) CALL stop_program('extended_system_type_destroy', &
!       'failed to allocate lnhc % dt_yosh' )
!      NULLIFY ( lnhc % dt_yosh )
!    END IF
!    IF ( ASSOCIATED ( lnhc % nvt ) ) THEN
!      DEALLOCATE ( lnhc % nvt , STAT = isos )
!      IF (isos/=0) CALL stop_program('extended_system_type_destroy', &
!       'failed to allocate lnhc % nvt' )
!      NULLIFY ( lnhc % nvt )
!    END IF
!    IF ( ASSOCIATED ( lnhc % p_kin ) ) THEN
!      DEALLOCATE ( lnhc % p_kin , STAT = isos )
!      IF (isos/=0) CALL stop_program('extended_system_type_destroy', &
!       'failed to allocate lnhc % p_kin' )
!      NULLIFY ( lnhc % p_kin )
!    END IF
!    IF ( ASSOCIATED ( lnhc % p_scale ) ) THEN
!      DEALLOCATE ( lnhc % p_scale , STAT = isos )
!      IF (isos/=0) CALL stop_program('extended_system_type_destroy', &
!       'failed to allocate lnhc % p_scale' )
!      NULLIFY ( lnhc % p_scale )
!    END IF
!    IF ( ASSOCIATED ( lnhc % v_scale ) ) THEN
!      DEALLOCATE ( lnhc % v_scale , STAT = isos )
!      IF (isos/=0) CALL stop_program('extended_system_type_destroy', &
!       'failed to allocate lnhc % v_scale' )
!      NULLIFY ( lnhc % v_scale )
!    END IF
!    IF ( ASSOCIATED ( lnhc % s_kin ) ) THEN
!      DEALLOCATE ( lnhc % s_kin , STAT = isos )
!      IF (isos/=0) CALL stop_program('extended_system_type_destroy', &
!       'failed to allocate lnhc % s_kin' )
!      NULLIFY ( lnhc % s_kin )
!    END IF
!    NULLIFY ( lnhc )
!  END IF

!  IF ( PRESENT ( npt_info ) ) THEN
!    IF ( ASSOCIATED ( npt_info ) ) THEN
!      DEALLOCATE ( npt_info , STAT = isos )
!      IF (isos/=0) CALL stop_program('extended_system_type_destroy', &
!       'failed to allocate npt_info' )
!      NULLIFY ( npt_info )
!    END IF
!  END IF
  
!END SUBROUTINE extended_system_type_destroy

!!*****
!******************************************************************************
END MODULE extended_system_types
!******************************************************************************
