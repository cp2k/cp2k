!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
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
!!
!!   SOURCE
!******************************************************************************

MODULE extended_system_types
  USE kinds,                           ONLY: dp
  USE termination,                     ONLY: stop_program

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: extended_system_type, extended_type_init, &
            lnhc_dealloc, lnhc_init, npt_dealloc, &
            lnhc_parameters_type, npt_info_type, nhc_info_type, &
            particle_kinetic_type, particle_scale_type

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
     CHARACTER ( LEN = 10 ) :: dis_type
     INTEGER :: nyosh, nc, nhc_len, num_nhc
     REAL (KIND=dp), POINTER :: dt_yosh ( : )
     TYPE ( nhc_info_type ), POINTER :: nvt ( :, : )
     TYPE ( particle_kinetic_type ), POINTER :: p_kin ( :, : )
     TYPE ( particle_scale_type ), POINTER :: p_scale ( :, : )
     REAL (KIND=dp), POINTER :: v_scale ( : )
     REAL (KIND=dp), POINTER :: s_kin ( : )
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

  CONTAINS

!***********
!  nullify all pointers 
!***********
  SUBROUTINE extended_type_init(extended_type)
    TYPE(extended_system_type), INTENT(OUT) :: extended_type
    NULLIFY(extended_type%nhc_part)
    NULLIFY(extended_type%nhc_baro)
    NULLIFY(extended_type%nhc_coef)
    NULLIFY(extended_type%npt_info)
  END SUBROUTINE

  SUBROUTINE lnhc_init(lnhc)
    TYPE(lnhc_parameters_type), INTENT(OUT) :: lnhc
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
  END SUBROUTINE

!******
!  deallocate lnhc_dealloc
!******
SUBROUTINE lnhc_dealloc ( lnhc )
   TYPE ( lnhc_parameters_type ), POINTER, DIMENSION(:) :: lnhc

   INTEGER :: isos,i

   IF (ASSOCIATED(lnhc)) THEN
    DO I=1,SIZE(lnhc)
      IF ( ASSOCIATED ( lnhc(I) % dt_yosh ) ) THEN
        DEALLOCATE ( lnhc(I) % dt_yosh , STAT = isos )
        IF (isos/=0) CALL stop_program('extended_system_type_destroy', &
         'failed to allocate lnhc % dt_yosh' )
      END IF
      IF ( ASSOCIATED ( lnhc(I) % nvt ) ) THEN
        DEALLOCATE ( lnhc(I) % nvt , STAT = isos )
        IF (isos/=0) CALL stop_program('extended_system_type_destroy', &
         'failed to allocate lnhc % nvt' )
      ENDIF
      IF ( ASSOCIATED ( lnhc(I) % p_kin ) ) THEN
        DEALLOCATE ( lnhc(I) % p_kin , STAT = isos )
        IF (isos/=0) CALL stop_program('extended_system_type_destroy', &
         'failed to allocate lnhc % p_kin' )
      END IF
      IF ( ASSOCIATED ( lnhc(I) % p_scale ) ) THEN
        DEALLOCATE ( lnhc(I) % p_scale , STAT = isos )
        IF (isos/=0) CALL stop_program('extended_system_type_destroy', &
         'failed to allocate lnhc % p_scale' )
      END IF
      IF ( ASSOCIATED ( lnhc(I) % v_scale ) ) THEN
        DEALLOCATE ( lnhc(I) % v_scale , STAT = isos )
        IF (isos/=0) CALL stop_program('extended_system_type_destroy', &
         'failed to allocate lnhc % v_scale' )
      END IF
      IF ( ASSOCIATED ( lnhc(I) % s_kin ) ) THEN
        DEALLOCATE ( lnhc(I) % s_kin , STAT = isos )
        IF (isos/=0) CALL stop_program('extended_system_type_destroy', &
         'failed to allocate lnhc % s_kin' )
      END IF
      IF ( ASSOCIATED ( lnhc(I) % index ) ) THEN
        DEALLOCATE ( lnhc(I) % index , STAT = isos )
        IF (isos/=0) CALL stop_program('extended_system_type_destroy', &
         'failed to allocate lnhc % index' )
      END IF
    ENDDO
    DEALLOCATE ( lnhc , STAT = isos )
    IF (isos/=0) CALL stop_program('extended_system_type_destroy', &
         'failed to allocate lnhc ' )
   ENDIF
 
 END SUBROUTINE lnhc_dealloc

 SUBROUTINE npt_dealloc ( npt_info )
   TYPE ( npt_info_type ), POINTER :: npt_info(:,:)
   INTEGER :: isos

   IF ( ASSOCIATED ( npt_info ) ) THEN
        DEALLOCATE ( npt_info , STAT = isos )
        IF (isos/=0) CALL stop_program('extended_system_type_destroy', &
         'failed to allocate npt_info' )
  END IF

 END SUBROUTINE npt_dealloc

!!*****
!******************************************************************************
END MODULE extended_system_types
!******************************************************************************
