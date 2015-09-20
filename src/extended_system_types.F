!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Lumps all possible extended system variables into one
!>    type for easy access and passing
!> \par History
!>      Teodoro Laino - 09.2007 - University of Zurich
!>        Cleaned the typo.. no need to have an extended
!>        type. Thermostat and Barostat type have been created
!> \author CJM
! *****************************************************************************
MODULE extended_system_types
  USE bibliography,                    ONLY: Nose1984a,&
                                             Nose1984b,&
                                             cite_reference
  USE input_constants,                 ONLY: do_thermo_no_communication
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: lnhc_dealloc,&
            lnhc_init,&
            lnhc_parameters_type,&
            map_info_type,&
            create_map_info_type,&
            release_map_info_type,&
            npt_info_type

  LOGICAL, PARAMETER, PUBLIC :: debug_isotropic_limit = .FALSE.
  LOGICAL, PARAMETER, PUBLIC :: debug_uniaxial_limit  = .FALSE.

! *****************************************************************************
  TYPE npt_info_type
     REAL (KIND=dp) :: eps
     REAL (KIND=dp) :: v
     REAL (KIND=dp) :: f
     REAL (KIND=dp) :: mass
  END TYPE npt_info_type

! *****************************************************************************
  TYPE nhc_info_type
     INTEGER        :: degrees_of_freedom
     REAL (KIND=dp) :: eta
     REAL (KIND=dp) :: v
     REAL (KIND=dp) :: f
     REAL (KIND=dp) :: nkt
     REAL (KIND=dp) :: mass
  END TYPE nhc_info_type

! *****************************************************************************
  TYPE point_info_type
     REAL (KIND=dp), POINTER :: point
  END TYPE point_info_type

! *****************************************************************************
  TYPE map_info_type
     INTEGER                                         :: dis_type
     INTEGER, POINTER, DIMENSION(:)                  :: index, map_index
     REAL (KIND=dp), POINTER, DIMENSION(:)           :: v_scale
     REAL (KIND=dp), POINTER, DIMENSION(:)           :: s_kin
     TYPE (point_info_type), POINTER, DIMENSION(:,:) :: p_scale
     TYPE (point_info_type), POINTER, DIMENSION(:,:) :: p_kin
  END TYPE map_info_type

! *****************************************************************************
  TYPE lnhc_parameters_type
     INTEGER                                 :: nyosh, nc, nhc_len
     INTEGER                                 :: glob_num_nhc, loc_num_nhc, region
     REAL(KIND=dp)                           :: tau_nhc, dt_fact
     REAL(KIND=dp), POINTER                  :: dt_yosh(:)
     TYPE(nhc_info_type), POINTER            :: nvt(:,:)
     TYPE(map_info_type), POINTER            :: map_info
  END TYPE lnhc_parameters_type

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'extended_system_types'

CONTAINS

! *****************************************************************************
!> \brief Initialize type for Nose-Hoover thermostat
!> \param lnhc ...
!> \param section ...
! *****************************************************************************
  SUBROUTINE lnhc_init(lnhc, section)
    TYPE(lnhc_parameters_type), POINTER      :: lnhc
    TYPE(section_vals_type), POINTER         :: section

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lnhc_init', &
      routineP = moduleN//':'//routineN

    NULLIFY(lnhc%dt_yosh)
    NULLIFY(lnhc%nvt)
    NULLIFY(lnhc%map_info)
    lnhc%loc_num_nhc=0
    lnhc%glob_num_nhc=0
    lnhc%dt_fact=1.0_dp
    CALL cite_reference(Nose1984a)
    CALL cite_reference(Nose1984b)
    CALL section_vals_val_get(section,"LENGTH",i_val=lnhc%nhc_len)
    CALL section_vals_val_get(section,"YOSHIDA",i_val=lnhc%nyosh)
    CALL section_vals_val_get(section,"TIMECON",r_val=lnhc%tau_nhc)
    CALL section_vals_val_get(section,"MTS",i_val=lnhc%nc)
    CALL create_map_info_type(lnhc%map_info)

  END SUBROUTINE lnhc_init

! *****************************************************************************
!> \brief create the map_info type
!> \param map_info ...
! *****************************************************************************
  SUBROUTINE create_map_info_type(map_info)
    TYPE(map_info_type), POINTER             :: map_info

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_map_info_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    ALLOCATE(map_info, stat=stat)
    CPASSERT(stat==0)
    NULLIFY(map_info%index, map_info%map_index)
    NULLIFY(map_info%v_scale)
    NULLIFY(map_info%p_scale)
    NULLIFY(map_info%s_kin)
    NULLIFY(map_info%p_kin)
    map_info%dis_type=do_thermo_no_communication

  END SUBROUTINE create_map_info_type

! *****************************************************************************
!> \brief release the map_info type
!> \param map_info ...
! *****************************************************************************
  SUBROUTINE release_map_info_type(map_info)
    TYPE(map_info_type), POINTER             :: map_info

    CHARACTER(LEN=*), PARAMETER :: routineN = 'release_map_info_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.

    IF (ASSOCIATED(map_info)) THEN
       IF (ASSOCIATED (map_info%p_kin)) THEN
          DEALLOCATE (map_info%p_kin, STAT=stat)
          CPASSERT(stat==0)
       END IF
       IF (ASSOCIATED (map_info%p_scale)) THEN
          DEALLOCATE (map_info%p_scale, STAT=stat)
          CPASSERT(stat==0)
       END IF
       IF (ASSOCIATED (map_info%v_scale)) THEN
          DEALLOCATE (map_info%v_scale, STAT=stat)
          CPASSERT(stat==0)
       END IF
       IF (ASSOCIATED (map_info%s_kin)) THEN
          DEALLOCATE (map_info%s_kin, STAT=stat)
          CPASSERT(stat==0)
       END IF
       IF (ASSOCIATED (map_info%index)) THEN
          DEALLOCATE (map_info%index, STAT=stat)
          CPASSERT(stat==0)
       END IF
       IF (ASSOCIATED (map_info%map_index)) THEN
          DEALLOCATE (map_info%map_index, STAT=stat)
          CPASSERT(stat==0)
       END IF

       DEALLOCATE (map_info, STAT=stat)
       CPASSERT(stat==0)
    END IF

  END SUBROUTINE release_map_info_type

! *****************************************************************************
!> \brief Deallocate type for Nose-Hoover thermostat
!> \param lnhc ...
! *****************************************************************************
  SUBROUTINE lnhc_dealloc ( lnhc)
    TYPE(lnhc_parameters_type), POINTER      :: lnhc

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lnhc_dealloc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.

    IF (ASSOCIATED(lnhc)) THEN
       IF ( ASSOCIATED(lnhc%dt_yosh)) THEN
          DEALLOCATE (lnhc%dt_yosh, STAT=stat)
          CPASSERT(stat==0)
       END IF
       IF (ASSOCIATED (lnhc%nvt)) THEN
          DEALLOCATE (lnhc%nvt, STAT=stat)
          CPASSERT(stat==0)
       ENDIF
       CALL release_map_info_type(lnhc%map_info)
       DEALLOCATE (lnhc, STAT=stat)
       CPASSERT(stat==0)
    ENDIF

  END SUBROUTINE lnhc_dealloc

END MODULE extended_system_types
