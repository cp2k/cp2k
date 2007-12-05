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
  USE bibliography,                    ONLY: Nose1984a,&
                                             Nose1984b,&
                                             cite_reference
  USE input_constants,                 ONLY: do_thermo_no_communication
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE f77_blas
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: lnhc_dealloc,&
            lnhc_init,&
            lnhc_parameters_type,&
            map_info_type,&
            create_map_info_type,&
            release_map_info_type,&
            npt_info_type,&
            nhc_info_type,&
            point_info_type,&
            extended_energy_type

  LOGICAL, PARAMETER, PUBLIC :: debug_isotropic_limit = .FALSE.
  LOGICAL, PARAMETER, PUBLIC :: debug_uniaxial_limit  = .FALSE.

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

  TYPE point_info_type
     REAL (KIND=dp), POINTER :: point
  END TYPE point_info_type

  TYPE map_info_type
     INTEGER                                         :: dis_type
     INTEGER, POINTER, DIMENSION(:)                  :: index, map_index
     REAL (KIND=dp), POINTER, DIMENSION(:)           :: v_scale
     REAL (KIND=dp), POINTER, DIMENSION(:)           :: s_kin
     TYPE (point_info_type), POINTER, DIMENSION(:,:) :: p_scale
     TYPE (point_info_type), POINTER, DIMENSION(:,:) :: p_kin
  END TYPE map_info_type

  TYPE lnhc_parameters_type
     INTEGER                                 :: nyosh, nc, nhc_len
     INTEGER                                 :: glob_num_nhc, loc_num_nhc, region
     REAL(KIND=dp)                           :: tau_nhc
     REAL(KIND=dp), POINTER                  :: dt_yosh(:)
     TYPE(nhc_info_type), POINTER            :: nvt(:,:)
     TYPE(map_info_type), POINTER            :: map_info
  END TYPE lnhc_parameters_type

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
  SUBROUTINE lnhc_init(lnhc, section, error)
    TYPE(lnhc_parameters_type), POINTER      :: lnhc
    TYPE(section_vals_type), POINTER         :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lnhc_init', &
      routineP = moduleN//':'//routineN

    NULLIFY(lnhc%dt_yosh)
    NULLIFY(lnhc%nvt)
    NULLIFY(lnhc%map_info)
    lnhc%loc_num_nhc=0
    lnhc%glob_num_nhc=0
    CALL cite_reference(Nose1984a)
    CALL cite_reference(Nose1984b)
    CALL section_vals_val_get(section,"LENGTH",i_val=lnhc%nhc_len,error=error)
    CALL section_vals_val_get(section,"YOSHIDA",i_val=lnhc%nyosh,error=error)
    CALL section_vals_val_get(section,"TIMECON",r_val=lnhc%tau_nhc,error=error)
    CALL section_vals_val_get(section,"MTS",i_val=lnhc%nc,error=error)
    CALL create_map_info_type(lnhc%map_info, error)

  END SUBROUTINE lnhc_init

! *****************************************************************************
!!****f* extended_system_types/create_map_info_type [1.0] *
!!
!!   NAME
!!     create_map_info_type
!!
!!   FUNCTION
!!     create the map_info type
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
  SUBROUTINE create_map_info_type(map_info, error)
    TYPE(map_info_type), POINTER             :: map_info
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_map_info_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    ALLOCATE(map_info, stat=stat)
    CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
    NULLIFY(map_info%index, map_info%map_index)
    NULLIFY(map_info%v_scale)
    NULLIFY(map_info%p_scale)
    NULLIFY(map_info%s_kin)
    NULLIFY(map_info%p_kin)
    map_info%dis_type=do_thermo_no_communication

  END SUBROUTINE create_map_info_type

! *****************************************************************************
!!****f* extended_system_types/release_map_info_type [1.0] *
!!
!!   NAME
!!     release_map_info_type
!!
!!   FUNCTION
!!     release the map_info type
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
  SUBROUTINE release_map_info_type(map_info, error)
    TYPE(map_info_type), POINTER             :: map_info
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'release_map_info_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.

    IF (ASSOCIATED(map_info)) THEN
       IF (ASSOCIATED (map_info%p_kin)) THEN
          DEALLOCATE (map_info%p_kin, STAT=stat)
          CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED (map_info%p_scale)) THEN
          DEALLOCATE (map_info%p_scale, STAT=stat)
          CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED (map_info%v_scale)) THEN
          DEALLOCATE (map_info%v_scale, STAT=stat)
          CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED (map_info%s_kin)) THEN
          DEALLOCATE (map_info%s_kin, STAT=stat)
          CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED (map_info%index)) THEN
          DEALLOCATE (map_info%index, STAT=stat)
          CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED (map_info%map_index)) THEN
          DEALLOCATE (map_info%map_index, STAT=stat)
          CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
       END IF
       
       DEALLOCATE (map_info, STAT=stat)
       CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
    END IF
    
  END SUBROUTINE release_map_info_type

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
    TYPE(lnhc_parameters_type), POINTER      :: lnhc
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lnhc_dealloc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.

    IF (ASSOCIATED(lnhc)) THEN
       IF ( ASSOCIATED(lnhc%dt_yosh)) THEN
          DEALLOCATE (lnhc%dt_yosh, STAT=stat)
          CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED (lnhc%nvt)) THEN
          DEALLOCATE (lnhc%nvt, STAT=stat)
          CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
       ENDIF
       CALL release_map_info_type(lnhc%map_info, error)
       DEALLOCATE (lnhc, STAT=stat)
       CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
    ENDIF

  END SUBROUTINE lnhc_dealloc

END MODULE extended_system_types
!******************************************************************************
