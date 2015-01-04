!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Type for the canonical sampling through velocity rescaling     
!> \author Teodoro Laino - 09.2007 University of Zurich [tlaino]
! *****************************************************************************
MODULE al_system_types
  USE bibliography,                    ONLY: Jones2011,&
                                             cite_reference
  USE extended_system_types,           ONLY: create_map_info_type,&
                                             map_info_type,&
                                             release_map_info_type
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE simpar_types,                    ONLY: simpar_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: al_system_type,&
            al_init,&
            al_dealloc,&
            al_thermo_create

! *****************************************************************************
  TYPE al_thermo_type
     INTEGER                                 :: degrees_of_freedom
     REAL(KIND=dp)                           :: nkt
     REAL(KIND=dp)                           :: chi
     REAL(KIND=dp)                           :: mass
     REAL(KIND=dp)                           :: region_kin_energy
  END TYPE al_thermo_type

! *****************************************************************************
  TYPE al_system_type
     INTEGER                                 :: region, glob_num_al, loc_num_al
     REAL(KIND=dp)                           :: tau_nh, tau_langevin, dt_fact
     REAL(KIND=dp)                           :: dt
     TYPE(al_thermo_type), POINTER           :: nvt(:)
     TYPE(map_info_type), POINTER            :: map_info
  END TYPE al_system_type

! *** Global parameters ***
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'al_system_types'

CONTAINS

! *****************************************************************************
!> \brief Initialize type for Adaptive Langevin (AD_LANGEVIN)
!> \param al ...
!> \param simpar ...
!> \param section ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Noam Bernstein [noamb] 02.2012
! *****************************************************************************
  SUBROUTINE al_init(al, simpar, section, error)
    TYPE(al_system_type), POINTER            :: al
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(section_vals_type), POINTER         :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'al_init', &
      routineP = moduleN//':'//routineN

    NULLIFY(al%nvt)
    NULLIFY(al%map_info)
    al%loc_num_al=0
    al%glob_num_al=0
    al%dt_fact=1.0_dp
    al%dt=simpar%dt
    CALL cite_reference(Jones2011)
    CALL section_vals_val_get(section,"TIMECON_NH",r_val=al%tau_nh,error=error)
    CALL section_vals_val_get(section,"TIMECON_LANGEVIN",r_val=al%tau_langevin,error=error)
    CALL create_map_info_type(al%map_info, error)
    
  END SUBROUTINE al_init

! *****************************************************************************
!> \brief Initialize NVT type for AD_LANGEVIN thermostat
!> \param al ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Noam Bernstein [noamb]  02.2012
! *****************************************************************************
  SUBROUTINE al_thermo_create(al, error)
    TYPE(al_system_type), POINTER            :: al
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'al_thermo_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: seed

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(al),cp_fatal_level,routineP,error,failure)
    CPPrecondition(.NOT.ASSOCIATED(al%nvt),cp_fatal_level,routineP,error,failure)
    
    ALLOCATE ( al%nvt(al%loc_num_al),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO i = 1, al%loc_num_al
       al%nvt(i)%chi = 0.0_dp
    END DO
    ! Initialize the gaussian stream random number
    ALLOCATE (seed(3,2,al%glob_num_al),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE al_thermo_create

! *****************************************************************************
!> \brief Deallocate type for AD_LANGEVIN thermostat
!> \param al ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Noam Bernstein [noamb] 02.2012
! *****************************************************************************
  SUBROUTINE al_dealloc ( al, error )
    TYPE(al_system_type), POINTER            :: al
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'al_dealloc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(al)) THEN
       CALL al_thermo_dealloc(al%nvt, error)
       CALL release_map_info_type(al%map_info, error)
       DEALLOCATE (al, STAT=stat)
       CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
    ENDIF

  END SUBROUTINE al_dealloc

! *****************************************************************************
!> \brief Deallocate NVT type for AD_LANGEVIN thermostat
!> \param nvt ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Noam Bernstein [noamb] 02.2012
! *****************************************************************************
  SUBROUTINE al_thermo_dealloc ( nvt, error )
    TYPE(al_thermo_type), DIMENSION(:), &
      POINTER                                :: nvt
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'al_thermo_dealloc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(nvt)) THEN
       DEALLOCATE (nvt, STAT=stat)
       CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
    ENDIF
  END SUBROUTINE al_thermo_dealloc
  
END MODULE al_system_types

