!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Thermal regions type: to initialize and control the temperature of
!>        different regions
!> \par History
!>   - Added support for langevin regions (2014/01/08, LT)
!> \author MI
! *****************************************************************************
MODULE thermal_region_types

  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: thermal_regions_type,&
            thermal_region_type,&
            allocate_thermal_regions,&
            release_thermal_regions,&
            retain_thermal_regions

  TYPE thermal_regions_type
     INTEGER :: id_nr, ref_count, nregions
     LOGICAL :: force_rescaling
     REAL(KIND=dp) :: temp_reg0
     LOGICAL, DIMENSION(:), POINTER                   :: do_langevin
     TYPE(section_vals_type), POINTER                 :: section
     TYPE(thermal_region_type), DIMENSION(:), POINTER :: thermal_region
  END TYPE thermal_regions_type

  TYPE thermal_region_type
     INTEGER :: region_index, npart
     INTEGER, DIMENSION(:), POINTER :: part_index
     REAL(KIND=dp) :: ekin, temperature, temp_expected, temp_tol
  END TYPE thermal_region_type

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'thermal_region_types'
CONTAINS

! *****************************************************************************
!> \brief allocate thermal_regions
!> \param thermal_regions ...
!> \param error ...
!> \author
! *****************************************************************************
  SUBROUTINE allocate_thermal_regions(thermal_regions,error)
    TYPE(thermal_regions_type), POINTER      :: thermal_regions
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_thermal_regions', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: check, failure

    failure = .FALSE.
    check = .NOT.ASSOCIATED(thermal_regions)
    CPPrecondition(check,cp_failure_level,routineP,error,failure)

    ALLOCATE(thermal_regions,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    thermal_regions%ref_count =  1
    thermal_regions%nregions = 0
    NULLIFY(thermal_regions%thermal_region)
    NULLIFY(thermal_regions%do_langevin)

  END SUBROUTINE allocate_thermal_regions
! *****************************************************************************
!> \brief retains  thermal_regions
!> \param thermal_regions ...
!> \param error ...
!> \author
! *****************************************************************************
  SUBROUTINE retain_thermal_regions(thermal_regions,error)

    TYPE(thermal_regions_type), POINTER      :: thermal_regions
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'retain_thermal_regions', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    IF (.NOT. failure) THEN
       IF (ASSOCIATED(thermal_regions)) THEN
          CPPrecondition(thermal_regions%ref_count>0,cp_failure_level,routineP,error,failure)
          thermal_regions%ref_count=thermal_regions%ref_count+1
       END IF
    END IF

  END SUBROUTINE retain_thermal_regions

! *****************************************************************************
!> \brief release thermal_regions
!> \param thermal_regions ...
!> \param error ...
!> \author
! *****************************************************************************
  SUBROUTINE release_thermal_regions(thermal_regions, error)

    TYPE(thermal_regions_type), POINTER      :: thermal_regions
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'release_thermal_regions', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ireg, stat
    LOGICAL                                  :: check, failure

    failure =.FALSE.
    check = ASSOCIATED(thermal_regions)
    IF (check) THEN
       check = thermal_regions%ref_count>0
       CPPrecondition(check,cp_failure_level,routineP,error,failure)
       thermal_regions%ref_count=thermal_regions%ref_count-1
       IF (thermal_regions%ref_count<1) THEN
          IF (ASSOCIATED(thermal_regions%thermal_region)) THEN
            DO ireg = 1,SIZE(thermal_regions%thermal_region)
              DEALLOCATE(thermal_regions%thermal_region(ireg)%part_index,STAT=stat)
              CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
            END DO
            DEALLOCATE(thermal_regions%thermal_region,STAT=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(thermal_regions%do_langevin)) THEN
             DEALLOCATE(thermal_regions%do_langevin, STAT=stat)
             CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          END IF
          DEALLOCATE(thermal_regions, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF

  END SUBROUTINE release_thermal_regions

END MODULE thermal_region_types
