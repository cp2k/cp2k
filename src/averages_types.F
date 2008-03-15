!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Handles the type to compute averages during an MD
!> \author Teodoro Laino [tlaino] - 03.2008 - University of Zurich
! *****************************************************************************
MODULE averages_types
  USE cell_types,                      ONLY: cell_type
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE md_ener_types,                   ONLY: md_ener_type
  USE virial_types,                    ONLY: virial_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *****************************************************************************
  TYPE average_quantities_type
     INTEGER       :: id_nr, ref_count
     REAL(KIND=dp) :: avetemp, avepot, avekin, avevol, aveca, avecb, avecc
     REAL(KIND=dp) :: avetemp_baro, avehugoniot, avecpu
     REAL(KIND=dp) :: aveal, avebe, avega, avepress, avekinc, avetempc, avepxx
     REAL(KIND=dp) :: avetemp_qm, avekin_qm, econs
  END TYPE average_quantities_type


! *** Public subroutines and data types ***
  PUBLIC :: average_quantities_type, create_averages,  release_averages,&
       retain_averages, compute_averages

! *** Global parameters ***
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'averages_types'
  INTEGER, SAVE, PRIVATE :: last_avg_env_id=0

CONTAINS

! *****************************************************************************
!> \brief Creates averages environment
!> \author Teodoro Laino [tlaino] - 03.2008 - University of Zurich
! *****************************************************************************
  SUBROUTINE create_averages ( averages, error )
    TYPE(average_quantities_type), POINTER   :: averages
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_averages', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    ALLOCATE (averages, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    last_avg_env_id=last_avg_env_id+1
    averages%id_nr=last_avg_env_id
    averages%ref_count=1
    averages%avetemp=0.0_dp
    averages%avepot=0.0_dp
    averages%avekin=0.0_dp
    averages%avevol=0.0_dp
    averages%aveca=0.0_dp
    averages%avecb=0.0_dp
    averages%avecc=0.0_dp
    averages%avetemp_baro=0.0_dp
    averages%avehugoniot=0.0_dp
    averages%avecpu=0.0_dp
    averages%aveal=0.0_dp
    averages%avebe=0.0_dp
    averages%avega=0.0_dp
    averages%avepress=0.0_dp
    averages%avekinc=0.0_dp
    averages%avetempc=0.0_dp
    averages%avepxx=0.0_dp
    averages%avetemp_qm=0.0_dp
    averages%avekin_qm=0.0_dp
    averages%econs=0.0_dp

  END SUBROUTINE create_averages

! *****************************************************************************
!> \brief retains the given averages env
!> \author Teodoro Laino [tlaino] - 03.2008 - University of Zurich
! *****************************************************************************
  SUBROUTINE retain_averages(averages, error)
    TYPE(average_quantities_type), POINTER   :: averages
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'retain_averages', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(averages),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPrecondition(averages%ref_count>0,cp_failure_level,routineP,error,failure)
       averages%ref_count=averages%ref_count+1
    END IF
  END SUBROUTINE retain_averages

! *****************************************************************************
!> \brief releases the given averages env
!> \author Teodoro Laino [tlaino] - 03.2008 - University of Zurich
! *****************************************************************************
  SUBROUTINE release_averages(averages, error)
    TYPE(average_quantities_type), POINTER   :: averages
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'release_averages', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(averages)) THEN
       CPPrecondition(averages%ref_count>0,cp_failure_level,routineP,error,failure)
       averages%ref_count=averages%ref_count-1
       IF (averages%ref_count==0) THEN
          DEALLOCATE(averages,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF
    
  END SUBROUTINE release_averages

! *****************************************************************************
!> \brief computes the averages
!> \author Teodoro Laino [tlaino] - 03.2008 - University of Zurich
! *****************************************************************************
  SUBROUTINE compute_averages(averages, md_ener, cell,  virial, pv_scalar, pv_xx,&
       used_time, hugoniot, abc, cell_angle, delta_t, error)

    TYPE(average_quantities_type), POINTER   :: averages
    TYPE(md_ener_type), POINTER              :: md_ener
    TYPE(cell_type), POINTER                 :: cell
    TYPE(virial_type), POINTER               :: virial
    REAL(dp)                                 :: pv_scalar, pv_xx
    REAL(KIND=dp), POINTER                   :: used_time
    REAL(KIND=dp)                            :: hugoniot
    REAL(dp)                                 :: abc( 3 ), cell_angle( 3 )
    INTEGER                                  :: delta_t
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'compute_averages', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    SELECT CASE (delta_t)
    CASE (1)
       IF ( virial%pv_availability ) averages%avepress = pv_scalar
       IF ( virial%pv_availability ) averages%avepxx = pv_xx
       averages%avecpu = used_time
       averages%avehugoniot = hugoniot
       averages%avetemp_baro = md_ener%temp_baro
       averages%avepot = md_ener%epot
       averages%avekin = md_ener%ekin
       averages%avetemp = md_ener%temp_part
       averages%avekin_qm  = md_ener%ekin_qm
       averages%avetemp_qm = md_ener%temp_qm
       averages%avevol = cell%deth
       averages%aveca = abc ( 1 )
       averages%avecb = abc ( 2 )
       averages%avecc = abc ( 3 )
       averages%aveal = cell_angle ( 3 )
       averages%avebe = cell_angle ( 2 )
       averages%avega = cell_angle ( 1 )
       averages%econs =  0._dp
    CASE DEFAULT
       IF (virial%pv_availability) averages%avepress = (averages%avepress*REAL(delta_t-1,dp)+pv_scalar)/REAL(delta_t, dp)
       IF (virial%pv_availability) averages%avepxx = (averages%avepxx*REAL(delta_t-1,dp)+pv_xx)/REAL(delta_t, dp)
       averages%avecpu = (averages%avecpu*REAL(delta_t-1,dp)+used_time)/REAL(delta_t,dp)
       averages%avehugoniot = (averages%avehugoniot*REAL(delta_t-1,dp)+hugoniot)/REAL(delta_t,dp)
       averages%avetemp_baro= (averages%avetemp_baro*REAL(delta_t-1,dp)+md_ener%temp_baro)/REAL(delta_t,dp)
       averages%avepot = (averages%avepot*REAL(delta_t-1,dp)+md_ener%epot)/REAL(delta_t,dp)
       averages%avekin = (averages%avekin*REAL(delta_t-1,dp)+md_ener%ekin)/REAL(delta_t,dp)
       averages%avetemp = (averages%avetemp*REAL(delta_t-1,dp)+md_ener%temp_part)/REAL(delta_t,dp)
       averages%avekin_qm = (averages%avekin_qm*REAL(delta_t-1,dp)+md_ener%ekin_qm)/REAL(delta_t,dp)
       averages%avetemp_qm = (averages%avetemp_qm*REAL(delta_t-1,dp)+md_ener%temp_qm)/REAL(delta_t,dp)
       averages%avevol = (averages%avevol*REAL(delta_t-1,dp)+cell%deth)/REAL(delta_t,dp)
       averages%aveca = (averages%aveca*REAL(delta_t-1,dp)+abc(1))/REAL(delta_t,dp)
       averages%avecb = (averages%avecb*REAL(delta_t-1,dp)+abc(2))/REAL(delta_t,dp)
       averages%avecc = (averages%avecc*REAL(delta_t-1,dp)+abc(3))/REAL(delta_t,dp)
       averages%aveal = (averages%aveal*REAL(delta_t-1,dp)+cell_angle(3))/REAL(delta_t,dp)
       averages%avebe = (averages%avebe*REAL(delta_t-1,dp)+cell_angle(2))/REAL(delta_t,dp)
       averages%avega = (averages%avega*REAL(delta_t-1,dp)+cell_angle(1))/REAL(delta_t,dp)
       averages%econs = (averages%econs*REAL(delta_t-1,dp) + md_ener%delta_cons)/REAL(delta_t,dp)
    END SELECT

  END SUBROUTINE compute_averages

END MODULE averages_types
