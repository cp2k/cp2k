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
  USE virial_types,                    ONLY: virial_create,&
                                             virial_release,&
                                             virial_type
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
     TYPE(virial_type), POINTER :: virial
  END TYPE average_quantities_type

! *****************************************************************************
  INTERFACE get_averages
     MODULE PROCEDURE get_averages_1r, get_averages_3x3r
  END INTERFACE get_averages

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
  SUBROUTINE create_averages ( averages, virial_avg, error )
    TYPE(average_quantities_type), POINTER   :: averages
    LOGICAL, INTENT(IN), OPTIONAL            :: virial_avg
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
    NULLIFY(averages%virial)
    IF (PRESENT(virial_avg)) THEN
       IF (virial_avg) CALL virial_create(averages%virial, error)
    END IF
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
          CALL virial_release(averages%virial, error)
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
       IF (ASSOCIATED(averages%virial).AND.virial%pv_availability) THEN
          averages%virial%pv_total      = virial%pv_total
          averages%virial%pv_virial     = virial%pv_virial
          averages%virial%pv_kinetic    = virial%pv_kinetic
          averages%virial%pv_constraint = virial%pv_constraint
          averages%virial%pv_xc         = virial%pv_xc
       END IF
    CASE DEFAULT
       IF (virial%pv_availability) averages%avepress = get_averages(averages%avepress,pv_scalar,delta_t)
       IF (virial%pv_availability) averages%avepxx   = get_averages(averages%avepxx,pv_xx,delta_t)
       averages%avecpu       = get_averages(averages%avecpu,used_time,delta_t)
       averages%avehugoniot  = get_averages(averages%avehugoniot,hugoniot,delta_t)
       averages%avetemp_baro = get_averages(averages%avetemp_baro,md_ener%temp_baro,delta_t)
       averages%avepot       = get_averages(averages%avepot,md_ener%epot,delta_t)
       averages%avekin       = get_averages(averages%avekin,md_ener%ekin,delta_t)
       averages%avetemp      = get_averages(averages%avetemp,md_ener%temp_part,delta_t)
       averages%avekin_qm    = get_averages(averages%avekin_qm,md_ener%ekin_qm,delta_t)
       averages%avetemp_qm   = get_averages(averages%avetemp_qm,md_ener%temp_qm,delta_t)
       averages%avevol       = get_averages(averages%avevol,cell%deth,delta_t)
       averages%aveca        = get_averages(averages%aveca,abc(1),delta_t)
       averages%avecb        = get_averages(averages%avecb,abc(2),delta_t)
       averages%avecc        = get_averages(averages%avecc,abc(3),delta_t)
       averages%aveal        = get_averages(averages%aveal,cell_angle(3),delta_t) 
       averages%avebe        = get_averages(averages%avebe,cell_angle(2),delta_t)
       averages%avega        = get_averages(averages%avega,cell_angle(1),delta_t)
       averages%econs        = get_averages(averages%econs,md_ener%delta_cons,delta_t)
       IF (ASSOCIATED(averages%virial).AND.virial%pv_availability) THEN
          averages%virial%pv_total      = get_averages(averages%virial%pv_total,virial%pv_total,delta_t)
          averages%virial%pv_virial     = get_averages(averages%virial%pv_virial,virial%pv_virial,delta_t)
          averages%virial%pv_kinetic    = get_averages(averages%virial%pv_kinetic,virial%pv_kinetic,delta_t)
          averages%virial%pv_constraint = get_averages(averages%virial%pv_constraint,virial%pv_constraint,delta_t)
          averages%virial%pv_xc         = get_averages(averages%virial%pv_xc,virial%pv_xc,delta_t)
       END IF
    END SELECT
  END SUBROUTINE compute_averages

! *****************************************************************************
!> \brief computes the averages - low level for REAL
!> \author Teodoro Laino [tlaino] - 03.2008 - University of Zurich
! *****************************************************************************
  FUNCTION get_averages_1r(avg, add, delta_t) RESULT(res)
    REAL(KIND=dp), INTENT(INOUT)             :: avg
    REAL(KIND=dp), INTENT(IN)                :: add
    INTEGER, INTENT(IN)                      :: delta_t
    REAL(KIND=dp)                            :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'get_averages_1r', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    res = (avg*REAL(delta_t-1,dp) + add)/REAL(delta_t,dp)
  END FUNCTION get_averages_1r

! *****************************************************************************
!> \brief computes the averages - low level for REAL (3,3)
!> \author Teodoro Laino [tlaino] - 03.2008 - University of Zurich
! *****************************************************************************
  FUNCTION get_averages_3x3r(avg, add, delta_t) RESULT(res)
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(INOUT)                          :: avg
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN)                             :: add
    INTEGER, INTENT(IN)                      :: delta_t
    REAL(KIND=dp), DIMENSION(3, 3)           :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'get_averages_3x3r', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    res = (avg*REAL(delta_t-1,dp) + add)/REAL(delta_t,dp)
  END FUNCTION get_averages_3x3r

END MODULE averages_types
