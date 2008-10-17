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
  USE colvar_utils,                    ONLY: get_clv_force,&
                                             number_of_colvar
  USE f77_blas
  USE force_env_types,                 ONLY: force_env_type
  USE input_section_types,             ONLY: section_vals_type
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
     INTEGER                                       :: id_nr, ref_count
     TYPE(section_vals_type), POINTER              :: averages_section
     REAL(KIND=dp)                                 :: avetemp, avepot, avekin,&
                                                      avevol, aveca, avecb, avecc
     REAL(KIND=dp)                                 :: avetemp_baro, avehugoniot, avecpu
     REAL(KIND=dp)                                 :: aveal, avebe, avega, avepress,&
                                                      avekinc, avetempc, avepxx
     REAL(KIND=dp)                                 :: avetemp_qm, avekin_qm, econs
     TYPE(virial_type), POINTER                    :: virial
     ! Colvar 
     REAL(KIND=dp), POINTER, DIMENSION(:)          :: avecolvar
     REAL(KIND=dp), POINTER, DIMENSION(:)          :: aveMmatrix
  END TYPE average_quantities_type

! *****************************************************************************
  INTERFACE get_averages
     MODULE PROCEDURE get_averages_rs, get_averages_rv, get_averages_rm
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
  SUBROUTINE create_averages(averages, averages_section, virial_avg, force_env, error)
    TYPE(average_quantities_type), POINTER   :: averages
    TYPE(section_vals_type), POINTER         :: averages_section
    LOGICAL, INTENT(IN), OPTIONAL            :: virial_avg
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_averages', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, nint, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    ALLOCATE (averages, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    NULLIFY(averages%virial)
    NULLIFY(averages%avecolvar)
    NULLIFY(averages%aveMmatrix)
    ! Point to the averages section
    averages%averages_section => averages_section
    ! Initialize averages
    last_avg_env_id=last_avg_env_id+1
    averages%id_nr=last_avg_env_id
    averages%ref_count      = 1
    averages%avetemp        = 0.0_dp
    averages%avepot         = 0.0_dp
    averages%avekin         = 0.0_dp
    averages%avevol         = 0.0_dp
    averages%aveca          = 0.0_dp
    averages%avecb          = 0.0_dp
    averages%avecc          = 0.0_dp
    averages%avetemp_baro   = 0.0_dp
    averages%avehugoniot    = 0.0_dp
    averages%avecpu         = 0.0_dp
    averages%aveal          = 0.0_dp
    averages%avebe          = 0.0_dp
    averages%avega          = 0.0_dp
    averages%avepress       = 0.0_dp
    averages%avekinc        = 0.0_dp
    averages%avetempc       = 0.0_dp
    averages%avepxx         = 0.0_dp
    averages%avetemp_qm     = 0.0_dp
    averages%avekin_qm      = 0.0_dp
    averages%econs          = 0.0_dp
    IF (PRESENT(virial_avg)) THEN
       IF (virial_avg) CALL virial_create(averages%virial, error)
    END IF
    ! Total number of COLVARs
    nint = number_of_colvar(force_env, error=error)
    ALLOCATE(averages%avecolvar(nint), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(averages%aveMmatrix(nint*nint), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO i = 1, nint
       averages%avecolvar(i)         = 0.0_dp
    END DO
    DO i = 1, nint*nint
       averages%aveMmatrix(i)        = 0.0_dp
    END DO
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
          DEALLOCATE(averages%avecolvar, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(averages%aveMmatrix, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(averages,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF
    
  END SUBROUTINE release_averages

! *****************************************************************************
!> \brief computes the averages
!> \author Teodoro Laino [tlaino] - 03.2008 - University of Zurich
! *****************************************************************************
  SUBROUTINE compute_averages (averages, force_env, md_ener, cell,  virial, &
       pv_scalar, pv_xx, used_time, hugoniot, abc, cell_angle, delta_t, nat,&
       error)
    TYPE(average_quantities_type), POINTER   :: averages
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(md_ener_type), POINTER              :: md_ener
    TYPE(cell_type), POINTER                 :: cell
    TYPE(virial_type), POINTER               :: virial
    REAL(dp)                                 :: pv_scalar, pv_xx
    REAL(KIND=dp), POINTER                   :: used_time
    REAL(KIND=dp)                            :: hugoniot
    REAL(dp)                                 :: abc( 3 ), cell_angle( 3 )
    INTEGER, INTENT(IN)                      :: delta_t, nat
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'compute_averages', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nint, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: cvalues, Mmatrix

    failure = .FALSE.
    NULLIFY(cvalues, Mmatrix)
    nint = SIZE(averages%avecolvar)
    SELECT CASE (delta_t)
    CASE (1)
       IF (virial%pv_availability) averages%avepress = pv_scalar
       IF (virial%pv_availability) averages%avepxx = pv_xx
       averages%avecpu        = used_time
       averages%avehugoniot   = hugoniot
       averages%avetemp_baro  = md_ener%temp_baro
       averages%avepot        = md_ener%epot
       averages%avekin        = md_ener%ekin
       averages%avetemp       = md_ener%temp_part
       averages%avekin_qm     = md_ener%ekin_qm
       averages%avetemp_qm    = md_ener%temp_qm
       averages%avevol        = cell%deth
       averages%aveca         = abc(1)
       averages%avecb         = abc(2)
       averages%avecc         = abc(3)
       averages%aveal         = cell_angle(3)
       averages%avebe         = cell_angle(2)
       averages%avega         = cell_angle(1)
       averages%econs         = 0._dp
       ! Virial
       IF (ASSOCIATED(averages%virial).AND.virial%pv_availability) THEN
          averages%virial%pv_total      = virial%pv_total
          averages%virial%pv_virial     = virial%pv_virial
          averages%virial%pv_kinetic    = virial%pv_kinetic
          averages%virial%pv_constraint = virial%pv_constraint
          averages%virial%pv_xc         = virial%pv_xc
       END IF
       ! Colvars
       IF (nint>0) THEN
          CALL get_clv_force(force_env, nsize_xyz=nat*3, nsize_int=nint, cvalues=averages%avecolvar,&
                             Mmatrix=averages%aveMmatrix, error=error)
       END IF
    CASE DEFAULT
       IF (virial%pv_availability) CALL get_averages(averages%avepress,pv_scalar,delta_t,error)
       IF (virial%pv_availability) CALL get_averages(averages%avepxx,pv_xx,delta_t,error)
       CALL get_averages(averages%avecpu,used_time,delta_t,error)
       CALL get_averages(averages%avehugoniot,hugoniot,delta_t,error)
       CALL get_averages(averages%avetemp_baro,md_ener%temp_baro,delta_t,error)
       CALL get_averages(averages%avepot,md_ener%epot,delta_t,error)
       CALL get_averages(averages%avekin,md_ener%ekin,delta_t,error)
       CALL get_averages(averages%avetemp,md_ener%temp_part,delta_t,error)
       CALL get_averages(averages%avekin_qm,md_ener%ekin_qm,delta_t,error)
       CALL get_averages(averages%avetemp_qm,md_ener%temp_qm,delta_t,error)
       CALL get_averages(averages%avevol,cell%deth,delta_t,error)
       CALL get_averages(averages%aveca,abc(1),delta_t,error)
       CALL get_averages(averages%avecb,abc(2),delta_t,error)
       CALL get_averages(averages%avecc,abc(3),delta_t,error)
       CALL get_averages(averages%aveal,cell_angle(3),delta_t,error) 
       CALL get_averages(averages%avebe,cell_angle(2),delta_t,error)
       CALL get_averages(averages%avega,cell_angle(1),delta_t,error)
       CALL get_averages(averages%econs,md_ener%delta_cons,delta_t,error)
       ! Virial
       IF (ASSOCIATED(averages%virial).AND.virial%pv_availability) THEN
          CALL get_averages(averages%virial%pv_total,virial%pv_total,delta_t,error)
          CALL get_averages(averages%virial%pv_virial,virial%pv_virial,delta_t,error)
          CALL get_averages(averages%virial%pv_kinetic,virial%pv_kinetic,delta_t,error)
          CALL get_averages(averages%virial%pv_constraint,virial%pv_constraint,delta_t,error)
          CALL get_averages(averages%virial%pv_xc,virial%pv_xc,delta_t,error)
       END IF
       ! Colvars
       IF (nint>0) THEN
          ALLOCATE(cvalues(nint),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(Mmatrix(nint*nint),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          CALL get_clv_force(force_env, nsize_xyz=nat*3, nsize_int=nint, cvalues=cvalues,&
                             Mmatrix=Mmatrix, error=error)
          CALL get_averages(averages%avecolvar, cvalues, delta_t,error)
          CALL get_averages(averages%aveMmatrix, Mmatrix, delta_t,error)
          DEALLOCATE(cvalues,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(Mmatrix,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END SELECT
  END SUBROUTINE compute_averages

! *****************************************************************************
!> \brief computes the averages - low level for REAL
!> \author Teodoro Laino [tlaino] - 03.2008 - University of Zurich
! *****************************************************************************
  SUBROUTINE get_averages_rs(avg, add, delta_t, error)
    REAL(KIND=dp), INTENT(INOUT)             :: avg
    REAL(KIND=dp), INTENT(IN)                :: add
    INTEGER, INTENT(IN)                      :: delta_t
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_averages_rs', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    avg = (avg*REAL(delta_t-1,dp) + add)/REAL(delta_t,dp)
  END SUBROUTINE get_averages_rs

! *****************************************************************************
!> \brief computes the averages - low level for REAL vector
!> \author Teodoro Laino [tlaino] - 10.2008 - University of Zurich
! *****************************************************************************
  SUBROUTINE get_averages_rv(avg, add, delta_t, error)
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: avg
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: add
    INTEGER, INTENT(IN)                      :: delta_t
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_averages_rv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: check, failure

    failure = .FALSE.
    check = SIZE(avg)==SIZE(add)
    CPPostcondition(check,cp_failure_level,routineP,error,failure)
    DO i = 1, SIZE(avg)
       avg(i) = (avg(i)*REAL(delta_t-1,dp) + add(i))/REAL(delta_t,dp)
    END DO
  END SUBROUTINE get_averages_rv

! *****************************************************************************
!> \brief computes the averages - low level for REAL matrix
!> \author Teodoro Laino [tlaino] - 10.2008 - University of Zurich
! *****************************************************************************
  SUBROUTINE get_averages_rm(avg, add, delta_t, error)
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: avg
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: add
    INTEGER, INTENT(IN)                      :: delta_t
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_averages_rm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j
    LOGICAL                                  :: check, failure

    failure = .FALSE.
    check = SIZE(avg,1)==SIZE(add,1)
    CPPostcondition(check,cp_failure_level,routineP,error,failure)
    check = SIZE(avg,2)==SIZE(add,2)
    CPPostcondition(check,cp_failure_level,routineP,error,failure)
    DO i = 1, SIZE(avg,2)
       DO j = 1, SIZE(avg,1)
          avg(j,i) = (avg(j,i)*REAL(delta_t-1,dp) + add(j,i))/REAL(delta_t,dp)
       END DO
    END DO
  END SUBROUTINE get_averages_rm

END MODULE averages_types
