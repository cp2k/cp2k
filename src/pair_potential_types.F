!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****** cp2k/pair_potential_types [1.0] *
!!
!!   NAME
!!     pair_potential_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE pair_potential_types

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
  USE kinds,                           ONLY: dp,&
                                             default_string_length
  USE splines,                         ONLY: spline_data_type, &
                                             spline_data_release
  USE termination,                     ONLY: stop_memory

  IMPLICIT NONE

  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='pair_potential_types'

  PRIVATE
  PUBLIC :: pair_potential_type, pot_parameters_release
  PUBLIC :: pot_parameters_create
  PUBLIC :: lj_pot_type, williams_pot_type, goodwin_pot_type, ipbv_pot_type
  PUBLIC :: eam_pot_type, clean_pot_parameter

  PUBLIC :: reallocate_potparm

  TYPE ipbv_pot_type
    REAL (KIND=dp), DIMENSION (2:15) :: a
    REAL (KIND=dp) :: rcore
    REAL (KIND=dp) :: m
    REAL (KIND=dp) :: b
  END TYPE ipbv_pot_type

  TYPE lj_pot_type
    REAL (KIND=dp) :: epsilon
    REAL (KIND=dp) :: sigma6
    REAL (KIND=dp) :: sigma12
  END TYPE lj_pot_type

  TYPE williams_pot_type
    REAL (KIND=dp) :: a
    REAL (KIND=dp) :: b
    REAL (KIND=dp) :: c
  END TYPE williams_pot_type

  TYPE goodwin_pot_type
    REAL (KIND=dp) :: vr0
    REAL (KIND=dp) :: m, mc
    REAL (KIND=dp) :: d, dc
  END TYPE goodwin_pot_type

  TYPE eam_pot_type
    CHARACTER ( len = default_string_length ) :: eam_file_name
    REAL ( KIND = dp ) drar, drhoar, acutal
    REAL ( KIND = dp ), DIMENSION ( 5000 ) :: rho, phi, frho, rhoval, rval
    REAL ( KIND = dp ), DIMENSION ( 5000 ) :: rhop, phip, frhop
  END TYPE eam_pot_type

  TYPE pair_potential_type
    CHARACTER ( LEN = default_string_length ) :: type
    CHARACTER ( LEN = default_string_length ) :: at1
    CHARACTER ( LEN = default_string_length ) :: at2
    TYPE ( ipbv_pot_type ) :: ipbv
    TYPE ( lj_pot_type ) :: lj
    TYPE ( williams_pot_type ) :: willis
    TYPE ( goodwin_pot_type ) :: goodwin
    TYPE ( eam_pot_type ) :: eam
    REAL (KIND=dp) :: rcutsq
    REAL (KIND=dp) :: energy_cutoff
    REAL (KIND=dp) :: e_cutoff_coul
    TYPE ( spline_data_type ), POINTER :: pair_spline_data
! for EAM type potentials
    TYPE ( spline_data_type ), POINTER :: rho_spline_data
  END TYPE pair_potential_type

!!*****
!******************************************************************************

CONTAINS
!******************************************************************************
!!****
!!****s* pair_potential_types/pot_parameters_release [1.0] *
!!
!!   NAME
!!     pot_parameters_release
!!
!!   FUNCTION
!!     Release Data-structure that constains potential parameters
!!
!!   AUTHOR
!!     unknown
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE pot_parameters_release ( potparm, error )
  TYPE(pair_potential_type), POINTER :: potparm ( :, : )
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  CHARACTER(len=*), PARAMETER :: routineN='pot_parameters_release',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat, i, j
  
  IF ( ASSOCIATED ( potparm ) ) THEN
     DO i = 1, SIZE ( potparm, 1 )
        DO j = 1, SIZE ( potparm, 2 )
           CALL spline_data_release ( potparm ( i, j ) % pair_spline_data )
           CALL spline_data_release ( potparm ( i, j ) % rho_spline_data )
        END DO
     END DO
     DEALLOCATE ( potparm, stat = stat )
     CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)        
  END IF
  NULLIFY ( potparm )
END SUBROUTINE pot_parameters_release
!******************************************************************************
!!****
!!****s* pair_potential_types/pot_parameters_create [1.0] *
!!
!!   NAME
!!     pot_parameters_create
!!
!!   FUNCTION
!!     Data-structure that constains potential parameters
!!
!!   AUTHOR
!!     unknown
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE pot_parameters_create ( potparm, nkinds, error )
  TYPE(pair_potential_type), POINTER :: potparm ( :, : )
  INTEGER, INTENT ( IN ) :: nkinds
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='pot_parameters_create',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat, i, j

  failure=.FALSE.
  
  ALLOCATE ( potparm ( nkinds, nkinds ), stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
    potparm ( :, : ) % type = 'NULL'
    potparm ( :, : ) % rcutsq = 0._dp
    potparm ( :, : ) % energy_cutoff = 0._dp
    potparm ( :, : ) % e_cutoff_coul = 0._dp
    potparm ( :, : ) % ipbv % a ( 2 ) = 0.0_dp
    potparm ( :, : ) % ipbv % a ( 3 ) = 0.0_dp
    potparm ( :, : ) % ipbv % a ( 4 ) = 0.0_dp
    potparm ( :, : ) % ipbv % a ( 5 ) = 0.0_dp
    potparm ( :, : ) % ipbv % a ( 6 ) = 0.0_dp
    potparm ( :, : ) % ipbv % a ( 7 ) = 0.0_dp
    potparm ( :, : ) % ipbv % a ( 8 ) = 0.0_dp
    potparm ( :, : ) % ipbv % a ( 9 ) = 0.0_dp
    potparm ( :, : ) % ipbv % a ( 10 ) = 0.0_dp
    potparm ( :, : ) % ipbv % a ( 11 ) = 0.0_dp
    potparm ( :, : ) % ipbv % a ( 12 ) = 0.0_dp
    potparm ( :, : ) % ipbv % a ( 13 ) = 0.0_dp
    potparm ( :, : ) % ipbv % a ( 14 ) = 0.0_dp
    potparm ( :, : ) % ipbv % a ( 15 ) = 0.0_dp
    potparm ( :, : ) % ipbv % rcore = 0.0_dp
    potparm ( :, : ) % ipbv % b = 0.0_dp
    potparm ( :, : ) % ipbv % m = 0.0_dp
    potparm ( :, : ) % lj % epsilon = 0.0_dp
    potparm ( :, : ) % lj % sigma6 = 0.0_dp
    potparm ( :, : ) % lj % sigma12 = 0.0_dp
    potparm ( :, : ) % willis % a = 0.0_dp
    potparm ( :, : ) % willis % b = 0.0_dp
    potparm ( :, : ) % willis % c = 0.0_dp
    potparm ( :, : ) % goodwin % vr0 = 0.0_dp
    potparm ( :, : ) % goodwin % d = 0.0_dp
    potparm ( :, : ) % goodwin % dc = 0.0_dp
    potparm ( :, : ) % goodwin % m = 0.0_dp
    potparm ( :, : ) % goodwin % mc = 0.0_dp
    potparm ( :, : ) % eam % eam_file_name = 'NULL'
    DO i = 1, nkinds
     DO j = 1, nkinds
      NULLIFY ( potparm ( i, j ) % pair_spline_data )
      NULLIFY ( potparm ( i, j ) % rho_spline_data )
     END DO
    END DO
    potparm ( :, : ) % eam % drar = 0._dp
    potparm ( :, : ) % eam % drhoar = 0._dp
    potparm ( :, : ) % eam % acutal = 0._dp
    DO i = 1, 5000
      potparm ( :, : ) % eam % rho ( i ) = 0._dp
      potparm ( :, : ) % eam % phi ( i ) = 0._dp
      potparm ( :, : ) % eam % frho ( i ) = 0._dp
      potparm ( :, : ) % eam % rhoval ( i ) = 0._dp
      potparm ( :, : ) % eam % rval ( i ) = 0._dp
      potparm ( :, : ) % eam % rhop ( i ) = 0._dp
      potparm ( :, : ) % eam % phip ( i ) = 0._dp
      potparm ( :, : ) % eam % frhop ( i ) = 0._dp
    END DO
  END IF
END SUBROUTINE pot_parameters_create

!******************************************************************************
!!****
!!****s* pair_potential_types/clean_pot_parameter [1.0] *
!!
!!   NAME
!!     clean_pot_parameter
!!
!!   FUNCTION
!!     Cleans the potential parameter type
!!
!!   AUTHOR
!!     unknown
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE clean_pot_parameter(potparm)
  IMPLICIT NONE
  TYPE(pair_potential_type) :: potparm
  INTEGER :: i

  potparm % TYPE = 'NULL'
  potparm % rcutsq = 0._dp
  potparm % energy_cutoff = 0._dp
  potparm % e_cutoff_coul = 0._dp
  potparm % ipbv % a ( 2 ) = 0.0_dp
  potparm % ipbv % a ( 3 ) = 0.0_dp
  potparm % ipbv % a ( 4 ) = 0.0_dp
  potparm % ipbv % a ( 5 ) = 0.0_dp
  potparm % ipbv % a ( 6 ) = 0.0_dp
  potparm % ipbv % a ( 7 ) = 0.0_dp
  potparm % ipbv % a ( 8 ) = 0.0_dp
  potparm % ipbv % a ( 9 ) = 0.0_dp
  potparm % ipbv % a ( 10 ) = 0.0_dp
  potparm % ipbv % a ( 11 ) = 0.0_dp
  potparm % ipbv % a ( 12 ) = 0.0_dp
  potparm % ipbv % a ( 13 ) = 0.0_dp
  potparm % ipbv % a ( 14 ) = 0.0_dp
  potparm % ipbv % a ( 15 ) = 0.0_dp
  potparm % ipbv % rcore = 0.0_dp
  potparm % ipbv % b = 0.0_dp
  potparm % ipbv % m = 0.0_dp
  potparm % lj % epsilon = 0.0_dp
  potparm % lj % sigma6 = 0.0_dp
  potparm % lj % sigma12 = 0.0_dp
  potparm % willis % a = 0.0_dp
  potparm % willis % b = 0.0_dp
  potparm % willis % c = 0.0_dp
  potparm % goodwin % vr0 = 0.0_dp
  potparm % goodwin % d = 0.0_dp
  potparm % goodwin % dc = 0.0_dp
  potparm % goodwin % m = 0.0_dp
  potparm % goodwin % mc = 0.0_dp
  potparm % eam % eam_file_name = 'NULL'
  IF (ASSOCIATED( potparm % pair_spline_data ))&
       CALL spline_data_release ( potparm % pair_spline_data )
  IF (ASSOCIATED( potparm % rho_spline_data )) &
       CALL spline_data_release ( potparm % rho_spline_data )
  potparm % eam % drar = 0._dp
  potparm % eam % drhoar = 0._dp
  potparm % eam % acutal = 0._dp
  DO i = 1, 5000
     potparm % eam % rho ( i ) = 0._dp
     potparm % eam % phi ( i ) = 0._dp
     potparm % eam % frho ( i ) = 0._dp
     potparm % eam % rhoval ( i ) = 0._dp
     potparm % eam % rval ( i ) = 0._dp
     potparm % eam % rhop ( i ) = 0._dp
     potparm % eam % phip ( i ) = 0._dp
     potparm % eam % frhop ( i ) = 0._dp
  END DO

END SUBROUTINE clean_pot_parameter

!******************************************************************************
  SUBROUTINE reallocate_potparm(p,lb1_new,ub1_new)
    TYPE (pair_potential_type), POINTER       :: p(:)
    INTEGER, INTENT(IN)                :: lb1_new, ub1_new

    CHARACTER(LEN=*), PARAMETER        :: routine_name = "reallocate_potparm"

    INTEGER                            :: istat, lb1, lb1_old, &
                                          new_size, old_size, ub1, &
                                          ub1_old, i
    TYPE (pair_potential_type), ALLOCATABLE :: work(:)

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(p)) THEN
      lb1_old = LBOUND(p,1)
      ub1_old = UBOUND(p,1)
      ALLOCATE (work(lb1_old:ub1_old),STAT=istat)
      IF (istat /= 0) THEN
      END IF
      work(:) = p(:)
      DEALLOCATE (p,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,moduleN,__LINE__,"p")
    END IF

    ALLOCATE (p(lb1_new:ub1_new),STAT=istat)
    IF (istat /= 0) THEN
      new_size = (ub1_new - lb1_new + 1)
      CALL stop_memory(routine_name,moduleN,__LINE__,"p",new_size)
    END IF

    p ( : ) % type = 'NULL'
    p ( : ) % ipbv % a ( 2 ) = 0.0_dp
    p ( : ) % ipbv % a ( 3 ) = 0.0_dp
    p ( : ) % ipbv % a ( 4 ) = 0.0_dp
    p ( : ) % ipbv % a ( 5 ) = 0.0_dp
    p ( : ) % ipbv % a ( 6 ) = 0.0_dp
    p ( : ) % ipbv % a ( 7 ) = 0.0_dp
    p ( : ) % ipbv % a ( 8 ) = 0.0_dp
    p ( : ) % ipbv % a ( 9 ) = 0.0_dp
    p ( : ) % ipbv % a ( 10 ) = 0.0_dp
    p ( : ) % ipbv % a ( 11 ) = 0.0_dp
    p ( : ) % ipbv % a ( 12 ) = 0.0_dp
    p ( : ) % ipbv % a ( 13 ) = 0.0_dp
    p ( : ) % ipbv % a ( 14 ) = 0.0_dp
    p ( : ) % ipbv % a ( 15 ) = 0.0_dp
    p ( : ) % ipbv % rcore = 0.0_dp
    p ( : ) % ipbv % b = 0.0_dp
    p ( : ) % ipbv % m = 0.0_dp
    p ( : ) % lj % epsilon = 0.0_dp
    p ( : ) % lj % sigma6 = 0.0_dp
    p ( : ) % lj % sigma12 = 0.0_dp
    p ( : ) % willis % a = 0.0_dp
    p ( : ) % willis % b = 0.0_dp
    p ( : ) % willis % c = 0.0_dp
    p ( : ) % goodwin % vr0 = 0.0_dp
    p ( : ) % goodwin % d = 0.0_dp
    p ( : ) % goodwin % dc = 0.0_dp
    p ( : ) % goodwin % m = 0.0_dp
    p ( : ) % goodwin % mc = 0.0_dp
    p ( : ) % eam % eam_file_name = 'NULL'
    DO i = 1, SIZE ( p )
      NULLIFY ( p ( i ) % pair_spline_data )
      NULLIFY ( p ( i ) % rho_spline_data )
    END DO

    IF (ALLOCATED(work)) THEN
      lb1 = MAX(lb1_new,lb1_old)
      ub1 = MIN(ub1_new,ub1_old)
      p(lb1:ub1) = work(lb1:ub1)
      DEALLOCATE (work,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine_name,moduleN,__LINE__,"work")
      END IF
    END IF

  END SUBROUTINE reallocate_potparm

END MODULE pair_potential_types
