!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C ) 2002 - 2003  CP2K developers group                         !
!-----------------------------------------------------------------------------!
!!****** cp2k/kg_rho_methods [1.0] *
!!
!!   NAME
!!     kg_rho_methods
!!
!!   FUNCTION
!!     rho_methods    
!!
!!   AUTHOR
!!     gloria,30.09.2002, previously in kg_force 
!!
!!   MODIFICATION HISTORY
!!     JGH (22-Feb-03) PW grid options added
!!     gt 16-nov-03 moved initialization in this new module
!!     13-10-04 only density: pw stuff moved in kg_rspw_methods
!!
!!   SOURCE
!******************************************************************************

MODULE kg_rho_methods
  USE coefficient_types,               ONLY: coeff_deallocate
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kg_environment_types,            ONLY: get_kg_env,&
                                             kg_environment_type
  USE kg_rho_types,                    ONLY: kg_rho_type
  USE kg_rspw_types,                   ONLY: kg_rspw_type,&
                                             kg_rspw_get
  USE kinds,                           ONLY: dp
  USE pw_pool_types,                   ONLY: pw_pool_init_coeff,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='kg_rho_methods'
  INTEGER, PRIVATE, SAVE :: last_rho_id_nr=0

  PUBLIC :: kg_rho_create

!-----------------------------------------------------------------------------!

CONTAINS

!!*************************************************************************
SUBROUTINE kg_rho_create(kg_rho,kg_env,error)
    TYPE(kg_rho_type), POINTER               :: kg_rho
    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(kg_rspw_type), POINTER              :: kg_rspw
    TYPE(dft_control_type), POINTER          :: dft_control
    CHARACTER(len=*), PARAMETER :: routineN = 'kg_rho_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle,stat
    LOGICAL                                  :: failure,kgpol
!-----------------------------------------------------------------------------!

  CALL timeset(routineN,'I',' ',handle)

  NULLIFY(kg_rspw,auxbas_pw_pool,dft_control)
  failure=.FALSE.
  kgpol=.FALSE.
  CALL get_kg_env(kg_env,rspw=kg_rspw,dft_control=dft_control,error=error)
  kgpol = dft_control % qs_control % polarization
  CALL kg_rspw_get(kg_rspw,auxbas_pw_pool=auxbas_pw_pool,error=error)
 
  IF (.NOT.ASSOCIATED(kg_rho)) THEN
    ALLOCATE(kg_rho, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
      kg_rho%ref_count=1
      last_rho_id_nr=last_rho_id_nr+1
      kg_rho%id_nr=last_rho_id_nr
      NULLIFY(kg_rho%rho_r, kg_rho%rho_g, kg_rho%rhop_r, kg_rho%rhop_g,&
          kg_rho%rho_core)
    END IF
  END IF

!rho_r
  IF(.NOT.failure) THEN
    IF(ASSOCIATED(kg_rho%rho_r)) THEN
      CALL coeff_deallocate(kg_rho%rho_r)
      DEALLOCATE(kg_rho%rho_r,stat=stat)
      CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
    ALLOCATE(kg_rho%rho_r,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(.NOT.failure) THEN
    CALL pw_pool_init_coeff(auxbas_pw_pool,kg_rho%rho_r,&
         use_data=REALDATA3D,in_space=REALSPACE,error=error)
  END IF

!rho_g
  IF(.NOT.failure) THEN
    IF(ASSOCIATED(kg_rho%rho_g)) THEN
      CALL coeff_deallocate(kg_rho%rho_g)
      DEALLOCATE(kg_rho%rho_g,stat=stat)
      CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
    ALLOCATE(kg_rho%rho_g,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(.NOT.failure) THEN
    CALL pw_pool_init_coeff(auxbas_pw_pool,kg_rho%rho_g,&
         use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,error=error)
  END IF

!rho_core
  IF(.NOT.failure) THEN
    IF(ASSOCIATED(kg_rho%rho_core)) THEN
      CALL coeff_deallocate(kg_rho%rho_core)
      DEALLOCATE(kg_rho%rho_core,stat=stat)
      CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
    ALLOCATE(kg_rho%rho_core,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(.NOT.failure) THEN
    CALL pw_pool_init_coeff(auxbas_pw_pool,kg_rho%rho_core,&
         use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,error=error)
  END IF

  IF (kgpol) THEN
!rhop_r
    IF(.NOT.failure) THEN
      IF(ASSOCIATED(kg_rho%rhop_r)) THEN
        CALL coeff_deallocate(kg_rho%rhop_r)
        DEALLOCATE(kg_rho%rhop_r,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF
      ALLOCATE(kg_rho%rhop_r,stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(.NOT.failure) THEN
      CALL pw_pool_init_coeff(auxbas_pw_pool,kg_rho%rhop_r,&
           use_data=REALDATA3D,in_space=REALSPACE,error=error)
    END IF
!rhop_g
    IF(.NOT.failure) THEN
      IF(ASSOCIATED(kg_rho%rhop_g)) THEN
        CALL coeff_deallocate(kg_rho%rhop_g)
        DEALLOCATE(kg_rho%rhop_g,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF
      ALLOCATE(kg_rho%rhop_g,stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(.NOT.failure) THEN
      CALL pw_pool_init_coeff(auxbas_pw_pool,kg_rho%rhop_g,&
           use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,error=error)
    END IF

  END IF

  CALL timestop(0.0_dp,handle)

END SUBROUTINE kg_rho_create
!!*** **********************************************************************

END MODULE kg_rho_methods

!******************************************************************************
