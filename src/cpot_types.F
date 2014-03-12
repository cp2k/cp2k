!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief makes a run that valuates the chemical potential difference between two
!>      different two box situations
!> \par History
!>      06.2006 created
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE cpot_types
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_release,&
                                             section_vals_retain,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             dp
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
#include "f77_int.h"

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cpot_types'
  INTEGER, PRIVATE, SAVE :: last_cpot_id=0

  PUBLIC :: cpot_type, cpot_release, cpot_retain, cpot_calc

! *****************************************************************************
!> \brief type for the chemical potential difference
!> \author fawzi
! *****************************************************************************
TYPE cpot_type
   INTEGER :: id_nr, ref_count, s_1_2_id, s_2_1_id, s_2_2_id
   LOGICAL :: active
   TYPE(section_vals_type), POINTER :: input
   REAL(kind=dp) :: sum_cpot, sum_cpotsq,&
        e_s_1_1,e_s_1_2,e_s_2_1,e_s_2_2, beta
   INTEGER :: npot
END TYPE cpot_type

CONTAINS

! *****************************************************************************
!> \brief create a cpot_env
!> \param cpot_env the cpot_env to create
!> \param cpot_section the cpot_input_section
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE cpot_create(cpot_env,cpot_section,para_env,error)
    TYPE(cpot_type), POINTER                 :: cpot_env
    TYPE(section_vals_type), POINTER         :: cpot_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cpot_create', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_path_length)       :: input_path
    INTEGER                                  :: handle, ierr, stat
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: temp

    failure=.FALSE.

    CALL timeset(routineN,handle)
    CPPrecondition(.NOT.ASSOCIATED(cpot_env),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       ALLOCATE(cpot_env,stat=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
       last_cpot_id=last_cpot_id+1
       cpot_env%id_nr=last_cpot_id
       cpot_env%ref_count=1
       cpot_env%input => cpot_section
       CALL section_vals_retain(cpot_env%input,error=error)
       cpot_env%s_1_2_id=-HUGE(0)
       cpot_env%s_2_1_id=-HUGE(0)
       cpot_env%s_2_2_id=-HUGE(0)

       CALL section_vals_val_get(cpot_env%input,"_DEFAULT_KEYWORD_",&
            l_val=cpot_env%active,error=error)
       CALL section_vals_val_get(cpot_env%input,"temp",&
            r_val=temp,error=error)
       cpot_env%beta=1/temp
       IF (cpot_env%active) THEN
          CALL section_vals_val_get(cpot_env%input,"state_1_2_input",&
               c_val=input_path,error=error)
          CALL cp_create_fenv_comm(cpot_env%s_1_2_id,input_path,&
               input_path//".out",para_env%group,ierr)
          CALL section_vals_val_get(cpot_env%input,"state_2_1_input",&
               c_val=input_path,error=error)
          CALL cp_create_fenv_comm(cpot_env%s_2_1_id,input_path,&
               input_path//".out",para_env%group,ierr)
          CALL section_vals_val_get(cpot_env%input,"state_2_2_input",&
               c_val=input_path,error=error)
          CALL cp_create_fenv_comm(cpot_env%s_2_2_id,input_path,&
               input_path//".out",para_env%group,ierr)
       END IF
    END IF
    CALL timestop(handle)
  END SUBROUTINE cpot_create

! *****************************************************************************
!> \brief retains the cpot
!> \param cpot_env the object to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
SUBROUTINE cpot_retain(cpot_env,error)
    TYPE(cpot_type), POINTER                 :: cpot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cpot_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  CPPrecondition(ASSOCIATED(cpot_env),cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     CPPrecondition(cpot_env%ref_count>0,cp_failure_level,routineP,error,failure)
     cpot_env%ref_count=cpot_env%ref_count+1
  END IF
END SUBROUTINE cpot_retain

! *****************************************************************************
!> \brief releases the cpot
!> \param cpot_env the object to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
SUBROUTINE cpot_release(cpot_env,error)
    TYPE(cpot_type), POINTER                 :: cpot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cpot_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  IF (ASSOCIATED(cpot_env)) THEN
     CPPrecondition(cpot_env%ref_count>0,cp_failure_level,routineP,error,failure)
     cpot_env%ref_count=cpot_env%ref_count-1
     IF (cpot_env%ref_count==0) THEN
        CALL section_vals_release(cpot_env%input,error=error)
        DEALLOCATE(cpot_env,stat=stat)
        CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
     END IF
  END IF
  NULLIFY(cpot_env)
END SUBROUTINE cpot_release

! *****************************************************************************
!> \brief calculates the chemical potential
!> \param cpot_env the environment in which to calculate the chemical
!>        potential
!> \param e_s_ 1_1: potential energy of the first box of the first system
!>        (the driving system)
!> \param pos_s_ 1_1: positions of the first box
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
SUBROUTINE cpot_calc(cpot_env, e_s_1_1,pos_s_1_1,error)
    TYPE(cpot_type), POINTER                 :: cpot_env
    REAL(kind=dp), INTENT(in)                :: e_s_1_1
    REAL(kind=dp), DIMENSION(*), INTENT(in)  :: pos_s_1_1
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cpot_calc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, natom
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: dcpot, e1, e2

  failure=.FALSE.

  CALL timeset(routineN,handle)
  IF (.NOT.failure) THEN
     cpot_env%e_s_1_1=e_s_1_1
     CALL cp_get_natom(cpot_env%s_1_2_id,natom,ierr)
     CPAssert(ierr/=0,cp_failure_level,routineP,error,failure)
     CALL cp_set_pos(cpot_env%s_1_2_id,pos_s_1_1,3*natom,ierr)
     CPAssert(ierr/=0,cp_failure_level,routineP,error,failure)
     CALL cp_do_shake(cpot_env%s_1_2_id,1._dp,1.e-10_dp,ierr)
     CPAssert(ierr/=0,cp_failure_level,routineP,error,failure)
     CALL cp_calc_energy_force(cpot_env%s_1_2_id,0,ierr)
     CPAssert(ierr/=0,cp_failure_level,routineP,error,failure)
     CALL cp_get_energy(cpot_env%s_1_2_id,cpot_env%e_s_1_2,ierr)
     CPAssert(ierr/=0,cp_failure_level,routineP,error,failure)

     CALL cp_get_natom(cpot_env%s_2_1_id,natom,ierr)
     CPAssert(ierr/=0,cp_failure_level,routineP,error,failure)
     CALL cp_set_pos(cpot_env%s_2_1_id,pos_s_1_1,3*natom,ierr)
     CPAssert(ierr/=0,cp_failure_level,routineP,error,failure)
     CALL cp_do_shake(cpot_env%s_2_1_id,1._dp,1.e-10_dp,ierr)
     CPAssert(ierr/=0,cp_failure_level,routineP,error,failure)
     CALL cp_calc_energy_force(cpot_env%s_2_1_id,0,ierr)
     CPAssert(ierr/=0,cp_failure_level,routineP,error,failure)
     CALL cp_get_energy(cpot_env%s_2_1_id,cpot_env%e_s_2_1,ierr)
     CPAssert(ierr/=0,cp_failure_level,routineP,error,failure)

     CALL cp_get_natom(cpot_env%s_2_2_id,natom,ierr)
     CPAssert(ierr/=0,cp_failure_level,routineP,error,failure)
     CALL cp_set_pos(cpot_env%s_2_2_id,pos_s_1_1,3*natom,ierr)
     CPAssert(ierr/=0,cp_failure_level,routineP,error,failure)
     CALL cp_do_shake(cpot_env%s_2_2_id,1._dp,1.e-10_dp,ierr)
     CPAssert(ierr/=0,cp_failure_level,routineP,error,failure)
     CALL cp_calc_energy_force(cpot_env%s_2_2_id,0,ierr)
     CPAssert(ierr/=0,cp_failure_level,routineP,error,failure)
     CALL cp_get_energy(cpot_env%s_2_2_id,cpot_env%e_s_2_2,ierr)
     CPAssert(ierr/=0,cp_failure_level,routineP,error,failure)

     e1=cpot_env%e_s_1_1+cpot_env%e_s_2_1
     e2=cpot_env%e_s_2_1+cpot_env%e_s_2_2
     IF (e1>e2) THEN
        dcpot=EXP(-cpot_env%beta*(e1-e2))
     ELSE
        dcpot=-EXP(-cpot_env%beta*(e2-e1))
     END IF
     cpot_env%sum_cpot=cpot_env%sum_cpot+dcpot
     cpot_env%sum_cpotsq=cpot_env%sum_cpotsq+dcpot*dcpot
  END IF
  CALL timestop(handle)
END SUBROUTINE cpot_calc

END MODULE cpot_types
