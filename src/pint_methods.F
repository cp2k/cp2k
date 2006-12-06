!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2006  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****h* cp2k/pint_methods *
!!
!!   NAME
!!     pint_methods
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2005 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE pint_methods
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_constants,                 ONLY: do_qs,&
                                             mol_dyn_run
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_set
  USE kinds,                           ONLY: dp
  USE pint_types,                      ONLY: pint_calc_uf_h,&
                                             pint_do_run,&
                                             pint_env_create,&
                                             pint_env_type,&
                                             pint_f2uf,&
                                             pint_release,&
                                             pint_start,&
                                             pint_u2x,&
                                             pint_x2u
  USE replica_methods,                 ONLY: rep_env_test
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pint_methods'

  PUBLIC :: do_pint_run
!!***
!****************************************************************************
CONTAINS

!!****f* pint_methods/pint_env_test *
!!
!!   NAME
!!     pint_env_test
!!
!!   FUNCTION
!!     tests the path integral methods
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - para_env: parallel environment
!!     - input: the input to test
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE pint_env_test(para_env,input,error)
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: input
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_env_test', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ib, idim, stat
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: c, e_h, err
    REAL(kind=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: x1
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(pint_env_type), POINTER             :: pint_env

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(para_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(input),cp_failure_level,routineP,error,failure)
  CPPrecondition(para_env%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(input%ref_count>0,cp_failure_level,routineP,error,failure)
  NULLIFY(pint_env)
  IF (.NOT. failure) THEN
     logger => cp_error_get_logger(error)
     CALL pint_env_create(pint_env,input,para_env,error=error)
     IF (ASSOCIATED(pint_env)) THEN
        ALLOCATE(x1(pint_env%ndim,pint_env%p),stat=stat)
        CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
        x1=pint_env%x
        CALL pint_x2u(pint_env,error=error)
        pint_env%x=0._dp
        CALL pint_u2x(pint_env,error=error)
        err=0._dp
        DO i=1,pint_env%ndim
           err=MAX(err,ABS(x1(1,i)-pint_env%x(1,i)))
        END DO
        CALL cp_log(logger,cp_note_level+1,routineP,"diff_r1="//cp_to_string(err),&
             local=.FALSE.)

        CALL pint_calc_uf_h(pint_env,e_h=e_h,error=error)
        c=-pint_env%staging_env%w_p**2
        pint_env%f=0._dp
        DO idim=1,pint_env%ndim
           DO ib=1,pint_env%p
              pint_env%f(ib,idim)=pint_env%f(ib,idim)+&
                   c*(2._dp*pint_env%x(ib,idim)&
                   -pint_env%x(MODULO(ib-2,pint_env%p)+1,idim)&
                   -pint_env%x(MODULO(ib,pint_env%p)+1,idim))
           END DO
        END DO
        CALL pint_f2uf(pint_env,error=error)
        err=0._dp
        DO idim=1,pint_env%ndim
           DO ib=1,pint_env%p
              err=MAX(err,ABS(pint_env%uf(ib,idim)-pint_env%uf_h(ib,idim)))
           END DO
        END DO
        CALL cp_log(logger,cp_note_level+1,routineP,"diff_f_h="//cp_to_string(err),&
             local=.FALSE.)
     END IF
  END IF
END SUBROUTINE pint_env_test
!***************************************************************************

!!****f* pint_methods/do_pint_run [1.0] *
!!
!!   NAME
!!     do_pint_run
!!
!!   FUNCTION
!!     perform a path integral simulation
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - para_env: parallel environment
!!     - input: the input to test
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE do_pint_run(para_env,input,error)
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: input
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error

    TYPE(pint_env_type), POINTER             :: pint_env

    CHARACTER(len=*), PARAMETER :: routineN = 'do_pint_run', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(para_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(input),cp_failure_level,routineP,error,failure)
  CPPrecondition(para_env%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(input%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL section_vals_val_set(input,"GLOBAL%RUN_TYPE",i_val=mol_dyn_run,&
          error=error)
     CALL section_vals_val_set(input,"FORCE_EVAL%METHOD",i_val=do_qs,&
          error=error)
!x     CALL rep_env_test(para_env,input,error=error)
!x     CALL pint_env_test(para_env,input,error=error)
!x

     NULLIFY(pint_env)
     CALL pint_env_create(pint_env,input,para_env,error=error)
     CALL pint_start(pint_env,error=error)
     CALL pint_do_run(pint_env,error=error)
     CALL pint_release(pint_env,error=error)
  END IF
END SUBROUTINE do_pint_run
!***************************************************************************

END MODULE pint_methods
