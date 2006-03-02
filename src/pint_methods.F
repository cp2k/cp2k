!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
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
                                             pint_env_create,&
                                             pint_env_type,&
                                             pint_f2uf,&
                                             pint_run,&
                                             pint_start,&
                                             pint_u2x,&
                                             pint_x2u
  USE replica_methods,                 ONLY: rep_env_calc_e_f,&
                                             rep_env_create
  USE replica_types,                   ONLY: rep_env_release,&
                                             replica_env_type
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


!!****f* pint_methods/rep_env_test *
!!
!!   NAME
!!     rep_env_test
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE rep_env_test(para_env,input,error)
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: input
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rep_env_test', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: err
    REAL(kind=dp), DIMENSION(:), POINTER     :: e1, e2, e3
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: f1, f2, f3
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(replica_env_type), POINTER          :: rep_env

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(para_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(input),cp_failure_level,routineP,error,failure)
  CPPrecondition(para_env%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(input%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     logger => cp_error_get_logger(error)
     NULLIFY(rep_env)
     CALL rep_env_create(rep_env, para_env=para_env, input=input,&
          nrep=1,prep=1,error=error)
     IF (ASSOCIATED(rep_env)) THEN
        CALL rep_env_calc_e_f(rep_env,error=error)
        ALLOCATE(f1(rep_env%ndim,rep_env%nrep),e1(rep_env%nrep),stat=stat)
        f1=rep_env%f(:rep_env%ndim,:)
        e1=rep_env%f(rep_env%ndim+1,:)
     END IF
     CALL rep_env_release(rep_env,error=error)

     CALL rep_env_create(rep_env, para_env=para_env, input=input,&
          nrep=1,prep=2,error=error)
     IF (ASSOCIATED(rep_env)) THEN
        CALL rep_env_calc_e_f(rep_env,error=error)
        ALLOCATE(f2(rep_env%ndim,rep_env%nrep),e2(rep_env%nrep),stat=stat)
        f2=rep_env%f(:rep_env%ndim,:)
        e2=rep_env%f(rep_env%ndim+1,:)
     END IF
     CALL rep_env_release(rep_env,error=error)

     CALL rep_env_create(rep_env, para_env=para_env, input=input,&
          nrep=2,prep=1,error=error)
     IF (ASSOCIATED(rep_env)) THEN
        rep_env%r(:,2)=rep_env%r(:,1)+0.5_dp
        CALL rep_env_calc_e_f(rep_env,error=error)
        ALLOCATE(f3(rep_env%ndim,rep_env%nrep),e3(rep_env%nrep),stat=stat)
        f3=rep_env%f(:rep_env%ndim,:)
        e3=rep_env%f(rep_env%ndim+1,:)
     END IF
     CALL rep_env_release(rep_env,error=error)

     IF (ASSOCIATED(f1).AND.ASSOCIATED(f2)) THEN
        err=0._dp
        DO i=1,SIZE(f1,1)
           err=MAX(err,ABS(f1(i,1)-f2(i,1)))
        END DO
        CALL cp_log(logger,cp_note_level+1,routineP,"diff_f12="//cp_to_string(err),&
             local=.FALSE.)
        CALL cp_log(logger,cp_note_level+1,routineP,"diff_e12="//&
             cp_to_string(ABS(e1(1)-e2(1)))//cp_to_string(e1(1))//&
             cp_to_string(e2(1)), local=.FALSE.)
     END IF
     IF (ASSOCIATED(f1).AND.ASSOCIATED(f3)) THEN
        err=0._dp
        DO i=1,SIZE(f1,1)
           err=MAX(err,ABS(f1(i,1)-f3(i,1)))
        END DO
        CALL cp_log(logger,cp_note_level+1,routineP,"diff_f13="//cp_to_string(err),&
             local=.FALSE.)
        CALL cp_log(logger,cp_note_level+1,routineP,"diff_e13="//&
             cp_to_string(ABS(e1(1)-e3(1)))//cp_to_string(e1(1))//&
             cp_to_string(e3(1)),&
             local=.FALSE.)
     END IF
     IF (ASSOCIATED(f3)) THEN
        err=0._dp
        DO i=1,SIZE(f1,1)
           err=MAX(err,ABS(f3(i,1)-f3(i,2)))
        END DO
        CALL cp_log(logger,cp_note_level+1,routineP,"diff_f33="//cp_to_string(err),&
             local=.FALSE.)
        CALL cp_log(logger,cp_note_level+1,routineP,"diff_e33="//&
             cp_to_string(ABS(e3(1)-e3(2)))//cp_to_string(e3(1))//&
             cp_to_string(e3(2)), local=.FALSE.)
        DEALLOCATE(f3,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
     IF (ASSOCIATED(f1)) THEN
        DEALLOCATE(f1,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
     IF (ASSOCIATED(f2)) THEN
        DEALLOCATE(f2,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
END SUBROUTINE rep_env_test
!****************************************************************************

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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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

        e_h=pint_calc_uf_h(pint_env,error=error)
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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'do_pint_run', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(pint_env_type), POINTER             :: pint_env

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(para_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(input),cp_failure_level,routineP,error,failure)
  CPPrecondition(para_env%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(input%ref_count>0,cp_failure_level,routineP,error,failure)
  NULLIFY(pint_env)
  IF (.NOT. failure) THEN
     CALL section_vals_val_set(input,"GLOBAL%RUN_TYPE",i_val=mol_dyn_run,&
          error=error)
     CALL section_vals_val_set(input,"GLOBAL%PROGRAM",i_val=do_qs,&
          error=error)
     CALL rep_env_test(para_env,input,error=error)
     CALL pint_env_test(para_env,input,error=error)

     CALL pint_env_create(pint_env,input,para_env,error=error)
     CALL pint_start(pint_env,error=error)
     CALL pint_run(pint_env,error=error)
  END IF
END SUBROUTINE do_pint_run
!***************************************************************************

END MODULE pint_methods
