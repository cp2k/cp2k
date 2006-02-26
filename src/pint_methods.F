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
  use kinds, only: dp
  use cp_para_types, only: cp_para_env_type
  use timings, only: timeset, timestop
  USE pint_types, ONLY: pint_env_type, pint_env_type, staging_env_type,&
       pint_env_create, pint_retain, pint_release,&
       pint_x2u, pint_u2x, pint_f2uf, pint_calc_uf_h,&
       pint_calc_f, pint_init_v, pint_write_centroids, pint_write_ener,&
       pint_write_restart, pint_step, pint_start, pint_run,&
       staging_env_create, staging_x2u, staging_u2x, staging_f2uf,&
       staging_init_masses
  USE input_section_types, ONLY: section_vals_type, section_vals_val_get, &
       section_vals_val_set, section_vals_write, section_vals_get_subs_vals,&
       section_vals_release, section_vals_get
  USE replica_methods, ONLY: rep_env_create, rep_env_calc_e_f
  USE replica_types, ONLY: replica_env_type, rep_env_retain, rep_env_release
  USE input_constants, ONLY: do_qs, mol_dyn_run
#include "cp_common_uses.h"

  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='pint_methods'

  public :: do_pint_run
!!***
!****************************************************************************
contains


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
  TYPE(cp_para_env_type), pointer :: para_env
  TYPE(section_vals_type), pointer :: input
  type(cp_error_type), optional, intent(inout) :: error

  logical :: failure
  character(len=*), parameter :: routineN='rep_env_test',&
       routineP=moduleN//':'//routineN
  TYPE(replica_env_type), POINTER :: rep_env
  INTEGER :: i, stat
  REAL(kind=dp), dimension(:,:), pointer :: f1,f2,f3
  REAL(kind=dp), dimension(:), pointer :: e1,e2,e3
  real(kind=dp) :: err
  TYPE(cp_logger_type), pointer :: logger
  
  failure=.FALSE.

  CPPrecondition(ASSOCIATED(para_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(input),cp_failure_level,routineP,error,failure)
  CPPrecondition(para_env%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(input%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     logger => cp_error_get_logger(error)
     nullify(rep_env)
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
  TYPE(section_vals_type), pointer :: input
  TYPE(cp_para_env_type), pointer :: para_env
  type(cp_error_type), optional, intent(inout) :: error

  logical :: failure
  character(len=*), parameter :: routineN='pint_env_test',&
       routineP=moduleN//':'//routineN
  TYPE(pint_env_type), POINTER :: pint_env
  REAL(kind=dp) :: err, c, e_h
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: x1
  INTEGER :: i, stat,ib,idim
  TYPE(cp_logger_type), pointer :: logger
  failure=.false.

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
  TYPE(cp_para_env_type), pointer :: para_env
  TYPE(section_vals_type), pointer :: input
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='do_pint_run',&
        routineP=moduleN//':'//routineN
  TYPE(pint_env_type), pointer :: pint_env

  failure=.false.
  
  CPPrecondition(ASSOCIATED(para_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(input),cp_failure_level,routineP,error,failure)
  CPPrecondition(para_env%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(input%ref_count>0,cp_failure_level,routineP,error,failure)
  nullify(pint_env)
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
