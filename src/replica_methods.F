!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****h* cp2k/replica_methods *
!!
!!   NAME
!!     replica_methods
!!
!!   FUNCTION
!!     methods to setup replicas of the same system differing only by atom
!!     positions and velocities (as used in path integral or nudged elastic
!!     band for example)
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     09.2005 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE replica_methods
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_output_handling,              ONLY: cp_add_iter_level
  USE cp_para_env,                     ONLY: cp_cart_create,&
                                             cp_para_env_create,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_cart_type,&
                                             cp_para_env_type
  USE f77_interface,                   ONLY: &
       calc_force, create_force_env, f_env_add_defaults, f_env_rm_defaults, &
       f_env_type, get_force, get_natom, get_pos, set_vel, destroy_force_env
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set,&
                                             section_vals_write
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE message_passing,                 ONLY: MPI_COMM_NULL,&
                                             mp_cart_create,&
                                             mp_cart_sub,&
                                             mp_sum
  USE replica_types,                   ONLY: replica_env_type,&
       rep_env_sync, rep_env_write, rep_env_create,&
       rep_env_calc_e_f, rep_env_release, rep_envs_get_rep_env
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE cp_control_types, ONLY: dft_control_type
  USE qs_environment_types, ONLY: qs_environment_type, get_qs_env, set_qs_env
  USE force_env_types, ONLY: force_env_get,use_qs_force
  USE qs_wf_history_types, ONLY: qs_wf_history_type,wfi_release,wfi_retain
  USE qs_wf_history_methods, only: wfi_create
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'replica_methods'

  PUBLIC :: rep_env_init_low,rep_env_test, rep_env_calc_e_f_low
!!***
!****************************************************************************
CONTAINS

!!****f* replica_methods/rep_env_init_low *
!!
!!   NAME
!!     rep_env_init_low
!!
!!   FUNCTION
!!     finishes the low level initialization of the replica env
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - rep_env_id: id_nr of the replica environment that should be initialized
!!     - ierr: will be non zero if there is an initialization error
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE rep_env_init_low(rep_env_id, ierr)
  INTEGER, INTENT(in) :: rep_env_id
  INTEGER, intent(out) :: ierr

    TYPE(replica_env_type), POINTER          :: rep_env

    CHARACTER(len=*), PARAMETER :: routineN = 'rep_env_init_low', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_path_length)       :: input_file_path, &
                                                output_file_path
    INTEGER                                  :: comm_cart, comm_f, &
                                                comm_inter_rep, i, i0, &
                                                ip, ir, lp, my_prep, natom, &
                                                new_env_id, stat, unit_nr,&
                                                in_use
    INTEGER, DIMENSION(2)                    :: dims, pos
    LOGICAL                                  :: failure,first_rep
    LOGICAL, DIMENSION(2)                    :: rdim
    TYPE(cp_error_type)                      :: error
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_cart_type), POINTER         :: cart
    TYPE(cp_para_env_type), POINTER          :: para_env_f, para_env_inter_rep
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(qs_environment_type), pointer :: qs_env
    TYPE(qs_wf_history_type), pointer :: wf_history
    TYPE(dft_control_type), pointer :: dft_control

  failure=.FALSE.
  rep_env => rep_envs_get_rep_env(rep_env_id)
  CALL cp_assert(ASSOCIATED(rep_env),cp_failure_level,cp_assertion_failed,&
       routineP,"could not find rep_env with id_nr"//cp_to_string(rep_env_id),&
       failure=failure)
  NULLIFY(wf_history,qs_env,dft_control)
  NULLIFY(cart,para_env_f,para_env_inter_rep, f_env)
  IF (.NOT. failure) THEN
     CALL f_env_add_defaults(f_env_id=rep_env%f_env_id,f_env=f_env,&
          new_error=error, failure=failure)
     IF (.NOT.failure) THEN
        logger => cp_error_get_logger(error)
        cart => rep_env%cart
        para_env_f => rep_env%para_env_f
        para_env_inter_rep => rep_env%para_env_inter_rep
        logger%iter_info%iteration(1)=cart%mepos(2)
        CALL cp_add_iter_level(iteration_info=logger%iter_info,&
             level_name="REP_EVAL",error=error)
        !wf interp
        IF (rep_env%keep_wf_history) THEN
           CALL force_env_get(f_env%force_env,in_use=in_use,error=error)
           IF (in_use==use_qs_force) THEN
              CALL force_env_get(f_env%force_env,qs_env=qs_env,&
                   error=error)
              CALL get_qs_env(qs_env,dft_control=dft_control,error=error)
              ALLOCATE(rep_env%wf_history(size(rep_env%local_rep_indices)),stat=stat)
              CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
              DO i=1,SIZE(rep_env%wf_history)
                 NULLIFY(rep_env%wf_history(i)%wf_history)
                 IF (i==1) THEN
                    CALL get_qs_env(qs_env,&
                         wf_history=rep_env%wf_history(i)%wf_history,&
                         error=error)
                    CALL wfi_retain(rep_env%wf_history(i)%wf_history,&
                         error=error)
                 ELSE
                    CALL wfi_create(rep_env%wf_history(i)%wf_history,&
                         interpolation_method_nr=&
                         dft_control%qs_control%wf_interpolation_method_nr,&
                         extrapolation_order = dft_control%qs_control%wf_extrapolation_order,&
                         error=error)
                 END IF
              END DO
           ELSE
              rep_env%keep_wf_history=.FALSE.
           END IF
        END IF
     END IF
     IF (.NOT.failure) THEN
        CALL rep_env_sync(rep_env,rep_env%r,error=error)
        CALL rep_env_sync(rep_env,rep_env%v,error=error)
        CALL rep_env_sync(rep_env,rep_env%f,error=error)
     END IF

     CALL f_env_rm_defaults(f_env,error,ierr)
     CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
  END IF
END SUBROUTINE rep_env_init_low
!***************************************************************************


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
    REAL(kind=dp)                            :: err,max_err
    REAL(kind=dp), DIMENSION(:), POINTER     :: e1, e2, e3
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: f1, f2, f3
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(replica_env_type), POINTER          :: rep_env
    CHARACTER(len=default_string_length)     :: project_name

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(para_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(input),cp_failure_level,routineP,error,failure)
  CPPrecondition(para_env%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(input%ref_count>0,cp_failure_level,routineP,error,failure)
  max_err=0._dp
  IF (.NOT. failure) THEN
     logger => cp_error_get_logger(error)
     NULLIFY(rep_env,e1,e2,e3,f1,f2,f3)
!FM     CALL cp_create_fenv_comm(i,"Ar-sys.inp","Ar-sys.out",&
!FM     0,stat)
!FM     CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
!FM     CALL cp_destroy_fenv(i,stat)
!FM     CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
     CALL section_vals_val_get(input,"GLOBAL%PROJECT_NAME",&
          c_val=project_name,error=error)

     CALL section_vals_val_set(input,"GLOBAL%PROJECT_NAME",&
          c_val=TRIM(project_name)//"t11",error=error)
     CALL rep_env_create(rep_env, para_env=para_env, input=input,&
          nrep=1,prep=1,error=error)
     CALL rep_env_write(rep_env,cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
          error=error)
     IF (ASSOCIATED(rep_env)) THEN
        CALL rep_env_calc_e_f(rep_env,error=error)
        ALLOCATE(f1(rep_env%ndim,rep_env%nrep),e1(rep_env%nrep),stat=stat)
        f1=rep_env%f(:rep_env%ndim,:)
        e1=rep_env%f(rep_env%ndim+1,:)
     END IF
     CALL rep_env_release(rep_env,error=error)

     CALL section_vals_val_set(input,"GLOBAL%PROJECT_NAME",&
          c_val=TRIM(project_name)//"t12",error=error)
     CALL rep_env_create(rep_env, para_env=para_env, input=input,&
          nrep=1,prep=2,error=error)
     CALL rep_env_write(rep_env,cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
          error=error)
     IF (ASSOCIATED(rep_env)) THEN
        CALL rep_env_calc_e_f(rep_env,error=error)
        ALLOCATE(f2(rep_env%ndim,rep_env%nrep),e2(rep_env%nrep),stat=stat)
        f2=rep_env%f(:rep_env%ndim,:)
        e2=rep_env%f(rep_env%ndim+1,:)
     END IF
     CALL rep_env_release(rep_env,error=error)

     CALL section_vals_val_set(input,"GLOBAL%PROJECT_NAME",&
          c_val=TRIM(project_name)//"t21",error=error)
     CALL rep_env_create(rep_env, para_env=para_env, input=input,&
          nrep=2,prep=1,error=error)
     CALL rep_env_write(rep_env,cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
          error=error)
     IF (ASSOCIATED(rep_env)) THEN
        rep_env%r(:,2)=rep_env%r(:,1)+0.5_dp
        CALL rep_env_calc_e_f(rep_env,error=error)
        ALLOCATE(f3(rep_env%ndim,rep_env%nrep),e3(rep_env%nrep),stat=stat)
        f3=rep_env%f(:rep_env%ndim,:)
        e3=rep_env%f(rep_env%ndim+1,:)
     END IF
     CALL rep_env_release(rep_env,error=error)

     CALL section_vals_val_set(input,"GLOBAL%PROJECT_NAME",&
          c_val=TRIM(project_name),error=error)

     IF (ASSOCIATED(f1).AND.ASSOCIATED(f2)) THEN
        err=0._dp
        DO i=1,SIZE(f1,1)
           err=MAX(err,ABS(f1(i,1)-f2(i,1)))
        END DO
        CALL cp_log(logger,cp_note_level+1,routineP,"diff_f12="//cp_to_string(err),&
             local=.FALSE.)
        max_err=MAX(err,max_err)
        err=ABS(e1(1)-e2(1))
        max_err=MAX(err,max_err)
        CALL cp_log(logger,cp_note_level+1,routineP,"diff_e12="//&
             cp_to_string(err)//cp_to_string(e1(1))//&
             cp_to_string(e2(1)), local=.FALSE.)
     END IF
     IF (ASSOCIATED(f1).AND.ASSOCIATED(f3)) THEN
        err=0._dp
        DO i=1,SIZE(f1,1)
           err=MAX(err,ABS(f1(i,1)-f3(i,1)))
        END DO
        CALL cp_log(logger,cp_note_level+1,routineP,"diff_f13="//cp_to_string(err),&
             local=.FALSE.)
        max_err=MAX(err,max_err)
        err=ABS(e1(1)-e3(1))
        max_err=MAX(err,max_err)
        CALL cp_log(logger,cp_note_level+1,routineP,"diff_e13="//&
             cp_to_string(err)//cp_to_string(e1(1))//&
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
        max_err=MAX(err,max_err)
        err=ABS(e3(1)-e3(2))
        max_err=MAX(err,max_err)
        CALL cp_log(logger,cp_note_level+1,routineP,"diff_e33="//&
             cp_to_string(err)//cp_to_string(e3(1))//&
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
     CALL cp_log(logger,cp_note_level+1,routineP,"CUMULATIVE ERROR"//&
          cp_to_string(max_err), local=.FALSE.)
     CPAssert(max_err<1.e-6,cp_failure_level,routineP,error,failure)
  END IF
END SUBROUTINE rep_env_test
!****************************************************************************

!!****f* replica_methods/rep_env_calc_e_f_int *
!!
!!   NAME
!!     rep_env_calc_e_f_int
!!
!!   FUNCTION
!!     calculates energy and force, internal private method
!!
!!   NOTES
!!     this is the where the real work is done
!!
!!   ARGUMENTS
!!     - rep_env: the replica env to update
!!     - if the force should be calculated as well (defaults to true)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE rep_env_calc_e_f_int(rep_env,calc_f,error)
    TYPE(replica_env_type), POINTER          :: rep_env
    LOGICAL, OPTIONAL                        :: calc_f
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rep_env_calc_e_f_int', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, irep, md_iter, my_calc_f,i
    LOGICAL                                  :: failure
    TYPE(cp_error_type)                      :: new_error
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(f_env_type), pointer :: f_env
    TYPE(qs_environment_type), pointer :: qs_env

  failure=.false.
  NULLIFY(f_env,qs_env)
  CPPrecondition(ASSOCIATED(rep_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(rep_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     my_calc_f=3*rep_env%nat
     IF (PRESENT(calc_f)) THEN
        IF (.NOT.calc_f) my_calc_f=0
     END IF
     
     CALL f_env_add_defaults(f_env_id=rep_env%f_env_id,f_env=f_env,&
          new_error=new_error, failure=failure)
     logger => cp_error_get_logger(new_error)
     md_iter=logger%iter_info%iteration(2)+1
     call f_env_rm_defaults(f_env,new_error,ierr)
     CPAssert(ierr==0,cp_failure_level,routineP,error,failure)

     DO i=1,SIZE(rep_env%local_rep_indices)
        irep=rep_env%local_rep_indices(i)
        IF (rep_env%sync_v) THEN
           CALL cp_set_vel(rep_env%f_env_id,rep_env%v(:,irep),&
                3*rep_env%nat,ierr)
           CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
           IF (failure) EXIT
        END IF
        
        logger%iter_info%iteration(1)=irep
        logger%iter_info%iteration(2)=md_iter
        
        IF (rep_env%keep_wf_history) THEN
           CALL f_env_add_defaults(f_env_id=rep_env%f_env_id,f_env=f_env,&
                new_error=new_error, failure=failure)
           CALL force_env_get(f_env%force_env,qs_env=qs_env,error=error)
           CALL set_qs_env(qs_env,&
                wf_history=rep_env%wf_history(i)%wf_history,&
                error=error)
           CALL f_env_rm_defaults(f_env,new_error,ierr)
           CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
        END IF
        
        CALL cp_calc_force(rep_env%f_env_id,rep_env%r,3*rep_env%nat,&
             rep_env%f(3*rep_env%nat+1,irep),rep_env%f(:3*rep_env%nat,irep),&
             my_calc_f,ierr)
        CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
        IF (failure) EXIT
     END DO
     IF (.NOT.failure) THEN
        CALL rep_env_sync(rep_env,rep_env%f,error=error)
     END IF
     CALL cp_error_synchronize_error(error,rep_env%para_env)
  END IF
     
END SUBROUTINE rep_env_calc_e_f_int
!***************************************************************************

!!****f* replica_methods/rep_env_calc_e_f_int *
!!
!!   NAME
!!     rep_env_calc_e_f_int
!!
!!   FUNCTION
!!     calculates energy and force, internal private method
!!
!!   NOTES
!!     low level wrapper to export this function in f77_int_low and work
!!     around the handling of circular dependecies in fortran 
!!
!!   ARGUMENTS
!!     - rep_env_id: the id if the replica environment in which energy and
!!       forces have to be evaluated
!!     - calc_f: if nonzero calculates also the forces along with the 
!!       energy
!!     - ierr: if an error happens this will be nonzero
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
RECURSIVE SUBROUTINE rep_env_calc_e_f_low(rep_env_id,calc_f,ierr)
  INTEGER, INTENT(in) :: rep_env_id, calc_f
  INTEGER, intent(out) :: ierr

  CHARACTER(len=*), PARAMETER :: routineN = 'rep_env_calc_e_f_int', &
       routineP = moduleN//':'//routineN

  logical :: failure
  TYPE(cp_error_type)                      :: new_error
  TYPE(cp_logger_type), POINTER            :: logger
  TYPE(replica_env_type), pointer :: rep_env
  TYPE(f_env_type), pointer :: f_env

  failure=.false.
  rep_env => rep_envs_get_rep_env(rep_env_id)
  IF (ASSOCIATED(rep_env)) THEN
     CALL f_env_add_defaults(f_env_id=rep_env%f_env_id,f_env=f_env,&
          new_error=new_error, failure=failure)
     CALL rep_env_calc_e_f_int(rep_env,calc_f/=0,new_error)
     CALL f_env_rm_defaults(f_env,new_error,ierr)
  ELSE
     ierr=111
  END IF
END SUBROUTINE rep_env_calc_e_f_low
!***************************************************************************

END MODULE replica_methods
