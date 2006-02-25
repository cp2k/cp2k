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
       f_env_type, get_force, get_natom, get_pos, set_vel
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set,&
                                             section_vals_write
  USE kinds,                           ONLY: default_path_length,&
                                             dp
  USE message_passing,                 ONLY: MPI_COMM_NULL,&
                                             mp_cart_create,&
                                             mp_cart_sub,&
                                             mp_sum
  USE replica_types,                   ONLY: replica_env_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'replica_methods'

  INTEGER, SAVE, PRIVATE :: last_rep_env_id=0

  PUBLIC :: rep_env_create, rep_env_sync, rep_env_calc_e_f
!!***
!****************************************************************************
CONTAINS

!!****f* replica_methods/rep_env_create *
!!
!!   NAME
!!     rep_env_create
!!
!!   FUNCTION
!!     creates a replica environment together with its force environment
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - rep_env: the replica environment that will be created
!!     - para_env: the parallel enviroment that will contain the replicas
!!     - input: the input used to initialize the force environment
!!     - nrep: the number of replicas to calculate
!!     - prep: the number of processors for each replica
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE rep_env_create(rep_env, para_env, input, nrep, prep, error)
    TYPE(replica_env_type), POINTER          :: rep_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: input
    INTEGER                                  :: nrep, prep
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rep_env_create', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_path_length)       :: input_file_path, &
                                                output_file_path
    INTEGER                                  :: comm_cart, comm_f, &
                                                comm_inter_rep, i, i0, ierr, &
                                                ip, ir, lp, my_prep, natom, &
                                                new_env_id, stat, unit_nr
    INTEGER, DIMENSION(2)                    :: dims, pos
    LOGICAL                                  :: failure
    LOGICAL, DIMENSION(2)                    :: rdim
    TYPE(cp_error_type)                      :: new_error
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_cart_type), POINTER         :: cart
    TYPE(cp_para_env_type), POINTER          :: para_env_f, para_env_inter_rep
    TYPE(f_env_type), POINTER                :: f_env

  failure=.FALSE.
  
  CPPrecondition(.NOT.ASSOCIATED(rep_env),cp_failure_level,routineP,error,failure)
  NULLIFY(cart,para_env_f,para_env_inter_rep, f_env)
  IF (.NOT. failure) THEN
     logger => cp_error_get_logger(error)
     new_env_id=-1
     my_prep=MIN(prep,para_env%num_pe)
     dims(1)=my_prep
     dims(2)=para_env%num_pe/my_prep
     IF (dims(1)*dims(2)/=para_env%num_pe) THEN
        unit_nr=cp_logger_get_default_unit_nr(logger,local=.FALSE.)
        WRITE(unit_nr,"(/,' WARNING',/,' number of processors not divisible in replica creation,')")
        WRITE(unit_nr,"(' wasting ',i4,' processors')") para_env%num_pe-dims(1)*dims(2)
     END IF
     CALL mp_cart_create ( comm_old=para_env%group, ndims=2, dims=dims, pos=pos, comm_cart=comm_cart)
     IF (comm_cart/=MPI_COMM_NULL) THEN
        CALL cp_cart_create(cart,comm_cart,ndims=2,owns_group=.TRUE.,error=error)
        rdim(1)=.TRUE.
        rdim(2)=.FALSE.
        CALL mp_cart_sub( comm=comm_cart, rdim=rdim , sub_comm=comm_f)
        CALL cp_para_env_create(para_env_f,comm_f,owns_group=.TRUE.,error=error)
        rdim(1)=.FALSE.
        rdim(2)=.TRUE.
        CALL mp_cart_sub( comm=comm_cart, rdim=rdim , sub_comm=comm_inter_rep)
        CALL cp_para_env_create(para_env_inter_rep,comm_inter_rep,&
             owns_group=.TRUE.,error=error)
        ALLOCATE(rep_env,stat=stat)
        CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
     END IF
  END IF
  IF (.not.failure.and.ASSOCIATED(rep_env)) THEN
     last_rep_env_id=last_rep_env_id+1
     rep_env%id_nr=last_rep_env_id
     rep_env%ref_count=1
     rep_env%nrep=nrep
     rep_env%sync_v=.TRUE.
     
     CALL section_vals_val_get(input,"GLOBAL%PROJECT_NAME",&
          c_val=input_file_path,error=error)
     lp=LEN_TRIM(input_file_path)
     input_file_path(lp+1:LEN(input_file_path))="-r-"//&
          ADJUSTL(cp_to_string(cart%mepos(2)))
     lp=LEN_TRIM(input_file_path)
     output_file_path=input_file_path(1:lp)//".out"
     input_file_path(lp+1:LEN(input_file_path))=".inp"
     CALL section_vals_val_set(input,"GLOBAL%OUTPUT_FILE_NAME",&
          c_val=TRIM(output_file_path),error=error)
     IF (para_env_f%source==para_env_f%mepos) THEN
        CALL open_file(file_name=TRIM(input_file_path),file_status="UNKNOWN",&
             file_form="FORMATTED",file_action="WRITE",&
             unit_number=unit_nr)
        CALL section_vals_write(input,unit_nr,hide_root=.TRUE.,error=error)
        CALL close_file(unit_nr)
     END IF
     CALL create_force_env(new_env_id=new_env_id,&
          input_file_path=input_file_path,&
          output_file_path=output_file_path,mpi_comm=para_env_f%group,ierr=ierr)
     CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
     IF (.NOT.failure) THEN
        CALL f_env_add_defaults(f_env_id=new_env_id,f_env=f_env,&
             new_error=new_error, failure=failure)
        IF (.NOT.failure) THEN
           logger => cp_error_get_logger(new_error)
           logger%iter_info%iteration(1)=cart%mepos(2)
           CALL cp_add_iter_level(iteration_info=logger%iter_info,&
                level_name="PINT_MD",error=error)
        END IF
        CALL f_env_rm_defaults(f_env,new_error,ierr)
        CPAssert(ierr==0,cp_failure_level,routineP,error,failure)

        rep_env%f_env_id=new_env_id
     END IF
     
     IF (.NOT.failure) THEN
        CALL get_natom(new_env_id,natom,ierr=ierr)
        CPPostcondition(ierr==0,cp_fatal_level,routineP,error,failure)
        rep_env%nat=natom
        rep_env%ndim=3*natom
        ALLOCATE(rep_env%replica_owner(nrep),stat=stat)
        CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
        
        i0=nrep/cart%num_pe(2)
        ir=MODULO(nrep,cart%num_pe(2))
        DO ip=1,cart%num_pe(2)
           DO i=i0*(ip-1)+MIN(ip,ir)+1,i0*ip+MIN(ip,ir)
              rep_env%replica_owner(i)=ip
           END DO
        END DO
        
        rep_env%cart => cart
        rep_env%para_env => para_env
        CALL cp_para_env_retain(rep_env%para_env,error=error)
        rep_env%para_env_f => para_env_f
        rep_env%para_env_inter_rep => para_env_inter_rep

        ALLOCATE(rep_env%r(rep_env%ndim,nrep),rep_env%v(rep_env%ndim,nrep),&
             rep_env%f(rep_env%ndim+1,nrep), stat=stat)
        CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)

        rep_env%f=0._dp
        rep_env%v=0._dp
        CALL set_vel(env_id=rep_env%f_env_id, &
             new_vel=rep_env%v(:,1), n_el=rep_env%ndim, ierr=ierr)
        CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
        DO i=1,nrep
           IF (cart%mepos(2)==rep_env%replica_owner(i)) THEN
              CALL get_pos(env_id=rep_env%f_env_id, &
                   pos=rep_env%r(:,i), n_el=rep_env%ndim, ierr=ierr)
              CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
           END IF
        END DO
        CALL rep_env_sync(rep_env,rep_env%r,error=error)
        CALL rep_env_sync(rep_env,rep_env%v,error=error)
        CALL rep_env_sync(rep_env,rep_env%f,error=error)
     END IF
  END IF
END SUBROUTINE rep_env_create
!***************************************************************************

!!****f* rep_env_methods/rep_env_sync *
!!
!!   NAME
!!     rep_env_sync
!!
!!   FUNCTION
!!     sends the data from each replica to all the other
!!     on replica j/=i data from replica i overwrites val(:,i)
!!
!!   NOTES
!!     could be optimized: all2all or shift vs sum
!!
!!   ARGUMENTS
!!     - rep_env: replica environment
!!     - vals: the values to synchronize (second index runs over replicas)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE rep_env_sync(rep_env,vals,error)
    TYPE(replica_env_type), POINTER          :: rep_env
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(inout)                          :: vals
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rep_env_sync', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: irep
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(rep_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(rep_env%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(SIZE(vals,2)==rep_env%nrep,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     DO irep=1,rep_env%nrep
        IF (rep_env%replica_owner(irep)/=rep_env%cart%mepos(2)) THEN
           vals(:,irep)=0._dp
        END IF
     END DO
     CALL mp_sum(vals,rep_env%para_env_inter_rep%group)
  END IF
END SUBROUTINE rep_env_sync
!***************************************************************************

!!****f* replica_methods/rep_env_evalf *
!!
!!   NAME
!!     rep_env_evalf
!!
!!   FUNCTION
!!     evaluates the forces
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - rep_env: the replica environment on which you want to evaluate the
!!       forces
!!     - calc_f: if true calculates also the forces, if false only the
!!       energy
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE rep_env_calc_e_f(rep_env,calc_f,error)
    TYPE(replica_env_type), POINTER          :: rep_env
    LOGICAL, OPTIONAL                        :: calc_f
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rep_env_calc_e_f', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, irep, md_iter, my_calc_f
    LOGICAL                                  :: failure
    TYPE(cp_error_type)                      :: new_error
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(f_env_type), POINTER                :: f_env

  failure=.FALSE.
  NULLIFY(f_env)

  CPPrecondition(ASSOCIATED(rep_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(rep_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     CALL f_env_add_defaults(f_env_id=rep_env%f_env_id,f_env=f_env,&
          new_error=new_error, failure=failure)
  END IF
  IF (.NOT.failure) THEN
     logger => cp_error_get_logger(new_error)
     md_iter=logger%iter_info%iteration(2)+1
     logger%iter_info%iteration(2)=md_iter
     CALL f_env_rm_defaults(f_env,new_error,ierr)
     CPAssert(ierr==0,cp_failure_level,routineP,error,failure)

     IF (.NOT. failure) THEN
        my_calc_f=3*rep_env%nat
        IF (PRESENT(calc_f)) THEN
           IF (.NOT.calc_f) my_calc_f=0
        END IF
        
        DO irep=1,rep_env%nrep
           IF (rep_env%replica_owner(irep)==rep_env%cart%mepos(2)) THEN
              IF (rep_env%sync_v) THEN
                 CALL set_vel(rep_env%f_env_id,rep_env%v(:,irep),&
                      3*rep_env%nat,ierr=ierr)
                 CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
                 IF (failure) EXIT
              END IF

              logger%iter_info%iteration(1)=irep
              logger%iter_info%iteration(2)=md_iter

              CALL calc_force(rep_env%f_env_id,rep_env%r,3*rep_env%nat,&
                   rep_env%f(3*rep_env%nat+1,irep),rep_env%f(:3*rep_env%nat,irep),&
                   3*rep_env%nat*my_calc_f,ierr)
              CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
              IF (failure) EXIT
              CALL get_force(rep_env%f_env_id,rep_env%f(:3*rep_env%nat,irep),&
                   my_calc_f,ierr=ierr)
              CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
              IF (failure) EXIT
           END IF
        END DO
        IF (.NOT.failure) THEN
           CALL rep_env_sync(rep_env,rep_env%f,error=error)
        END IF
     END IF
  END IF
END SUBROUTINE rep_env_calc_e_f
!***************************************************************************

END MODULE replica_methods
