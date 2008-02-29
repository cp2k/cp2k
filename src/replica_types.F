!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief types used to handle many replica of the same system that differ only
!>      in atom positions, and velocity.
!>      This is useful for things like path integrals or nudged elastic band
!> \note
!>      this is a stupid implementation that replicates all the information
!>      about the replicas, if you really want to do a *lot* of replicas on
!>      a lot of processors you should think about distributiong also that
!>      information
!> \par History
!>      09.2005 created [fawzi]
!> \author fawzi
! *****************************************************************************
MODULE replica_types
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_para_env,                     ONLY: cp_cart_create,&
                                             cp_cart_release,&
                                             cp_cart_write,&
                                             cp_para_env_create,&
                                             cp_para_env_release,&
                                             cp_para_env_write
  USE cp_para_types,                   ONLY: cp_para_cart_type,&
                                             cp_para_env_type
  USE cp_result_types,                 ONLY: cp_result_p_type,&
                                             cp_result_release,&
                                             cp_result_type,&
                                             put_results
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set,&
                                             section_vals_write
  USE kinds,                           ONLY: default_path_length,&
                                             dp
  USE message_passing,                 ONLY: MPI_COMM_NULL,&
                                             mp_alltoall,&
                                             mp_cart_create,&
                                             mp_cart_sub,&
                                             mp_sum,&
                                             mp_sync
  USE qs_wf_history_types,             ONLY: qs_wf_history_p_type,&
                                             wfi_release
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  LOGICAL, SAVE, PRIVATE :: module_initialized=.FALSE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'replica_types'
  INTEGER, SAVE, PRIVATE :: last_rep_env_id=0

  PUBLIC :: replica_env_type,replica_env_p_type
  PUBLIC :: rep_env_retain, rep_env_release,rep_env_write
  PUBLIC :: rep_env_sync, rep_env_calc_e_f,rep_env_sync_type
  PUBLIC :: rep_envs_get_rep_env, rep_env_create,rep_env_local_index

! *****************************************************************************
!> \brief keeps replicated information about the replicas
!> \param ref_count reference count
!> \param id_nr identity number (unique or each replica_env)
!> \param nrep number of replicas
!> \param nat number of atoms (in each replica)
!> \param ndim number of dimensions (3*nat)
!> \param f_env_id id of the force env that will do the calculations for the
!>        replicas owned by this processor
!> \param r ,v,f: positions, velocities and forces of the replicas.
!>        the indexing is as follow (idir,iat,irep)
!> \param replica_owner which processor column owns the replica irep
!> \param cart 2d distribution of the processors for the replicas,
!>        a column  work together on the same force_env (i.e. changing the
!>        row you stay in the same replica), rows have different replicas
!> \param para_env the global para env that contains all the replicas,
!>        this is just the cart as para_env
!> \param para_env_f parallel environment of the underlying force
!>        environment
!> \param row_rank mapping row -> rank in para_env_f
!> \param para_env_inter_rep parallel environment between processors
!>        on the same row of cart
!> \param col_rank mapping column -> rank in para_env_inter_rep
!> \param local_rep_indices indices of the local replicas, starting at 1
!> \param wf_history wavefunction history for the owned replicas
!> \param keep_wf_history if the wavefunction history for the owned replicas
!>        should be kept
!> \author fawzi
! *****************************************************************************
  TYPE replica_env_type
     INTEGER :: ref_count, id_nr, f_env_id, nrep,nat,ndim
     REAL(kind=dp), DIMENSION(:,:), POINTER :: r,v,f
     LOGICAL :: sync_v, keep_wf_history
     TYPE(qs_wf_history_p_type), DIMENSION(:), POINTER :: wf_history
     TYPE(cp_result_p_type),DIMENSION(:),POINTER :: results
     INTEGER, DIMENSION(:), POINTER :: local_rep_indices
     INTEGER, DIMENSION(:), POINTER :: replica_owner,row_rank,col_rank
     TYPE(cp_para_cart_type), POINTER :: cart
     TYPE(cp_para_env_type), POINTER :: para_env, para_env_f,para_env_inter_rep
  END TYPE replica_env_type

! *****************************************************************************
!> \brief ****s* replica_types/replica_env_p_type *
!> 
!>      to build arrays of pointers to a replica_env_type
!> \param rep_env the pointer to the replica_env
!> \author fawzi
! *****************************************************************************
  TYPE replica_env_p_type
     TYPE(replica_env_type), POINTER :: rep_env
  END TYPE replica_env_p_type

  TYPE(replica_env_p_type), POINTER, DIMENSION(:), PRIVATE :: rep_envs

CONTAINS

! *****************************************************************************
!> \brief creates a replica environment together with its force environment
!> \param rep_env the replica environment that will be created
!> \param para_env the parallel enviroment that will contain the replicas
!> \param input the input used to initialize the force environment
!> \param nrep the number of replicas to calculate
!> \param prep the number of processors for each replica
!> \param sync_v if the volocity should be synchronized (defaults to false)
!> \param keep_wf_history if wf history should be kept on a per replica
!>        basis (defaults to true for QS jobs)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE rep_env_create(rep_env, para_env, input, nrep, prep,&
       sync_v,keep_wf_history,error)
    TYPE(replica_env_type), POINTER          :: rep_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: input
    INTEGER                                  :: nrep, prep
    LOGICAL, INTENT(in), OPTIONAL            :: sync_v, keep_wf_history
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rep_env_create', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_path_length)       :: input_file_path, &
                                                output_file_path
    INTEGER :: comm_cart, comm_f, comm_inter_rep, i, i0, ierr, ip, ir, irep, &
      lp, my_prep, natom, new_env_id, nrep_local, stat, unit_nr
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: gridinfo
    INTEGER, DIMENSION(2)                    :: dims, pos
    LOGICAL                                  :: failure, ionode
    LOGICAL, DIMENSION(2)                    :: rdim
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_cart_type), POINTER         :: cart
    TYPE(cp_para_env_type), POINTER          :: para_env_f, para_env_full, &
                                                para_env_inter_rep

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(rep_env),cp_failure_level,routineP,error,failure)
    NULLIFY(cart,para_env_f,para_env_inter_rep)
    IF (.NOT. failure) THEN
       logger => cp_error_get_logger(error)
       ionode = para_env%mepos==para_env%source
       IF (ionode) THEN
          unit_nr=cp_logger_get_default_unit_nr(logger)
       ELSE
          unit_nr=-1
       END IF
       new_env_id=-1
       my_prep=MIN(prep,para_env%num_pe)
       dims(1)=my_prep
       dims(2)=MIN(para_env%num_pe/my_prep,nrep)
       IF ((dims(1)*dims(2)/=para_env%num_pe).AND.(unit_nr>0)) THEN
          WRITE(unit_nr,FMT="(T2,A)") "REPLICA| WARNING: number of processors is not divisible by the number of replicas"
          WRITE(unit_nr,FMT="(T2,A,I0,A)") "REPLICA| ",para_env%num_pe-dims(1)*dims(2)," MPI process(es) will be idle"
       END IF
       CALL mp_cart_create ( comm_old=para_env%group, ndims=2, dims=dims, pos=pos, comm_cart=comm_cart)
       IF (comm_cart/=MPI_COMM_NULL) THEN
          CALL cp_cart_create(cart,comm_cart,ndims=2,owns_group=.TRUE.,error=error)
          NULLIFY(para_env_full)
          CALL cp_para_env_create(para_env_full,comm_cart,owns_group=.FALSE.,error=error)
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
    IF (.NOT. failure) THEN
       ALLOCATE(gridinfo(2,0:para_env%num_pe-1))
       gridinfo=0
       gridinfo(:,para_env%mepos)=pos
       CALL mp_sum(gridinfo,para_env%group)
       IF (unit_nr>0) THEN
          WRITE(unit_nr,FMT="(T2,A,T71,I10)") "REPLICA| layout of the replica grid, number of groups ",dims(2)
          WRITE(unit_nr,FMT="(T2,A,T71,I10)") "REPLICA| layout of the replica grid, size of each group",dims(1)
          WRITE(unit_nr,FMT="(T2,A)",ADVANCE="NO") "REPLICA| MPI process to grid (rank,group) correspondence:"
          DO i=0,para_env%num_pe-1
             IF (MODULO(i,4)==0) WRITE(unit_nr,*)
             WRITE(unit_nr,FMT='(A3,I4,A3,I4,A1,I4,A1)',ADVANCE="NO")&
                  "  (",i," : ",gridinfo(1,i),",",gridinfo(2,i),")"
          END DO
          WRITE(unit_nr,*)
       ENDIF
       DEALLOCATE(gridinfo)
    ENDIF
    IF (.NOT.failure.AND.ASSOCIATED(rep_env)) THEN
       last_rep_env_id=last_rep_env_id+1
       rep_env%id_nr=last_rep_env_id
       rep_env%ref_count=1
       rep_env%nrep=nrep
       rep_env%sync_v=.FALSE.
       IF (PRESENT(sync_v)) rep_env%sync_v=sync_v
       rep_env%keep_wf_history=.TRUE.
       IF (PRESENT(keep_wf_history)) rep_env%keep_wf_history=keep_wf_history
       NULLIFY(rep_env%wf_history)
       NULLIFY(rep_env%results)

       ALLOCATE(rep_env%row_rank(0:cart%num_pe(1)-1),rep_env%col_rank(0:cart%num_pe(2)-1),&
            stat=stat)
       rep_env%row_rank=0
       rep_env%row_rank(cart%mepos(1))=para_env_f%mepos
       CALL mp_sum(rep_env%row_rank,para_env_f%group)
       rep_env%col_rank=0
       rep_env%col_rank(cart%mepos(2))=para_env_inter_rep%mepos
       CALL mp_sum(rep_env%col_rank,para_env_inter_rep%group)

       CALL section_vals_val_get(input,"GLOBAL%PROJECT_NAME",&
            c_val=input_file_path,error=error)
       lp=LEN_TRIM(input_file_path)
       input_file_path(lp+1:LEN(input_file_path))="-r-"//&
            ADJUSTL(cp_to_string(cart%mepos(2)))
       lp=LEN_TRIM(input_file_path)
            CALL section_vals_val_set(input,"GLOBAL%PROJECT_NAME",&
                 c_val=input_file_path,error=error) !needed, or using a differen iteration is enough?
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
       CALL cp_create_fenv_comm(new_env_id,input_file_path,output_file_path,&
            para_env_f%group,ierr)
       CPAssert(ierr==0,cp_failure_level,routineP,error,failure)

       IF (.NOT.failure) THEN
          rep_env%f_env_id=new_env_id
          CALL cp_get_natom(new_env_id,natom,ierr)
          CPPostcondition(ierr==0,cp_fatal_level,routineP,error,failure)
          rep_env%nat=natom
          rep_env%ndim=3*natom
          ALLOCATE(rep_env%replica_owner(nrep),stat=stat)
          CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)

          i0=nrep/cart%num_pe(2)
          ir=MODULO(nrep,cart%num_pe(2))
          DO ip=0,cart%num_pe(2)-1
             DO i=i0*ip+MIN(ip,ir)+1,i0*(ip+1)+MIN(ip+1,ir)
                rep_env%replica_owner(i)=ip
             END DO
          END DO

          nrep_local=i0
          IF (cart%mepos(2)<ir) nrep_local=nrep_local+1
          ALLOCATE(rep_env%local_rep_indices(nrep_local),stat=stat)
          CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
          nrep_local=0
          DO irep=1,nrep
             IF (rep_env%replica_owner(irep)==cart%mepos(2)) THEN
                nrep_local=nrep_local+1
                rep_env%local_rep_indices(nrep_local)=irep
             END IF
          END DO
          CPPostcondition(nrep_local==SIZE(rep_env%local_rep_indices),cp_failure_level,routineP,error,failure)

          rep_env%cart => cart
          rep_env%para_env => para_env_full
          rep_env%para_env_f => para_env_f
          rep_env%para_env_inter_rep => para_env_inter_rep

          ALLOCATE(rep_env%r(rep_env%ndim,nrep),rep_env%v(rep_env%ndim,nrep),&
               rep_env%f(rep_env%ndim+1,nrep), stat=stat)
          CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)

          rep_env%r=0._dp
          rep_env%f=0._dp
          rep_env%v=0._dp
          CALL cp_set_vel(rep_env%f_env_id, rep_env%v(:,1), rep_env%ndim, ierr)
          CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
          DO i=1,nrep
             IF (cart%mepos(2)==rep_env%replica_owner(i)) THEN
                CALL cp_get_pos(rep_env%f_env_id, &
                     rep_env%r(:,i), rep_env%ndim, ierr)
                CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
             END IF
          END DO
       END IF
    END IF
    IF (.NOT.failure.AND.ASSOCIATED(rep_env)) THEN
       CALL rep_envs_add_rep_env(rep_env,error=error)
       CALL cp_rep_init(rep_env%id_nr,ierr)
       CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
  END SUBROUTINE rep_env_create

! *****************************************************************************
!> \brief releases the given replica environment
!> \param rep_env the replica environment to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      here and not in replica_types to allow the use of replica_env_type
!>      in a force_env (call to destroy_force_env gives circular dep)
!> \author fawzi
! *****************************************************************************
  SUBROUTINE rep_env_release(rep_env,error)
    TYPE(replica_env_type), POINTER          :: rep_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rep_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ierr, stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    CALL timeset(routineN,"I","",handle)
    IF (ASSOCIATED(rep_env)) THEN
       CPPrecondition(rep_env%ref_count>0,cp_failure_level,routineP,error,failure)
       rep_env%ref_count=rep_env%ref_count-1
       IF (rep_env%ref_count==0) THEN
          CALL cp_rep_destroy(rep_env%id_nr,ierr)
          IF (rep_env%f_env_id>0) THEN
             CALL cp_destroy_fenv(rep_env%f_env_id,ierr)
             CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(rep_env%r)) THEN
             DEALLOCATE(rep_env%r,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (ASSOCIATED(rep_env%v)) THEN
             DEALLOCATE(rep_env%v,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (ASSOCIATED(rep_env%f)) THEN
             DEALLOCATE(rep_env%f,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (ASSOCIATED(rep_env%wf_history)) THEN
             DO i=1,SIZE(rep_env%wf_history)
                CALL wfi_release(rep_env%wf_history(i)%wf_history,error=error)
             END DO
             DEALLOCATE(rep_env%wf_history,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (ASSOCIATED(rep_env%results)) THEN
             DO i=1,SIZE(rep_env%results)
                CALL cp_result_release(rep_env%results(i)%results,error=error)
             END DO
             DEALLOCATE(rep_env%results,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          DEALLOCATE(rep_env%local_rep_indices,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          IF (ASSOCIATED(rep_env%replica_owner)) THEN
             DEALLOCATE(rep_env%replica_owner,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          DEALLOCATE(rep_env%col_rank,rep_env%row_rank,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          CALL cp_cart_release(rep_env%cart,error=error)
          CALL cp_para_env_release(rep_env%para_env,error=error)
          CALL cp_para_env_release(rep_env%para_env_f,error=error)
          CALL cp_para_env_release(rep_env%para_env_inter_rep,error=error)
          CALL rep_envs_rm_rep_env(rep_env,error=error)
          DEALLOCATE(rep_env,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF
    NULLIFY(rep_env)
    CALL timestop(handle)
  END SUBROUTINE rep_env_release

! *****************************************************************************
!> \brief retains the given replica environment
!> \param rep_env the replica environment to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE rep_env_retain(rep_env,error)
    TYPE(replica_env_type), POINTER          :: rep_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rep_env_retain', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    failure=.FALSE.
    CALL timeset(routineN,"I","",handle)
    CPPrecondition(ASSOCIATED(rep_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPrecondition(rep_env%ref_count>0,cp_failure_level,routineP,error,failure)
       rep_env%ref_count=rep_env%ref_count+1
    END IF
    CALL timestop(handle)
  END SUBROUTINE rep_env_retain

! *****************************************************************************
!> \brief writes out information about the rep_env
!> \param rep_env the replica env to describe
!> \param unit_nr the unit to write to
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE rep_env_write(rep_env, unit_nr, error)
    TYPE(replica_env_type), POINTER          :: rep_env
    INTEGER, INTENT(in)                      :: unit_nr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rep_env_write', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    failure=.FALSE.

    CALL timeset(routineN,"I","",handle)
    IF (ASSOCIATED(rep_env)) THEN
       WRITE(unit_nr,"('<replica_env, id_nr=',i6,' ref_count=',i6,'> sync_v=',l1,',')")&
            rep_env%id_nr,rep_env%ref_count, rep_env%sync_v
       WRITE (unit_nr,"(' f_env_id=',i6,', nrep=',i6,', nat=',i8,', ndim=',i8,',')")&
            rep_env%f_env_id, rep_env%nrep,rep_env%nat, rep_env%ndim
       WRITE (unit_nr,"(' replica_owner=')", advance="no")
       WRITE (unit_nr,"(10i6)") rep_env%replica_owner
       WRITE (unit_nr,"(' cart=')", advance="no")
       CALL cp_cart_write(rep_env%cart,unit_nr,error=error)
       WRITE (unit_nr,"(' para_env=')", advance="no")
       CALL cp_para_env_write(rep_env%para_env,unit_nr,error=error)
       WRITE (unit_nr,"(' para_env_f=')", advance="no")
       CALL cp_para_env_write(rep_env%para_env_f,unit_nr,error=error)
       WRITE (unit_nr,"(' para_env_inter_rep=')", advance="no")
       CALL cp_para_env_write(rep_env%para_env_inter_rep,unit_nr,error=error)
       WRITE (unit_nr,"(' col_rank=(')", advance="no")
       WRITE (unit_nr,"(10i6)") rep_env%col_rank
       WRITE (unit_nr,"(')')")
       WRITE (unit_nr,"(' row_rank=(')", advance="no")
       WRITE (unit_nr,"(10i6)") rep_env%row_rank
       WRITE (unit_nr,"(')')")
       WRITE (unit_nr,"(' r=(')", advance="no")
       WRITE (unit_nr,"(3es12.5)") rep_env%r
       WRITE (unit_nr,"(')')")
       WRITE (unit_nr,"(' v=(')", advance="no")
       WRITE (unit_nr,"(3es12.5)") rep_env%v
       WRITE (unit_nr,"(')')")
       WRITE (unit_nr,"(' v=(')", advance="no")
       WRITE (unit_nr,"(3es12.5)") rep_env%f
       WRITE (unit_nr,"(')')")
       WRITE (unit_nr,"(' keep_wf_history=',l1,', associated(wf_history)=',l1,',')")&
            rep_env%keep_wf_history,ASSOCIATED(rep_env%wf_history)
       WRITE (unit_nr,"('</replica_env>')")
    ELSE
       WRITE(unit_nr,"('<replica_env=*null*/>')")
    END IF
    CALL timestop(handle)
  END SUBROUTINE rep_env_write

! *****************************************************************************
!> \brief sends the data from each replica to all the other
!>      on replica j/=i data from replica i overwrites val(:,i)
!> \param rep_env replica environment
!> \param vals the values to synchronize (second index runs over replicas)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      could be optimized: bcast in inter_rep, all2all or shift vs sum
!> \author fawzi
! *****************************************************************************
  SUBROUTINE rep_env_sync(rep_env,vals,error)
    TYPE(replica_env_type), POINTER          :: rep_env
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(inout)                          :: vals
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rep_env_sync', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, irep
    LOGICAL                                  :: failure

    failure=.FALSE.
    CALL timeset(routineN,"I","",handle)
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
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE rep_env_sync

! *****************************************************************************
!> \brief sends the data from each replica to all the other
!>      in this case the result type is passed
!> \param rep_env replica environment
!> \param results is an array of result_types
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fschiff
! *****************************************************************************
  SUBROUTINE rep_env_sync_type(rep_env,results,error)
    TYPE(replica_env_type), POINTER          :: rep_env
    TYPE(cp_result_p_type), DIMENSION(:), &
      POINTER                                :: results
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rep_env_sync_type', &
      routineP = moduleN//':'//routineN

    CHARACTER(30), ALLOCATABLE, DIMENSION(:) :: rbuf_list
    INTEGER                                  :: handle, irep, nrep, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: rbuf_map, rcount, &
                                                rdispl_list, rdispl_map, &
                                                scount, sdispl, size_list, &
                                                size_map
    LOGICAL                                  :: failure
    TYPE(cp_result_type), POINTER            :: loc_res

    failure=.FALSE.
    CALL timeset(routineN,"I","",handle)
    nrep=rep_env%nrep
    NULLIFY(loc_res)
    CPPrecondition(ASSOCIATED(rep_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(rep_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(results)==rep_env%nrep,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       ALLOCATE(size_list(nrep),stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ALLOCATE(size_map(nrep),stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ALLOCATE(scount(nrep),stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ALLOCATE(rcount(nrep),stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)     
       ALLOCATE(sdispl(nrep),stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ALLOCATE(rdispl_map(nrep),stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error) 
       ALLOCATE(rdispl_list(nrep),stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error) 
       size_list=0
       size_map=0
       DO irep=1,nrep
          IF (rep_env%replica_owner(irep)==rep_env%cart%mepos(2)) THEN
             size_list(irep)=SIZE(results(irep)%results%result_list)
             size_map(irep)=SIZE(results(irep)%results%desc_map)        
             loc_res=>results(irep)%results
          END IF
       END DO
       CALL mp_sync(rep_env%para_env_inter_rep%group)
       CALL mp_sum(size_list,rep_env%para_env_inter_rep%group)
       CALL mp_sum(size_map,rep_env%para_env_inter_rep%group)
       CALL mp_sync(rep_env%para_env_inter_rep%group)
       ALLOCATE(rbuf_map(SUM(size_map)),stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ALLOCATE(rbuf_list(SUM(size_list)),stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       scount(:)=SIZE(loc_res%desc_map)
       sdispl(:)=0
       rcount(:)=size_map(:)
       rdispl_map(1)=0
       DO irep=2,nrep
          rdispl_map(irep)=rdispl_map(irep-1)+rcount(irep-1)
       END DO

       CALL mp_sync(rep_env%para_env_inter_rep%group)
       CALL mp_alltoall(loc_res%desc_map,scount,sdispl,rbuf_map,rcount,rdispl_map,&
            rep_env%para_env_inter_rep%group)

       scount(:)=SIZE(loc_res%result_list)
       rcount(:)=size_list(:)
       rdispl_list(1)=0

       DO irep=2,nrep
          rdispl_list(irep)=rdispl_list(irep-1)+rcount(irep-1)
       END DO
       CALL mp_sync(rep_env%para_env_inter_rep%group)
       CALL mp_alltoall(loc_res%result_list,scount,sdispl,rbuf_list,rcount,rdispl_list,&
            LEN(loc_res%result_list(1)),rep_env%para_env_inter_rep%group)     

       DO irep=1,nrep
          CALL put_results(results(irep)%results,rbuf_list,rbuf_map,size_list(irep),size_map(irep),&
               rdispl_list(irep),rdispl_map(irep),error)
       END DO

       DEALLOCATE(size_list,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       DEALLOCATE(size_map,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       DEALLOCATE(scount,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       DEALLOCATE(rcount,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)     
       DEALLOCATE(sdispl,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       DEALLOCATE(rdispl_list,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       DEALLOCATE(rdispl_map,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)     
       DEALLOCATE(rbuf_map,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       DEALLOCATE(rbuf_list,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE rep_env_sync_type

! *****************************************************************************
!> \brief evaluates the forces
!> \param rep_env the replica environment on which you want to evaluate the
!>        forces
!> \param calc_f if true calculates also the forces, if false only the
!>        energy
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      indirect through f77_int_low to work around fortran madness
!> \author fawzi
! *****************************************************************************
  SUBROUTINE rep_env_calc_e_f(rep_env,calc_f,error)
    TYPE(replica_env_type), POINTER          :: rep_env
    LOGICAL, OPTIONAL                        :: calc_f
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rep_env_calc_e_f', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, my_calc_f
    LOGICAL                                  :: failure

    failure=.FALSE.

    CALL timeset(routineN,"I","",handle)
    CPPrecondition(ASSOCIATED(rep_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(rep_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       my_calc_f=0
       IF (PRESENT(calc_f)) THEN
          IF (calc_f) my_calc_f=1
       END IF
       CALL cp_rep_calc_e_f(rep_env%id_nr,my_calc_f,ierr)
       CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(handle)
  END SUBROUTINE rep_env_calc_e_f

! *****************************************************************************
!> \brief returns the replica environment with the given id_nr
!> \param id_nr the id_nr of the requested rep_envs
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  FUNCTION rep_envs_get_rep_env(id_nr,error) RESULT(res)
    INTEGER, INTENT(in)                      :: id_nr
    TYPE(cp_error_type), INTENT(inout)       :: error
    TYPE(replica_env_type), POINTER          :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'rep_envs_get_rep_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

    NULLIFY(res)
    IF (module_initialized) THEN
       IF (ASSOCIATED(rep_envs)) THEN
          DO i=1,SIZE(rep_envs)
             IF (rep_envs(i)%rep_env%id_nr==id_nr) THEN
                res => rep_envs(i)%rep_env
                EXIT
             END IF
          END DO
       END IF
    END IF
  END FUNCTION rep_envs_get_rep_env

! *****************************************************************************
!> \brief adds the given rep_env to the list of controlled rep_envs.
!> \param rep_env the rep_env to add
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE rep_envs_add_rep_env(rep_env,error)
    TYPE(replica_env_type), POINTER          :: rep_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rep_envs_add_rep_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure
    TYPE(replica_env_p_type), DIMENSION(:), &
      POINTER                                :: new_rep_envs
    TYPE(replica_env_type), POINTER          :: rep_env2

    failure=.FALSE.

    IF (ASSOCIATED(rep_env)) THEN
       rep_env2 => rep_envs_get_rep_env(rep_env%id_nr,error=error)
       IF (.NOT.ASSOCIATED(rep_env2)) THEN
          IF (module_initialized) THEN
             IF (.NOT.ASSOCIATED(rep_envs)) THEN
                ALLOCATE(rep_envs(1),stat=stat)
                CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
             ELSE
                ALLOCATE(new_rep_envs(SIZE(rep_envs)+1),stat=stat)
                CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
                DO i=1,SIZE(rep_envs)
                   new_rep_envs(i)%rep_env => rep_envs(i)%rep_env
                END DO
                DEALLOCATE(rep_envs,stat=stat)
                CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
                rep_envs => new_rep_envs
             END IF
          ELSE
             ALLOCATE(rep_envs(1),stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          rep_envs(SIZE(rep_envs))%rep_env => rep_env
          module_initialized=.TRUE.
       END IF
    END IF
  END SUBROUTINE rep_envs_add_rep_env

! *****************************************************************************
!> \brief removes the given rep_env to the list of controlled rep_envs.
!> \param rep_env the rep_env to remove
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE rep_envs_rm_rep_env(rep_env,error)
    TYPE(replica_env_type), POINTER          :: rep_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rep_envs_rm_rep_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ii, stat
    LOGICAL                                  :: failure
    TYPE(replica_env_p_type), DIMENSION(:), &
      POINTER                                :: new_rep_envs

    failure=.FALSE.

    IF (ASSOCIATED(rep_env)) THEN
       CPPrecondition(module_initialized,cp_failure_level,routineP,error,failure)
       ALLOCATE(new_rep_envs(SIZE(rep_envs)-1),stat=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
       ii=0
       DO i=1,SIZE(rep_envs)
          IF (rep_envs(i)%rep_env%id_nr/=rep_env%id_nr) THEN
             ii=ii+1
             new_rep_envs(ii)%rep_env => rep_envs(i)%rep_env
          END IF
       END DO
       CPPostcondition(ii==SIZE(new_rep_envs),cp_failure_level,routineP,error,failure)
       DEALLOCATE(rep_envs,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       rep_envs => new_rep_envs
       IF (SIZE(rep_envs)==0) THEN
          DEALLOCATE(rep_envs,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF
  END SUBROUTINE rep_envs_rm_rep_env

! *****************************************************************************
!> \brief returns the local index of the replica (-1 if it is not a local replica)
!> \param rep_env the replica env
!> \param global_index the global replica index
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  FUNCTION rep_env_local_index(rep_env,global_index,error) RESULT(res)
    TYPE(replica_env_type), POINTER          :: rep_env
    INTEGER, INTENT(in)                      :: global_index
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'rep_env_local_index', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i
    LOGICAL                                  :: failure

    failure=.FALSE.

    CALL timeset(routineN,handle)
    CPPrecondition(ASSOCIATED(rep_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(rep_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       res=-1
       DO i=1,SIZE(rep_env%local_rep_indices)
          IF (rep_env%local_rep_indices(i)==global_index) THEN
             res=i
             EXIT
          END IF
       END DO
       IF (res==-1) THEN
          PRINT *,routineP," ",global_index," not in ",rep_env%local_rep_indices
       END IF
    END IF
    CALL timestop(handle)
  END FUNCTION rep_env_local_index

END MODULE replica_types
