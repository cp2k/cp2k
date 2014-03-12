!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
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
  USE cp_result_methods,               ONLY: cp_results_mp_bcast
  USE cp_result_types,                 ONLY: cp_result_p_type,&
                                             cp_result_release
  USE f77_blas
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
  USE qs_wf_history_types,             ONLY: qs_wf_history_p_type,&
                                             wfi_release
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  LOGICAL, SAVE, PRIVATE :: module_initialized=.FALSE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'replica_types'
  INTEGER, SAVE, PRIVATE :: last_rep_env_id=0

  PUBLIC :: replica_env_type,replica_env_p_type
  PUBLIC :: rep_env_retain, rep_env_release,rep_env_write
  PUBLIC :: rep_env_sync, rep_env_calc_e_f,rep_env_sync_results
  PUBLIC :: rep_envs_get_rep_env, rep_env_create,rep_env_local_index

! *****************************************************************************
!> \brief keeps replicated information about the replicas
!> \param ref_count reference count
!> \param id_nr identity number (unique or each replica_env)
!> \param nrep number of replicas
!> \param nparticle number of particles (usually atoms) in each replica
!> \param ndim = 3*nparticle
!> \param f_env_id id of the force env that will do the calculations for the
!>        replicas owned by this processor
!> \param r ,v,f: positions, velocities and forces of the replicas.
!>        the indexing is as follow (idir,iat,irep)
!> \param replica_owner which replica group number owns the replica irep
!> \param cart 2d distribution of the processors for the replicas,
!>        a column (or row if row_force was true in the rep_env_create call)
!>        work together on the same force_env (i.e. changing the
!>        row (column) you stay in the same replica), rows (columns) have
!>        different replicas
!> \param force_dim which dimension of cart works on forces together
!>        used to be hardcoded to 1. Default is still 1, will
!>        be 2 if row_force is true in the rep_env_create call.
!> \param para_env the global para env that contains all the replicas,
!>        this is just the cart as para_env
!> \param para_env_f parallel environment of the underlying force
!>        environment
!> \param inter_rep_rank mapping replica group number -> rank in para_env_inter_rep
!>        (this used to be col_rank)
!> \param para_env_inter_rep parallel environment between replica
!> \param force_rank mapping number of processor in force env -> rank in para_env_f
!>        (this used to be row_rank)
!> \param local_rep_indices indices of the local replicas, starting at 1
!> \param rep_is_local logical if specific replica is a local one.
!> \param my_rep_group which replica group number this process belongs to
!>        (this used to be just cart%mepos(2) but with transposing the cart
!>        (row_force=.true.) became cart%mepos(1), and to generalize this it
!>        is now a separate variable, so one does not need to know
!>        which way the cart is mapped.)
!> \param wf_history wavefunction history for the owned replicas
!> \param keep_wf_history if the wavefunction history for the owned replicas
!>        should be kept
!> \author fawzi
! *****************************************************************************
  TYPE replica_env_type
     INTEGER                                           :: ref_count, id_nr, f_env_id,&
                                                          nrep, ndim, nparticle,&
                                                          my_rep_group, force_dim
     REAL(kind=dp), DIMENSION(:,:), POINTER            :: r,v,f
     LOGICAL                                           :: sync_v, keep_wf_history
     CHARACTER(LEN=default_path_length)                :: original_project_name
     TYPE(qs_wf_history_p_type), DIMENSION(:), POINTER :: wf_history
     TYPE(cp_result_p_type),DIMENSION(:),POINTER       :: results
     INTEGER, DIMENSION(:), POINTER                    :: local_rep_indices
     INTEGER, DIMENSION(:), POINTER                    :: replica_owner, force_rank, &
                                                          inter_rep_rank
     LOGICAL, DIMENSION(:), POINTER                    :: rep_is_local
     TYPE(cp_para_cart_type), POINTER                  :: cart
     TYPE(cp_para_env_type), POINTER                   :: para_env, para_env_f,&
                                                          para_env_inter_rep
  END TYPE replica_env_type

! *****************************************************************************
!> \brief ****s* replica_types/replica_env_p_type *
!>
!>      to build arrays of pointers to a replica_env_type
!> \param rep_env the pointer to the replica_env
!> \author fawzi
! *****************************************************************************
  TYPE replica_env_p_type
     TYPE(replica_env_type), POINTER                   :: rep_env
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
!> \param row_force to use the new mapping to the cart with rows
!>        working on force instead of columns.
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE rep_env_create(rep_env, para_env, input, nrep, prep,&
       sync_v,keep_wf_history,row_force,error)
    TYPE(replica_env_type), POINTER          :: rep_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: input
    INTEGER                                  :: nrep, prep
    LOGICAL, INTENT(in), OPTIONAL            :: sync_v, keep_wf_history, &
                                                row_force
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rep_env_create', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_path_length)       :: input_file_path, &
                                                output_file_path
    INTEGER :: comm_cart, comm_f, comm_inter_rep, forcedim, i, i0, ierr, ip, &
      ir, irep, lp, my_prep, new_env_id, nparticle, nrep_local, stat, unit_nr
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: gridinfo
    INTEGER, DIMENSION(2)                    :: dims, pos
    LOGICAL                                  :: failure
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
       unit_nr=cp_logger_get_default_io_unit(logger)
       new_env_id=-1
       forcedim = 1
       IF (PRESENT(row_force)) THEN
          IF (row_force) forcedim = 2
       END IF
       my_prep=MIN(prep,para_env%num_pe)
       dims(3-forcedim)=MIN(para_env%num_pe/my_prep,nrep)
       dims(forcedim)=my_prep
       IF ((dims(1)*dims(2)/=para_env%num_pe).AND.(unit_nr>0)) THEN
          WRITE(unit_nr,FMT="(T2,A)") "REPLICA| WARNING: number of processors is not divisible by the number of replicas"
          WRITE(unit_nr,FMT="(T2,A,I0,A)") "REPLICA| ",para_env%num_pe-dims(1)*dims(2)," MPI process(es) will be idle"
       END IF
       CALL mp_cart_create ( comm_old=para_env%group, ndims=2, dims=dims, pos=pos, comm_cart=comm_cart)
       IF (comm_cart/=MPI_COMM_NULL) THEN
          CALL cp_cart_create(cart,comm_cart,ndims=2,owns_group=.TRUE.,error=error)
          NULLIFY(para_env_full)
          CALL cp_para_env_create(para_env_full,comm_cart,owns_group=.FALSE.,error=error)
          rdim(3-forcedim)=.FALSE.
          rdim(forcedim)=.TRUE.
          CALL mp_cart_sub( comm=comm_cart, rdim=rdim , sub_comm=comm_f)
          CALL cp_para_env_create(para_env_f,comm_f,owns_group=.TRUE.,error=error)
          rdim(3-forcedim)=.TRUE.
          rdim(forcedim)=.FALSE.
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
          WRITE(unit_nr,FMT="(T2,A,T71,I10)") "REPLICA| layout of the replica grid, number of groups ",para_env_inter_rep%num_pe
          WRITE(unit_nr,FMT="(T2,A,T71,I10)") "REPLICA| layout of the replica grid, size of each group",para_env_f%num_pe
          WRITE(unit_nr,FMT="(T2,A)",ADVANCE="NO") "REPLICA| MPI process to grid (group,rank) correspondence:"
          DO i=0,para_env%num_pe-1
             IF (MODULO(i,4)==0) WRITE(unit_nr,*)
             WRITE(unit_nr,FMT='(A3,I4,A3,I4,A1,I4,A1)',ADVANCE="NO")&
                  "  (",i," : ",gridinfo(3-forcedim,i),",",&
                  gridinfo(forcedim,i),")"
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

       rep_env%force_dim = forcedim
       rep_env%my_rep_group = cart%mepos(3-forcedim)
       ALLOCATE(rep_env%inter_rep_rank(0:para_env_inter_rep%num_pe-1),&
                rep_env%force_rank(0:para_env_f%num_pe-1),stat=stat)
       rep_env%inter_rep_rank=0
       rep_env%inter_rep_rank(rep_env%my_rep_group)=para_env_inter_rep%mepos
       CALL mp_sum(rep_env%inter_rep_rank,para_env_inter_rep%group)
       rep_env%force_rank=0
       rep_env%force_rank(cart%mepos(forcedim))=para_env_f%mepos
       CALL mp_sum(rep_env%force_rank,para_env_f%group)


       CALL section_vals_val_get(input,"GLOBAL%PROJECT_NAME",&
            c_val=input_file_path,error=error)
       rep_env%original_project_name = input_file_path
       ! By default replica_env handles files for each replica
       ! with the structure PROJECT_NAME-r-N where N is the
       ! number of the local replica..
       lp=LEN_TRIM(input_file_path)
       input_file_path(lp+1:LEN(input_file_path))="-r-"//&
            ADJUSTL(cp_to_string(rep_env%my_rep_group))
       lp=LEN_TRIM(input_file_path)
       ! Setup new project name
       CALL section_vals_val_set(input,"GLOBAL%PROJECT_NAME",&
            c_val=input_file_path,error=error)
       ! Redirect the output of each replica on a same local file
       output_file_path=input_file_path(1:lp)//".out"
       CALL section_vals_val_set(input,"GLOBAL%OUTPUT_FILE_NAME",&
            c_val=TRIM(output_file_path),error=error)

       ! Dump an input file to warm-up new force_eval structures and
       ! delete them immediately afterwards..
       input_file_path(lp+1:LEN(input_file_path))=".inp"
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

       ! Delete input files..
       IF (para_env_f%source==para_env_f%mepos) THEN
          CALL open_file(file_name=TRIM(input_file_path),file_status="OLD",&
               file_form="FORMATTED",file_action="READ",unit_number=unit_nr)
          CALL close_file(unit_number=unit_nr,file_status="DELETE")
       END IF

       IF (.NOT.failure) THEN
          rep_env%f_env_id=new_env_id
          CALL cp_get_nparticle(new_env_id,nparticle,ierr)
          CPPostcondition(ierr==0,cp_fatal_level,routineP,error,failure)
          rep_env%nparticle = nparticle
          rep_env%ndim = 3*nparticle
          ALLOCATE(rep_env%replica_owner(nrep),stat=stat)
          CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)

          i0=nrep/para_env_inter_rep%num_pe
          ir=MODULO(nrep,para_env_inter_rep%num_pe)
          DO ip=0,para_env_inter_rep%num_pe-1
             DO i=i0*ip+MIN(ip,ir)+1,i0*(ip+1)+MIN(ip+1,ir)
                rep_env%replica_owner(i)=ip
             END DO
          END DO

          nrep_local=i0
          IF (rep_env%my_rep_group<ir) nrep_local=nrep_local+1
          ALLOCATE(rep_env%local_rep_indices(nrep_local),&
                   rep_env%rep_is_local(nrep),stat=stat)
          CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
          nrep_local=0
          rep_env%rep_is_local = .FALSE.
          DO irep=1,nrep
             IF (rep_env%replica_owner(irep)==rep_env%my_rep_group) THEN
                nrep_local=nrep_local+1
                rep_env%local_rep_indices(nrep_local)=irep
                rep_env%rep_is_local(irep) = .TRUE.
             END IF
          END DO
          CPPostcondition(nrep_local==SIZE(rep_env%local_rep_indices),cp_failure_level,routineP,error,failure)

          rep_env%cart => cart
          rep_env%para_env => para_env_full
          rep_env%para_env_f => para_env_f
          rep_env%para_env_inter_rep => para_env_inter_rep

          ALLOCATE (rep_env%r(rep_env%ndim,nrep),rep_env%v(rep_env%ndim,nrep),&
                    rep_env%f(rep_env%ndim+1,nrep), stat=stat)
          CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)

          rep_env%r=0._dp
          rep_env%f=0._dp
          rep_env%v=0._dp
          CALL cp_set_vel(rep_env%f_env_id, rep_env%v(:,1), rep_env%ndim, ierr)
          CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
          DO i=1,nrep
             IF (rep_env%rep_is_local(i)) THEN
                CALL cp_get_pos(rep_env%f_env_id,rep_env%r(:,i),rep_env%ndim,ierr)
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
    CALL timeset(routineN,handle)
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
          DEALLOCATE(rep_env%rep_is_local,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          IF (ASSOCIATED(rep_env%replica_owner)) THEN
             DEALLOCATE(rep_env%replica_owner,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          DEALLOCATE(rep_env%inter_rep_rank,rep_env%force_rank,stat=stat)
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
    CALL timeset(routineN,handle)
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

    CALL timeset(routineN,handle)
    IF (ASSOCIATED(rep_env)) THEN
       WRITE(unit_nr,"('<replica_env, id_nr=',i6,' ref_count=',i6,'> sync_v=',l1,',')")&
            rep_env%id_nr,rep_env%ref_count, rep_env%sync_v
       WRITE (unit_nr,"(' f_env_id=',i6,', nrep=',i6,', nparticle=',i8,', ndim=',i8,',')")&
            rep_env%f_env_id, rep_env%nrep,rep_env%nparticle, rep_env%ndim
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
       WRITE (unit_nr,"(' force_rank=(')", advance="no")
       WRITE (unit_nr,"(10i6)") rep_env%force_rank
       WRITE (unit_nr,"(')')")
       WRITE (unit_nr,"(' inter_rep_rank=(')", advance="no")
       WRITE (unit_nr,"(10i6)") rep_env%inter_rep_rank
       WRITE (unit_nr,"(')')")
       WRITE (unit_nr,"(' force_dim=(')", advance="no")
       WRITE (unit_nr,"(10i6)") rep_env%force_dim
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
    CALL timeset(routineN,handle)
    CPPrecondition(ASSOCIATED(rep_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(rep_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(vals,2)==rep_env%nrep,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       DO irep=1,rep_env%nrep
          IF (.NOT.rep_env%rep_is_local(irep)) THEN
             vals(:,irep)=0._dp
          END IF
       END DO
       CALL mp_sum(vals,rep_env%para_env_inter_rep%group)
    END IF
    CALL timestop(handle)
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
  SUBROUTINE rep_env_sync_results(rep_env,results,error)
    TYPE(replica_env_type), POINTER          :: rep_env
    TYPE(cp_result_p_type), DIMENSION(:), &
      POINTER                                :: results
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rep_env_sync_results', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, irep, nrep, source
    LOGICAL                                  :: failure

    failure=.FALSE.
    CALL timeset(routineN,handle)
    nrep=rep_env%nrep
    CPPrecondition(ASSOCIATED(rep_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(rep_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(results)==rep_env%nrep,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       DO irep=1, nrep
!         source = 0
!         IF (rep_env%rep_is_local(irep)) source = rep_env%para_env_inter_rep%mepos
!         CALL mp_sum(source, rep_env%para_env_inter_rep%group)
! above three lines should be the same as just:
          source = rep_env%inter_rep_rank(rep_env%replica_owner(irep))
          CALL cp_results_mp_bcast(results(irep)%results, source, rep_env%para_env_inter_rep, error)
       END DO
    END IF
    CALL timestop(handle)
  END SUBROUTINE rep_env_sync_results

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

    CALL timeset(routineN,handle)
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
  FUNCTION rep_envs_get_rep_env(id_nr,ierr) RESULT(res)
    INTEGER, INTENT(in)                      :: id_nr
    INTEGER, INTENT(OUT)                     :: ierr
    TYPE(replica_env_type), POINTER          :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'rep_envs_get_rep_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

    NULLIFY(res)
    ierr = -1
    IF (module_initialized) THEN
       IF (ASSOCIATED(rep_envs)) THEN
          DO i=1,SIZE(rep_envs)
             IF (rep_envs(i)%rep_env%id_nr==id_nr) THEN
                res => rep_envs(i)%rep_env
                ierr = 0
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
       rep_env2 => rep_envs_get_rep_env(rep_env%id_nr,ierr=stat)
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
