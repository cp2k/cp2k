!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
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
  USE cp_output_handling,              ONLY: cp_rm_iter_level
  USE cp_para_env,                     ONLY: cp_cart_release,&
                                             cp_cart_write,&
                                             cp_para_env_release,&
                                             cp_para_env_write
  USE cp_para_types,                   ONLY: cp_para_cart_type,&
                                             cp_para_env_type
  USE cp_result_methods,               ONLY: cp_results_mp_bcast
  USE cp_result_types,                 ONLY: cp_result_p_type,&
                                             cp_result_release
  USE f77_interface,                   ONLY: destroy_force_env,&
                                             f_env_add_defaults,&
                                             f_env_rm_defaults,&
                                             f_env_type
  USE kinds,                           ONLY: default_path_length,&
                                             dp
  USE message_passing,                 ONLY: mp_sum
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

  PUBLIC :: replica_env_type
  PUBLIC :: rep_env_release
  PUBLIC :: rep_env_sync, rep_env_sync_results, rep_envs_add_rep_env
  PUBLIC :: rep_envs_get_rep_env

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
!> \brief releases the given replica environment
!> \param rep_env the replica environment to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
!> \note
!>      here and not in replica_types to allow the use of replica_env_type
!>      in a force_env (call to destroy_force_env gives circular dep)
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
          CALL rep_env_destroy_low(rep_env%id_nr,ierr)
          IF (rep_env%f_env_id>0) THEN
             CALL destroy_force_env(rep_env%f_env_id,ierr)
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
!> \brief initializes the destruction of the replica_env
!> \param rep_env_id id_nr of the replica environment that should be initialized
!> \param ierr will be non zero if there is an initialization error
!> \author fawzi
! *****************************************************************************
  SUBROUTINE rep_env_destroy_low(rep_env_id, ierr)
    INTEGER, INTENT(in)                      :: rep_env_id
    INTEGER, INTENT(out)                     :: ierr

    CHARACTER(len=*), PARAMETER :: routineN = 'rep_env_destroy_low', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure
    TYPE(cp_error_type)                      :: error
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(replica_env_type), POINTER          :: rep_env

    failure=.FALSE.
    rep_env => rep_envs_get_rep_env(rep_env_id,ierr=stat)
    CALL cp_assert(ASSOCIATED(rep_env),cp_failure_level,cp_assertion_failed,&
         routineP,"could not find rep_env with id_nr"//cp_to_string(rep_env_id),&
         failure=failure)
    CALL f_env_add_defaults(f_env_id=rep_env%f_env_id,f_env=f_env,&
         new_error=error, failure=failure)
    logger => cp_error_get_logger(error)
    CALL cp_rm_iter_level(iteration_info=logger%iter_info,&
         level_name="REPLICA_EVAL",error=error)
    CALL f_env_rm_defaults(f_env,error,ierr)
    CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE rep_env_destroy_low

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
    CPPrecondition(rep_env%ref_count>0,cp_failure_level,routineP,error,failure)
    rep_env%ref_count=rep_env%ref_count+1
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
!> \author fawzi
!> \note
!>      could be optimized: bcast in inter_rep, all2all or shift vs sum
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
    DO irep=1,rep_env%nrep
       IF (.NOT.rep_env%rep_is_local(irep)) THEN
          vals(:,irep)=0._dp
       END IF
    END DO
    CALL mp_sum(vals,rep_env%para_env_inter_rep%group)
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
    DO irep=1, nrep
!         source = 0
!         IF (rep_env%rep_is_local(irep)) source = rep_env%para_env_inter_rep%mepos
!         CALL mp_sum(source, rep_env%para_env_inter_rep%group)
! above three lines should be the same as just:
       source = rep_env%inter_rep_rank(rep_env%replica_owner(irep))
       CALL cp_results_mp_bcast(results(irep)%results, source, rep_env%para_env_inter_rep, error)
    END DO
    CALL timestop(handle)
  END SUBROUTINE rep_env_sync_results


! *****************************************************************************
!> \brief returns the replica environment with the given id_nr
!> \param id_nr the id_nr of the requested rep_envs
!> \param ierr ...
!> \retval res ...
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
!> \retval res ...
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
    CALL timestop(handle)
  END FUNCTION rep_env_local_index

END MODULE replica_types
