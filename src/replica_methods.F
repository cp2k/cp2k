!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief methods to setup replicas of the same system differing only by atom
!>      positions and velocities (as used in path integral or nudged elastic
!>      band for example)
!> \par History
!>      09.2005 created [fawzi]
!> \author fawzi
! *****************************************************************************
MODULE replica_methods
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_output_handling,              ONLY: cp_add_iter_level
  USE cp_para_env,                     ONLY: cp_cart_create,&
                                             cp_para_env_create
  USE cp_para_types,                   ONLY: cp_para_cart_type,&
                                             cp_para_env_type
  USE cp_result_types,                 ONLY: cp_result_create,&
                                             cp_result_retain
  USE f77_blas
  USE f77_interface,                   ONLY: calc_force,&
                                             create_force_env,&
                                             f_env_add_defaults,&
                                             f_env_rm_defaults,&
                                             f_env_type,&
                                             get_nparticle,&
                                             get_pos,&
                                             set_vel
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_set,&
                                             use_qs_force
  USE input_section_types,             ONLY: section_type,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set,&
                                             section_vals_write
  USE iso_c_binding
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE message_passing,                 ONLY: MPI_COMM_NULL,&
                                             mp_cart_create,&
                                             mp_cart_sub,&
                                             mp_sum
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_wf_history_methods,           ONLY: wfi_create
  USE qs_wf_history_types,             ONLY: wfi_retain
  USE replica_types,                   ONLY: rep_env_release,&
                                             rep_env_sync,&
                                             rep_env_sync_results,&
                                             rep_env_write,&
                                             rep_envs_add_rep_env,&
                                             rep_envs_get_rep_env,&
                                             replica_env_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'replica_methods'
  INTEGER, SAVE, PRIVATE :: last_rep_env_id=0

  PUBLIC :: rep_env_create, rep_env_calc_e_f

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
  SUBROUTINE rep_env_create(rep_env, para_env, input, input_declaration, nrep, prep,&
       sync_v,keep_wf_history,row_force,error)
    TYPE(replica_env_type), POINTER          :: rep_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: input
    TYPE(section_type), POINTER              :: input_declaration
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
    CPPrecondition(ASSOCIATED(input_declaration),cp_failure_level,routineP,error,failure)

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
       CALL create_force_env(new_env_id,input_declaration,input_file_path,&
             output_file_path,para_env_f%group,ierr=ierr)
       CPAssert(ierr==0,cp_failure_level,routineP,error,failure)

       ! Delete input files..
       IF (para_env_f%source==para_env_f%mepos) THEN
          CALL open_file(file_name=TRIM(input_file_path),file_status="OLD",&
               file_form="FORMATTED",file_action="READ",unit_number=unit_nr)
          CALL close_file(unit_number=unit_nr,file_status="DELETE")
       END IF

       IF (.NOT.failure) THEN
          rep_env%f_env_id=new_env_id
          CALL get_nparticle(new_env_id,nparticle,ierr)
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
          CALL set_vel(rep_env%f_env_id, rep_env%v(:,1), rep_env%ndim, ierr)
          CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
          DO i=1,nrep
             IF (rep_env%rep_is_local(i)) THEN
                CALL get_pos(rep_env%f_env_id,rep_env%r(:,i),rep_env%ndim,ierr)
                CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
             END IF
          END DO
       END IF
    END IF
    IF (.NOT.failure.AND.ASSOCIATED(rep_env)) THEN
       CALL rep_envs_add_rep_env(rep_env,error=error)
       CALL rep_env_init_low(rep_env%id_nr,ierr)
       CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
  END SUBROUTINE rep_env_create


! *****************************************************************************
!> \brief finishes the low level initialization of the replica env
!> \param rep_env_id id_nr of the replica environment that should be initialized
!> \param ierr will be non zero if there is an initialization error
!> \author fawzi
! *****************************************************************************
  SUBROUTINE rep_env_init_low(rep_env_id, ierr)
    INTEGER, INTENT(in)                      :: rep_env_id
    INTEGER, INTENT(out)                     :: ierr

    CHARACTER(len=*), PARAMETER :: routineN = 'rep_env_init_low', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, in_use, stat
    LOGICAL                                  :: failure, has_unit_metric
    TYPE(cp_error_type)                      :: error
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(replica_env_type), POINTER          :: rep_env

    failure=.FALSE.
    rep_env => rep_envs_get_rep_env(rep_env_id,ierr=stat)
    CALL cp_assert(ASSOCIATED(rep_env),cp_failure_level,cp_assertion_failed,&
         routineP,"could not find rep_env with id_nr"//cp_to_string(rep_env_id),&
         failure=failure)
    NULLIFY(qs_env,dft_control)
    IF (.NOT. failure) THEN
       CALL f_env_add_defaults(f_env_id=rep_env%f_env_id,f_env=f_env,&
            new_error=error, failure=failure)
       IF (.NOT.failure) THEN
          logger => cp_error_get_logger(error)
          logger%iter_info%iteration(1)=rep_env%my_rep_group
          CALL cp_add_iter_level(iteration_info=logger%iter_info,&
               level_name="REPLICA_EVAL",error=error)
          !wf interp
          IF (rep_env%keep_wf_history) THEN
             CALL force_env_get(f_env%force_env,in_use=in_use,error=error)
             IF (in_use==use_qs_force) THEN
                CALL force_env_get(f_env%force_env,qs_env=qs_env,&
                     error=error)
                CALL get_qs_env(qs_env,dft_control=dft_control,error=error)
                ALLOCATE(rep_env%wf_history(SIZE(rep_env%local_rep_indices)),stat=stat)
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
                      CALL get_qs_env(qs_env,has_unit_metric=has_unit_metric,error=error)
                      CALL wfi_create(rep_env%wf_history(i)%wf_history,&
                           interpolation_method_nr=&
                           dft_control%qs_control%wf_interpolation_method_nr,&
                           extrapolation_order = dft_control%qs_control%wf_extrapolation_order,&
                           has_unit_metric = has_unit_metric, &
                           error=error)
                   END IF
                END DO
             ELSE
                rep_env%keep_wf_history=.FALSE.
             END IF
          END IF
          ALLOCATE(rep_env%results(rep_env%nrep),stat=stat)
          CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
          DO i=1, rep_env%nrep
             NULLIFY(rep_env%results(i)%results)
             IF (i==1) THEN
                CALL force_env_get(f_env%force_env,&
                     results=rep_env%results(i)%results,&
                     error=error)
                CALL cp_result_retain(rep_env%results(i)%results,error)
             ELSE
                CALL cp_result_create(rep_env%results(i)%results,error)
             END IF
          END DO
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
       CALL rep_env_calc_e_f_low(rep_env%id_nr,my_calc_f,ierr)
       CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(handle)
  END SUBROUTINE rep_env_calc_e_f


! *****************************************************************************
!> \brief calculates energy and force, internal private method
!> \param rep_env_id the id if the replica environment in which energy and
!>        forces have to be evaluated
!> \param calc_f if nonzero calculates also the forces along with the
!>        energy
!> \param ierr if an error happens this will be nonzero
!> \note
!>      low level wrapper to export this function in f77_int_low and work
!>      around the handling of circular dependecies in fortran
!> \author fawzi
! *****************************************************************************
  RECURSIVE SUBROUTINE rep_env_calc_e_f_low(rep_env_id,calc_f,ierr)
    INTEGER, INTENT(in)                      :: rep_env_id, calc_f
    INTEGER, INTENT(out)                     :: ierr

    CHARACTER(len=*), PARAMETER :: routineN = 'rep_env_calc_e_f_low', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(cp_error_type)                      :: new_error
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(replica_env_type), POINTER          :: rep_env

    failure=.FALSE.
    rep_env => rep_envs_get_rep_env(rep_env_id, ierr)
    IF (ASSOCIATED(rep_env)) THEN
       CALL f_env_add_defaults(f_env_id=rep_env%f_env_id,f_env=f_env,&
            new_error=new_error, failure=failure)
       CALL rep_env_calc_e_f_int(rep_env,calc_f/=0,new_error)
       CALL f_env_rm_defaults(f_env,new_error,ierr)
    ELSE
       ierr=111
    END IF
  END SUBROUTINE rep_env_calc_e_f_low


! *****************************************************************************
!> \brief calculates energy and force, internal private method
!> \param rep_env the replica env to update
!> \param if the force should be calculated as well (defaults to true)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      this is the where the real work is done
!> \author fawzi
! *****************************************************************************
  SUBROUTINE rep_env_calc_e_f_int(rep_env,calc_f,error)
    TYPE(replica_env_type), POINTER          :: rep_env
    LOGICAL, OPTIONAL                        :: calc_f
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rep_env_calc_e_f_int', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ierr, irep, md_iter, &
                                                my_calc_f, ndim
    LOGICAL                                  :: failure
    TYPE(cp_error_type)                      :: new_error
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(qs_environment_type), POINTER       :: qs_env

    failure=.FALSE.
    NULLIFY(f_env,qs_env)
    CPPrecondition(ASSOCIATED(rep_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(rep_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       my_calc_f=3*rep_env%nparticle
       IF (PRESENT(calc_f)) THEN
          IF (.NOT.calc_f) my_calc_f=0
       END IF

       CALL f_env_add_defaults(f_env_id=rep_env%f_env_id,f_env=f_env,&
            new_error=new_error, failure=failure)
       logger => cp_error_get_logger(new_error)
       !     md_iter=logger%iter_info%iteration(2)+1
       md_iter=logger%iter_info%iteration(2)
       CALL f_env_rm_defaults(f_env,new_error,ierr)
       CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
       DO i=1,SIZE(rep_env%local_rep_indices)
          irep = rep_env%local_rep_indices(i)
          ndim = 3*rep_env%nparticle
          IF (rep_env%sync_v) THEN
             CALL set_vel(rep_env%f_env_id,rep_env%v(:,irep),ndim,ierr)
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

          CALL f_env_add_defaults(f_env_id=rep_env%f_env_id,f_env=f_env,&
               new_error=new_error, failure=failure)
          CALL force_env_set(force_env=f_env%force_env,&
               results=rep_env%results(irep)%results,&
               error=error)
          CALL f_env_rm_defaults(f_env,new_error,ierr)
          CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
          CALL calc_force(rep_env%f_env_id,rep_env%r(:,irep),ndim,&
                             rep_env%f(ndim+1,irep),rep_env%f(:ndim,irep),&
                             my_calc_f,ierr)
          CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
          IF (failure) EXIT
       END DO
       IF (.NOT.failure) THEN
          CALL rep_env_sync(rep_env,rep_env%f,error=error)
          CALL rep_env_sync_results(rep_env,rep_env%results,error)
       END IF
       CALL cp_error_synchronize_error(error,rep_env%para_env)
    END IF

  END SUBROUTINE rep_env_calc_e_f_int


! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE rep_env_test(para_env,input,input_declaration,error)
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: input
    TYPE(section_type), POINTER              :: input_declaration
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rep_env_test', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: project_name
    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: err, max_err
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: e1, e2, e3
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: f1, f2, f3
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(replica_env_type), POINTER          :: rep_env

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(para_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(input),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(input_declaration),cp_failure_level,routineP,error,failure)
    CPPrecondition(para_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(input%ref_count>0,cp_failure_level,routineP,error,failure)
    max_err=0._dp
    IF (.NOT. failure) THEN
       logger => cp_error_get_logger(error)
       NULLIFY(rep_env,f1,f2,f3)
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
            input_declaration=input_declaration,nrep=1,prep=1,error=error)
       CALL rep_env_write(rep_env,cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
            error=error)
       IF (ASSOCIATED(rep_env)) THEN
          CALL rep_env_calc_e_f(rep_env,error=error)
          ALLOCATE(f1(rep_env%ndim,rep_env%nrep),e1(rep_env%nrep),stat=stat)
          CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
          f1=rep_env%f(:rep_env%ndim,:)
          e1(:)=rep_env%f(rep_env%ndim+1,:)
       END IF
       CALL rep_env_release(rep_env,error=error)

       CALL section_vals_val_set(input,"GLOBAL%PROJECT_NAME",&
            c_val=TRIM(project_name)//"t12",error=error)
       CALL rep_env_create(rep_env, para_env=para_env, input=input,&
             input_declaration=input_declaration,nrep=1,prep=2,error=error)
       CALL rep_env_write(rep_env,cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
            error=error)
       IF (ASSOCIATED(rep_env)) THEN
          CALL rep_env_calc_e_f(rep_env,error=error)
          ALLOCATE(f2(rep_env%ndim,rep_env%nrep),e2(rep_env%nrep),stat=stat)
          CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
          f2=rep_env%f(:rep_env%ndim,:)
          e2(:)=rep_env%f(rep_env%ndim+1,:)
       END IF
       CALL rep_env_release(rep_env,error=error)

       CALL section_vals_val_set(input,"GLOBAL%PROJECT_NAME",&
            c_val=TRIM(project_name)//"t21",error=error)
       CALL rep_env_create(rep_env, para_env=para_env, input=input,&
             input_declaration=input_declaration,nrep=2,prep=1,error=error)
       CALL rep_env_write(rep_env,cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
            error=error)
       IF (ASSOCIATED(rep_env)) THEN
          rep_env%r(:,2)=rep_env%r(:,1)+0.5_dp
          CALL rep_env_calc_e_f(rep_env,error=error)
          ALLOCATE(f3(rep_env%ndim,rep_env%nrep),e3(rep_env%nrep),stat=stat)
          CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
          f3=rep_env%f(:rep_env%ndim,:)
          e3(:)=rep_env%f(rep_env%ndim+1,:)
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
          DEALLOCATE(f3,e3,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
       IF (ASSOCIATED(f1)) THEN
          DEALLOCATE(f1,e1,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
       IF (ASSOCIATED(f2)) THEN
          DEALLOCATE(f2,e2,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
       CALL cp_log(logger,cp_note_level+1,routineP,"CUMULATIVE ERROR"//&
            cp_to_string(max_err), local=.FALSE.)
       CPAssert(max_err<1.e-6,cp_failure_level,routineP,error,failure)
    END IF
  END SUBROUTINE rep_env_test

END MODULE replica_methods
