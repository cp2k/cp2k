!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2002  CP2K developers group
!-----------------------------------------------------------------------------!
!!****** cp2k/message_passing [1.0] *
!!
!!   NAME
!!     message_passing
!!
!!   FUNCTION
!!     Interface to the message passing library MPI
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     JGH (02-Jan-2001): New error handling
!!                        Performance tools
!!     JGH (14-Jan-2001): New routines mp_comm_compare, mp_cart_coords, 
!!                                     mp_rank_compare, mp_alltoall
!!     JGH (06-Feb-2001): New routines mp_comm_free
!!     JGH (22-Mar-2001): New routines mp_comm_dup
!!     fawzi (04-NOV-2004): storable performance info (for f77 interface)
!!
!!   SOURCE
!******************************************************************************

MODULE message_passing
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_cputime
  USE parallel_include

  IMPLICIT NONE

  PRIVATE

#ifdef __parallel
  LOGICAL, PARAMETER :: cp2k_is_parallel=.TRUE.
#else
  LOGICAL, PARAMETER :: cp2k_is_parallel=.FALSE.
#endif

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN="message_passing"
  ! init and error
  PUBLIC :: mp_world_init, mp_world_finalize
  PUBLIC :: mp_stop, mp_abort

  ! performance gathering
  PUBLIC :: mp_perf_env_type, mp_perf_env_p_type,mp_perf_type
  PUBLIC :: mp_perf_env_create, mp_perf_env_retain, mp_perf_env_release,&
       mp_perf_env_describe
  PUBLIC :: add_mp_perf_env, rm_mp_perf_env, get_mp_perf_env, describe_mp_perf_env

  ! informational / generation of sub comms
  PUBLIC :: mp_environ, mp_comm_compare, mp_cart_coords, mp_rank_compare
  PUBLIC :: mp_cart_create, mp_dims_create, mp_cart_rank, mp_cart_sub, mp_comm_free
  PUBLIC :: mp_comm_dup, mp_comm_split
  PUBLIC :: cp2k_is_parallel

  ! message passing

  PUBLIC :: mp_bcast, mp_sum, mp_max, mp_maxloc, mp_minloc, mp_min, mp_sync
  PUBLIC :: mp_gather, mp_alltoall, mp_sendrecv, mp_allgather
  PUBLIC :: mp_cart_shift, mp_isend, mp_irecv
  PUBLIC :: mp_sum_scatter, mp_shift, mp_isendrecv, mp_wait, mp_waitall

  ! currently not used
  PUBLIC :: mp_group

  ! default communicators
  PUBLIC :: MPI_COMM_SELF, MPI_COMM_WORLD
#ifndef __parallel
  INTEGER, PARAMETER :: MPI_COMM_SELF=0,MPI_COMM_WORLD=0
#endif

  !
  ! interfaces to deal easily with scalars / vectors / matrice / ...
  ! of the different types (integers, doubles, logicals, characters)
  !
  INTERFACE mp_minloc
     MODULE PROCEDURE mp_minloc_rv
  END INTERFACE

  INTERFACE mp_maxloc
     MODULE PROCEDURE mp_maxloc_i
  END INTERFACE

  INTERFACE mp_shift
     MODULE PROCEDURE mp_shift_im, mp_shift_r1, mp_shift_rm
  END INTERFACE
  
  INTERFACE mp_bcast
     MODULE PROCEDURE mp_bcast_i1, mp_bcast_r1, mp_bcast_c1, mp_bcast_z, &
          mp_bcast_iv, mp_bcast_rv, mp_bcast_cv, mp_bcast_l, mp_bcast_rm, &
          mp_bcast_cm, mp_bcast_im, mp_bcast_rm3
  END INTERFACE

  INTERFACE mp_sum
     MODULE PROCEDURE mp_sum_i1, mp_sum_r1, mp_sum_c1, mp_sum_iv, &
          mp_sum_rv, mp_sum_cv, mp_sum_im, mp_sum_rm, mp_sum_cm, &
          mp_sum_im3, mp_sum_rm3, mp_sum_cm3, mp_sum_root_rv, &
          mp_sum_root_rm, mp_sum_root_rm3
  END INTERFACE

  INTERFACE mp_max
     MODULE PROCEDURE mp_max_i, mp_max_iv, mp_max_r, mp_max_rv
  END INTERFACE

  INTERFACE mp_min
     MODULE PROCEDURE mp_min_i, mp_min_iv, mp_min_r, mp_min_rv
  END INTERFACE

  INTERFACE mp_gather
     MODULE PROCEDURE mp_gather_i, mp_gather_iv, mp_gather_r, mp_gather_rv
  END INTERFACE

  INTERFACE mp_allgather
     MODULE PROCEDURE mp_allgather_i1, mp_allgatherv_rv
  END INTERFACE

  INTERFACE mp_sum_scatter
     MODULE PROCEDURE mp_sum_scatter_rv
  END INTERFACE

  INTERFACE mp_environ
     MODULE PROCEDURE mp_environ_l, mp_environ_c
  END INTERFACE

  INTERFACE mp_alltoall
     MODULE PROCEDURE mp_alltoall_c22v, mp_alltoall_i, mp_alltoall_i11v, mp_alltoall_r11v
  END INTERFACE

  INTERFACE mp_sendrecv
     MODULE PROCEDURE mp_sendrecv_rm3, mp_sendrecv_rm2, mp_sendrecv_rv
  END INTERFACE

  INTERFACE mp_isendrecv
     MODULE PROCEDURE mp_isendrecv_rm2,mp_isendrecv_rv
  END INTERFACE

  INTERFACE mp_isend
     MODULE PROCEDURE mp_isend_rm2
  END INTERFACE

  INTERFACE mp_irecv
     MODULE PROCEDURE mp_irecv_rm2
  END INTERFACE

  ! type internally used to store message passing performance indicators
  TYPE mp_perf_type
    CHARACTER ( LEN = 20 ) :: name
    INTEGER :: count
    REAL (KIND=dp) :: msg_size
    REAL (KIND=dp) :: time
  END TYPE mp_perf_type

  INTEGER, PARAMETER :: MAX_PERF = 12

  TYPE mp_perf_env_type
     !private
     INTEGER :: ref_count, id_nr
     TYPE ( mp_perf_type ), DIMENSION ( MAX_PERF ) :: mp_perfs
  END TYPE mp_perf_env_type

  TYPE mp_perf_env_p_type
     TYPE(mp_perf_env_type), POINTER :: mp_perf_env
  END TYPE mp_perf_env_p_type
     
  ! introduce a stack of mp_perfs, first index is the stack pointer, for convience is replacing
  INTEGER, PARAMETER :: max_stack_size = 10
  INTEGER            :: stack_pointer = 0
  ! target attribute needed as a hack around ifc 7.1 bug 
  TYPE ( mp_perf_env_p_type ), DIMENSION (max_stack_size), TARGET :: mp_perf_stack

  CHARACTER ( LEN = 20 ), PARAMETER :: sname ( MAX_PERF ) =  &
   (/"MP_Group            ", "MP_Bcast            ", "MP_Allreduce        ", &
     "MP_Gather           ", "MP_Sync             ", "MP_Alltoall         ", &
     "MP_SendRecv         ", "MP_ISendRecv        ", "MP_Wait             ", &
     "MP_comm_split       ", "MP_ISend            ", "MP_IRecv            "/)
  REAL (KIND=dp) :: t_start, t_end

  ! we make some assumptions on the length of INTEGERS, REALS and LOGICALS
  INTEGER, PARAMETER :: intlen=BIT_SIZE ( 0 ) / 8, reallen=8, loglen=BIT_SIZE ( 0 ) / 8, charlen=1
  INTEGER, SAVE, PRIVATE :: last_mp_perf_env_id=0

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************

!!****f* message_passing/mp_world_init *
!!
!!   NAME
!!     mp_world_init
!!
!!   FUNCTION
!!     initializes the system default communicator
!!
!!   NOTES
!!     should only be called once 
!!
!!   INPUTS
!!    -mp_comm [output] : handle of the default communicator
!!
!!   MODIFICATION HISTORY
!!     2.2004 created [Joost VandeVondele ]
!!
!!   SOURCE
!!*** **********************************************************************
SUBROUTINE mp_world_init(mp_comm)
  INTEGER, INTENT(OUT) :: mp_comm
#if defined(__parallel)
  INTEGER :: ierr
  CALL mpi_init ( ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_init @ mp_world_init" )
  CALL mpi_errhandler_set ( MPI_COMM_WORLD, MPI_ERRORS_RETURN, ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_errhandler_set @ mp_world_init" )
  mp_comm = MPI_COMM_WORLD
#else
  mp_comm = 0
#endif
END SUBROUTINE mp_world_init

!!****f* message_passing/mp_world_finalize *
!!
!!   NAME
!!     mp_world_finalize
!!
!!   FUNCTION
!!     finalizes the system default communicator
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     2.2004 created [Joost VandeVondele]
!!
!!   SOURCE
!!*** **********************************************************************
SUBROUTINE mp_world_finalize()
#if defined(__parallel)
  INTEGER :: ierr
  CALL mpi_barrier ( MPI_COMM_WORLD,ierr ) ! call mpi directly to avoid 0 stack pointer
  IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_barrier @ mp_world_finalize" )
  CALL mpi_finalize ( ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_finalize @ mp_world_finalize" )
#endif
END SUBROUTINE

!******************************************************************************
! all the following routines should work for a given communicator, not MPI_WORLD
!******************************************************************************

!!****f* message_passing/mp_perf_start/mp_perf_stop *
!!
!!   NAME
!!     mp_perf_start  / mp_perf_stop
!!
!!   FUNCTION
!!     start and stop the performance indicators
!!     for every call to start there has to be (exactly) one call to stop
!!
!!   NOTES
!!     can be used to measure performance of a sub-part of a program. 
!!     timings measured here will not show up in the outer start/stops
!!     Doesn't need a fresh communicator
!!
!!   INPUTS
!!     - mp_comm : the mpi communicator
!!     - scr     : output unit
!!
!!   MODIFICATION HISTORY
!!     2.2004 created [Joost VandeVondele]
!!
!!   SOURCE
!!*** **********************************************************************
SUBROUTINE add_mp_perf_env(perf_env)
    TYPE(mp_perf_env_type), OPTIONAL, &
      POINTER                                :: perf_env

  stack_pointer = stack_pointer + 1
  IF (stack_pointer > max_stack_size) THEN
     CALL mp_stop ( stack_pointer, "stack_pointer too large : message_passing @ add_mp_perf_env" )   
  ENDIF
  NULLIFY(mp_perf_stack(stack_pointer)%mp_perf_env)
  IF (PRESENT(perf_env)) THEN
     mp_perf_stack(stack_pointer)%mp_perf_env => perf_env
     IF (ASSOCIATED(perf_env)) CALL mp_perf_env_retain(perf_env)
  END IF
  IF (.NOT.ASSOCIATED(mp_perf_stack(stack_pointer)%mp_perf_env)) THEN
     CALL mp_perf_env_create(mp_perf_stack(stack_pointer)%mp_perf_env)
  END IF
END SUBROUTINE add_mp_perf_env

SUBROUTINE mp_perf_env_create(perf_env)
    TYPE(mp_perf_env_type), OPTIONAL, &
      POINTER                                :: perf_env

    INTEGER                                  :: i, stat

    NULLIFY(perf_env)
    ALLOCATE(perf_env,stat=stat)
    IF (stat/=0) THEN
       CALL mp_stop ( stat, "stack_pointer too large : message_passing @ mp_perf_create" )   
    ENDIF
    last_mp_perf_env_id=last_mp_perf_env_id+1
    perf_env%id_nr=last_mp_perf_env_id
    perf_env%ref_count=1
    DO i = 1, MAX_PERF
       perf_env%mp_perfs(i) % name = sname ( i )
       perf_env%mp_perfs(i) % count = 0
       perf_env%mp_perfs(i) % msg_size = 0.0_dp
       perf_env%mp_perfs(i) % time = 0.0_dp
    END DO

END SUBROUTINE mp_perf_env_create

SUBROUTINE mp_perf_env_release(perf_env)
    TYPE(mp_perf_env_type), POINTER          :: perf_env

    INTEGER                                  :: stat

  IF (ASSOCIATED(perf_env)) THEN
     IF (perf_env%ref_count<1) THEN
        CALL mp_stop(perf_env%ref_count,&
             "invalid ref_count: message_passing @ mp_perf_env_release")
     END IF
     perf_env%ref_count=perf_env%ref_count-1
     IF (perf_env%ref_count==0) THEN
        DEALLOCATE(perf_env,stat=stat)
        IF (stat/=0) THEN
           CALL mp_stop(stat,&
                "deallocation error: message_passing @ mp_perf_env_release")
        END IF
     END IF
  END IF
  NULLIFY(perf_env)
END SUBROUTINE mp_perf_env_release

SUBROUTINE mp_perf_env_retain(perf_env)
    TYPE(mp_perf_env_type), POINTER          :: perf_env

  IF (.NOT.ASSOCIATED(perf_env)) THEN
     CALL mp_stop(-1,&
          "unassociated perf_env: message_passing @ mp_perf_env_retain")
  END IF
  IF (perf_env%ref_count<1) THEN
     CALL mp_stop(perf_env%ref_count,&
          "invalid ref_count: message_passing @ mp_perf_env_retain")
  END IF
  perf_env%ref_count=perf_env%ref_count+1
END SUBROUTINE mp_perf_env_retain

!.. reports the performance counters for the MPI run
SUBROUTINE mp_perf_env_describe ( perf_env, mp_comm, scr )
  TYPE(mp_perf_env_type), POINTER :: perf_env
    INTEGER, INTENT(IN)                      :: mp_comm, scr

    INTEGER                                  :: i, ierr, iw, numtask, taskid
    REAL(KIND=dp)                                :: per, vol

  IF (.not.ASSOCIATED(perf_env)) THEN
     CALL mp_stop ( 1, "unassociated perf_env : message_passing @ mp_perf_env_describe" )   
  ENDIF
  IF (perf_env%ref_count<1) THEN
     CALL mp_stop ( perf_env%ref_count, "invalid perf_env%ref_count : message_passing @ mp_perf_env_describe" )   
  ENDIF
#if defined(__parallel)
  CALL mp_environ ( numtask, taskid, mp_comm )
  iw = scr
  IF ( taskid == 0 ) THEN
    WRITE ( iw, '( /, 1X, 79("-") )' )
    WRITE ( iw, '( " -", 77X, "-" )' )
    WRITE ( iw, '( " -", 24X, A, 24X, "-" )' ) ' MESSAGE PASSING PERFORMANCE '
    WRITE ( iw, '( " -", 77X, "-" )' )
    WRITE ( iw, '( 1X, 79("-"), / )' )
    WRITE ( iw, '( A, A, A )' ) ' ROUTINE', '             CALLS ', &
            ' TOT TIME [s]  AVE VOLUME [Bytes]  PERFORMANCE [MB/s]'
    DO i = 1, MAX_PERF

      IF ( perf_env%mp_perfs( i ) % count > 0 ) THEN
        vol = perf_env%mp_perfs( i ) % msg_size / REAL ( perf_env%mp_perfs( i ) % count,KIND=dp)
        IF ( perf_env%mp_perfs( i ) % time > 0.0_dp ) THEN
          per = perf_env%mp_perfs( i ) % msg_size / perf_env%mp_perfs (i) % time * 1.e-6_dp
        ELSE
          per = 0.0_dp
        ENDIF
        IF ( vol < 1.0_dp ) THEN
          WRITE ( iw, '(1X,A15,T17,I10,T27,F14.3)' ) &
            ADJUSTL ( perf_env%mp_perfs( i ) % name ), perf_env%mp_perfs( i ) % count, &
            perf_env%mp_perfs( i ) % time
        ELSE
          WRITE ( iw, '(1X,A15,T17,I10,T27,F14.3,T50,F11.0,T69,F12.2)' ) &
            ADJUSTL ( perf_env%mp_perfs( i ) % name ), perf_env%mp_perfs( i ) % count, &
            perf_env%mp_perfs( i ) % time, vol, per
        END IF
      ENDIF

    END DO
    WRITE ( iw, '( 1X, 79("-"), / )' )
  END IF
#endif
END SUBROUTINE mp_perf_env_describe

SUBROUTINE rm_mp_perf_env ()
  IF (stack_pointer<1) THEN
     CALL mp_stop ( stack_pointer, "no perf_env in the stack : message_passing @ rm_mp_perf_env" )   
  ENDIF
  CALL mp_perf_env_release(mp_perf_stack(stack_pointer)%mp_perf_env)
  stack_pointer = stack_pointer - 1
END SUBROUTINE rm_mp_perf_env

FUNCTION get_mp_perf_env () RESULT(res)
    TYPE(mp_perf_env_type), POINTER          :: res

  IF (stack_pointer<1) THEN
     CALL mp_stop ( stack_pointer, "no perf_env in the stack : message_passing @ get_mp_perf_env" )   
  ENDIF
  res => mp_perf_stack(stack_pointer)%mp_perf_env

END FUNCTION get_mp_perf_env

SUBROUTINE describe_mp_perf_env(mp_comm, scr)
    INTEGER, INTENT(in)                      :: mp_comm, scr

    TYPE(mp_perf_env_type), POINTER          :: perf_env

  perf_env => get_mp_perf_env()
  CALL mp_perf_env_describe(perf_env,mp_comm,scr)
END SUBROUTINE describe_mp_perf_env

!!****f* message_passing/add_perf *
!!
!!   NAME
!!     add_perf
!!   FUNCTION
!!     adds the performance informations of one call
!!   NOTES
!!     -
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** *********************************************************************
SUBROUTINE add_perf(perf_id,count,time,msg_size)
  INTEGER, INTENT(in) :: perf_id
  INTEGER, INTENT(in), OPTIONAL :: count,msg_size
  REAL(KIND=dp), INTENT(in), OPTIONAL :: time
  TYPE(mp_perf_type), POINTER :: mp_perf

#if defined(__parallel)
  
  mp_perf => mp_perf_stack (stack_pointer)%mp_perf_env%mp_perfs( perf_id )
  IF (PRESENT(count)) THEN
      mp_perf%count = mp_perf%count + count
  END IF
  IF (PRESENT(time)) THEN
     mp_perf%time = mp_perf%time + time
  END IF
  IF (PRESENT(msg_size)) THEN
     mp_perf%msg_size = mp_perf%msg_size+REAL(msg_size,dp)
  END IF
#endif

END SUBROUTINE add_perf

!!****f* message_passing/mp_stop / mp_abort *
!!
!!   NAME
!!      mp_stop / mp_abort
!!   FUNCTION
!!      stops all tasks of MPI_COMM_WORLD
!!   NOTES
!!      we should learn to error out more gracefully (i.e. for farming) 
!!   INPUTS
!!    -
!!    -
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!!*** **********************************************************************
SUBROUTINE mp_stop ( ierr, prg_code )

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: ierr
  CHARACTER ( LEN=* ) :: prg_code

! Locals
  INTEGER :: numtask, taskid, groupid, len
#if defined(__parallel)
  CHARACTER ( LEN = MPI_MAX_ERROR_STRING ) :: error_string
  CALL mp_environ ( numtask, taskid, MPI_COMM_WORLD )
#else
  numtask=1
  taskid=0
#endif
  WRITE(*,'(/,A,T71,I10)') ' CP2K| Stopped by processor number',taskid
  WRITE(*,'(A," ",A)') ' CP2K| ', TRIM(prg_code)
  WRITE(*,'(A,T71,I10,/)') ' CP2K| Error number was ',ierr

#if defined(__parallel)
  CALL mpi_error_string ( ierr, error_string, len )
  WRITE(*,'(A,A)') ' CP2K| ', ADJUSTR ( error_string )
  CALL mpi_abort ( MPI_COMM_WORLD, prg_code )
#endif

  STOP "mp_stop"


END SUBROUTINE mp_stop

!******************************************************************************

!..mp_abort
SUBROUTINE mp_abort ( )

    INTEGER                                  :: numtask, taskid
#if defined(__parallel)
  CALL mp_environ ( numtask, taskid, MPI_COMM_WORLD )
#else
  numtask=1
  taskid=0
#endif
  WRITE(*,'(/,A,T71,I10)') ' CP2K| Stopped by process number',taskid
  WRITE(*,'(A,/)')         ' CP2K| Abnormal program termination'

!FM  CALL m_abort() ! uncomment if you want nice core dumps
#if defined(__parallel)
  CALL mpi_abort ( MPI_COMM_WORLD, 0 )
#else
  STOP
#endif

END SUBROUTINE mp_abort



!!****f* message_passing/mp_sync *
!!
!!   NAME
!!      mp_sync
!!   FUNCTION
!!      synchronizes with a barrier a given group of mpi tasks
!!   NOTES
!!
!!   INPUTS
!!    - group : mpi communicator
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!!*** **********************************************************************
SUBROUTINE mp_sync ( group )

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: group

! Locals
  INTEGER :: ierr

#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_barrier ( group, ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_barrier @ mp_sync" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=4,count=1,time=t_end-t_start)
#endif

END SUBROUTINE mp_sync


!!****f* message_passing/mp_environ *
!!
!!   NAME
!!      mp_environ_l / mp_environ_c
!!   FUNCTION
!!      returns number of tasks and task id for a given mpi group
!!      simple and cartesian version
!!   NOTES
!!        ..mp_world_setup is gone, use mp_environ instead (i.e. give a groupid explicitly)
!!   INPUTS
!!    - groupid : mpi communicator
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!!*** **********************************************************************
SUBROUTINE mp_environ_l ( numtask, taskid, groupid )

    INTEGER, INTENT(OUT)                     :: numtask, taskid
    INTEGER, INTENT(IN)                      :: groupid

    INTEGER                                  :: ierr

  ierr = 0
  numtask = 1
  taskid = 0
#if defined(__parallel)
  CALL mpi_comm_rank ( groupid, taskid, ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_rank @ mp_environ_l" )

  CALL mpi_comm_size ( groupid, numtask, ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_size @ mp_environ_l" )
#endif

END SUBROUTINE mp_environ_l

SUBROUTINE mp_environ_c ( numtask, dims, task_coor, groupid )


    INTEGER, INTENT(OUT)                     :: numtask, dims( 2 ), &
                                                task_coor( 2 )
    INTEGER, INTENT(IN)                      :: groupid

    INTEGER                                  :: ierr
    LOGICAL                                  :: periods( 2 )

  ierr = 0
  numtask = 1
  task_coor = 0
  dims = 1
#if defined(__parallel)
  CALL mpi_comm_size ( groupid, numtask, ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_size @ mp_environ_c" )

  CALL mpi_cart_get ( groupid, 2, dims, periods, task_coor, ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_cart_get @ mp_environ_c" )
#endif

END SUBROUTINE mp_environ_c

!******************************************************************************


!..mp_cart_create
SUBROUTINE mp_cart_create ( comm_old, ndims, dims, pos, comm_cart )


    INTEGER, INTENT(IN)                      :: comm_old, ndims
    INTEGER, INTENT(INOUT)                   :: dims( : )
    INTEGER, INTENT(OUT)                     :: pos( : ), comm_cart

    INTEGER                                  :: ierr, nodes
    LOGICAL                                  :: period(1:ndims), reorder

  ierr = 0
  pos ( 1:ndims ) = 1
  comm_cart = comm_old
#if defined(__parallel)

  t_start = m_cputime ( )
  CALL mpi_comm_size ( comm_old, nodes, ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_size @ mp_cart_create" )

  CALL mpi_dims_create ( nodes, ndims, dims, ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_dims_create @ mp_cart_create" )

  reorder = .TRUE.
  period = .TRUE.
  CALL mpi_cart_create ( comm_old, ndims, dims, period, reorder, comm_cart, &
       ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_cart_create @ mp_cart_create" )

  CALL mpi_cart_get ( comm_cart, ndims, dims, period, pos, ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_cart_get @ mp_cart_create" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=1,count=1)
  CALL add_perf(perf_id=2,time=t_end-t_start)
#else
  dims = 1
  comm_cart = 0
#endif
  
END SUBROUTINE mp_cart_create

!******************************************************************************

!..mp_cart_coords 
SUBROUTINE mp_cart_coords ( comm, rank, coords)


    INTEGER, INTENT(IN)                      :: comm, rank
    INTEGER, DIMENSION(:), INTENT(OUT)       :: coords

    INTEGER                                  :: ierr, m

  ierr = 0
  m = SIZE ( coords )
#if defined(__parallel)
  CALL mpi_cart_coords ( comm, rank, m, coords, ierr )
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_cart_coords @ mp_cart_coords" )
#else
  coords = 0
#endif

END SUBROUTINE mp_cart_coords

!******************************************************************************

!..mp_cart_shift
SUBROUTINE mp_cart_shift ( comm, dir, disp, source, dest )


    INTEGER, INTENT(IN)                      :: comm, dir, disp
    INTEGER, INTENT(OUT)                     :: source, dest

    INTEGER                                  :: ierr

  ierr = 0
#if defined(__parallel)
  CALL mpi_cart_shift ( comm, dir, disp, source, dest, ierr )
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_cart_shift @ mp_cart_shift" )
#else
  source = 0
  dest = 0
#endif

END SUBROUTINE mp_cart_shift

!******************************************************************************

!..mp_comm_compare
SUBROUTINE mp_comm_compare ( comm1, comm2, RESULT)


    INTEGER, INTENT(IN)                      :: comm1, comm2
    INTEGER, INTENT(OUT)                     :: RESULT

    INTEGER                                  :: ierr, iout

  ierr = 0
  RESULT = 0
#if defined(__parallel)
  CALL mpi_comm_compare ( comm1, comm2, iout, ierr )
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_comm_compare @ mp_comm_compare" )
  IF ( iout == MPI_UNEQUAL ) RESULT = 1
#endif

END SUBROUTINE mp_comm_compare

!******************************************************************************

!..mp_cart_sub
SUBROUTINE mp_cart_sub ( comm, rdim, sub_comm )


    INTEGER, INTENT(IN)                      :: comm
    LOGICAL, DIMENSION(:), INTENT(IN)        :: rdim
    INTEGER, INTENT(OUT)                     :: sub_comm

    INTEGER                                  :: ierr

  ierr = 0
  sub_comm = 0
#if defined(__parallel)
  CALL mpi_cart_sub ( comm, rdim, sub_comm, ierr )
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_cart_sub @ mp_cart_sub" )
#endif

END SUBROUTINE mp_cart_sub

!******************************************************************************

!..mp_comm_free
SUBROUTINE mp_comm_free ( comm )


    INTEGER, INTENT(IN)                      :: comm

    INTEGER                                  :: ierr

  ierr = 0
#if defined(__parallel)
  CALL mpi_comm_free ( comm, ierr )
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_comm_free @ mp_comm_free" )
#endif

END SUBROUTINE mp_comm_free

!******************************************************************************

!..mp_comm_dup
SUBROUTINE mp_comm_dup ( comm1, comm2 )


    INTEGER, INTENT(IN)                      :: comm1
    INTEGER, INTENT(OUT)                     :: comm2

    INTEGER                                  :: ierr

  ierr = 0
#if defined(__parallel)
  CALL mpi_comm_dup ( comm1, comm2, ierr )
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_comm_dup @ mp_comm_dup" )
#else
  comm2=comm1
#endif

END SUBROUTINE mp_comm_dup

!******************************************************************************

!..mp_rank_compare
SUBROUTINE mp_rank_compare ( comm1, comm2, rank )


    INTEGER, INTENT(IN)                      :: comm1, comm2
    INTEGER, DIMENSION(:), INTENT(OUT)       :: rank

    INTEGER                                  :: g1, g2, i, ierr, n, n1, n2
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: rin

  ierr = 0
  rank = 0
#if defined(__parallel)
  CALL mpi_comm_size ( comm1, n1, ierr )
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_comm_size @ mp_rank_compare" )
  CALL mpi_comm_size ( comm2, n2, ierr )
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_comm_size @ mp_rank_compare" )
  n = MAX ( n1, n2 )
  CALL mpi_comm_group ( comm1, g1, ierr )
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_comm_group @ mp_rank_compare" )
  CALL mpi_comm_group ( comm2, g2, ierr )
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_comm_group @ mp_rank_compare" )
  ALLOCATE ( rin ( 0 : n - 1 ), STAT = ierr )
  IF ( ierr /= 0 ) CALL mp_stop( 0, "allocate @ mp_rank_compare" )
  DO i = 0, n-1
    rin ( i ) = i
  END DO
  CALL mpi_group_translate_ranks ( g1, n, rin, g2, rank, ierr )
  IF ( ierr /= 0 ) CALL mp_stop( ierr, &
      "mpi_group_translate_rank @ mp_rank_compare" )
  DEALLOCATE ( rin, STAT = ierr )
  IF ( ierr /= 0 ) CALL mp_stop( 0, "deallocate @ mp_rank_compare" )
#endif

END SUBROUTINE mp_rank_compare

!******************************************************************************

!..mp_dims_create
SUBROUTINE mp_dims_create ( nodes, dims )

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: nodes
  INTEGER, DIMENSION ( : ), INTENT ( INOUT ) :: dims

! Locals
  INTEGER :: ndim, ierr

#if defined(__parallel)
  ndim = SIZE ( dims )
  CALL mpi_dims_create ( nodes, ndim, dims ,ierr)
  IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_dims_create @ mp_dims_create" )
#else
  dims = 1
#endif

END SUBROUTINE mp_dims_create

!******************************************************************************

!******************************************************************************
! this routine will shift around the data in msg
! i.e. msg will be moved from rank to rank+displ_in (in a circular way)
! displ_in will be 1 by default (others not tested)
! the message array needs to be the same size on all processes
!*****************************************************************************
SUBROUTINE mp_shift_im( msg, group, displ_in)

  IMPLICIT NONE
  INTEGER, INTENT ( IN )               :: group
  INTEGER, INTENT(INOUT)               :: msg ( :, : )
  INTEGER, INTENT ( IN ), OPTIONAL     :: displ_in
!
  INTEGER, DIMENSION(:), ALLOCATABLE ::  status
  INTEGER displ,myrank,nprocs,left,right,ierror,msglen,tag,elen
  CHARACTER (LEN=80) :: estring

#if defined(__parallel)
  ALLOCATE(status(MPI_STATUS_SIZE))
  CALL mpi_comm_rank(group,myrank,ierror)
  IF ( ierror /= 0 ) CALL mp_stop ( ierror, "mpi_comm_rank @ mp_shift" )
  CALL mpi_comm_size(group,nprocs,ierror)
  IF ( ierror /= 0 ) CALL mp_stop ( ierror, "mpi_comm_size @ mp_shift" )
  IF (PRESENT(displ_in)) THEN
     displ=displ_in
  ELSE
     displ=1
  ENDIF
  right=MODULO(myrank+displ,nprocs)
  left =MODULO(myrank-displ,nprocs)
  tag=17
  msglen = SIZE(msg)
  t_start = m_cputime ( )
  CALL mpi_sendrecv_replace(msg,msglen,MPI_INTEGER,right,tag,left,tag, &
     group,status(1),ierror)
  t_end = m_cputime ( )
  IF ( ierror /= 0 ) THEN
! this seems not to work (on SUN)
!    CALL mpi_error_string ( ierror, estring, elen )
!    WRITE (*,*) estring(1:elen)
     WRITE (*,*) msglen,status
     CALL mp_stop ( ierror, "mpi_sendrecv_replace @ mp_shift_im" )
  END IF
  CALL add_perf(perf_id=7,count=1,time=t_end-t_start,msg_size=msglen*intlen)
  DEALLOCATE(status)
#endif

END SUBROUTINE mp_shift_im

SUBROUTINE mp_shift_rm( msg, group, displ_in)

  IMPLICIT NONE
  INTEGER, INTENT ( IN )               :: group
  REAL (KIND=dp) , INTENT(INOUT)         :: msg ( :, : )
  INTEGER, INTENT ( IN ), OPTIONAL     :: displ_in
! 
  INTEGER, DIMENSION(:), ALLOCATABLE ::  status
  INTEGER displ,myrank,nprocs,left,right,ierror,msglen,tag,elen
  CHARACTER (LEN=80) :: estring

#if defined(__parallel)
  ALLOCATE(status(MPI_STATUS_SIZE))
  CALL mpi_comm_rank(group,myrank,ierror)
  IF ( ierror /= 0 ) CALL mp_stop ( ierror, "mpi_comm_rank @ mp_shift" )
  CALL mpi_comm_size(group,nprocs,ierror)
  IF ( ierror /= 0 ) CALL mp_stop ( ierror, "mpi_comm_size @ mp_shift" )
  IF (PRESENT(displ_in)) THEN
     displ=displ_in 
  ELSE
     displ=1
  ENDIF
  right=MODULO(myrank+displ,nprocs) 
  left =MODULO(myrank-displ,nprocs) 
  tag=18
  msglen = SIZE(msg)
  t_start = m_cputime ( )
  CALL mpi_sendrecv_replace(msg,msglen,MPI_DOUBLE_PRECISION,right,tag,left,&
                            tag,group,status(1),ierror)
  t_end = m_cputime ( )
  IF ( ierror /= 0 ) THEN
! this seems not to work (on SUN)
!    CALL mpi_error_string ( ierror, estring, elen )
!    WRITE (*,*) estring(1:elen)
     CALL mp_stop ( ierror, "mpi_sendrecv_replace @ mp_shift_rm" )
  END IF
  CALL add_perf(perf_id=7,count=1,time=t_end-t_start,msg_size=msglen*reallen)
  DEALLOCATE(status)
#endif

END SUBROUTINE mp_shift_rm

SUBROUTINE mp_shift_r1( msg, group, displ_in)

  IMPLICIT NONE
  INTEGER, INTENT ( IN )               :: group
  REAL (KIND=dp) , INTENT(INOUT)         :: msg ( : )
  INTEGER, INTENT ( IN ), OPTIONAL     :: displ_in
!
  INTEGER, DIMENSION(:), ALLOCATABLE ::  status
  INTEGER displ,myrank,nprocs,left,right,ierror,msglen,tag,elen
  CHARACTER (LEN=80) :: estring

#if defined(__parallel)
  ALLOCATE(status(MPI_STATUS_SIZE))
  CALL mpi_comm_rank(group,myrank,ierror)
  IF ( ierror /= 0 ) CALL mp_stop ( ierror, "mpi_comm_rank @ mp_shift" )
  CALL mpi_comm_size(group,nprocs,ierror)
  IF ( ierror /= 0 ) CALL mp_stop ( ierror, "mpi_comm_size @ mp_shift" )
  IF (PRESENT(displ_in)) THEN
     displ=displ_in
  ELSE
     displ=1
  ENDIF
  right=MODULO(myrank+displ,nprocs)
  left =MODULO(myrank-displ,nprocs)
  tag=19
  msglen = SIZE(msg)
  t_start = m_cputime ( )
  CALL mpi_sendrecv_replace(msg,msglen,MPI_DOUBLE_PRECISION,right,tag,left,&
                            tag,group,status(1),ierror)
  t_end = m_cputime ( )
  IF ( ierror /= 0 ) THEN
! this seems not to work (on SUN)
!    CALL mpi_error_string ( ierror, estring, elen )
!    WRITE (*,*) estring(1:elen)
     CALL mp_stop ( ierror, "mpi_sendrecv_replace @ mp_shift_r1" )
  END IF
  CALL add_perf(perf_id=7,count=1,time=t_end-t_start,msg_size=msglen*reallen)
  DEALLOCATE(status)
#endif

END SUBROUTINE mp_shift_r1

!******************************************************************************

!..mp_cart_rank
SUBROUTINE mp_cart_rank ( group, pos, rank )

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: group
  INTEGER, DIMENSION ( : ), INTENT ( IN ) :: pos
  INTEGER, INTENT ( OUT ) :: rank

! Locals
  INTEGER :: ierr

#if defined(__parallel)
  CALL mpi_cart_rank ( group, pos, rank, ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_cart_rank @ mp_cart_rank" )
#else
  rank = 1
#endif

END SUBROUTINE mp_cart_rank
!******************************************************************************
!..mp_alltoall

SUBROUTINE mp_alltoall_i11v ( sb, scount, sdispl, rb, rcount, rdispl, group )


    INTEGER, DIMENSION(:), INTENT(IN)        :: sb, scount, sdispl
    INTEGER, DIMENSION(:), INTENT(OUT)       :: rb
    INTEGER, DIMENSION(:), INTENT(IN)        :: rcount, rdispl
    INTEGER, INTENT(IN)                      :: group

    INTEGER                                  :: ierr, msglen

  ierr = 0
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_alltoallv ( sb, scount, sdispl, MPI_INTEGER, &
                       rb, rcount, rdispl, MPI_INTEGER, group, ierr )
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_alltoallv @ mp_alltoall_i11v" )
  t_end = m_cputime ( )
  msglen = SUM ( scount ) + SUM ( rcount )
  CALL add_perf(perf_id=6,count=1,time=t_end-t_start,msg_size=msglen*intlen)
#else
  rb(rdispl(1)+1:rdispl(1)+rcount(1))=sb(sdispl(1)+1:sdispl(1)+scount(1))
#endif

END SUBROUTINE mp_alltoall_i11v

SUBROUTINE mp_alltoall_r11v ( sb, scount, sdispl, rb, rcount, rdispl, group )


    REAL(KIND=dp), DIMENSION(:), INTENT(IN)      :: sb
    INTEGER, DIMENSION(:), INTENT(IN)        :: scount, sdispl
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT)     :: rb
    INTEGER, DIMENSION(:), INTENT(IN)        :: rcount, rdispl
    INTEGER, INTENT(IN)                      :: group

    INTEGER                                  :: ierr, msglen

  ierr = 0
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_alltoallv ( sb, scount, sdispl, MPI_DOUBLE_PRECISION, &
                       rb, rcount, rdispl, MPI_DOUBLE_PRECISION, group, ierr )
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_alltoallv @ mp_alltoall_r11v" )
  msglen = SUM ( scount ) + SUM ( rcount )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=6,count=1,time=t_end-t_start,msg_size=msglen*reallen)
#else
  rb(rdispl(1)+1:rdispl(1)+rcount(1))=sb(sdispl(1)+1:sdispl(1)+scount(1))
#endif

END SUBROUTINE mp_alltoall_r11v


!..mp_alltoall
SUBROUTINE mp_alltoall_c22v ( sb, scount, sdispl, rb, rcount, rdispl, group )


    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: sb
    INTEGER, DIMENSION(:), INTENT(IN)        :: scount, sdispl
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: rb
    INTEGER, DIMENSION(:), INTENT(IN)        :: rcount, rdispl
    INTEGER, INTENT(IN)                      :: group

    INTEGER                                  :: ierr, msglen

  ierr = 0
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_alltoallv ( sb, scount, sdispl, MPI_DOUBLE_COMPLEX, &
                       rb, rcount, rdispl, MPI_DOUBLE_COMPLEX, group, ierr )
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_alltoallv @ mp_alltoall_c22v" )
  msglen = SUM ( scount ) + SUM ( rcount )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=6,count=1,time=t_end-t_start,msg_size=msglen*2*reallen)
#else
  rb=sb
#endif

END SUBROUTINE mp_alltoall_c22v

SUBROUTINE mp_alltoall_i ( sb, rb, count, group )


    INTEGER, DIMENSION(:), INTENT(IN)        :: sb
    INTEGER, DIMENSION(:), INTENT(OUT)       :: rb
    INTEGER, INTENT(IN)                      :: count, group

    INTEGER                                  :: ierr, msglen, np

  ierr = 0
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_alltoall ( sb, count, MPI_INTEGER, &
                      rb, count, MPI_INTEGER, group, ierr )
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_alltoall @ mp_alltoall_i" )
  CALL mpi_comm_size ( group, np, ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_size @ mp_alltoall_i" )
  msglen = 2 * count * np
  t_end = m_cputime ( )
  CALL add_perf(perf_id=6,count=1,time=t_end-t_start,msg_size=msglen*intlen)
#else
  rb=sb
#endif

END SUBROUTINE mp_alltoall_i

!******************************************************************************

!..mp_bcast
SUBROUTINE mp_bcast_i1(msg,source,gid)
    INTEGER                                  :: msg, source, gid

    INTEGER                                  :: ierr, msglen

  ierr = 0
  msglen = 1
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_bcast(msg,msglen,MPI_INTEGER,source,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ mp_bcast_i1" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=2,count=1,time=t_end-t_start,msg_size=msglen*intlen)
#endif
END SUBROUTINE mp_bcast_i1
SUBROUTINE mp_bcast_iv(msg,source,gid)
    INTEGER                                  :: msg( : ), source, gid

    INTEGER                                  :: ierr, msglen

  ierr = 0
  msglen = SIZE(msg)
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_bcast(msg,msglen,MPI_INTEGER,source,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ mp_bcast_iv" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=2,count=1,time=t_end-t_start,msg_size=msglen*intlen)
#endif
END SUBROUTINE mp_bcast_iv
SUBROUTINE mp_bcast_im(msg,source,gid)
    INTEGER                                  :: msg( :, : ), source, gid

    INTEGER                                  :: ierr, msglen

  ierr = 0
  msglen = SIZE(msg)
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_bcast(msg,msglen,MPI_INTEGER,source,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ mp_bcast_im" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=2,count=1,time=t_end-t_start,msg_size=msglen*intlen)
#endif
END SUBROUTINE mp_bcast_im
SUBROUTINE mp_bcast_r1(msg,source,gid)
    REAL(KIND=dp)                                :: msg
    INTEGER                                  :: source, gid

    INTEGER                                  :: ierr, msglen

  ierr = 0
  msglen = 1
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_bcast(msg,msglen,MPI_DOUBLE_PRECISION,source,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ mp_bcast_r1" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=2,count=1,time=t_end-t_start,msg_size=msglen*reallen)
#endif
END SUBROUTINE mp_bcast_r1
SUBROUTINE mp_bcast_rv(msg,source,gid)
    REAL(KIND=dp)                                :: msg( : )
    INTEGER                                  :: source, gid

    INTEGER                                  :: ierr, msglen

  ierr = 0
  msglen = SIZE(msg)
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_bcast(msg,msglen,MPI_DOUBLE_PRECISION,source,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ mp_bcast_rv" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=2,count=1,time=t_end-t_start,msg_size=msglen*reallen)
#endif
END SUBROUTINE mp_bcast_rv
SUBROUTINE mp_bcast_rm(msg,source,gid)
    REAL(KIND=dp)                                :: msg( :, : )
    INTEGER                                  :: source, gid

    INTEGER                                  :: ierr, msglen

  ierr = 0
  msglen = SIZE(msg)
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_bcast(msg,msglen,MPI_DOUBLE_PRECISION,source,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ mp_bcast_rm" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=2,count=1,time=t_end-t_start,msg_size=msglen*reallen)
#endif
END SUBROUTINE mp_bcast_rm
SUBROUTINE mp_bcast_rm3(msg,source,gid)
    REAL(KIND=dp)                                :: msg( :, :, : )
    INTEGER                                  :: source, gid

    INTEGER                                  :: ierr, msglen

  ierr = 0
  msglen = SIZE(msg)
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_bcast(msg,msglen,MPI_DOUBLE_PRECISION,source,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ mp_bcast_rm" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=2,count=1,time=t_end-t_start,msg_size=msglen*reallen)
#endif
END SUBROUTINE mp_bcast_rm3
SUBROUTINE mp_bcast_c1(msg,source,gid)
    COMPLEX(KIND=dp)                             :: msg
    INTEGER                                  :: source, gid

    INTEGER                                  :: ierr, msglen

  ierr = 0
  msglen = 1
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_bcast(msg,msglen,MPI_DOUBLE_COMPLEX,source,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ mp_bcast_c1" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=2,count=1,time=t_end-t_start,msg_size=msglen*2*reallen)
#endif
END SUBROUTINE mp_bcast_c1
SUBROUTINE mp_bcast_cv(msg,source,gid)
    COMPLEX(KIND=dp)                             :: msg( : )
    INTEGER                                  :: source, gid

    INTEGER                                  :: ierr, msglen

  ierr = 0
  msglen = SIZE(msg)
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_bcast(msg,msglen,MPI_DOUBLE_COMPLEX,source,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ mp_bcast_cv" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=2,count=1,time=t_end-t_start,msg_size=msglen*2*reallen)
#endif
END SUBROUTINE mp_bcast_cv
SUBROUTINE mp_bcast_cm(msg,source,gid)
    COMPLEX(KIND=dp)                             :: msg( :, : )
    INTEGER                                  :: source, gid

    INTEGER                                  :: ierr, msglen

  ierr = 0
  msglen = SIZE(msg)
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_bcast(msg,msglen,MPI_DOUBLE_COMPLEX,source,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ mp_bcast_cm")
  t_end = m_cputime ( )
  CALL add_perf(perf_id=2,count=1,time=t_end-t_start,msg_size=msglen*2*reallen)
#endif
END SUBROUTINE mp_bcast_cm
SUBROUTINE mp_bcast_l(msg,source,gid)
    LOGICAL                                  :: msg
    INTEGER                                  :: source, gid

    INTEGER                                  :: ierr, msglen

  ierr = 0
  msglen = 1
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_bcast(msg,msglen,MPI_LOGICAL,source,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ mp_bcast_l" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=2,count=1,time=t_end-t_start,msg_size=msglen*loglen)
#endif
END SUBROUTINE mp_bcast_l
SUBROUTINE mp_bcast_z(msg,source,gid)
    CHARACTER(LEN=*)                         :: msg
    INTEGER                                  :: source, gid

    INTEGER                                  :: i, ierr, msglen
    INTEGER, ALLOCATABLE                     :: imsg( : )

  ierr = 0
  msglen = LEN(msg)
#if defined(__parallel)
  t_start = m_cputime ( )
! this is a workaround to avoid problems on the T3E
! at the moment we have a data alignment error when trying to
! broadcats characters on the T3E (not always!)
! JH 19/3/99 on galileo
! CALL mpi_bcast(msg,msglen,MPI_CHARACTER,source,gid,ierr)
  ALLOCATE (imsg(1:msglen))
  DO i = 1, msglen
     imsg(i) = ICHAR(msg(i:i))
  END DO
  CALL mpi_bcast(imsg,msglen,MPI_INTEGER,source,gid,ierr)
  DO i = 1, msglen
     msg(i:i) = CHAR(imsg(i))
  END DO
  DEALLOCATE (imsg)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ mp_bcast_z" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=2,count=1,time=t_end-t_start,msg_size=msglen*charlen)
#endif
END SUBROUTINE mp_bcast_z

!******************************************************************************

!..mp_sum
SUBROUTINE mp_sum_i1(msg,gid)
    INTEGER, INTENT(INOUT)                   :: msg
    INTEGER, INTENT(IN)                      :: gid

    INTEGER                                  :: ierr, msglen, res

  msglen = 1
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_allreduce(msg,res,msglen,MPI_INTEGER,MPI_SUM,gid,ierr)
  msg = res
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_sum_i1" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*intlen)
#endif
END SUBROUTINE mp_sum_i1
SUBROUTINE mp_sum_iv(msg,gid)
  IMPLICIT NONE
  INTEGER, INTENT ( INOUT ) :: msg ( : )
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
  INTEGER, ALLOCATABLE :: res ( : )
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = SIZE(msg)
  ALLOCATE (res(1:msglen),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop( 0, "allocate @ mp_sum_iv" )
  CALL mpi_allreduce(msg,res,msglen,MPI_INTEGER,MPI_SUM,gid,ierr)
  msg = res
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_sum_iv" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*intlen)
#endif
END SUBROUTINE mp_sum_iv
SUBROUTINE mp_sum_im(msg,gid)
  IMPLICIT NONE
  INTEGER, INTENT ( INOUT ) :: msg ( :, : )
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, m1, m2, ierr
  INTEGER, ALLOCATABLE :: res ( :, : )
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = SIZE(msg)
  m1 = SIZE(msg,1)
  m2 = SIZE(msg,2)
  ALLOCATE (res(m1,m2),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop( 0, "allocate @ mp_sum_im" )
  CALL mpi_allreduce(msg,res,msglen,MPI_INTEGER,MPI_SUM,gid,ierr)
  msg = res
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_sum_im" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*intlen)
#endif
END SUBROUTINE mp_sum_im
SUBROUTINE mp_sum_im3(msg,gid)
  IMPLICIT NONE
  INTEGER, INTENT ( INOUT ) :: msg ( :, :, : )
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, m1, m2, m3, ierr
  INTEGER, ALLOCATABLE :: res ( :, :, : )
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = SIZE(msg)
  m1 = SIZE(msg,1)
  m2 = SIZE(msg,2)
  m3 = SIZE(msg,3)
  ALLOCATE (res(m1,m2,m3),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop( 0, "allocate @ mp_sum_im3" )
  CALL mpi_allreduce(msg,res,msglen,MPI_INTEGER,MPI_SUM,gid,ierr)
  msg = res
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_sum_im3" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*intlen)
#endif
END SUBROUTINE mp_sum_im3
SUBROUTINE mp_sum_r1(msg,gid)
    REAL(KIND=dp), INTENT(INOUT)                 :: msg
    INTEGER, INTENT(IN)                      :: gid

    INTEGER                                  :: ierr, msglen
    REAL(KIND=dp)                                :: res

  msglen = 1
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_allreduce(msg,res,msglen,MPI_DOUBLE_PRECISION,MPI_SUM,gid, &
       ierr)
  msg = res
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_sum_r1" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*reallen)
#endif
END SUBROUTINE mp_sum_r1
SUBROUTINE mp_sum_rv(msg,gid)
  IMPLICIT NONE
  REAL (KIND=dp), INTENT ( INOUT ) :: msg ( : )
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
  REAL (KIND=dp), ALLOCATABLE :: res ( : )
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = SIZE(msg)
  ALLOCATE (res(1:msglen),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop( 0, "allocate @ mp_sum_rv" )
  CALL mpi_allreduce(msg,res,msglen,MPI_DOUBLE_PRECISION,MPI_SUM,gid, &
       ierr)
  msg = res
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_sum_rv" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*reallen)
#endif
END SUBROUTINE mp_sum_rv
SUBROUTINE mp_sum_root_rv(msg,root,gid)
  IMPLICIT NONE
  REAL (KIND=dp), INTENT ( INOUT ) :: msg ( : )
  INTEGER, INTENT ( IN ) :: root,gid
  INTEGER :: msglen, m1, ierr, taskid
  REAL (KIND=dp), ALLOCATABLE :: res ( : )
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_comm_rank ( gid, taskid, ierr )
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_comm_rank @ mp_sum_root_rm" )
  msglen = SIZE(msg)
  m1 = SIZE(msg,1)
  ALLOCATE (res(m1),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop( 0, "allocate @ mp_sum_root_rv" )
  CALL mpi_reduce(msg,res,msglen,MPI_DOUBLE_PRECISION,MPI_SUM,&
       root,gid,ierr)
  IF ( taskid == root ) THEN
    msg = res
  END IF
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_reduce @ mp_sum_root_rv" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*reallen)
#endif
END SUBROUTINE mp_sum_root_rv
SUBROUTINE mp_sum_rm(msg,gid)
  IMPLICIT NONE
  REAL (KIND=dp), INTENT ( INOUT ) :: msg ( :, : )
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, m1, m2, ierr
  REAL (KIND=dp), ALLOCATABLE :: res ( :, : )
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = SIZE(msg)
  m1 = SIZE(msg,1)
  m2 = SIZE(msg,2)
  ALLOCATE (res(m1,m2),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop( 0, "allocate @ mp_sum_rm" )
  CALL mpi_allreduce(msg,res,msglen,MPI_DOUBLE_PRECISION,MPI_SUM,gid, &
       ierr)
  msg = res
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_sum_rm" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*reallen)
#endif
END SUBROUTINE mp_sum_rm
SUBROUTINE mp_sum_root_rm(msg,root,gid)
  IMPLICIT NONE
  REAL (KIND=dp), INTENT ( INOUT ) :: msg ( :, : )
  INTEGER, INTENT ( IN ) :: root,gid
  INTEGER :: msglen, m1, m2, ierr, taskid
  REAL (KIND=dp), ALLOCATABLE :: res ( :, : )
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_comm_rank ( gid, taskid, ierr )
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_comm_rank @ mp_sum_root_rm" )
  msglen = SIZE(msg)
  m1 = SIZE(msg,1)
  m2 = SIZE(msg,2)
  ALLOCATE (res(m1,m2),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop( 0, "allocate @ mp_sum_root_rm" )
  CALL mpi_reduce(msg,res,msglen,MPI_DOUBLE_PRECISION,MPI_SUM,&
       root,gid,ierr)
  IF ( taskid == root ) THEN
    msg = res
  END IF
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_reduce @ mp_sum_root_rm" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*reallen)
#endif
END SUBROUTINE mp_sum_root_rm
SUBROUTINE mp_sum_rm3(msg,gid)
  IMPLICIT NONE
  REAL (KIND=dp), INTENT ( INOUT ) :: msg ( :, :, : )
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, m1, m2, m3, ierr
  REAL (KIND=dp), ALLOCATABLE :: res ( :, :, : )
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = SIZE(msg)
  m1 = SIZE(msg,1)
  m2 = SIZE(msg,2)
  m3 = SIZE(msg,3)
  ALLOCATE (res(m1,m2,m3),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop( 0, "allocate @ mp_sum_rm3" )
  CALL mpi_allreduce(msg,res,msglen,MPI_DOUBLE_PRECISION,MPI_SUM,gid, &
       ierr)
  msg = res
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_sum_rm3" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*reallen)
#endif
END SUBROUTINE mp_sum_rm3
SUBROUTINE mp_sum_root_rm3(msg,root,gid)
  IMPLICIT NONE
  REAL (KIND=dp), INTENT ( INOUT ) :: msg ( :, :, : )
  INTEGER, INTENT ( IN ) :: root,gid
  INTEGER :: msglen, m1, m2, m3, ierr, taskid
  REAL (KIND=dp), ALLOCATABLE :: res ( :, :, : )
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_comm_rank ( gid, taskid, ierr )
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_comm_rank @ mp_sum_root_rm3" )
  msglen = SIZE(msg)
  m1 = SIZE(msg,1)
  m2 = SIZE(msg,2)
  m3 = SIZE(msg,3)
  ALLOCATE (res(m1,m2,m3),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop( 0, "allocate @ mp_sum_root_rm3" )
  CALL mpi_reduce(msg,res,msglen,MPI_DOUBLE_PRECISION,MPI_SUM,&
       root,gid,ierr)
  IF ( taskid == root ) THEN
    msg = res
  END IF
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_reduce @ mp_sum_root_rm3" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*reallen)
#endif
END SUBROUTINE mp_sum_root_rm3
SUBROUTINE mp_sum_c1(msg,gid)
    COMPLEX(KIND=dp), INTENT(INOUT)              :: msg
    INTEGER, INTENT(IN)                      :: gid

    COMPLEX(KIND=dp)                             :: res
    INTEGER                                  :: ierr, msglen

  msglen = 2
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_allreduce(msg,res,msglen,MPI_DOUBLE_PRECISION,MPI_SUM,gid, &
       ierr)
  msg = res
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_sum_c1" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*2*reallen)
#endif
END SUBROUTINE mp_sum_c1
SUBROUTINE mp_sum_cv(msg,gid)
  IMPLICIT NONE
  COMPLEX (KIND=dp), INTENT ( INOUT ) :: msg ( : )
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
  COMPLEX (KIND=dp), ALLOCATABLE :: res ( : )
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = SIZE(msg)
  ALLOCATE (res(1:msglen),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop( 0, "allocate @ mp_sum_cv" )
  CALL mpi_allreduce(msg,res,2*msglen,MPI_DOUBLE_PRECISION,MPI_SUM,gid, &
       ierr)
  msg = res
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_sum_cv" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*2*reallen)
#endif
END SUBROUTINE mp_sum_cv
SUBROUTINE mp_sum_cm(msg,gid)
  IMPLICIT NONE
  COMPLEX (KIND=dp), INTENT ( INOUT ) :: msg ( :, : )
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, m1, m2, ierr
  COMPLEX (KIND=dp), ALLOCATABLE :: res ( :, : )
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = 2*SIZE(msg)
  m1 = SIZE(msg,1)
  m2 = SIZE(msg,2)
  ALLOCATE (res(m1,m2),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop( 0, "allocate @ mp_sum_cm" )
  CALL mpi_allreduce(msg,res,msglen,MPI_DOUBLE_PRECISION,MPI_SUM,gid, &
       ierr)
  msg = res
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_sum_cm" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*2*reallen)
#endif
END SUBROUTINE mp_sum_cm
SUBROUTINE mp_sum_cm3(msg,gid)
  IMPLICIT NONE
  COMPLEX (KIND=dp), INTENT ( INOUT ) :: msg ( :, :, : )
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, m1, m2, m3, ierr
  COMPLEX (KIND=dp), ALLOCATABLE :: res ( :, :, : )
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = 2*SIZE(msg)
  m1 = SIZE(msg,1)
  m2 = SIZE(msg,2)
  m3 = SIZE(msg,3)
  ALLOCATE (res(m1,m2,m3),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop( 0, "allocate @ mp_sum_cm3" )
  CALL mpi_allreduce(msg,res,msglen,MPI_DOUBLE_PRECISION,MPI_SUM,gid, &
       ierr)
  msg = res
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_sum_cm3" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*2*reallen)
#endif
END SUBROUTINE mp_sum_cm3

!******************************************************************************

SUBROUTINE mp_max_i(msg,gid)
  IMPLICIT NONE
  INTEGER, INTENT ( INOUT ) :: msg
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
  INTEGER :: res
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = 1
  CALL mpi_allreduce(msg,res,msglen,MPI_INTEGER,MPI_MAX,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_max_i" )
  msg = res
  t_end = m_cputime ( )
  CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*intlen)
#endif
END SUBROUTINE mp_max_i
SUBROUTINE mp_max_iv(msg,gid)
  IMPLICIT NONE
  INTEGER, INTENT ( INOUT ) :: msg ( : )
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
  INTEGER, ALLOCATABLE :: res ( : )
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = SIZE(msg)
  ALLOCATE (res(1:msglen),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop( 0, "allocate @ mp_max_iv" )
  CALL mpi_allreduce(msg,res,msglen,MPI_INTEGER,MPI_MAX,gid,ierr)
  msg = res
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_max_iv" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*intlen)
#endif
END SUBROUTINE mp_max_iv
SUBROUTINE mp_max_r(msg,gid)
  IMPLICIT NONE
  REAL (KIND=dp), INTENT ( INOUT ) :: msg
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
  REAL (KIND=dp) :: res
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = 1
  CALL mpi_allreduce(msg,res,msglen,MPI_DOUBLE_PRECISION,MPI_MAX,gid, &
       ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_max_r" )
  msg = res
  t_end = m_cputime ( )
  CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*reallen)
#endif
END SUBROUTINE mp_max_r

SUBROUTINE mp_maxloc_i(msg,gid)
  IMPLICIT NONE
  INTEGER, INTENT ( INOUT ) :: msg ( : )
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
  INTEGER, ALLOCATABLE :: res ( : )
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = SIZE(msg)
  ALLOCATE (res(1:msglen),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop( 0, "allocate @ mp_maxloc_i" )
  CALL mpi_allreduce(msg,res,1,MPI_2INTEGER,MPI_MAXLOC,gid, &
       ierr)
  msg = res
  DEALLOCATE (res,STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_maxloc_i" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*intlen)
#endif
END SUBROUTINE mp_maxloc_i

SUBROUTINE mp_minloc_rv(msg,gid)
  IMPLICIT NONE
  REAL (KIND=dp), INTENT ( INOUT ) :: msg ( : )
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
  REAL (KIND=dp), ALLOCATABLE :: res ( : )
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = SIZE(msg)
  ALLOCATE (res(1:msglen),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop( 0, "allocate @ mp_min_rv" )
  CALL mpi_allreduce(msg,res,1,MPI_2DOUBLE_PRECISION,MPI_MINLOC,gid, &
       ierr)
  msg = res
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_min_rv" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*reallen)
#endif
END SUBROUTINE mp_minloc_rv

SUBROUTINE mp_max_rv(msg,gid)
  IMPLICIT NONE
  REAL (KIND=dp), INTENT ( INOUT ) :: msg ( : )
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
  REAL (KIND=dp), ALLOCATABLE :: res ( : )
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = SIZE(msg)
  ALLOCATE (res(1:msglen),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop( 0, "allocate @ mp_max_rv" )
  CALL mpi_allreduce(msg,res,msglen,MPI_DOUBLE_PRECISION,MPI_MAX,gid, &
       ierr)
  msg = res
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_max_rv" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*reallen)
#endif
END SUBROUTINE mp_max_rv

!******************************************************************************

SUBROUTINE mp_min_i(msg,gid)
  IMPLICIT NONE
  INTEGER, INTENT ( INOUT ) :: msg
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
  INTEGER :: res
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = 1
  CALL mpi_allreduce(msg,res,msglen,MPI_INTEGER,MPI_MIN,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_min_i" )
  msg = res
  t_end = m_cputime ( )
  CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*intlen)
#endif
END SUBROUTINE mp_min_i
SUBROUTINE mp_min_iv(msg,gid)
  IMPLICIT NONE
  INTEGER, INTENT ( INOUT ) :: msg ( : )
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
  INTEGER, ALLOCATABLE :: res ( : )
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = SIZE(msg)
  ALLOCATE (res(1:msglen),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop( 0, "allocate @ mp_min_iv" )
  CALL mpi_allreduce(msg,res,msglen,MPI_INTEGER,MPI_MIN,gid,ierr)
  msg = res
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_min_iv" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*intlen)
#endif
END SUBROUTINE mp_min_iv
SUBROUTINE mp_min_r(msg,gid)
  IMPLICIT NONE
  REAL (KIND=dp), INTENT ( INOUT ) :: msg
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
  REAL (KIND=dp) :: res
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = 1
  CALL mpi_allreduce(msg,res,msglen,MPI_DOUBLE_PRECISION,MPI_MIN,gid, &
       ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_min_r" )
  msg = res
  t_end = m_cputime ( )
  CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*reallen)
#endif
END SUBROUTINE mp_min_r
SUBROUTINE mp_min_rv(msg,gid)
  IMPLICIT NONE
  REAL (KIND=dp), INTENT ( INOUT ) :: msg ( : )
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
  REAL (KIND=dp), ALLOCATABLE :: res ( : )
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = SIZE(msg)
  ALLOCATE (res(1:msglen),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop( 0, "allocate @ mp_min_rv" )
  CALL mpi_allreduce(msg,res,msglen,MPI_DOUBLE_PRECISION,MPI_MIN,gid, &
       ierr)
  msg = res
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_min_rv" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*reallen)
#endif
END SUBROUTINE mp_min_rv

!******************************************************************************

SUBROUTINE mp_gather_i(msg,msg_gather,root,gid)
  IMPLICIT NONE
  INTEGER, INTENT ( IN ) :: msg
  INTEGER, INTENT ( OUT ) :: msg_gather ( : )
  INTEGER, INTENT ( IN ) :: root
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = 1
  CALL mpi_gather(msg,msglen,MPI_INTEGER,msg_gather,&
                  msglen,MPI_INTEGER,root,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_gather @ mp_gather_i" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=4,count=1,time=t_end-t_start,msg_size=msglen*intlen)
#else
  msg_gather = msg
#endif
END SUBROUTINE mp_gather_i

SUBROUTINE mp_gather_iv(msg,msg_gather,root,gid)
  IMPLICIT NONE
  INTEGER, INTENT ( IN ) :: msg ( : )
  INTEGER, INTENT ( OUT ) :: msg_gather ( : )
  INTEGER, INTENT ( IN ) :: root
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = SIZE(msg)
  CALL mpi_gather(msg,msglen,MPI_INTEGER,msg_gather,&
                  msglen,MPI_INTEGER,root,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_gather @ mp_gather_iv" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=4,count=1,time=t_end-t_start,msg_size=msglen*intlen)
#else
  msg_gather = msg
#endif
END SUBROUTINE mp_gather_iv

SUBROUTINE mp_gather_r(msg,msg_gather,root,gid)
  IMPLICIT NONE
  REAL (KIND=dp), INTENT ( IN ) :: msg
  REAL (KIND=dp), INTENT ( OUT ) :: msg_gather ( : )
  INTEGER, INTENT ( IN ) :: root
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = 1
  CALL mpi_gather(msg,msglen,MPI_DOUBLE_PRECISION,msg_gather,&
                  msglen,MPI_DOUBLE_PRECISION,root,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_gather @ mp_gather_r" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=4,count=1,time=t_end-t_start,msg_size=msglen*reallen)
#else
  msg_gather = msg
#endif
END SUBROUTINE mp_gather_r

SUBROUTINE mp_gather_rv(msg,msg_gather,root,gid)
  IMPLICIT NONE
  REAL (KIND=dp), INTENT ( IN ) :: msg ( : )
  REAL (KIND=dp), INTENT ( OUT ) :: msg_gather ( : )
  INTEGER, INTENT ( IN ) :: root
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = SIZE(msg)
  CALL mpi_gather(msg,msglen,MPI_DOUBLE_PRECISION,msg_gather,&
                  msglen,MPI_DOUBLE_PRECISION,root,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_gather @ mp_gather_rv" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=4,count=1,time=t_end-t_start,msg_size=msglen*reallen)
#else
  msg_gather = msg
#endif
END SUBROUTINE mp_gather_rv

! wrapper for MPI_ALLGATHERV subroutine
! one dimensional real arrays

SUBROUTINE mp_allgather_i1(msgout,msgin,gid)
  IMPLICIT NONE
  INTEGER, INTENT ( IN ) :: msgout
  INTEGER, INTENT ( OUT ) :: msgin ( : )
  INTEGER, INTENT ( IN ) :: gid

  INTEGER :: ierr, scount, rcount

#if defined(__parallel)
  scount = 1
  rcount = 1
  CALL MPI_ALLGATHER(msgout, scount, MPI_INTEGER, &
                     msgin , rcount, MPI_INTEGER, &
                     gid, ierr )
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allgather @ mp_allgather_i1" )
#else
  msgin = msgout
  ierr = 0
#endif
END SUBROUTINE mp_allgather_i1

SUBROUTINE mp_allgatherv_rv(msgout,msgin,rcount,rdispl,gid)
  IMPLICIT NONE
  REAL (KIND=dp), INTENT ( IN ) :: msgout ( : )
  REAL (KIND=dp), INTENT ( OUT ) :: msgin ( : )
  INTEGER, INTENT ( IN ) :: rcount ( : ), rdispl ( : )
  INTEGER, INTENT ( IN ) :: gid

  INTEGER :: ierr, scount

#if defined(__parallel)
  scount = SIZE ( msgout )
  CALL MPI_ALLGATHERV(msgout, scount, MPI_DOUBLE_PRECISION, msgin, rcount, &
                      rdispl, MPI_DOUBLE_PRECISION, gid, ierr )
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allgatherv @ mp_allgatherv_rv" )
#else
  msgin = msgout
  ierr = 0
#endif
END SUBROUTINE mp_allgatherv_rv

!******************************************************************************
! wrapper for MPI_REDUCE_SCATTER subroutine
! one dimensional real arrays
SUBROUTINE mp_sum_scatter_rv(msgout,msgin,rcount,gid)
  IMPLICIT NONE
  REAL (KIND=dp), INTENT ( IN ) :: msgout ( : )
  REAL (KIND=dp), INTENT ( OUT ) :: msgin ( : )
  INTEGER, INTENT ( IN ) :: rcount ( : )
  INTEGER, INTENT ( IN ) :: gid

  INTEGER :: ierr

#if defined(__parallel)
  t_start = m_cputime ( )
  CALL MPI_REDUCE_SCATTER(msgout, msgin, rcount, MPI_DOUBLE_PRECISION, MPI_SUM, &
                          gid, ierr )
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_reduce_scatter @ mp_sum_scatter_rv" )

  t_end = m_cputime ( )
  CALL add_perf(perf_id=3,count=1,time=t_end-t_start,&
       msg_size=rcount(1)*2*reallen)
#else
  msgin = msgout
  ierr = 0
#endif
END SUBROUTINE mp_sum_scatter_rv

!******************************************************************************

SUBROUTINE mp_sendrecv_rm3(msgin,dest,msgout,source,comm)
  IMPLICIT NONE
  REAL (KIND=dp), INTENT ( IN ) :: msgin ( :, :, : )
  REAL (KIND=dp), INTENT ( OUT ) :: msgout ( :, :, : )
  INTEGER, INTENT ( IN ) :: dest, source, comm
  INTEGER :: msglen, ierr, send_tag, recv_tag
  INTEGER, DIMENSION(:), ALLOCATABLE ::  status
#if defined(__parallel)
  ALLOCATE(status(MPI_STATUS_SIZE))
  t_start = m_cputime ( )
  msglen = SIZE(msgin)
  send_tag = 0 ! cannot think of something better here, this might be dangerous
  recv_tag = 0 ! cannot think of something better here, this might be dangerous
  CALL mpi_sendrecv(msgin,msglen,MPI_DOUBLE_PRECISION,dest,send_tag,msgout,&
       msglen,MPI_DOUBLE_PRECISION,source,recv_tag,comm,status(1),ierr)
  ! we do not check the status
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_sendrecv @ mp_send_recv_rm3" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=7,count=1,time=t_end-t_start,msg_size=msglen*reallen)
  DEALLOCATE(status)
#else
  msgout = msgin
#endif
END SUBROUTINE mp_sendrecv_rm3

!******************************************************************************

SUBROUTINE mp_sendrecv_rm2(msgin,dest,msgout,source,comm)
  IMPLICIT NONE
  REAL (KIND=dp), INTENT ( IN ) :: msgin ( :, : )
  REAL (KIND=dp), INTENT ( OUT ) :: msgout ( :, : )
  INTEGER, INTENT ( IN ) :: dest, source, comm
  INTEGER :: msglen, ierr, send_tag, recv_tag
  INTEGER, DIMENSION(:), ALLOCATABLE ::  status
#if defined(__parallel)
  ALLOCATE(status(MPI_STATUS_SIZE))
  t_start = m_cputime ( )
  msglen = SIZE(msgin,1)*SIZE(msgin,2)
  send_tag = 0 ! cannot think of something better here, this might be dangerous
  recv_tag = 0 ! cannot think of something better here, this might be dangerous
  CALL mpi_sendrecv(msgin,msglen,MPI_DOUBLE_PRECISION,dest,send_tag,msgout,&
       msglen,MPI_DOUBLE_PRECISION,source,recv_tag,comm,status(1),ierr)
  ! we do not check the status
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_sendrecv @ mp_send_recv_rm2" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=7,count=1,time=t_end-t_start,msg_size=msglen*reallen)
  DEALLOCATE(status)
#else
  msgout = msgin
#endif
END SUBROUTINE mp_sendrecv_rm2
!******************************************************************************

SUBROUTINE mp_sendrecv_rv(msgin,dest,msgout,source,comm)
  IMPLICIT NONE
  REAL (KIND=dp), INTENT ( IN ) :: msgin ( : )
  REAL (KIND=dp), INTENT ( OUT ) :: msgout ( : )
  INTEGER, INTENT ( IN ) :: dest, source, comm
  INTEGER :: msglen_in, msglen_out, ierr, send_tag, recv_tag
  INTEGER, DIMENSION(:), ALLOCATABLE ::  status
#if defined(__parallel)
  ALLOCATE(status(MPI_STATUS_SIZE))
  t_start = m_cputime ( )
  msglen_in = SIZE(msgin)
  msglen_out = SIZE(msgout)
  send_tag = 0 ! cannot think of something better here, this might be dangerous
  recv_tag = 0 ! cannot think of something better here, this might be dangerous
  CALL mpi_sendrecv(msgin,msglen_in,MPI_DOUBLE_PRECISION,dest,send_tag,msgout,&
       msglen_out,MPI_DOUBLE_PRECISION,source,recv_tag,comm,status(1),ierr)
  ! we do not check the status
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_sendrecv @ mp_send_recv_rv" )
  t_end = m_cputime ( )
  CALL add_perf(perf_id=7,count=1,time=t_end-t_start,&
       msg_size=(msglen_in+msglen_out)*reallen/2)
  DEALLOCATE(status)
#else
  msgout = msgin
#endif
END SUBROUTINE mp_sendrecv_rv

!******************************************************************************

!!****f* message_passing/mp_isendrecv_rv [1.0] *
!!
!!   NAME
!!     mp_isendrecv_rv
!!
!!   FUNCTION
!!     non blocking send of real arrays
!!
!!   NOTES
!!     the arguments must be pointers to be sure that we do not get a temporaries.
!!     they must point to contiguous memory
!!
!!   MODIFICATION HISTORY
!!     11.2004 created [Joost VandeVondele]
!!
!!*** **********************************************************************
SUBROUTINE mp_isendrecv_rv(msgin,dest,msgout,source,comm,send_request,&
                           recv_request,tag)
  IMPLICIT NONE
  REAL (KIND=dp), POINTER, DIMENSION (:) :: msgin,msgout
  INTEGER, INTENT ( IN ) :: dest, source, comm
  INTEGER, INTENT(out) :: send_request, recv_request
  INTEGER, OPTIONAL, INTENT(in) :: tag
  INTEGER :: msglen, ierr, my_tag
  INTEGER :: lower1
  REAL(KIND=dp) :: foo
#if defined(__parallel)
  t_start = m_cputime ( )
  my_tag = 0
  IF (PRESENT(tag)) my_tag=tag

  msglen = SIZE(msgout,1)
  IF (msglen>0) THEN
     lower1=LBOUND(msgout,1)
     CALL mpi_irecv(msgout(lower1),msglen,MPI_DOUBLE_PRECISION,source, my_tag,&
          comm,recv_request,ierr)
  ELSE
     CALL mpi_irecv(foo,msglen,MPI_DOUBLE_PRECISION,source, my_tag,&
          comm,recv_request,ierr)
  END IF
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_ircv @ mp_isendrecv_rm2" )

  msglen = SIZE(msgin,1)
  IF (msglen>0) THEN
     lower1=LBOUND(msgin,1)
     CALL mpi_isend(msgin(lower1),msglen,MPI_DOUBLE_PRECISION,dest,my_tag,&
          comm,send_request,ierr)
  ELSE
     CALL mpi_isend(foo,msglen,MPI_DOUBLE_PRECISION,dest,my_tag,&
          comm,send_request,ierr)
  END IF
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_isend @ mp_isendrecv_rm2" )

  msglen = (msglen+SIZE(msgout,1))/2.0_dp
  t_end = m_cputime ( )
  CALL add_perf(perf_id=8,count=1,time=t_end-t_start,msg_size=msglen*reallen)
#else
  send_request=0
  recv_request=0
  msgout = msgin
#endif
END SUBROUTINE mp_isendrecv_rv


!!****f* message_passing/mp_isendrecv_rm2 [1.0] *
!!
!!   NAME
!!     mp_isendrecv_rm2
!!
!!   FUNCTION
!!     non blocking send of real arrays 
!!
!!   NOTES
!!     the arguments must be pointers to be sure that we do not get a temporaries.
!!     they must point to contiguous memory
!!     
!!   AUTHOR
!!     joost & fawzi
!!
!!   MODIFICATION HISTORY
!!     08.2003 created [f&j]
!!
!!*** **********************************************************************
SUBROUTINE mp_isendrecv_rm2(msgin,dest,msgout,source,comm,send_request,&
     recv_request,tag)
  IMPLICIT NONE
  REAL (KIND=dp), POINTER, DIMENSION (:,:) :: msgin,msgout
  INTEGER, INTENT ( IN ) :: dest, source, comm
  INTEGER, INTENT(out) :: send_request, recv_request
  INTEGER, OPTIONAL, INTENT(in) :: tag
  INTEGER :: msglen, ierr, my_tag
  INTEGER :: lower1,lower2
  REAL(KIND=dp) :: foo

#if defined(__parallel)
  t_start = m_cputime ( )
  my_tag = 0
  IF (PRESENT(tag)) my_tag=tag

  msglen = SIZE(msgout,1)*SIZE(msgout,2)
  IF (msglen>0) THEN
     lower1=LBOUND(msgout,1)
     lower2=LBOUND(msgout,2)
     CALL mpi_irecv(msgout(lower1,lower2),msglen,MPI_DOUBLE_PRECISION,source, my_tag,&
          comm,recv_request,ierr)
  ELSE
     CALL mpi_irecv(foo,msglen,MPI_DOUBLE_PRECISION,source, my_tag,&
          comm,recv_request,ierr)
  END IF
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_ircv @ mp_isendrecv_rm2" )

  msglen = SIZE(msgin,1)*SIZE(msgin,2)
  IF (msglen>0) THEN
     lower1=LBOUND(msgin,1)
     lower2=LBOUND(msgin,2)
     CALL mpi_isend(msgin(lower1,lower2),msglen,MPI_DOUBLE_PRECISION,dest,my_tag,&
          comm,send_request,ierr)
  ELSE
     CALL mpi_isend(foo,msglen,MPI_DOUBLE_PRECISION,dest,my_tag,&
          comm,send_request,ierr)
  END IF
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_isend @ mp_isendrecv_rm2" )

  msglen = (msglen+SIZE(msgout,1)*SIZE(msgout,2))/2.0_dp
  t_end = m_cputime ( )
  CALL add_perf(perf_id=8,count=1,time=t_end-t_start,msg_size=msglen*reallen)
#else
  send_request=0
  recv_request=0 
  msgout = msgin
#endif
END SUBROUTINE mp_isendrecv_rm2
!***************************************************************************

!!****f* message_passing/mp_isend_rm2 [1.0] *
!!
!!   NAME
!!     mp_isend_rm2
!!
!!   FUNCTION
!!     non blocking send of real arrays 
!!
!!   NOTES
!!     the arguments must be pointers to be sure that we do not get a temporaries.
!!     
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE mp_isend_rm2(msgin,dest,comm,request,tag)
  IMPLICIT NONE
  REAL (KIND=dp), POINTER, DIMENSION (:,:) :: msgin
  INTEGER, INTENT ( IN ) :: dest, comm
  INTEGER, INTENT(out) :: request
  INTEGER, OPTIONAL, INTENT(in) :: tag

  INTEGER :: msglen, ierr, my_tag
  INTEGER :: lower1,lower2
  REAL(KIND=dp) :: foo(1)

#if defined(__parallel)
  t_start = m_cputime ( )
  my_tag = 0
  IF (PRESENT(tag)) my_tag=tag

  msglen = SIZE(msgin,1)*SIZE(msgin,2)
  IF (msglen>0) THEN
     lower1=LBOUND(msgin,1)
     lower2=LBOUND(msgin,2)
     CALL mpi_isend(msgin(lower1,lower2),msglen,MPI_DOUBLE_PRECISION,dest,my_tag,&
          comm,request,ierr)
  ELSE
     CALL mpi_isend(foo,msglen,MPI_DOUBLE_PRECISION,dest,my_tag,&
          comm,request,ierr)
  END IF
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_isend @ mp_isend_rm2" )

  t_end = m_cputime ( )
  CALL add_perf(perf_id=11,count=1,time=t_end-t_start,msg_size=msglen*reallen)
#else
  CALL mp_stop( ierr, "mp_isend called in non parallel case" )
#endif
END SUBROUTINE mp_isend_rm2
!***************************************************************************

!!****f* message_passing/mp_irecv_rm2 [1.0] *
!!
!!   NAME
!!     mp_irecv_rm2
!!
!!   FUNCTION
!!     non blocking recieve of real arrays 
!!
!!   NOTES
!!     the arguments must be pointers to be sure that we do not get a temporaries.
!!     
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE mp_irecv_rm2(msgout,source,comm,request,tag)
  IMPLICIT NONE
  REAL (KIND=dp), POINTER, DIMENSION (:,:) :: msgout
  INTEGER, INTENT ( IN ) :: source, comm
  INTEGER, INTENT(out) :: request
  INTEGER, OPTIONAL, INTENT(in) :: tag
  INTEGER :: msglen, ierr, my_tag
  INTEGER :: lower1,lower2
  REAL(KIND=dp) :: foo(1)

#if defined(__parallel)
  t_start = m_cputime ( )
  my_tag = 0
  IF (PRESENT(tag)) my_tag=tag

  msglen = SIZE(msgout,1)*SIZE(msgout,2)
  IF (msglen>0) THEN
     lower1=LBOUND(msgout,1)
     lower2=LBOUND(msgout,2)
     CALL mpi_irecv(msgout(lower1,lower2),msglen,MPI_DOUBLE_PRECISION,source, my_tag,&
          comm,request,ierr)
  ELSE
     CALL mpi_irecv(foo,msglen,MPI_DOUBLE_PRECISION,source, my_tag,&
          comm,request,ierr)
  END IF
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_ircv @ mp_irecv_rm2" )

  t_end = m_cputime ( )
  CALL add_perf(perf_id=12,count=1,time=t_end-t_start,msg_size=msglen*reallen)
#else
  CALL mp_stop( ierr, "mp_irecv called in non parallel case" )
#endif
END SUBROUTINE mp_irecv_rm2
!***************************************************************************

!!****f* message_passing/mp_wait [1.0] *
!!
!!   NAME
!!     mp_wait
!!
!!   FUNCTION
!!     waits for completion of the given request
!!
!!   NOTES
!!     see isendrecv
!!
!!   AUTHOR
!!     joost & fawzi
!!
!!   MODIFICATION HISTORY
!!     08.2003 created [f&j]
!!
!!*** **********************************************************************
SUBROUTINE mp_wait(request)
  INTEGER, INTENT(inout) :: request

  INTEGER, DIMENSION(:), ALLOCATABLE :: status
  INTEGER :: ierr

#if defined(__parallel)
  ALLOCATE(status(MPI_STATUS_SIZE))
  t_start = m_cputime ( )

  CALL mpi_wait(request,status,ierr)
  ! we do not check the status
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_wait @ mp_wait" )

  t_end = m_cputime ( )
  CALL add_perf(perf_id=9,count=1,time=t_end-t_start)
  DEALLOCATE(status)
#endif
END SUBROUTINE mp_wait
!***************************************************************************

!!****f* message_passing/mp_waitall [1.0] *
!!
!!   NAME
!!     mp_waitall
!!
!!   FUNCTION
!!     waits for completion of the given requests
!!
!!   NOTES
!!     see isendrecv
!!
!!   AUTHOR
!!     joost & fawzi
!!
!!   MODIFICATION HISTORY
!!     08.2003 created [f&j]
!!
!!*** **********************************************************************
SUBROUTINE mp_waitall(requests)
  INTEGER, INTENT(inout), DIMENSION(:) :: requests

  INTEGER, DIMENSION(:,:), ALLOCATABLE :: status
  INTEGER :: ierr,count

#if defined(__parallel)
  count=SIZE(requests)
  ALLOCATE(status(MPI_STATUS_SIZE,count))
  t_start = m_cputime ( )

  CALL mpi_waitall(count,requests,status,ierr)
  ! we do not check the status
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_waitall @ mp_waitall" )

  t_end = m_cputime ( )
  CALL add_perf(perf_id=9,count=1,time=t_end-t_start)
  DEALLOCATE(status)
#endif
END SUBROUTINE mp_waitall
!***************************************************************************

!!****f* message_passing/mp_comm_split [1.0] *
!!
!!   NAME
!!     mp_comm_split
!!
!!   FUNCTION
!!     splits the given communicator in group in subgroups trying to organize
!!     them in a way that the communication within each subgroup is
!!     efficent (but not necessarily the comunication between subgroups)
!!
!!   NOTES
!!     at least one of subgroup_min_size and n_subgroups is needed,
!!     the other default to the value needed to use most processors.
!!     if less cpus are present than needed for subgroup min size, n_subgroups,
!!     just one comm is created that contains all cpus
!!
!!   INPUTS
!!     - comm: the mpi communicator that you want to split
!!     - sub_comm: the communicator for the subgroup (created, needs to be freed later)
!!     - ngroups: actual number of groups
!!     - group_distribution  input  : allocated with array with the nprocs entries (0 .. nprocs-1)
!!                           output : the number of the group the proc belongs to (0..ngroups-1)
!!     - subgroup_min_size: the minimum size of the subgroup
!!     - n_subgroups: the number of subgroups wanted
!!     - group_partition: n_subgroups sized array containing the number of cpus wanted per group.
!!                        should match the total number of cpus (only used if present and associated) (0..ngroups-1)
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     10.2003 created [fawzi]
!!     02.2004 modified [Joost VandeVondele]
!!
!!*** **********************************************************************
SUBROUTINE mp_comm_split(comm,sub_comm,ngroups, group_distribution, &
                    subgroup_min_size, n_subgroups, group_partition)
    INTEGER, INTENT(in)                      :: comm
    INTEGER, INTENT(out)                     :: sub_comm, ngroups
    INTEGER, DIMENSION(:), POINTER           :: group_distribution
    INTEGER, INTENT(in), OPTIONAL            :: subgroup_min_size, n_subgroups
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: group_partition

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_comm_split', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: color, i, ierr, j, k, mepos, &
                                                my_subgroup_min_size, nnodes

! actual number of groups

    IF (.NOT. PRESENT(subgroup_min_size) .AND. .NOT. PRESENT(n_subgroups)) THEN
       CALL mp_stop(1,routineP//" missing arguments ")
    ENDIF
    IF (PRESENT(subgroup_min_size) .AND. PRESENT(n_subgroups)) THEN
       CALL mp_stop(1,routineP//" too many arguments ")
    ENDIF

    CALL mp_environ(nnodes,mepos,comm)

    IF (.NOT. ASSOCIATED(group_distribution)) THEN
       CALL mp_stop(1,routineP//" group_distribution not associated")
    ENDIF
    IF (LBOUND(group_distribution,1) .NE. 0 .OR. &
        UBOUND(group_distribution,1).NE.nnodes-1) THEN
       CALL mp_stop(1,routineP//" group_distribution wrong bounds")
    ENDIF
  
#if defined(__parallel)
    t_start = m_cputime ( )
    IF (PRESENT(subgroup_min_size)) THEN
       IF (subgroup_min_size<0) THEN
          CALL mp_stop(mpi_err_comm,routineP//" subgroup_min_size too small")
       ENDIF
       ngroups= nnodes / subgroup_min_size
       my_subgroup_min_size = subgroup_min_size
    ELSE  ! n_subgroups
       IF (n_subgroups<=0) THEN
          CALL mp_stop(mpi_err_comm,routineP//" n_subgroups too small")
       ENDIF
       IF (nnodes/n_subgroups > 0) THEN ! we have a least one cpu per group
          ngroups = n_subgroups 
       ELSE ! well, only one group then
          ngroups = 1
       ENDIF
       my_subgroup_min_size = nnodes / ngroups
    ENDIF
    DO i=0,nnodes-1
       group_distribution(i)=i / my_subgroup_min_size 
       ! if part of the rest, join the last group
       IF ( group_distribution(i) >= ngroups ) group_distribution(i)=ngroups-1
    ENDDO
    ! even the user gave a partition, see if we can use it to overwrite this choice
    IF (PRESENT(group_partition)) THEN
       IF (ASSOCIATED(group_partition)) THEN
          IF (ALL(group_partition>0) .AND. SUM(group_partition).EQ.nnodes .AND. ngroups == SIZE(group_partition)) THEN
             k=0
             DO i=0,SIZE(group_partition)-1
                DO j=1,group_partition(i)
                   group_distribution(k)=i
                   k=k+1
                ENDDO
             ENDDO
          ELSE
             ! just ignore silently as we have reasonable defaults. Probably a warning would not be to bad         
          ENDIF 
       ENDIF
    ENDIF
    color=group_distribution(mepos)
    CALL mpi_comm_split(comm,color,0,sub_comm,ierr)
    IF (ierr/=mpi_success) CALL mp_stop(ierr,"in "//routineP//" split")
     
    t_end = m_cputime ( )
    CALL add_perf(perf_id=10,count=1,time=t_end-t_start)
#else
    CALL mp_comm_dup(comm,sub_comm) 
    group_distribution(0)=0
    ngroups=1
#endif

END SUBROUTINE mp_comm_split

!***************************************************************************
! a bunch of currently unused features
!***************************************************************************
!..mp_group
SUBROUTINE mp_group(group_list,group_size,base_group,groupid)


    INTEGER, INTENT(IN)                      :: group_list( : ), group_size, &
                                                base_group
    INTEGER, INTENT(OUT)                     :: groupid

    INTEGER                                  :: base, ierr, newgroup

  ierr = 0
  groupid = base_group
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_comm_group ( base_group, base, ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_group @ mp_group" )

  CALL mpi_group_incl ( base, group_size, group_list, newgroup, ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_group_incl @ mp_group" )

  CALL mpi_comm_create ( base_group, newgroup, groupid, ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_create @ mp_group" )
  CALL add_perf(perf_id=2,count=1)
  t_end = m_cputime ( )
  CALL add_perf(perf_id=1,time=t_end-t_start)
#endif

END SUBROUTINE mp_group

!******************************************************************************

END MODULE message_passing

!******************************************************************************
